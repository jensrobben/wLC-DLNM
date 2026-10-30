###### (1) Main function ----

aLL.NB <- function(xv, isoyv, isowv, etx, dtx, eps, verbose = FALSE, 
                   max.iter = 50, update.dispersion = TRUE){
  
  # Notation
  n   <- length(xv)	
  mY  <- length(unique(isoyv))
  m   <- nrow(etx)
  R   <- dim(dtx)[3]
  
  ### Stage 1: common trend, fitted on the aggregated data -----
  dtx.agg <- array(rowSums(dtx, dims = 2), dim = c(m, n, 1))
  etx.agg <- array(rowSums(etx, dims = 2), dim = c(m, n, 1))
  mtx.agg <- dtx.agg/etx.agg
  mtx.agg[mtx.agg == 0] <- 10^(-8)
  
  Av <- matrix(apply(log(mtx.agg), c(2, 3), mean, na.rm = TRUE), n, 1)
  Bv <- rep(1/n, n)
  Kv <- sapply(seq_len(mY), function(t) 
    sum(apply(log(mtx.agg[isoyv == t,,,drop = FALSE]), c(2, 3), mean) - Av))
  ds.x.agg <- rep(log(5), n) - mean(rep(log(5), n))
  ds.r.agg <- log(50)
  exp.agg  <- aperm(replicate(m, exp(outer(ds.x.agg, ds.r.agg, "+")), 
                              simplify = "array"), c(3, 1, 2))
  
  l0    <- -100
  l1    <- -99
  iter1 <- 0
  
  while(abs((l1-l0)/l0) > eps & iter1 <= max.iter)  {
    iter1 <- iter1 + 1
    l0    <- l1
    
    Bv <- llmax.B.LL(a = Av, B = Bv, K = Kv, dv = dtx.agg, 
                     ev = etx.agg, phi = exp.agg)
    Kv <- llmax.K.LL(a = Av, B = Bv, K = Kv, dv = dtx.agg, 
                     ev = etx.agg, phi = exp.agg)
    Av <- llmax.alpha.LL(a = Av, B = Bv, K = Kv, dv = dtx.agg, 
                         ev = etx.agg, phi = exp.agg)
    
    mhat.agg <- exp(lp.LL(Av, B = Bv, K = Kv))
    
    if(update.dispersion){
      ds.r.agg <- llmax.phir.LL(mu = mhat.agg, dv = dtx.agg, ev = etx.agg, 
                                phix = ds.x.agg, phir = ds.r.agg)
      ds.x.agg <- llmax.phix.LL(mu = mhat.agg, dv = dtx.agg, ev = etx.agg, 
                                phix = ds.x.agg, phir = ds.r.agg)
    }
    
    # Constraints: sum_x B_x^2 = 1, sum_t K_t = 0
    fac1 <- sum(Bv)
    Bv   <- Bv/fac1
    Kv   <- Kv*fac1
    
    fac2 <- mean(Kv)
    Kv   <- Kv - fac2
    Av   <- Av + Bv*fac2
    
    fac7     <- mean(ds.x.agg)
    ds.r.agg <- ds.r.agg + fac7
    ds.x.agg <- ds.x.agg - fac7
    exp.agg  <- aperm(replicate(m, exp(outer(ds.x.agg, ds.r.agg, "+")), simplify = "array"), c(3, 1, 2))
    
    mhat.agg <- exp(lp.LL(Av, B = Bv, K = Kv))
    
    l1 <- sum(lgamma(exp.agg + dtx.agg) - lgamma(exp.agg) - lgamma(dtx.agg + 1) + 
                dtx.agg * log(etx.agg * mhat.agg) - 
                dtx.agg * log(exp.agg + etx.agg * mhat.agg) + 
                exp.agg * log(exp.agg) - 
                exp.agg * log(exp.agg + etx.agg * mhat.agg), na.rm = TRUE)
    
    if(verbose) cat(l1, ' -> ')
  }
  
  Av <- as.vector(Av)
  if(verbose) cat('\n  Stage 1 (common trend) LogLik:', l1, '\n')
  
  ### Stage 2: regional deviation, with exposures adjusted by the common trend -----
  off     <- array(rep(outer(Kv[isoyv], Bv) + rep(Av, each = m), R), 
                   dim = c(m, n, R))
  mtx     <- dtx/etx
  mtx[mtx == 0] <- 10^(-8)
  
  beta1v  <- apply(log(mtx) - off, c(2, 3), mean, na.rm = TRUE)    
  beta2v  <- matrix(1/n, nrow = n, ncol = R)                      
  kappa2v <- matrix(0, nrow = mY, ncol = R)                       
  
  ds.x   <- rep(log(5), n) - mean(rep(log(5), n))
  ds.r   <- rep(log(50), R)
  exp.ds <- aperm(replicate(m, exp(outer(ds.x, ds.r, "+")), simplify = "array"), 
                  c(3, 1, 2))
  
  l0    <- -100
  l1    <- -99
  iter2 <- 0
  
  while(abs((l1-l0)/l0) > eps & iter2 <= max.iter)  {
    iter2 <- iter2 + 1
    l0    <- l1
    
    kappa2v <- llmax.kappa.LL(a = beta1v, b = beta2v, k = kappa2v,
                              dv = dtx, ev = etx, phi = exp.ds, off = off)
    beta2v  <- llmax.beta.LL(a = beta1v, b = beta2v, k = kappa2v,
                             dv = dtx, ev = etx, phi = exp.ds, off = off)
    beta1v  <- llmax.alpha.LL(a = beta1v, b = beta2v, k = kappa2v,
                              dv = dtx, ev = etx, phi = exp.ds, off = off)
    
    mhat <- exp(lp.LL(beta1v, b = beta2v, k = kappa2v, off = off))
    
    if(update.dispersion){
      ds.r   <- llmax.phir.LL(mu = mhat, dv = dtx, ev = etx, phix = ds.x, phir = ds.r)
      ds.x   <- llmax.phix.LL(mu = mhat, dv = dtx, ev = etx, phix = ds.x, phir = ds.r)
    }
    
    # Constraints: sum_x b_{x,r}^2 = 1, sum_t kappa_{t,r} = 0, sum_x phi_x = 0
    fac3    <- colSums(beta2v)
    beta2v  <- sweep(beta2v, 2, fac3, '/')
    kappa2v <- sweep(kappa2v, 2, fac3, '*')
    
    fac4    <- colMeans(kappa2v)
    kappa2v <- sweep(kappa2v, 2, fac4, '-')
    beta1v  <- beta1v + sweep(beta2v, 2, fac4, '*')
    
    fac7 <- mean(ds.x)
    ds.r <- ds.r + fac7
    ds.x <- ds.x - fac7
    exp.ds  <- aperm(replicate(m, exp(outer(ds.x, ds.r, "+")), 
                               simplify = "array"), c(3, 1, 2))
    
    mhat <- exp(lp.LL(beta1v, b = beta2v, k = kappa2v, off = off))
    dimnames(mhat) <- dimnames(dtx)
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat + (etx*mhat)^2/exp.ds)
    l1 <- sum(lgamma(exp.ds + dtx) - lgamma(exp.ds) - lgamma(dtx + 1) + 
                dtx * log(etx * mhat) - dtx * log(exp.ds + etx * mhat) + 
                exp.ds * log(exp.ds) - exp.ds * log(exp.ds + etx * mhat), na.rm = TRUE)
    if(verbose) cat(l1, ' -> ')
  }
  if(verbose) cat('\n  Stage 2 (regional deviation) LogLik:', l1, '\n')
  
  npar <- length(Av) + length(Bv) + length(Kv) + length(beta1v) + length(beta2v) + 
    length(kappa2v) + length(ds.x) + length(ds.r) - (2 + 2*R + 1)
  
  # Output: beta1 collects A_x + a_{x,r}, so that log mu = beta1 + B K + beta2 kappa
  list(A      = Av,      B      = Bv,      K      = Kv,
       a      = beta1v,  beta1  = sweep(beta1v, 1, Av, '+'),
       beta2  = beta2v,  kappa2 = kappa2v, 
       phix   = ds.x,    phir   = ds.r, 
       phi.common = list(phix = ds.x.agg, phir = ds.r.agg),   
       epsilon = epsilon, mhat  = mhat,  
       ll     = l1,      npar   = npar,
       itr    = c('stage1' = iter1, 'stage2' = iter2),
       cvg    = c(iter1, iter2) <= max.iter)	
}

###### (2) Linear predictor ----
lp.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, off = NULL){
  m   <- length(isoyv)
  out <- array(0, dim = c(m, nrow(a), ncol(a)))
  for(r in seq_len(ncol(a))){
    out[,,r] <- rep(a[, r], each = m)
    if(!is.null(B)) out[,,r] <- out[,,r] + outer(K[isoyv], B)
    if(!is.null(b)) out[,,r] <- out[,,r] + outer(k[isoyv, r], b[, r])
  }
  if(!is.null(off)) out <- out + off
  out
}

###### (3) Newton Raphson iterations -----
nb.grad.LL <- function(dv, ev, phi, lp){
  thetat.exp <- ev * exp(lp)
  list(f = (dv + phi) * (1 - phi/(phi + thetat.exp)) - dv,          
       h = (dv + phi) * (phi * thetat.exp)/((phi + thetat.exp)^2))   
}
llmax.alpha.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, dv, ev, phi, off = NULL){
  a1 <- a
  a0 <- a1 - 1
  while(max(abs(a1 - a0)) > 0.0001){
    a0  <- a1
    g   <- nb.grad.LL(dv, ev, phi, lp.LL(a0, B, K, b, k, off))
    f0  <- colSums(g$f, dims = 1, na.rm = TRUE)
    df0 <- colSums(g$h, dims = 1, na.rm = TRUE)
    a1  <- a0 - f0/df0
  }
  a1
}
llmax.B.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, dv, ev, phi, off = NULL){
  B1 <- B
  B0 <- B1 - 1
  z  <- K[isoyv]
  while(max(abs(B1 - B0)) > 0.0001){
    B0  <- B1
    g   <- nb.grad.LL(dv, ev, phi, lp.LL(a, B0, K, b, k, off))
    f0  <- rowSums(colSums(g$f * z,   dims = 1, na.rm = TRUE))
    df0 <- rowSums(colSums(g$h * z^2, dims = 1, na.rm = TRUE))
    B1  <- B0 - f0/df0
  }
  B1
}
llmax.K.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, dv, ev, phi, off = NULL){
  K1 <- K
  K0 <- K1 - 1
  while(max(abs(K1 - K0)) > 0.0001) {
    K0  <- K1
    g   <- nb.grad.LL(dv, ev, phi, lp.LL(a, B, K0, b, k, off))
    f0  <- rowsum(rowSums(sweep(g$f, 2, B,   '*'), dims = 1, na.rm = TRUE), isoyv)[,1]
    df0 <- rowsum(rowSums(sweep(g$h, 2, B^2, '*'), dims = 1, na.rm = TRUE), isoyv)[,1]
    K1  <- K0 - f0/df0
  }
  K1
}
llmax.beta.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, dv, ev, phi, off = NULL){
  b1 <- b
  b0 <- b1 - 1
  z  <- k[isoyv, , drop = FALSE]
  while(max(abs(b1 - b0)) > 0.0001){
    b0  <- b1
    g   <- nb.grad.LL(dv, ev, phi, lp.LL(a, B, K, b0, k, off))
    f0  <- colSums(sweep(g$f, c(1, 3), z,   '*'), dims = 1, na.rm = TRUE)
    df0 <- colSums(sweep(g$h, c(1, 3), z^2, '*'), dims = 1, na.rm = TRUE)
    b1  <- b0 - f0/df0
  }
  b1
}
llmax.kappa.LL <- function(a, B = NULL, K = NULL, b = NULL, k = NULL, dv, ev, phi, off = NULL){
  k1 <- k
  k0 <- k1 - 1
  while(max(abs(k1 - k0)) > 0.0001) {
    k0  <- k1
    g   <- nb.grad.LL(dv, ev, phi, lp.LL(a, B, K, b, k0, off))
    f0  <- rowsum(apply(sweep(g$f, c(2, 3), b,   '*'), c(1, 3), sum, na.rm = TRUE), isoyv)
    df0 <- rowsum(apply(sweep(g$h, c(2, 3), b^2, '*'), c(1, 3), sum, na.rm = TRUE), isoyv)
    k1  <- k0 - f0/df0
  }
  unname(k1)
}
llmax.phir.LL <- function(mu, dv, ev, phix, phir){
  phir1     <- phir
  phir0     <- phir - 0.1
  it        <- 0
  
  while(max(abs(phir1 - phir0)) > 0.001 & it < 100) {
    it        <- it + 1
    phir0     <- phir1
    exp.phi0  <- exp(outer(phix, phir0, "+"))
    exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  
    exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))
    
    dg1 <- digamma(dv + exp.phi0)
    tg1 <- trigamma(dv + exp.phi0)
    dg2 <- digamma(exp.phi0)
    tg2 <- trigamma(exp.phi0)
    c1  <- dv + exp.phi0
    c2  <- ev * mu + exp.phi0
    lc2 <- log(c2)
    
    f0  <- apply((dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
                 3, sum, na.rm = TRUE)
    df0 <- apply(((tg1 - tg2 + 1/exp.phi0 - 2/c2 + c1 * 1/(c2^2)) * exp.phi0^2) +
                   (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0, 
                 3, sum, na.rm = TRUE)
    
    # Safeguarded Newton step: Newton where the log-likelihood is concave in phi,
    # a small ascent step elsewhere, and at most 1 per iteration
    step  <- ifelse(df0 < 0, f0/df0, -0.1 * sign(f0))
    step  <- pmax(pmin(step, 1), -1)
    phir1 <- phir0 - step
  }
  phir1
}
llmax.phix.LL <- function(mu, dv, ev, phix, phir){
  phix1     <- phix
  phix0     <- phix - 0.1
  it        <- 0
  
  while(max(abs(phix1 - phix0)) > 0.001 & it < 100) {
    it        <- it + 1
    phix0     <- phix1
    exp.phi0  <- exp(outer(phix0, phir, "+"))
    exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  
    exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))
    
    dg1 <- digamma(dv + exp.phi0)
    tg1 <- trigamma(dv + exp.phi0)
    dg2 <- digamma(exp.phi0)
    tg2 <- trigamma(exp.phi0)
    c1  <- dv + exp.phi0
    c2  <- ev * mu + exp.phi0
    lc2 <- log(c2)
    
    f0  <- apply((dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
                 2, sum, na.rm = TRUE)
    df0 <- apply(((tg1 - tg2 + 1/exp.phi0 - 2/c2 + c1 * 1/(c2^2)) * exp.phi0^2) +
                   (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0, 
                 2, sum, na.rm = TRUE)
    
    # Safeguarded Newton step: Newton where the log-likelihood is concave in phi,
    # a small ascent step elsewhere, and at most 1 per iteration
    step  <- ifelse(df0 < 0, f0/df0, -0.1 * sign(f0))
    step  <- pmax(pmin(step, 1), -1)
    phix1 <- phix0 - step
  }
  phix1
}