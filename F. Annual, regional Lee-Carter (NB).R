###### (1) Main function ----
aLC.NB <- function(xv, isoyv, isowv, etx, dtx, int = xv * 0, constraints = 'ALL', 
                   eps, adjmat, initialize = NULL, prev.dlnm = NULL, 
                   verbose = FALSE, max.iter = 50, update.dispersion = TRUE){
  
  # Notation
  n   <- length(xv)	
  mY  <- length(unique(isoyv))
  mW  <- length(unique(isowv)) 
  mWv <- sapply(split(isowv, isoyv), max)
  m   <- nrow(etx)
  R   <- dim(dtx)[3]
  
  # Death rates
  mtx <- dtx/etx	    
  mtx[mtx==0] <- 10^(-8)
  
  # Initialization
  beta1v   <- matrix(rep(int, each = R), nrow = n, ncol = R)
  beta2v   = beta3v <- (1:n)*0
  kappa2v  <- matrix(0, nrow = mY, ncol = R)
  kappa3v  <- matrix(0, nrow = mW, ncol = R)
      
  for(j in 1:n) {
    if(all(int == 0)) 
      beta1v[j,] <- colMeans(log(mtx[,j,]), na.rm = TRUE)
    beta2v[j] <- 1
  }
      
  for(r in 1:R){
    kappa2v[,r] <- sapply(1:mY, function(t) 
      (sum(log(mtx[isoyv==t,,r]), na.rm = TRUE)/mWv[t] - 
         sum(beta1v[,r], na.rm = TRUE)))
  }
      
  # Dispersion
  ds.x   <- rep(log(5), n) - mean(rep(log(5), n))
  ds.r   <- rep(log(50), R)
  exp.ds <- exp(outer(ds.x, ds.r, "+"))
  exp.ds <- replicate(m, exp.ds, simplify = "array")  
  exp.ds <- aperm(exp.ds, c(3, 1, 2))
    
  # Empty objects
  betal  <- NULL
  gammal <- NULL
  cb     <- NULL
  lambda <- NULL
    
  # Iterate
  l0   <- -100
  l1   <- -99
  iter <- 0
  
  while(abs((l1-l0)/l0) > eps & iter <= max.iter)  {
    iter <- iter + 1
    l0   <- l1
    
    ### (1) Stage 1B2 optimise over the beta2(x) ----
    beta2v <- sapply(seq_len(n), function(j) {
      dv  <- dtx[,j,]	
      ev  <- etx[,j,]	
      phi <- exp.ds[,j,]
      llmax.beta2(b1  = beta1v[j,], 
                  b2  = beta2v[j-1*(iter==1 & is.null(initialize) & j!=1)],
                  b3  = beta3v[j],  
                  bl  = if (!is.null(prev.dlnm)) sapply(betal, `[[`, j),
                  k2  = kappa2v,    
                  k3  = kappa3v,    
                  gl  = if (!is.null(prev.dlnm)) gammal,
                  cb  = if (!is.null(prev.dlnm)) cb,
                  dv  = dv,         
                  ev  = ev,         
                  phi = phi)
    })
    
    ### (2) Stage 1K2 optimise over the kappa2(t,r) ----
    sy      <- isoyear(as.Date(rownames(etx)[1]))
    
    kappa2v <- t(sapply(seq_len(mY), function(i) {
      dv  <- dtx[isoyv == i,,]
      ev  <- etx[isoyv == i,,]
      phi <- exp.ds[isoyv == i,,]
      llmax.kappa2(b1 = beta1v,  
                   b2 = beta2v,
                   b3 = beta3v,  
                   bl  = if (!is.null(prev.dlnm)) betal,
                   k2 = kappa2v[i-1*(iter==1 & is.null(initialize) & i!=1),], 
                   k3 = kappa3v, 
                   gl  = if (!is.null(prev.dlnm)) gammal,
                   cb  = if (!is.null(prev.dlnm)) lapply(cb, function(x) 
                     x[isoyv == i, , ]),
                   dv = dv, 
                   ev = ev,      
                   phi = phi)
    }))
    
    ### (3) Stage 1A  optimise over the beta1(x) ----
    beta1v <- t(sapply(seq_len(n), function(j) {
      dv  <- dtx[,j,]	
      ev  <- etx[,j,]	
      phi <- exp.ds[,j,]
      llmax.beta1(b1  = beta1v[j,], 
                  b2  = beta2v[j],
                  b3  = beta3v[j],  
                  bl  = if (!is.null(prev.dlnm)) sapply(betal, `[[`, j),
                  k2  = kappa2v,    
                  k3  = kappa3v,  
                  gl  = if (!is.null(prev.dlnm)) gammal,
                  cb  = if (!is.null(prev.dlnm)) cb,
                  dv  = dv,         
                  ev  = ev,         
                  phi = phi)
    }))
    
    mhat <- mtx * NA
    dlnm <- Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, cb, gammal)
    
    for(j in 1:n){
      term.j    <- if (length(betal) > 0) 
        Reduce(`+`, Map(`*`, sapply(betal, `[[`, j), dlnm)) else 0
      mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] + 
                               beta3v[j] * kappa3v[isowv,] + term.j, 2, beta1v[j,], '+'))
    }
    
    ### (4) Stage 1F  optimize over the phi(r) ----
    if(update.dispersion){
      ds.r <- llmax.phir(mu = mhat, dv = dtx,
                         ev = etx, phix = ds.x, phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))   
    }
    
    ### (5) Stage 1F  optimize over the phi(x) ----
    if(update.dispersion){
      ds.x <- llmax.phix(mu   = mhat[,,], dv   = dtx[,,],
                         ev   = etx[,,],  phix = ds.x[], 
                         phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))  
    }
    
    ### (6) Constraints ----
    if(constraints == "ALL"){
      fac21   <- colMeans(kappa2v)
      fac22   <- sum(beta2v)
      fac71   <- mean(ds.x)
      kappa2v <- fac22 * sweep(kappa2v, 2, fac21, '-')   
      beta2v  <- beta2v/fac22     
      beta1v  <- beta1v + 
        (beta2v * fac22) %o% fac21  
      ds.r    <- ds.r + fac71
      ds.x    <- ds.x - fac71
    }
    
    exp.ds <- exp(outer(ds.x, ds.r, "+"))
    exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
    exp.ds <- aperm(exp.ds, c(3, 1, 2))  
    
    ### (9) Log-likelihood ----
    mhat <- mtx * NA
    dlnm <- Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, cb, gammal)
    
    for(j in 1:n){
      term.j    <- if (length(betal) > 0) 
        Reduce(`+`, Map(`*`, sapply(betal, `[[`, j), dlnm)) else 0
      mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] + 
                               beta3v[j] * kappa3v[isowv,] + term.j, 2, beta1v[j,], '+'))
    }
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat + (etx*mhat)^2/exp.ds)
    l1      <- sum(lgamma(exp.ds + dtx) - lgamma(exp.ds) - 
                     lgamma(dtx + 1) + dtx * log(etx * mhat) - 
                     dtx * log(exp.ds + etx * mhat) + 
                     exp.ds * log(exp.ds) - 
                     exp.ds * log(exp.ds + etx * mhat), na.rm = TRUE)
    pen     <- if(!is.null(prev.dlnm)) 
      Reduce(`+`, Map(function(lam, gam)
        lam/2 * sum(diag(gam %*% adjmat %*% t(gam))),
        lambda, gammal)) else 0
    l1.pen  <- l1 - pen
    
    l1      <- l1.pen
    
    if(verbose) cat(l1, ' -> ')
  }	

  # Output prev dlnm 
  list(beta1  = beta1v,  beta2   = beta2v,  beta3 = beta3v,
       kappa2 = kappa2v, kappa3  = kappa3v, 
       betal  = betal,   gammal  = gammal,
       phix   = ds.x,    phir    = ds.r, 
       lambda = lambda,  epsilon = epsilon, 
       mhat   = mhat,    ll      = l1,
       cb     = cb)	
}

###### (2) Newton Raphson iterations baseline model -----
llmax.beta1 <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi){
  b11    <- b1
  b10    <- b11 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    Reduce(`+`, Map(`*`, bl, dlnm)) else 0
  thetat <- ev * 
    exp(b2 * k2[isoyv,] + b3 * k3[isowv,] + term)
  s1     <- colSums(dv, na.rm = TRUE)
  
  while(max(abs(b11 - b10)) > 0.0001){
    b10        <- b11
    thetat.exp <- sweep(thetat, 2, exp(b10), '*')
    f0         <- colSums((dv + phi) * (1 - phi/(phi + thetat.exp)), 
                          na.rm = TRUE) - s1
    df0        <- colSums((dv + phi) * (phi * thetat.exp)/((phi + thetat.exp)^2), 
                          na.rm = TRUE)
    b11        <- b10 - f0/df0
    b11
  }
  b11
}
llmax.beta2 <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi){
  b21    <- b2
  b20    <- b21 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    Reduce(`+`, Map(`*`, bl, dlnm)) else 0
  thetat <- ev * 
    exp(sweep(b3 * k3[isowv,] + term, 2, b1, '+'))
  s1     <- sum(dv * k2[isoyv,], na.rm = TRUE)
  
  while(abs(b21 - b20) > 0.0001){
    b20        <- b21
    thetat.exp <- thetat * exp(b20 * k2[isoyv,])
    f0         <- sum((dv + phi) * (1 - phi/(phi + thetat.exp)) * 
                        k2[isoyv,], na.rm = TRUE) - s1
    df0        <- sum((dv + phi) * (phi * thetat.exp)/((phi + thetat.exp)^2) * 
                        k2[isoyv,]^2, na.rm = TRUE)
    b21        <- b20 - f0/df0
    b21
  }
  b21
}
llmax.kappa2  <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi){
  k21    <- k2
  k20    <- k21 - 1
  k3     <- k3[1:dim(dv)[1],]
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    aperm(Reduce(`+`, Map(`%o%`, bl, dlnm)), c(2, 1, 3)) else 0
  thetat <- ev * exp(sweep(aperm(b3 %o% k3, c(2, 1, 3)) + term, c(2,3), b1, '+'))
  s1     <- apply(sweep(dv, 2, b2, '*'), 3, sum, na.rm = TRUE)
  
  while(max(abs(k21 - k20)) > 0.0001) {
    k20        <- k21
    thetat.exp <- sweep(thetat, c(2,3), exp(b2 %o% k20), '*')
    f0         <- apply(sweep((dv + phi) * (1 - phi/(phi + thetat.exp)), 2, b2, '*'), 
                        3, sum, na.rm = TRUE) - s1
    df0        <- apply(sweep((dv + phi) * (phi * thetat.exp/((phi + thetat.exp)^2)),
                              2, b2^2, '*'), 3, sum, na.rm = TRUE)
    k21        <- k20 - f0/df0
    k21
  }
  k21
} 
llmax.phir <- function(mu, dv, ev, phix, phir){
  phir1     <- phir
  phir0     <- phir - 0.1
  exp.phi0  <- exp(outer(phix, phir0, "+"))
  exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  # gives an n x R x T array
  exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))  
  
  while(max(abs(phir1 - phir0)) > 0.001) {
    phir0     <- phir1
    exp.phi0  <- exp(outer(phix, phir0, "+"))
    exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  # gives an n x R x T array
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
    
    phir1 <- phir0 - f0/df0
    
    phir1
  }
  phir1
}
llmax.phix <- function(mu, dv, ev, phix, phir){
  phix1     <- phix
  phix0     <- phix - 0.1
  exp.phi0  <- exp(outer(phix0, phir, "+"))
  exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  # gives an n x R x T array
  exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))  
  
  while(max(abs(phix1 - phix0)) > 0.001) {
    phix0     <- phix1
    exp.phi0  <- exp(outer(phix0, phir, "+"))
    exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  # gives an n x R x T array
    exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))
    
    dg1 <- digamma(dv + exp.phi0)
    tg1 <- trigamma(dv + exp.phi0)
    dg2 <- digamma(exp.phi0)
    tg2 <- trigamma(exp.phi0)
    c1  <- dv + exp.phi0
    c2  <- ev * mu + exp.phi0
    lc2 <- log(c2)
    
    f0  <- apply((dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
                 2, sum)
    df0 <- apply(((tg1 - tg2 + 1/exp.phi0 - 2/c2 + c1 * 1/(c2^2)) * exp.phi0^2) +
                   (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0, 
                 2, sum)
    
    phix1 <- phix0 - f0/df0
    
    phix1
  }
  # fit <- lm(phix1[1:7] ~ I(1:(length(phix1) - 1)))
  # phix1[length(phix1)] <- fit$coefficients[1] + fit$coefficients[2]*length(phix1)
  phix1
}
