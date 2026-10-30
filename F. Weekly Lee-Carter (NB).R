###### (1) Main function ----
wLC.NB <- function(xv, isoyv, isowv, etx, dtx, int, constraints, eps, adjmat,
                   initialize = NULL, prev.dlnm = NULL, verbose = FALSE,
                   max.iter = 50, update.dispersion = TRUE){
  
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
  if(is.null(prev.dlnm)){ 
    
    if(is.null(initialize)){
      beta1v   <- matrix(rep(int, each = R), nrow = n, ncol = R)
      beta2v   = beta3v = beta4v <- (1:n)*0
      kappa2v  <- matrix(0, nrow = mY, ncol = R)
      kappa3v  <- matrix(0, nrow = mW, ncol = R)
      
      for(j in 1:n) {
        if(all(int == 0)) 
          beta1v[j,] <- colMeans(log(mtx[,j,, drop = FALSE]), na.rm = TRUE)
        beta2v[j] <- 1
        beta3v[j] <- 1
      }
      
      for(r in 1:R){
        kappa2v[,r] <- sapply(1:mY, function(t) 
          (sum(log(mtx[isoyv==t,,r]), na.rm = TRUE)/mWv[t] - 
             sum(beta1v[,r], na.rm = TRUE)))
        kappa3v[,r] <- sapply(1:mW, function(w){
          N <- sum(isowv == w)
          sum(log(mtx[isowv==w,,r]), na.rm = TRUE)/sum(isowv == w, na.rm = TRUE) - 
            sum(beta1v[,r], na.rm = TRUE)
        })
      }
      
      # Dispersion
      ds.x   <- rep(log(5), n) - mean(rep(log(5), n))
      ds.r   <- rep(log(50), R)
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  
      exp.ds <- aperm(exp.ds, c(3, 1, 2))} else {
        
        # Initialize
        beta1v  <- initialize$beta1
        beta2v  <- initialize$beta2
        beta3v  <- initialize$beta3
        kappa2v <- initialize$kappa2
        kappa3v <- initialize$kappa3
        ds.x    <- initialize$phix
        ds.r    <- initialize$phir
        exp.ds  <- exp(outer(ds.x, ds.r, "+"))
        exp.ds  <- replicate(m, exp.ds, simplify = "array")  
        exp.ds  <- aperm(exp.ds, c(3, 1, 2))
      }
    
    # Empty objects
    betal  <- NULL
    gammal <- NULL
    cb     <- NULL
    lambda <- NULL
  } else {
    
    # Initialization
    beta1v  <- prev.dlnm$beta1
    beta2v  <- prev.dlnm$beta2
    beta3v  <- prev.dlnm$beta3
    kappa2v <- prev.dlnm$kappa2
    kappa3v <- prev.dlnm$kappa3
    betal   <- prev.dlnm$betal
    gammal  <- prev.dlnm$gammal
    ds.x    <- prev.dlnm$phix
    ds.r    <- prev.dlnm$phir
    exp.ds  <- exp(outer(ds.x, ds.r, "+"))
    exp.ds  <- replicate(m, exp.ds, simplify = "array") 
    exp.ds  <- aperm(exp.ds, c(3, 1, 2))
    lambda  <- prev.dlnm$lambda
    cb      <- prev.dlnm$cb
  }    
  
  # Iterate
  l0   <- -100
  l1   <- -99
  iter <- 0
  
  while(abs((l1-l0)/l0) > eps & iter <= max.iter)  {
    iter <- iter + 1
    l0   <- l1
    
    ### (1) Stage 1B2 optimise over the beta2(x) ----
    beta2v <- sapply(seq_len(n), function(j) {
      dv  <- abind::adrop(dtx[,j,,drop = FALSE],2)	
      ev  <- abind::adrop(etx[,j,,drop = FALSE],2)	
      phi <- abind::adrop(exp.ds[,j,,drop = FALSE],2)	
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
    
    kappa2v <- do.call('rbind', lapply(seq_len(mY), function(i) {
      dv  <- dtx[isoyv == i,,,drop=FALSE]
      ev  <- etx[isoyv == i,,,drop=FALSE]
      phi <- exp.ds[isoyv == i,,,drop=FALSE]
      llmax.kappa2(b1 = beta1v,  
                   b2 = beta2v,
                   b3 = beta3v,  
                   bl  = if (!is.null(prev.dlnm)) betal,
                   k2 = kappa2v[i-1*(iter==1 & is.null(initialize) & i!=1),, drop = FALSE], 
                   k3 = kappa3v, 
                   gl  = if (!is.null(prev.dlnm)) gammal,
                   cb  = if (!is.null(prev.dlnm)) lapply(cb, function(x) 
                     x[isoyv == i, , ]),
                   dv = dv, 
                   ev = ev,      
                   phi = phi)
    }))
    
    ### (3) Stage 1B3 optimize over the beta3(x) ----
    beta3v <- sapply(seq_len(n), function(j) {
      dv  <- abind::adrop(dtx[,j,,drop = FALSE],2)	
      ev  <- abind::adrop(etx[,j,,drop = FALSE],2)	
      phi <- abind::adrop(exp.ds[,j,,drop = FALSE],2)	
      llmax.beta3(b1  = beta1v[j,], 
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
    })
    
    ### (4) Stage 1K3 optimise over the kappa3(t,r) ----
    kappa3v <- do.call('rbind', lapply(seq_len(mW), function(i) {
      dv  <- dtx[isowv == i,,,drop=FALSE]	
      ev  <- etx[isowv == i,,,drop=FALSE]
      phi <- exp.ds[isowv == i,,,drop=FALSE]
      llmax.kappa3(b1  = beta1v,      
                   b2  = beta2v, 
                   b3  = beta3v,  
                   bl  = if (!is.null(prev.dlnm)) betal,
                   k2  = kappa2v, 
                   k3  = kappa3v[i,,drop=FALSE], 
                   gl  = if (!is.null(prev.dlnm)) gammal,
                   cb  = if (!is.null(prev.dlnm)) lapply(cb, function(x) 
                     x[isowv == i, , ]),
                   dv  = dv, 
                   ev  = ev,          
                   phi = phi)
    }))
    
    ### (5) Stage 1A  optimise over the beta1(x) ----
    beta1v <- do.call('rbind', lapply(seq_len(n), function(j) {
      dv  <- abind::adrop(dtx[,j,,drop = FALSE],2)	
      ev  <- abind::adrop(etx[,j,,drop = FALSE],2)	
      phi <- abind::adrop(exp.ds[,j,,drop = FALSE],2)	
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
      mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,,drop = FALSE] + 
                               beta3v[j] * kappa3v[isowv,,drop = FALSE] + 
                               term.j, 2, beta1v[j,], '+'))
    }
    
    ### (7) Stage 1F  optimize over the phi(x) ----
    if(update.dispersion){
      ds.x <- llmax.phix(mu   = mhat, dv   = dtx,
                         ev   = etx,  phix = ds.x, 
                         phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))  
    }
    
    ### (6) Stage 1F  optimize over the phi(r) ----
    if(update.dispersion){
      ds.r <- llmax.phir(mu = mhat, dv = dtx,
                         ev = etx, phix = ds.x, phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))   
    }
    
    ### (8) Constraints ----
    if(constraints == "ALL"){
      fac21   <- colMeans(kappa2v) # kappa2v[1,]
      fac31   <- colMeans(kappa3v) # kappa3v[1,]
      fac22   <- sum(beta2v)       # beta2v[1] 
      fac32   <- sum(beta3v)       # beta3v[1]
      fac71   <- mean(ds.x)        # ds.x[1]
      kappa2v <- fac22 * sweep(kappa2v, 2, fac21, '-')   
      kappa3v <- fac32 * sweep(kappa3v, 2, fac31, '-')
      beta2v  <- beta2v/fac22     
      beta3v  <- beta3v/fac32
      beta1v  <- beta1v + 
        (beta2v * fac22) %o% fac21 + 
        (beta3v * fac32) %o% fac31
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
      mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,,drop=FALSE] + 
                               beta3v[j] * kappa3v[isowv,,drop=FALSE] + 
                               term.j, 2, beta1v[j,], '+'))
    }
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat + (etx*mhat)^2/exp.ds)
    l1      <- sum(lgamma(exp.ds + dtx) - lgamma(exp.ds) - 
                     lgamma(dtx + 1) + dtx * log(etx * mhat) - 
                     dtx * log(exp.ds + etx * mhat) + 
                     exp.ds * log(exp.ds) - 
                     exp.ds * log(exp.ds + etx * mhat), na.rm = TRUE)
    pen     <- if(!is.null(prev.dlnm)) 
      dlnm_penalty(lambda, betal, gammal, adjmat) else 0
    l1.up   <- l1
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
       ll.up  = l1.up,   cb      = cb)	
  
}

###### (2) Compute Hessian wLC-NB ----
hessian.wLC.NB <- function(etx, dtx, fit){
  H.beta1v  <- hessian.beta1v(dv = dtx,  ev = etx, lca = fit)
  H.beta2v  <- hessian.beta2v(dv = dtx,  ev = etx, lca = fit)
  H.kappa2v <- hessian.kappa2v(dv = dtx, ev = etx, lca = fit)
  H.beta3v  <- hessian.beta3v(dv = dtx,  ev = etx, lca = fit)
  H.kappa3v <- hessian.kappa3v(dv = dtx, ev = etx, lca = fit)
  H.base    <- cbind(H.beta1v, H.beta2v, H.kappa2v, H.beta3v, H.kappa3v)
  
  colnames(H.base) = rownames(H.base) <- 
    c(outer(dimnames(dtx)[[2]], dimnames(dtx)[[3]], 
            FUN = function(x,y) paste0('Alpha_Age', x, "_Region", y)),
      paste0('Beta2_Age', dimnames(dtx)[[2]]),
      outer(unique(isoyv), dimnames(dtx)[[3]], 
            FUN = function(x,y) paste0('Kappa2_Year', x, "_Region", y)),
      paste0('Beta3_Age', dimnames(dtx)[[2]]),
      outer(unique(isowv), dimnames(dtx)[[3]], 
            FUN = function(x,y) paste0('Kappa3_Week', x, "_Region", y)))
  
  # # Remove some columns due to ID constraints
  # ind1.rm <- c(length(fit$beta1) + 1,
  #              length(fit$beta1) + length(fit$beta2) + (0:(R-1))*max(isoyv) + 1,
  #              length(fit$beta1) + length(fit$beta2) + length(fit$kappa2) + 1,
  #              length(fit$beta1) + length(fit$beta2) + length(fit$kappa2) + 
  #                length(fit$beta3) + (0:(R-1))*max(isowv) + 1)
  # 
  # H.base[-ind1.rm, -ind1.rm]
  
  H.base
  
}

###### (3) Newton Raphson iterations baseline model -----
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
    exp(b2 * k2[isoyv,,drop = FALSE] + b3 * k3[isowv,] + term)
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
    exp(sweep(b3 * k3[isowv,,drop = FALSE] + term, 2, b1, '+'))
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
  k3     <- k3[1:dim(dv)[1],, drop = FALSE]
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
llmax.beta3 <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi){
  b31    <- b3
  b30    <- b31 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    Reduce(`+`, Map(`*`, bl, dlnm)) else 0
  thetat <- ev * 
    exp(sweep(b2 * k2[isoyv,,drop = FALSE] + term, 2, b1, '+'))
  s1     <- sum(dv * k3[isowv,], na.rm = TRUE)
  
  while(abs(b31 - b30) > 0.0001) {
    b30        <- b31
    thetat.exp <- thetat * exp(b30 * k3[isowv,])
    f0         <- sum((dv + phi) * (1 - phi/(phi + thetat.exp)) * 
                        k3[isowv,], na.rm = TRUE) - s1
    df0        <- sum((dv + phi) * (phi * thetat.exp)/((phi + thetat.exp)^2) * 
                        k3[isowv,]^2, na.rm = TRUE)
    b31        <- b30 - f0/df0
  }
  b31
}
llmax.kappa3 <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi){
  sy     <- 1990
  k31    <- k3
  k30    <- k31 - 1
  k2     <- k2[isoyear(as.Date(rownames(dv))) - sy + 1, ,drop=FALSE]
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    aperm(Reduce(`+`, Map(`%o%`, bl, dlnm)), c(2, 1, 3)) else 0
  thetat <- ev * exp(sweep(aperm(b2 %o% k2, c(2, 1, 3)) + term, c(2,3), b1, '+'))
  s1     <- apply(sweep(dv, 2, b3, '*'), 3, sum, na.rm = TRUE)
  
  while(max(abs(k31 - k30)) > 0.0001) {
    k30        <- k31
    thetat.exp <- sweep(thetat, c(2,3), exp(b3 %o% k30), '*')
    f0         <- apply(sweep((dv + phi) * (1 - phi/(phi + thetat.exp)), 2, b3, '*'), 
                        3, sum, na.rm = TRUE) - s1
    df0        <- apply(sweep((dv + phi) * (phi * thetat.exp/((phi + thetat.exp)^2)),
                              2, b3^2, '*'), 3, sum, na.rm = TRUE)
    k31        <- k30 - f0/df0
  }
  k31
}
llmax.phir <- function(mu, dv, ev, phix, phir){
  phir1     <- phir
  phir0     <- phir - 0.1
  exp.phi0  <- exp(outer(phix, phir1, "+"))
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
  exp.phi0  <- exp(outer(phix1, phir, "+"))
  exp.phi0  <- replicate(nrow(dv), exp.phi0, simplify = "array")  
  exp.phi0  <- aperm(exp.phi0, c(3, 1, 2))  
  
  while(max(abs(phix1 - phix0)) > 0.001) {
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

###### (4) Hessian computations baseline model ----
hessian.beta1v <- function(dv, ev, lca){
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array")  # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat    <- lca$mhat
  R       <- dim(dv)[3]
  na      <- dim(dv)[2]
  isoyv   <- isoyear(as.Date(rownames(dv))) - 
    min(isoyear(as.Date(rownames(dv)))) + 1
  isowv   <- isoweek(as.Date(rownames(dv)))
  ny      <- length(unique(isoyv))
  nw      <- length(unique(isowv))
  
  # Define grid
  grid <- expand.grid(Age = 1:na, Region = 1:R)
  
  # Calculation
  H.beta1v <- sapply(1:nrow(grid), function(i){
    # Values
    x = grid$Age[i]
    r = grid$Region[i]
    
    # Second order derivatives
    base <- - (dv[,x,r] + exp.phi[,x,r]) * exp.phi[,x,r] * ev[,x,r] * mhat[,x,r] *
      (ev[,x,r] * mhat[,x,r] + exp.phi[,x,r])^(-2)
    
    v1 <- sapply(1:R, function(s)
      c(rep(0, (x - 1)), sum(base)*(r==s), rep(0, (na - x))))
    v2 <- c(rep(0, (x - 1)), sum(base * lca$kappa2[isoyv,r]), rep(0, (na - x)))
    v3 <- sapply(1:R, function(s) sapply(1:ny, function(t)
      sum(base[isoyv == t] * lca$beta2[x]) * (r==s)))
    v4 <- unname(c(rep(0, (x - 1)), sum(base * lca$kappa3[isowv,r]), rep(0, (na - x))))
    v5 <- sapply(1:R, function(s) sapply(1:nw, function(w)
      sum(base[isowv == w] * lca$beta3[x]) * (r==s)))
    
    
    # Return
    
    c(v1, v2, v3, v4, v5)
  })
  
  H.beta1v
}
hessian.beta2v <- function(dv, ev, lca){
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array")  # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat    <- lca$mhat
  R       <- dim(dv)[3]
  na      <- dim(dv)[2]
  isoyv   <- isoyear(as.Date(rownames(dv))) - 
    min(isoyear(as.Date(rownames(dv)))) + 1
  isowv   <- isoweek(as.Date(rownames(dv)))
  ny      <- length(unique(isoyv))
  nw      <- length(unique(isowv))
  
  # Calculation
  H.beta2v <- sapply(1:ncol(dv), function(x){
    deriv <- (dv[,x,] - (dv[,x,] + exp.phi[,x,]) * 
                (1 - (exp.phi[,x,]/(ev[,x,] * mhat[,x,] + exp.phi[,x,]))))
    base <- - (dv[,x,] + exp.phi[,x,]) * exp.phi[,x,] * ev[,x,] * mhat[,x,] *
      (ev[,x,] * mhat[,x,] + exp.phi[,x,])^(-2)
    base <- sweep(base, c(1,2), lca$kappa2[isoyv,], '*')
    
    v1 <- sapply(1:R, function(r)
      unname(c(rep(0, (x - 1)), sum(base[,r]), rep(0, (na - x)))))
    v2 <- unname(c(rep(0, (x - 1)),
                   sum(sweep(base, c(1,2), lca$kappa2[isoyv,], '*')),
                   rep(0, (na - x))))
    v3 <- t(sapply(1:ny, function(t)
      colSums(base[isoyv == t,] * lca$beta2[x]) + colSums(deriv[isoyv==t,]))) 
    v4 <- unname(c(rep(0, (x - 1)),
                   sum(sweep(base, c(1,2), lca$kappa3[isowv,], '*')),
                   rep(0, (na - x))))
    v5 <- t(sapply(1:nw, function(w)
      colSums(base[isowv == w,] * lca$beta3[x])))
    
    c(v1, v2, v3, v4, v5)
  })
  
  H.beta2v
}
hessian.kappa2v <- function(dv, ev, lca){
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array")  # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat    <- lca$mhat
  R       <- dim(dv)[3]
  na      <- dim(dv)[2]
  isoyv   <- isoyear(as.Date(rownames(dv))) - 
    min(isoyear(as.Date(rownames(dv)))) + 1
  isowv   <- isoweek(as.Date(rownames(dv)))
  ny      <- length(unique(isoyv))
  nw      <- length(unique(isowv))
  
  # Define grid
  grid <- expand.grid(Year = 1:length(unique(isoyv)), Region = 1:R)
  
  # Calculation 
  H.kappa2v <- sapply(1:nrow(grid), function(j) {
    # Grid values
    t <- grid$Year[j]
    r <- grid$Region[j]
    
    # First order
    deriv <- (dv[isoyv==t,,r] - (dv[isoyv==t,,r] + exp.phi[isoyv==t,,r]) * 
                (1 - exp.phi[isoyv==t,,r]/(ev[isoyv==t,,r] * mhat[isoyv==t,,r] + 
                                             exp.phi[isoyv==t,,r])))
    
    # Second order derivatives
    base <- - (dv[isoyv==t,,r] + exp.phi[isoyv==t,,r]) * exp.phi[isoyv==t,,r] *
      ev[isoyv==t,,r] * mhat[isoyv==t,,r] *
      (ev[isoyv==t,,r] * mhat[isoyv==t,,r] + exp.phi[isoyv==t,,r])^(-2)
    base <- sweep(base, 2, lca$beta2, '*')
    
    v1 <- sapply(1:R, function(s)
      colSums(base)*(r==s))
    v2 <- colSums(base * lca$kappa2[t, r]) + colSums(deriv)
    v3 <- sapply(1:R, function(s)
      unname(c(rep(0, (t - 1)),
               sum(sweep(base, 2, lca$beta2, '*'))*(r==s),
               rep(0, (nrow(lca$kappa2) - t)))))
    v4 <- colSums(sweep(base, 1, lca$kappa3[1:nrow(base),r], '*'))
    v5 <- sapply(1:R, function(s)
      c(rowSums(sweep(base, 2, lca$beta3, '*'))*(r==s), rep(0, nw - nrow(base))))
    
    # Return
    c(v1, v2, v3, v4, v5)
  })
  
  H.kappa2v
}
hessian.beta3v <- function(dv, ev, lca){
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array")  # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat    <- lca$mhat
  R       <- dim(dv)[3]
  na      <- dim(dv)[2]
  isoyv   <- isoyear(as.Date(rownames(dv))) - 
    min(isoyear(as.Date(rownames(dv)))) + 1
  isowv   <- isoweek(as.Date(rownames(dv)))
  ny      <- length(unique(isoyv))
  nw      <- length(unique(isowv))
  
  # Calculation
  H.beta3v <- sapply(1:ncol(dv), function(x){
    
    # First order
    deriv <- (dv[,x,] - (dv[,x,] + exp.phi[,x,]) * 
                (1 - exp.phi[,x,]/(ev[,x,] * mhat[,x,] + exp.phi[,x,])))
    
    # Second order derivatives
    base <- - (dv[,x,] + exp.phi[,x,]) * exp.phi[,x,] * ev[,x,] * mhat[,x,] *
      (ev[,x,] * mhat[,x,] + exp.phi[,x,])^(-2)
    base <- sweep(base, c(1,2), lca$kappa3[isowv,], '*')
    
    v1 <- sapply(1:R, function(r)
      unname(c(rep(0, (x - 1)), sum(base[,r]), rep(0, (na - x)))))
    v2 <- unname(c(rep(0, (x - 1)),
                   sum(sweep(base, c(1,2), lca$kappa2[isoyv,], '*')),
                   rep(0, (na - x))))
    v3 <- t(sapply(1:nrow(lca$kappa2), function(t)
      colSums(base[isoyv == t,] * lca$beta2[x])))
    v4 <- unname(c(rep(0, (x - 1)),
                   sum(sweep(base, c(1,2), lca$kappa3[isowv,], '*')),
                   rep(0, (na - x))))
    v5 <- t(sapply(1:nrow(lca$kappa3), function(w)
      colSums(base[isowv == w,] * lca$beta3[x])  + colSums(deriv[isowv==w,])))
    
    # Return
    c(v1, v2, v3, v4, v5)
  })
  
  H.beta3v
}
hessian.kappa3v <- function(dv, ev, lca){
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array")  # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat    <- lca$mhat
  R       <- dim(dv)[3]
  na      <- dim(dv)[2]
  isoyv   <- isoyear(as.Date(rownames(dv))) - 
    min(isoyear(as.Date(rownames(dv)))) + 1
  isowv   <- isoweek(as.Date(rownames(dv)))
  ny      <- length(unique(isoyv))
  nw      <- length(unique(isowv))
  
  # Define grid
  grid <- expand.grid(Week = 1:length(unique(isowv)), Region = 1:R)
  
  # Calculation 
  H.kappa3v <- sapply(1:nrow(grid), function(j) {  
    # Grid values
    w <- grid$Week[j]
    r <- grid$Region[j]
    
    # First order
    deriv <- (dv[isowv==w,,r] - (dv[isowv==w,,r] + exp.phi[isowv==w,,r]) * 
                (1 - exp.phi[isowv==w,,r]/(ev[isowv==w,,r] * mhat[isowv==w,,r] + 
                                             exp.phi[isowv==w,,r])))
    
    # Calculation
    base <- - (dv[isowv==w,,r] + exp.phi[isowv==w,,r]) * exp.phi[isowv==w,,r] *
      ev[isowv==w,,r] * mhat[isowv==w,,r] *
      (ev[isowv==w,,r] * mhat[isowv==w,,r] + exp.phi[isowv==w,,r])^(-2)
    base <- sweep(base, 2, lca$beta3, '*')
    
    v1 <- sapply(1:R, function(s)
      colSums(base)*(r==s))
    v2 <- colSums(sweep(base, 1,
                        lca$kappa2[isoyear(as.Date(rownames(base))) -
                                     min(isoyear(as.Date(rownames(dv)))) + 1,r], '*'))
    vv3 <- sapply(1:R, function(s)
      c(rowSums(sweep(base, 2, lca$beta2, '*'))*(r==s)))
    dimnames(vv3) <- list(isoyear(as.Date(rownames(vv3))) -
                            min(isoyear(as.Date(rownames(dv)))) + 1,
                          dimnames(dv)[[3]])
    
    v3 <- matrix(0, nrow = length(unique(isoyv)), ncol = R,
                 dimnames = list(unique(isoyv), dimnames(dv)[[3]]))
    v3[rownames(vv3),colnames(vv3)] <- vv3
    
    v4 <- colSums(base * lca$kappa3[w, r]) + colSums(deriv)
    v5 <- sapply(1:R, function(s)
      unname(c(rep(0, (w - 1)),
               sum(sweep(base, 2, lca$beta3, '*'))*(r==s),
               rep(0, (nrow(lca$kappa3) - w)))))
    
    # Return
    c(v1, v2, v3, v4, v5)
  })
}

###### (5) Contrast function ----
contrast_base <- function(fit){
  nA <- length(fit$beta2)
  R  <- ncol(fit$kappa2)
  nT <- nrow(fit$kappa2)
  nW <- nrow(fit$kappa3)
  as.matrix(Matrix::bdiag(diag(nA * R), 
                          sumzero(nA), 
                          kronecker(diag(R), sumzero(nT)),
                          sumzero(nA), 
                          kronecker(diag(R), sumzero(nW))))
}