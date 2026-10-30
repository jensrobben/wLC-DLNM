###### (1) Main function ----
wLC.POI <- function(xv, isoyv, isowv, etx, dtx, int, constraints, eps, adjmat,
                    initialize = NULL, prev.dlnm = NULL, verbose = FALSE,
                    max.iter = 50){
  
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
        })}} else {
        
        # Initialize
        beta1v  <- initialize$beta1
        beta2v  <- initialize$beta2
        beta3v  <- initialize$beta3
        kappa2v <- initialize$kappa2
        kappa3v <- initialize$kappa3
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
      llmax.beta2p(b1  = beta1v[j,], 
                   b2  = beta2v[j-1*(iter==1 & is.null(initialize) & j!=1)],
                   b3  = beta3v[j],  
                   bl  = if (!is.null(prev.dlnm)) sapply(betal, `[[`, j),
                   k2  = kappa2v,    
                   k3  = kappa3v,    
                   gl  = if (!is.null(prev.dlnm)) gammal,
                   cb  = if (!is.null(prev.dlnm)) cb,
                   dv  = dv,         
                   ev  = ev)
    })
    
    ### (2) Stage 1K2 optimise over the kappa2(t,r) ----
    sy      <- isoyear(as.Date(rownames(etx)[1]))
    
    kappa2v <- do.call('rbind', lapply(seq_len(mY), function(i) {
      dv  <- dtx[isoyv == i,,,drop=FALSE]
      ev  <- etx[isoyv == i,,,drop=FALSE]
      llmax.kappa2p(b1 = beta1v,  
                    b2 = beta2v,
                    b3 = beta3v,  
                    bl  = if (!is.null(prev.dlnm)) betal,
                    k2 = kappa2v[i-1*(iter==1 & is.null(initialize) & i!=1),,drop = FALSE], 
                    k3 = kappa3v, 
                    gl  = if (!is.null(prev.dlnm)) gammal,
                    cb  = if (!is.null(prev.dlnm)) lapply(cb, function(x) 
                      x[isoyv == i, , ]),
                    dv = dv, 
                    ev = ev)
    }))
    ### (3) Stage 1B3 optimize over the beta3(x) ----
    beta3v <- sapply(seq_len(n), function(j) {
      dv  <- abind::adrop(dtx[,j,,drop = FALSE],2)	
      ev  <- abind::adrop(etx[,j,,drop = FALSE],2)	
      llmax.beta3p(b1  = beta1v[j,], 
                   b2  = beta2v[j],
                   b3  = beta3v[j],  
                   bl  = if (!is.null(prev.dlnm)) sapply(betal, `[[`, j),
                   k2  = kappa2v,   
                   k3  = kappa3v, 
                   gl  = if (!is.null(prev.dlnm)) gammal,
                   cb  = if (!is.null(prev.dlnm)) cb,
                   dv  = dv,         
                   ev  = ev)
    })
    
    
    ### (4) Stage 1K3 optimise over the kappa3(t,r) ----
    kappa3v <- do.call('rbind', lapply(seq_len(mW), function(i) {
      dv  <- dtx[isowv == i,,,drop=FALSE]	
      ev  <- etx[isowv == i,,,drop=FALSE]
      llmax.kappa3p(b1  = beta1v,      
                    b2  = beta2v, 
                    b3  = beta3v,  
                    bl  = if (!is.null(prev.dlnm)) betal,
                    k2  = kappa2v, 
                    k3  = kappa3v[i,], 
                    gl  = if (!is.null(prev.dlnm)) gammal,
                    cb  = if (!is.null(prev.dlnm)) lapply(cb, function(x) 
                      x[isowv == i, , ]),
                    dv  = dv, 
                    ev  = ev)
    }))
    
    
    ### (5) Stage 1A  optimise over the beta1(x) ----
    if(all(int == 0)){
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n) {		 		 
        s1     <- colSums(abind::adrop(dtx[,j,,drop=FALSE],2))
        term.j <- if (length(betal) > 0) 
          Reduce(`+`, Map(`*`, sapply(betal, `[[`, j), dlnm)) else 0
        s2 <- colSums(etx[,j,] * exp(beta2v[j] * kappa2v[isoyv,,drop=FALSE] + 
                                       beta3v[j] * kappa3v[isowv,,drop=FALSE] + term.j))
        beta1v[j,] <- log(s1) - log(s2)
      }
    }
    ### (6) Constraints ----
    if(constraints == "ALL"){
      fac21   <- colMeans(kappa2v) # kappa2v[1,]
      fac31   <- colMeans(kappa3v) # kappa3v[1,]
      fac22   <- sum(beta2v)       # beta2v[1] 
      fac32   <- sum(beta3v)       # beta3v[1]
      kappa2v <- fac22 * sweep(kappa2v, 2, fac21, '-')   
      kappa3v <- fac32 * sweep(kappa3v, 2, fac31, '-')
      beta2v  <- beta2v/fac22     
      beta3v  <- beta3v/fac32
      beta1v  <- beta1v + 
        (beta2v * fac22) %o% fac21 + 
        (beta3v * fac32) %o% fac31
    }
    
    ### (7) Log-likelihood ----
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
                               term.j, 2, beta1v[j,,drop=FALSE], '+'))
    }
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - 
                      etx*mhat - 
                      lgamma(dtx+1))) 
    pen     <- if(!is.null(prev.dlnm)) 
      dlnm_penalty(lambda, betal, gammal, adjmat) else 0
    l1.up   <- l1
    l1.pen  <- l1 - pen
    
    l1      <- l1.pen
    
    if(verbose) cat(l1, ' -> ')
  }
  
  # Calculate number of parameters and deduct the number of constraints
  if(constraints=="ALL"){
    npar <- length(beta1v) + length(beta2v) + length(kappa2v) + 
      length(beta3v) + length(kappa3v) - R*2 - 2
  }
  
  # Output prev dlnm 
  list(beta1  = beta1v,  beta2   = beta2v,  beta3 = beta3v,
       kappa2 = kappa2v, kappa3  = kappa3v, 
       betal  = betal,   gammal  = gammal,
       lambda = lambda,  epsilon = epsilon, 
       mhat   = mhat,    ll      = l1, 
       ll.up  = l1.up,   cb      = cb)	
}

###### (2) Newton Raphson iterations baseline model -----
llmax.beta2p  <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev){
  b21    <- b2
  b20    <- b21 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    Reduce(`+`, Map(`*`, bl, dlnm)) else 0
  thetat <- k2[isoyv,] * ev * 
    exp(sweep(b3 * k3[isowv,,drop = FALSE] + term, 2, b1, '+'))
  s1     <- sum(dv * k2[isoyv,])
  
  while(abs(b21 - b20) > 0.0001){
    b20 <- b21
    f0  <- sum(exp(b20*k2[isoyv,]) * thetat) - s1
    df0 <- sum(exp(b20*k2[isoyv,]) * k2[isoyv,] * thetat)
    b21 <- b20 - f0/df0
  }
  b21
}
llmax.kappa2p <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev){
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
  thetat <- sweep(ev, 2, b2, '*') * 
    exp(sweep(aperm(b3 %o% k3, c(2, 1, 3)) + term, c(2,3), b1, '+'))
  s1     <- apply(sweep(dv, 2, b2, '*'), 3, sum)
  
  while(max(abs(k21 - k20)) > 0.0001) {
    k20 <- k21
    f0  <- apply(sweep(thetat, c(2,3), exp(b2 %o% k20), '*'), 3, sum) - s1
    df0 <- apply(sweep(thetat, c(2,3), sweep(exp(b2 %o% k20), 1, b2, '*'), '*'), 3, sum)
    k21 <- k20 - f0/df0
  }
  k21
} 
llmax.beta3p  <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev){
  b31    <- b3
  b30    <- b31 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    Reduce(`+`, Map(`*`, bl, dlnm)) else 0
  thetat <- k3[isowv,] * ev * 
    exp(sweep(b2 * k2[isoyv,, drop = FALSE] + term, 2, b1, '+'))
  s1     <- sum(dv * k3[isowv,])
  
  while(abs(b31 - b30) > 0.0001) {
    b30 <- b31
    f0  <- sum(exp(b30 * k3[isowv,]) * thetat) - s1
    df0 <- sum(exp(b30 * k3[isowv,]) * k3[isowv,] * thetat)
    b31 <- b30 - f0/df0
  }
  b31
}
llmax.kappa3p <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev){
  sy     <- 1990
  k31    <- k3
  k30    <- k31 - 1
  k2     <- k2[isoyear(as.Date(rownames(dv))) - sy + 1,, drop = FALSE]
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl) > 0) 
    aperm(Reduce(`+`, Map(`%o%`, bl, dlnm)), c(2, 1, 3)) else 0
  thetat <- sweep(ev, 2, b3, '*') * 
    exp(sweep(aperm(b2 %o% k2, c(2, 1, 3)) + term, c(2,3), b1, '+'))
  s1     <- apply(sweep(dv, 2, b3, '*'), 3, sum)
  
  while(max(abs(k31 - k30)) > 0.0001) {
    k30 <- k31
    f0  <- apply(sweep(thetat, c(2,3), exp(b3 %o% k30), '*'), 3, sum) - s1
    df0 <- apply(sweep(thetat, c(2,3), sweep(exp(b3 %o% k30), 1, b3, '*'), '*'), 3, sum)
    k31 <- k30 - f0/df0
  }
  k31
}
