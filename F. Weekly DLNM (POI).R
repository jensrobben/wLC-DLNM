###### (1) Main function ----
wLC.DLNM.POI <-  function(xv, isoyv, isowv, etx, dtx, cb, int, adjmat,
                          constraints, lambda, eps, initialize = NULL, 
                          prev.base = NULL, verbose = FALSE, max.iter = 50){
  
  # Notation
  n   <- length(xv)
  mY  <- length(unique(isoyv))
  mW  <- length(unique(isowv)) 
  mWv <- sapply(split(isowv, isoyv), max)
  m   <- nrow(etx)
  R   <- dim(dtx)[3]
  Q   <- unlist(lapply(cb, function(x) dim(x)[2]))
  
  # Death rates
  mtx <- dtx/etx	 
  mtx[mtx==0] <- 10^(-8)
  
  # Initialization
  if(is.null(prev.base$betal)){
    if(is.null(initialize)){
      # Parameters
      betal  <- lapply(Q, function(q) (1:n)*0 + 1)
      gammal <- lapply(Q, function(q) matrix(10^(-5), nrow = q, ncol = R))
      } else {
        # Parameters
        betal  <- initialize$betal
        gammal <- initialize$gammal
      }
  } else {
    # Parameters
    betal  <- prev.base$betal
    gammal <- prev.base$gammal
  }
  
  # Previous Lee-Carter estimates
  beta1v  <- prev.base$beta1
  beta2v  <- prev.base$beta2
  beta3v  <- prev.base$beta3
  kappa2v <- prev.base$kappa2
  kappa3v <- prev.base$kappa3
  
  # Indicators
  sel1 <- which('Temp' == names(Q))
  sel2 <- which('Infl' == names(Q))
  sel3 <- which('Int'  == names(Q))
  if (is.null(names(lambda))) names(lambda) <- names(cb)
  
  # Stage 1: iterate
  l0   <- -100
  l1   <- -99
  iter <- 0
  
  while(abs((l1-l0)/l0) > eps & iter <= max.iter)  {
    iter <- iter + 1
    l0   <- l1
    
    ### (1)  Optimize over beta4 ----
    if(length(sel1) > 0){
      betal[['Temp']] <- sapply(seq_len(n), function(j) {
        dv  <- dtx[,j,]
        ev  <- etx[,j,]
        llmax.betap.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                         b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                         k2  = kappa2v,    k3  = kappa3v,
                         gl  = gammal,     cb  = cb,
                         dv  = dv,         ev  = ev,
                         sel = sel1,       lambda = lambda,  
                         adjmat = adjmat)
      })
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    ### (2)  Optimize over gamma4 ----
    if(length(sel1) > 0){
      gammal[['Temp']] <- llmax.gammap.dlnm(b1  = beta1v,   b2     = beta2v,
                                            b3  = beta3v,   bl     = betal,
                                            k2  = kappa2v,  k3     = kappa3v,
                                            gl  = gammal,   cb     = cb,
                                            dv  = dtx,      ev     = etx,
                                            sel = sel1,     lambda = lambda,
                                            adjmat = adjmat)
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    
    
    
    ### (3)  Optimize over beta5 ----
    if(length(sel2) > 0){
      betal[['Infl']] <- sapply(seq_len(n), function(j) {
        dv  <- dtx[,j,]
        ev  <- etx[,j,]
        llmax.betap.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                         b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                         k2  = kappa2v,    k3  = kappa3v,
                         gl  = gammal,     cb  = cb,
                         dv  = dv,         ev  = ev,
                         sel = sel2,       lambda = lambda,  
                         adjmat = adjmat)
      })
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    
    ### (4)  Optimize over gamma5 ----
    if(length(sel2) > 0){
      gammal[['Infl']] <- llmax.gammap.dlnm(b1  = beta1v,   b2     = beta2v,
                                            b3  = beta3v,   bl     = betal,
                                            k2  = kappa2v,  k3     = kappa3v,
                                            gl  = gammal,   cb     = cb,
                                            dv  = dtx,      ev     = etx,
                                            sel = sel2,     lambda = lambda,
                                            adjmat = adjmat)
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    
    ### (5)  Optimize over beta5 ----
    if(length(sel3) > 0){
      betal[['Int']] <- sapply(seq_len(n), function(j) {
        dv  <- dtx[,j,]
        ev  <- etx[,j,]
        llmax.betap.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                         b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                         k2  = kappa2v,    k3  = kappa3v,
                         gl  = gammal,     cb  = cb,
                         dv  = dv,         ev  = ev,
                         sel = sel3,       lambda = lambda,  
                         adjmat = adjmat)
      })
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    
    ### (6)  Optimize over gamma5 ----
    if(length(sel3) > 0){
      gammal[['Int']] <- llmax.gammap.dlnm(b1  = beta1v,   b2     = beta2v,
                                           b3  = beta3v,   bl     = betal,
                                           k2  = kappa2v,  k3     = kappa3v,
                                           gl  = gammal,   cb     = cb,
                                           dv  = dtx,      ev     = etx,
                                           sel = sel3,   lambda = lambda,
                                           adjmat = adjmat)
      
      mhat   <- mtx*NA
      dlnm <- Map(function(c, g){
        out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
        out[is.na(out)] <- 0
        out
      }, cb, gammal)
      for(j in 1:n){
        mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                                 beta3v[j] * kappa3v[isowv,] +
                                 Reduce(`+`, Map(function(b, d) 
                                   b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
      }
    }
    
    ### (7)  Constraints ----
    
    if(constraints == "ALL"){
      fac1   <- sapply(betal, function(x) unname(sum(x)))
      betal  <- Map(function(b, f) unname(b / f), betal, fac1)
      gammal <- Map(function(g, f) unname(g * f), gammal, fac1)
    }
    ### (8)  Log-likelihood ----
    
    mhat   <- mtx*NA
    dlnm <- Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, cb, gammal)
    for(j in 1:n){
      mhat[,j,] <- exp(sweep(beta2v[j] * kappa2v[isoyv,] +
                               beta3v[j] * kappa3v[isowv,] +
                               Reduce(`+`, Map(function(b, d) 
                                 b[j] * d, betal, dlnm)), 2, beta1v[j,], '+'))
    }
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - 
                      etx*mhat - 
                      lgamma(dtx+1))) 
    pen     <- dlnm_penalty(lambda, betal, gammal, adjmat)
    l1.up   <- l1
    l1.pen  <- l1 - pen
    l1      <- l1.pen
    
    if(verbose) cat(l1, ' -> ')
  }
  
  # Set some names
  names(betal)  <- names(cb)
  names(gammal) <- names(cb)
  
  # Output
  list(beta1   = beta1v,  beta2   = beta2v,  beta3 = beta3v,
       kappa2  = kappa2v, kappa3  = kappa3v, 
       betal   = betal,   gammal  = gammal,
       lambda  = lambda,  epsilon = epsilon, 
       mhat    = mhat,    ll      = l1,
       ll.up   = l1.up,   cb      = cb,
       adjmat  = adjmat)	
}

  

###### (3) Newton Raphson iterations DLNM -----
llmax.betap.dlnm  <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, sel, 
                              lambda,  adjmat){
  lam.k  <- lambda[[names(cb)[sel]]]
  P.k    <- sum(diag(gl[[sel]] %*% adjmat %*% t(gl[[sel]])))
  b41    <- bl[sel]
  b40    <- b41 - 1
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb, gl)
  term   <- if (length(bl[-sel]) > 0) 
    Reduce(`+`, Map(`*`, bl[-sel], dlnm[-sel])) else 0
  thetat <- dlnm[[sel]] * ev * 
    exp(sweep(b2 * k2[isoyv,] + b3 * k3[isowv,] + term, 2, b1, '+'))
  s1     <- sum(dv * dlnm[[sel]])
  
  while(abs(b41 - b40) > 0.001) {
    b40 <- b41
    f0  <- sum(exp(b40 * dlnm[[sel]]) * thetat) - s1 + lam.k * b40 * P.k
    df0 <- sum(exp(b40 * dlnm[[sel]]) * dlnm[[sel]] * thetat) + lam.k * P.k
    b41 <- b40 - f0/df0
  }
  b41
}
llmax.gammap.dlnm <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, lambda, adjmat, sel){
  R      <- ncol(gl[[1]])
  Q      <- unlist(lapply(cb, function(x) dim(x)[2]))
  g11    <- gl[[sel]]
  g10    <- g11 - 1
  cb     <- lapply(cb, function(x) { x[is.na(x)] <- 0; x })
  dlnm <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, cb[-sel], gl[-sel])
  term <- if (length(bl[-sel]) > 0)
    aperm(Reduce(`+`, Map(`%o%`, bl[-sel], dlnm)), c(2, 1, 3)) else 0
  thetat <- sweep(ev, 2, bl[[sel]], '*') * 
    exp(sweep(aperm(b2 %o% k2[isoyv,], c(2, 1, 3)), c(2,3), b1, '+') +
          aperm(b3 %o% k3[isowv,], c(2,1,3)) +
          term)
  thetat <- aperm(simplify2array(lapply(1:Q[[sel]], function(q) 
    sweep(thetat, c(1,3), cb[[sel]][,q,], '*'))), c(1,2,4,3))
  
  s1     <-  sapply(1:R, function(r) 
    colSums(sweep(cb[[sel]][,,r], 1, dv[,,r] %*% bl[[sel]], '*'), na.rm = TRUE))
  
  while(max(abs(g11 - g10)) > 0.001) {
    g10    <- g11
    dlnm1  <- sapply(1:dim(cb[[sel]])[3], function(r) 
      cb[[sel]][,,r] %*% g10[,r])
    f0     <- sapply(1:R, function(r) {
      mu <- aperm(exp(bl[[sel]] %o% dlnm1[,r]), c(2,1))
      v  <- sweep(thetat[,,,r], c(1,2), mu, '*')
      apply(v, 3, sum, na.rm = TRUE) - s1[,r]
    })
    lam.eff <- lambda[[names(cb)[sel]]] * sum(bl[[sel]]^2)
    f0.pen <- f0 + lam.eff * 
      matrix(kronecker(adjmat, diag(Q[sel])) %*% as.vector(g10), ncol = R)
    df0    <- lapply(1:R, function(r) {
      mu <- aperm(exp(bl[[sel]] %o% dlnm1[,r]), c(2,1))
      v1 <- sweep(mu, 2, bl[[sel]], '*')
      v2 <- apply(sweep(thetat[,,,r], c(1,2), v1, '*'), c(1,3), sum)
      apply(abind::abind(lapply(1:Q[sel], function(q) sweep(v2, 1, cb[[sel]][,q,r], '*')),
                         along = 3), c(2,3), sum, na.rm = TRUE)
    })
    df0 <- abind::abind(df0, along = 3)
    HH <- as.matrix(Matrix::bdiag(lapply(1:R, function(r) df0[,,r]))) +
      lam.eff * kronecker(adjmat, diag(Q[sel]))
    DD   <- as.vector(f0.pen)
    g11 <- matrix(as.vector(g10) - solve(HH) %*% DD, 
                  nrow = Q[sel], ncol = R)
  }
  g11
}

