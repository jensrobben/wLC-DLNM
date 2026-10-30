###### (1) Main function ----
wLC.DLNM.NB <-  function(xv, isoyv, isowv, etx, dtx, cb, int, adjmat,
                         constraints, lambda, eps, update.dispersion = TRUE,
                         initialize = NULL, prev.base = NULL, verbose = FALSE,
                         max.iter = 50){
  
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
      # Dispersion values
      ds.x   <- prev.base$phix
      ds.r   <- prev.base$phir
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  
      exp.ds <- aperm(exp.ds, c(3, 1, 2))    
      
      # Parameters
      betal  <- lapply(Q, function(q) (1:n)*0 + 1)
      gammal <- lapply(Q, function(q) matrix(10^(-5), nrow = q, ncol = R))} else {
        # Dispersion values
        ds.x   <- initialize$phix
        ds.r   <- initialize$phir
        exp.ds <- exp(outer(ds.x, ds.r, "+"))
        exp.ds <- replicate(m, exp.ds, simplify = "array")  
        exp.ds <- aperm(exp.ds, c(3, 1, 2))    
        
        # Parameters
        betal  <- initialize$betal
        gammal <- initialize$gammal
      }
  } else {
    # Dispersion values
    ds.x   <- prev.base$phix
    ds.r   <- prev.base$phir
    exp.ds <- exp(outer(ds.x, ds.r, "+"))
    exp.ds <- replicate(m, exp.ds, simplify = "array")  
    exp.ds <- aperm(exp.ds, c(3, 1, 2))    
    
    # Parameters
    betal  <- prev.base$betal
    gammal <- prev.base$gammal
    # lambda <- prev.base$lambda
  }
  
  # Previous Lee-Carter estimates
  beta1v  <- prev.base$beta1
  beta2v  <- prev.base$beta2
  beta3v  <- prev.base$beta3
  kappa2v <- prev.base$kappa2
  kappa3v <- prev.base$kappa3
  
  # Indicators
  if (is.null(names(lambda))) 
    names(lambda) <- names(cb)
  sel1 <- which('Temp' == names(Q))
  sel2 <- which('Infl' == names(Q))
  sel3 <- which('Int'  == names(Q))
  
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
        phi <- exp.ds[,j,]
        llmax.beta.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                        b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                        k2  = kappa2v,    k3  = kappa3v,
                        gl  = gammal,     cb  = cb,
                        dv  = dv,         ev  = ev,
                        phi = phi,        sel = sel1,
                        lambda = lambda,  adjmat = adjmat)
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
      gammal[['Temp']] <- llmax.gamma.dlnm(b1  = beta1v,   b2     = beta2v,
                                           b3  = beta3v,   bl     = betal,
                                           k2  = kappa2v,  k3     = kappa3v,
                                           gl  = gammal,   cb     = cb,
                                           dv  = dtx,      ev     = etx,
                                           phi = exp.ds,   lambda = lambda,
                                           sel = sel1,     adjmat = adjmat)
      
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
        phi <- exp.ds[,j,]
        llmax.beta.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                        b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                        k2  = kappa2v,    k3  = kappa3v,
                        gl  = gammal,     cb  = cb,
                        dv  = dv,         ev  = ev,
                        phi = phi,        sel = sel2,
                        lambda = lambda,  adjmat = adjmat)
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
      gammal[['Infl']] <- llmax.gamma.dlnm(b1  = beta1v,   b2     = beta2v,
                                           b3  = beta3v,   bl     = betal,
                                           k2  = kappa2v,  k3     = kappa3v,
                                           gl  = gammal,   cb     = cb,
                                           dv  = dtx,      ev     = etx,
                                           phi = exp.ds,   lambda = lambda,
                                           sel = sel2,     adjmat = adjmat)
      
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
        phi <- exp.ds[,j,]
        llmax.beta.dlnm(b1  = beta1v[j,], b2  = beta2v[j],
                        b3  = beta3v[j],  bl  = sapply(betal, `[[`, j),
                        k2  = kappa2v,    k3  = kappa3v,
                        gl  = gammal,     cb  = cb,
                        dv  = dv,         ev  = ev,
                        phi = phi,        sel = sel3,
                        lambda = lambda,  adjmat = adjmat)
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
      gammal[['Int']] <- llmax.gamma.dlnm(b1  = beta1v,   b2     = beta2v,
                                          b3  = beta3v,   bl     = betal,
                                          k2  = kappa2v,  k3     = kappa3v,
                                          gl  = gammal,   cb     = cb,
                                          dv  = dtx,      ev     = etx,
                                          phi = exp.ds,   lambda = lambda,
                                          sel = sel3,     adjmat = adjmat)
      
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
    
    ### (8)  Optimize over Phix ----
    if(update.dispersion){
      ds.x <- llmax.phix(mu   = mhat[,,], dv   = dtx[,,],
                         ev   = etx[,,],  phix = ds.x[], 
                         phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))
      
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
    
    ### (7)  Optimize over Phir ----
    if(update.dispersion){
      ds.r <- llmax.phir(mu = mhat, dv = dtx,
                         ev = etx, phix = ds.x, phir = ds.r)
      
      exp.ds <- exp(outer(ds.x, ds.r, "+"))
      exp.ds <- replicate(m, exp.ds, simplify = "array")  # gives an n x R x T array
      exp.ds <- aperm(exp.ds, c(3, 1, 2))  
      
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
    
    ### (9)  Constraints ----
    
    # Now apply the constraints
    if(constraints == "ALL"){
      # fac1   <- sapply(betal, function(x) unname(x[1]))
      fac1   <- sapply(betal, function(x) unname(sum(x)))  
      betal  <- Map(function(b, f) unname(b / f), betal, fac1)
      gammal <- Map(function(g, f) unname(g * f), gammal, fac1)
      if(update.dispersion){
        fac71   <- mean(ds.x)
        ds.r    <- ds.r + fac71
        ds.x    <- ds.x - fac71
      }
      
    }
    
    ### (10)  Log-likelihood ----
    
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
    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat + (etx*mhat)^2/exp.ds)
    l1      <- sum(lgamma(exp.ds + dtx) - lgamma(exp.ds) -
                     lgamma(dtx + 1) + dtx * log(etx * mhat) -
                     dtx * log(exp.ds + etx * mhat) +
                     exp.ds * log(exp.ds) -
                     exp.ds * log(exp.ds + etx * mhat))
    
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
  list(beta1  = beta1v,  beta2   = beta2v,  beta3 = beta3v,
       kappa2 = kappa2v, kappa3  = kappa3v, 
       betal  = betal,   gammal  = gammal,
       phix   = ds.x,    phir    = ds.r, 
       lambda = lambda,  epsilon = epsilon, 
       mhat   = mhat,    ll      = l1,
       cb     = cb,      adjmat  = adjmat,
       ll.up  = l1.up)	
}

### DLNM penalty
dlnm_penalty <- function(lambda, betal, gammal, adjmat){
  if (length(gammal) == 0) return(0)
  sum(sapply(names(gammal), function(k)
    lambda[[k]] / 2 * sum(betal[[k]]^2) *
      sum(diag(gammal[[k]] %*% adjmat %*% t(gammal[[k]])))))
}

###### (2) Compute Hessian wLC-DLNM ----
hessian.wLC.DLNM <- function(etx, dtx, fit){
  # Dimensions
  Q <- unlist(lapply(fit$cb, function(x) dim(x)[2]))
  
  # Penalized Hessian
  H.dlnm <- Reduce('cbind', lapply(1:length(fit$cb), function(j) {
    H.betav  <- hessian.beta(dv = dtx, ev = etx, lca = fit, sel = j) 
    H.gammav <- hessian.gamma(dv = dtx, ev = etx, lca = fit, sel = j)
    cbind(H.betav, H.gammav)
  }))
  
  # set dimnames
  colnames(H.dlnm) = rownames(H.dlnm) <- 
    unlist(lapply(1:length(fit$cb), function(j) {
      selj <- which(names(Q[j]) == c('Temp', 'Infl', 'Int'))
      c(paste0('Beta',3 + selj, '_Age', dimnames(dtx)[[2]]),
        outer(1:Q[j], dimnames(dtx)[[3]], 
              FUN = function(x,y) paste0('Gamma', selj, '_cb', x, "_Region", y)))
    }))
  
  # ind.rm    <- cumsum(c(1, head(Q*R + ncol(dtx), -1)))
  # 
  # H.dlnm[-ind.rm, -ind.rm]
  
  H.dlnm
}
hessian.wLC.DLNM.up <- function(etx, dtx, fit){
  # Dimensions
  Q <- unlist(lapply(fit$cb, function(x) dim(x)[2]))
  
  # Unpenalized hessian
  fit.up <- fit
  fit.up$lambda <- fit$lambda * 0
  H.dlnm.up <- Reduce('cbind', lapply(1:length(fit$cb), function(j) {
    H.betav  <- hessian.beta(dv = dtx, ev = etx, lca = fit.up, sel = j) 
    H.gammav <- hessian.gamma(dv = dtx, ev = etx, lca = fit.up, sel = j)
    cbind(H.betav, H.gammav)
  }))
  
  # set dimnames
  colnames(H.dlnm.up) = rownames(H.dlnm.up) <- 
    unlist(lapply(1:length(fit$cb), function(j) {
      selj <- which(names(Q[j]) == c('Temp', 'Infl', 'Int'))
      c(paste0('Beta',3 + selj, '_Age', dimnames(dtx)[[2]]),
        outer(1:Q[j], dimnames(dtx)[[3]], 
              FUN = function(x,y) paste0('Gamma', selj, '_cb', x, "_Region", y)))
    }))
  
  # ind.rm    <- cumsum(c(1, head(Q*R + ncol(dtx), -1)))
  # 
  # H.dlnm.up[-ind.rm, -ind.rm]
  
  H.dlnm.up
}

###### (3) Newton Raphson iterations DLNM -----
llmax.beta.dlnm  <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, 
                             phi, sel, lambda, adjmat){
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
  thetat <- ev * exp(sweep(b2 * k2[isoyv,] + b3 * k3[isowv,] + term, 2, b1, '+'))
  s1     <- sum(dv * dlnm[[sel]])
  
  while(abs(b41 - b40) > 0.001) {
    b40        <- b41
    thetat.exp <- thetat * exp(b40 * dlnm[[sel]])
    f0         <- sum((dv + phi) * (1 - phi/(phi + thetat.exp)) * 
                        dlnm[[sel]]) + lam.k * b40 * P.k - s1
    df0        <- sum((dv + phi) * (phi * thetat.exp)/((phi + thetat.exp)^2) * 
                        dlnm[[sel]]^2) + lam.k * P.k
    b41 <- b40 - f0/df0
  }
  b41
}

llmax.gamma.dlnm <- function(b1, b2, b3, bl, k2, k3, gl, cb, dv, ev, phi, lambda, adjmat, sel){
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
  thetat <- exp(sweep(aperm(b2 %o% k2[isoyv,], c(2, 1, 3)), c(2,3), b1, '+') +
                  aperm(b3 %o% k3[isowv,], c(2,1,3)) +
                  term) * ev
  s1     <-  sapply(1:R, function(r) 
    colSums(sweep(cb[[sel]][,,r], 1, (dv[,,r]) %*% bl[[sel]], '*')))
  
  while(max(abs(g11 - g10)) > 0.001) {
    g10    <- g11
    dlnm1  <- sapply(1:dim(cb[[sel]])[3], function(r) cb[[sel]][,,r] %*% g10[,r])
    thetat.exp <- thetat * exp(aperm(bl[[sel]] %o% dlnm1, c(2, 1, 3))) 
    
    f0 <- apply(sweep((dv + phi) * (1 - phi/(phi + thetat.exp)), 
                      2, bl[[sel]], '*'), c(1,3), sum)
    f0 <- t(apply(simplify2array(lapply(1:Q[sel], function(q) cb[[sel]][,q,] * f0)), 
                  c(2,3), sum)) - s1
    lam.eff <- lambda[[names(cb)[sel]]] * sum(bl[[sel]]^2)
    f0.pen  <- f0 + lam.eff * 
      matrix(kronecker(adjmat, diag(Q[sel])) %*% as.vector(g10), ncol = R)
    
    df0 <- apply(sweep((dv + phi) * (phi * thetat.exp / ((phi + thetat.exp)^2)), 
                       2, bl[[sel]]^2, '*'), c(1,3), sum)
    df0  <- aperm(simplify2array(lapply(1:Q[sel], function(q) cb[[sel]][,q,] * df0)), c(1,3,2))
    df0    <- lapply(1:R, function(r) {
      apply(abind::abind(lapply(1:Q[sel], function(q) sweep(df0[,,r], 1, cb[[sel]][,q,r], '*')),
                         along = 3), c(2,3), sum, na.rm = TRUE)
    })
    df0 <- abind::abind(df0, along = 3)
    
    HH <- as.matrix(Matrix::bdiag(lapply(1:R, function(r) df0[,,r]))) +
      lam.eff * kronecker(adjmat, diag(Q[sel]))
    DD   <- as.vector(f0.pen)
    step <- tryCatch(solve(HH, DD), error = function(e) MASS::ginv(HH) %*% DD)
    g11  <- matrix(as.vector(g10) - step, nrow = Q[sel], ncol = R)
  }
  g11
}

###### (4) Hessian computations  ----
hessian.beta  <- function(dv, ev, lca, sel){
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
  cb      <- lca$cb
  L       <- lca$adjmat
  lambda  <- lca$lambda
  Q       <- unlist(lapply(cb, function(x) dim(x)[2]))
  
  # DLNM structure
  cb   <- lapply(cb, function(x) {
    x[is.na(x)] <- 0
    x})
  dlnm <- Map(function(cb, g) {
    sapply(seq_len(dim(cb)[3]), function(r) cb[,,r] %*% g[,r])
  }, cb, lca$gammal)
  
  # Calculation
  H.beta <- sapply(1:ncol(dv), function(x){
    deriv <- (dv[,x,] - (dv[,x,] + exp.phi[,x,]) * 
                (1 - (exp.phi[,x,]/(ev[,x,] * mhat[,x,] + exp.phi[,x,]))))
    base <- - (dv[,x,] + exp.phi[,x,]) * exp.phi[,x,] * ev[,x,] * mhat[,x,] *
      (ev[,x,] * mhat[,x,] + exp.phi[,x,])^(-2)
    base <- sweep(base, c(1,2), dlnm[[sel]], '*')
    
    # Vector to store
    vec <- rep(NA, sum(Q*R) + na * length(Q))
    
    # Indices
    lens    <- as.vector(rbind(na, Q*R))           
    starts  <- cumsum(c(1, lens))        
    ind.lst <- Map(seq, starts[-length(starts)], starts[-1]-1)
    ind     <- 1:(2*length(Q))
    
    # First set (selection)
    vec[ind.lst[[2*(sel - 1) + 1]]] <- 
      c(rep(0, (x - 1)), sum(sweep(base, c(1,2), dlnm[[sel]], '*')),
        rep(0, (na - x)))
    vec[ind.lst[[2*sel]]] <- t(sapply(1:Q[sel], function(q)
      colSums(base * cb[[sel]][,q,] * lca$betal[[sel]][x]) + 
        colSums(deriv * cb[[sel]][,q,]))) 
    lam.k <- lambda[[names(Q)[sel]]]; G.k <- lca$gammal[[sel]]
    vec[ind.lst[[2*(sel-1)+1]]][x] <- vec[ind.lst[[2*(sel-1)+1]]][x] - 
      lam.k * sum(diag(G.k %*% L %*% t(G.k)))
    vec[ind.lst[[2*sel]]] <- vec[ind.lst[[2*sel]]] - 
      2 * lam.k * lca$betal[[sel]][x] * as.vector(G.k %*% L)
    ind <- ind[-c(2*(sel-1)+1,2*sel)]
    
    if(length(ind) > 0){
      ijk <- ceiling(ind[1]/2)
      vec[ind.lst[[ind[1]]]] <- 
        c(rep(0, (x - 1)), sum(sweep(base, c(1,2), dlnm[[ijk]], '*')),
          rep(0, (na - x)))
      vec[ind.lst[[ind[2]]]] <- t(sapply(1:Q[ijk], function(q)
        colSums(base * cb[[ijk]][,q,] * lca$betal[[ijk]][x])))
      ind <- ind[-c(1,2)]
    }
    
    if(length(ind) > 0){
      ijk <- ceiling(ind[1]/2)
      vec[ind.lst[[ind[1]]]] <- 
        c(rep(0, (x - 1)), sum(sweep(base, c(1,2), dlnm[[ijk]], '*')),
          rep(0, (na - x)))
      vec[ind.lst[[ind[2]]]] <- t(sapply(1:Q[ijk], function(q)
        colSums(base * cb[[ijk]][,q,] * lca$betal[[ijk]][x])))
      ind <- ind[-c(1,2)]
    }
    
    vec
  })
  
  H.beta
}
hessian.gamma <- function(dv, ev, lca, sel){
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
  cb      <- lca$cb
  L       <- lca$adjmat
  lambda  <- lca$lambda
  Q       <- unlist(lapply(cb, function(x) dim(x)[2]))
  
  # DLNM structure
  cb   <- lapply(cb, function(x) {
    x[is.na(x)] <- 0
    x})
  dlnm <- Map(function(cb, g) {
    sapply(seq_len(dim(cb)[3]), function(r) cb[,,r] %*% g[,r])
  }, cb, lca$gammal)
  
  # Define grid
  grid <- expand.grid(qq = 1:Q[sel], Region = 1:R)
  
  # Calculation
  H.gamma <- sapply(1:nrow(grid), function(j){
    # Grid values
    q <- grid$qq[j]
    r <- grid$Region[j]
    
    deriv <- (dv[,,r] - (dv[,,r] + exp.phi[,,r]) * 
                (1 - (exp.phi[,,r]/(ev[,,r] * mhat[,,r] + exp.phi[,,r]))))
    base <- - (dv[,,r] + exp.phi[,,r]) * exp.phi[,,r] * ev[,,r] * mhat[,,r] *
      (ev[,,r] * mhat[,,r] + exp.phi[,,r])^(-2)
    base <- sweep(base, c(1,2), cb[[sel]][,q,r] %o% lca$betal[[sel]], '*')
    
    # Vector to store
    vec <- rep(NA, sum(Q*R) + na * length(Q))
    
    # Indices
    lens    <- as.vector(rbind(na, Q*R))           
    starts  <- cumsum(c(1, lens))        
    ind.lst <- Map(seq, starts[-length(starts)], starts[-1]-1)
    ind     <- 1:(2*length(Q))
    
    lam.k <- lambda[[names(Q)[sel]]]
    S.k   <- sum(lca$betal[[sel]]^2)
    GL.k  <- lca$gammal[[sel]] %*% L
    
    vec[ind.lst[[2*(sel - 1) + 1]]] <- 
      colSums(sweep(base, 1, dlnm[[sel]][,r], '*')) + 
      colSums(sweep(deriv, 1, cb[[sel]][,q,r], '*')) -
      2 * lam.k * lca$betal[[sel]] * GL.k[q, r]
    vec[ind.lst[[2*sel]]] <- sapply(1:R, function(s) {
      if(r == s)
        sapply(1:Q[sel], function(qp) 
          sum(sweep(base, c(1,2), cb[[sel]][,qp,s] %o% lca$betal[[sel]], '*')) - 
            lam.k * S.k * L[r,s] * (q == qp)) else
              sapply(1:Q[sel], function(qp) - lam.k * S.k * L[r,s] * (q == qp))
    })
    ind <- ind[-c(2*(sel-1)+1,2*sel)]
    
    if(length(ind) > 0){
      ijk <- ceiling(ind[1]/2)
      vec[ind.lst[[ind[1]]]] <- colSums(sweep(base, 1, dlnm[[ijk]][,r], '*'))
      vec[ind.lst[[ind[2]]]] <- sapply(1:R, function(s) {
        if(r == s)
          sapply(1:Q[ijk], function(qp)
            sum(sweep(base, c(1,2), cb[[ijk]][,qp,s] %o% lca$betal[[ijk]], '*'))) else
              rep(0, Q[ijk])
      })
      ind <- ind[-c(1,2)]
    }
    
    if(length(ind) > 0){
      ijk <- ceiling(ind[1]/2)
      vec[ind.lst[[ind[1]]]] <- colSums(sweep(base, 1, dlnm[[ijk]][,r], '*'))
      vec[ind.lst[[ind[2]]]] <- sapply(1:R, function(s) {
        if(r == s)
          sapply(1:Q[ijk], function(qp)
            sum(sweep(base, c(1,2), cb[[ijk]][,qp,s] %o% lca$betal[[ijk]], '*'))) else
              rep(0, Q[ijk])
      })
      ind <- ind[-c(1,2)]
    }
    
    vec
  })
  
  H.gamma
}


###### (5) Adopt centering in cross-basis
center_cb <- function(cb, isowv, means = NULL){
  if (is.null(means)) {
    means <- simplify2array(lapply(seq_len(dim(cb)[3]), function(r)
      apply(cb[,,r, drop = FALSE][,,1], 2, function(z) 
        tapply(z, isowv, mean, na.rm = TRUE))))
    dimnames(means) <- list(sort(unique(isowv)), NULL, dimnames(cb)[[3]])
  }
  out <- cb
  for (r in seq_len(dim(cb)[3]))
    out[,,r] <- cb[,,r] - means[as.character(isowv),,r, drop = FALSE][,,1]
  attr(out, "means") <- means
  out
}
###### (5) Additional functions ----
sumzero <- function(k) rbind(diag(k - 1), -1)

contrast_dlnm <- function(fit){
  nA <- length(fit$betal[[1]])
  R <- ncol(fit$gammal[[1]])
  bl <- list()
  for (k in names(fit$gammal))
    bl <- c(bl, list(sumzero(nA), diag(nrow(fit$gammal[[k]]) * R)))
  as.matrix(Matrix::bdiag(bl))
}

cov_constrained <- function(H, C){
  V <- C %*% solve(t(C) %*% (-H) %*% C) %*% t(C)
  dimnames(V) <- dimnames(H); V
}

edf_constrained <- function(H.pen, H.up, C)
  sum(diag((t(C) %*% H.up %*% C) %*% solve(t(C) %*% H.pen %*% C)))