###### (1) Lee-Carter by SVD ----
LC.SVD <- function(lm.mat){
  a    <- rowMeans(lm.mat)
  mbar <- lm.mat - a
  
  sv   <- svd(mbar)
  u    <- sv$u[, 1]
  v    <- sv$v[, 1]
  d1   <- sv$d[1]
  
  b    <- u / sum(u)
  k    <- d1 * v * sum(u)
  
  list(a = a, b = b, k = k)
}

###### (2) Li-Lee by the product-ratio method ----
LL.SVD <- function(lm.arr){
  N  <- dim(lm.arr)[1]
  Tn <- dim(lm.arr)[2]
  R  <- dim(lm.arr)[3]
  
  # Product term: geometric mean across regions (mean on the log scale)
  lp  <- apply(lm.arr, c(1, 2), mean)
  fit <- LC.SVD(lp)
  Ap  <- fit$a
  B   <- fit$b
  K   <- fit$k
  
  # Ratio term: region relative to the product term
  A <- matrix(0, N, R)
  b <- matrix(0, N, R)
  k <- matrix(0, Tn, R)
  
  for(r in seq_len(R)){
    fr     <- LC.SVD(lm.arr[,,r] - lp)
    A[, r] <- Ap + fr$a
    b[, r] <- fr$b
    k[, r] <- fr$k
  }
  
  # Constraints of Min et al.: sum B^2 = 1, sum K = 0, sum b^2 = 1, sum kappa = 0
  s  <- sqrt(sum(B^2)); B <- B/s; K <- K*s
  mK <- mean(K); K <- K - mK; A <- A + outer(B, rep(mK, R))
  for(r in seq_len(R)){
    sr     <- sqrt(sum(b[, r]^2)); b[, r] <- b[, r]/sr; k[, r] <- k[, r]*sr
    mk     <- mean(k[, r]); k[, r] <- k[, r] - mk; A[, r] <- A[, r] + b[, r] * mk
  }
  list(A = A, B = B, K = K, b = b, kappa = k)
}

###### (3) Cross-basis and wave variables ----
cb_min <- function(tg.daily, ref.idx, lag = 21, df.var = 4, df.lag = 4,
                   ref = tg.daily){
  bs.ref <- splines::ns(ref, df = df.var)
  cb <- dlnm::crossbasis(
    tg.daily, lag = lag,
    argvar = list(fun = 'ns',
                  knots          = attr(bs.ref, 'knots'),
                  Boundary.knots = attr(bs.ref, 'Boundary.knots')),
    arglag = list(fun = 'ns', knots = dlnm::logknots(lag, df.lag - 1)))
  matrix(cb[ref.idx, ], nrow = length(ref.idx), dimnames = list(NULL, colnames(cb)))
}

# Heat/cold wave days per ISO week (percentile-based definition, Min et al. Suppl. E)
wave_days <- function(tg.daily, week.idx, probs = c(0.025, 0.975), ref = tg.daily){
  qs  <- quantile(ref, probs = probs, na.rm = TRUE)
  hot <- tg.daily > qs[2]; cld <- tg.daily < qs[1]
  spell <- function(z) { z & c(FALSE, head(z, -1)) & c(FALSE, FALSE, head(z, -2)) }
  data.frame(HWD = as.vector(tapply(spell(hot), week.idx, sum)),
             CWD = as.vector(tapply(spell(cld), week.idx, sum)))
}

###### (4) Backfitting DLNM-LL ----
fit_wDLNM.LL.min <- function(lm.arr, cb.lst, wave.lst, J = 20, delta = 10^(-2),
                             verbose = FALSE){
  N  <- dim(lm.arr)[1]
  Tn <- dim(lm.arr)[2]
  R  <- dim(lm.arr)[3]
  
  lmj  <- lm.arr                                     
  coef <- array(0, dim = c(ncol(cb.lst[[1]]) + 3, N, R))   
  vcv  <- vector('list', N * R)
  j    <- 0
  
  repeat {
    j  <- j + 1
    
    ### (a) Li-Lee step on the current log mortality array
    ll <- LL.SVD(lmj)
    
    ### (b) DLNM step on the partial residuals, per age group and region
    ehat <- array(0, dim = c(N, Tn, R))
    for(r in seq_len(R)){
      X <- cbind(cb.lst[[r]], HWD = wave.lst[[r]]$HWD, CWD = wave.lst[[r]]$CWD)
      for(x in seq_len(N)){
        e   <- lmj[x, , r] - ll$A[x, r]
        fit <- lm(e ~ X)
        cf  <- coef(fit); cf[is.na(cf)] <- 0
        ehat[x, , r] <- cbind(1, X) %*% cf
        coef[, x, r] <- coef[, x, r] + cf
        vcv[[(r - 1) * N + x]] <- vcov(fit)
      }
    }
    
    ### (c) Remove the fitted DLNM component
    lmj <- lmj - aperm(ehat, c(1, 2, 3))
    
    ### (d) Convergence
    if(j > 1){
      dev <- max(abs(ll$A - prev$A), abs(ll$B - prev$B), abs(ll$K - prev$K),
                 abs(ll$b - prev$b), abs(ll$kappa - prev$kappa))
      if(verbose) cat('iteration', j, ': max parameter change =', signif(dev, 3), '\n')
      if(dev < delta) break
    }
    prev <- ll
    if(j >= J) break
  }
  
  fitted <- array(0, dim = c(N, Tn, R))
  for(r in seq_len(R)){
    X  <- cbind(cb.lst[[r]], HWD = wave.lst[[r]]$HWD, CWD = wave.lst[[r]]$CWD)
    for(x in seq_len(N)){
      e   <- lm.arr[x, , r] - (ll$A[x, r] + ll$B[x] * ll$K + ll$b[x, r] * ll$kappa[, r])
      fit <- lm(e ~ X)
      cf  <- coef(fit); cf[is.na(cf)] <- 0
      coef[, x, r] <- cf
      vcv[[(r - 1) * N + x]] <- vcov(fit)
      fitted[x, , r] <- ll$A[x, r] + ll$B[x] * ll$K + ll$b[x, r] * ll$kappa[, r] +
        cbind(1, X) %*% cf
    }
  }
  sigma <- apply(lm.arr - fitted, c(1, 3), sd)
  
  list(A = ll$A, B = ll$B, K = ll$K, b = ll$b, kappa = ll$kappa,
       coef = coef, vcov = vcv, fitted = fitted, sigma = sigma, itr = j)
}