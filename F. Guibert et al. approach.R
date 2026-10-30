###### (1) Cross-basis matrices ----

# Info basis functions
basis_info <- function(tg.daily, lag = 21) {
  list(knots = quantile(tg.daily, c(0.10, 0.75, 0.90), na.rm = TRUE),
       bk    = range(tg.daily, na.rm = TRUE),
       lag   = lag)
}

# Cross-basis Guibert et al. paper
cb_guibert <- function(tg.daily, ref.idx, info){
  cb <- dlnm::crossbasis(tg.daily,
                         lag = info$lag,
                         argvar = list(fun = 'ns', knots = info$knots,
                                       Boundary.knots = info$bk),
                         arglag = list(fun = 'ns', intercept = TRUE,
                                       knots = dlnm::logknots(info$lag, 3)))
  
  matrix(cb[ref.idx, ], nrow = length(ref.idx),
         dimnames = list(NULL, colnames(cb)))
}

# Cumulative log relative risk
cum_logRR <- function(v, theta, info){
  Bv <- dlnm::onebasis(v, fun = 'ns', 
                       knots = info$knots, 
                       Boundary.knots = info$bk)
  Bl <- dlnm::onebasis(0:info$lag, 
                       fun = 'ns', 
                       intercept = TRUE,
                       knots = dlnm::logknots(info$lag, 3))
  
  as.vector(mgcv::tensor.prod.model.matrix(
    list(Bv, matrix(colSums(Bl), nrow = 1)[
      rep(1, length(v)), , drop = FALSE])) %*% theta)
}

###### (2) Stage 1: DLNM and temperature adjustment ----
stage1_guibert <- function(dtx, cb.lst, tg.daily.lst, week.lst, info.lst, 
                           time.idx, isoyv, strata, df.year = 8, n.year, 
                           verbose = FALSE){
  
  Tn     <- dim(dtx)[1]
  N      <- dim(dtx)[2]
  R      <- dim(dtx)[3]
  ns.t   <- splines::ns(time.idx, df = df.year * n.year)
  strata <- as.factor(strata)
  S      <- levels(strata)
  
  theta  <- vcv <- vector('list', N * R)
  Tadj.w <- array(NA, dim = c(Tn, N, R))
  Tadj.y <- array(NA, dim = c(length(unique(isoyv)), N, R))
  mmt    <- numeric(R)
  
  for(r in seq_len(R)){
    
    # Grid over the daily temperature distribution 
    grd <- seq(quantile(tg.daily.lst[[r]], 0.01), 
               quantile(tg.daily.lst[[r]], 0.99),
               length.out = 200)
    
    # Minimum mortality temperature, from a DLNM on the deaths of all ages together
    fit.all <- glm(rowSums(dtx[, , r]) ~ cb.lst[[r]] + ns.t, 
                   family = quasipoisson())
    cf.all  <- coef(fit.all)
    cf.all[is.na(cf.all)] <- 0
    mmt[r]  <- grd[which.min(cum_logRR(grd, 
                                       cf.all[grep('^cb.lst', names(cf.all))],
                                       info.lst[[r]]))]
    
    for(k in S){
      ag  <- which(strata == k)                        
      fit <- glm(rowSums(dtx[, ag, r, drop = FALSE], dims = 1) ~ 
                   cb.lst[[r]] + ns.t,
                 family = quasipoisson())
      cf  <- coef(fit)
      cf[is.na(cf)] <- 0
      sel <- grep('^cb.lst', names(cf))
      
      # Daily attributable fractions, averaged within each ISO week
      AF.d <- 1 - exp(-(cum_logRR(tg.daily.lst[[r]], cf[sel], info.lst[[r]]) -
                          cum_logRR(mmt[r], cf[sel], info.lst[[r]])))
      AF.w <- as.vector(tapply(AF.d, week.lst[[r]], mean))
      Dhat <- fitted(fit)
      T.y  <- (1 - rowsum(AF.w * Dhat, isoyv)[, 1]/rowsum(Dhat, isoyv)[, 1])^(-1)
      
      # The stratum's estimates apply to every age group inside it
      for(x in ag){
        theta[[(r - 1) * N + x]] <- cf[sel]
        vcv[[(r - 1) * N + x]]   <- vcov(fit)[sel, sel]
        Tadj.w[, x, r]           <- (1 - AF.w)^(-1)
        Tadj.y[, x, r]           <- T.y
      }
    }
    
    if(verbose) cat('Stage 1 - region', regions[r], 
                    ': MMT =', round(mmt[r], 1), '\n')
  }
  
  list(theta = theta, vcov = vcv, Tadj.w = Tadj.w, Tadj.y = Tadj.y,
       mmt = mmt, strata = strata, info = info.lst)
}

###### (3) Temperature adjustment for new temperature values ----
Tadj_guibert <- function(tg.daily.new.lst, week.new.lst, st1, N, R){
  out <- array(NA, dim = c(length(unique(week.new.lst[[1]])), N, R))
  
  for(r in seq_len(R)) for(x in seq_len(N)){
    th   <- st1$theta[[(r - 1) * N + x]]
    AF.d <- 1 - exp(-(cum_logRR(tg.daily.new.lst[[r]], th, st1$info[[r]]) -
                        cum_logRR(st1$mmt[r], th, st1$info[[r]])))
    out[, x, r] <- (1 - as.vector(tapply(AF.d, week.new.lst[[r]], mean)))^(-1)
  }
  
  out
}