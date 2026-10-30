###### (1) Features: excess heat and cold indices, influenza anomalies ----

excess_index <- function(tg.daily, week.idx, ref = tg.daily, probs = c(0.05, 0.95)){
  qs  <- quantile(ref, probs = probs, na.rm = TRUE)
  
  data.frame(EHI = as.vector(tapply(pmax(0, tg.daily - qs[2]), week.idx, mean)),
             ECI = as.vector(tapply(pmax(0, qs[1] - tg.daily), week.idx, mean)))
}

influenza_index <- function(ia, isowv, ref.idx = seq_along(ia)){
  avg <- tapply(ia[ref.idx], isowv[ref.idx], mean)
  ano <- ia - avg[as.character(isowv)]
  
  pmax(ano - quantile(ano[ref.idx], 0.75), 0)
}

add_lags <- function(z){
  cbind('l0'   = z,
        'l1'   = c(NA, head(z, -1)),
        'l2.3' = (c(NA, NA, head(z, -2)) + c(NA, NA, NA, head(z, -3)))/2)
}

###### (2) Baseline model ----

fit_baseline_rs <- function(df, adjmat, regions, ages, verbose = FALSE){
  df <- as.data.frame(df) %>%
    dplyr::mutate('fsin52' = sin(2*pi*ISOWeek/52.1775),
                  'fcos52' = cos(2*pi*ISOWeek/52.1775),
                  'fsin26' = sin(4*pi*ISOWeek/52.1775),
                  'fcos26' = cos(4*pi*ISOWeek/52.1775),
                  Region   = factor(Region, levels = regions),
                  Age      = factor(Age, levels = ages),
                  Time     = ISOYear - min(ISOYear))
  
  fits <- vector('list', length(ages))
  names(fits) <- ages
  
  for(a in seq_along(ages)){
    dfa <- df[df$Age == ages[a], ]
    
    dfa$Deaths_f <- NA
    for(r in regions){
      sel <- dfa$Region == r
      dfa[sel, 'Deaths_f'] <- round(as.numeric(forecast::tsclean(
        ts(dfa$Deaths[sel], frequency = 52))))
    }
    
    # Model matrices, one block per baseline parameter
    dfa$M0 <- model.matrix(~ -1 + Region, data = dfa)
    dfa$M1 <- model.matrix(~ -1 + Time:Region, data = dfa)
    dfa$M2 <- model.matrix(~ -1 + fsin52:Region, data = dfa)
    dfa$M3 <- model.matrix(~ -1 + fcos52:Region, data = dfa)
    dfa$M4 <- model.matrix(~ -1 + fsin26:Region, data = dfa)
    dfa$M5 <- model.matrix(~ -1 + fcos26:Region, data = dfa)
    
    nR <- length(regions)
    fits[[a]] <- mgcv::gam(Deaths_f ~ -1 + M0 + M1 + M2 + M3 + M4 + M5,
                           offset = log(dfa$Expo), data = dfa,
                           family = poisson(link = 'log'),
                           control = list(epsilon = 10^(-5), trace = verbose,
                                          newton = list(conv.tol = 1e-5)),
                           drop.intercept = TRUE,
                           paraPen = list(M0 = list(adjmat, rank = nR, sp = -1),
                                          M1 = list(adjmat, rank = nR, sp = -1),
                                          M2 = list(adjmat, rank = nR, sp = -1),
                                          M3 = list(adjmat, rank = nR, sp = -1),
                                          M4 = list(adjmat, rank = nR, sp = -1),
                                          M5 = list(adjmat, rank = nR, sp = -1)))
    if(verbose) cat('baseline: age group', ages[a], 'done\n')
  }
  
  list(fits = fits, ages = ages, regions = regions, year.min = min(df$ISOYear))
}

# Baseline death counts for a (new) data frame with the same structure
predict_baseline_rs <- function(obj, newdf, regions, ages){
  newdf <- as.data.frame(newdf) %>%
    dplyr::mutate('fsin52' = sin(2*pi*ISOWeek/52.1775),
                  'fcos52' = cos(2*pi*ISOWeek/52.1775),
                  'fsin26' = sin(4*pi*ISOWeek/52.1775),
                  'fcos26' = cos(4*pi*ISOWeek/52.1775),
                  Region   = factor(Region, levels = regions),
                  Age      = factor(Age, levels = ages),
                  Time     = ISOYear - obj$year.min)
  
  out <- rep(NA_real_, nrow(newdf))
  for(a in seq_along(ages)){
    sel <- which(newdf$Age == ages[a])
    dfa <- newdf[sel, ]
    dfa$M0 <- model.matrix(~ -1 + Region, data = dfa)
    dfa$M1 <- model.matrix(~ -1 + Time:Region, data = dfa)
    dfa$M2 <- model.matrix(~ -1 + fsin52:Region, data = dfa)
    dfa$M3 <- model.matrix(~ -1 + fcos52:Region, data = dfa)
    dfa$M4 <- model.matrix(~ -1 + fsin26:Region, data = dfa)
    dfa$M5 <- model.matrix(~ -1 + fcos26:Region, data = dfa)
    
    out[sel] <- as.numeric(predict(obj$fits[[a]], dfa, type = 'response')) * dfa$Expo
  }
  out
}

###### (3) Objective functions for the M-step ----
# Column-wise products 
colprods <- function(m) exp(colSums(log(m)))

# Expected Poisson log-likelihood in the shock states 
obj.alpha <- function(alpha, msprobs, xmat, dtr, btr, nA){
  Pt  <- do.call('rbind', replicate(nA, msprobs, simplify = FALSE))
  xa  <- sapply(1:dim(xmat)[3], function(j) t(t(xmat[,,j])) %*% alpha)
  
  - sum(Pt * (-btr*exp(xa) + dtr*xa + dtr*log(btr) - lgamma(dtr+1)))
}

# Expected log-likelihood of the transition probabilities
obj.beta <- function(beta, jsprobs, xmat01, xmat02, xmat11, xmat22,
                     pena01, pena02, pena11, pena22, R){
  i01 <- dim(xmat01)[2]
  i02 <- dim(xmat02)[2]
  i11 <- dim(xmat11)[2]
  i22 <- dim(xmat22)[2]
  b01 <- beta[1:i01]
  b02 <- beta[(i01+1):(i01+i02)]
  b11 <- beta[(i01+i02+1):(i01+i02+i11)]
  b22 <- beta[(i01+i02+i11+1):(i01+i02+i11+i22)]
  
  x1b01 <- pena01 + sapply(1:R, function(r) t(t(xmat01[,,r])) %*% b01)
  x2b02 <- pena02 + sapply(1:R, function(r) t(t(xmat02[,,r])) %*% b02)
  x1b11 <- pena11 + sapply(1:R, function(r) t(t(xmat11[,,r])) %*% b11)
  x2b22 <- pena22 + sapply(1:R, function(r) t(t(xmat22[,,r])) %*% b22)
  
  exp.x1b01 <- exp(x1b01)
  exp.x2b02 <- exp(x2b02)
  exp.x1b11 <- exp(x1b11)
  exp.x2b22 <- exp(x2b22)
  
  p00 <- 1/(1 + exp.x1b01 + exp.x2b02)
  p01 <- exp.x1b01/(1 + exp.x1b01 + exp.x2b02)
  p02 <- exp.x2b02/(1 + exp.x1b01 + exp.x2b02)
  p10 <- 1/(1 + exp.x1b11)
  p11 <- exp.x1b11/(1 + exp.x1b11)
  p20 <- 1/(1 + exp.x2b22)
  p22 <- exp.x2b22/(1 + exp.x2b22)
  
  # avoid numerical issues
  for(p in c('p00','p01','p02','p10','p11','p20','p22')){
    z <- get(p)
    z[which(z == 0)] <- 10^(-300)
    assign(p, z)
  }
  
  - sum(jsprobs[-1,1:R]             * log(p00[-1,]) +
          jsprobs[-1,(R+1):(2*R)]   * log(p01[-1,]) +
          jsprobs[-1,(2*R+1):(3*R)] * log(p02[-1,]) +
          jsprobs[-1,(3*R+1):(4*R)] * log(p10[-1,]) +
          jsprobs[-1,(4*R+1):(5*R)] * log(p11[-1,]) +
          jsprobs[-1,(6*R+1):(7*R)] * log(p20[-1,]) +
          jsprobs[-1,(8*R+1):(9*R)] * log(p22[-1,]))
}

obj.gamma <- function(gamma, msprobs, R){
  gvec <- rep(c(gamma, 0, 1 - gamma), each = R)
  
  - sum(msprobs * log(gvec + 10^(-300)))
}

###### (4) EM algorithm ----
fit_regime_switching <- function(df, form1, form2, regions, ages,
                                 eps = 10^(-7), max.iter = 100, verbose = TRUE){
  
  df <- as.data.frame(df)
  R  <- length(regions)
  nA <- length(ages)
  nT <- nrow(df)/(nA*R)
  
  # Design matrices of the state-specific Poisson models
  fit1.0 <- glm(form1, data = df, family = poisson(link = 'log'))
  m10    <- model.matrix(fit1.0)
  x1     <- simplify2array(lapply(1:R, function(j)
    matrix(m10[(nA*nT*(j-1) + 1):(nA*nT*j), ], ncol = ncol(m10))))
  alpha.1 <- unname(fit1.0$coefficients)
  
  fit2.0 <- glm(form2, data = df, family = poisson(link = 'log'))
  m20    <- model.matrix(fit2.0)
  x2     <- simplify2array(lapply(1:R, function(j)
    matrix(m20[(nA*nT*(j-1) + 1):(nA*nT*j), ], ncol = ncol(m20))))
  alpha.2 <- unname(fit2.0$coefficients)
  
  # Design matrices of the transition probabilities (age invariant)
  dfa  <- df %>% dplyr::filter(Age == ages[nA])
  mmat <- function(f) {
    m <- model.matrix(f, data = dfa)
    simplify2array(lapply(1:R, function(j)
      matrix(m[(nT*(j-1) + 1):(nT*j), ], ncol = ncol(m))))
  }
  xx01 <- mmat(~ w_avg_ehi)
  xx11 <- mmat(~ w_avg_ehi)
  xx02 <- mmat(~ w_avg_eci)
  xx22 <- mmat(~ w_avg_eci_l2.3)
  
  beta.01 <- rep(0, dim(xx01)[2])
  beta.02 <- rep(0, dim(xx02)[2])
  beta.11 <- rep(0, dim(xx11)[2])
  beta.22 <- rep(0, dim(xx22)[2])
  beta.0  <- c(beta.01, beta.02)
  gamma   <- 0.3
  
  # Data in matrix form
  dt    <- matrix(df$Deaths,    nrow = nT*nA, ncol = R, 
                  dimnames = list(NULL, regions))
  bt    <- matrix(df$bDeaths,   nrow = nT*nA, ncol = R, 
                  dimnames = list(NULL, regions))
  pen01 <- matrix(df$penalty01, nrow = nT*nA, ncol = R)[1:nT, ]
  pen02 <- matrix(df$penalty02, nrow = nT*nA, ncol = R)[1:nT, ]
  pen11 <- pen01*0
  pen22 <- pen02*0
  
  create_ind_lst <- function(i)
    unlist(rep(list((3*R*i+1):(3*R*i+R), (3*R*i+R+1):(3*R*i+2*R),
                    (3*R*i+2*R+1):(3*R*i+3*R)), each = 3))
  indl0 <- create_ind_lst(0)
  indl1 <- create_ind_lst(1)
  indl2 <- create_ind_lst(2)
  
  ll0  <- 20*10^8
  ll1  <- 10*10^8
  iter <- 0
  best <- NULL
  
  while(abs(ll0 - ll1)/abs(ll0) > eps & iter < max.iter){
    ll0  <- ll1
    iter <- iter + 1
    
    ### E-step
    
    # Conditional densities
    x1a1  <- sapply(1:dim(x1)[3], function(j) t(t(x1[,,j])) %*% alpha.1)
    x2a2  <- sapply(1:dim(x2)[3], function(j) t(t(x2[,,j])) %*% alpha.2)
    lfd   <- cbind(dpois(dt, lambda = bt, log = TRUE),
                   dpois(dt, lambda = bt * exp(x1a1), log = TRUE),
                   dpois(dt, lambda = bt * exp(x2a2), log = TRUE))
    
    # Transition probabilities
    x1b01 <- sapply(1:R, function(j) t(t(xx01[,,j])) %*% beta.01)
    x2b02 <- sapply(1:R, function(j) t(t(xx02[,,j])) %*% beta.02)
    x1b11 <- sapply(1:R, function(j) t(t(xx11[,,j])) %*% beta.11)
    x2b22 <- sapply(1:R, function(j) t(t(xx22[,,j])) %*% beta.22)
    
    pt.01 <- exp(pen01 + x1b01)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
    pt.02 <- exp(pen02 + x2b02)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
    pt.00 <- 1 - pt.01 - pt.02
    
    pt.11 <- exp(pen11 + x1b11)/(1 + exp(pen11 + x1b11))
    pt.12 <- matrix(0, nrow = nT, ncol = R)
    pt.10 <- 1 - pt.11
    
    pt.21 <- matrix(0, nrow = nT, ncol = R)
    pt.22 <- exp(pen22 + x2b22)/(1 + exp(pen22 + x2b22))
    pt.20 <- 1 - pt.22
    
    P <- unname(cbind(pt.00, pt.01, pt.02, pt.10, pt.11, pt.12, 
                      pt.20, pt.21, pt.22))
    P[1, ] <- NA
    
    gvec <- rep(c(gamma, 0, 1 - gamma), each = R)
    lst  <- rep(list(list('a' = rep(NA, R*9), 'b' = rep(NA, R), 
                          'c' = rep(NA, R*9))), nT)
    llobs <- 0
    
    for(t in 2:nT){
      # log density of the nA age groups, summed, per region and state
      lp <- colSums(matrix(lfd[t + nT*(1:nA-1), ], nrow = nA))
      kk <- pmax(lp[1:R], lp[(R+1):(2*R)], lp[(2*R+1):(3*R)])      # scaling per region
      fw <- exp(lp - rep(kk, times = 3))
      
      if(t == 2){
        a <- rep(fw, times = 3) * P[t, ] *
          c(rep(gvec[1:R], times = 3), rep(gvec[(R+1):(2*R)], times = 3),
            rep(gvec[(2*R+1):(3*R)], times = 3))
      } else {
        sl <- lst[[t-1]][['c']]
        a  <- rep(fw, times = 3) * P[t, ] * (sl[indl0] + sl[indl1] + sl[indl2])
      }
      b <- rowSums(matrix(a, nrow = R))
      llobs <- llobs + sum(log(b) + kk)
      b[b == 0] <- 10^(-300)
      lst[[t]][['a']] <- a
      lst[[t]][['b']] <- b
      lst[[t]][['c']] <- a/rep(b, times = 9)
    }
    
    js.prob <- do.call('rbind', lapply(1:length(lst), function(x) lst[[x]][['c']]))
    ms.prob <- rbind(c(sapply(1:R, function(j) sum(js.prob[2, c(j, R+j, 2*R+j)])),
                       sapply(1:R, function(j) sum(js.prob[2, c(3*R+j, 4*R+j, 5*R+j)])),
                       sapply(1:R, function(j) sum(js.prob[2, c(6*R+j, 7*R+j, 8*R+j)]))),
                     cbind(sapply(1:R, function(j) rowSums(js.prob[-1, c(j, 3*R+j, 6*R+j)])),
                           sapply(1:R, function(j) rowSums(js.prob[-1, c(R+j, 4*R+j, 7*R+j)])),
                           sapply(1:R, function(j) rowSums(js.prob[-1, c(2*R+j, 5*R+j, 8*R+j)]))))
    
    ### M-step
    
    # Optimize gamma
    fit0  <- optim(par     = gamma, 
                   fn      = obj.gamma, 
                   msprobs = ms.prob[1, ], 
                   R       = R,
                   method  = 'L-BFGS-B', 
                   lower   = 0.001, 
                   upper   = 0.999,
                   control = list(factr = 10^(-14)))
    gamma <- fit0$par
    
    # Define weights
    df$.w1 <- as.vector(do.call('rbind', replicate(nA, ms.prob[,(R+1):(2*R)], simplify = FALSE)))
    df$.w2 <- as.vector(do.call('rbind', replicate(nA, ms.prob[,(2*R+1):(3*R)], simplify = FALSE)))
    
    # Optimize alpha_1
    alpha.1 <- glm(form1, 
                   family  = poisson(link = 'log'), 
                   data    = df,
                   start   = alpha.1,
                   control = list(epsilon = 10^(-12)), 
                   weights = .w1)$coefficients
    
    # Optimize alpha_2
    alpha.2 <- glm(form2, 
                   family  = poisson(link = 'log'), 
                   data    = df, 
                   start   = alpha.2,
                   control = list(epsilon = 10^(-12)), 
                   weights = .w2)$coefficients
    
    # Transition parameters
    fit3 <- optim(par     = c(beta.0, beta.11, beta.22), 
                  fn      = obj.beta, 
                  jsprobs = js.prob,
                  xmat01  = xx01, 
                  xmat02  = xx02, 
                  xmat11  = xx11, 
                  xmat22  = xx22,
                  pena01  = pen01, 
                  pena02  = pen02, 
                  pena11  = pen11, 
                  pena22  = pen22, 
                  R       = R, 
                  method  = 'BFGS', 
                  control = list(reltol = 10^(-12)))
    
    i01     <- length(beta.01)
    i02     <- length(beta.02)
    i11     <- length(beta.11)
    
    beta.01 <- fit3$par[1:i01]
    beta.02 <- fit3$par[(i01+1):(i01+i02)]
    beta.11 <- fit3$par[(i01+i02+1):(i01+i02+i11)]
    beta.22 <- fit3$par[(i01+i02+i11+1):length(fit3$par)]
    beta.0  <- c(beta.01, beta.02)
    
    # Log-likelihood
    ll1 <- - llobs
    if(verbose) cat('iteration', iter, ': -loglik =', round(ll1, 4),
                    '| change =', signif(ll0 - ll1, 3), '\n')
    
    if(!is.null(best) && ll1 > best$ll){
      if(verbose) cat('no further improvement: keeping the previous iterate\n')
      alpha.1 <- best$alpha.1; alpha.2 <- best$alpha.2; gamma <- best$gamma
      beta.01 <- best$beta.01; beta.02 <- best$beta.02
      beta.11 <- best$beta.11; beta.22 <- best$beta.22
      ms.prob <- best$ms.prob; js.prob <- best$js.prob; ll1 <- best$ll
      break
    }
    best <- list(alpha.1 = alpha.1, alpha.2 = alpha.2, gamma = gamma,
                 beta.01 = beta.01, beta.02 = beta.02, beta.11 = beta.11,
                 beta.22 = beta.22, ms.prob = ms.prob, js.prob = js.prob, ll = ll1)
  }
  
  npar <- dim(x1)[2] + dim(x2)[2] + dim(xx01)[2] + dim(xx02)[2] +
    dim(xx11)[2] + dim(xx22)[2] + 1
  
  list(alpha.1 = alpha.1, 
       alpha.2 = alpha.2, 
       gamma   = gamma,
       beta.01 = beta.01, 
       beta.02 = beta.02, 
       beta.11 = beta.11, 
       beta.22 = beta.22,
       ms.prob = ms.prob, 
       js.prob = js.prob, 
       ll      = ll1, 
       npar    = npar,
       bic     = 2*ll1 + log(nrow(df))*npar, 
       itr     = iter, 
       cvg     = iter < max.iter)
}

###### (5) Simulate regime trajectories and death counts ----
# Their Section 5: state uncertainty (simulated trajectories) and Poisson uncertainty.
simulate_rs <- function(obj, xnew1, xnew2, xx01n, xx02n, xx11n, xx22n,
                        pen01n, pen02n, bnew, s0, nT, R, nA){
  x1b01 <- pen01n + sapply(1:R, function(j) t(t(xx01n[,,j])) %*% obj$beta.01)
  x2b02 <- pen02n + sapply(1:R, function(j) t(t(xx02n[,,j])) %*% obj$beta.02)
  x1b11 <- sapply(1:R, function(j) t(t(xx11n[,,j])) %*% obj$beta.11)
  x2b22 <- sapply(1:R, function(j) t(t(xx22n[,,j])) %*% obj$beta.22)
  
  p01 <- exp(x1b01)/(1 + exp(x1b01) + exp(x2b02))
  p02 <- exp(x2b02)/(1 + exp(x1b01) + exp(x2b02))
  p11 <- exp(x1b11)/(1 + exp(x1b11))
  p22 <- exp(x2b22)/(1 + exp(x2b22))
  
  # Regime trajectory per region
  s <- matrix(NA, nT, R)
  for(r in 1:R){
    st <- s0[r]
    for(t in 1:nT){
      pr <- switch(st + 1,
                   c(1 - p01[t,r] - p02[t,r], p01[t,r], p02[t,r]),
                   c(1 - p11[t,r], p11[t,r], 0),
                   c(1 - p22[t,r], 0, p22[t,r]))
      st      <- sample(0:2, 1, prob = pmax(pr, 0))
      s[t, r] <- st
    }
  }
  
  # Expected deaths given the sampled states, then Poisson noise
  x1a1 <- sapply(1:R, function(j) t(t(xnew1[,,j])) %*% obj$alpha.1)
  x2a2 <- sapply(1:R, function(j) t(t(xnew2[,,j])) %*% obj$alpha.2)
  amp  <- matrix(0, nT*nA, R)
  for(r in 1:R){
    sr <- rep(s[, r], times = nA)
    amp[sr == 1, r] <- x1a1[sr == 1, r]
    amp[sr == 2, r] <- x2a2[sr == 2, r]
  }
  mu <- bnew * exp(amp)
  
  list(states = s, mu = mu, deaths = matrix(rpois(length(mu), mu), nT*nA, R))
}