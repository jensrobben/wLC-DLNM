##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'spdep',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'data.table', 'xgboost',
              'parallel', 'doParallel', 'readxl', 'tidyr', 'purrr', 'imputeTS',
              'dlnm','Matrix', 'metR', 'pbapply', 'scoringRules', 'corpcor')
invisible(sapply(packages, require, character.only = TRUE))
`%notin%` <- Negate(`%in%`)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Colour scales
colors <- c('11' = "#A6CEE3", '84' = "#1F78B4", '75' = "#B2DF8A", '76' = "#33A02C", 
            '32' = "#FB9A99", '44' = "#E31A1C", '93' = "#FDBF6F", '53' = "#FF7F00", 
            '28' = "#CAB2D6", '52' = "#6A3D9A", '27' = "#FFFF99", '24' = "#B15928")

# Source
source('F. Weekly Lee-Carter (NB).R')
source('F. Weekly DLNM (NB).R')
source('F. Weekly Lee-Carter + DLNM (NB).R')

##### 1) Load data ----
# Load data
Df <- fread(file = '../Results updated/Df_var.txt')
Df$ia100 <- Df$ia100/1000

# Read shape file - downloaded from Eurostat
shapefile <- read_sf('../NUTS/NUTS_RG_03M_2024_3035.shp')
nuts.spec <- 1
shapef    <- shapefile[shapefile$LEVL_CODE == 1 & 
                         shapefile$CNTR_CODE %in% c('FR') &
                         ! shapefile$NUTS_ID %in% c('FRY', 'FRM'),] %>%
  dplyr::arrange(NUTS_NAME)
shapef    <- shapef %>% 
  dplyr::left_join(unique(Df[,c('region', 'NUTS_ID')])) %>%
  arrange(region)

# Select gender, ages, and time period
DfM <- Df %>% 
  dplyr::filter(gender == 2,
                age >= 50, 
                ISOYear >= 1990,
                ISOYear <= 2014) %>% 
  dplyr::arrange(age, ISODate)

# Avoid numerical issues
DfM$deaths[DfM$deaths == 0] <- 0.001

# Age range, year range, week range
xv      <- as.integer(sort(unique(DfM$age)))
isowv   <- isoweek(as.Date(sort(unique(DfM$ISODate))))
isoyv   <- isoyear(as.Date(sort(unique(DfM$ISODate)))) - min(DfM$ISOYear) + 1
regions <- as.character(sort(unique(DfM$region)))
R       <- length(regions)

# Adjacency matrix and penalty matrix
W     <- nb2mat(spdep::poly2nb(shapef, queen = TRUE), 
                style = 'B', zero.policy = TRUE)
D     <- diag(rowSums(W))
L     <- D - W

# Anomalies
sum <- DfM %>%
  dplyr::group_by(ISOWeek, region) %>%
  reframe('mtg' = quantile(tg, 0.5),
          'mia100' = quantile(ia100,0.90)) %>%
  as.data.table()
DfM <- DfM %>%
  dplyr::left_join(sum, by = c('ISOWeek', 'region'))
DfM$ia100 <- DfM$ia100 - DfM$mia100

##### 2) Cross-basis matrices ----
# Data
tavg <- unique(DfM[,c('ISODate','region', 'tg')]) %>%
  dplyr::arrange(ISODate) %>%
  as.data.table()
infl <- unique(DfM[,c('ISODate','region', 'ia100')]) %>%
  dplyr::arrange(ISODate) %>%
  as.data.table()

# Global bounds
bnd.tavg <- tavg %>%
  dplyr::summarise(
    min = floor(min(tg)),
    max = ceiling(max(tg))
  ) %>%
  unlist(use.names = FALSE)

# Maximum lags
mxlag.tavg <- 4
mxlag.infl <- 8

# Create arguments per region
args <- lapply(regions, function(r) {
  # Temperature
  argvar.tavg <- list(fun = 'bs', 
                      degree = 3,
                      knots = quantile(tavg[, tg], c(0.1,0.9)),
                      Boundary.knots = bnd.tavg)
  arglag.tavg <- list(fun   = "ns",
                      intercept = TRUE,
                      knots = dlnm::logknots(mxlag.tavg, nk = 2),
                      Boundary.knots = c(0,mxlag.tavg))
  
  # Influenza
  argvar.infl <- list(fun = 'poly', 
                      degree = 1,
                      scale = max(infl[region == r,ia100]))
  arglag.infl <- list(fun = "ns",
                      intercept = TRUE,
                      knots = dlnm::logknots(mxlag.infl, nk = 2),
                      Boundary.knots = c(0,mxlag.infl))
  
  list('Temp' = list(argvar = argvar.tavg, arglag = arglag.tavg),
       'Infl' = list(argvar = argvar.infl, arglag = arglag.infl))
})
names(args) <- regions

# Function to create cross-basis matrix
cb_fun <- function(dat, lags, argvar, arglag){
  # Append data to argument lists
  argvar <- append(argvar, list('x' = dat))
  arglag <- append(arglag, list('x' = lags))
  
  # One-basis functions
  Z <- do.call(onebasis, argvar)
  C <- do.call(onebasis, arglag)
  
  # Construct cross-basis matrix
  lag <- c(min(lags), max(lags))
  cb  <- matrix(0, nrow = dim(Z)[1], ncol = ncol(Z) * ncol(C))
  
  for (v in seq(length = ncol(Z))) {
    mat <- as.matrix(tsModel::Lag(Z[, v], dlnm:::seqlag(lag)))
    
    for (l in seq(length = ncol(C))) {
      cb[, ncol(C) * (v - 1) + l] <- mat %*% (C[, l])
    }
  }
  
  cb
}
cb_int <- function(dat1, dat2, lags, argvar1, argvar2, arglag){
  # Append data to argument lists
  argvar1 <- append(argvar1, list('x' = dat1))
  argvar2 <- append(argvar2, list('x' = dat2))
  arglag  <- append(arglag,  list('x' = lags))
  
  # One-basis functions
  Z1  <- do.call(onebasis, argvar1)
  Z2  <- do.call(onebasis, argvar2)
  Z12 <- t(sapply(1:dim(Z1)[1], function(t) 
    as.vector(outer(Z1[t, ], Z2[t, ]))))
  C   <- do.call(onebasis, arglag)
  
  # Construct cross-basis matrix
  lag <- c(min(lags), max(lags))
  cb  <- matrix(0, nrow = dim(Z12)[1], ncol = ncol(Z12) * ncol(C))
  
  for (v in seq(length = ncol(Z12))) {
    mat <- as.matrix(tsModel::Lag(Z12[, v], dlnm:::seqlag(lag)))
    
    for (l in seq(length = ncol(C))) {
      cb[, ncol(C) * (v - 1) + l] <- mat %*% (C[, l])
    }
  }
  
  cb
}

# Result
cb1  <- simplify2array(lapply(regions, function(r) {
  argv <- args[[r]]$Temp$argvar
  argl <- args[[r]]$Temp$arglag
  
  cb_fun(dat = tavg[region == r, tg], lags = 0:mxlag.tavg,
         argvar = argv, arglag = argl)
}))
cb2  <- simplify2array(lapply(regions, function(r) {
  argv <- args[[r]]$Infl$argvar
  argl <- args[[r]]$Infl$arglag
  
  cb_fun(dat = infl[region == r, ia100], lags = 0:mxlag.infl,
         argvar = argv, arglag = argl)
}))
cb12 <- simplify2array(lapply(regions, function(r) {
  argv1 <- args[[r]]$Temp$argvar
  argv2 <- args[[r]]$Infl$argvar
  argl  <- args[[r]]$Temp$arglag
  
  cb_int(dat1 = tavg[region == r, tg], 
         dat2 = infl[region == r, ia100], 
         lags = 0:min(c(mxlag.tavg,mxlag.infl)),
         argvar1 = argv1,
         argvar2 = argv2,
         arglag = argl)
}))

cb1    <- center_cb(cb1, isowv)
cb2    <- center_cb(cb2, isowv)
cb12   <- center_cb(cb12, isowv)
means1 <- attr(cb1, 'means')
means2 <- attr(cb2, 'means')

##### 3) New data ----
# Updated dataset
Dfnew <- fread(file = '../Results updated/Df_var.txt') %>%
  dplyr::filter(ISOYear >= 2015,
                ISOYear <= 2019,
                gender == 2,
                age >= 50)

# Exposure
expo.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'expo')) %>%
  unique() %>%
  dplyr::arrange(age, region, ISODate)
expo.obs <- aperm(simplify2array(lapply(unique(expo.obs$age), function(a)
  matrix(expo.obs[expo.obs$age == a,]$expo, 
         ncol = length(unique(expo.obs$region)), 
         nrow = length(unique(expo.obs$ISODate)),
         dimnames = list(as.character(unique(expo.obs$ISODate)), 
                         unique(expo.obs$region))))), c(1,3,2))
dimnames(expo.obs)[[2]] <- unique(Dfnew$age)

# Deaths
deaths.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'deaths')) %>%
  unique() %>%
  dplyr::arrange(age, region, ISODate)
deaths.obs <- aperm(simplify2array(lapply(unique(deaths.obs$age), function(a)
  matrix(deaths.obs[deaths.obs$age == a,]$deaths, 
         ncol = length(unique(deaths.obs$region)), 
         nrow = length(unique(deaths.obs$ISODate)),
         dimnames = list(as.character(unique(deaths.obs$ISODate)), 
                         unique(deaths.obs$region))))),
  c(1,3,2))
dimnames(deaths.obs)[[2]] <- unique(Dfnew$age)


##### 4) Run models ----

# Data preparation
dtxr <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM)
etxr <- stats::xtabs(expo ~ ISODate + age + region, data = DfM)

dtxr_ft  <- deaths.obs
etxr_ft  <- expo.obs
dimnames(dtxr_ft) = dimnames(etxr_ft) <- 
  list(as.character(sort(unique(Dfnew$ISODate))),
       unique(DfM$age), unique(DfM$region))

# Optimal lambda
lambda.opt <- c('Temp' = 10^(1.75), 'Infl' = 10^(3.75))

# Fit
# lca1 <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                         etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                         cb = list('Temp' = cb1, 'Infl' = cb2)[1],
#                         lambda = lambda.opt[1])
# lca2  <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                          etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                          cb = list('Temp' = cb1, 'Infl' = cb2)[2],
#                          lambda = lambda.opt[2])
# lca12 <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                          etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                          cb = list('Temp' = cb1, 'Infl' = cb2),
#                          lambda = lambda.opt)
# saveRDS(lca1,  file = '../Results updated/lca_Temp_benchmark.rds')
# saveRDS(lca2,  file = '../Results updated/lca_Infl_benchmark.rds')
# saveRDS(lca12, file = '../Results updated/lca_Temp_Infl_benchmark.rds')

lca1 <- readRDS('../Results updated/lca_Temp_benchmark.rds')
lca2 <- readRDS('../Results updated/lca_Infl_benchmark.rds')
lca3 <- readRDS('../Results updated/lca_Temp_Infl_benchmark.rds')

##### 5) Point Forecast -----

# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Recover original temperature/influenza values
hist.tg <- sapply(1:length(regions), function(r) {
  tavg <- unique(DfM[,c('ISODate','region', 'tg')]) %>%
    dplyr::arrange(ISODate) %>%
    as.data.table()
  tavg[region == regions[r], tg]
}) 
hist.ia <- sapply(1:length(regions), function(r) {
  infl <- unique(DfM[,c('ISODate','region', 'ia100')]) %>%
    dplyr::arrange(ISODate) %>%
    as.data.table()
  infl[region == regions[r], ia100]
}) 

# Cross-basis for new covariate values: prepend the last weeks of 1990-2014 (lags),
# drop them again, and centre with the calibration means
cb_forecast <- function(new, hist, var, mxlag, means){
  cb <- simplify2array(lapply(seq_along(regions), function(j) {
    a <- args[[regions[j]]][[var]]
    x <- c(tail(hist[,j], mxlag), new[,j])
    cb_fun(dat = x, lags = 0:mxlag, argvar = a$argvar,
           arglag = a$arglag)[-(1:mxlag),,drop = FALSE]
  }))
  dimnames(cb)[[3]] <- regions
  center_cb(cb, isowv_fc, means = means)
}

# Temperature
tg.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'tg')) %>%
  unique() %>%
  dplyr::arrange(region, ISODate)
tg.obs <- matrix(tg.obs$tg, 
                 ncol = length(regions), 
                 nrow = length(unique(tg.obs$ISODate)),
                 dimnames = list(as.character(unique(tg.obs$ISODate)), 
                                 unique(tg.obs$region)))

# Influenza
ia.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'ia100')) %>%
  unique() %>%
  dplyr::arrange(region, ISODate)
ia.obs$ia100 <- ia.obs$ia100/1000
ia.obs <- ia.obs %>%
  dplyr::left_join(sum[,c('ISOWeek', 'region', 'mia100')], by = c('ISOWeek', 'region'))
ia.obs$ia100 <- ia.obs$ia100 - ia.obs$mia100
ia.obs <- matrix(ia.obs$ia100, 
                 ncol = length(regions), 
                 nrow = length(unique(ia.obs$ISODate)),
                 dimnames = list(as.character(unique(ia.obs$ISODate)), 
                                 unique(ia.obs$region)))

# Fit ARIMAX for period effect
obj.kt <- lapply(1:length(regions), function(r){
  list('Temp' = forecast::Arima(ts(lca1$lcd$kappa2[,r]), 
                                order = c(0,1,1),
                                include.drift = TRUE),
       'Infl' = forecast::Arima(ts(lca2$lcd$kappa2[,r]), 
                                order = c(0,1,1),
                                include.drift = TRUE),
       'Temp.Infl' = forecast::Arima(ts(lca3$lcd$kappa2[,r]), 
                                     order = c(0,1,1),
                                     include.drift = TRUE))
})

# Point forecasts
tg_fc <- tg.obs
ia_fc <- ia.obs
kt_fc <- list(
  'Temp' = sapply(1:dim(dtxr)[3], function(r) {
    c(forecast::forecast(obj.kt[[r]]$Temp, nY)$mean)}),
  'Infl' = sapply(1:dim(dtxr)[3], function(r) {
    c(forecast::forecast(obj.kt[[r]]$Infl, nY)$mean)}),
  'Temp.Infl' = sapply(1:dim(dtxr)[3], function(r) {
    c(forecast::forecast(obj.kt[[r]]$Temp.Infl, nY)$mean)})
)

# Cross-bases
cb1_fc <- cb_forecast(tg_fc, hist.tg, 'Temp', mxlag.tavg, means1)
cb2_fc <- cb_forecast(ia_fc, hist.ia, 'Infl', mxlag.infl, means2)

# DLNM
dlnm <- list(
  'Temp' = Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out}, list('Temp' = cb1_fc), lca1$lcd$gammal),
  'Infl' = Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out}, list('Infl' = cb2_fc), lca2$lcd$gammal),
  'Temp.Infl' = Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out}, list('Temp' = cb1_fc, 'Infl' = cb2_fc), lca3$lcd$gammal)
)

# Force of mortality
mhat1 <- exp(sweep(aperm(
  lca1$lcd$beta2 %o% kt_fc$Temp[isoyv_fc - min(isoyv_fc) + 1,] +
    lca1$lcd$beta3 %o% lca1$lcd$kappa3[isowv_fc,] +
    Reduce(`+`, Map(function(b, d) b %o% d, lca1$lcd$betal, dlnm$Temp)),
                        c(2,1,3)), c(2,3), lca1$lcd$beta1, '+'))
mhat2 <- exp(sweep(aperm(
  lca2$lcd$beta2 %o% kt_fc$Infl[isoyv_fc - min(isoyv_fc) + 1,] +
    lca2$lcd$beta3 %o% lca2$lcd$kappa3[isowv_fc,] +
    Reduce(`+`, Map(function(b, d) b %o% d, lca2$lcd$betal, dlnm$Infl)),
  c(2,1,3)), c(2,3), lca2$lcd$beta1, '+'))

mhat3 <- exp(sweep(aperm(
  lca3$lcd$beta2 %o% kt_fc$Temp.Infl[isoyv_fc - min(isoyv_fc) + 1,] +
    lca3$lcd$beta3 %o% lca3$lcd$kappa3[isowv_fc,] +
    Reduce(`+`, Map(function(b, d) b %o% d, lca3$lcd$betal, dlnm$Temp.Infl)),
  c(2,1,3)), c(2,3), lca3$lcd$beta1, '+'))

dimnames(mhat1) = dimnames(mhat2) = dimnames(mhat3) <- 
  list(c(as.character(sort(unique(Dfnew$ISODate)))),
       dimnames(lca1$lcd$mhat)[[2]],
       dimnames(lca1$lcd$mhat)[[3]])

# Error metrics
dt <- dtxr_ft
et <- etxr_ft
mhat_rep1 <- mhat1
mhat_rep2 <- mhat2
mhat_rep3 <- mhat3
diff1 <- dt - mhat1 * et
diff2 <- dt - mhat2 * et
diff3 <- dt - mhat3 * et

plot(as.Date(rownames(dt)), (dt/et)[,9,1], type = 'l')
lines(as.Date(rownames(mhat1)), mhat1[,9,1], col = 'orange', lwd = 2)
lines(as.Date(rownames(mhat2)), mhat2[,9,1], col = 'steelblue', lwd = 2)
lines(as.Date(rownames(mhat3)), mhat3[,9,1], col = 'forestgreen', lwd = 2)

rmse <- c('Temp'      = sqrt(1/length(diff1) * sum(((diff1)^2))),
          'Infl'      = sqrt(1/length(diff2) * sum(((diff2)^2))),
          'Temp.Infl' = sqrt(1/length(diff3) * sum(((diff3)^2))))
mae  <- c('Temp'      = 1/length(diff1) * sum(abs(diff1)),
          'Infl'      = 1/length(diff2) * sum(abs(diff2)),
          'Temp.Infl' = 1/length(diff3) * sum(abs(diff3)))

rmse_age <- data.frame(
  'Temp'      = sapply(1:dim(diff1)[2], function(a) 
    sqrt(1/length(diff1[,a,]) * sum(((diff1[,a,])^2)))),
  'Infl'      = sapply(1:dim(diff2)[2], function(a) 
    sqrt(1/length(diff2[,a,]) * sum(((diff2[,a,])^2)))),
  'Temp.Infl' = sapply(1:dim(diff3)[2], function(a) 
    sqrt(1/length(diff3[,a,]) * sum(((diff3[,a,])^2)))))
mae_age <- data.frame(
  'Temp'      = sapply(1:dim(diff1)[2], function(a) 
    1/length(diff1[,a,]) * sum(abs(diff1[,a,]))),
  'Infl'      = sapply(1:dim(diff2)[2], function(a) 
    1/length(diff2[,a,]) * sum(abs(diff2[,a,]))),
  'Temp.Infl' = sapply(1:dim(diff3)[2], function(a) 
    1/length(diff3[,a,]) * sum(abs(diff3[,a,]))))

# saveRDS(rmse_age, '../Results updated/RMSE_age_Models4-6.rds')
# saveRDS(mae_age, '../Results updated/MAE_age_Models4-6.rds')
# saveRDS(rmse, '../Results updated/RMSE_overall_Models4-6.rds')
# saveRDS(mae, '../Results updated/MAE_overall_Models4-6.rds')

##### 6) Simulations covariates/kappa ----

### Define time series 

# Quantiles influenza
ia.quant <- sum[,c('ISOWeek', 'region', 'mia100')] %>%
  dplyr::arrange(region, ISOWeek)
ia.quant <- matrix(ia.quant$mia100, ncol = R, nrow = length(unique(ia.quant$ISOWeek)),
                   dimnames = list(as.character(unique(ia.quant$ISOWeek)), 
                                   unique(ia.quant$region)))
ia.quant <- ia.quant[isowv_fc,]

# Time series per region - transform ia
wv <- sapply(split(isowv, isoyv), function(x) max(x))

lst.ts <- lapply(1:length(regions), function(r){
  data.frame(
    'kappa1' = rep(lca1$lcd$kappa2[,r], times = wv), 
    'kappa2' = rep(lca2$lcd$kappa2[,r], times = wv),
    'kappa3' = rep(lca3$lcd$kappa2[,r], times = wv),
    'Temp'   = tavg %>% dplyr::filter(region == regions[r]) %>% pull(tg),
    'Infl'   = log(infl %>% dplyr::filter(region == regions[r]) %>% pull(ia100) + 
                     DfM %>% dplyr::filter(region == regions[r],
                                           age == '50') %>% pull(mia100) + 0.01))
})
names(lst.ts) <- regions

# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

### Shrinkage-based copulas -> towards constant correlation
shrink_corr <- function(u) {
  d <- ncol(u)
  n <- nrow(u)
  Rt <- sin(pi / 2 * cor(u, method = "kendall"))

  X  <- scale(qnorm(u))
  S  <- cor(X)
  sb <- mean(S[upper.tri(S)])
  
  num <- den <- 0
  for (i in 1:(d - 1)) for (j in (i + 1):d) {
    w   <- X[, i] * X[, j]
    num <- num + n / (n - 1)^3 * sum((w - mean(w))^2)
    den <- den + (S[i, j] - sb)^2
  }
  
  lam <- min(1, max(0, num / den))
  Tm  <- matrix(mean(Rt[upper.tri(Rt)]), d, d); diag(Tm) <- 1
  Rs  <- lam * Tm + (1 - lam) * Rt
  list(R = as.matrix(Matrix::nearPD(Rs, corr = TRUE)$mat), lambda = lam)
}

fit_shrink_tcop <- function(u, df.start = 4) {
  s   <- shrink_corr(u)
  d   <- ncol(u)
  cop <- copula::tCopula(copula::P2p(s$R), dim = d, dispstr = "un", df = df.start)
  copula::fixedParam(cop) <- c(rep(TRUE, d * (d - 1) / 2), FALSE)
  fit <- copula::fitCopula(cop, data = u, method = "ml", start = df.start)
  list(fit = fit, lambda = s$lambda)
}

fit_shrink_normalcop <- function(u) {
  s   <- shrink_corr(u)
  d   <- ncol(u)
  cop <- copula::normalCopula(copula::P2p(s$R), dim = d, dispstr = "un")
  fit <- new("fitCopula", copula = cop, estimate = numeric(0), 
             var.est = matrix(numeric(0), 0, 0),
             loglik = sum(copula::dCopula(u, cop, log = TRUE)), nsample = nrow(u),
             method = "shrinkage", call = match.call(), fitting.stats = list(convergence = 0L))
  list(fit = fit, lambda = s$lambda)
}

### Temperature models

# Individual SARIMA - models
sarima.tg = res.tg <- list()
for(r in 1:length(regions)){
  cat(regions[r], ' ... ')
  sarima.tg[[r]] <- forecast::Arima(ts(lst.ts[[r]]$Temp, frequency = 52), 
                                    order = c(1,0,1),
                                    seasonal = c(1,1,1),
                                    include.drift = TRUE)
  res.tg[[r]] <- residuals(sarima.tg[[r]], type = 'innovation')/
    sqrt(sarima.tg[[r]]$sigma2)
}
names(sarima.tg) <- regions

# Pseudo-observations of standardized residuals
res.mat.tg <- do.call(cbind, res.tg)[-c(1:52),]
u.mat.tg   <- copula::pobs(res.mat.tg)

# Fit copulas
copt.tg  <- copula::tCopula(dim = ncol(u.mat.tg), dispstr = 'ex')
copn.tg  <- copula::normalCopula(dim = ncol(u.mat.tg), dispstr = 'ex')
copts.tg <- fit_shrink_tcop(u.mat.tg)
copns.tg <- fit_shrink_normalcop(u.mat.tg)

fit.copt.tg  <- copula::fitCopula(copt.tg, data = u.mat.tg, method = "ml")
fit.copn.tg  <- copula::fitCopula(copn.tg, data = u.mat.tg, method = "ml")
fit.copts.tg <- copts.tg$fit
fit.copns.tg <- copns.tg$fit

# BIC comparison
d       <- ncol(u.mat.tg)
pt.eff   <- copts.tg$lambda + (1 - copts.tg$lambda) * (d * (d - 1) / 2) + 
  length(fit.copts.tg@estimate)
pn.eff   <- copns.tg$lambda + (1 - copns.tg$lambda) * (d * (d - 1) / 2) + 
  length(fit.copns.tg@estimate)

bic.t  <- -2 * fit.copt.tg@loglik + 
  log(nrow(u.mat.tg)) * length(fit.copt.tg@estimate)
bic.n  <- -2 * fit.copn.tg@loglik + 
  log(nrow(u.mat.tg)) * length(fit.copn.tg@estimate)
bic.ts <- -2 * fit.copts.tg@loglik + 
  log(nrow(u.mat.tg)) * pt.eff
bic.ns <- -2 * fit.copns.tg@loglik + 
  log(nrow(u.mat.tg)) * pn.eff

c(bic.t, bic.n, bic.ts, bic.ns)

# View summary of fitted copula
summary(fit.copts.tg)

### Influenza models

# Individual SARIMA - models 
sarima.ia = sarima.ext.ia = res.ia <- list()
for(r in 1:length(regions)){
  cat(regions[r], ' ... ')
  sarima.ia[[r]] <- forecast::Arima(ts(lst.ts[[r]]$Infl, frequency = 52), 
                                    order = c(1,0,1),
                                    seasonal = c(1,1,1), 
                                    include.drift = TRUE)
  
  sarima.ext.ia[[r]] <- forecast::Arima(ts(lst.ts[[r]]$Infl, frequency = 52), 
                                    order = c(1,0,1),
                                    seasonal = c(1,1,1), 
                                    xreg = lst.ts[[r]]$Temp,
                                    include.drift = TRUE)

  res.ia[[r]] <- residuals(sarima.ia[[r]], type = 'innovation')/
    sqrt(sarima.ia[[r]]$sigma2)
}
names(sarima.ia) = names(sarima.ext.ia) <- regions

# Pseudo-observations of standardized residuals
res.mat.ia <- do.call(cbind, res.ia)[-c(1:52),]
u.mat.ia   <- copula::pobs(res.mat.ia)

# Fit copulas
copt.ia  <- copula::tCopula(dim = ncol(u.mat.ia), dispstr = 'ex')
copn.ia  <- copula::normalCopula(dim = ncol(u.mat.ia), dispstr = 'ex')
copts.ia <- fit_shrink_tcop(u.mat.ia)
copns.ia <- fit_shrink_normalcop(u.mat.ia)

fit.copt.ia  <- copula::fitCopula(copt.ia, data = u.mat.ia, method = "ml")
fit.copn.ia  <- copula::fitCopula(copn.ia, data = u.mat.ia, method = "ml")
fit.copts.ia <- copts.ia$fit
fit.copns.ia <- copns.ia$fit

# BIC comparison
d       <- ncol(u.mat.ia)
pt.eff   <- copts.ia$lambda + (1 - copts.ia$lambda) * (d * (d - 1) / 2) + 
  length(fit.copts.ia@estimate)
pn.eff   <- copns.ia$lambda + (1 - copns.ia$lambda) * (d * (d - 1) / 2) + 
  length(fit.copns.ia@estimate)

bic.t  <- -2 * fit.copt.ia@loglik + 
  log(nrow(u.mat.ia)) * length(fit.copt.ia@estimate)
bic.n  <- -2 * fit.copn.ia@loglik + 
  log(nrow(u.mat.ia)) * length(fit.copn.ia@estimate)
bic.ts <- -2 * fit.copts.ia@loglik + 
  log(nrow(u.mat.ia)) * pt.eff
bic.ns <- -2 * fit.copns.ia@loglik + 
  log(nrow(u.mat.ia)) * pn.eff

c(bic.t, bic.n, bic.ts, bic.ns)

# View summary of fitted copula
summary(fit.copt.ia)


### Kappa_t models 

# Individual ARIMA - models (use tg, and ia as external regressor)
sarima.kt1 = sarima.kt2 = sarima.kt3 = sarima.ext.kt <- list()
res.kt1 = res.kt2 = res.kt3 <- list()

for(r in 1:length(regions)){
  cat(regions[r], ' ... ')
  ext <- cbind('tg' = tapply(lst.ts[[r]]$Temp, isoyv, mean),
               'ia' = tapply(exp(lst.ts[[r]]$Infl) - 0.01, isoyv, mean))
  
  sarima.kt1[[r]] <- forecast::Arima(ts(lca1$lcd$kappa2[,r]), 
                                     order = c(0,1,1),
                                     include.drift = TRUE)
  sarima.kt2[[r]] <- forecast::Arima(ts(lca2$lcd$kappa2[,r]), 
                                     order = c(0,1,1),
                                     include.drift = TRUE)
  sarima.kt3[[r]] <- forecast::Arima(ts(lca3$lcd$kappa2[,r]), 
                                     order = c(0,1,1),
                                     include.drift = TRUE)
  
  sarima.ext.kt[[r]] <- forecast::Arima(ts(lca3$lcd$kappa2[,r]), 
                                        order = c(0,1,1),
                                        include.drift = TRUE,
                                        xreg = ext)
  
  res.kt1[[r]] <- residuals(sarima.kt1[[r]], type = 'innovation')/
    sqrt(sarima.kt1[[r]]$sigma2)
  res.kt2[[r]] <- residuals(sarima.kt2[[r]], type = 'innovation')/
    sqrt(sarima.kt2[[r]]$sigma2)
  res.kt3[[r]] <- residuals(sarima.kt3[[r]], type = 'innovation')/
    sqrt(sarima.kt3[[r]]$sigma2)
}
names(sarima.kt1) = names(sarima.kt2) = names(sarima.kt3) <- regions
names(sarima.ext.kt) <- regions

# Pseudo-observations of standardized residuals
res.mat.kt1 <- do.call(cbind, res.kt1)[-c(1),]
res.mat.kt2 <- do.call(cbind, res.kt2)[-c(1),]
res.mat.kt3 <- do.call(cbind, res.kt3)[-c(1),]

u.mat.kt1   <- copula::pobs(res.mat.kt1)
u.mat.kt2   <- copula::pobs(res.mat.kt2)
u.mat.kt3   <- copula::pobs(res.mat.kt3)

# Fit t-copula
cop.kt1     <- copula::normalCopula(dim = ncol(u.mat.kt1))
cop.kt2     <- copula::normalCopula(dim = ncol(u.mat.kt2))
cop.kt3     <- copula::normalCopula(dim = ncol(u.mat.kt3))

fit.cop.kt1 <- copula::fitCopula(cop.kt1, data = u.mat.kt1, method = "ml")
fit.cop.kt2 <- copula::fitCopula(cop.kt2, data = u.mat.kt2, method = "ml")
fit.cop.kt3 <- copula::fitCopula(cop.kt3, data = u.mat.kt3, method = "ml")

# View summary of fitted copula
summary(fit.cop.kt1)
summary(fit.cop.kt2)
summary(fit.cop.kt3)

### Comparison SARIMA vs SARIMAX
a1 <- stack(lapply(sarima.ia, function(s) BIC(s)))
a2 <- stack(lapply(sarima.ext.ia, function(s) BIC(s)))
names(a1) <- c('V1', 'code')
names(a2) <- c('V2', 'code')

b1 <- stack(lapply(sarima.kt3, function(s) BIC(s)))
b2 <- stack(lapply(sarima.ext.kt, function(s) BIC(s)))
names(b1) <- c('V3', 'code')
names(b2) <- c('V4', 'code')

ddd <- list(a1, a2, b1, b2) %>%
  reduce(full_join, by = "code") %>%
  mutate(code = as.integer(as.character(code)), 
         .before = V1) %>%
  dplyr::mutate(region = shapef$NUTS_NAME[which(shapef$region %in% code)], 
                .before = 1)%>%
  dplyr::select(-c('code'))
head(ddd)
rownames(ddd) <- NULL


print(xtable::xtable(ddd), include.rownames = FALSE)


##### 7) Simulate trajectories ----
n.sim <- 5000
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl,c('n.sim', 'n.hor', 'regions', 'R', 'expo.obs', 'infl', 
                   'sarima.tg', 'sarima.ia', 'sarima.kt1', 'sarima.kt2', 
                   'sarima.kt3', 'ia.quant', 'tavg', 'args', 'lca1', 'lca2',
                   'lca3', 'mxlag.infl', 'mxlag.tavg', 'isowv_fc', 'isoyv_fc',
                   'fit.copts.tg', 'fit.copt.ia', 'fit.cop.kt1', 'cb_forecast', 
                   'center_cb', 'hist.tg', 'hist.ia', 'means1', 'means2', 
                   'cb1_fc', 'cb2_fc','regions',
                   'fit.cop.kt2', 'fit.cop.kt3', 'nY', 'cb_fun', 'xv'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(copula)
  library(dlnm)
  library(forecast)
  library(data.table)
})
results <- pblapply(cl=cl, X=1:n.sim, FUN = function(s) {
  ## Simulate temperature
  sim.u.tg <- copula::rCopula(n.hor, fit.copts.tg@copula)
  sim.z.tg <- qnorm(sim.u.tg)
  sim.tg   <- sapply(1:length(regions), function(r) {
    simulate(object = sarima.tg[[r]], nsim = n.hor,
             innov = sim.z.tg[,r] * sqrt(sarima.tg[[r]]$sigma2)) 
  })
  dimnames(sim.tg) <- list(as.character(rownames(expo.obs)), regions)
  
  ## Simulate influenza
  sim.u.ia <- copula::rCopula(n.hor, fit.copt.ia@copula)
  sim.z.ia <- qnorm(sim.u.ia)
  sim.ia   <- sapply(1:length(regions), function(r) {
    exp(simulate(object = sarima.ia[[r]], nsim = n.hor,
                 innov  = sim.z.ia[,r] * sqrt(sarima.ia[[r]]$sigma2))) - 0.01
  }) - ia.quant
  dimnames(sim.ia) <- list(as.character(rownames(expo.obs)), regions)
  
  ## Simulate kt-process - temperature
  sim.u.kt1 <- copula::rCopula(nY, fit.cop.kt1@copula)
  sim.z.kt1 <- qnorm(sim.u.kt1)
  sim.kt1   <- sapply(1:length(regions), function(r) {
    simulate(object = sarima.kt1[[r]], nsim = nY,
             innov = sim.z.kt1[,r] * sqrt(sarima.kt1[[r]]$sigma2))
  })
  
  ## Simulate kt-process - influenza
  sim.u.kt2 <- copula::rCopula(nY, fit.cop.kt2@copula)
  sim.z.kt2 <- qnorm(sim.u.kt2)
  sim.kt2   <- sapply(1:length(regions), function(r) {
    simulate(object = sarima.kt2[[r]], nsim = nY,
             innov = sim.z.kt2[,r] * sqrt(sarima.kt2[[r]]$sigma2))
  })
  
  ## Simulate kt-process - temperature + influenza
  sim.u.kt3 <- copula::rCopula(nY, fit.cop.kt3@copula)
  sim.z.kt3 <- qnorm(sim.u.kt3)
  sim.kt3   <- sapply(1:length(regions), function(r) {
    simulate(object = sarima.kt3[[r]], nsim = nY,
             innov = sim.z.kt3[,r] * sqrt(sarima.kt3[[r]]$sigma2))
  })
  
  ## Cross-basis temperature & influenza - observed
  obs.cb1 <- cb1_fc
  obs.cb2 <- cb2_fc
  
  ## Cross-basis temperature & influenza - simulated
  sim.cb1 <- cb_forecast(sim.tg, hist.tg, 'Temp', mxlag.tavg, means1)
  sim.cb2 <- cb_forecast(sim.ia, hist.ia, 'Infl', mxlag.infl, means2)
  
  # DLNM - observed
  obs.dlnm <- list(
    'Temp' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Temp' = obs.cb1), lca1$lcd$gammal),
    'Infl' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Infl' = obs.cb2), lca2$lcd$gammal),
    'Temp.Infl' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Temp' = obs.cb1, 'Infl' = obs.cb2), lca3$lcd$gammal)
  )
  
  # DLNM - simulated
  sim.dlnm <- list(
    'Temp' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Temp' = sim.cb1), lca1$lcd$gammal),
    'Infl' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Infl' = sim.cb2), lca2$lcd$gammal),
    'Temp.Infl' = Map(function(c, g){
      out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
      out[is.na(out)] <- 0
      out
    }, list('Temp' = sim.cb1, 'Infl' = sim.cb2), lca3$lcd$gammal)
  )
  
  # Force of mortality - observed
  ijk <- isoyv_fc - min(isoyv_fc) + 1
  
  obs.mu <- list(
    'Temp' = exp(sweep(aperm(lca1$lcd$beta2 %o% sim.kt1[ijk,] + 
                               lca1$lcd$beta3 %o% lca1$lcd$kappa3[isowv_fc,] + 
                               Reduce(`+`, Map(function(b, d) 
                                 b %o% d, lca1$lcd$betal, obs.dlnm$Temp)),
                             c(2,1,3)), c(2,3), lca1$lcd$beta1, '+')),
    'Infl' = exp(sweep(aperm(lca2$lcd$beta2 %o% sim.kt2[ijk,] + 
                               lca2$lcd$beta3 %o% lca2$lcd$kappa3[isowv_fc,] + 
                               Reduce(`+`, Map(function(b, d) 
                                 b %o% d, lca2$lcd$betal, obs.dlnm$Infl)),
                             c(2,1,3)), c(2,3), lca2$lcd$beta1, '+')),
    'Temp.Infl' = exp(sweep(aperm(lca3$lcd$beta2 %o% sim.kt3[ijk,] + 
                                    lca3$lcd$beta3 %o% lca3$lcd$kappa3[isowv_fc,] +
                                    Reduce(`+`, Map(function(b, d) 
                                      b %o% d, lca3$lcd$betal, obs.dlnm$Temp.Infl)),
                                  c(2,1,3)), c(2,3), lca3$lcd$beta1, '+'))
  )
  
  # Force of mortality - simulated
  ijk <- isoyv_fc - min(isoyv_fc) + 1
  
  sim.mu <- list(
    'Temp' = exp(sweep(aperm(lca1$lcd$beta2 %o% sim.kt1[ijk,] + 
                               lca1$lcd$beta3 %o% lca1$lcd$kappa3[isowv_fc,] + 
                               Reduce(`+`, Map(function(b, d) 
                                 b %o% d, lca1$lcd$betal, sim.dlnm$Temp)),
                             c(2,1,3)), c(2,3), lca1$lcd$beta1, '+')),
    'Infl' = exp(sweep(aperm(lca2$lcd$beta2 %o% sim.kt2[ijk,] + 
                               lca2$lcd$beta3 %o% lca2$lcd$kappa3[isowv_fc,] + 
                               Reduce(`+`, Map(function(b, d) 
                                 b %o% d, lca2$lcd$betal, sim.dlnm$Infl)),
                             c(2,1,3)), c(2,3), lca2$lcd$beta1, '+')),
    'Temp.Infl' = exp(sweep(aperm(lca3$lcd$beta2 %o% sim.kt3[ijk,] + 
                                    lca3$lcd$beta3 %o% lca3$lcd$kappa3[isowv_fc,] +
                                    Reduce(`+`, Map(function(b, d) 
                                      b %o% d, lca3$lcd$betal, sim.dlnm$Temp.Infl)),
                             c(2,1,3)), c(2,3), lca3$lcd$beta1, '+'))
  )
  
  # Dispersion 
  exp.ds <- list(
    'Temp' = aperm(replicate(n.hor, 
                             exp(outer(lca1$lcd$phix, lca1$lcd$phir, "+")),
                             simplify = 'array'), c(3, 1, 2)),
    'Infl' = aperm(replicate(n.hor, 
                             exp(outer(lca2$lcd$phix, lca2$lcd$phir, "+")),
                             simplify = 'array'), c(3, 1, 2)),
    'Temp.Infl' = aperm(replicate(n.hor, 
                                  exp(outer(lca3$lcd$phix, lca3$lcd$phir, "+")),
                                  simplify = 'array'), c(3, 1, 2))
  )
  
  len <- n.hor * length(xv) * R
  dms <- c(n.hor, length(xv), R)
  
  list('O.A1' = obs.mu$Temp * expo.obs,
       'O.A2' = obs.mu$Infl * expo.obs,
       'O.A3' = obs.mu$Temp.Infl * expo.obs,
       'S.A1' = sim.mu$Temp * expo.obs,
       'S.A2' = sim.mu$Infl * expo.obs,
       'S.A3' = sim.mu$Temp.Infl * expo.obs,
       'O.B1' = array(rnbinom(len, mu = obs.mu$Temp * expo.obs, 
                              size = exp.ds$Temp), dim = dms),
       'O.B2' = array(rnbinom(len, mu = obs.mu$Infl * expo.obs, 
                              size = exp.ds$Infl), dim = dms),
       'O.B3' = array(rnbinom(len, mu = obs.mu$Temp.Infl * expo.obs, 
                              size = exp.ds$Temp.Infl), dim = dms),
       'S.B1' = array(rnbinom(len, mu = sim.mu$Temp * expo.obs, 
                              size = exp.ds$Temp), dim = dms),
       'S.B2' = array(rnbinom(len, mu = sim.mu$Infl * expo.obs, 
                              size = exp.ds$Infl), dim = dms),
       'S.B3' = array(rnbinom(len, mu = sim.mu$Temp.Infl * expo.obs, 
                              size = exp.ds$Temp.Infl), dim = dms)
  )
})
stopCluster(cl)

# Estimated force of mortality - observed covariates - NB noise
mhat.obs1 <- simplify2array(lapply(results, function(x) x$O.B1))
mhat.obs2 <- simplify2array(lapply(results, function(x) x$O.B2))
mhat.obs3 <- simplify2array(lapply(results, function(x) x$O.B3))

# Estimated force of mortality - observed covariates - no NB noise
mhat.obsA1 <- simplify2array(lapply(results, function(x) x$O.A1))
mhat.obsA2 <- simplify2array(lapply(results, function(x) x$O.A2))
mhat.obsA3 <- simplify2array(lapply(results, function(x) x$O.A3))

# Estimated force of mortality - simulated covariates - NB noise
mhat.sim1 <- simplify2array(lapply(results, function(x) x$S.B1))
mhat.sim2 <- simplify2array(lapply(results, function(x) x$S.B2))
mhat.sim3 <- simplify2array(lapply(results, function(x) x$S.B3))

# Estimated force of mortality - simulated covariates - no NB noise
mhat.simA1 <- simplify2array(lapply(results, function(x) x$S.A1))
mhat.simA2 <- simplify2array(lapply(results, function(x) x$S.A2))
mhat.simA3 <- simplify2array(lapply(results, function(x) x$S.A3))

# Compute quantiles - observed covariates - NB noise
quant.mhat.obs1 <- apply(mhat.obs1, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.obs2 <- apply(mhat.obs2, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.obs3 <- apply(mhat.obs3, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()

# Compute quantiles - observed covariates - no NB noise
quant.mhat.obsA1 <- apply(mhat.obsA1, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.obsA2 <- apply(mhat.obsA2, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.obsA3 <- apply(mhat.obsA3, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()

# Compute quantiles - simulated covariates - NB noise
quant.mhat.sim1 <- apply(mhat.sim1, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.sim2 <- apply(mhat.sim2, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.sim3 <- apply(mhat.sim3, MARGIN = c(1, 2, 3), 
                         FUN = function(x) 
                           quantile(x, probs = c(0.025, 0.5, 0.975))); gc()

# Compute quantiles - simulated covariates - no NB noise
quant.mhat.simA1 <- apply(mhat.simA1, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.simA2 <- apply(mhat.simA2, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()
quant.mhat.simA3 <- apply(mhat.simA3, MARGIN = c(1, 2, 3), 
                          FUN = function(x) 
                            quantile(x, probs = c(0.025, 0.5, 0.975))); gc()

# Set dim names
dimnames(mhat.obs1) = dimnames(mhat.obs2) = dimnames(mhat.obs3) = 
  dimnames(mhat.sim1) = dimnames(mhat.sim2) = dimnames(mhat.sim3) = 
  dimnames(mhat.obsA1) = dimnames(mhat.obsA2) = dimnames(mhat.obsA3) = 
  dimnames(mhat.simA1) = dimnames(mhat.simA2) = dimnames(mhat.simA3) <- 
  list(as.character(unique(Dfnew$ISODate)), unique(Dfnew$age), 
       unique(Dfnew$region), 1:n.sim)

dimnames(quant.mhat.obs1) = dimnames(quant.mhat.obs2) = 
  dimnames(quant.mhat.obs3) =  dimnames(quant.mhat.sim1) = 
  dimnames(quant.mhat.sim2) = dimnames(quant.mhat.sim3) = 
  dimnames(quant.mhat.obsA1) = dimnames(quant.mhat.obsA2) = 
  dimnames(quant.mhat.obsA3) =  dimnames(quant.mhat.simA1) = 
  dimnames(quant.mhat.simA2) = dimnames(quant.mhat.simA3) <- 
  list(rownames(quant.mhat.sim1), as.character(unique(Dfnew$ISODate)),
       unique(Dfnew$age), unique(Dfnew$region))

# Plot
df.mhat.obs1  <- as.data.frame.table(quant.mhat.obs1, responseName = "value") 
df.mhat.obs2  <- as.data.frame.table(quant.mhat.obs2, responseName = "value") 
df.mhat.obs3  <- as.data.frame.table(quant.mhat.obs3, responseName = "value") 

df.mhat.obsA1 <- as.data.frame.table(quant.mhat.obsA1, responseName = "value") 
df.mhat.obsA2 <- as.data.frame.table(quant.mhat.obsA2, responseName = "value") 
df.mhat.obsA3 <- as.data.frame.table(quant.mhat.obsA3, responseName = "value") 

df.mhat.sim1  <- as.data.frame.table(quant.mhat.sim1, responseName = "value") 
df.mhat.sim2  <- as.data.frame.table(quant.mhat.sim2, responseName = "value") 
df.mhat.sim3  <- as.data.frame.table(quant.mhat.sim3, responseName = "value") 

df.mhat.simA1 <- as.data.frame.table(quant.mhat.simA1, responseName = "value") 
df.mhat.simA2 <- as.data.frame.table(quant.mhat.simA2, responseName = "value") 
df.mhat.simA3 <- as.data.frame.table(quant.mhat.simA3, responseName = "value") 

df.deaths <- as.data.frame.table(dtxr_ft)
df.expo   <- as.data.frame.table(etxr_ft)
colnames(df.mhat.obs1) = colnames(df.mhat.obs2) = colnames(df.mhat.obs3) = 
  colnames(df.mhat.obsA1) = colnames(df.mhat.obsA2) = colnames(df.mhat.obsA3) = 
  colnames(df.mhat.sim1) = colnames(df.mhat.sim2) = colnames(df.mhat.sim3) = 
  colnames(df.mhat.simA1) = colnames(df.mhat.simA2) = colnames(df.mhat.simA3) <- 
  c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo) <- c('ISODate', 'age', 'region', 'expo')

df.tot.mhat.obs <- lapply(list(df.mhat.obs1, df.mhat.obs2, df.mhat.obs3,
                               df.mhat.obsA1, df.mhat.obsA2, df.mhat.obsA3),
                          function(x) plyr::join_all(list(x, df.deaths, df.expo),
                                                     by = c('ISODate', 'age', 'region')) %>%
                            tidyr::pivot_wider(names_from = Quantile,
                                               values_from = sim,
                                               names_prefix = "q"))
names(df.tot.mhat.obs) <- c('Temp', 'Infl', 'Temp.Infl', 
                            'TempA', 'InflA', 'Temp.InflA')

df.tot.mhat.sim <- lapply(list(df.mhat.sim1, df.mhat.sim2, df.mhat.sim3,
                               df.mhat.simA1, df.mhat.simA2, df.mhat.simA3),
                          function(x) plyr::join_all(list(x, df.deaths, df.expo),
                                                     by = c('ISODate', 'age', 'region')) %>%
                            tidyr::pivot_wider(names_from = Quantile,
                                               values_from = sim,
                                               names_prefix = "q"))
names(df.tot.mhat.sim) <- c('Temp', 'Infl', 'Temp.Infl', 
                            'TempA', 'InflA', 'Temp.InflA')


age_vec <- c('50', '70', '90')
reg     <- regions[c(5)]
dfs <- df.tot.mhat.sim$Temp.InflA %>% 
  dplyr::filter(age %in% age_vec, region == reg)
DfM$age <- as.factor(DfM$age)
dfs$age <- as.factor(dfs$age)

ggplot(dfs) + 
  theme_bw(base_size = 15) + 
  geom_line(data = DfM %>% dplyr::filter(age %in% age_vec, region == reg),
            aes(x = as.Date(ISODate), y = log(deaths), group = age, col = age))  + 
  geom_point(aes(x = as.Date(ISODate), y = log(obs), group = age, col = age),
             shape = 3, col = 'black', cex = 0.2) + 
  geom_ribbon(aes(x = as.Date(ISODate), ymin = log(`q2.5%`), ymax = log(`q97.5%`), 
                  group = age, fill = age, col = age), alpha = 0.2) +
  geom_line(aes(x = as.Date(ISODate), y = log(`q50%`), group = age, col = age)) +
  coord_cartesian(xlim= c(as.Date('2010-01-01'), as.Date('2019-12-31')),
                  ylim = c(pmax(min(log(dfs$obs)),0), max(log(dfs$obs)))) + 
  scale_fill_manual(name = 'Age', 
                    values = c("#33A02C", "#1F78B4", "#E31A1C"),
                    labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
  scale_color_manual(name = 'Age', 
                     values = c("#33A02C", "#1F78B4", "#E31A1C"),
                     labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) 

##### 8) Scoring rules ----
grid   <- expand.grid(age = unique(DfM$age), region = unique(DfM$region))
crps_obs_all <- lapply(list(mhat.obs1, mhat.obs2, mhat.obs3), function(m) 
  sapply(1:nrow(grid), function(i){
    x <- as.character(grid[i,1])
    r <- as.character(grid[i,2])
    mean(crps_sample(y = dt[,x,r], dat = m[,x,r,]))
    }))
crps_sim_all <- lapply(list(mhat.sim1, mhat.sim2, mhat.sim3), function(m) 
  sapply(1:nrow(grid), function(i){
    x <- as.character(grid[i,1])
    r <- as.character(grid[i,2])
    mean(crps_sample(y = dt[,x,r], dat = m[,x,r,]))
  }))
logs_obs_all <- lapply(list(mhat.obs1, mhat.obs2, mhat.obs3), function(m) 
  sapply(1:nrow(grid), function(i){
    x <- as.character(grid[i,1])
    r <- as.character(grid[i,2])
    mean(logs_sample(y = dt[,x,r], dat = m[,x,r,]))
    }))
logs_sim_all <- lapply(list(mhat.sim1, mhat.sim2, mhat.sim3), function(m) 
  sapply(1:nrow(grid), function(i){
    x <- as.character(grid[i,1])
    r <- as.character(grid[i,2])
    mean(logs_sample(y = dt[,x,r], dat = m[,x,r,]))
  }))
grid[,'crps.obs.Temp']      <- crps_obs_all[[1]]
grid[,'crps.obs.Infl']      <- crps_obs_all[[2]]
grid[,'crps.obs.Temp.Infl'] <- crps_obs_all[[3]]
grid[,'crps.sim.Temp']      <- crps_sim_all[[1]]
grid[,'crps.sim.Infl']      <- crps_sim_all[[2]]
grid[,'crps.sim.Temp.Infl'] <- crps_sim_all[[3]]

grid[,'logs.obs.Temp']      <- logs_obs_all[[1]]
grid[,'logs.obs.Infl']      <- logs_obs_all[[2]]
grid[,'logs.obs.Temp.Infl'] <- logs_obs_all[[3]]
grid[,'logs.sim.Temp']      <- logs_sim_all[[1]]
grid[,'logs.sim.Infl']      <- logs_sim_all[[2]]
grid[,'logs.sim.Temp.Infl'] <- logs_sim_all[[3]]

crps_obs <- sapply(c('crps.obs.Temp', 'crps.obs.Infl', 'crps.obs.Temp.Infl'),
                   function(c) mean(grid[,c]))
crps_sim <- sapply(c('crps.sim.Temp', 'crps.sim.Infl', 'crps.sim.Temp.Infl'),
                   function(c) mean(grid[,c]))
logs_obs <- sapply(c('logs.obs.Temp', 'logs.obs.Infl', 'logs.obs.Temp.Infl'),
                   function(c) mean(grid[,c]))
logs_sim <- sapply(c('logs.sim.Temp', 'logs.sim.Infl', 'logs.sim.Temp.Infl'),
                   function(c) mean(grid[,c]))

crps_age_obs <- sapply(c('crps.obs.Temp', 'crps.obs.Infl', 'crps.obs.Temp.Infl'),
                       function(c) grid %>% dplyr::group_by(age) %>% 
                         dplyr::reframe(m = mean(!!sym(c))) %>% pull(m))
crps_age_sim <- sapply(c('crps.sim.Temp', 'crps.sim.Infl', 'crps.sim.Temp.Infl'),
                      function(c) grid %>% dplyr::group_by(age) %>% 
                        dplyr::reframe(m = mean(!!sym(c))) %>% pull(m))
logs_age_obs <- sapply(c('logs.obs.Temp', 'logs.obs.Infl', 'logs.obs.Temp.Infl'),
                       function(c) grid %>% dplyr::group_by(age) %>% 
                         dplyr::reframe(m = mean(!!sym(c))) %>% pull(m))
logs_age_sim <- sapply(c('logs.sim.Temp', 'logs.sim.Infl', 'logs.sim.Temp.Infl'),
                      function(c) grid %>% dplyr::group_by(age) %>% 
                        dplyr::reframe(m = mean(!!sym(c))) %>% pull(m))

# saveRDS(crps_age_obs, '../Results updated/CRPS_age_obs_Models4-6.rds')
# saveRDS(crps_age_sim, '../Results updated/CRPS_age_sim_Models4-6.rds')
# saveRDS(crps_obs, '../Results updated/CRPS_overall_obs_Models4-6.rds')
# saveRDS(crps_sim, '../Results updated/CRPS_overall_sim_Models4-6.rds')
# saveRDS(logs_age_obs, '../Results updated/LogS_age_obs_Models4-6.rds')
# saveRDS(logs_age_sim, '../Results updated/LogS_age_sim_Models4-6.rds')
# saveRDS(logs_obs, '../Results updated/LogS_overall_obs_Models4-6.rds')
# saveRDS(logs_sim, '../Results updated/LogS_overall_sim_Models4-6.rds')

##### 9) Coverage ratio and interval score ----

for(c in c('Temp', 'Infl', 'Temp.Infl')) {
  df.tot.mhat.obs[[c]]$coverage <- df.tot.mhat.obs[[c]]$obs < df.tot.mhat.obs[[c]]$`q97.5%` & 
    df.tot.mhat.obs[[c]]$obs > df.tot.mhat.obs[[c]]$`q2.5%`
  df.tot.mhat.obs[[c]]$iscore <- with(
    df.tot.mhat.obs[[c]], (`q97.5%` - `q2.5%`) + 
      (2 / 0.05) * (`q2.5%` - obs) * (obs < `q2.5%`) +
      (2 / 0.05) * (obs - `q97.5%`) * (obs > `q97.5%`))
  
  df.tot.mhat.sim[[c]]$coverage <- df.tot.mhat.sim[[c]]$obs < df.tot.mhat.sim[[c]]$`q97.5%` & 
    df.tot.mhat.sim[[c]]$obs > df.tot.mhat.sim[[c]]$`q2.5%`
  df.tot.mhat.sim[[c]]$iscore <- with(
    df.tot.mhat.sim[[c]], (`q97.5%` - `q2.5%`) + 
      (2 / 0.05) * (`q2.5%` - obs) * (obs < `q2.5%`) +
      (2 / 0.05) * (obs - `q97.5%`) * (obs > `q97.5%`))
}

cov_obs    <- sapply(c('Temp', 'Infl', 'Temp.Infl'), 
                  function(c) mean(df.tot.mhat.obs[[c]]$coverage))
cov_sim    <- sapply(c('Temp', 'Infl', 'Temp.Infl'), 
                  function(c) mean(df.tot.mhat.sim[[c]]$coverage))
iscore_obs <- sapply(c('Temp', 'Infl', 'Temp.Infl'), 
                  function(c) mean(df.tot.mhat.obs[[c]]$iscore))
iscore_sim <- sapply(c('Temp', 'Infl', 'Temp.Infl'), 
                  function(c) mean(df.tot.mhat.sim[[c]]$iscore))

cov_age_obs    <- sapply(c('Temp', 'Infl', 'Temp.Infl'), function(c) 
  df.tot.mhat.obs[[c]] %>% dplyr::group_by(age) %>% 
    dplyr::reframe(m = mean(coverage)) %>% pull(m))
cov_age_sim    <- sapply(c('Temp', 'Infl', 'Temp.Infl'), function(c) 
  df.tot.mhat.sim[[c]] %>% dplyr::group_by(age) %>% 
    dplyr::reframe(m = mean(coverage)) %>% pull(m))
iscore_age_obs <- sapply(c('Temp', 'Infl', 'Temp.Infl'), function(c) 
  df.tot.mhat.obs[[c]] %>% dplyr::group_by(age) %>% 
    dplyr::reframe(m = mean(iscore)) %>% pull(m))
iscore_age_sim <- sapply(c('Temp', 'Infl', 'Temp.Infl'), function(c) 
  df.tot.mhat.sim[[c]] %>% dplyr::group_by(age) %>% 
    dplyr::reframe(m = mean(iscore)) %>% pull(m))

# saveRDS(cov_age_obs, '../Results updated/COV_age_obs_Models4-6.rds')
# saveRDS(cov_age_sim, '../Results updated/COV_age_sim_Models4-6.rds')
# saveRDS(cov_obs, '../Results updated/COV_overall_obs_Models4-6.rds')
# saveRDS(cov_sim, '../Results updated/COV_overall_sim_Models4-6.rds')
# saveRDS(iscore_age_obs, '../Results updated/IScore_age_obs_Models4-6.rds')
# saveRDS(iscore_age_sim, '../Results updated/IScore_age_sim_Models4-6.rds')
# saveRDS(iscore_obs, '../Results updated/IScore_overall_obs_Models4-6.rds')
# saveRDS(iscore_sim, '../Results updated/IScore_overall_sim_Models4-6.rds')

##### 10) Uncertainty quantification -----

dt.sim.TempInfl <- df.tot.mhat.sim$Temp.InflA %>%
  dplyr::left_join(df.tot.mhat.sim$Temp.Infl, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.obs.TempInfl <- df.tot.mhat.obs$Temp.InflA %>%
  dplyr::left_join(df.tot.mhat.obs$Temp.Infl, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.sim.Temp     <- df.tot.mhat.sim$TempA %>%
  dplyr::left_join(df.tot.mhat.sim$Temp, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.obs.Temp     <- df.tot.mhat.obs$TempA %>%
  dplyr::left_join(df.tot.mhat.obs$Temp, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.sim.Infl     <- df.tot.mhat.sim$InflA %>%
  dplyr::left_join(df.tot.mhat.sim$Infl, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.obs.Infl     <- df.tot.mhat.obs$InflA %>%
  dplyr::left_join(df.tot.mhat.obs$Infl, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))

sum.sim.TempInfl <- dt.sim.TempInfl %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum.obs.TempInfl <- dt.obs.TempInfl %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 

sum.sim.Temp <- dt.sim.Temp %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum.obs.Temp <- dt.obs.Temp %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 

sum.sim.Infl <- dt.sim.Infl %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum.obs.Infl <- dt.obs.Infl %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 


# saveRDS(sum.sim.Temp, '../Results updated/UDecom_sim_Model4.rds')
# saveRDS(sum.obs.Temp, '../Results updated/UDecom_obs_Model4.rds')
# saveRDS(sum.sim.Infl, '../Results updated/UDecom_sim_Model5.rds')
# saveRDS(sum.obs.Infl, '../Results updated/UDecom_obs_Model5.rds')
# saveRDS(sum.sim.TempInfl, '../Results updated/UDecom_sim_Model6.rds')
# saveRDS(sum.obs.TempInfl, '../Results updated/UDecom_obs_Model6.rds')

##### 11) Figure ----

pvis <- list()
reg_vec <- c('32', '53', '84')
age_vec <- c('50', '70', '90')

for(r in 1:length(reg_vec)){
  dfs <- df.tot.mhat.sim$Temp.Infl %>% 
    dplyr::filter(region == reg_vec[r],
                  age %in% age_vec)
  
  pvis[[r]] <- ggplot(dfs) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Date','')) + 
    ylab(if(r==1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(data = DfM %>% dplyr::filter(region == reg_vec[r], age %in% c('50', '70', '90')),
              aes(x = as.Date(ISODate), y = log(deaths/expo), group = age, col = age))  + 
    geom_point(aes(x = as.Date(ISODate), y = log(obs/expo), group = age, col = age),
               shape = 3, col = 'black', cex = 0.15) + 
    geom_ribbon(aes(x = as.Date(ISODate), ymin = log(`q2.5%`/expo), ymax = log(`q97.5%`/expo), 
                    group = age, fill = age, col = age), alpha = 0.2) +
    geom_line(aes(x = as.Date(ISODate), y = log(`q50%`/expo), group = age, col = age)) +
    coord_cartesian(xlim= c(as.Date('2012-01-01'), as.Date('2020-01-01')),
                    ylim = c(-7.5,-1)) + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) 
}

plot.forecast <- ggpubr::ggarrange(plotlist = pvis, common.legend = TRUE, 
                                   legend = 'bottom', nrow = 1, ncol = 3, 
                                   align = 'h') 

# ggsave('../Figures/predperform.pdf',
#        plot = plot.forecast, height = 12, width = 30, units = 'cm')
