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

# Lee-Carter model
# lca1 <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                         etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                         cb = list('Temp' = cb1, 'Infl' = cb2),
#                         lambda = c('Temp' = 0, 'Infl' = 0))
# saveRDS(lca1, file = '../Results updated/lca_Temp_Infl_NoPenal_benchmark.rds')
# lca1 <- readRDS('../Results updated/lca_Temp_Infl_NoPenal_benchmark.rds')

lc0 <- readRDS('../Results updated/lca_Temp_Infl_benchmark.rds')
lc1 <- readRDS('../Results updated/lca_Temp_Infl_NoPenal_benchmark.rds')

lc0$lcd$ll
lc1$lcd$ll

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

# Fit ARIMA for period effect
obj.kt <- lapply(1:length(regions), function(r){
  list('Temp.Infl0' = forecast::Arima(ts(lc0$lcd$kappa2[,r]), 
                                      order = c(0,1,1),
                                      include.drift = TRUE),
       'Temp.Infl1' = forecast::Arima(ts(lc1$lcd$kappa2[,r]), 
                                      order = c(0,1,1),
                                      include.drift = TRUE))
})

# Point forecasts
tg_fc <- tg.obs
ia_fc <- ia.obs
kt_fc <- list(
  'Temp.Infl0' = sapply(1:dim(dtxr)[3], function(r) {
    c(forecast::forecast(obj.kt[[r]]$Temp.Infl0, nY)$mean)}),
  'Temp.Infl1' = sapply(1:dim(dtxr)[3], function(r) {
    c(forecast::forecast(obj.kt[[r]]$Temp.Infl1, nY)$mean)})
)

# Cross-bases
cb1_fc <- cb_forecast(tg_fc, hist.tg, 'Temp', mxlag.tavg, means1)
cb2_fc <- cb_forecast(ia_fc, hist.ia, 'Infl', mxlag.infl, means2)

# DLNM
dlnm <- list(
  'Temp.Infl0' = Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out}, list('Temp' = cb1_fc, 'Infl' = cb2_fc), lc0$lcd$gammal),
  'Temp.Infl1' = Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out}, list('Temp' = cb1_fc, 'Infl' = cb2_fc), lc1$lcd$gammal)
)

# Force of mortality
mhat0 <- exp(sweep(aperm(
  lc0$lcd$beta2 %o% kt_fc$Temp.Infl0[isoyv_fc - min(isoyv_fc) + 1,] +
    lc0$lcd$beta3 %o% lc0$lcd$kappa3[isowv_fc,] +
    Reduce(`+`, Map(function(b, d) b %o% d, lc0$lcd$betal, dlnm$Temp.Infl0)),
  c(2,1,3)), c(2,3), lc0$lcd$beta1, '+'))
mhat1 <- exp(sweep(aperm(
  lc1$lcd$beta2 %o% kt_fc$Temp.Infl1[isoyv_fc - min(isoyv_fc) + 1,] +
    lc1$lcd$beta3 %o% lc1$lcd$kappa3[isowv_fc,] +
    Reduce(`+`, Map(function(b, d) b %o% d, lc1$lcd$betal, dlnm$Temp.Infl1)),
  c(2,1,3)), c(2,3), lc1$lcd$beta1, '+'))

dimnames(mhat0) = dimnames(mhat1) <- 
  list(c(as.character(sort(unique(Dfnew$ISODate)))),
       dimnames(lc1$lcd$mhat)[[2]],
       dimnames(lc1$lcd$mhat)[[3]])

# Error metrics
dt <- dtxr_ft
et <- etxr_ft
mhat_rep0 <- mhat0
mhat_rep1 <- mhat1
diff0 <- dt - mhat0 * et
diff1 <- dt - mhat1 * et

plot(as.Date(rownames(dt)), (dt/et)[,9,1], type = 'l')
lines(as.Date(rownames(mhat0)), mhat0[,9,1], col = 'orange', lwd = 2)
lines(as.Date(rownames(mhat1)), mhat1[,9,1], col = 'purple4', lwd = 2)

rmse <- c('Temp.Infl0' = sqrt(1/length(diff0) * sum(((diff0)^2))),
          'Temp.Infl1' = sqrt(1/length(diff1) * sum(((diff1)^2))))
mae  <- c('Temp.Infl0' = 1/length(diff0) * sum(abs(diff0)),
          'Temp.Infl1' = 1/length(diff1) * sum(abs(diff1)))

## Table C.1
by_reg <- function(d, f) sapply(seq_len(dim(d)[3]), function(r) f(d[,,r]))

tabC1 <- data.frame(
  Region     = shapef$NUTS_NAME[match(regions, shapef$region)],
  RMSE_unpen = by_reg(diff1, function(z) sqrt(mean(z^2))),
  RMSE_pen   = by_reg(diff0, function(z) sqrt(mean(z^2))),
  MAE_unpen  = by_reg(diff1, function(z) mean(abs(z))),
  MAE_pen    = by_reg(diff0, function(z) mean(abs(z))))

apply(tabC1[,c(2,3)], 1, which.min)
apply(tabC1[,c(4,5)], 1, which.min)

print(tabC1, digits = 4)

print(xtable::xtable(tabC1, digits = c(0, 0, 2, 2, 2, 2),
                     align = 'llSSSS'),
      include.rownames = FALSE, booktabs = TRUE,
      sanitize.colnames.function = function(x)
        c('Region', 'RMSE: Unpenalized', 'RMSE: Penalized',
          'MAE: Unpenalized', 'MAE: Penalized'))

