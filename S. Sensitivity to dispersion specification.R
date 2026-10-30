##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'spdep',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'data.table', 'xgboost',
              'parallel', 'doParallel', 'readxl', 'tidyr', 'purrr', 'imputeTS',
              'dlnm','Matrix', 'metR', 'pbapply', 'snow', 'doSNOW', 'cli',
              'abind', 'collapse')
invisible(sapply(packages, require, character.only = TRUE))
`%notin%` <- Negate(`%in%`)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Colour scales
colors <- c('11' = "#A6CEE3", '84' = "#1F78B4", '75' = "#B2DF8A", '76' = "#33A02C", 
            '32' = "#FB9A99", '44' = "#E31A1C", '93' = "#FDBF6F", '53' = "#FF7F00", 
            '28' = "#CAB2D6", '52' = "#6A3D9A", '27' = "#FFFF99", '24' = "#B15928")

# Source
source('F. Weekly Lee-Carter + time disp (NB).R')
source('F. Weekly DLNM + time disp (NB).R')
source('F. Weekly Lee-Carter + DLNM + time disp (NB).R')


##### 1) Fit age-enhanced DLNM ----
### 1.1) Make everything ready to fit ----
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
                ISOYear <= 2019,
                ISOYear >= 1990) %>% 
  dplyr::arrange(age, ISODate)

# Death counts and exposures - matrix format
dtxr  <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM)
etxr  <- stats::xtabs(expo ~ ISODate + age + region, data = DfM)

# Avoid numerical issues
dtxr[dtxr == 0] <- 0.001

# Age range, year range, week range
xv    <- as.integer(colnames(dtxr))
isoyv <- isoyear(as.Date(rownames(dtxr))) - min(DfM$ISOYear) + 1
isowv <- isoweek(as.Date(rownames(dtxr)))
R     <- dim(dtxr)[3]

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

### 1.2) Construct cross-basis matrices ----

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
args <- lapply(dimnames(dtxr)[[3]], function(r) {
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
                      scale = max(infl[,ia100]))
  arglag.infl <- list(fun = "ns",
                      intercept = TRUE,
                      knots = dlnm::logknots(mxlag.infl, nk = 2),
                      Boundary.knots = c(0,mxlag.infl))
  
  list('Temp' = list(argvar = argvar.tavg, arglag = arglag.tavg),
       'Infl' = list(argvar = argvar.infl, arglag = arglag.infl))
})
names(args) <- dimnames(dtxr)[[3]]

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
cb1  <- simplify2array(lapply(dimnames(dtxr)[[3]], function(r) {
  argv <- args[[r]]$Temp$argvar
  argl <- args[[r]]$Temp$arglag
  
  cb_fun(dat = tavg[region == r, tg], lags = 0:mxlag.tavg,
         argvar = argv, arglag = argl)
}))
cb2  <- simplify2array(lapply(dimnames(dtxr)[[3]], function(r) {
  argv <- args[[r]]$Infl$argvar
  argl <- args[[r]]$Infl$arglag
  
  cb_fun(dat = infl[region == r, ia100], lags = 0:mxlag.infl,
         argvar = argv, arglag = argl)
}))
cb12 <- simplify2array(lapply(dimnames(dtxr)[[3]], function(r) {
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

cb1  <- center_cb(cb1, isowv)
cb2  <- center_cb(cb2, isowv)
cb12 <- center_cb(cb12, isowv)

### 1.3) Fit model ----

# Fit model
obj0 <- readRDS('../Results updated/lca_orig.rds')
obj1 <- fit_wLC.DLNM.NB.ext(xv = xv, isoyv = isoyv, isowv = isowv,
                            etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
                            cb = list('Temp' = cb1, 'Infl' = cb2),
                            lambda = c('Temp' = 0, 'Infl' = 0))

c(obj0$lcd$ll, obj1$lcd$ll)

# BIC
npar0 <- length(obj0$lcd$beta1) + length(obj0$lcd$beta2) + length(obj0$lcd$beta3) +
  length(unlist(obj0$lcd$betal)) + length(obj0$lcd$kappa2) +
  length(obj0$lcd$kappa3) + length(unlist(obj0$lcd$gammal)) +
  length(obj0$lcd$phix) + length(obj0$lcd$phir) -
  (2 + 2*R + 1 + length(obj0$lcd$betal))
npar1 <- npar0 + length(obj1$lcd$phiw) - 1

c(-2 * obj0$lcd$ll + log(length(dtxr)) * npar0,
  -2 * obj1$lcd$ll + log(length(dtxr)) * npar1)
