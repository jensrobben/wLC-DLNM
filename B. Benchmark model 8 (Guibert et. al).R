##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'spdep',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'data.table', 'xgboost',
              'parallel', 'doParallel', 'readxl', 'tidyr', 'purrr', 'imputeTS',
              'dlnm','Matrix', 'metR', 'pbapply', 'scoringRules', 'copula')
invisible(sapply(packages, require, character.only = TRUE))
`%notin%` <- Negate(`%in%`)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Colour scales
colors <- c('11' = "#A6CEE3", '84' = "#1F78B4", '75' = "#B2DF8A", '76' = "#33A02C", 
            '32' = "#FB9A99", '44' = "#E31A1C", '93' = "#FDBF6F", '53' = "#FF7F00", 
            '28' = "#CAB2D6", '52' = "#6A3D9A", '27' = "#FFFF99", '24' = "#B15928")

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

# Death counts and exposures - matrix format
dtxr     <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM)
etxr     <- stats::xtabs(expo ~ ISODate + age + region, data = DfM)

##### 2) New data ----
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

dtxr_ft <- deaths.obs
etxr_ft <- expo.obs


##### 3) Model 8: two-stage DLNM + Li-Lee of Guibert et al. (2026) ----

### 3.1) Adapt data ----

# Source
source('F. Annual, regional Li-Lee (NB).R')
source('F. Guibert et al. approach.R')

# Daily temperature (NUTS 1), mapped to the INSEE region codes
nuts.tab <- data.frame(Region = c('FR1','FRB','FRC','FRD','FRE','FRF',
                                  'FRG','FRH','FRI','FRJ','FRK','FRL'),
                       region = c('11','24','27','28','32','44',
                                  '52','53','75','76','84','93'))

df.daily <- fread('../Temperature/FR_NUTS1_tg_daily.txt') %>%
  dplyr::left_join(nuts.tab, by = 'Region') %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(ISOYear = isoyear(Date), ISOWeek = isoweek(Date)) %>%
  dplyr::arrange(region, Date) %>%
  as.data.table()

# Weeks of the calibration and forecast period
weeks.cal <- sort(unique(DfM$ISODate))
weeks.fc  <- sort(unique(Dfnew$ISODate))

# Year index of the forecast period
isoyv_fc <- isoyear(as.Date(weeks.fc)) - min(DfM$ISOYear) + 1

# Cross-basis matrices and reference-day temperatures per region
cb.cal <- tg.cal <- tg.fc <- wk.cal <- 
  wk.fc <- info.lst <- vector('list', length(regions))

for(j in seq_along(regions)){
  dj  <- df.daily[region == regions[j]]
  wk  <- ISOweek::ISOweek2date(paste0(sprintf("%04d-W%02d", dj$ISOYear, dj$ISOWeek), "-1"))
  ref <- which(!duplicated(wk, fromLast = TRUE))
  idx.cal <- ref[wk[ref] %in% weeks.cal]
  idx.fc  <- ref[wk[ref] %in% weeks.fc]
  
  cal.days      <- wk %in% weeks.cal  
  fc.days       <- wk %in% weeks.fc
  
  info.lst[[j]] <- basis_info(dj$tg[cal.days])
  
  cb.cal[[j]]   <- cb_guibert(tg.daily = dj$tg, 
                              ref.idx = idx.cal, 
                              info = info.lst[[j]])
  
  tg.cal[[j]]   <- dj$tg[cal.days]
  wk.cal[[j]]   <- wk[cal.days]
  tg.fc[[j]]    <- dj$tg[fc.days]
  wk.fc[[j]]    <- wk[fc.days]
}

### 3.2) Stage 1: DLNM and temperature adjustment ----

# Age strata to stabilize the MMTs -> similar to Guibert approach
strata_M8 <- cut(xv, breaks = c(50, 65, 75, 85, Inf), right = FALSE)

# Stage 1 estimation
st1_M8 <- stage1_guibert(dtx = dtxr, 
                         cb.lst = cb.cal, 
                         tg.daily.lst = tg.cal,
                         strata = strata_M8,
                         week.lst = wk.cal, 
                         info.lst = info.lst, 
                         time.idx = seq_along(weeks.cal),
                         isoyv = isoyv, 
                         df.year = 8,
                         n.year = length(unique(DfM$ISOYear)), 
                         verbose = TRUE)

# Weekly adjustment over the forecast period
Tadj_fc_M8 <- Tadj_guibert(tg.fc, wk.fc, st1_M8, N = length(xv), R = length(regions))

### 3.3) Stage 2: Li-Lee on temperature-adjusted exposures ----

lca_M8 <- aLL.NB(xv = xv, isoyv = isoyv, isowv = isowv,
                 etx = etxr * st1_M8$Tadj.y[isoyv, , ], 
                 dtx = dtxr, eps = 10^(-9))

### 3.4) Point forecast ----
nY       <- 5
n.hor    <- 52*nY + 1

# Fit time series
obj.K_M8  <- forecast::Arima(ts(lca_M8$K), order = c(0,1,0), include.drift = TRUE)
obj.kt_M8 <- lapply(seq_along(regions), function(r)
  forecast::Arima(ts(lca_M8$kappa2[, r]), order = c(1,0,0), include.mean = TRUE))

# Point forecast K and kappa
K_fc_M8  <- c(lca_M8$K, c(forecast::forecast(obj.K_M8, nY)$mean))
kt_fc_M8 <- rbind(lca_M8$kappa2, sapply(seq_along(regions), function(r)
  c(forecast::forecast(obj.kt_M8[[r]], nY)$mean)))

# Point forecast death rates: not attributable to temperature, times the adjustment
mhat_M8 <- array(NA, dim = c(length(weeks.fc), length(xv), length(regions)),
                 dimnames = list(as.character(weeks.fc), as.character(xv), regions))

for(r in seq_along(regions))
  mhat_M8[, , r] <- Tadj_fc_M8[, , r] * 
  exp(outer(K_fc_M8[isoyv_fc], lca_M8$B) +
        outer(kt_fc_M8[isoyv_fc, r], lca_M8$beta2[, r]) +
        rep(lca_M8$beta1[, r], each = length(weeks.fc))) 


# Error metrics
dt <- dtxr_ft
et <- etxr_ft

diff_M8 <- mhat_M8 * et - dt

rmse_M8 <- sqrt(1/length(diff_M8) * sum(((diff_M8)^2)))
mae_M8  <- 1/length(diff_M8) * sum(abs(diff_M8))

rmse_age_M8 <- sapply(1:dim(diff_M8)[2], function(a)
  sqrt(1/length(diff_M8[,a,]) * sum(((diff_M8[,a,])^2))))
mae_age_M8  <- sapply(1:dim(diff_M8)[2], function(a)
  1/length(diff_M8[,a,]) * sum(abs(diff_M8[,a,])))

### 3.5) Scenarios ----

n.sim <- 5000

cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl, c('n.sim', 'nY', 'regions', 'xv', 'isoyv_fc', 'etxr_ft', 
                    'weeks.fc', 'lca_M8', 'st1_M8', 'obj.K_M8', 'obj.kt_M8', 
                    'tg.fc', 'wk.fc', 'cum_logRR'))
clusterEvalQ(cl, expr = {
  library(dplyr)
  library(forecast)
  library(dlnm)
  library(MASS)
})
mlst.sim_M8 <- pblapply(cl = cl, X = 1:n.sim, FUN = function(s) {
  
  # Simulate K-process (common) and kt-processes (regional)
  sim.K  <- c(lca_M8$K, simulate(obj.K_M8, nsim = nY))
  sim.kt <- rbind(lca_M8$kappa2, sapply(seq_along(regions), function(r)
    simulate(obj.kt_M8[[r]], nsim = nY)))
  
  # Force of mortality with bootstrapped DLNM coefficients
  mu <- array(NA, dim = c(length(weeks.fc), length(xv), length(regions)))
  
  for(r in seq_along(regions)) for(x in seq_along(xv)){
    th <- MASS::mvrnorm(1, st1_M8$theta[[(r - 1) * length(xv) + x]],
                        st1_M8$vcov[[(r - 1) * length(xv) + x]])
    AF.d <- 1 - exp(-(cum_logRR(tg.fc[[r]], th, st1_M8$info[[r]]) -
                        cum_logRR(st1_M8$mmt[r], th, st1_M8$info[[r]])))
    Tb   <- (1 - as.vector(tapply(AF.d, wk.fc[[r]], mean)))^(-1)
    mu[, x, r] <- exp(lca_M8$beta1[x, r] + lca_M8$B[x] * sim.K[isoyv_fc] +
                        lca_M8$beta2[x, r] * sim.kt[isoyv_fc, r]) * Tb
  }
  
  # Dispersion
  exp.ds <- exp(outer(lca_M8$phix, lca_M8$phir, "+"))
  exp.ds <- aperm(replicate(length(weeks.fc), exp.ds, simplify = "array"), c(3, 1, 2))
  
  list('A' = mu * etxr_ft,
       'B' = array(rnbinom(length(mu), mu = mu * etxr_ft, size = exp.ds), dim = dim(mu)))
})
stopCluster(cl)

mhat.sim_M8 <- simplify2array(lapply(mlst.sim_M8, function(x) x$B))
dimnames(mhat.sim_M8) <- list(as.character(weeks.fc), as.character(xv), regions, NULL)

# Quantiles, coverage ratio and interval score
quant_M8 <- apply(mhat.sim_M8, c(1, 2, 3), quantile, probs = c(0.025, 0.5, 0.975))
lo_M8    <- quant_M8[1,,,]
hi_M8    <- quant_M8[3,,,]

cov_M8    <- mean(dt < hi_M8 & dt > lo_M8)
iscore_M8 <- mean((hi_M8 - lo_M8) +
                    (2/0.05) * (lo_M8 - dt) * (dt < lo_M8) +
                    (2/0.05) * (dt - hi_M8) * (dt > hi_M8))

cov_age_M8    <- sapply(1:dim(dt)[2], function(a)
  mean(dt[,a,] < hi_M8[,a,] & dt[,a,] > lo_M8[,a,]))
iscore_age_M8 <- sapply(1:dim(dt)[2], function(a)
  mean((hi_M8[,a,] - lo_M8[,a,]) +
         (2/0.05) * (lo_M8[,a,] - dt[,a,]) * (dt[,a,] < lo_M8[,a,]) +
         (2/0.05) * (dt[,a,] - hi_M8[,a,]) * (dt[,a,] > hi_M8[,a,])))

# Scoring rules
grid <- expand.grid(age = as.character(xv), region = regions, stringsAsFactors = FALSE)
grid[, 'crps'] <- sapply(1:nrow(grid), function(i)
  mean(crps_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M8[, grid$age[i], grid$region[i], ])))
grid[, 'logs'] <- sapply(1:nrow(grid), function(i)
  mean(logs_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M8[, grid$age[i], grid$region[i], ])))

crps_M8 <- mean(grid$crps)
logs_M8 <- mean(grid$logs)

crps_age_M8 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_age_M8 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)
