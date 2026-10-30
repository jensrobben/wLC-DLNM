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


##### 3) Model 7: weekly DLNM-Li-Lee of Min et al. (2026) ----

### 3.1) Adapt data ----

# Sourcef
source('F. Min et. al approach.R')

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
weeks.all <- c(weeks.cal, weeks.fc)

# Log death rates (age x week x region), zero counts already replaced by 0.001
mtx_M7 <- stats::xtabs(deaths ~ age + ISODate + region, data = DfM) /
  stats::xtabs(expo ~ age + ISODate + region, data = DfM)
lm_M7  <- log(mtx_M7)

# Cross-basis and wave variables per region, over calibration and forecast period
cb.cal <- cb.all <- wave.cal <- wave.all <- vector('list', length(regions))

for(j in seq_along(regions)){
  dj  <- df.daily[region == regions[j]]
  wk  <- ISOweek::ISOweek2date(paste0(sprintf("%04d-W%02d", dj$ISOYear, dj$ISOWeek), "-1"))
  ref <- which(!duplicated(wk, fromLast = TRUE))                
  idx.cal <- ref[wk[ref] %in% weeks.cal]
  idx.all <- ref[wk[ref] %in% weeks.all]
  
  cal.days      <- wk %in% weeks.cal                     
  info          <- dj$tg[cal.days]                      
  
  cb.cal[[j]]   <- cb_min(dj$tg, idx.cal, ref = info)
  cb.all[[j]]   <- cb_min(dj$tg, idx.all, ref = info)
  
  wv            <- wave_days(dj$tg, wk, ref = info)
  rownames(wv)  <- as.character(unique(wk))
  wave.cal[[j]] <- wv[as.character(weeks.cal), ]
  wave.all[[j]] <- wv[as.character(weeks.all), ]
}

### 3.2) Fit ----
lca_M7 <- fit_wDLNM.LL.min(lm.arr = lm_M7, cb.lst = cb.cal, wave.lst = wave.cal,
                           J = 20, delta = 10^(-2), verbose = TRUE)

### 3.3) Point forecast ----
nY       <- 5
n.hor    <- 52*nY + 1

# Time series models for the common and regional factors (Min et al. use auto.arima
# with a seasonal lag of 52)
obj.K_M7  <- forecast::auto.arima(ts(lca_M7$K, frequency = 52), 
                                  max.p = 1, max.q = 1, max.P = 1, max.Q = 1, 
                                  approximation = TRUE, stepwise = TRUE)
cl <- makeCluster(13)
clusterExport(cl, c('lca_M7', 'regions'), envir = environment())
clusterEvalQ(cl, library(forecast))
obj.kt_M7 <- parLapply(cl, seq_along(regions), function(r)
  forecast::auto.arima(ts(lca_M7$kappa[, r], frequency = 52), 
                          max.p = 1, max.q = 1, max.P = 1, max.Q = 1, 
                          approximation = TRUE, stepwise = TRUE))
stopCluster(cl)

K_fc_M7  <- c(forecast::forecast(obj.K_M7, n.hor)$mean)
kt_fc_M7 <- sapply(seq_along(regions), function(r)
  c(forecast::forecast(obj.kt_M7[[r]], n.hor)$mean))

# Point forecast of the log death rates
mhat_M7 <- array(NA, dim = c(n.hor, length(xv), length(regions)),
                 dimnames = list(as.character(weeks.fc), as.character(xv), regions))
i.fc <- (length(weeks.cal) + 1):length(weeks.all)
for(r in seq_along(regions)){
  X <- cbind(1, cb.all[[r]][i.fc, ], wave.all[[r]]$HWD[i.fc], wave.all[[r]]$CWD[i.fc])
  for(x in seq_along(xv))
    mhat_M7[, x, r] <- exp(lca_M7$A[x, r] + lca_M7$B[x] * K_fc_M7 +
                             lca_M7$b[x, r] * kt_fc_M7[, r] + X %*% lca_M7$coef[, x, r])
}

# Error metrics
dt <- dtxr_ft
et <- etxr_ft

diff_M7 <- mhat_M7 * et - dt

rmse_M7 <- sqrt(1/length(diff_M7) * sum(((diff_M7)^2)))
mae_M7  <- 1/length(diff_M7) * sum(abs(diff_M7))

rmse_age_M7 <- sapply(1:dim(diff_M7)[2], function(a)
  sqrt(1/length(diff_M7[,a,]) * sum(((diff_M7[,a,])^2))))
mae_age_M7  <- sapply(1:dim(diff_M7)[2], function(a)
  1/length(diff_M7[,a,]) * sum(abs(diff_M7[,a,])))

### 3.4) Scenarios ----

n.sim <- 5000

cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl, c('n.sim', 'n.hor', 'regions', 'xv', 'etxr_ft', 'lca_M7',
                    'obj.K_M7', 'obj.kt_M7', 'cb.all', 'wave.all', 'i.fc'))
clusterEvalQ(cl, expr = {
  library(dplyr)
  library(forecast)
  library(MASS)
})
mlst.sim_M7 <- pblapply(cl = cl, X = 1:n.sim, FUN = function(s) {
  
  # Simulate the common and regional factors
  sim.K  <- simulate(obj.K_M7, nsim = n.hor)
  sim.kt <- sapply(seq_along(regions), function(r) simulate(obj.kt_M7[[r]], nsim = n.hor))
  
  # Force of mortality with bootstrapped DLNM coefficients and Gaussian error
  mu.e <- mu <- array(NA, dim = c(n.hor, length(xv), length(regions)))
  
  for(r in seq_along(regions)){
    X <- cbind(1, cb.all[[r]][i.fc, ], wave.all[[r]]$HWD[i.fc], wave.all[[r]]$CWD[i.fc])
    for(x in seq_along(xv)){
      cf <- MASS::mvrnorm(1, lca_M7$coef[, x, r],
                          lca_M7$vcov[[(r - 1) * length(xv) + x]])
      lp <- lca_M7$A[x, r] + lca_M7$B[x] * sim.K +
        lca_M7$b[x, r] * sim.kt[, r] + X %*% cf
      mu[, x, r]   <- exp(lp)
      mu.e[, x, r] <- exp(lp + rnorm(n.hor, 0, lca_M7$sigma[x, r]))
    }
  }
  
  # Return: 'A' without and 'B' with negative binomial observation noise
  list('A' = mu * etxr_ft,
       'B' = mu.e * etxr_ft)
})
stopCluster(cl)

mhat.sim_M7 <- simplify2array(lapply(mlst.sim_M7, function(x) x$B))
dimnames(mhat.sim_M7) <- list(as.character(weeks.fc), as.character(xv), regions, NULL)

# Quantiles, coverage and interval score
quant_M7 <- apply(mhat.sim_M7, c(1, 2, 3), quantile, probs = c(0.025, 0.5, 0.975))
lo_M7    <- quant_M7[1,,,]
hi_M7 <- quant_M7[3,,,]

cov_M7    <- mean(dt < hi_M7 & dt > lo_M7)
iscore_M7 <- mean((hi_M7 - lo_M7) +
                    (2/0.05) * (lo_M7 - dt) * (dt < lo_M7) +
                    (2/0.05) * (dt - hi_M7) * (dt > hi_M7))

cov_age_M7    <- sapply(1:dim(dt)[2], function(a) mean(dt[,a,] < hi_M7[,a,] & dt[,a,] > lo_M7[,a,]))
iscore_age_M7 <- sapply(1:dim(dt)[2], function(a)
  mean((hi_M7[,a,] - lo_M7[,a,]) +
         (2/0.05) * (lo_M7[,a,] - dt[,a,]) * (dt[,a,] < lo_M7[,a,]) +
         (2/0.05) * (dt[,a,] - hi_M7[,a,]) * (dt[,a,] > hi_M7[,a,])))

# Scoring rules
grid <- expand.grid(age = as.character(xv), region = regions, stringsAsFactors = FALSE)
grid[, 'crps'] <- sapply(1:nrow(grid), function(i)
  mean(crps_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M7[, grid$age[i], grid$region[i], ])))
grid[, 'logs'] <- sapply(1:nrow(grid), function(i)
  mean(logs_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M7[, grid$age[i], grid$region[i], ])))

crps_M7 <- mean(grid$crps)
logs_M7 <- mean(grid$logs)

crps_age_M7 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_age_M7 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)

### Plot
df.sim    <- as.data.frame.table(quant_M7, responseName = "value") 
df.deaths <- as.data.frame.table(dtxr_ft)
df.expo   <- as.data.frame.table(etxr_ft)
colnames(df.sim) <- c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo) <- c('ISODate', 'age', 'region', 'expo')

df.tot <- plyr::join_all(list(df.sim, df.deaths, df.expo)) %>%
  tidyr::pivot_wider(names_from = Quantile,
                     values_from = sim,
                     names_prefix = "q")

age_vec <- c('50', '70', '90')
reg     <- c('11')
dfs <- df.tot %>% 
  dplyr::filter(age %in% age_vec, region == reg)
DfM$age <- as.factor(DfM$age)
dfs$age    <- as.factor(dfs$age)

ggplot(dfs) + 
  theme_bw(base_size = 15) + 
  ggtitle('Country-level') +
  geom_line(data = DfM %>% dplyr::filter(age %in% age_vec, region == reg),
            aes(x = as.Date(ISODate), y = log(deaths), group = age, col = age))  + 
  geom_point(aes(x = as.Date(ISODate), y = log(obs), group = age, col = age),
             shape = 3, col = 'black', cex = 0.2) + 
  geom_ribbon(aes(x = as.Date(ISODate), ymin = log(`q2.5%`), ymax = log(`q97.5%`), 
                  group = age, fill = age, col = age), alpha = 0.2) +
  geom_line(aes(x = as.Date(ISODate), y = log(`q50%`), group = age, col = age)) +
  coord_cartesian(xlim= c(as.Date('2010-01-01'), as.Date('2019-12-31'))) + 
  scale_fill_manual(name = 'Age', 
                    values = c("#33A02C", "#1F78B4", "#E31A1C"),
                    labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
  scale_color_manual(name = 'Age', 
                     values = c("#33A02C", "#1F78B4", "#E31A1C"),
                     labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) 
