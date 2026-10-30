##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'spdep',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'data.table', 'xgboost',
              'parallel', 'doParallel', 'readxl', 'tidyr', 'purrr', 'imputeTS',
              'dlnm','Matrix', 'metR', 'pbapply', 'scoringRules')
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


##### 3) Model 0: annual, regional Lee-Carter ----
### 3.1) Adapt data ----

# Data
DfM_M0 <- DfM 

# Death counts and exposures - matrix format
dtxr_M0     <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM_M0)
etxr_M0     <- stats::xtabs(expo ~ ISODate + age + region, data = DfM_M0)

dtxr_ft_M0  <- deaths.obs
etxr_ft_M0  <- expo.obs

### 3.2) Fit ----
# Source
source('F. Annual, regional Lee-Carter (NB).R')

# Fit
lca_M0 <- aLC.NB(xv = xv, isoyv = isoyv, isowv = isowv, etx = etxr_M0, 
                 adjmat = L, dtx = dtxr_M0, eps = 10^(-9), verbose = TRUE)


### 3.3) Point forecast ----
# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Fit time series
obj.kt_M0  <- lapply(1:dim(dtxr_M0)[3], function(r)
  forecast::Arima(ts(lca_M0$kappa2[,r]), 
                  order = c(0,1,1),
                  include.drift = TRUE))

# Point forecast kappa
kt_fc_M0 <- lapply(1:dim(dtxr_M0)[3], function(r) 
  c(forecast::forecast(obj.kt_M0[[r]], nY)$mean))
kt_fc_M0 <- matrix(unlist(kt_fc_M0), ncol = dim(dtxr_M0)[3])
kt_fc_M0 <- rbind(lca_M0$kappa2, kt_fc_M0)

# Point forecast death rates
mhat_M0 <- array(NA, dim = c(nrow(lca_M0$mhat) + n.hor, ncol(lca_M0$mhat), 12),
                 dimnames = list(c(dimnames(lca_M0$mhat)[[1]], 
                                   as.character(sort(unique(Dfnew$ISODate)))),
                                 dimnames(lca_M0$mhat)[[2]],
                                 dimnames(lca_M0$mhat)[[3]]))

for(j in 1:ncol(mhat_M0)){
  mhat_M0[,j,] <- exp(sweep(lca_M0$beta2[j] * kt_fc_M0[c(isoyv,isoyv_fc),], 2, 
                            lca_M0$beta1[j,], '+'))
}
mhat_M0 <- mhat_M0[-c(1:nrow(dtxr_M0)),,]

# Error metrics
dt <- dtxr_ft_M0
et <- etxr_ft_M0

diff_M0 <- mhat_M0*et - dt

rmse_M0 <- sqrt(1/length(diff_M0) * sum(((diff_M0)^2)))
mae_M0  <- 1/length(diff_M0) * sum(abs(diff_M0))

rmse_reg_M0 <- sapply(1:dim(diff_M0)[3], function(r) 
  sqrt(1/length(diff_M0[,,r]) * sum(((diff_M0[,,r])^2))))
rmse_age_M0 <- sapply(1:dim(diff_M0)[2], function(a) 
  sqrt(1/length(diff_M0[,a,]) * sum(((diff_M0[,a,])^2))))
mae_reg_M0 <- sapply(1:dim(diff_M0)[3], function(r) 
  1/length(diff_M0[,,r]) * sum(abs(diff_M0[,,r])))
mae_age_M0 <- sapply(1:dim(diff_M0)[2], function(a) 
  1/length(diff_M0[,a,]) * sum(abs(diff_M0[,a,])))

### 3.4) Scenarios ----

# Number of scenarios
n.sim <- 5000

# Do the simulations
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl,c('n.sim', 'n.hor', 'regions', 'etxr_ft_M0',
                   'dtxr_ft_M0', 'obj.kt_M0', 'lca_M0', 'nY',
                   'isoyv', 'isoyv_fc'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(copula)
  library(dlnm)
  library(forecast)
})
mlst.sim_M0 <- pblapply(cl=cl, X=1:n.sim, FUN = function(s) {
  
  # Simulate kt-process
  sim.kt   <- sapply(1:length(regions), function(r) {
    simulate(object = obj.kt_M0[[r]], nsim = nY)
  })
  
  # Force of mortality
  mu <- exp(sweep(aperm(lca_M0$beta2 %o% sim.kt[isoyv_fc - min(isoyv_fc) + 1,,drop=FALSE],
                        c(2, 1, 3)), c(2,3), lca_M0$beta1, '+'))

  # Return
  exp.ds <- exp(outer(lca_M0$phix, lca_M0$phir, "+"))
  exp.ds <- replicate(nrow(mu), exp.ds, simplify = "array")  # gives an n x R x T array
  exp.ds <- aperm(exp.ds, c(3, 1, 2))

  list('A' = mu * etxr_ft_M0,
       'B' = array(rnbinom(length(mu), mu = mu * etxr_ft_M0, size = exp.ds), 
                   dim = dim(mu)))
})
stopCluster(cl)
mhat.sim_M0 <- aperm(simplify2array(lapply(mlst.sim_M0, function(x) x$B)), c(4,1,2,3))

# Compute quantiles
quant.sim_M0 <- apply(mhat.sim_M0, MARGIN = c(2, 3, 4), 
                      FUN = function(x) 
                        quantile(x, probs = c(0.025, 0.5, 0.975)))
# Set dimnames
dimnames(mhat.sim_M0) <- list(1:n.sim,
                              as.character(unique(Dfnew$ISODate)),
                              unique(DfM_M0$age), 
                              unique(Dfnew$region))
dimnames(quant.sim_M0) <- list(rownames(quant.sim_M0), 
                               as.character(unique(Dfnew$ISODate)),
                               unique(DfM_M0$age), 
                               unique(Dfnew$region))

# Plot
df.sim    <- as.data.frame.table(quant.sim_M0, responseName = "value") 
df.deaths <- as.data.frame.table(dtxr_ft_M0)
df.expo   <- as.data.frame.table(etxr_ft_M0)
colnames(df.sim) <- c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo) <- c('ISODate', 'age', 'region', 'expo')

df.tot_M0 <- plyr::join_all(list(df.sim, df.deaths, df.expo)) %>%
  tidyr::pivot_wider(names_from = Quantile,
                     values_from = sim,
                     names_prefix = "q")

age_vec <- c('50', '70', '90')
reg     <- c('11')
dfs <- df.tot_M0 %>% 
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

# Coverage ratio
df.tot_M0$coverage <- df.tot_M0$obs < df.tot_M0$`q97.5%` & 
  df.tot_M0$obs > df.tot_M0$`q2.5%`
df.tot_M0$iscore <- with(df.tot_M0,
                         (`q97.5%` - `q2.5%`) +
                           (2 / 0.05) * (`q2.5%` - obs) * (obs < `q2.5%`) +
                           (2 / 0.05) * (obs - `q97.5%`) * (obs > `q97.5%`))

cov_M0    <- mean(df.tot_M0$coverage)
iscore_M0 <- mean(df.tot_M0$iscore)

cov_age_M0    <- df.tot_M0 %>% dplyr::group_by(age) %>% 
  dplyr::reframe(m = mean(coverage)) %>% pull(m)
iscore_age_M0 <- df.tot_M0 %>% dplyr::group_by(age) %>% 
  dplyr::reframe(m = mean(iscore)) %>% pull(m)

# Scoring rules
grid   <- expand.grid(age = unique(DfM$age), region = unique(DfM$region))
crps_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(crps_sample(y = dt[,x,r], dat = t(mhat.sim_M0[,,x,r])))
})
logs_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(logs_sample(y = dt[,x,r], dat = t(mhat.sim_M0[,,x,r])))
})
grid[,'crps'] <- crps_all
grid[,'logs'] <- logs_all

crps_M0 <- mean(grid$crps)
logs_M0 <- mean(grid$logs)

crps_reg_M0 <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
crps_age_M0 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_reg_M0  <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)
logs_age_M0  <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)


##### 5) Model 1: weekly, country-level Lee-Carter ----
### 5.1) Adapt data ----

# Data
DfM_M1 <- DfM %>%
  dplyr::group_by(ISODate, ISOYear, ISOWeek, age, gender) %>%
  dplyr::reframe(deaths = sum(deaths),
                 expo   = sum(expo)) %>%
  dplyr::mutate(region = '00')

# Death counts and exposures - matrix format
dtxr_M1     <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM_M1)
etxr_M1     <- stats::xtabs(expo ~ ISODate + age + region, data = DfM_M1)

dtxr_ft_M1  <- deaths.obs
etxr_ft_M1  <- expo.obs
dtx_ft_M1   <- apply(deaths.obs, c(1,2), sum)
dtx_ft_M1   <- array(dtx_ft_M1, dim = c(dim(dtx_ft_M1), 1))
etx_ft_M1   <- apply(expo.obs, c(1,2), sum)
etx_ft_M1   <- array(etx_ft_M1, dim = c(dim(etx_ft_M1), 1))

dimnames(dtx_ft_M1) = dimnames(etx_ft_M1) <- 
  list(as.character(sort(unique(Dfnew$ISODate))),
       unique(DfM_M1$age), 
       unique(DfM_M1$region))

### 5.2) Fit -----
# Source
source('F. Weekly Lee-Carter (NB).R')

# Fit
lca_M1 <- wLC.NB(xv = xv, isoyv = isoyv, isowv = isowv, etx = etxr_M1, 
                 dtx = dtxr_M1, int = xv*0, constraints = 'ALL', 
                 adjmat = L, eps = 10^(-9), verbose = TRUE)

### 5.3) Point forecast ----
# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Fit time series
obj.kt_M1  <- lapply(1:dim(dtxr_M1)[3], function(r)
  forecast::Arima(ts(lca_M1$kappa2[,r]), 
                  order = c(0,1,1),
                  include.drift = TRUE))

# Point forecast kappa
kt_fc_M1 <- lapply(1:dim(dtxr_M1)[3], function(r) 
  c(forecast::forecast(obj.kt_M1[[r]], nY)$mean))
kt_fc_M1 <- matrix(unlist(kt_fc_M1), ncol = dim(dtxr_M1)[3])
kt_fc_M1 <- rbind(lca_M1$kappa2, kt_fc_M1)

# Point forecast death rates
mhat_M1 <- array(NA, dim = c(nrow(lca_M0$mhat) + n.hor, ncol(lca_M0$mhat), 12),
                 dimnames = list(c(dimnames(lca_M0$mhat)[[1]], 
                                   as.character(sort(unique(Dfnew$ISODate)))),
                                 dimnames(lca_M0$mhat)[[2]],
                                 dimnames(lca_M0$mhat)[[3]]))

for(j in 1:ncol(mhat_M1)){
  mhat_M1[,j,] <- exp(lca_M1$beta2[j] * kt_fc_M1[c(isoyv, isoyv_fc),] + 
                        lca_M1$beta3[j] * lca_M1$kappa3[c(isowv, isowv_fc),] +
                        lca_M1$beta1[j,])
}

# Error metrics
dt <- dtxr_ft_M1
et <- etxr_ft_M1
mhat_rep_M1 <- array(rep(mhat_M1, dim(dt)[3]), 
                     dim = c(nrow(mhat_M1), ncol(mhat_M1), dim(dt)[3]))[-c(1:nrow(dtxr_M1)),,]
dimnames(mhat_rep_M1) <- list(rownames(dtxr_ft_M1), colnames(mhat_M1), unique(Dfnew$region))
diff_M1 <- mhat_rep_M1*et - dt

rmse_M1 <- sqrt(1/length(diff_M1) * sum(((diff_M1)^2)))
mae_M1  <- 1/length(diff_M1) * sum(abs(diff_M1))

rmse_reg_M1 <- sapply(1:dim(diff_M1)[3], function(r) 
  sqrt(1/length(diff_M1[,,r]) * sum(((diff_M1[,,r])^2))))
rmse_age_M1 <- sapply(1:dim(diff_M1)[2], function(a) 
  sqrt(1/length(diff_M1[,a,]) * sum(((diff_M1[,a,])^2))))
mae_reg_M1 <- sapply(1:dim(diff_M1)[3], function(r) 
  1/length(diff_M1[,,r]) * sum(abs(diff_M1[,,r])))
mae_age_M1 <- sapply(1:dim(diff_M1)[2], function(a) 
  1/length(diff_M1[,a,]) * sum(abs(diff_M1[,a,])))

### 5.4) Scenarios ----

# Number of scenarios
n.sim <- 5000

# Do the simulations
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl,c('n.sim', 'n.hor', 'regions', 'etx_ft_M1', 'etxr_ft_M1',
                   'dtx_ft_M1', 'obj.kt_M1', 'lca_M1', 'nY', 'isoyv', 'isowv',
                   'isoyv_fc', 'isowv_fc'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(copula)
  library(dlnm)
  library(forecast)
})
mlst.sim_M1 <- pblapply(cl=cl, X=1:n.sim, FUN = function(s) {
  
  # Simulate kt-process
  sim.kt   <- sapply(dim(etx_ft_M1)[3], function(r) {
    simulate(object = obj.kt_M1[[r]], nsim = nY)
  })
  
  # Force of mortality
  mu <- exp(sweep(aperm(lca_M1$beta2 %o% sim.kt[isoyv_fc - min(isoyv_fc) + 1,,drop=FALSE],
                        c(2, 1, 3)), c(2,3), lca_M1$beta1, '+') +
              aperm(lca_M1$beta3 %o% lca_M1$kappa3[isowv_fc,,drop=FALSE], 
                    c(2, 1, 3)))
  mu <- array(rep(mu, 12), dim = dim(etxr_ft_M1))
  

  # Return
  exp.ds <- exp(outer(lca_M1$phix, lca_M1$phir, "+"))
  exp.ds <- replicate(nrow(mu), exp.ds, simplify = "array")  # gives an n x R x T array
  exp.ds <- aperm(exp.ds, c(3, 1, 2))
  exp.ds <- array(rep(exp.ds, 12), dim = dim(etxr_ft_M1))

  list('A' = mu * etxr_ft_M1,
       'B' = array(rnbinom(length(mu), mu = mu * etxr_ft_M1, size = exp.ds), 
                   dim = dim(mu)))
})
stopCluster(cl)
mhat.sim_M1 <- aperm(simplify2array(lapply(mlst.sim_M1, function(x) x$B)), c(4,1,2,3))

# Compute quantiles
quant.sim_M1 <- apply(mhat.sim_M1, MARGIN = c(2, 3, 4), 
                      FUN = function(x) 
                        quantile(x, probs = c(0.025, 0.5, 0.975)))
# Set dimnames
dimnames(mhat.sim_M1) <- list(1:n.sim,
                              as.character(unique(Dfnew$ISODate)),
                              unique(DfM_M1$age), 
                              unique(Dfnew$region))
dimnames(quant.sim_M1) <- list(rownames(quant.sim_M1), 
                               as.character(unique(Dfnew$ISODate)),
                               unique(DfM_M1$age), 
                               unique(Dfnew$region))

# Plot
df.sim    <- as.data.frame.table(quant.sim_M1, responseName = "value") 
df.deaths <- as.data.frame.table(dtxr_ft_M1)
df.expo   <- as.data.frame.table(etxr_ft_M1)
colnames(df.sim) <- c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo) <- c('ISODate', 'age', 'region', 'expo')

df.tot_M1 <- plyr::join_all(list(df.sim, df.deaths, df.expo)) %>%
  tidyr::pivot_wider(names_from = Quantile,
                     values_from = sim,
                     names_prefix = "q")

age_vec <- c('50', '70', '90')
reg     <- c('11')
dfs <- df.tot_M1 %>% 
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

# Coverage ratio
df.tot_M1$coverage <- df.tot_M1$obs < df.tot_M1$`q97.5%` & 
  df.tot_M1$obs > df.tot_M1$`q2.5%`
df.tot_M1$iscore <- with(df.tot_M1,
                         (`q97.5%` - `q2.5%`) +
                           (2 / 0.05) * (`q2.5%` - obs) * (obs < `q2.5%`) +
                           (2 / 0.05) * (obs - `q97.5%`) * (obs > `q97.5%`))

cov_M1    <- mean(df.tot_M1$coverage)
iscore_M1 <- mean(df.tot_M1$iscore)

cov_age_M1    <- df.tot_M1 %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(coverage)) %>% pull(m)
iscore_age_M1 <- df.tot_M1 %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(iscore)) %>% pull(m)

# Scoring rules
grid   <- expand.grid(age = unique(DfM$age), region = unique(DfM$region))
crps_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(crps_sample(y = dt[,x,r], dat = t(mhat.sim_M1[,,x,r])))
})
logs_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(logs_sample(y = dt[,x,r], dat = t(mhat.sim_M1[,,x,r])))
})
grid[,'crps'] <- crps_all
grid[,'logs'] <- logs_all

crps_M1 <- mean(grid$crps)
logs_M1 <- mean(grid$logs)

crps_reg_M1 <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
crps_age_M1 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_reg_M1  <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)
logs_age_M1  <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)

##### 6) Model 2: weekly, region-specific Lee-Carter ----
### 6.1) Adapt data ----

# Data
DfM_M2 <- DfM 

# Death counts and exposures - matrix format
dtxr_M2     <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM_M2)
etxr_M2     <- stats::xtabs(expo ~ ISODate + age + region, data = DfM_M2)

dtxr_ft_M2  <- deaths.obs
etxr_ft_M2  <- expo.obs
dimnames(dtxr_ft_M2) = dimnames(etxr_ft_M2) <- 
  list(as.character(sort(unique(Dfnew$ISODate))),
       unique(DfM_M2$age), unique(DfM_M2$region))

### 6.2) Fit -----
# Source
source('F. Weekly Lee-Carter (NB).R')

# Fit
lca_M2 <- wLC.NB(xv = xv, isoyv = isoyv, isowv = isowv, etx = etxr_M2, 
                 dtx = dtxr_M2, int = xv*0, constraints = 'ALL', 
                 adjmat = L, eps = 10^(-9), verbose = TRUE)

### 6.3) Point forecast ----
# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Fit time series
obj.kt_M2  <- lapply(1:dim(dtxr_M2)[3], function(r)
  forecast::Arima(ts(lca_M2$kappa2[,r]), 
                  order = c(0,1,1),
                  include.drift = TRUE))

# Point forecast kappa
kt_fc_M2 <- lapply(1:dim(dtxr_M2)[3], function(r) 
  c(forecast::forecast(obj.kt_M2[[r]], nY)$mean))
kt_fc_M2 <- matrix(unlist(kt_fc_M2), ncol = dim(dtxr_M2)[3])
kt_fc_M2 <- rbind(lca_M2$kappa2, kt_fc_M2)

# Point forecast death rates
mhat_M2 <- exp(sweep(aperm(lca_M2$beta2 %o% kt_fc_M2[c(isoyv, isoyv_fc),],
                           c(2, 1, 3)), c(2,3), lca_M2$beta1, '+') +
                 aperm(lca_M2$beta3 %o% lca_M2$kappa3[c(isowv, isowv_fc),], 
                       c(2, 1, 3)))
dimnames(mhat_M2) <- list(c(dimnames(lca_M2$mhat)[[1]], 
                            as.character(sort(unique(Dfnew$ISODate)))),
                          dimnames(lca_M2$mhat)[[2]],
                          dimnames(lca_M2$mhat)[[3]])

# Error metrics
dt <- dtxr_ft_M2
et <- etxr_ft_M2
mhat_rep_M2 <- mhat_M2[-c(1:nrow(dtxr_M2)),,]

diff_M2 <- mhat_rep_M2*et - dt

rmse_M2 <- sqrt(1/length(diff_M2) * sum(((diff_M2)^2)))
mae_M2  <- 1/length(diff_M2) * sum(abs(diff_M2))

rmse_reg_M2 <- sapply(1:dim(diff_M2)[3], function(r) 
  sqrt(1/length(diff_M2[,,r]) * sum(((diff_M2[,,r])^2))))
rmse_age_M2 <- sapply(1:dim(diff_M2)[2], function(a) 
  sqrt(1/length(diff_M2[,a,]) * sum(((diff_M2[,a,])^2))))
mae_reg_M2 <- sapply(1:dim(diff_M2)[3], function(r) 
  1/length(diff_M2[,,r]) * sum(abs(diff_M2[,,r])))
mae_age_M2 <- sapply(1:dim(diff_M2)[2], function(a) 
  1/length(diff_M2[,a,]) * sum(abs(diff_M2[,a,])))

### 6.4) Scenarios ----

# Number of scenarios
n.sim <- 5000

# Do the simulations
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl,c('n.sim', 'n.hor', 'regions', 'etxr_ft_M2',
                   'dtxr_ft_M2', 'obj.kt_M2', 'lca_M2', 'nY',
                   'isoyv_fc', 'isowv_fc'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(copula)
  library(dlnm)
  library(forecast)
})
mlst.sim_M2 <- pblapply(cl=cl, X=1:n.sim, FUN = function(s) {
  
  # Simulate kt-process
  sim.kt   <- sapply(1:length(regions), function(r) {
    simulate(object = obj.kt_M2[[r]], nsim = nY)
  })
  
  # Force of mortality
  mu <- exp(sweep(aperm(lca_M2$beta2 %o% sim.kt[isoyv_fc - min(isoyv_fc) + 1,],
                        c(2, 1, 3)), c(2,3), lca_M2$beta1, '+') +
              aperm(lca_M2$beta3 %o% lca_M2$kappa3[isowv_fc,], 
                    c(2, 1, 3)))
  
  # Return
  exp.ds <- exp(outer(lca_M2$phix, lca_M2$phir, "+"))
  exp.ds <- replicate(nrow(mu), exp.ds, simplify = "array")  # gives an n x R x T array
  exp.ds <- aperm(exp.ds, c(3, 1, 2))

  list('A' = mu * etxr_ft_M2,
       'B' = array(rnbinom(length(mu), mu = mu * etxr_ft_M2, size = exp.ds), 
                   dim = dim(mu)))
})
stopCluster(cl)
mhat.sim_M2 <- aperm(simplify2array(lapply(mlst.sim_M2, function(x) x$B)), 
                     c(4,1,2,3))

# Compute quantiles
quant.sim_M2 <- apply(mhat.sim_M2, MARGIN = c(2, 3, 4), 
                      FUN = function(x) 
                        quantile(x, probs = c(0.025, 0.5, 0.975)))

# Set dim names
dimnames(mhat.sim_M2) <- list(1:n.sim,
                              as.character(unique(Dfnew$ISODate)),
                              unique(Dfnew$age), 
                              unique(Dfnew$region))
dimnames(quant.sim_M2) <- list(rownames(quant.sim_M2), 
                               as.character(unique(Dfnew$ISODate)),
                               unique(Dfnew$age), 
                               unique(Dfnew$region))

# Plot
df.sim    <- as.data.frame.table(quant.sim_M2, responseName = "value") 
df.deaths <- as.data.frame.table(dtxr_ft_M2)
df.expo   <- as.data.frame.table(etxr_ft_M2)
colnames(df.sim) <- c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo) <- c('ISODate', 'age', 'region', 'expo')

df.tot_M2 <- plyr::join_all(list(df.sim, df.deaths, df.expo)) %>%
  tidyr::pivot_wider(names_from = Quantile,
                     values_from = sim,
                     names_prefix = "q")

age_vec <- c('50', '70', '90')
reg     <- '11'
dfs <- df.tot_M2 %>% 
  dplyr::filter(age %in% age_vec, region == reg)
DfM_M2$age <- as.factor(DfM_M2$age)
dfs$age    <- as.factor(dfs$age)

ggplot(dfs) + 
  theme_bw(base_size = 15) + 
  ggtitle('Country-level') +
  geom_line(data = DfM_M2 %>% dplyr::filter(age %in% age_vec, region == reg),
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

# Coverage ratio
df.tot_M2$coverage <- df.tot_M2$obs < df.tot_M2$`q97.5%` & 
  df.tot_M2$obs > df.tot_M2$`q2.5%`
df.tot_M2$iscore <- with(df.tot_M2,
                         (`q97.5%` - `q2.5%`) +
                           (2 / 0.05) * (`q2.5%` - obs) * (obs < `q2.5%`) +
                           (2 / 0.05) * (obs - `q97.5%`) * (obs > `q97.5%`))

cov_M2    <- mean(df.tot_M2$coverage)
iscore_M2 <- mean(df.tot_M2$iscore)

cov_age_M2    <- df.tot_M2 %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(coverage)) %>% pull(m)
iscore_age_M2 <- df.tot_M2 %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(iscore)) %>% pull(m)

# Scoring rules
grid   <- expand.grid(age = unique(DfM$age), region = unique(DfM$region))
crps_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(crps_sample(y = dt[,x,r], dat = t(mhat.sim_M2[,,x,r])))
})
logs_all <- sapply(1:nrow(grid), function(i){
  x <- as.character(grid[i,1])
  r <- as.character(grid[i,2])
  mean(logs_sample(y = dt[,x,r], dat = t(mhat.sim_M2[,,x,r])))
})
grid[,'crps'] <- crps_all
grid[,'logs'] <- logs_all

crps_M2 <- mean(grid$crps)
logs_M2 <- mean(grid$logs)

crps_reg_M2 <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
crps_age_M2 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_reg_M2  <- grid %>% dplyr::group_by(region) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)
logs_age_M2  <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)

##### 7) Model 3: weekly, region-specific Lee-Carter with DLNM ----

# RMSE & MAE
rmse_age_M3 <- readRDS('../Results updated/RMSE_age_Models4-6.rds')
mae_age_M3  <- readRDS('../Results updated/MAE_age_Models4-6.rds')
rmse_M3     <- readRDS('../Results updated/RMSE_overall_Models4-6.rds')
mae_M3      <- readRDS('../Results updated/MAE_overall_Models4-6.rds')

# CRPS
crps_age_obs_M3 <- readRDS('../Results updated/CRPS_age_obs_Models4-6.rds')
crps_age_sim_M3 <- readRDS('../Results updated/CRPS_age_sim_Models4-6.rds')
crps_obs_M3     <- readRDS('../Results updated/CRPS_overall_obs_Models4-6.rds')
crps_sim_M3     <- readRDS('../Results updated/CRPS_overall_sim_Models4-6.rds')

# LogS
logs_age_obs_M3 <- readRDS('../Results updated/LogS_age_obs_Models4-6.rds')
logs_age_sim_M3 <- readRDS('../Results updated/LogS_age_sim_Models4-6.rds')
logs_obs_M3     <- readRDS('../Results updated/LogS_overall_obs_Models4-6.rds')
logs_sim_M3     <- readRDS('../Results updated/LogS_overall_sim_Models4-6.rds')

# Coverage ratio
cov_age_obs_M3 <- readRDS('../Results updated/COV_age_obs_Models4-6.rds')
cov_age_sim_M3 <- readRDS('../Results updated/COV_age_sim_Models4-6.rds')
cov_obs_M3     <- readRDS('../Results updated/COV_overall_obs_Models4-6.rds')
cov_sim_M3     <- readRDS('../Results updated/COV_overall_sim_Models4-6.rds')

# Interval score
iscore_age_obs_M3 <- readRDS('../Results updated/IScore_age_obs_Models4-6.rds')
iscore_age_sim_M3 <- readRDS('../Results updated/IScore_age_sim_Models4-6.rds')
iscore_obs_M3     <- readRDS('../Results updated/IScore_overall_obs_Models4-6.rds')
iscore_sim_M3     <- readRDS('../Results updated/IScore_overall_sim_Models4-6.rds')

##### 9) Make table -----
library(xtable)

### RMSE & MAE
dt <- cbind(do.call('cbind', list(rmse_age_M0, rmse_age_M1,
                                  rmse_age_M2, rmse_age_M3)),
            do.call('cbind', list(mae_age_M0, mae_age_M1,
                                  mae_age_M2, mae_age_M3)))
dt <- rbind(dt, c(rmse_M0, rmse_M1, rmse_M2, rmse_M3,
                  mae_M0, mae_M1, mae_M2, mae_M3))
rownames(dt) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                  '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(dt)
apply(dt[,1:6], 1, which.min)
apply(dt[,7:12], 1, which.min)

### CRPS
dt <- do.call('cbind', list(crps_age_M0, crps_age_M1,
                            crps_age_M2, crps_age_obs_M3,
                            crps_age_sim_M3))
dt <- rbind(dt, c(crps_M0, crps_M1, crps_M2, crps_obs_M3, crps_sim_M3))
rownames(dt) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                  '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(dt)
apply(dt, 1, which.min)

### LogS
dt <- do.call('cbind', list(logs_age_M0, logs_age_M1,
                            logs_age_M2, logs_age_obs_M3,
                            logs_age_sim_M3))
dt <- rbind(dt, c(logs_M0, logs_M1, logs_M2, logs_obs_M3, logs_sim_M3))
rownames(dt) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                  '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(dt)
apply(dt, 1, which.min)

### Coverage ratio
dt <- do.call('cbind', list(cov_age_M0, cov_age_M1,
                            cov_age_M2, cov_age_obs_M3,
                            cov_age_sim_M3))
dt <- rbind(dt, c(cov_M0, cov_M1, cov_M2, cov_obs_M3, cov_sim_M3))
rownames(dt) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                  '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(dt)
apply(dt, 1, which.max)

### Interval score
dt <- do.call('cbind', list(iscore_age_M0, iscore_age_M1,
                            iscore_age_M2, iscore_age_obs_M3,
                            iscore_age_sim_M3))
dt <- rbind(dt, c(iscore_M0, iscore_M1, iscore_M2, iscore_obs_M3, iscore_sim_M3))
rownames(dt) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                  '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(dt)
apply(dt, 1, which.min)

##### 10) Make figure ----

dt0 <- as.data.frame.table(mhat_M0) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Model 1')) 
dt1 <- as.data.frame.table(mhat_rep_M1) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Model 2')) 
dt2 <- as.data.frame.table(mhat_M2) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Model 3')) 
dt3 <- as.data.frame.table(mhat_M3) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Model 4')) 
dt4 <- as.data.frame.table(mhat_M4) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Model 5')) 
dto <- as.data.frame.table(dtxr_ft_M0/etxr_ft_M0) %>%
  `colnames<-`(c('ISODate', 'age', 'region', 'Observed')) 

dt <- plyr::join_all(list(dt0, dt1, dt2, dt3, dt4, dto),
                     by = c('ISODate', 'age', 'region'))
dt_sub <- dt %>% 
  dplyr::filter(region == '11', age == '90')

ggplot(dt_sub) + 
  theme_bw(base_size = 15) + 
  xlab('Date') + 
  ylab(bquote('log('*mu['x,t,w,r']*')'))  +
  geom_line(aes(x = as.Date(ISODate), y = log(`Observed`), col = 'Obs'), 
            linewidth = 1, alpha = 0.5) + 
  ggtitle(shapef$NUTS_NAME[shapef$region == '11']) +
  geom_line(aes(x = as.Date(ISODate), y = log(`Model 1`), col = 'M1')) + 
  geom_line(aes(x = as.Date(ISODate), y = log(`Model 2`), col = 'M2')) + 
  geom_line(aes(x = as.Date(ISODate), y = log(`Model 3`), col = 'M3')) + 
  geom_line(aes(x = as.Date(ISODate), y = log(`Model 4`), col = 'M4')) + 
  geom_line(aes(x = as.Date(ISODate), y = log(`Model 5`), col = 'M5')) + 
  coord_cartesian(xlim= c(as.Date('2012-01-01'), as.Date('2020-01-01'))) 

pvis <- list()
reg_vec <- c('32', '53', '84')
age_vec <- c('50', '70', '90')

for(r in 1:length(reg_vec)){
  dfs <- df.tot_M3 %>% dplyr::filter(region == reg_vec[r],
                                     age %in% age_vec)
  
  pvis[[r]] <- ggplot(dfs) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Date','')) + 
    ylab(if(r==1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(data = DfM_M4 %>% dplyr::filter(region == reg_vec[r], age %in% c('50', '70', '90')),
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

# ggsave('../Figures updated/predperform.pdf',
#        plot = plot.forecast, height = 12, width = 30, units = 'cm')

##### 11) Uncertainty decomposition ----

### Models 1-3

# Extract
mhat.simA_M0 <- aperm(simplify2array(lapply(mlst.sim_M0, function(x) x$A)), 
                      c(4,1,2,3))
mhat.simA_M1 <- aperm(simplify2array(lapply(mlst.sim_M1, function(x) x$A)),
                      c(4,1,2,3))
mhat.simA_M2 <- aperm(simplify2array(lapply(mlst.sim_M2, function(x) x$A)),
                      c(4,1,2,3))

# Quantiles
quant.simA_M0 <- apply(mhat.simA_M0, MARGIN = c(2, 3, 4), 
                       FUN = function(x) 
                         quantile(x, probs = c(0.025, 0.5, 0.975)))
quant.simA_M1 <- apply(mhat.simA_M1, MARGIN = c(2, 3, 4), 
                       FUN = function(x) 
                         quantile(x, probs = c(0.025, 0.5, 0.975)))
quant.simA_M2 <- apply(mhat.simA_M2, MARGIN = c(2, 3, 4), 
                       FUN = function(x) 
                         quantile(x, probs = c(0.025, 0.5, 0.975)))

# Dim names
dimnames(mhat.simA_M0) = dimnames(mhat.simA_M1) = dimnames(mhat.simA_M2) <- 
  list(1:n.sim, as.character(unique(Dfnew$ISODate)), 
       unique(Dfnew$age), unique(Dfnew$region))
dimnames(quant.simA_M0) = dimnames(quant.simA_M1) = dimnames(quant.simA_M2) <- 
  list(rownames(quant.sim_M2), as.character(unique(Dfnew$ISODate)),
       unique(Dfnew$age), unique(Dfnew$region))

# Data frame
df.sim_M0  <- as.data.frame.table(quant.sim_M0, responseName = "value") 
df.simA_M0 <- as.data.frame.table(quant.simA_M0, responseName = "value") 
df.sim_M1  <- as.data.frame.table(quant.sim_M1, responseName = "value")  
df.simA_M1 <- as.data.frame.table(quant.simA_M1, responseName = "value")
df.sim_M2  <- as.data.frame.table(quant.sim_M2, responseName = "value")
df.simA_M2 <- as.data.frame.table(quant.simA_M2, responseName = "value")

colnames(df.sim_M0) = colnames(df.simA_M0) = colnames(df.sim_M1) = 
  colnames(df.simA_M1) = colnames(df.sim_M2) = colnames(df.simA_M2) <- 
  c('Quantile', 'ISODate', 'age', 'region', 'sim')


# Join
df.sim_all <- lapply(list(df.sim_M0, df.simA_M0, df.sim_M1, df.simA_M1,
                          df.sim_M2, df.simA_M2), 
                     function(x) plyr::join_all(list(x, df.deaths, df.expo),
                                                by = c('ISODate', 'age', 'region')) %>%
                            tidyr::pivot_wider(names_from = Quantile,
                                               values_from = sim,
                                               names_prefix = "q"))
names(df.sim_all) <- c('Model1', 'Model1A', 'Model2', 
                       'Model2A', 'Model3', 'Model3A')

# Join
dt.sim_M1 <- df.sim_all$Model1A %>%
  dplyr::left_join(df.sim_all$Model1, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.sim_M2 <- df.sim_all$Model2A %>%
  dplyr::left_join(df.sim_all$Model2, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))
dt.sim_M3 <- df.sim_all$Model3A %>%
  dplyr::left_join(df.sim_all$Model3, 
                   by = c('ISODate', 'age', 'region', 'obs', 'expo'))

# Proportions
sum_M1 <- dt.sim_M1 %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum_M2 <- dt.sim_M2 %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum_M3 <- dt.sim_M3 %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(prop1 = sum(`q97.5%.x` - `q2.5%.x`)/sum(`q97.5%.y` - `q2.5%.y`)) 
sum.obs_M4 <- readRDS('../Results updated/UDecom_obs_Model4.rds')
sum.obs_M5 <- readRDS('../Results updated/UDecom_obs_Model5.rds')
sum.obs_M6 <- readRDS('../Results updated/UDecom_obs_Model6.rds')
sum.sim_M4 <- readRDS('../Results updated/UDecom_sim_Model4.rds')
sum.sim_M5 <- readRDS('../Results updated/UDecom_sim_Model5.rds')
sum.sim_M6 <- readRDS('../Results updated/UDecom_sim_Model6.rds')


sum <- list(sum_M1, sum_M2, sum_M3, sum.obs_M4, sum.obs_M5, sum.obs_M6,
            sum.sim_M4, sum.sim_M5, sum.sim_M6) %>%
  plyr::join_all(by = 'age') %>%
  as.data.frame() 
rownames(sum) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                   '75--79', '80--84', '85--89', '90--94', '95+')
xtable(sum[,-1])

