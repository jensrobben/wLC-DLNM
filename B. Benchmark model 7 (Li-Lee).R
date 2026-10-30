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

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

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


##### 3) Model 0: annual, regional Li-Lee ----
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
source('F. Annual, regional Li-Lee (NB).R')

# Fit
lca_M0 <- aLL.NB(xv = xv, isoyv = isoyv, isowv = isowv, etx = etxr_M0, 
                 dtx = dtxr_M0, eps = 10^(-9), verbose = TRUE)

# In-sample fit
bic_M0 <- -2 * lca_M0$ll + log(length(dtxr_M0)) * lca_M0$npar

### 3.3) Point forecast ----
# Forecast horizon
nY       <- 5
n.hor    <- 52*nY + 1
nTT      <- max(isoyv) + nY
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Fit time series: common factor (random walk with drift), 
# regional factors (mean-reverting AR(1), Li & Lee, 2005)
obj.K_M0   <- forecast::Arima(ts(lca_M0$K), 
                              order = c(0,1,0),
                              include.drift = TRUE)
obj.kt_M0  <- lapply(1:dim(dtxr_M0)[3], function(r)
  forecast::Arima(ts(lca_M0$kappa2[,r]), 
                  order = c(1,0,0),
                  include.mean = FALSE))

# Point forecast K and kappa
K_fc_M0  <- c(lca_M0$K, c(forecast::forecast(obj.K_M0, nY)$mean))
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
  mhat_M0[,j,] <- exp(sweep(lca_M0$B[j] * K_fc_M0[c(isoyv,isoyv_fc)] + 
                              sweep(kt_fc_M0[c(isoyv,isoyv_fc),], 2, 
                                    lca_M0$beta2[j,], '*'), 
                            2, lca_M0$beta1[j,], '+'))
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
                   'dtxr_ft_M0', 'obj.K_M0', 'obj.kt_M0', 'lca_M0', 'nY',
                   'isoyv', 'isoyv_fc'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(copula)
  library(dlnm)
  library(forecast)
})
mlst.sim_M0 <- pblapply(cl=cl, X=1:n.sim, FUN = function(s) {
  
  # Simulate K-process (common) and kt-processes (regional)
  sim.K    <- simulate(object = obj.K_M0, nsim = nY)
  sim.kt   <- sapply(1:length(regions), function(r) {
    simulate(object = obj.kt_M0[[r]], nsim = nY)
  })
  
  # Force of mortality
  idx <- isoyv_fc - min(isoyv_fc) + 1
  mu  <- simplify2array(lapply(1:length(regions), function(r) 
    exp(outer(sim.K[idx], lca_M0$B) + outer(sim.kt[idx, r], lca_M0$beta2[,r]) + 
          rep(lca_M0$beta1[,r], each = length(idx)))))
  
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