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
adjmat             <- nb2mat(spdep::poly2nb(shapef, queen = TRUE), 
                             style = 'B', zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat)   <- list(regions, regions)
diag(adjmat)       <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))

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
Dfnew <- Df %>%
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


##### 3) Model 9: three-state regime-switching model of Robben et al. (2026) ----

### 3.1) Adapt data ----

# Source
source('F. Robben et al. approach.R')

# Daily temperature (NUTS 1)
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
weeks.cal <- as.character(sort(unique(DfM$ISODate)))
weeks.fc  <- as.character(sort(unique(Dfnew$ISODate)))
weeks.all <- c(weeks.cal, weeks.fc)

# Excess heat and cold indices per region
feat.lst <- vector('list', length(regions))

for(j in seq_along(regions)){
  dj <- df.daily[region == regions[j]]
  wk <- ISOweek::ISOweek2date(paste0(sprintf("%04d-W%02d", dj$ISOYear, dj$ISOWeek), "-1"))
  
  ex <- excess_index(dj$tg, wk, ref = dj$tg[as.character(wk) %in% weeks.cal])
  rownames(ex) <- as.character(sort(unique(wk)))
  ex <- ex[weeks.all, ]
  
  # Lags 0, 1 and the average of lags 2-3
  feat.lst[[j]] <- data.frame(
    Date   = weeks.all,
    Region = as.character(regions[j]),
    w_avg_ehi      = add_lags(ex$EHI)[, 'l0'],
    w_avg_ehi_l1   = add_lags(ex$EHI)[, 'l1'],
    w_avg_eci      = add_lags(ex$ECI)[, 'l0'],
    w_avg_eci_l1   = add_lags(ex$ECI)[, 'l1'],
    w_avg_eci_l2.3 = add_lags(ex$ECI)[, 'l2.3'])
}
df.feat <- do.call('rbind', feat.lst)

# Influenza activity index
df.ia <- Df %>%
  dplyr::filter(gender == 2, age == 50, ISOYear >= 1990, ISOYear <= 2019) %>%
  dplyr::select(ISODate, ISOWeek, region, ia100) %>%
  dplyr::arrange(region, ISODate) %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(ia100q = influenza_index(ia100, ISOWeek,
                                         ref.idx = which(as.character(ISODate) 
                                                         %in% weeks.cal))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(ia100q_l1   = dplyr::lag(ia100q, 1),
                ia100q_l2.3 = (dplyr::lag(ia100q, 2) + dplyr::lag(ia100q, 3))/2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Region = as.character(region), Date = as.character(ISODate)) %>%
  dplyr::select(Date, Region, ia100q, ia100q_l1, ia100q_l2.3)

# Death counts and exposures in long format
df.rs <- rbind(DfM[, c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'deaths', 'expo')],
               Dfnew[, c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'deaths', 'expo')]) %>%
  `colnames<-`(c('Date', 'ISOYear', 'ISOWeek', 'Region', 'Age', 'Deaths', 'Expo')) %>%
  dplyr::mutate(Region = as.character(Region), Date = as.character(Date)) %>%
  dplyr::left_join(df.feat, by = c('Date', 'Region')) %>%
  dplyr::left_join(df.ia,   by = c('Date', 'Region')) %>%
  dplyr::mutate(Age = factor(Age, levels = xv)) %>%
  dplyr::arrange(Region, Age, Date) %>%
  na.omit()

ages <- levels(df.rs$Age)

### 3.2) Baseline model ----
base_M9 <- fit_baseline_rs(df = df.rs %>% dplyr::filter(Date %in% weeks.cal),
                           adjmat = adjmat, regions = regions, ages = ages,
                           verbose = TRUE)

# Baseline deaths over the calibration and the forecast period
df.rs$bDeaths <- predict_baseline_rs(base_M9, df.rs, regions, ages)

### 3.3) Regime-switching model ----

# State restrictions: a transition is only possible when the driving index is non-zero
df.rs$penalty01 <- ifelse(df.rs$w_avg_ehi > 0, 0, -1000)
df.rs$penalty02 <- ifelse(df.rs$ia100q > 0 | df.rs$w_avg_eci > 0 |
                            df.rs$ia100q_l1 > 0 | df.rs$w_avg_eci_l1 > 0, 0, -1000)

# State-specific Poisson models (no COVID-19 terms: the data end in 2019)
form1_M9 <- as.formula(Deaths ~ offset(log(bDeaths)) +
                         w_avg_ehi:Age +
                         w_avg_ehi_l1:Age)

form2_M9 <- as.formula(Deaths ~ offset(log(bDeaths)) +
                         w_avg_eci:Age +
                         ia100q:Age +
                         w_avg_eci_l1:Age +
                         ia100q_l1:Age +
                         w_avg_eci_l2.3:Age +
                         ia100q_l2.3:Age)

df.rs$Deaths[df.rs$Deaths < 1] <- 0

lca_M9 <- fit_regime_switching(df = df.rs %>% dplyr::filter(Date %in% weeks.cal),
                               form1 = form1_M9, form2 = form2_M9,
                               regions = regions, ages = ages, verbose = TRUE)

### 3.4) Point forecast and scenarios ----

n.sim <- 5000

df.fc <- df.rs %>% dplyr::filter(Date %in% weeks.fc)
nT.fc <- length(unique(df.fc$Date))
nA    <- length(ages)
R     <- length(regions)

# Design matrices over the forecast period
mm <- function(f, d) {
  m <- model.matrix(f, data = d)
  n <- nrow(m)/R
  simplify2array(lapply(1:R, function(j) matrix(m[(n*(j-1) + 1):(n*j), ], ncol = ncol(m))))
}
xnew1 <- mm(form1_M9, df.fc)
xnew2 <- mm(form2_M9, df.fc)
dfa   <- df.fc %>% dplyr::filter(Age == ages[nA])
xx01n <- mm(~ w_avg_ehi, dfa)
xx11n <- mm(~ w_avg_ehi, dfa)
xx02n <- mm(~ w_avg_eci, dfa)
xx22n <- mm(~ w_avg_eci_l2.3, dfa)

pen01n <- matrix(dfa$penalty01, nrow = nT.fc, ncol = R)
pen02n <- matrix(dfa$penalty02, nrow = nT.fc, ncol = R)
bnew   <- matrix(df.fc$bDeaths, nrow = nT.fc*nA, ncol = R)

# Most probable state at the end of the calibration period
s0_M9 <- apply(matrix(lca_M9$ms.prob[nrow(lca_M9$ms.prob), ], ncol = 3), 1, which.max) - 1

cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl, c('n.sim', 'lca_M9', 'xnew1', 'xnew2', 'xx01n', 'xx02n', 'xx11n',
                    'xx22n', 'pen01n', 'pen02n', 'bnew', 's0_M9', 'nT.fc', 'R', 'nA',
                    'simulate_rs'))
mlst.sim_M9 <- pblapply(cl = cl, X = 1:n.sim, FUN = function(s)
  simulate_rs(lca_M9, xnew1, xnew2, xx01n, xx02n, xx11n, xx22n,
              pen01n, pen02n, bnew, s0_M9, nT.fc, R, nA))
stopCluster(cl)

# Reshape to (week, age, region, simulation)
to.array <- function(m) aperm(array(m, dim = c(nT.fc, nA, R)), c(1, 2, 3))
mhat.sim_M9 <- simplify2array(lapply(mlst.sim_M9, function(x) to.array(x$deaths)))
dimnames(mhat.sim_M9) <- list(as.character(weeks.fc), as.character(xv), regions, NULL)

# Best estimate: 50% quantile of the simulated trajectories (their Tables 6 and 7)
dt <- dtxr_ft
et <- etxr_ft

mhat_M9 <- apply(mhat.sim_M9, c(1, 2, 3), median)
diff_M9 <- mhat_M9 - dt

rmse_M9 <- sqrt(1/length(diff_M9) * sum(((diff_M9)^2)))
mae_M9  <- 1/length(diff_M9) * sum(abs(diff_M9))

rmse_age_M9 <- sapply(1:dim(diff_M9)[2], function(a)
  sqrt(1/length(diff_M9[,a,]) * sum(((diff_M9[,a,])^2))))
mae_age_M9  <- sapply(1:dim(diff_M9)[2], function(a)
  1/length(diff_M9[,a,]) * sum(abs(diff_M9[,a,])))

# Quantiles, coverage ratio and interval score
quant_M9 <- apply(mhat.sim_M9, c(1, 2, 3), quantile, probs = c(0.025, 0.5, 0.975))
lo_M9    <- quant_M9[1,,,]
hi_M9    <- quant_M9[3,,,]

cov_M9    <- mean(dt < hi_M9 & dt > lo_M9)
iscore_M9 <- mean((hi_M9 - lo_M9) +
                    (2/0.05) * (lo_M9 - dt) * (dt < lo_M9) +
                    (2/0.05) * (dt - hi_M9) * (dt > hi_M9))

cov_age_M9    <- sapply(1:dim(dt)[2], function(a)
  mean(dt[,a,] < hi_M9[,a,] & dt[,a,] > lo_M9[,a,]))
iscore_age_M9 <- sapply(1:dim(dt)[2], function(a)
  mean((hi_M9[,a,] - lo_M9[,a,]) +
         (2/0.05) * (lo_M9[,a,] - dt[,a,]) * (dt[,a,] < lo_M9[,a,]) +
         (2/0.05) * (dt[,a,] - hi_M9[,a,]) * (dt[,a,] > hi_M9[,a,])))

# Scoring rules
grid <- expand.grid(age = as.character(xv), region = regions, stringsAsFactors = FALSE)
grid[, 'crps'] <- sapply(1:nrow(grid), function(i)
  mean(crps_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M9[, grid$age[i], grid$region[i], ])))
grid[, 'logs'] <- sapply(1:nrow(grid), function(i)
  mean(logs_sample(y = dt[, grid$age[i], grid$region[i]],
                   dat = mhat.sim_M9[, grid$age[i], grid$region[i], ])))

crps_M9 <- mean(grid$crps)
logs_M9 <- mean(grid$logs)

crps_age_M9 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(crps)) %>% pull(m)
logs_age_M9 <- grid %>% dplyr::group_by(age) %>% dplyr::reframe(m = mean(logs)) %>% pull(m)


### Plot
df.sim    <- as.data.frame.table(quant_M9, responseName = "value") 
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
reg     <- c('75')
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
  coord_cartesian(xlim = c(as.Date('2010-01-01'), as.Date('2019-12-31')),
                  ylim = c(1,5.75)) + 
  scale_fill_manual(name = 'Age', 
                    values = c("#33A02C", "#1F78B4", "#E31A1C"),
                    labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
  scale_color_manual(name = 'Age', 
                     values = c("#33A02C", "#1F78B4", "#E31A1C"),
                     labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94'))

