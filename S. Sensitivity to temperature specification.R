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
source('F. Weekly Lee-Carter (NB).R')
source('F. Weekly DLNM (NB).R')
source('F. Weekly Lee-Carter + DLNM (NB).R')

##### 1) Influenza data ----

# Read (no data for overseas region 'OUTRE-MER)
df_ia <- fread('../Sentinel/inc-3-RDD.csv', sep = ",", skip = 1) %>%
  mutate('Year' = as.integer(substr(week, 1, 4)),
         'Week' = as.integer(substr(week, 5, 6)),
         'Date' = ISOweek2date(sprintf("%d-W%02d-%d", Year, Week, 1)),
         'ia100' = inc100,
         .before = 1) %>%
  dplyr::filter(! geo_name %in% c('OUTRE-MER','CORSE')) %>% 
  dplyr::select(c(Year, Week, Date, ia100, geo_name)) %>%
  dplyr::arrange(Date, geo_name)

# Shapefile
shapefile <- read_sf('../NUTS/NUTS_RG_03M_2024_3035.shp') %>%
  dplyr::filter(LEVL_CODE == 1, 
                CNTR_CODE == "FR",
                ! NUTS_ID %in% c('FRY', 'FRM'))
plot(shapefile[,c('NUTS_ID')])

# Match with NUTS 1
df_ia$Region <- plyr::mapvalues(df_ia$geo_name, 
                                from = sort(unique(df_ia$geo_name)),
                                to = sort(unique(shapefile$NUTS_NAME)))
df_ia <- df_ia %>%
  dplyr::left_join(shapefile[,c('NUTS_ID','NUTS_NAME')], 
                   by = c('Region' = 'NUTS_NAME')) %>%
  arrange(Year, Week) %>%
  dplyr::filter(Year >= 1985 & Year <= 2024) %>%
  mutate(ia100 = as.integer(ia100))

# Example data for one region
region <- c('FRC')
dfr    <- df_ia %>% dplyr::filter(NUTS_ID == region)
plot(dfr$Date, as.numeric(dfr$ia100), type = 'l')

# Check if data for every isoyear, isoweek
check1 <- df_ia %>% group_by(Year, Week) %>%
  summarize(l = length(ia100))
all(check1$l == 12) #21 is number of regions we consider

# Missing value imputation per region
check2 <- df_ia %>%
  group_by(Region) %>%
  dplyr::summarize(sum(is.na(ia100))*1) %>%
  as.data.frame()
check2 

df_ia <- df_ia %>%
  group_by(Region) %>%
  arrange(Date) %>%
  mutate(ia100 = na_seasplit(ts(ia100), algorithm = "interpolation")) %>%
  ungroup()

# Recode
df_ia <- df_ia %>% 
  dplyr::mutate(Region = plyr::mapvalues(geo_name, 
                                         from = sort(unique(df_ia$geo_name)),
                                         to = gsub("-", " ", sort(unique(df_ia$geo_name))))) %>%
  dplyr::select(-c(geo_name))

# Plot
dfiar <- df_ia %>%
  dplyr::filter(NUTS_ID %in% c('FRE', 'FRH', 'FRK'),
                Year %in% 1990:2019)

pia <- ggplot(dfiar) +
  theme_bw(base_size = 15) + xlab('Date') + 
  ylab('ILI cases / 100 inhabitants') + 
  geom_line(aes(x = Date, y = ia100/1000, group = NUTS_ID, col = NUTS_ID)) + 
  scale_color_manual(name = '',
                     values = c('FRE' = "#FB9A99",
                                'FRH' = "#FF7F00",
                                'FRK' = "#1F78B4"),
                     labels = c('FRE' = '32',
                                'FRH' = '53',
                                'FRK' = '84')) +
  guides(col = guide_legend(override.aes = list(linewidth = 1.5, size = 5)))

# ggsave('../Figures/ia100.pdf', 
#        plot = pia, height = 12.5, width = 25, units = 'cm')

##### 2) Temperature data ----

# Load daily data
df.temp <- fread('../Temperature/FR_NUTS1_tg_daily.txt')

# Add ISOyear, ISOWeek, ISODate
df.temp <- df.temp %>%
  dplyr::mutate(ISOYear = isoyear(Date),
                ISOWeek = isoweek(Date),
                ISODate = ISOweek2date(paste0(sprintf("%04d-W%02d", ISOYear, ISOWeek), "-1")))

# Aggregate on week basis
df.temp <- df.temp %>%
  dplyr::group_by(ISODate, ISOYear, ISOWeek, Region) %>%
  dplyr::arrange(Date) %>%
  dplyr::reframe(T1 = max(tg),
                 T2 = min(tg), 
                 T3 = max((tg + lag(tg, 1) + lag(tg, 2)) / 3, na.rm = TRUE),
                 T4 = mean(tg)) %>%
  ungroup()


##### 3) Construct final data set ----

# Read data
Df <- fread('../Results updated/Df_final.txt')

# Join temp and influenza
df.var <-  df_ia %>%
  dplyr::left_join(df.temp, by = c('Date' = 'ISODate', 'Year' = 'ISOYear', 
                                   'Week' = 'ISOWeek', 'NUTS_ID' = 'Region'))

# Change dates
df.var$Date <- as.IDate(df.var$Date)
df.var      <- df.var %>% as.data.table()

# Add temp + influenza
Df <- Df %>% 
  dplyr::left_join(df.var, by = c('ISODate' = 'Date', 'ISOYear' = 'Year', 
                                  'ISOWeek' = 'Week', 'region_name' = 'Region')) %>%
  na.omit()

# Remove geometry
Df <- Df %>% dplyr::select(-c('geometry'))

# Interpolate 1997 June July data - misreporting
Df[, is_target := (ISOYear == 1997 & ISOWeek %in% 25:32)]  

interp_fun <- function(x) {
  ts_data <- ts(x, frequency = 52)
  as.numeric(pmax(round(forecast::na.interp(ts_data)), 0.001))
}

Df[, deaths := {
  x <- deaths
  x[is_target] <- NA
  interp_fun(x)
}, by = .(region, gender, age)] %>%
  dplyr::select(-c('is_target'))

# Save
fwrite(Df, file = '../Results updated/Df_var_new.txt')

##### 4) Fit age-enhanced DLNM ----
### 4.1) Make everything ready to fit ----
# Load data
Df <- fread(file = '../Results updated/Df_var_new.txt')
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
  reframe('mia100' = quantile(ia100,0.9)) %>%
  as.data.table()
DfM <- DfM %>%
  dplyr::left_join(sum, by = c('ISOWeek', 'region'))

DfM$ia100 <- DfM$ia100 - DfM$mia100

### 4.2) Construct cross-basis matrices ----

DfM$tg <- DfM$T3

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
                      scale = max(infl[region == r,ia100]))
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

cb1  <- center_cb(cb1, isowv)
cb2  <- center_cb(cb2, isowv)

### 4.3) Fit model ----

## Optimal lambda
lambda.opt <- c('Temp' = 10^(1.75), 'Infl' = 10^(3.75)) 
gc()

# Fit model
# obj <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                        etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                        cb = list('Temp' = cb1, 'Infl' = cb2),
#                        lambda = lambda.opt)

# saveRDS(object = obj, file = '../Results updated/lca_MaxTemp.rds')
# saveRDS(object = obj, file = '../Results updated/lca_MinTemp.rds')
# saveRDS(object = obj, file = '../Results updated/lca_3DayMax.rds')

# Load models
lca.Tavg <- readRDS(file = '../Results updated/lca_final.rds')
lca.Tmax <- readRDS(file = '../Results updated/lca_MaxTemp.rds')
lca.Tmin <- readRDS(file = '../Results updated/lca_MinTemp.rds')
lca.3Day <- readRDS(file = '../Results updated/lca_3DayMax.rds')

# Log-likelihoods
c(lca.Tavg$lcd$ll.up, lca.Tmax$lcd$ll.up, 
  lca.Tmin$lcd$ll.up, lca.3Day$lcd$ll.up)

# In-sample plot
reg_vec <- dimnames(dtxr)[[3]][c(5,8,11)]
age_vec <- colnames(dtxr)[c(5,9)]

Df <- data.frame(
  expand.grid('Date'   = as.Date(dimnames(dtxr)[[1]]), 
              'Age'    = dimnames(dtxr)[[2]], 
              'Region' = dimnames(dtxr)[[3]]),
  'deaths'   = as.vector(dtxr),
  'expo'     = as.vector(etxr),
  'rate'     = as.vector(dtxr)/as.vector(etxr),
  'fit.Tavg' = as.vector(lca.Tavg$lcd$mhat),
  'fit.Tmax' = as.vector(lca.Tmax$lcd$mhat),
  'fit.Tmin' = as.vector(lca.Tmin$lcd$mhat),
  'fit.3Day' = as.vector(lca.3Day$lcd$mhat)) %>%
  pivot_longer(
    cols = starts_with("fit."),
    names_to = "model",
    values_to = "fit"
  )

Df$model <- factor(Df$model, levels = c('fit.3Day','fit.Tmax',
                                        'fit.Tavg', 'fit.Tmin'))

pvis1 <- list()
for(a in 1:length(age_vec)){
  dfs1 <- Df %>% dplyr::filter(Region == reg_vec[1],
                                   Age == age_vec[a])
  
  pvis1[[3*(a-1) + 1]] <- ggplot(dfs1 %>% dplyr::filter(Date >= as.Date('2015-01-01'))) + 
    theme_bw(base_size = 15) + 
    xlab('') + 
    ylab(bquote(log(mu[.(age_vec[a])*',t,w,r']))) + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'rate'), linewidth = 1.2) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = model, group = model)) +
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[1]], '')) + 
    scale_color_manual(name = '', 
                       values = c('rate' = "gray90", 'fit.Tavg' = agreen,
                                  'fit.Tmax' = ared, 'fit.Tmin' = ablue,
                                  'fit.3Day' = aorange),
                       labels = c('rate' = 'Observed', 'fit.Tavg' = 'Tavg',
                                  'fit.Tmax' = 'Tmax', 'fit.Tmin' = 'Tmin',
                                  'fit.3Day' = '3Day'),
                       breaks = c('rate', 'fit.Tmax', 'fit.Tmin', 'fit.3Day', 'fit.Tavg')) + 
    coord_cartesian(ylim = quantile(Df %>% 
                                      dplyr::filter(Age == age_vec[a], 
                                                    Region %in% reg_vec) %>% 
                                      pull(rate) %>% log(), c(0.005,1))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  dfs2 <- Df %>% dplyr::filter(Region == reg_vec[2],
                               Age == age_vec[a])
  pvis1[[3*(a-1) + 2]] <- ggplot(dfs2 %>% dplyr::filter(Date >= as.Date('2015-01-01'))) + 
    theme_bw(base_size = 15) + 
    xlab('') + 
    ylab(bquote(log(mu[.(age_vec[a])*',t,w,r']))) + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'rate'), linewidth = 1.2) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = model, group = model)) +
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[2]], '')) + 
    scale_color_manual(name = '', 
                       values = c('rate' = "gray90", 'fit.Tavg' = agreen,
                                  'fit.Tmax' = ared, 'fit.Tmin' = ablue,
                                  'fit.3Day' = aorange),
                       labels = c('rate' = 'Observed', 'fit.Tavg' = 'Tavg',
                                  'fit.Tmax' = 'Tmax', 'fit.Tmin' = 'Tmin',
                                  'fit.3Day' = '3Day'),
                       breaks = c('rate', 'fit.Tmax', 'fit.Tmin', 'fit.3Day', 'fit.Tavg')) + 
    coord_cartesian(ylim = quantile(Df %>% 
                                      dplyr::filter(Age == age_vec[a], 
                                                    Region %in% reg_vec) %>% 
                                      pull(rate) %>% log(), c(0.005,1))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  
  dfs3 <- Df %>% dplyr::filter(Region == reg_vec[3],
                               Age == age_vec[a])
  
  pvis1[[3*(a-1) + 3]] <- ggplot(dfs3 %>% dplyr::filter(Date >= as.Date('2015-01-01'))) + 
    theme_bw(base_size = 15) + 
    xlab('') + 
    ylab(bquote(log(mu[.(age_vec[a])*',t,w,r']))) + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'rate'), linewidth = 1.2) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = model, group = model)) +
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[3]], '')) + 
    scale_color_manual(name = '', 
                       values = c('rate' = "gray90", 'fit.Tavg' = agreen,
                                  'fit.Tmax' = ared, 'fit.Tmin' = ablue,
                                  'fit.3Day' = aorange),
                       labels = c('rate' = 'Observed', 'fit.Tavg' = 'Tavg',
                                  'fit.Tmax' = 'Tmax', 'fit.Tmin' = 'Tmin',
                                  'fit.3Day' = '3Day'),
                       breaks = c('rate', 'fit.Tmax', 'fit.Tmin', 'fit.3Day', 'fit.Tavg')) + 
    coord_cartesian(ylim = quantile(Df %>% 
                                      dplyr::filter(Age == age_vec[a], 
                                                    Region %in% reg_vec) %>% 
                                      pull(rate) %>% log(), c(0.005,1))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
}

plot.in.sample <- ggpubr::ggarrange(plotlist = pvis1, common.legend = TRUE, 
                                    legend = 'bottom', nrow = 2, ncol = 3, 
                                    align = 'h', widths = c(1,1,1)) 
plot.in.sample


# Check RMSE at extreme temperatures, per region and age
grid <- expand.grid(region = dimnames(dtxr)[[3]],
                    age    = dimnames(dtxr)[[2]])

temp <- df.temp %>%
  dplyr::filter(ISODate >= as.Date('1990-01-01'),
                ISODate <= as.Date('2019-12-29'),
                Region %in% shapef$NUTS_ID) %>%
  dplyr::left_join(shapef[,c('NUTS_ID', 'region')],
                   by = c('Region' = 'NUTS_ID'))

res <- lst()

for(j in 1:nrow(grid)){
  tavgr  <- temp %>% dplyr::filter(region == grid$region[j])
  
  dates1 <- unique(tavgr$ISODate[tavgr$T1 > quantile(tavgr$T1, 0.95) | 
                                   tavgr$T2 > quantile(tavgr$T2, 0.95) |
                                   tavgr$T3 > quantile(tavgr$T3, 0.95) |
                                   tavgr$T4 > quantile(tavgr$T4, 0.95)])
  dates2 <- unique(tavgr$ISODate[tavgr$T1 > quantile(tavgr$T1, 0.975) | 
                                   tavgr$T2 > quantile(tavgr$T2, 0.975) |
                                   tavgr$T3 > quantile(tavgr$T3, 0.975) |
                                   tavgr$T4 > quantile(tavgr$T4, 0.975)])
  dates3 <- unique(tavgr$ISODate[tavgr$T1 < quantile(tavgr$T1, 0.05) | 
                                   tavgr$T2 < quantile(tavgr$T2, 0.05) |
                                   tavgr$T3 < quantile(tavgr$T3, 0.05) |
                                   tavgr$T4 < quantile(tavgr$T4, 0.05)])
  dates4 <- unique(tavgr$ISODate[tavgr$T1 < quantile(tavgr$T1, 0.025) | 
                                   tavgr$T2 < quantile(tavgr$T2, 0.025) |
                                   tavgr$T3 < quantile(tavgr$T3, 0.025) |
                                   tavgr$T4 < quantile(tavgr$T4, 0.025)])
  
  Dfra1 <- Df %>%
    dplyr::filter(Age == grid$age[j],
                  Region == grid$region[j],
                  Date %in% dates1)
  Dfra2 <- Df %>%
    dplyr::filter(Age == grid$age[j],
                  Region == grid$region[j],
                  Date %in% dates2)
  Dfra3 <- Df %>%
    dplyr::filter(Age == grid$age[j],
                  Region == grid$region[j],
                  Date %in% dates3)
  
  Dfra4 <- Df %>%
    dplyr::filter(Age == grid$age[j],
                  Region == grid$region[j],
                  Date %in% dates4)
  
  vec1 <- Dfra1 %>%
    dplyr::group_by(model) %>%
    dplyr::reframe(rmse = sqrt(1/n() * 
                                 sum((deaths - expo * fit)^2)),
                   mae = 1/n() * sum(abs(deaths - expo * fit))) 
  vec2 <- Dfra2 %>%
    dplyr::group_by(model) %>%
    dplyr::reframe(rmse = sqrt(1/n() * 
                                 sum((deaths - expo * fit)^2)),
                   mae = 1/n() * sum(abs(deaths - expo * fit))) 
  vec3 <- Dfra3 %>%
    dplyr::group_by(model) %>%
    dplyr::reframe(rmse = sqrt(1/n() * 
                                 sum((deaths - expo * fit)^2)),
                   mae = 1/n() * sum(abs(deaths - expo * fit))) 
  vec4 <- Dfra4 %>%
    dplyr::group_by(model) %>%
    dplyr::reframe(rmse = sqrt(1/n() * 
                                 sum((deaths - expo * fit)^2)),
                   mae = 1/n() * sum(abs(deaths - expo * fit))) 
  
  rmse_vec1 <- setNames(vec1$rmse, gsub("\\.fit", "", 
                                        paste0('RMSE.', vec1$model, '.q1')))
  mae_vec1  <- setNames(vec1$mae, gsub("\\.fit", "", 
                                       paste0('MAE.', vec1$model, '.q1')))
  rmse_vec2 <- setNames(vec2$rmse, gsub("\\.fit", "", 
                                        paste0('RMSE.', vec2$model, '.q2')))
  mae_vec2  <- setNames(vec2$mae, gsub("\\.fit", "", 
                                       paste0('MAE.', vec2$model, '.q2')))
  rmse_vec3 <- setNames(vec3$rmse, gsub("\\.fit", "", 
                                        paste0('RMSE.', vec3$model, '.q3')))
  mae_vec3  <- setNames(vec3$mae, gsub("\\.fit", "", 
                                       paste0('MAE.', vec3$model, '.q3')))
  rmse_vec4 <- setNames(vec4$rmse, gsub("\\.fit", "", 
                                        paste0('RMSE.', vec4$model, '.q4')))
  mae_vec4  <- setNames(vec4$mae, gsub("\\.fit", "", 
                                       paste0('MAE.', vec4$model, '.q4')))
  
  res[[j]] <- c(rmse_vec1, mae_vec1, rmse_vec2, mae_vec2, 
                rmse_vec3, mae_vec3, rmse_vec4, mae_vec4)
}

Df.res <- bind_rows(res)
Dff    <- bind_cols(grid, Df.res)

tab1 <- Dff %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(mean(RMSE.Tavg.q1),
                 mean(RMSE.Tmax.q1),
                 mean(RMSE.Tmin.q1),
                 mean(RMSE.3Day.q1),
                 mean(MAE.Tavg.q1),
                 mean(MAE.Tmax.q1),
                 mean(MAE.Tmin.q1),
                 mean(MAE.3Day.q1)) %>%
  dplyr::select(-c('age'))

tab1 <- rbind(tab1,
              c(mean(Dff$RMSE.Tavg.q1), mean(Dff$RMSE.Tmax.q1),
                mean(Dff$RMSE.Tmin.q1), mean(Dff$RMSE.3Day.q1),
                mean(Dff$MAE.Tavg.q1), mean(Dff$MAE.Tmax.q1),
                mean(Dff$MAE.Tmin.q1), mean(Dff$MAE.3Day.q1))) %>%
  as.data.frame()
rownames(tab1) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                    '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
colnames(tab1) <- NULL
xtable(tab1)
apply(tab1[,1:4], 1, which.min)
apply(tab1[,5:8], 1, which.min)

tab2 <- Dff %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(mean(RMSE.Tavg.q2),
                 mean(RMSE.Tmax.q2),
                 mean(RMSE.Tmin.q2),
                 mean(RMSE.3Day.q2),
                 mean(MAE.Tavg.q2),
                 mean(MAE.Tmax.q2),
                 mean(MAE.Tmin.q2),
                 mean(MAE.3Day.q2)) %>%
  dplyr::select(-c('age'))

tab2 <- rbind(tab2,
              c(mean(Dff$RMSE.Tavg.q2), mean(Dff$RMSE.Tmax.q2),
                mean(Dff$RMSE.Tmin.q2), mean(Dff$RMSE.3Day.q2),
                mean(Dff$MAE.Tavg.q2), mean(Dff$MAE.Tmax.q2),
                mean(Dff$MAE.Tmin.q2), mean(Dff$MAE.3Day.q2))) %>%
  as.data.frame()
rownames(tab2) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                    '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(tab2)
apply(tab2[,1:4], 1, which.min)
apply(tab2[,5:8], 1, which.min)

tab3 <- Dff %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(mean(RMSE.Tavg.q3),
                 mean(RMSE.Tmax.q3),
                 mean(RMSE.Tmin.q3),
                 mean(RMSE.3Day.q3),
                 mean(MAE.Tavg.q3),
                 mean(MAE.Tmax.q3),
                 mean(MAE.Tmin.q3),
                 mean(MAE.3Day.q3)) %>%
  dplyr::select(-c('age'))

tab3 <- rbind(tab3,
              c(mean(Dff$RMSE.Tavg.q3), mean(Dff$RMSE.Tmax.q3),
                mean(Dff$RMSE.Tmin.q3), mean(Dff$RMSE.3Day.q3),
                mean(Dff$MAE.Tavg.q3), mean(Dff$MAE.Tmax.q3),
                mean(Dff$MAE.Tmin.q3), mean(Dff$MAE.3Day.q3))) %>%
  as.data.frame()
rownames(tab3) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                    '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')
xtable(tab3)
apply(tab3[,1:4], 1, which.min)
apply(tab3[,5:8], 1, which.min)

tab4 <- Dff %>%
  dplyr::group_by(age) %>%
  dplyr::reframe(mean(RMSE.Tavg.q4),
                 mean(RMSE.Tmax.q4),
                 mean(RMSE.Tmin.q4),
                 mean(RMSE.3Day.q4),
                 mean(MAE.Tavg.q4),
                 mean(MAE.Tmax.q4),
                 mean(MAE.Tmin.q4),
                 mean(MAE.3Day.q4)) %>%
  dplyr::select(-c('age'))

tab4 <- rbind(tab4,
              c(mean(Dff$RMSE.Tavg.q4), mean(Dff$RMSE.Tmax.q4),
                mean(Dff$RMSE.Tmin.q4), mean(Dff$RMSE.3Day.q4),
                mean(Dff$MAE.Tavg.q4), mean(Dff$MAE.Tmax.q4),
                mean(Dff$MAE.Tmin.q4), mean(Dff$MAE.3Day.q4))) %>%
  as.data.frame()
rownames(tab4) <- c('50--54', '55--59', '60--64', '65--69', '70--74',
                    '75--79', '80--84', '85--89', '90--94', '95+', 'Overall')

xtable(tab4)
apply(tab4[,1:4], 1, which.min)
apply(tab4[,5:8], 1, which.min)

grid$best_model <- fit_cols[max.col(-as.matrix(grid[, fit_cols]))]

prop.table(table(grid$best_model))

mean(grid$fit.3Day)
mean(grid$fit.Tavg)
mean(grid$fit.Tmax)
mean(grid$fit.Tmin)
