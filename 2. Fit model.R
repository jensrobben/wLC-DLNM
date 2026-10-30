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

# Color scales
colors <- c('11' = "#A6CEE3", '84' = "#1F78B4", '75' = "#B2DF8A", '76' = "#33A02C", 
            '32' = "#FB9A99", '44' = "#E31A1C", '93' = "#FDBF6F", '53' = "#FF7F00", 
            '28' = "#CAB2D6", '52' = "#6A3D9A", '27' = "#FFFF99", '24' = "#B15928")

# Source
source('F. Weekly Lee-Carter (NB).R')
source('F. Weekly DLNM (NB).R')
source('F. Weekly Lee-Carter + DLNM (NB).R')
source('F. Weekly Lee-Carter (POI).R')
source('F. Weekly DLNM (POI).R')
source('F. Weekly Lee-Carter + DLNM (POI).R')
source('F. Bootstrap CI.R')


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

# Show overdispersion
disp <- DfM %>% 
  dplyr::group_by(region, age) %>%
  dplyr::reframe('Ratio' = var(deaths)/mean(deaths))
mat <- matrix(disp$Ratio, 
              nrow = length(unique(DfM$age)), 
              ncol = length(unique(DfM$region)))
dimnames(mat) <- list(unique(DfM$age), unique(DfM$region))

library(xtable)
print(xtable(mat, caption = "Descriptive Table", align = "ccccccccccccc", 
             label = "tab:mytable"), 
      booktabs = TRUE, 
      include.rownames = TRUE,
      sanitize.text.function = identity)

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

### 1.3) Find optimal values for lambda ----

# # Grid
# lgrid <- expand.grid('Temp' = 10^seq(0, 4, by = 0.25),
#                      'Infl' = 10^seq(0, 4, by = 0.25))
# 
# # run parallel
# library(pbapply)
# cl <- snow::makeCluster(15)
# doSNOW::registerDoSNOW(cl)
# snow::clusterEvalQ(cl,  c(library(dplyr), library(lubridate),
#                           library(data.table), library(MASS)))
# snow::clusterExport(cl, list = c('xv', 'isoyv', 'isowv', 'etxr', 'dtxr',
#                                  'cb1', 'cb2', 'cb12', 'lgrid', 'L',
#                                  'llmax.beta1', 'llmax.beta2',
#                                  'llmax.beta3', 'llmax.beta.dlnm',
#                                  'llmax.kappa2', 'llmax.kappa3',
#                                  'llmax.gamma.dlnm', 'llmax.phix',
#                                  'llmax.phir', 'fit_wLC.DLNM.NB',
#                                  'wLC.DLNM.NB', 'wLC.NB', 'dlnm_penalty'),
#                     envir = environment())
# 
# pb <- progress::progress_bar$new(format = "Grid search = :letter [:bar] :elapsed | eta: :eta",
#                                  total = nrow(lgrid), width = 60)
# progress_letter <- 1:nrow(lgrid)
# progress <- function(n){
#   pb$tick(tokens = list(letter = progress_letter[n]))
# }
# opts <- list(progress = progress)
# 
# results <- foreach(k = 1:nrow(lgrid), .options.snow = opts) %dopar% {
#   lca <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                          etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                          cb = list('Temp' = cb1, 'Infl' = cb2),
#                          lambda = unlist(lgrid[k,]))
#   lca$lcd
#   }
# stopCluster(cl)
# 
# saveRDS(results, '../Results updated/modelsel1.rds')
# results <- readRDS('../Results updated/modelsel1.rds')
# 
# # Optimal values
# npar.base <- length(results[[1]]$beta1) + length(results[[1]]$beta2) +
#   length(results[[1]]$beta3) + length(results[[1]]$kappa2) +
#   length(results[[1]]$kappa3) + length(results[[1]]$phix) +
#   length(results[[1]]$phir) - 2*R - 2 - 1
# 
# library(pbapply)
# cl <- snow::makeCluster(15)
# doSNOW::registerDoSNOW(cl)
# snow::clusterEvalQ(cl,  c(library(dplyr), library(lubridate),
#                           library(data.table), library(MASS)))
# snow::clusterExport(cl, list = c('npar.base', 'hessian.wLC.DLNM',
#                                  'etxr', 'dtxr', 'results',
#                                  'sumzero', 'contrast_dlnm',
#                                  'edf_constrained', 'hessian.wLC.DLNM.up',
#                                  'hessian.beta', 'hessian.gamma'),
#                     envir = environment())
# pb <- progress::progress_bar$new(format = "Grid search = :letter [:bar] :elapsed | eta: :eta",
#                                  total = nrow(lgrid), width = 60)
# progress_letter <- 1:nrow(lgrid)
# progress <- function(n){
#   pb$tick(tokens = list(letter = progress_letter[n]))
# }
# opts <- list(progress = progress)
# 
# val <-  foreach(k = 1:length(results), .options.snow = opts) %dopar% {
#   H.dlnm    <- hessian.wLC.DLNM(etx = etxr, dtx = dtxr, fit = results[[k]])
#   H.dlnm.up <- hessian.wLC.DLNM.up(etx = etxr, dtx = dtxr, fit = results[[k]])
#   edf       <- edf_constrained(H.dlnm, H.dlnm.up, contrast_dlnm(results[[k]]))
#   -2 * results[[k]]$ll.up + log(length(results[[k]]$mhat)) *
#     (edf + npar.base)
# }
# stopCluster(cl)
# 
# lgrid$val <- unlist(val)
# 
# ijk <- which(lgrid$Temp == lgrid$Infl)
# lst <- results[ijk]
# beta.Temp <- data.frame('Age' = rep(xv, times = length(ijk)),
#                         'Val' = unlist(lapply(lst, function(x) x$betal[['Temp']])),
#                         'Phi' = as.factor(rep(lgrid$Temp[ijk], each = length(xv))))
# 
# # Surface plot
# space_comma <- function(x) {
#   format(x, big.mark = " ", scientific = FALSE, trim = TRUE)
# }
# lgrid$log_Temp <- log10(lgrid$Temp)
# lgrid$log_Infl <- log10(lgrid$Infl)
# min_point <- lgrid[which.min(lgrid$val), ]
# 
# ps <- ggplot(lgrid,
#              aes(x = log_Temp, y = log_Infl, fill = val)) +
#   geom_tile(color = "white", size = 0.5) +
#   scale_fill_gradientn(
#     colors = c("gray30", "gray50", "gray70", "tomato", "red"),
#     values = scales::rescale(c(
#       max(lgrid$val), quantile(lgrid$val, 0.75),
#       quantile(lgrid$val, 0.5), quantile(lgrid$val, 0.25),min(lgrid$val))),
#     name = "BIC",
#     labels = space_comma) +
#   geom_point(data = min_point,
#              aes(x = log_Temp, y = log_Infl),
#              color = "yellow", size = 5, shape = 8,stroke = 1.5) +
#   geom_contour(aes(z = val), color = "white", alpha = 0.3,
#                linetype = "dashed", bins = 8) +
#   labs(x = expression(log[10](psi[1])~"(Temperature smoothing)"),
#        y = expression(log[10](psi[2])~"(Influenza smoothing)"),
#        title = "BIC Surface",
#        subtitle = paste0("Optimal: (",
#                         round(min_point$log_Temp, 2), ", ",
#                         round(min_point$log_Infl, 2), ")")) +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         legend.position = "right",
#         legend.key.width = unit(0.8, "cm"),
#         legend.key.height = unit(2, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill = NA, color = "gray80",
#                                     linewidth = 0.5)) +
#   coord_fixed(ratio = 1, expand = TRUE)
# ps
# 
# ggsave('../Figures updated/smoothing.pdf',
#        plot = ps, height = 18, width = 18, units = 'cm')

### 1.4) Fit model ----

# Optimal lambda
lambda.opt <- c('Temp' = 10^(1.75), 'Infl' = 10^(3.75)) 
gc()

# Fit model
# obj <- fit_wLC.DLNM.NB(xv = xv, isoyv = isoyv, isowv = isowv,
#                        etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
#                        cb = list('Temp' = cb1, 'Infl' = cb2),
#                        lambda = lambda.opt)
# saveRDS(obj, '../Results updated/lca_final.rds')
obj  <- readRDS('../Results updated/lca_final.rds')
obj0 <- readRDS('../Results updated/lca_orig.rds')

### 1.5) CIs ----

## Notation
Q       <- unlist(lapply(obj$lcd$cb, function(x) dim(x)[2]))
regions <- dimnames(dtxr)[[3]]
yv      <- unique(isoyv) + min(DfM$ISOYear) - 1
wv      <- unique(isowv)
nA      <- dim(dtxr)[2]
nT      <- length(unique(isoyv))
nW      <- length(unique(isowv))
nQ      <- length(Q)

## Parameter estimates
Df.par <- data.frame(
  'par' = c(obj$lcd$beta1, obj$lcd$beta2, obj$lcd$kappa2,
            obj$lcd$beta3, obj$lcd$kappa3, unlist(obj$lcd$betal)),
  'type' = c(rep('alpha', nA * R), rep('beta2', nA),  rep('kappa2', nT * R),
             rep('beta3', nA), rep('kappa3', nW * R),
             sapply(seq_along(Q), function(x) 
               rep(paste0('beta',3 + which(names(Q[x]) == c('Temp', 'Infl', 'Int'))), nA))),
  'age' = c(rep(xv, times = R), xv, rep(yv * NA, times = R), xv,
            rep(wv * NA, times = R),  rep(xv, times = nQ)),
  'year' = c(rep(xv * NA, times = R), xv * NA, rep(yv, times = R), xv * NA,
             rep(wv * NA, times = R), rep(xv * NA, times = nQ)),
  'week' = c(rep(xv * NA, times = R), xv * NA,  rep(yv * NA, times = R), xv * NA,
             rep(wv, times = R),  rep(xv * NA, times = nQ)),
  'region' = c(rep(regions, each = nA), xv * NA, rep(regions, each = nT), xv * NA,
               rep(regions, each = nW), rep(xv * NA, times = nQ)))

## Hessians
H.base    <- hessian.wLC.NB(etx = etxr, dtx = dtxr, fit = obj$lcb)
H.dlnm    <- hessian.wLC.DLNM(etx = etxr, dtx = dtxr, fit = obj$lcd)
H.dlnm.up <- hessian.wLC.DLNM.up(etx = etxr, dtx = dtxr, fit = obj$lcd)
C.base    <- contrast_base(obj$lcb)
C.dlnm    <- contrast_dlnm(obj$lcd)
edf.dlnm  <- edf_constrained(H.dlnm, H.dlnm.up, C.dlnm)
npar0     <- ncol(C.base) + ncol(C.dlnm) + 10 + 12 - 1
npar1     <- ncol(C.base) + edf.dlnm + 10 + 12 - 1
bic0      <- -2 * obj0$lcd$ll.up + log(length(dtxr)) * npar0  
bic1      <- -2 * obj$lcd$ll.up + log(length(dtxr)) * npar1  

I1  <- cov_constrained(H.base, C.base)
I2  <- cov_constrained(H.dlnm, C.dlnm)
sd1 <- sqrt(diag(I1))
sd2 <- sqrt(diag(I2))

Df.sd <- data.frame('sd' = unname(c(sd1, sd2)),
                    'Name' = c(names(sd1), names(sd2))) %>%
  dplyr::mutate(type = tolower(sub("_.*", "", Name)),
                age  = as.integer(sub(".*Age([0-9]+).*", "\\1", Name)),
                year = as.integer(sub(".*Year([0-9]+).*", "\\1", 
                                      Name)) + min(DfM$ISOYear) - 1,
                week = as.integer(sub(".*Week([0-9]+).*", "\\1", Name)),
                region = ifelse(grepl("Region", Name),
                                sub(".*Region([0-9]+).*", "\\1", Name), NA)) %>%
  dplyr::filter(! type %in% c('gamma1', 'gamma2', 'gamma3'))

Df.par <- Df.par %>%
  dplyr::left_join(Df.sd, by = c('type', 'age', 'year', 'week', 'region')) %>%
  dplyr::mutate('sd' = dplyr::coalesce(sd, 0),
                'lowerH' = par - 1.96 * sd,
                'upperH' = par + 1.96 * sd)

## Bootstrap CI
# btci <- bootstrap_nb(B = 5000, obj.unpen = obj0, obj.pen = obj, xv = xv,
#                      isoyv = isoyv, isowv = isowv, etx = etxr, adjmat = L)
# saveRDS(object = btci, file = '../Results updated/bootstrap_Temp_Infl.rds')

btci <- readRDS(file = '../Results updated/bootstrap_Temp_Infl.rds')

alpha.b <- sapply(btci, function(x) as.vector(x$beta1)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
beta2.b <- sapply(btci, function(x) as.vector(x$beta2)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
kappa2.b <- sapply(btci, function(x) as.vector(x$kappa2)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
beta3.b <- sapply(btci, function(x) as.vector(x$beta3)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
kappa3.b <- sapply(btci, function(x) as.vector(x$kappa3)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
betal.b <- sapply(btci, function(x) unlist(x$betal)) %>%
  matrixStats::rowQuantiles(probs = c(0.025,0.975))

par.b <- do.call('rbind', list(alpha.b, beta2.b, kappa2.b, 
                               beta3.b, kappa3.b, betal.b)) %>%
  as.data.frame() 
colnames(par.b) <- c('lowerB', 'upperB')

Df.par <- bind_cols(Df.par, par.b)

### 1.6) Plot parameter estimates and CI ----

# Function
plot_parameter <- function(Dfp, var, ci) {
  # Colours
  colors <- c('11' = "#A6CEE3", '84' = "#1F78B4", '75' = "#B2DF8A", '76' = "#33A02C", 
              '32' = "#FB9A99", '44' = "#E31A1C", '93' = "#FDBF6F", '53' = "#FF7F00", 
              '28' = "#CAB2D6", '52' = "#6A3D9A", '27' = "#FFFF99", '24' = "#B15928")
  
  # X-axis label
  xlab <- if(var == 'alpha' | grepl('beta', var, fixed = TRUE))
    'Age' else if(var == 'kappa2') 
      'Year' else if(var == 'kappa3') 'Week'
  xvar <- tolower(xlab)
  
  
  # Y-axis label
  ylab <- if(var == 'alpha') 
    bquote(alpha['x,r']) else if(var == 'beta2') 
      bquote(beta['x']) else if(var == 'beta3')
        bquote(gamma['x']) else if(var == 'beta4')
          bquote(delta['x']) else if(var == 'beta5')
            bquote(epsilon['x']) else if(var == 'beta6')
              bquote(zeta['x']) else if(var == 'kappa2')
                bquote(kappa['t,r']) else if(var == 'kappa3')
                  bquote(lambda['w,r'])
  
  # Group
  if(var %in% c('alpha', 'kappa2', 'kappa3')) {
    group <- sym('region')
    cls   <- colors } else {
      group <- NA
      cls   <- 'black'
    }
  
  # Type of confidence interval
  if(ci == 'boot') {
    lower <- 'lowerB'
    upper <- 'upperB'
  } else if (ci == 'hess') {
    lower <- 'lowerH'
    upper <- 'upperH'
  }
  
  # Bounds
  if(var %in% c('beta4', 'beta5', 'beta6')){
    Dfss <- Dfp %>% 
      dplyr::filter(type %in% c('beta4', 'beta5', 'beta6'))
    ll <- Dfss[[lower]]
    uu <- Dfss[[upper]]
    bnds <- range(c(ll, uu))
  } else bnds <- NULL
  
  # Check if present
  if(sum(Dfp$type == var) == 0) {
    message(paste0('The parameter ', var, ' is not present in "Dfp"'))
    return(NULL)
  }
  
  options(warn = -1)
  p1 <- ggplot(Dfp %>% dplyr::filter(type == !!var)) + 
    theme_bw(base_size = 15) + 
    geom_line(aes(x = !!sym(xvar), y = par, group = !!group, col = !!group),
              show.legend = if(is.na(group)) FALSE else TRUE) + 
    geom_ribbon(aes(x = !!sym(xvar), ymin = !!sym(lower), ymax = !!sym(upper), 
                    group = !!group, fill = !!group), 
                show.legend = if(is.na(group)) FALSE else TRUE,  alpha = 0.1) + 
    scale_fill_manual(name = '', values = cls) + 
    scale_color_manual(name = '', values = cls) + 
    xlab(xlab) + ylab(ylab) +
    theme(axis.title.y=element_text(size=20), legend.position = 'bottom') + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) + 
    coord_cartesian(ylim = bnds)
  options(warn = -1)
  
  p1
}

# Hessian
p.alpha.hess  <- plot_parameter(Dfp = Df.par, var = 'alpha',  ci = 'hess')
p.beta2.hess  <- plot_parameter(Dfp = Df.par, var = 'beta2',  ci = 'hess')
p.kappa2.hess <- plot_parameter(Dfp = Df.par, var = 'kappa2', ci = 'hess')
p.beta3.hess  <- plot_parameter(Dfp = Df.par, var = 'beta3',  ci = 'hess')
p.kappa3.hess <- plot_parameter(Dfp = Df.par, var = 'kappa3', ci = 'hess')
p.beta4.hess  <- plot_parameter(Dfp = Df.par, var = 'beta4',  ci = 'hess')
p.beta5.hess  <- plot_parameter(Dfp = Df.par, var = 'beta5',  ci = 'hess')
p.beta6.hess  <- plot_parameter(Dfp = Df.par, var = 'beta6',  ci = 'hess')

p.hess <- ggpubr::ggarrange(
  p.alpha.hess,               
  ggpubr::ggarrange(p.beta2.hess, p.kappa2.hess, ncol = 2, 
                    labels = c("B", "C"), legend = 'none'), 
  ggpubr::ggarrange(p.beta3.hess, p.kappa3.hess, ncol = 2, 
                    labels = c("D", "E"), legend = 'none'), 
  ggpubr::ggarrange(p.beta4.hess, p.beta5.hess, ncol = 2, 
                    labels = c("F", "G"), legend = 'none'), 
  nrow = 4, 
  labels = "A",
  common.legend = TRUE,
  legend = 'bottom'
) + guides(fill=guide_legend(nrow=2,byrow=TRUE))

# ggsave('../Figures updated/param_hess.pdf',
#        plot = p.hess, height = 30, width = 25, units = 'cm')

# Bootstrap
p.alpha.boot  <- plot_parameter(Dfp = Df.par, var = 'alpha',  ci = 'boot')
p.beta2.boot  <- plot_parameter(Dfp = Df.par, var = 'beta2',  ci = 'boot')
p.kappa2.boot <- plot_parameter(Dfp = Df.par, var = 'kappa2', ci = 'boot')
p.beta3.boot  <- plot_parameter(Dfp = Df.par, var = 'beta3',  ci = 'boot')
p.kappa3.boot <- plot_parameter(Dfp = Df.par, var = 'kappa3', ci = 'boot')
p.beta4.boot  <- plot_parameter(Dfp = Df.par, var = 'beta4',  ci = 'boot')
p.beta5.boot  <- plot_parameter(Dfp = Df.par, var = 'beta5',  ci = 'boot')
p.beta6.boot  <- plot_parameter(Dfp = Df.par, var = 'beta6',  ci = 'boot')


p.boot <- ggpubr::ggarrange(
  p.alpha.boot,               
  ggpubr::ggarrange(p.beta2.boot, p.kappa2.boot, ncol = 2, 
                    labels = c("B", "C"), legend = 'none'), 
  ggpubr::ggarrange(p.beta3.boot, p.kappa3.boot, ncol = 2, 
                    labels = c("D", "E"), legend = 'none'), 
  ggpubr::ggarrange(p.beta4.boot, p.beta5.boot, ncol = 2, 
                    labels = c("F", "G"), legend = 'none'), 
  nrow = 4, 
  labels = "A",
  common.legend = TRUE,
  legend = 'bottom' 
) + guides(fill=guide_legend(nrow=2,byrow=TRUE))

p.hess
p.boot

# ggsave('../Figures updated/param_boot.pdf',
#        plot = p.boot, height = 30, width = 25, units = 'cm')

##### 2) DLNM visualisations ----

### Pre-sets

# Presets
Q       <- unlist(lapply(obj$lcd$cb, function(x) dim(x)[2]))
regions <- dimnames(dtxr)[[3]]
ages    <- dimnames(dtxr)[[2]]

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# Temperature and influenza
objRR <- setNames(
  rep(list(setNames(rep(list(NA), length(regions)), regions)), length(Q)),
  names(Q)
)

# Covariance matrices
patterns <- list(Temp = "Beta4|Gamma1", Infl = "Beta5|Gamma2")
sel      <- lapply(patterns, function(p) grepl(p, rownames(I2)))
cov      <- lapply(names(sel), function(n) {
  s <- sel[[n]]
  I2[s, s]
})
names(cov) <- names(patterns)

### Calculate the MMT for every region
mmt_region <- function(model.fit, model.args, var, x.grid, region){
  argsv    <- model.args[[region]][[var]]$argvar
  argsl    <- model.args[[region]][[var]]$arglag
  lagvec   <- 0:argsl$Boundary.knots[2]
  gammav   <- model.fit$lcd$gammal[[var]][, which(region == names(model.args))]
  basislag <- do.call("onebasis", c(list(x = lagvec), argsl))
  
  Fcum <- vapply(x.grid, function(xx){
    basisvar <- do.call("onebasis", c(list(x = rep(xx, length(lagvec))), argsv))
    sum(colSums(mgcv::tensor.prod.model.matrix(list(basisvar, basislag))) * gammav)
  }, numeric(1))
  
  x.grid[which.min(Fcum)]
}

tgr_seq_lst <- lapply(regions, function(r){
  Dfr <- DfM %>% dplyr::filter(region == r)
  seq(ceiling_dec(min(Dfr$tg), 1), floor_dec(max(Dfr$tg), 1), by = 0.1)
})
iar_seq_lst <- lapply(regions, function(r){
  Dfr <- DfM %>% dplyr::filter(region == r)
  seq(ceiling_dec(min(Dfr$ia100), 2), floor_dec(max(Dfr$ia100), 2), by = 0.01)
})
names(tgr_seq_lst) <- names(iar_seq_lst) <- regions

ref_temp <- vapply(regions, function(r)
  mmt_region(obj, args, 'Temp', tgr_seq_lst[[r]], r), numeric(1))
ref_temp

### Functions to calculate relative risk
relative_risk <- function(model.fit, model.args, var, 
                          covmat, x, ref, region) {
  # Vector of lags and covariate points
  at       <- x
  lagvec   <- 0:model.args[[region]][[var]]$arglag$Boundary.knots[2]
  varvec   <- rep(at,length(lagvec))
  
  # Construct onebasis objects
  argsv    <- model.args[[region]][[var]]$argvar
  argsl    <- model.args[[region]][[var]]$arglag
  basisvar <- do.call("onebasis", c(list(x = varvec), argsv))
  basislag <- do.call("onebasis", c(list(x = lagvec), argsl))
  basiscen <- do.call("onebasis", c(list(x = ref), argsv))
  basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
  
  # Construct crossbasis as tensor product
  Zlq <- mgcv::tensor.prod.model.matrix(list(basisvar, basislag))
  
  # Region index and ages
  regions <- names(model.args)
  ages    <- colnames(model.fit$lcd$mhat)
  reg.j   <- which(region == regions)
  
  # Parameter estimates
  betav  <- model.fit$lcd$betal[[var]]
  gammav <- model.fit$lcd$gammal[[var]]
  
  # RR per lag
  RR1 <- as.vector(exp(Zlq %*% gammav[,reg.j]))
  RR1 <- unlist(lapply(1:length(ages), function(j) RR1^(betav[j])))
  
  # Overall RR
  RR2 <- as.vector(exp(colSums(Zlq) %*% gammav[,reg.j]))^(betav)
  
  # RR per lag - Bootstrap
  RR1.boot <- sapply(1:length(btci), function(b) {
    v <- as.vector(exp(Zlq %*% btci[[b]]$gammal[[var]][,reg.j]))
    unlist(lapply(1:length(ages), function(j) 
      v^(btci[[b]]$betal[[var]][j])))
  }) %>%
    matrixStats::rowQuantiles(probs = c(0.025,0.975))
  colnames(RR1.boot) <- c('lowerB', 'upperB')
  
  ## Overall RR - Bootstrap
  RR2.boot <- sapply(1:length(btci), function(b){
    vec <- as.vector(exp(colSums(Zlq) %*% btci[[b]]$gammal[[var]][,reg.j]))
    vec^(btci[[b]]$betal[[var]])
  }) %>%
    matrixStats::rowQuantiles(probs = c(0.025,0.975)) 
  colnames(RR2.boot) <- c('lowerB', 'upperB')
  
  # SD RR per lag - Hessian
  SD1.hess <- t(sapply(lagvec, function(l) {
    sqrt((betav * RR1[l+1])^2 * as.numeric((t(Zlq[l+1,]) %*% 
                                              covmat[[var]] %*% Zlq[l+1,])))})) %>%
    as.vector()
  RR1.hess <- cbind('lowerH' = RR1 - 1.96 * SD1.hess,
                    'upperH' = RR1 + 1.96 * SD1.hess)
  
  # SD Overall RR
  SD2.hess <- as.vector(sqrt((betav * RR2)^2 * as.numeric(
    (t(colSums(Zlq)) %*% covmat[[var]] %*% colSums(Zlq)))))
  RR2.hess <- cbind('lowerH'  = RR2 - 1.96 * SD2.hess, 
                    'upperH'  = RR2 + 1.96 * SD2.hess)
  
  # Combine
  RR1.ci <- bind_cols('value'  = x, 
                      'lag'    = rep(lagvec, times = length(ages)), 
                      'region' = region,
                      'age'    = rep(ages, each = length(lagvec)),
                      'RR.lag' = RR1, 
                      RR1.boot, RR1.hess)
  
  RR2.ci <- bind_cols('value'  = x, 
                      'region' = region, 
                      'age'    = ages,
                      'RR'     = RR2, 
                      RR2.boot, RR2.hess)
  
  list('RR' = RR2.ci, 'RR.lag' = RR1.ci)
}

### Calculate the relative risk for each region and age class
# for(r in regions){
#   message(paste0('Region ', r))
# 
#   # Covariance matrices
#   selr <- lapply(1:length(cov), function(j) {
#     rn <- rownames(cov[[j]])
#     grepl(paste0('_Region', r, '$'), rn)
#   })
#   covr <- lapply(1:length(selr), function(j)
#     cov[[j]][selr[[j]],selr[[j]]])
#   names(covr) <- names(patterns)
# 
#   # Sequences to evaluate RR at
#   tgr_seq     <- tgr_seq_lst[[r]]
#   iar_seq     <- iar_seq_lst[[r]]
# 
#   # Xpred for temperature
#   Xpred_tg <- if('Temp' %in% names(Q)) lapply(tgr_seq, function(x)
#     relative_risk(model.fit = obj, model.args = args, var = 'Temp',
#                   covmat = covr, x = x, ref = ref_temp[[r]], region = r))
# 
#   # Xpred for influenza
#   Xpred_ia <- if('Infl' %in% names(Q)) lapply(iar_seq, function(x)
#     relative_risk(model.fit = obj, model.args = args, var = 'Infl',
#                   covmat = covr, x = x, ref = 0, region = r))
# 
#   objRR[['Temp']][[r]] <- Xpred_tg
#   objRR[['Infl']][[r]] <- Xpred_ia
# }

# saveRDS(object = objRR, file = '../Results updated/relative_risk.rds')
objRR <- readRDS('../Results updated/relative_risk.rds')

tg_overall <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Temp[[r]], `[[`, "RR"))))
ia_overall <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Infl[[r]], `[[`, "RR"))))

### Overall RR
reg_vec <- regions[c(5,8,11)]
age_vec <- ages[c(1,5,9)]
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

plots_overall <- rep(list(), 6)
for(r in 1:length(reg_vec)){
  sub1 <- tg_overall %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  sub2 <- ia_overall %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  
  plots_overall[[r]] <- 
    ggplot(sub1 %>% dplyr::filter(region == reg_vec[r])) +
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Average temperature','')) + 
    ylab(ifelse(r==1,'Overall relative risk','')) + 
    geom_ribbon(aes(x = value, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = value, y = RR, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub1$RR, sub1$lowerB, sub1$upperB))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  plots_overall[[3 + r]] <- 
    ggplot(sub2 %>% dplyr::filter(region == reg_vec[r])) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'ILI exceedance','')) + 
    ylab(ifelse(r==1,'Overall relative risk','')) + 
    geom_ribbon(aes(x = value, ymin = lowerB, ymax = upperB, 
                    group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = value, y = RR, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub2$RR, sub2$lowerB, sub2$upperB))) 
  
  if(r == 1 | r == 3){
    plots_overall[[r]] <- plots_overall[[r]] + 
      theme(axis.title.x=element_blank())
    plots_overall[[3+r]] <- plots_overall[[3+r]] + 
      theme(axis.title.x=element_blank())
  }
  
  if(r == 2 | r == 3 | r == 1){
    plots_overall[[r]] <- plots_overall[[r]] + 
      theme(axis.title.y=element_blank())
    plots_overall[[3+r]] <- plots_overall[[3+r]] + 
      theme(axis.title.y=element_blank())
  }
}

fig1 <- ggpubr::ggarrange(plotlist = plots_overall, common.legend = TRUE, 
                          legend = 'right', nrow = 2, ncol = 3, align = 'h',
                          widths = c(1,1,1)) 

fig1 <- ggpubr::annotate_figure(fig1, 
                                left = ggpubr::text_grob("Overall relative risk", 
                                                         rot = 90, vjust = 0.5, 
                                                         hjust = 0.4, size = 17))
fig1

# ggsave('../Figures updated/RRoverall.pdf',
#        plot = fig1, height = 18, width = 30, units = 'cm')

## All regions
reg_vec <- regions
age_vec <- ages[c(1,5,9)]
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

plots_tg <- rep(list(), 12)
plots_ia <- rep(list(), 12)

for(r in 1:length(reg_vec)){
  sub1 <- tg_overall %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  sub2 <- ia_overall %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  
  plots_tg[[r]] <- 
    ggplot(sub1 %>% dplyr::filter(region == reg_vec[r])) +
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r%%3==2,'Average temperature','')) + 
    ylab(ifelse(r%%3==1,'Overall relative risk','')) + 
    geom_ribbon(aes(x = value, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = value, y = RR, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = c(0.8,5)) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  plots_ia[[r]] <- 
    ggplot(sub2 %>% dplyr::filter(region == reg_vec[r])) + 
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r%%3==2,'ILI exceedance','')) + 
    ylab(ifelse(r%%3==1,'Overall relative risk','')) + 
    geom_ribbon(aes(x = value, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = value, y = RR, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = c(0.8,2.5)) 
  
}

fig_tg <- ggpubr::ggarrange(plotlist = plots_tg, common.legend = TRUE, 
                            legend = 'bottom', nrow = 4, ncol = 3, align = 'h') 

fig_ia <- ggpubr::ggarrange(plotlist = plots_ia, common.legend = TRUE, 
                            legend = 'bottom', nrow = 4, ncol = 3, align = 'h') 

# ggsave('../Figures updated/RRoverall_tg_all.pdf',
#        plot = fig_tg, height = 50, width = 36, units = 'cm')
# ggsave('../Figures updated/RRoverall_ia_all.pdf',
#        plot = fig_ia, height = 50, width = 36, units = 'cm')

### Overview table
reg_vec <- regions
age_vec <- ages
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

# Quantiles
quants <- round(DfM %>% 
                  dplyr::group_by(region) %>% 
                  dplyr::reframe(q = quantile(tg, 0.995)) %>% 
                  dplyr::pull(q), 1)
names(quants) <- regions


xtable(tg_overall %>% 
         dplyr::group_by(region) %>%
         dplyr::mutate(QQ = quants[region[1]]) %>%
         dplyr::filter(round(value,1) == quants[region[1]]) %>%  
         ungroup() %>%
         dplyr::select(region, age, RR) %>%
         pivot_wider(
           id_cols = c(age),
           names_from = region,
           values_from = RR,
           names_prefix = "region_"
         ) %>%
         column_to_rownames('age'), 
       booktabs = TRUE, 
       include.rownames = FALSE,
       sanitize.text.function = identity)

tavg %>%
  dplyr::group_by(region) %>%
  dplyr::reframe(q = quantile(tg, 0.995)) %>%
  pull(q) %>%
  round(1) %>%
  paste(collapse = ') & (')

### RR at specific value over lags (99.5% quantile)

reg_vec <- regions[c(5,8,11)]
age_vec <- ages[c(1,5,9)]
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

# Quantiles
quants_tg <- round(DfM %>% 
                     dplyr::group_by(region) %>% 
                     dplyr::reframe(q = quantile(tg, 0.995)) %>% 
                     dplyr::pull(q), 1)
quants_ia <- round(DfM %>% 
                     dplyr::group_by(region) %>% 
                     dplyr::reframe(q = quantile(ia100, 0.995)) %>% 
                     dplyr::pull(q), 2)
names(quants_tg) <- regions
names(quants_ia) <- regions

# Data
tg_lag <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Temp[[r]], `[[`, "RR.lag")))) %>%
  dplyr::filter(round(value,1) == round(quants_tg[region],1))
ia_lag <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Infl[[r]], `[[`, "RR.lag")))) %>%
  dplyr::filter(round(value,2) == round(quants_ia[region],2))

plots_lag <- rep(list(), 6)
for(r in 1:length(reg_vec)){
  sub1 <- tg_lag %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  sub2 <- ia_lag %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  
  plots_lag[[r]] <- 
    ggplot(sub1 %>% dplyr::filter(region == reg_vec[r])) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Lag','')) + 
    ylab(ifelse(r==1,'RR','')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    geom_ribbon(aes(x = lag, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = lag, y = RR.lag, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub1$RR.lag, sub1$lowerB, sub1$upperB))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  plots_lag[[3 + r]] <- 
    ggplot(sub2 %>% dplyr::filter(region == reg_vec[r])) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Lag','')) + 
    ylab(ifelse(r==1,'RR','')) + 
    geom_ribbon(aes(x = lag, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = lag, y = RR.lag, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub2$RR.lag, sub2$lowerB, sub2$upperB))) + 
    theme(plot.title = element_text(size=12))
  
  
  if(r == 1 | r == 3 | r == 2){
    plots_lag[[r]] <- plots_lag[[r]] + 
      theme(axis.title.x=element_blank())
    plots_lag[[3+r]] <- plots_lag[[3+r]] + 
      theme(axis.title.x=element_blank())
  }
  
  if(r == 2 | r == 3 | r == 1){
    plots_lag[[r]] <- plots_lag[[r]] + 
      theme(axis.title.y=element_blank())
    plots_lag[[3+r]] <- plots_lag[[3+r]] + 
      theme(axis.title.y=element_blank())
  }
}

fig2 <- ggpubr::ggarrange(plotlist = plots_lag, common.legend = TRUE, 
                          legend = 'right', nrow = 2, ncol = 3, align = 'h',
                          widths = c(1,1,1)) + xlab('sdfl') 

fig2 <- ggpubr::annotate_figure(fig2, 
                                left = ggpubr::text_grob("Relative risk by lag", rot = 90, 
                                                         vjust = 0.5, hjust = 0.4, size = 17),
                                bottom = ggpubr::text_grob("Lag", size = 17, hjust = 1.4))
fig2
# ggsave('../Figures updated/RRlag.pdf',
#        plot = fig2, height = 18, width = 30, units = 'cm')

## All regions

reg_vec <- regions
age_vec <- ages[c(1,5,9)]
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

plots_lag_tg = plots_lag_ia <- rep(list(), 12)

for(r in 1:length(reg_vec)){
  sub1 <- tg_lag %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  sub2 <- ia_lag %>% 
    dplyr::filter(region %in% reg_vec,
                  age %in% age_vec)
  
  plots_lag_tg[[r]] <- 
    ggplot(sub1 %>% dplyr::filter(region == reg_vec[r])) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r%%3==2,'Lag','')) + 
    ylab(ifelse(r%%3==1,'RR','')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    geom_ribbon(aes(x = lag, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = lag, y = RR.lag, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub1$RR.lag, sub1$lowerB, sub1$upperB))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"),
          axis.title.y=element_blank())
  
  plots_lag_ia[[r]] <- 
    ggplot(sub2 %>% dplyr::filter(region == reg_vec[r])) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r%%3==2,'Lag','')) + 
    ylab(ifelse(r%%3==1,'RR','')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    geom_ribbon(aes(x = lag, ymin = lowerB, ymax = upperB, group = age, fill = age),  
                col = 'gray', alpha = 0.2) + 
    geom_line(aes(x = lag, y = RR.lag, group = age, col = age), linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = range(c(sub2$RR.lag, sub2$lowerB, sub2$upperB))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"),
          axis.title.y=element_blank())
}

fig_lag_tg <- ggpubr::ggarrange(plotlist = plots_lag_tg, common.legend = TRUE,
                                legend = 'right', nrow = 4, ncol = 3, align = 'h') 
fig_lag_ia <- ggpubr::ggarrange(plotlist = plots_lag_ia, common.legend = TRUE,
                                legend = 'right', nrow = 4, ncol = 3, align = 'h') 

# ggsave('../Figures updated/RRlag_tg_all.pdf',
#        plot = fig_lag_tg, height = 50, width = 36, units = 'cm')
# ggsave('../Figures updated/RRlag_ia_all.pdf',
#        plot = fig_lag_ia, height = 50, width = 36, units = 'cm')

### Surface plot
reg_vec <- regions[c(5,8,11)]
age_vec <- ages[c(1,5,9)]
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

tg_lag <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Temp[[r]], `[[`, "RR.lag")))) %>%
  dplyr::filter(region %in% reg_vec)
ia_lag <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Infl[[r]], `[[`, "RR.lag"))))  %>%
  dplyr::filter(region %in% reg_vec)

plots_surf <- rep(list(), 6)
for(r in 1:length(reg_vec)){
  
  df1 <- tg_lag %>% 
    dplyr::filter(region == reg_vec[r],
                  age    == age_vec[3])
  df2 <- ia_lag %>% 
    dplyr::filter(region == reg_vec[r],
                  age    == age_vec[3])
  
  plots_surf[[r]] <-  
    ggplot(df1, aes(x = value, y = lag, z = RR.lag)) + 
    theme_bw(base_size = 15) + 
    geom_contour_fill(bins = 30) +
    scale_fill_gradientn(
      name = '',
      colours = c("#33A02C", "white", "white", "#E31A1C", "#E31A1C"),
      values = scales::rescale(c(min(tg_lag$RR.lag), 0.99,1.01, 
                                 max(tg_lag$RR.lag) - 0.25, max(tg_lag$RR.lag))),
      breaks = pretty(range(tg_lag$RR.lag), 9),
      limits = range(tg_lag$RR.lag)
    ) +
    xlab(ifelse(r==2,'Average temperature','')) + 
    ylab(ifelse(r==1,'Relative risk surface','')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  plots_surf[[3 + r]] <- 
    ggplot(df2, aes(x = value, y = lag, z = RR.lag)) + 
    theme_bw(base_size = 15) + 
    geom_contour_fill(bins = 30) +
    scale_fill_gradientn(
      name = '',
      colours = c("#33A02C", "white", "#E31A1C"),
      values = scales::rescale(c(min(ia_lag$RR.lag), 1, max(ia_lag$RR.lag))),
      breaks = pretty(range(ia_lag$RR.lag), 8),
      limits = range(ia_lag$RR.lag)
    ) +
    xlab(ifelse(r==2,'ILI exceedance','')) + 
    ylab(ifelse(r==1,'Relative risk surface','')) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1, "cm"))
  
  
  if(r == 1 | r == 3 ){
    plots_surf[[r]] <- plots_surf[[r]] + 
      theme(axis.title.x=element_blank())
    plots_surf[[3+r]] <- plots_surf[[3+r]] + 
      theme(axis.title.x=element_blank())
  }
  
  if(r == 2 | r == 3 | r == 1){
    plots_surf[[r]] <- plots_surf[[r]] + 
      theme(axis.title.y=element_blank())
    plots_surf[[3+r]] <- plots_surf[[3+r]] + 
      theme(axis.title.y=element_blank())
  }
}


fig34 <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plots_surf[1:3], 
                                             common.legend = TRUE, align = 'h', 
                                             legend = 'right', nrow = 1),
                           ggpubr::ggarrange(plotlist = plots_surf[4:6], 
                                             common.legend = TRUE, align = 'h',
                                             legend = 'right', nrow = 1), 
                           nrow = 2, heights = c(1.15,1))

fig34 <- ggpubr::annotate_figure(fig34, 
                                 left = ggpubr::text_grob("Lag", rot = 90, 
                                                          vjust = 0.5, hjust = 0.4, 
                                                          size = 17))
fig34
# ggsave('../Figures updated/RRsurface.pdf',
#        plot = fig34, height = 18, width = 30, units = 'cm')

## All regions
reg_vec <- regions
age_vec <- ages
com_vec <- expand.grid(reg = reg_vec, age = age_vec)

tg_surf <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Temp[[r]], `[[`, "RR.lag")))) %>%
  dplyr::filter(region %in% reg_vec)
ia_surf <- bind_rows(lapply(regions, function(r) 
  bind_rows(lapply(objRR$Infl[[r]], `[[`, "RR.lag"))))  %>%
  dplyr::filter(region %in% reg_vec)


for(a in 1:length(age_vec)){
  plots_surf_tg = plots_surf_ia <- rep(list(), 12)
  
  for(r in 1:length(reg_vec)){
    
    df1 <- tg_surf %>% 
      dplyr::filter(region == reg_vec[r],
                    age    == age_vec[a])
    df2 <- ia_surf %>% 
      dplyr::filter(region == reg_vec[r],
                    age    == age_vec[a])
    
    dff1 <- tg_surf %>% 
      dplyr::filter(age == age_vec[a])
    dff2 <- ia_surf %>% 
      dplyr::filter(age == age_vec[a])
    
    plots_surf_tg[[r]] <-  
      ggplot(df1, aes(x = value, y = lag, z = RR.lag)) + 
      theme_bw(base_size = 15) + 
      geom_contour_fill(bins = 30) +
      scale_fill_gradientn(
        name = '',
        colours = c("#33A02C", "white", "#E31A1C", "#E31A1C"),
        values = scales::rescale(c(min(dff1$RR.lag), 1, max(dff1$RR.lag)-0.25, max(dff1$RR.lag))),
        breaks = pretty(range(dff1$RR.lag), 8),
        limits = range(dff1$RR.lag)
      ) +
      xlab(ifelse(r%%3==2,'Average temperature','')) + 
      ylab(ifelse(r%%3==1,'Relative risk surface','')) + 
      ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
      theme(legend.title = element_text(size = 20), 
            legend.text  = element_text(size = 20),
            legend.spacing.x = unit(10, 'cm'),
            legend.key.size = unit(1.25, "cm"),
            axis.title.y=element_blank())
    
    plots_surf_ia[[r]] <- 
      ggplot(df2, aes(x = value, y = lag, z = RR.lag)) + 
      theme_bw(base_size = 15) + 
      geom_contour_fill(bins = 30) +
      scale_fill_gradientn(
        name = '',
        colours = c("#33A02C", "white", "#E31A1C"),
        values = scales::rescale(c(min(dff2$RR.lag), 1, max(dff2$RR.lag))),
        breaks = pretty(range(dff2$RR.lag), 8),
        limits = range(dff2$RR.lag)
      ) +
      xlab(ifelse(r%%3==2,'ILI exceedance','')) + 
      ylab(ifelse(r%%3==1,'Relative risk surface','')) + 
      ggtitle(shapef$NUTS_NAME[shapef$region ==  reg_vec[r]]) + 
      theme(legend.title = element_text(size = 20), 
            legend.text  = element_text(size = 20),
            legend.spacing.x = unit(10, 'cm'),
            legend.key.size = unit(1, "cm"),
            axis.title.y=element_blank())
  }
  
  agegr <- ifelse(a != 10, paste0("Age group ", age_vec[a], '-', 
                                  as.integer(age_vec[a+1]) -1), 'Age group 95+')
  surf_tg_all <- ggpubr::ggarrange(plotlist = plots_surf_tg, 
                                   common.legend = TRUE, align = 'h', 
                                   legend = 'right', nrow = 4, ncol = 3) 
  surf_tg_all <- ggpubr::annotate_figure(surf_tg_all, 
                                         top = ggpubr::text_grob(
                                           agegr, size = 25, vjust = 0,
                                           hjust = 0.5)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  surf_ia_all <- ggpubr::ggarrange(plotlist = plots_surf_ia, 
                                   common.legend = TRUE, align = 'h', 
                                   legend = 'right', nrow = 4, ncol = 3) 
  surf_ia_all <- ggpubr::annotate_figure(surf_ia_all,
                                         top = ggpubr::text_grob(
                                           agegr, size = 25, vjust = 0,
                                           hjust = 0.5)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  
  # ggsave(paste0('../Figures updated/RRsurface_tg_Age', age_vec[a], '.pdf'),
  #        plot = surf_tg_all, height = 50, width = 36, units = 'cm')
  # ggsave(paste0('../Figures updated/RRsurface_ia_Age', age_vec[a], '.pdf'),
  #        plot = surf_ia_all, height = 50, width = 36, units = 'cm')
}

##### 3) In-sample fit - Hessian & Bootstrap ----

# Dimensions dlnm-part
Q <- unlist(lapply(obj$lcd$cb, function(x) dim(x)[2]))

# Store parameters and variance-covariance matrices
param.base <- c(obj$lcd$beta1, obj$lcd$beta2, obj$lcd$kappa2, 
                obj$lcd$beta3, obj$lcd$kappa3)
param.dlnm <- c(obj$lcd$betal[['Temp']], obj$lcd$gammal[['Temp']], 
                obj$lcd$betal[['Infl']], obj$lcd$gammal[['Infl']],
                obj$lcd$betal[['Int']], obj$lcd$gammal[['Int']])
names(param.base) <- colnames(I1)
names(param.dlnm) <- colnames(I2)

# Covariance matrices from Hessians
V.u.base <- solve(t(C.base) %*% (-H.base) %*% C.base)
V.u.dlnm <- solve(t(C.dlnm) %*% (-H.dlnm) %*% C.dlnm)

# List of indicators where each parameter is stored
sel <- list('b1' = which(grepl('Alpha',  rownames(H.base), fixed = TRUE)),
            'b2' = which(grepl('Beta2',  rownames(H.base), fixed = TRUE)),
            'k2' = which(grepl('Kappa2', rownames(H.base), fixed = TRUE)),
            'b3' = which(grepl('Beta3',  rownames(H.base), fixed = TRUE)),
            'k3' = which(grepl('Kappa3', rownames(H.base), fixed = TRUE)),
            'b4' = which(grepl('Beta4',  rownames(H.dlnm), fixed = TRUE)),
            'g1' = which(grepl('Gamma1', rownames(H.dlnm), fixed = TRUE)),
            'b5' = which(grepl('Beta5',  rownames(H.dlnm), fixed = TRUE)),
            'g2' = which(grepl('Gamma2', rownames(H.dlnm), fixed = TRUE)),
            'b6' = which(grepl('Beta6',  rownames(H.dlnm), fixed = TRUE)),
            'g3' = which(grepl('Gamma3', rownames(H.dlnm), fixed = TRUE)))

# Number of simulations
nsim <- 5000

# Do the simulations!
library(pbapply)
library(doParallel)
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl, c('obj', 'param.base', 'param.dlnm', 'C.base', 'V.u.base', 'C.dlnm', 'V.u.dlnm',
                    'sel', 'nsim', 'R', 'isoyv', 'isowv', 'etxr', 'btci'))
clusterEvalQ(cl,expr= {
  library(dplyr)
  library(mvnfast)
})
obj.sim <- pblapply(cl=cl, X = 1:nsim, FUN = function(i) {
  
  # Hessian-based
  set.seed(i)
  u.base    <- mvnfast::rmvn(1, rep(0, ncol(C.base)), V.u.base)
  boot.base <- param.base + drop(C.base %*% t(u.base))
  u.dlnm    <- mvnfast::rmvn(1, rep(0, ncol(C.dlnm)), V.u.dlnm)
  boot.dlnm <- param.dlnm + drop(C.dlnm %*% t(u.dlnm))
  
  b1.hess <- matrix(boot.base[sel$b1], ncol = R)
  b2.hess <- boot.base[sel$b2]
  k2.hess <- matrix(boot.base[sel$k2], ncol = R)
  b3.hess <- boot.base[sel$b3]
  k3.hess <- matrix(boot.base[sel$k3], ncol = R)
  bl.hess <- list('Temp' = if(length(sel$b4) > 0) boot.dlnm[sel$b4],
                  'Infl' = if(length(sel$b5) > 0) boot.dlnm[sel$b5],
                  'Int'  = if(length(sel$b6) > 0) boot.dlnm[sel$b6])
  gl.hess <- list('Temp' = if(length(sel$b4) > 0) matrix(boot.dlnm[sel$g1], ncol = R),
                  'Infl' = if(length(sel$b5) > 0) matrix(boot.dlnm[sel$g2], ncol = R),
                  'Int'  = if(length(sel$b6) > 0) matrix(boot.dlnm[sel$g3], ncol = R))
  
  common <- intersect(names(obj$lcd$cb), names(gl.hess))
  dlnm.hess <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, obj$lcd$cb[common], gl.hess[common])
  
  # Parametric bootstrap based
  b1.boot <- btci[[i]]$beta1
  b2.boot <- btci[[i]]$beta2
  k2.boot <- btci[[i]]$kappa2
  b3.boot <- btci[[i]]$beta3
  k3.boot <- btci[[i]]$kappa3
  bl.boot <- btci[[i]]$betal
  gl.boot <- btci[[i]]$gammal
  
  dlnm.boot <- Map(function(c, g){
    out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
    out[is.na(out)] <- 0
    out
  }, obj$lcd$cb[common], gl.boot[common])
  
  exp.ds <- exp(outer(obj$lcd$phix, obj$lcd$phir, "+"))
  exp.ds <- replicate(nrow(etxr), exp.ds, simplify = "array")  # gives an n x R x T array
  exp.ds <- aperm(exp.ds, c(3, 1, 2))    
  
  mutxr.hess <- exp(sweep(aperm(b2.hess %o% k2.hess[isoyv,] + 
                                  b3.hess %o% k3.hess[isowv,] + 
                                  Reduce(`+`, Map(function(b, d) 
                                    b %o% d, bl.hess[common], dlnm.hess[common])), 
                                c(2,1,3)), c(2,3), b1.hess, '+'))
  mutxr.boot <- exp(sweep(aperm(b2.boot %o% k2.boot[isoyv,] + 
                                  b3.boot %o% k3.boot[isowv,] + 
                                  Reduce(`+`, Map(function(b, d) 
                                    b %o% d, bl.boot[common], dlnm.boot[common])), 
                                c(2,1,3)), c(2,3), b1.boot, '+'))
  list('H1' = mutxr.hess  * etxr, 
       'H2' = array(rnbinom(length(mutxr.hess), mu = mutxr.hess * etxr, size = exp.ds),
                    dim = dim(mutxr.hess)),
       'B1' = mutxr.boot  * etxr, 
       'B2' = array(rnbinom(length(mutxr.boot), mu = mutxr.boot * etxr, size = exp.ds),
                    dim = dim(mutxr.boot)))
})
stopCluster(cl)

# Compute the quantiles
qarray <- function(lst, type, probs = c(0.025, 0.5, 0.975)){
  # Variables to extract
  var1 <- paste0(type, '1')
  var2 <- paste0(type, '2')
  
  # Do the extraction
  SA1  <- abind(lapply(lst, function(x) x[[var1]]), along = 4)
  SA2  <- abind(lapply(lst, function(x) x[[var2]]), along = 4)
  
  # Quantile arrays
  QA1 <- apply(SA1, c(1, 2, 3), collapse::fquantile, probs = probs)
  QA2 <- apply(SA2, c(1, 2, 3), collapse::fquantile, probs = probs)
  QA1 <- aperm(QA1, c(2, 3, 4, 1))
  QA2 <- aperm(QA2, c(2, 3, 4, 1))
  
  list(QA1, QA2)
}

gc()
QA.hess <- qarray(lst = obj.sim, type = 'H')
gc()
QA.boot <- qarray(lst = obj.sim, type = 'B')
gc()

base   <- obj$lcd$mhat*NA
for(j in 1:ncol(dtxr)){
  base[,j,] <- exp(sweep(obj$lcd$beta2[j] * obj$lcd$kappa2[isoyv,] +
                           obj$lcd$beta3[j] * obj$lcd$kappa3[isowv,], 
                         2, obj$lcd$beta1[j,], '+'))
}
df.sim <- data.frame(
  expand.grid('Date'   = dimnames(dtxr)[[1]], 
              'Age'    = dimnames(dtxr)[[2]], 
              'Region' = dimnames(dtxr)[[3]]),
  'deaths'   = as.vector(dtxr),
  'expo'     = as.vector(etxr),
  'rate'     = as.vector(dtxr)/as.vector(etxr),
  'fit'      = as.vector(obj$lcd$mhat),
  'base'     = as.vector(base),
  'Q11.hess' = as.vector(QA.hess[[1]][,,,1])/as.vector(etxr),
  'Q21.hess' = as.vector(QA.hess[[1]][,,,2])/as.vector(etxr),
  'Q31.hess' = as.vector(QA.hess[[1]][,,,3])/as.vector(etxr),
  'Q12.hess' = as.vector(QA.hess[[2]][,,,1])/as.vector(etxr),
  'Q22.hess' = as.vector(QA.hess[[2]][,,,2])/as.vector(etxr),
  'Q32.hess' = as.vector(QA.hess[[2]][,,,3])/as.vector(etxr),
  'Q11.boot' = as.vector(QA.boot[[1]][,,,1])/as.vector(etxr),
  'Q21.boot' = as.vector(QA.boot[[1]][,,,2])/as.vector(etxr),
  'Q31.boot' = as.vector(QA.boot[[1]][,,,3])/as.vector(etxr),
  'Q12.boot' = as.vector(QA.boot[[2]][,,,1])/as.vector(etxr),
  'Q22.boot' = as.vector(QA.boot[[2]][,,,2])/as.vector(etxr),
  'Q32.boot' = as.vector(QA.boot[[2]][,,,3])/as.vector(etxr)
)

## Plot for baseline and full force of mortality
reg_vec <- dimnames(dtxr)[[3]][c(5,8,11)]
age_vec <- colnames(dtxr)[c(5,9)]

pvis1 <- list()
for(a in 1:length(age_vec)){
  dfs1 <- df.sim %>% dplyr::filter(Region == reg_vec[1],
                                   Age == age_vec[a])
  pvis1[[3*(a-1) + 1]] <- ggplot(dfs1) + 
    theme_bw(base_size = 15) + 
    xlab('') + 
    ylab(bquote(log(mu[.(age_vec[a])*',t,w,r']))) + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'obs')) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = 'fit')) +
    geom_line(aes(x = as.Date(Date), y = log(base), col = 'base')) + 
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[1]], '')) + 
    scale_color_manual(name = '', 
                       values = c('obs' = "gray70", 'base' = "#1F78B4", 'fit' = "#E31A1C"),
                       labels = c('obs' = 'Observed', 'fit' = 'Full model', 'base' = 'Baseline'),
                       breaks = c('obs', 'base', 'fit')) + 
    coord_cartesian(ylim = quantile(df.sim %>% 
                                      dplyr::filter(Age == age_vec[a], 
                                                    Region %in% reg_vec) %>% 
                                      pull(rate) %>% log(), c(0.005,1))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  dfs2 <- df.sim %>% dplyr::filter(Region == reg_vec[2],
                                   Age == age_vec[a])
  pvis1[[3*(a-1) + 2]] <- ggplot(dfs2) + 
    theme_bw(base_size = 15) + 
    xlab('Date') + 
    ylab('') + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'obs')) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = 'fit')) +
    geom_line(aes(x = as.Date(Date), y = log(base), col = 'base')) + 
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[2]], '')) + 
    scale_color_manual(name = '', 
                       values = c('obs' = "gray70", 'base' = "#1F78B4", 'fit' = "#E31A1C"),
                       labels = c('obs' = 'Observed', 'fit' = 'Full model', 'base' = 'Baseline'),
                       breaks = c('obs', 'base', 'fit')) + 
    coord_cartesian(ylim = quantile(df.sim %>% 
                                      dplyr::filter(Age == age_vec[a], 
                                                    Region %in% reg_vec) %>% 
                                      pull(rate) %>% log(), c(0.005,1))) + 
    theme(legend.title = element_text(size = 20), 
          legend.text  = element_text(size = 20),
          legend.spacing.x = unit(10, 'cm'),
          legend.key.size = unit(1.25, "cm"))
  
  
  dfs3 <- df.sim %>% dplyr::filter(Region == reg_vec[3],
                                   Age == age_vec[a])
  pvis1[[3*(a-1) + 3]] <- ggplot(dfs3) + 
    theme_bw(base_size = 15) + 
    xlab('') + 
    ylab('') + 
    geom_line(aes(x = as.Date(Date), y = log(rate), col = 'obs')) +
    geom_line(aes(x = as.Date(Date), y = log(fit), col = 'fit')) +
    geom_line(aes(x = as.Date(Date), y = log(base), col = 'base')) + 
    ggtitle(ifelse(a==1, shapef$NUTS_NAME[shapef$region == reg_vec[3]], '')) + 
    scale_color_manual(name = '', 
                       values = c('obs' = "gray70", 'base' = "#1F78B4", 'fit' = "#E31A1C"),
                       labels = c('obs' = 'Observed', 'fit' = 'Full model', 'base' = 'Baseline'),
                       breaks = c('obs', 'base', 'fit')) + 
    coord_cartesian(ylim = quantile(df.sim %>% 
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
# ggsave('../Figures updated/InSampleFit.pdf',
#        plot = plot.in.sample, height = 22, width = 36, units = 'cm')

## Plot all regions
reg_vec <- dimnames(dtxr)[[3]]
age_vec <- colnames(dtxr)

for(a in 1:length(age_vec)){
  pvis1 <- list()
  for(r in 1:length(reg_vec)){
    dfs1 <- df.sim %>% dplyr::filter(Region == reg_vec[r],
                                     Age == age_vec[a])
    
    pvis1[[r]] <- ggplot(dfs1) + 
      theme_bw(base_size = 15) + 
      xlab(ifelse(r == 7, 'Date', '')) + 
      ylab(bquote(log(mu[.(age_vec[a])*',t,w,r']))) + 
      geom_line(aes(x = as.Date(Date), y = log(rate), col = 'obs')) +
      geom_line(aes(x = as.Date(Date), y = log(fit), col = 'fit')) +
      geom_line(aes(x = as.Date(Date), y = log(base), col = 'base')) + 
      ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) + 
      scale_color_manual(name = '', 
                         values = c('obs' = "gray70", 'base' = "#1F78B4", 'fit' = "#E31A1C"),
                         labels = c('obs' = 'Observed', 'fit' = 'Full model', 'base' = 'Baseline'),
                         breaks = c('obs', 'base', 'fit')) + 
      coord_cartesian(ylim = quantile(df.sim %>% 
                                        dplyr::filter(Age == age_vec[a], 
                                                      Region %in% reg_vec) %>% 
                                        pull(rate) %>% log(), c(0.005,1))) + 
      theme(legend.title = element_text(size = 20), 
            legend.text  = element_text(size = 20),
            legend.spacing.x = unit(10, 'cm'),
            legend.key.size = unit(1.25, "cm"))
    
    if(r %% 3 != 1) pvis1[[r]] <- pvis1[[r]] + ylab('') + 
      theme(axis.title.y=element_blank())
  }
  plot.in.sample <- ggpubr::ggarrange(plotlist = pvis1, common.legend = TRUE, 
                                      legend = 'bottom', nrow = 4, ncol = 3, 
                                      align = 'h') 
  
  agegr <- ifelse(a != 10, paste0("Age group ", age_vec[a], '-', 
                                  as.integer(age_vec[a+1]) -1), 'Age group 95+')
  plot.in.sample <- ggpubr::annotate_figure(plot.in.sample, 
                                            top = ggpubr::text_grob(
                                              agegr, size = 25, vjust = 0,
                                              hjust = 0.35)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  
  ggsave(paste0('../Figures updated/InSampleFit_Age', age_vec[a], '.pdf'),
         plot = plot.in.sample, height = 50, width = 36, units = 'cm')
}

### With uncertainty
reg_vec <- dimnames(dtxr)[[3]][c(5,8,11)]
age_vec <- colnames(dtxr)[c(1,5,9)]

pvis2 <- list()

for(r in 1:length(reg_vec)){
  dfs <- df.sim %>% dplyr::filter(Region == reg_vec[r],
                                  Age %in% age_vec)
  
  pvis2[[r]] <- ggplot(dfs) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r==2,'Date','')) + 
    ylab(if(r==1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(aes(x = as.Date(Date), y = log(rate), group = Age), col = 'gray70')  +
    geom_ribbon(aes(x = as.Date(Date), ymin = log(Q12.boot), ymax = log(Q32.boot),
                    group = Age, fill = Age), alpha = 0.2) + 
    geom_line(aes(x = as.Date(Date), y = log(Q32.boot), col = Age)) +
    geom_line(aes(x = as.Date(Date), y = log(Q12.boot), col = Age)) + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = c(-7.5,-0.8))
}

plot.uncertainty <- ggpubr::ggarrange(plotlist = pvis2, common.legend = TRUE, 
                                      legend = 'bottom', nrow = 1, ncol = 3, 
                                      align = 'h', widths = c(1,1,1)) 
plot.uncertainty

# ggsave('../Figures updated/Uncertainty.pdf',
#        plot = plot.uncertainty, height = 12, width = 30, units = 'cm')

## Plot all regions
reg_vec <- dimnames(dtxr)[[3]]
age_vec <- colnames(dtxr)[c(1,5,9)]

for(r in 1:length(reg_vec)){
  dfs <- df.sim %>% dplyr::filter(Region == reg_vec[r],
                                  Age %in% age_vec)
  
  pvis2[[r]] <- ggplot(dfs) + 
    theme_bw(base_size = 15) + 
    xlab(ifelse(r%%3==2,'Date','')) + 
    ylab(if(r%%3==1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(aes(x = as.Date(Date), y = log(rate), group = Age), col = 'gray70')  +
    geom_ribbon(aes(x = as.Date(Date), ymin = log(Q12.boot), ymax = log(Q32.boot),
                    group = Age, fill = Age), alpha = 0.2) + 
    geom_line(aes(x = as.Date(Date), y = log(Q32.boot), col = Age)) +
    geom_line(aes(x = as.Date(Date), y = log(Q12.boot), col = Age)) + 
    scale_fill_manual(name = 'Age', 
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    scale_color_manual(name = 'Age', 
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) + 
    coord_cartesian(ylim = c(-7.5,-0.5))
}

plot.uncertainty <- ggpubr::ggarrange(plotlist = pvis2, common.legend = TRUE, 
                                      legend = 'bottom', nrow = 4, ncol = 3, 
                                      align = 'h') 
plot.uncertainty

# ggsave('../Figures updated/Uncertainty_all.pdf',
#        plot = plot.uncertainty, height = 50, width = 36, units = 'cm')

##### 4) Check residuals ----

# Discrete scale colour gradient
cuts <- qchisq(c(0.75, 0.95, 0.995, 0.9995, 0.99995, 0.999995), df = 1)
vec  <- c(0, cuts, max(obj$lcd$epsilon^2))
vec <- (vec[-length(vec)] + vec[-1]) / 2
cols <- c("#E6F4E4","#B2DF8A", "#E31A1C", "#C31A1C", "#8F140F","#4F140F", "black")

plot.rho <- list()
reg_vec  <- dimnames(dtxr)[[3]][c(5,8,11)]

for(r in 1:length(reg_vec)){
  eps_long <- as.data.frame(obj$lcd$epsilon[,,reg_vec[r]])
  brks <- c(0, cuts)
  
  plot.rho[[r]] <- ggplot(eps_long) + 
    xlab(ifelse(r == 2, 'Date', '')) + 
    ylab(ifelse(r == 1, 'Age', '')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    geom_tile(aes(x = as.Date(ISODate), y = age, fill = Freq^2)) +
    scale_fill_stepsn(
      name    = '',
      breaks  = brks, 
      limits  = c(0, max(eps_long$Freq^2)),
      values  = scales::rescale(vec, from = c(0, max(eps_long$Freq^2))),
      colours = cols,
      labels  = function(x) sprintf("%.3f", x)) +
    theme(legend.position = 'bottom',
          legend.key.width = unit(3, 'cm'))
  
  if(r == 2 | r == 3| r == 1){
    plot.rho[[r]] <- plot.rho[[r]] + 
      theme(axis.title.y=element_blank())
  }
}

pearson <- ggpubr::ggarrange(plotlist = plot.rho, ncol = 3, nrow = 1,
                             align = 'h', common.legend = TRUE, legend = 'bottom')
pearson <- ggpubr::annotate_figure(pearson, 
                                   left = ggpubr::text_grob("Age group", rot = 90,
                                                            vjust = 0.5, hjust = 0.2, 
                                                            size = 15))
# ggsave('../Figures updated/Pearson.pdf',
#        plot = pearson, height = 12, width = 30, units = 'cm')

prop.table(table(obj$lcd$epsilon[,,reg_vec[1]]^2 > 3.841))['TRUE']
prop.table(table(obj$lcd$epsilon[,,reg_vec[2]]^2 > 3.841))['TRUE']
prop.table(table(obj$lcd$epsilon[,,reg_vec[3]]^2 > 3.841))['TRUE']
round(range(sapply(regions, function(r) 
  prop.table(table(obj$lcd$epsilon[,,r]^2 > 3.841))['TRUE'])), 4)

sum(obj$lcd$epsilon^2)
qchisq(0.95, df = length(dtxr) - (ncol(C.base) + edf.dlnm + 
                                    length(obj$lcd$phix) + length(obj$lcd$phir) - 1))

## All regions
plot.rho <- list()
reg_vec  <- dimnames(dtxr)[[3]]

for(r in 1:length(reg_vec)){
  eps_long <- as.data.frame(obj$lcd$epsilon[,,reg_vec[r]])
  
  plot.rho[[r]] <- ggplot(eps_long) + 
    xlab(ifelse(r%%3==2, 'Date', '')) + 
    ylab(ifelse(r%%3==1, 'Age', '')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    geom_tile(aes(x = as.Date(ISODate), y = age, fill = Freq^2)) +
    scale_fill_stepsn(
      name    = '',
      breaks  = c(0,cuts), 
      limits  = c(0, max(obj$lcd$epsilon^2)),
      values  = scales::rescale(vec, from = c(0, max(obj$lcd$epsilon^2))),
      colours = cols,
      labels  = function(x) sprintf("%.3f", x)) +
    theme(legend.position = 'bottom',
          legend.key.width = unit(3, 'cm'),
          axis.title.y=element_blank())
}

pearson <- ggpubr::ggarrange(plotlist = plot.rho, ncol = 3, nrow = 4,
                             align = 'h', common.legend = TRUE, legend = 'bottom')
pearson <- ggpubr::annotate_figure(pearson, 
                                   left = ggpubr::text_grob("Age group", rot = 90,
                                                            vjust = 0.5, hjust = 0.2, 
                                                            size = 15))
pearson

# ggsave('../Figures updated/Pearson_all.pdf',
#        plot = pearson, height = 50, width = 36, units = 'cm')

##### 5) Compare with Poisson case ----

# Number of parameters
npar.POI <- ncol(C.base) + ncol(C.dlnm)
npar.NB  <- ncol(C.base) + ncol(C.dlnm) + length(obj$lcd$phix) +
  length(obj$lcd$phir) - 1

# Fit the two models NB - and POI - model
obj.pois <- fit_wLC.DLNM.POI(xv = xv, isoyv = isoyv, isowv = isowv,
                             etx = etxr, dtx = dtxr, adjmat = L, max.iter = 50,
                             cb = list('Temp' = cb1, 'Infl' = cb2),
                             lambda = lambda.opt)

# Poisson DLNM Hessian 
pois.as.nb      <- obj.pois$lcd
pois.as.nb$phix <- rep(0, length(xv))
pois.as.nb$phir <- rep(log(1e20), R)

H.dlnm.P    <- hessian.wLC.DLNM(etx = etxr, dtx = dtxr, fit = pois.as.nb)
H.dlnm.up.P <- hessian.wLC.DLNM.up(etx = etxr, dtx = dtxr, fit = pois.as.nb)
edf.dlnm.P  <- edf_constrained(H.dlnm.P, H.dlnm.up.P, contrast_dlnm(obj.pois$lcd))

# EDF models
edf.POI  <- ncol(C.base) + edf.dlnm.P
edf.NB   <- ncol(C.base) + edf.dlnm + 
  length(obj$lcd$phix) + length(obj$lcd$phir) - 1

# Likelihood ratio test
stat <- as.numeric(2 * (obj$lcd$ll.up - obj.pois$lcd$ll.up))
pchisq(stat, df = length(c(obj$lcd$phir, obj$lcd$phix)), lower.tail = FALSE)/2

# Overdispersion check
sum((dtxr - obj.pois$lcd$mhat * etxr)^2/(obj.pois$lcd$mhat * etxr))/
  (length(dtxr) - npar.POI)

# Regions
reg_vec <- c('32', '53', '84')
age_vec <- c('50', '70', '90') 

# AIC/BIC
c('NB' = -2*obj$lcd$ll.up + 2 * edf.NB, 
  'POI' = -2 * obj.pois$lcd$ll.up + 2 * edf.POI)
c('NB' = -2*obj$lcd$ll.up + log(length(dtxr)) * edf.NB, 
  'POI' = -2 * obj.pois$lcd$ll.up + log(length(dtxr)) * edf.POI)

c('POI' = prop.table(table(obj.pois$lcd$epsilon[,,reg_vec[1]]^2 > 3.841))['TRUE'],
  'NB'  = prop.table(table(obj$lcd$epsilon[,,reg_vec[1]]^2 > 3.841))['TRUE'])
c('POI' = prop.table(table(obj.pois$lcd$epsilon[,,reg_vec[2]]^2 > 3.841))['TRUE'],
  'NB'  = prop.table(table(obj$lcd$epsilon[,,reg_vec[2]]^2 > 3.841))['TRUE'])
c('POI' = prop.table(table(obj.pois$lcd$epsilon[,,reg_vec[3]]^2 > 3.841))['TRUE'],
  'NB'  = prop.table(table(obj$lcd$epsilon[,,reg_vec[3]]^2 > 3.841))['TRUE'])  

range(sapply(dimnames(dtxr)[[3]], function(r) 
  prop.table(table(obj$lcd$epsilon[,,r]^2 > 3.841))['TRUE']))
range(sapply(dimnames(dtxr)[[3]], function(r) 
  prop.table(table(obj.pois$lcd$epsilon[,,r]^2 > 3.841))['TRUE']))

sum(obj$lcd$epsilon^2)
sum(obj.pois$lcd$epsilon^2)
qchisq(0.95, df = length(dtxr) - edf.POI)
qchisq(0.95, df = length(dtxr) - edf.NB)

# Discrete scale colour gradient
cuts <- qchisq(c(0.75, 0.95, 0.995, 0.9995, 0.99995, 0.999995), df = 1)
vec  <- c(0, cuts, max(obj.pois$lcd$epsilon^2))
vec <- (vec[-length(vec)] + vec[-1]) / 2
cols <- c("#E6F4E4","#B2DF8A", "#E31A1C", "#C31A1C", "#8F140F","#4F140F", "black")

plot.rho <- list()
reg_vec  <- dimnames(dtxr)[[3]][c(5,8,11)]

for(r in 1:length(reg_vec)){
  eps_long <- as.data.frame(obj.pois$lcd$epsilon[,,reg_vec[r]])
  
  plot.rho[[r]] <- ggplot(eps_long) + 
    xlab(ifelse(r == 2, 'Date', '')) + 
    ylab(ifelse(r == 1, 'Age', '')) + 
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) + 
    theme_bw(base_size = 15) + 
    geom_tile(aes(x = as.Date(ISODate), y = age, fill = Freq^2)) +
    scale_fill_stepsn(
      name    = '',
      breaks  = c(0,cuts), 
      limits  = c(0, max(eps_long$Freq^2)),
      values  = scales::rescale(vec, from = c(0, max(eps_long$Freq^2))),
      colours = cols,
      labels  = function(x) sprintf("%.3f", x)) +
    theme(legend.position = 'bottom',
          legend.key.width = unit(3, 'cm'))
  
  if(r == 2 | r == 3| r == 1){
    plot.rho[[r]] <- plot.rho[[r]] + 
      theme(axis.title.y=element_blank())
  }
}

pearson <- ggpubr::ggarrange(plotlist = plot.rho, ncol = 3, nrow = 1,
                             align = 'h', common.legend = TRUE, legend = 'bottom')
pearson <- ggpubr::annotate_figure(pearson, 
                                   left = ggpubr::text_grob("Age group", rot = 90,
                                                            vjust = 0.5, hjust = 0.2, 
                                                            size = 15))
pearson

# ggsave('../Figures updated/Pearson_poisson.pdf',
#        plot = pearson, height = 15, width = 36, units = 'cm')

##### 6) Forecasting ----

# Regions
regions <- dimnames(dtxr)[[3]]

# Number of simulations
nY    <- 5
n.sim <- 5000

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

# Individual SARIMA - models for latent mortality index
sarima.kt = res.kt <- list()

for(r in 1:length(regions)){
  cat(regions[r], ' ... ')
  sarima.kt[[r]] <- forecast::Arima(ts(obj$lcd$kappa2[,r]),
                                    order = c(0,1,1),
                                    include.drift = TRUE)
  res.kt[[r]] <- residuals(sarima.kt[[r]], type = 'innovation')/
    sqrt(sarima.kt[[r]]$sigma2)
}
names(sarima.kt) <- regions

# Pseudo-observations of standardized residuals
res.mat.kt <- do.call(cbind, res.kt)[-c(1),]
u.mat.kt   <- copula::pobs(res.mat.kt)
cop.kt     <- copula::normalCopula(dim = ncol(u.mat.kt), dispstr = 'ex')
fit.cop.kt <- copula::fitCopula(cop.kt, data = u.mat.kt, method = "ml")

# Updated dataset
Dfnew <- fread(file = '../Results updated/Df_var_2025.txt') %>%
  dplyr::filter(ISOYear >= 2020,
                ISOYear <= 2024,
                gender == 2,
                age >= 50)

# Temperature
tg.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'tg')) %>%
  unique() %>%
  dplyr::arrange(region, ISODate)
tg.obs <- matrix(tg.obs$tg, ncol = R, nrow = length(unique(tg.obs$ISODate)),
                 dimnames = list(as.character(unique(tg.obs$ISODate)),
                                 unique(tg.obs$region)))

# Influenza
ia.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'ia100'))  %>%
  unique() %>%
  dplyr::arrange(region, ISODate)
ia.obs$ia100 <- ia.obs$ia100/1000
ia.obs <- ia.obs %>%
  dplyr::left_join(sum[,c('ISOWeek', 'region', 'mia100')], by = c('ISOWeek', 'region'))
ia.obs$ia100 <- ia.obs$ia100 - ia.obs$mia100
ia.obs <- matrix(ia.obs$ia100, ncol = R, nrow = length(unique(ia.obs$ISODate)),
                 dimnames = list(as.character(unique(ia.obs$ISODate)),
                                 unique(ia.obs$region)))

# Exposure
expo.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'expo')) %>%
  unique() %>%
  dplyr::arrange(age, region, ISODate)
expo.obs <- aperm(simplify2array(lapply(unique(expo.obs$age), function(a)
  matrix(expo.obs[expo.obs$age == a,]$expo,
         ncol = R, nrow = length(unique(expo.obs$ISODate)),
         dimnames = list(as.character(unique(expo.obs$ISODate)),
                         unique(expo.obs$region))))),
  c(1,3,2))
dimnames(expo.obs)[[2]] <- unique(Dfnew$age)

# Deaths
deaths.obs <- Dfnew %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'region', 'age', 'deaths')) %>%
  unique() %>%
  dplyr::arrange(age, region, ISODate)
deaths.obs <- aperm(simplify2array(lapply(unique(deaths.obs$age), function(a)
  matrix(deaths.obs[deaths.obs$age == a,]$deaths,
         ncol = R, nrow = length(unique(deaths.obs$ISODate)),
         dimnames = list(as.character(unique(deaths.obs$ISODate)),
                         unique(deaths.obs$region))))),
  c(1,3,2))
dimnames(deaths.obs)[[2]] <- unique(Dfnew$age)

# ISO year and weeks for Dfnew
isowv_fc <- isoweek(as.Date(sort(unique(Dfnew$ISODate))))
isoyv_fc <- isoyear(as.Date(sort(unique(Dfnew$ISODate)))) - min(DfM$ISOYear) + 1

# Cross-basis over 2020-2024: prepend the last weeks of the calibration period for
# the lags, drop them again, and centre with the calibration week-of-year means
means1 <- attr(cb1, 'means')
means2 <- attr(cb2, 'means')

cb_forecast <- function(new, hist, var, mxlag, means){
  cb <- simplify2array(lapply(seq_along(regions), function(j) {
    a <- args[[regions[j]]][[var]]
    x <- c(tail(hist[,j], mxlag), new[,j])
    cb_fun(dat = x, lags = 0:mxlag, argvar = a$argvar,
           arglag = a$arglag)[-(1:mxlag), , drop = FALSE]
  }))
  dimnames(cb)[[3]] <- regions
  center_cb(cb, isowv_fc, means = means)
}

obs.cb1 <- cb_forecast(tg.obs, hist.tg, 'Temp', mxlag.tavg, means1)
obs.cb2 <- cb_forecast(ia.obs, hist.ia, 'Infl', mxlag.infl, means2)

# DLNM part - does not depend on the simulation, so computed once
obs.dlnm <- Map(function(c, g){
  out <- sapply(seq_len(dim(c)[3]), function(r) c[,,r] %*% g[,r])
  out[is.na(out)] <- 0
  out
}, list('Temp' = obs.cb1, 'Infl' = obs.cb2), obj$lcd$gammal)

# Run
cl <- makeCluster(15)
registerDoParallel(cl)
clusterExport(cl, c('n.sim', 'regions', 'R', 'isowv_fc', 'isoyv_fc',
                    'fit.cop.kt', 'sarima.kt', 'obj', 'expo.obs', 'nY',
                    'obs.dlnm'))
clusterEvalQ(cl, expr = {
  library(dplyr)
  library(copula)
  library(forecast)
})
mhat.sim <- pblapply(cl = cl, X = 1:n.sim, FUN = function(s) {
  
  # Simulate kt-process
  sim.u.kt <- copula::rCopula(nY, fit.cop.kt@copula)
  sim.z.kt <- qnorm(sim.u.kt)
  sim.kt   <- sapply(1:length(regions), function(r) {
    simulate(object = sarima.kt[[r]], nsim = nY,
             innov = sim.z.kt[,r] * sqrt(sarima.kt[[r]]$sigma2))
  })
  
  # Force of mortality
  mu <- exp(sweep(aperm(obj$lcd$beta2 %o% sim.kt[isoyv_fc - min(isoyv_fc) + 1,] +
                          obj$lcd$beta3 %o% obj$lcd$kappa3[isowv_fc,] +
                          Reduce(`+`, Map(function(b, d)
                            b %o% d, obj$lcd$betal, obs.dlnm)),
                        c(2,1,3)), c(2,3), obj$lcd$beta1, '+'))
  
  # Dispersion
  exp.ds <- exp(outer(obj$lcd$phix, obj$lcd$phir, "+"))
  exp.ds <- replicate(nrow(mu), exp.ds, simplify = "array")
  exp.ds <- aperm(exp.ds, c(3, 1, 2))
  
  # Return
  array(rnbinom(length(mu), mu = mu * expo.obs, size = exp.ds), dim = dim(mu))
})
stopCluster(cl)

mhat.sim <- aperm(simplify2array(mhat.sim), c(4,1,2,3))

## Quantiles
quant.sim1 <- apply(mhat.sim, MARGIN = c(2, 3, 4),
                    FUN = function(x)
                      quantile(x, probs = c(0.025, 0.5, 0.975)))
dimnames(quant.sim1) <- list(rownames(quant.sim1), as.character(unique(Dfnew$ISODate)),
                             unique(Dfnew$age), regions)

df.sim1    <- as.data.frame.table(quant.sim1, responseName = "value")
df.deaths1 <- as.data.frame.table(deaths.obs)
df.expo1   <- as.data.frame.table(expo.obs)
colnames(df.sim1)    <- c('Quantile', 'ISODate', 'age', 'region', 'sim')
colnames(df.deaths1) <- c('ISODate', 'age', 'region', 'obs')
colnames(df.expo1)   <- c('ISODate', 'age', 'region', 'expo')

df.tot1 <- plyr::join_all(list(df.sim1, df.deaths1, df.expo1)) %>%
  tidyr::pivot_wider(names_from = Quantile,
                     values_from = sim,
                     names_prefix = "q")

DfM$age <- as.factor(DfM$age)

## Three regions
pvis3   <- list()
reg_vec <- c('32', '53', '84')
age_vec <- c('50', '70', '90')

for(r in 1:length(reg_vec)){
  dfs <- df.tot1 %>% dplyr::filter(region == reg_vec[r],
                                   age %in% age_vec)
  
  pvis3[[r]] <- ggplot(dfs) +
    theme_bw(base_size = 15) +
    xlab(ifelse(r == 2, 'Date', '')) +
    ylab(if(r == 1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(data = DfM %>% dplyr::filter(region == reg_vec[r], age %in% age_vec),
              aes(x = as.Date(ISODate), y = log(deaths/expo), group = age, col = age)) +
    geom_point(aes(x = as.Date(ISODate), y = log(obs/expo), group = age, col = age),
               shape = 3, col = 'black', cex = 0.2) +
    geom_ribbon(aes(x = as.Date(ISODate), ymin = log(`q2.5%`/expo), ymax = log(`q97.5%`/expo),
                    group = age, fill = age, col = age), alpha = 0.2) +
    geom_line(aes(x = as.Date(ISODate), y = log(`q50%`/expo), group = age, col = age)) +
    coord_cartesian(xlim = c(as.Date('2015-01-01'), as.Date('2024-12-31'))) +
    scale_fill_manual(name = 'Age',
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) +
    scale_color_manual(name = 'Age',
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94'))
}

plot.forecast1 <- ggpubr::ggarrange(plotlist = pvis3, common.legend = TRUE,
                                    legend = 'bottom', nrow = 1, ncol = 3,
                                    align = 'h', widths = c(1,1,1))

# ggsave('../Figures updated/forecast1.pdf',
#        plot = plot.forecast1, height = 12, width = 30, units = 'cm')

## All regions
pvis4   <- list()
reg_vec <- dimnames(dtxr)[[3]]
age_vec <- colnames(dtxr)[c(1,5,9)]

for(r in 1:length(reg_vec)){
  dfs <- df.tot1 %>% dplyr::filter(region == reg_vec[r],
                                   age %in% age_vec)
  
  pvis4[[r]] <- ggplot(dfs) +
    theme_bw(base_size = 15) +
    xlab(ifelse(r == 2, 'Date', '')) +
    ylab(if(r == 1) bquote('log('*mu['x,t,w,r']*')') else '') +
    ggtitle(shapef$NUTS_NAME[shapef$region == reg_vec[r]]) +
    geom_line(data = DfM %>% dplyr::filter(region == reg_vec[r], age %in% age_vec),
              aes(x = as.Date(ISODate), y = log(deaths/expo), group = age, col = age)) +
    geom_point(aes(x = as.Date(ISODate), y = log(obs/expo), group = age, col = age),
               shape = 3, col = 'black', cex = 0.2) +
    geom_ribbon(aes(x = as.Date(ISODate), ymin = log(`q2.5%`/expo), ymax = log(`q97.5%`/expo),
                    group = age, fill = age, col = age), alpha = 0.2) +
    geom_line(aes(x = as.Date(ISODate), y = log(`q50%`/expo), group = age, col = age)) +
    coord_cartesian(xlim = c(as.Date('2015-01-01'), as.Date('2024-12-31'))) +
    scale_fill_manual(name = 'Age',
                      values = c("#33A02C", "#1F78B4", "#E31A1C"),
                      labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94')) +
    scale_color_manual(name = 'Age',
                       values = c("#33A02C", "#1F78B4", "#E31A1C"),
                       labels = c('50' = '50-54', '70' = '70-74', '90' = '90-94'))
}

plot.forecast1_all <- ggpubr::ggarrange(plotlist = pvis4, common.legend = TRUE,
                                        legend = 'bottom', nrow = 4, ncol = 3,
                                        align = 'h')

# ggsave('../Figures updated/forecast1_all.pdf',
#        plot = plot.forecast1_all, height = 50, width = 36, units = 'cm')