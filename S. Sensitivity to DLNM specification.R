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
source('F. Weekly Lee-Carter (POI).R')
source('F. Weekly DLNM (POI).R')
source('F. Weekly Lee-Carter + DLNM (POI).R')
source('F. Bootstrap CI.R')


##### 1) Set up everything ----
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
  reframe('mia100' = quantile(ia100,0.90)) %>%
  as.data.table()
DfM <- DfM %>%
  dplyr::left_join(sum, by = c('ISOWeek', 'region'))

DfM$ia100 <- DfM$ia100 - DfM$mia100

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

bnd.infl <- infl %>%
  dplyr::summarise(
    min = floor(min(ia100)),
    max = ceiling(max(ia100))
  ) %>%
  unlist(use.names = FALSE)

##### 2) Sensitivity ----

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

# Function for sensitivity specifications
fit_sens <- function(
    mxlag.tavg  = 4,
    mxlag.infl  = 8,
    argvar.tavg = list(fun = 'bs', degree = 3,
                       knots = quantile(tavg[, tg], c(0.1, 0.9))),
    argvar.infl = list(fun = 'poly', degree = 1),
    arglag.tavg = list(fun = 'ns', intercept = TRUE, nk = 2),
    arglag.infl = list(fun = 'ns', intercept = TRUE, nk = 2),
    label       = NA_character_,
    dat.tavg    = tavg, 
    dat.infl    = infl,
    bnd.t       = bnd.tavg, 
    bnd.i       = bnd.infl,
    dtx         = dtxr, 
    etx         = etxr,
    agev        = xv, 
    yearv       = isoyv, 
    weekv       = isowv,
    adjmat      = L,
    max.iter    = 50, 
    verbose     = TRUE) {
  
  regions <- dimnames(dtx)[[3]]
  
  ## --- lag bases: knots follow the maximum lag, boundaries stay at (0, L)
  mk_arglag <- function(a, mxlag) {
    if (is.null(a$knots) && !is.null(a$nk))
      a$knots <- dlnm::logknots(mxlag, nk = a$nk)
    a$nk <- NULL
    if (a$fun %in% c('ns', 'bs') && is.null(a$Boundary.knots))
      a$Boundary.knots <- c(0, mxlag)
    if (a$fun == 'lin') {
      a$knots <- NULL
      a$Boundary.knots <- NULL
    }
    a
  }
  arglag.tavg <- mk_arglag(arglag.tavg, mxlag.tavg)
  arglag.infl <- mk_arglag(arglag.infl, mxlag.infl)
  
  ## --- arguments per region (from 11.2) -------------------------------
  args <- lapply(regions, function(r) {
    av.t <- argvar.tavg
    if (av.t$fun %in% c('bs', 'ns') && is.null(av.t$Boundary.knots))
      av.t$Boundary.knots <- bnd.t
    
    av.i <- argvar.infl
    if (av.i$fun %in% c('bs', 'ns') && is.null(av.i$Boundary.knots))
      av.i$Boundary.knots <- bnd.i
    if (av.i$fun == 'poly' && is.null(av.i$scale))
      av.i$scale <- max(dat.infl[region == r, ia100])
    
    list('Temp' = list(argvar = av.t, arglag = arglag.tavg),
         'Infl' = list(argvar = av.i, arglag = arglag.infl))
  })
  names(args) <- regions
  
  ## --- cross-basis matrices (from 11.2) -------------------------------
  cb1 <- simplify2array(lapply(regions, function(r)
    cb_fun(dat    = dat.tavg[region == r, tg],
           lags   = 0:mxlag.tavg,
           argvar = args[[r]]$Temp$argvar,
           arglag = args[[r]]$Temp$arglag)))
  cb2 <- simplify2array(lapply(regions, function(r)
    cb_fun(dat    = dat.infl[region == r, ia100],
           lags   = 0:mxlag.infl,
           argvar = args[[r]]$Infl$argvar,
           arglag = args[[r]]$Infl$arglag)))
  
  cb1 <- center_cb(cb1, weekv)
  cb2 <- center_cb(cb2, weekv)
  
  ## --- fit, zero penalty ----------------------------------------------
  obj <- fit_wLC.DLNM.NB(xv = agev, isoyv = yearv, isowv = weekv,
                         etx = etx, dtx = dtx, adjmat = adjmat,
                         cb = list('Temp' = cb1, 'Infl' = cb2),
                         lambda   = c('Temp' = 0, 'Infl' = 0),
                         max.iter = max.iter, verbose = verbose)
  
  ## --- BIC (from 11.4), edf = npar because lambda = 0 -----------------
  Rn   <- dim(dtx)[3]
  npar <- length(obj$lcd$beta1) + length(obj$lcd$beta2) +
    length(obj$lcd$beta3) + length(unlist(obj$lcd$betal)) +
    length(obj$lcd$kappa2) + length(obj$lcd$kappa3) +
    length(unlist(obj$lcd$gammal)) + length(obj$lcd$phix) +
    length(obj$lcd$phir) - (4 + 2 * Rn + 1)
  nobs <- length(obj$lcd$mhat)
  ll   <- obj$lcd$ll.up
  
  dsc <- function(a) paste0(a$fun,
                            if (!is.null(a$degree)) paste0('-', a$degree) else '',
                            if (!is.null(a$knots)) paste0(' (', length(a$knots), ' knots)') else '')
  
  list(obj = obj,
       cb  = list('Temp' = cb1, 'Infl' = cb2),
       summary = data.frame(
         label      = label,
         mxlag.tavg = mxlag.tavg,
         mxlag.infl = mxlag.infl,
         var.tavg   = dsc(argvar.tavg),
         var.infl   = dsc(argvar.infl),
         lag.tavg   = dsc(arglag.tavg),
         lag.infl   = dsc(arglag.infl),
         Q1         = dim(cb1)[2],
         Q2         = dim(cb2)[2],
         npar       = npar,
         nobs       = nobs,
         ll         = ll,
         BIC        = -2 * ll + log(nobs) * npar,
         cvg        = obj$cvg,
         itr        = obj$itr,
         stringsAsFactors = FALSE))
}

### the specifications to compare 
specs <- list(
  list(label = 'Proposed'),
  
  ## maximum lag, temperature
  list(label = 'Lag temp 2', mxlag.tavg = 2),
  list(label = 'Lag temp 3', mxlag.tavg = 3),
  list(label = 'Lag temp 5', mxlag.tavg = 5),
  list(label = 'Lag temp 6', mxlag.tavg = 6),
  
  ## maximum lag, influenza
  list(label = 'Lag ILI 6',  mxlag.infl = 6),
  list(label = 'Lag ILI 7',  mxlag.infl = 7),
  list(label = 'Lag ILI 9',  mxlag.infl = 9),
  list(label = 'Lag ILI 10', mxlag.infl = 10),
  
  ## covariate basis, temperature
  list(label = 'Temperature: natural spline (same knots)',
       argvar.tavg = list(fun = 'ns',
                          knots = quantile(tavg[, tg], c(0.1, 0.9)))),
  list(label = 'Temperature: B-spline (quadratic)',
       argvar.tavg = list(fun = 'bs', degree = 2,
                          knots = quantile(tavg[, tg], c(0.1, 0.9)))),
  list(label = 'Temperature: B-spline (quartic)',
       argvar.tavg = list(fun = 'bs', degree = 4,
                          knots = quantile(tavg[, tg], c(0.1, 0.9)))),
  list(label = 'Temperature: B-spline (knots at 10th, 50th, 90th perc.)',
       argvar.tavg = list(fun = 'bs', degree = 3,
                          knots = quantile(tavg[, tg], c(0.1, 0.5, 0.9)))),
  list(label = 'Temperature: B-spline (knots at 10th, 75th, 90th perc.)',
       argvar.tavg = list(fun = 'bs', degree = 3,
                          knots = quantile(tavg[, tg], c(0.1, 0.75, 0.9)))),
  list(label = 'Temperature: P-spline (cubic, df = 5)',
       argvar.tavg = list(fun = 'ps', degree = 3, df = 5)),  
  
  ## ILI exposure-response
  list(label = 'ILI: quadratic (second-degree polynomial)',
       argvar.infl = list(fun = 'poly', degree = 2)),
  list(label = 'ILI: B-spline (linear, knots at 10th, 90th perc.)',
       argvar.infl = list(fun = 'bs', degree = 1,
                          knots = quantile(infl[, ia100], c(0.1, 0.9)))),
  list(label = 'ILI: B-spline (quadratic, knots at 10th, 90th perc.)',
       argvar.infl = list(fun = 'bs', degree = 2,
                          knots = quantile(infl[, ia100], c(0.1, 0.9)))),
  list(label = 'ILI: B-spline (cubic, knots at 10th, 90th perc.)',
       argvar.infl = list(fun = 'bs', degree = 3,
                          knots = quantile(infl[, ia100], c(0.1, 0.9)))),
  
  ## Lag-response
  list(label = 'Lag-response (temperature, ILI): B-spline',
       arglag.tavg = list(fun = 'bs', degree = 3, intercept = TRUE, nk = 2),
       arglag.infl = list(fun = 'bs', degree = 3, intercept = TRUE, nk = 2)),
  list(label = 'Lag-response (ILI): linear spline',
       arglag.infl = list(fun = 'lin', intercept = TRUE, knots = numeric(0)))
)

### run 
res <- vector('list', length(specs))
for (k in seq_along(specs)) {
  cat('\n----', specs[[k]]$label, '----\n')
  res[[k]] <- do.call(fit_sens, specs[[k]])
  res[[k]]$obj <- NULL   
  gc()
}

# saveRDS(res, '../Results updated/Sensitivity/results.rds')
res <- readRDS('../Results updated/Sensitivity/results.rds')

tab <- do.call(rbind, lapply(res, `[[`, 'summary'))
tab$dBIC <- round(tab$BIC - tab$BIC[1], 1)
tab

fmt <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = "\\\\,")
rows <- paste0(tab$label, ' & ', fmt(tab$npar), ' & ', fmt(round(tab$ll)), ' & ',
               fmt(round(tab$BIC)), ' & ', fmt(round(tab$dBIC)), ' \\\\')
writeLines(rows)

i = 1
cat(fmt(round(tab$ll[i])), "\n")

cat(tab$npar[i], ' & ', paste0(fmt(tab$ll[i]), ' & ', fmt(tab$BIC[i])), "\n")
