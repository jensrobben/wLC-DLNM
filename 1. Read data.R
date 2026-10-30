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

##### 1) Read data ----

### Death records 

# Vector of years
yv   <- 1990:2025

# Year 1970 
file <- paste0('../Deaths/deces-', yv[1],'.csv')
dt   <- data.table::fread(input = file)
dt$lieudeces <- stringr::str_pad(dt$lieudeces, width = 5, pad = "0")

# Loop over other years
for(k in 2:length(yv)){
  message(yv[k])
  file <- paste0('../Deaths/deces-', yv[k],'.csv')
  df   <- data.table::fread(input = file)
  df$lieudeces <- stringr::str_pad(df$lieudeces, width = 5, pad = "0")
  dt   <- rbindlist(list(dt, df))
}

### Some pre-processing

# Check for duplicates - and remove them
ind = duplicated(dt[,c('nomprenom', 'sexe', 'datenaiss', 'datedeces', 'lieunaiss')])
dt <- dt[-c(which(ind)),]

# Set column names
colnames(dt) <- c("fullname", "gender", "birthdate", "birthplace", 
                  "birthmunicipality",  "birthcountry", "deathdate", 
                  "deathplace", "deathrecord")

dt <- dt %>% 
  dplyr::select(-c(fullname, birthplace, birthmunicipality, birthcountry,
                   deathrecord)) 

# Check for missing death dates and remove 
dt <- dt[!is.na(dt$deathdate)]

# Delete observations with wrong age
dt <- dt[! (dt$birthdate == 0 | dt$deathdate == 0)]

### Deal with birth and death dates

# Remove observations with birthdate > deathdate
dt <- dt[! dt$birthdate > dt$deathdate]

# Extract birth year, month, and day
dt[, birthyear := as.numeric(substr(birthdate, 1, 4))]
dt[, birthmonth := as.numeric(substr(birthdate, 5, 6))]
dt[, birthday := as.numeric(substr(birthdate, 7, 8))]

# Extract death year, month, and day
dt[, deathyear := as.numeric(substr(deathdate, 1, 4))]
dt[, deathmonth := as.numeric(substr(deathdate, 5, 6))]
dt[, deathday := as.numeric(substr(deathdate, 7, 8))]

# Check whenever birth/death year, month or day are missing
# Remove when missing year
dt <- dt[dt$birthyear != 0] 
dt <- dt[!is.na(dt$birthyear)]
dt <- dt[dt$deathyear != 0] 
dt <- dt[!is.na(dt$deathyear)]

# draw randomly when missing/wrong month or day
i1b <- which(dt$birthmonth == 0 | dt$birthmonth > 12 | is.na(dt$birthmonth))
i1d <- which(dt$deathmonth == 0 | dt$deathmonth > 12 | is.na(dt$deathmonth))
i2b <- which(dt$birthday == 0 | dt$birthday > 31 | is.na(dt$birthday))
i2d <- which(dt$deathday == 0 | dt$deathday > 31 | is.na(dt$deathday))

# Analysis of missing death records
sum(table(dt$deathyear[i1d])[as.character(1990:2019)], na.rm = TRUE)
sum(table(dt$deathyear[i2d])[as.character(1990:2019)], na.rm = TRUE)
table(dt$deathplace[i1d])
table(dt$deathplace[i2d])

set.seed(1)
probsi1b <- prop.table(table(dt$birthmonth[-c(i1b)]))[as.character(1:12)]
probsi1b <- probsi1b/sum(probsi1b)
table(dt[i1b,deathyear] - dt[i1b,birthyear])
dt$birthmonth[i1b] <- sample(1:12, length(i1b), prob = probsi1b, replace = TRUE)

set.seed(2)
probsi1d <- prop.table(table(dt$deathmonth[-c(i1d)]))[as.character(1:12)]
probsi1d <- probsi1d/sum(probsi1d)
dt$deathmonth[i1d] <- sample(1:12, length(i1d), prob = probsi1d, replace = TRUE)

set.seed(3)
years.i2b  <- dt$birthyear[i2b]
months.i2b <- dt$birthmonth[i2b]
years.i2d  <- dt$deathyear[i2d]
months.i2d <- dt$deathmonth[i2d]
leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

sim_day <- function(x){
  # Extract year and month
  y <- x[1]
  m <- x[2]
  
  # Check if leap year
  if(m == 2){
    if(leap_year(y)){
      sample(1:29, 1)
    } else{
      sample(1:28, 1)
    }
  } else if(m %in% c(1,3,5,7,8,10,12)){
    sample(1:31, 1)
  } else {
    sample(1:30, 1)
  }
}
dt$birthday[i2b] <- sapply(1:length(years.i2b), function(k)
  sim_day(x = c(years.i2b[k], months.i2b[k])))
dt$deathday[i2d] <- sapply(1:length(years.i2d), function(k)
  sim_day(x = c(years.i2d[k], months.i2d[k])))

# Compute exact age at death in years (as a decimal)
dt <- dt[, age := deathyear - birthyear + 
           (deathmonth - birthmonth) / 12 + 
           (deathday - birthday) / 365.25]

dt <- dt[dt$age >= 0 & dt$age < 125]

##### 2) Map communes to departments  ----

# Departments
dt$department <- substr(dt$deathplace, 1, 2)

# Focus on France métropolitaine
dt <- dt %>% 
  dplyr::filter(! department %in% c('- ', '00', '96', '97', '98', '99'))

##### 3) Merge by 5-year age groups ----

# 5-year age groups
dt$age <- pmin(floor(dt$age/5) * 5, 95)
sort(unique(dt$age))

##### 4) Aggregate by week, age group, sex, and department ----

# Dates
dt$deathdate <- ymd(sprintf("%04d-%02d-%02d", dt$deathyear, dt$deathmonth, dt$deathday))

# ISOYear, ISOWeek
dt$ISOYear <- isoyear(dt$deathdate)
dt$ISOWeek <- isoweek(dt$deathdate)

# Agregate by week, age, sex, and department
Df <- dt %>%
  dplyr::group_by(ISOYear, ISOWeek, age, gender, department) %>%
  dplyr::reframe(Deaths = n()) %>%
  ungroup() %>%
  mutate('ISODate' = ISOweek2date(paste0(sprintf("%04d-W%02d", ISOYear, ISOWeek), "-1")), 
         .before = 1)

# Start in 1990, end in 2024
Df <- Df %>%
  dplyr::filter(ISOYear >= 1990, ISOYear <= 2024)

# Some observations are missing -> zero death counts
tab <- expand.grid(ISODate = sort(unique(Df$ISODate)),
                   age = sort(unique(Df$age)),
                   gender = c(1,2),
                   department = sort(unique(Df$department)))

Df <- tab %>%
  dplyr::left_join(Df, by = c('ISODate', 'age', 'gender', 'department')) %>%
  dplyr::mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths),
                ISOYear = isoyear(ISODate),
                ISOWeek = isoweek(ISODate)) %>%
  dplyr::select(c('ISODate', 'ISOYear', 'ISOWeek', 'age', 
                  'gender', 'department', 'Deaths'))

# Plot
sub <- Df %>% dplyr::filter(age == 95, gender == 2, department == '75')
plot(sub$ISODate, sub$Deaths, type = 'l')

rm(df)
rm(dt)
rm(tab)
gc()

# # Save data
# fwrite(Df, '../Results updated/Df_dept.txt')

##### 5) From department to region codes ----
gc()

# Read saved Df
Df_dept <- fread('../Results updated/Df_dept.txt')

# Read in transition table
trans.tab <- fread('../v_departement_2023.csv') %>%
  dplyr::select(c('DEP', 'REG'))
colnames(trans.tab) <- c('department', 'region')

# Check if it contains all our departments
ind <- which(! unique(Df_dept$department) %in% trans.tab$department)
unique(Df_dept$department)[ind]

# There is still the old code for Corsica used -> add it 
trans.tab <- rbind(trans.tab, 
                   data.table(department = 20,
                              region = 94))
# Add to original data frame
Df_dept <- Df_dept %>%
  dplyr::left_join(trans.tab, by = c('department'))


##### 6) Add also region names ----
# Read in transition table
trans.tab <- fread('../v_region_2023.csv') %>%
  dplyr::select(c('REG', 'NCC'))
colnames(trans.tab) <- c('region', 'region_name')

# Check if it contains all our departments
ind <- which(! unique(Df_dept$region) %in% trans.tab$region)
unique(Df_dept$region)[ind]

# Add to original data frame
Df_dept <- Df_dept %>%
  dplyr::left_join(trans.tab, by = c('region'))

# Aggegrate by region
Df_reg <- Df_dept %>%
  dplyr::group_by(ISODate, ISOYear, ISOWeek, age, gender, region, region_name) %>%
  dplyr::reframe(deaths = sum(Deaths)) 

# # Save
# fwrite(Df_reg, '../Results updated/Df_reg.txt')

##### 7) Population data and exposures ----

# Read
Df_reg <- fread('../Results updated/Df_reg.txt')

# Read population data
file_path <- "../Population/estim-pop-dep-sexe-aq-1975-2026.xlsx"
sheets    <- excel_sheets(file_path)

# Load all sheets, inspecting structure first
read_year_sheet <- function(sheet) {
  message(sheet)
  df <- read_excel(file_path, sheet = sheet, skip = 0, .name_repair = "minimal",
                   col_types = c("text", "text", rep("numeric", 63)),
                   range = 'A5:BM101')
  
  colnames(df) <- c('department', 'department_name',
                    paste0('comb', seq(0,95,5)),'comb_total',
                    paste0('male', seq(0,95,5)),'male_total',
                    paste0('female', seq(0,95,5)),'female_total')
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with(c("comb", "male", "female")),
      names_to = c("gender", "age"),
      names_pattern = "([a-zA-Z]+)(\\d+)",
      values_to = "pop") %>%
    mutate(
      year = as.numeric(sheet),
      age = as.integer(age),
      gender = recode(gender,
                      "comb" = 0,
                      "male" = 1,
                      "female" = 2))
  
  df_long <- df_long %>% 
    dplyr::filter(gender != 0,
                  (! department %in% c('- ', '00', '96', '971', '972', 
                                       '973', '974', '976'))) %>%
    na.omit()
  
  return(df_long)
}

obj <- lapply(2:length(sheets), function(j){
  read_year_sheet(sheets[j])
})

df.pop <- rbindlist(obj)

# Define date
df.pop$date <- as.IDate(ymd(paste0(df.pop$year, '01', '01')))

# Read in again transition table
trans.tab <- fread('../v_departement_2023.csv') %>%
  dplyr::select(c('DEP', 'REG'))
trans.tab <- rbind(trans.tab, 
                   data.table(DEP = 20,
                              REG = 94))
colnames(trans.tab) <- c('department', 'region')

df.pop <- df.pop %>%
  dplyr::left_join(trans.tab, by = 'department')

# Aggregate by region
df.pop <- df.pop %>%
  dplyr::group_by(date, year, age, gender, region) %>%
  dplyr::reframe(pop = sum(pop)) 

# Interpolate on every day
dfunq <- df.pop[,c('gender','age','region')] %>% unique()

obj <- lapply(1:nrow(dfunq), function(j){
  s <- dfunq$gender[j]
  a <- dfunq$age[j]
  r <- dfunq$region[j]
  
  sub <- df.pop %>% 
    dplyr::filter(age == a, gender == s, region == r)
  
  dfs <- sub %>%
    group_by(gender, region, age) %>%
    complete(date = seq(min(as.IDate('1974-12-30')), 
                        max(as.IDate('2028-12-28')), by = "week")) %>%
    arrange(gender, region, age, date) %>%
    mutate(pop = approx(date, pop, date, method = "linear")$y) %>%
    dplyr::select(-c(year)) %>%
    na.omit()
  
  dfs
})

# New df.pop
df.pop <- rbindlist(obj)

# Calculate exposure
df.pop <- df.pop %>%
  group_by(age, gender, region) %>%
  arrange(date) %>%
  mutate('expo' = 1/52.18 * c((pop[-1] + pop[-length(pop)])/2, pop[length(pop)])) %>%
  ungroup() %>%
  dplyr::select(-c('pop')) %>%
  as.data.table()

# Change type for matching
df.pop$gender  <- as.integer(df.pop$gender)
df.pop$region  <- as.integer(df.pop$region)
df.pop$date    <- as.Date(df.pop$date)
Df_reg$ISODate <- as.Date(Df_reg$ISODate)

# Add to Df
Df_reg <- Df_reg %>%
  dplyr::left_join(df.pop, by = c('ISODate' = 'date', 'gender' = 'gender', 
                                  'region' = 'region', 'age' = 'age'))

# Check for missing observations -> after 2023 all missing obs
ind.missing <- which(is.na(Df_reg))
range(Df_reg[!complete.cases(Df_reg)]$ISODate)

# Remove missing observations 
Df_reg <- Df_reg %>% 
  na.omit()

# Plot
ggplot(Df_reg %>% dplyr::filter(gender == 2, region == 11, age >= 50)) + 
  theme_bw(base_size = 15) + ylab(bquote(log('Death rate'))) + 
  geom_line(aes(x = ISODate, y = log(deaths/expo), group = as.factor(age), color = as.factor(age)))

# # Save
# fwrite(Df_reg, file = '../Results updated/Df_final.txt')

##### 8) Influenza data ----

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

# ggsave('../Figures updated/ia100.pdf',
#        plot = pia, height = 12.5, width = 25, units = 'cm')

##### 9) Temperature data ----

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
  dplyr::reframe(hi = max((tg + lag(tg, 1) + lag(tg, 2)) / 3,na.rm = TRUE),
                 tg = mean(tg)) %>%
  ungroup()

# Adaptive threshold
df.temp <- df.temp %>%
  arrange(Region, ISODate) %>%
  group_by(Region) %>%
  mutate(
    # Compute adaptive upper (hot) and lower (cold) thresholds
    past_q90 = map_dbl(seq_along(tg), function(i) {
      if (i == 1) return(NA_real_)
      quantile(tg[1:(i-1)], probs = 0.95, na.rm = TRUE)
    }),
    past_q10 = map_dbl(seq_along(tg), function(i) {
      if (i == 1) return(NA_real_)
      quantile(tg[1:(i-1)], probs = 0.05, na.rm = TRUE)
    }),
    
    # Compute adaptive exceedance variables
    exceed_heat = pmax(0, tg - past_q90),
    exceed_cold = pmax(0, past_q10 - tg),
    
    # Binary flags
    is_heat_extreme = tg > past_q90,
    is_cold_extreme = tg < past_q10
  ) %>%
  filter(ISOYear >= 1985) %>%
  ungroup()


##### 10) Construct final data set ----

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

# # Save
# fwrite(Df, file = '../Results updated/Df_var.txt')

### Some plots
Df <- fread(file = '../Results updated/Df_var.txt')

# Overview death counts and exposures
dff1 <- Df %>%
  dplyr::filter(ISOYear %in% 1990:2019,
                gender == 2, 
                age >= 50) %>%
  group_by(region, ISODate) %>%
  summarize(deaths = sum(deaths),
            expo = sum(expo))

regg <- dff1 %>%
  group_by(region) %>%
  summarize('TD' = sum(deaths)) %>%
  arrange(TD) %>%
  pull(region)

dff1$region <- factor(dff1$region, levels = regg)
cols <- colors[1:12]
names(cols) <- sort(unique(dff1$region))

p1 <- ggplot(dff1) + 
  theme_bw(base_size = 15) + xlab('Date') + 
  geom_bar(aes(x = ISODate, y = deaths, group = region, fill = region),
           stat = 'identity', width = 7) + 
  scale_fill_manual(name = 'Region',
                    values = colors[1:12],
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 

p2 <- ggplot(dff1) + 
  theme_bw(base_size = 15) + ylab('Population exposure') + xlab('Date') + 
  geom_bar(aes(x = ISODate, y = expo, group = region, fill = region),
           stat = 'identity', width = 7) + 
  scale_fill_manual(name = 'NUTS 2',
                    values = colors[1:12],
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 



p1p2 <- ggpubr::ggarrange(plotlist = list(p1, p2),
                          common.legend = TRUE, 
                          legend = 'bottom',
                          align = 'hv')

# ggsave('../Figures updated/overview.pdf',
#        plot = p1p2, height = 12.5, width = 25, units = 'cm')

# Plot death rates
df.reg11 <- Df %>% 
  dplyr::filter(age %in% c('50','60','70','80','90'),
                gender == 2,
                region %in% c(11),
                ISOYear %in% 1990:2019)
df.reg12 <- df.reg11 %>%
  group_by(ISOWeek, age) %>%
  dplyr::reframe('m' = mean(log(deaths/expo))) %>%
  ungroup() %>% 
  group_by(age) %>%
  dplyr::reframe(ISOWeek, m,
                 'mc' = m - mean(m))

p1 <- ggplot(df.reg11) + 
  theme_bw(base_size = 15) + xlab('Date') + ylab(bquote(log(m['x,t,w,r']))) + 
  geom_line(aes(x = ISODate, y = log(deaths/expo), group = as.factor(age), 
                col = as.factor(age)), linewidth = 0.2) + 
  scale_color_manual(name = 'Age group',
                     values = unname(colors[1:5]),
                     labels = c('50' = '50-54',
                                '60' = '60-64',
                                '70' = '70-74',
                                '80' = '80-84',
                                '90' = '90-94')) +
  guides(col = guide_legend(override.aes = list(linewidth = 1.5, size = 5)))

p2 <- ggplot(df.reg12) + 
  theme_bw(base_size = 15) + xlab('Week') + ylab(bquote('avg'*(log(m['x,t,w,r'])))) + 
  geom_line(aes(x = ISOWeek, y = mc, group = as.factor(age), 
                col = as.factor(age)), linewidth = 0.8) + 
  scale_color_manual(name = 'Age group',
                     values = unname(colors[1:5]),
                     labels = c('50' = '50-54',
                                '60' = '60-64',
                                '70' = '70-74',
                                '80' = '80-84',
                                '90' = '90-94'))

p12 <- ggpubr::ggarrange(plotlist = list(p1,p2), common.legend = TRUE, 
                         legend = 'bottom', align = 'hv')

# ggsave('../Figures updated/deathrate.pdf',
#        plot = p12, height = 12.5, width = 25, units = 'cm')
