##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'stats',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr','httr',
              'utils', 'data.table', 'viridis')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source function
source('F. Diverse.R')
source('F. Population Weights.R')

# Country/countries of interest
ctry.spec <- c("FR")

##### 1) Load shape file NUTS 3 European regions ----

# Read shape file - downloaded from Eurostat
shapefile <- read_sf('../NUTS/NUTS_RG_03M_2024_3035.shp')

# Extract shapefile for NUTS 2 regions in the countries of interest 
nuts.spec <- 1
shapef    <- shapefile[shapefile$LEVL_CODE == 1 & 
                         shapefile$CNTR_CODE %in% c('FR') &
                         ! shapefile$NUTS_ID %in% c('FRY'),] %>%
  dplyr::arrange(NUTS_NAME)

# Box long lat coordinates where country falls into
coords_ctry <- st_coordinates(st_transform(shapef, 4326)) %>% 
  as.data.frame() %>% dplyr::select('X','Y')
long.range <- range(coords_ctry[,1])
lat.range  <- range(coords_ctry[,2])

# plot
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
            "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
            "#FFFF99", "#B15928", "#66C2A5", "#FC8D62", "#8DA0CB", 
            "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
            "#FFB3BA", "#D4E157", "#D9B2E6", "#F1E8AE", "#A1B2A0", 
            "#FFC3A0", "#C2D6DA")

shapeff <- shapef %>% 
  dplyr::filter(NUTS_ID != 'FRM')

shapeff$regcode <- plyr::mapvalues(shapeff$NUTS_ID,
                                   from = unique(shapeff$NUTS_ID),
                                   to = c(84, 27, 53, 24, 44, 32,
                                          11, 28, 75, 76, 52, 93))
shapeff <- shapeff %>% 
  dplyr::arrange(regcode)

cols <- colors[c(1,12,11,10,5,6,9,8,3,4,2,7)]
names(cols) <- shapeff$regcode
shapeff$regcode <- as.factor(shapeff$regcode)

p <- ggplot(shapeff) + 
  theme_bw(base_size = 15) + 
  geom_sf(aes(fill = regcode), show.legend = TRUE) + 
  scale_fill_manual(values = cols, 
                    name = '',
                    labels = paste0(shapeff$regcode, ': ', shapeff$NUTS_NAME)) 

# ggsave('../Figures updated/regions.pdf',
#        plot = p, height = 10, width = 20, units = 'cm')


##### 2) Read E-OBS data ----
# Select weather variable
var <- 'tg'

file   <-  paste0('../Temperature/tg_ens_mean_0.1deg_reg_v31.0e.nc')
nc_ds  <- nc_open(file)

# Extract longitude, latitude coordinates and time object
long <- ncvar_get(nc_ds, "longitude") 
lat  <- ncvar_get(nc_ds, "latitude")
time <- ncvar_get(nc_ds, "time")

# Put as date object
t_units <- ncatt_get(nc_ds, "time", "units")
t_ustr  <- strsplit(t_units$value, " ")
t_dstr  <- strsplit(unlist(t_ustr)[3], "-")
date    <- ymd(t_dstr) + ddays(time)

# Indices to filter the spatial dataframes to the long/lat of interest
id.long <- which(nc_ds$dim$longitude$vals >= long.range[1] & 
                   nc_ds$dim$longitude$vals <= long.range[2])
id.lat  <- which(nc_ds$dim$latitude$vals >= lat.range[1] &
                   nc_ds$dim$latitude$vals <= lat.range[2])
id.time <- which(nc_ds$dim$time$vals >=  as.Date('1979-12-31') - 
                   as.Date('1950-01-01') & nc_ds$dim$time$vals <= 
                   as.Date('2024-12-31') - as.Date('1950-01-01'))
long    <- long[id.long]
lat     <- lat[id.lat]
date    <- date[id.time]

# Store the data in a 3-dimensional array 
var.array <- ncvar_get(nc_ds, 
                       start = c(id.long[1], id.lat[1], id.time[1]),
                       count = c(length(id.long), length(id.lat), length(id.time)),
                       verbose = TRUE) 
#var.array <- var.array[id.long, id.lat, id.time]
dimnames(var.array) <- list(long, lat, as.character(date))

# Plot the climate feature on world map at random time point date
pick.date <- "2003-07-20"
subarray  <- var.array[,,as.character(pick.date)] %>% as.data.frame() %>%
  rownames_to_column('long') %>%
  tidyr::gather(key = 'lat', value = val, - long) %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(long)) %>%
  na.omit()
subarray <- st_as_sf(subarray, coords = c(1,2), crs = 4326)

p.val <- ggplot() +
  geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
  geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
  geom_sf(data = subarray$geometry, aes(col = as.numeric(subarray$val))) +
  theme_bw(base_size = 15) + xlab('') + 
  theme(legend.position = 'bottom') + 
  scale_colour_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                         name = '') +
  geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
  geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
  ggtitle(paste0('Tavg',': ', pick.date))
p.val

# ggsave('../Figures updated/tempdata.pdf',
#        plot = p.val, height = 12.5, width = 12.5, units = 'cm')


# Plot population grid data
ncpop <- nc_open(filename = "../Population/PopulationCount.nc")
longp <- ncvar_get(ncpop, "longitude") 
latp  <- ncvar_get(ncpop, "latitude")
grpop <- ncvar_get(ncpop) 
dimnames(grpop) <- list(longp, latp, 1:20)
grpop <- grpop[longp >= min(long) & longp <= max(long),
               latp >= min(lat) & latp <= max(lat),5]

arrayp  <- grpop %>% as.data.frame() %>%
  rownames_to_column('long') %>%
  tidyr::gather(key = 'lat', value = val, - long) %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(long)) %>%
  na.omit()
subarrayp <- st_as_sf(arrayp, coords = c(1,2), crs = 4326)

p.valp <- ggplot() +
  geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
  geom_sf(data = subarrayp$geometry, aes(col = log(as.numeric(subarrayp$val)))) +
  theme_bw(base_size = 15) + xlab('Long') + 
  theme(legend.position = 'bottom') + 
  scale_color_gradient2(low = 'gray99', mid = 'gray80', high = 'black', midpoint = 3,
                        name = '') +
  geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
  geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
  ggtitle(paste0('Population count (log)',': ', 2020))
p.valp

# ggsave(plot = p.valp,
#        filename = paste0('../Figures updated/sedacfr.pdf'),
#        width = 4.5, height = 5)

##### 3) Intersection shape file with E-OBS grid ----

# All grid points in E-OBS data 
dat <- expand.grid(long, lat)
colnames(dat) <- c('Longitude', 'Latitude')

# Put in same CRS
eobs_sf  <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
nuts3_sf <- st_transform(shapef, 4326)     

# Make intersection
ints <- st_intersects(nuts3_sf, eobs_sf, sparse = TRUE)
for(k in 1:length(ints)){
  dat[ints[[k]],'region'] <- nuts3_sf$NUTS_ID[k]
}

# Keep E-OBS grid points in NUTS 3 regions
coords.ctry <- dat %>% na.omit() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 

# Visualize the intersection
n <- length(shapef$NUTS_ID)
palette <- createPalette(n,  c("#ff0000", "#00ff00", "#0000ff"))

ggplot(data = shapef) +
  geom_sf() + theme_bw() + 
  geom_sf(data = coords.ctry, aes(col = region), alpha = 0.5) + 
  scale_color_manual(values = unname(palette), name = '')  + 
  theme(legend.position = "none")

##### 4) Average climate feature over grid points in each NUTS region ----

# Remove missing data
df <- dat %>% na.omit()

# Calculate population weights
wpop <- get_pop_weights(ctry.spec = ctry.spec, long.ext = long, lat.ext = lat)
df   <- df %>% dplyr::left_join(wpop, by = c('Longitude' = 'Long.ext', 
                                             'Latitude'  = 'Lat.ext', 
                                             'region'    = 'Region'))

# If NA - zero population weight
df[which(is.na(df), arr.ind = TRUE)] <- 0

# Small check for population weights to sum up to 1 for each region
test <- df %>% dplyr::group_by(region) %>%
  summarize(sum(WPC2000, na.rm = TRUE), sum(WPC2005, na.rm = TRUE), 
            sum(WPC2010, na.rm = TRUE), sum(WPC2015, na.rm = TRUE), 
            sum(WPC2020, na.rm = TRUE))

# Split coordinates by region
coord.region <- df %>% split(f = df$region)

# Function to calculate population weighted average over each region
avg.clim.region <- function(df.sub){
  message(df.sub$region[1])
  lo.ind <- match(df.sub$Longitude, long)
  la.ind <- match(df.sub$Latitude, lat)
  obj <- sapply(1:nrow(df.sub), function(s) 
    var.array[lo.ind[s],la.ind[s],])
  
  weight.matrix <- matrix(rep(df.sub$WPC2015, times = nrow(obj)), 
                          nrow = nrow(obj), ncol = ncol(obj), byrow = TRUE)
  mm <- rowSums((!is.na(obj)) * weight.matrix)
  
  objw <- obj * weight.matrix
  
  v.objw <- rowSums(objw, na.rm = TRUE)/mm
  names(objw) <- dimnames(var.array)[[3]]
  v.objw
}

df.region.CLI <- sapply(coord.region, function(sub) 
  avg.clim.region(sub))

##### 5) Function to transform matrix to long data frame r (see ribiosUtils) ----
matrix2longdf <- function(mat,
                          row.names, col.names,
                          longdf.colnames=c("row","column","value")) {
  if(missing(row.names)) row.names <- rownames(mat)
  if(missing(col.names)) col.names <- colnames(mat)
  
  if(is.null(row.names)) row.names <- 1:nrow(mat)
  if(is.null(col.names)) col.names <- 1:ncol(mat)
  
  value <- as.vector(mat)
  if(length(row.names)!=nrow(mat))
    warning("row.names is inconsistent with the matrix dim")
  if(length(col.names)!=ncol(mat))
    warning("col.names is inconsistent with the matrix dim")
  
  rn <- rep(row.names, ncol(mat))
  cn <- rep(col.names, each=nrow(mat))
  res <- data.frame(row=rn,
                    column=cn,
                    value=value)
  colnames(res) <- longdf.colnames
  return(res)
}

df.region.CLI.daily  <- matrix2longdf(df.region.CLI,
                                      longdf.colnames = c('Date', 'Region', var))  

df.region.CLI.daily$Date <- as.Date(df.region.CLI.daily$Date)

# file.save <- paste0(paste0(ctry.spec, collapse = '_'),'_NUTS1_',var,'_daily','.txt')
# fwrite(df.region.CLI.daily, file = paste0('../Temperature/', file.save))

##### 6) Plot ----
file.save <- paste0(paste0("FR", collapse = '_'),'_NUTS1_',"tg",'_daily','.txt')
Df <- fread(file = paste0('../Temperature/', file.save)) %>%
  dplyr::left_join(shapeff, by = c('Region' = 'NUTS_ID'))
borders <- shapeff %>% group_by() %>% summarize('geom1' = st_union(geometry))

ppp <- ggplot(Df %>%
         dplyr::filter(Date == as.Date('2003-07-14'))) + 
  geom_sf(aes(fill = tg, geometry = geometry)) + 
  geom_sf(data = shapeff$geometry, aes(), col = 'gray70', fill = NA) + 
  scale_fill_viridis(name = '', option = 'rocket', direction = -1,
                     na.value = 'white') + 
  scale_color_viridis(name = '', option = 'rocket', direction = -1,
                      na.value = 'white') + 
  ggtitle(paste0('Tavg: ', 2003, '-',29)) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines")))
ppp <- ggpubr::ggarrange(plotlist = list(p.val, ppp), align = 'hv')

# ggsave('../Figures updated/tempaggr.pdf',
#        plot = ppp, height = 12.5, width = 25, units = 'cm')
