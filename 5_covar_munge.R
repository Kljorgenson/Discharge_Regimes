#### Covariate organization
# Compile data from csv file output from clipping ERA5 data to watersheds
library(tidyverse)
library(tidyhydat)
library(ggplot2)

filelist = list.files(path = "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/Region", pattern = "*snowfall*", full.names = TRUE, recursive = T)
filelist
#filelist <- filelist[c(1,4,6,9)]

df_list <- lapply(filelist, read.csv) 
names(df_list) <- basename( filelist )

df_list_named <- map2(df_list, names(df_list), ~mutate(.x, source = .y))

# Make all columns characters
df_list_named <- lapply(df_list_named, function(x) {x[] <- lapply(x, as.character); x})

cov.raw <- df_list_named %>% reduce(full_join)
head(cov.raw)

# Extract variable from file name
cov.raw$variable <- substr(cov.raw$source, start = 6, stop = 200)
cov.raw$variable <- gsub(".([0-9]+)|.nc.csv","",cov.raw$variable)

head(cov.raw)


# Remove duplicated time steps
covs.hrly <- cov.raw %>% dplyr::select(date_UTC, time, station, variable, value) 

# Check for duplicates in all but value
any(duplicated(covs.hrly[,1:4]))

names(covs.hrly) <- c("date_UTC", "hour_UTC", "station", "variable", "value")

head(covs.hrly)
tail(covs.hrly)
covs.hr <- covs.hrly

### Take monthly averages
# Based on timezone of discharge data

# Match stations and UTC-X
meta.ak <- read.csv("Data/river_metadata.csv") %>% dplyr::select(station) %>% unique()
meta.ak$UTC_minus <- 9
meta.ak$station <- as.character(meta.ak$station)
meta.can <- hy_stations() %>% filter(STATION_NUMBER %in% covs.hr$station) %>%
  rename(station = STATION_NUMBER) %>%
  mutate(UTC_minus = case_when(PROV_TERR_STATE_LOC %in% c("YT", "BC") ~ 8,
                               PROV_TERR_STATE_LOC %in% c("NT", "AB", "NU") ~ 7,
                               PROV_TERR_STATE_LOC %in% c("MB", "SK") ~ 6,
                               PROV_TERR_STATE_LOC %in% c("ON", "QC") ~ 5,
                               PROV_TERR_STATE_LOC %in% c("NL") ~ 4,
                               )) %>% dplyr::select(station, PROV_TERR_STATE_LOC, UTC_minus)
stations <- hy_stations()
write.csv(stations, "Data/can.stations.csv", row.names = F)
meta <- full_join(meta.ak, meta.can)
head(meta)         

# Adjust datetime
covs.hr$station <- as.character(covs.hr$station)
covs.hr <- left_join(covs.hr, meta)
covs.hr$time <- ifelse(nchar(covs.hr$hour_UTC) == 1, paste(0,covs.hr$hour_UTC, sep = ""), covs.hr$hour_UTC)

covs.hr$datetime_UTC <- as.POSIXct(paste(covs.hr$date_UTC, " ", covs.hr$time, ":00:00", sep = ""), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
covs.hr$datetime_adj <- covs.hr$datetime_UTC-covs.hr$UTC_minus* 3600
covs.hr$date_adj <- as.Date(covs.hr$datetime_adj)
covs.hr$year <- year(covs.hr$datetime_adj)
covs.hr$month <- month(covs.hr$datetime_adj)

tail(covs.hr)


### Calculate monthly average
covs.hr$value <- as.numeric(covs.hr$value)
covs.mo <- covs.hr %>% group_by(station, year, month, variable) %>%
  summarise(n = length(variable), value = if_else(first(variable) %in% c("temm", "temm-add.csv", "temm.nc-ad.csv"), mean(value), sum(value))) %>% 
  filter(n > 600) %>% dplyr::select(-c(n))

head(covs.mo)


# Export monthly data

write.csv(covs.mo, "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/Region/Precipitation_monthly.csv", row.names = F)

# Daily
covs.day <- covs.hr %>% mutate(jday = yday(datetime_adj)) %>% group_by(station, year, month, jday, variable) %>%
  summarise(n = length(variable), value = if_else(first(variable) %in% c("temm", "temm-add.csv", "temm.nc-ad.csv"), mean(value), sum(value))) %>% 
  filter(n == 24) %>% dplyr::select(-c(n))

tail(covs.day)
write.csv(covs.day, "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/Region/Precipitation_daily.csv", row.names = F)

################################################################

#### Combine all monthly data
filelist = list.files(path = "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/", pattern = "Covariates_m*", full.names = TRUE)
filelist


covars0 <- lapply(filelist, read.csv)
covars <- lapply(covars0, function(x) {x[] <- lapply(x, as.character); x}) %>% reduce(full_join) %>% unique() #%>% group_by(station, year, month, variable) %>% filter(values == max(values, na.rm = T))

covars$values <- ifelse(is.na(covars$values), covars$value, covars$values)
covars$values <- as.numeric(covars$values)
covrs <- covars %>% filter(is.na(values) == F,is.na(station) == F,is.na(year) == F,is.na(month) == F)

# Remove duplicated rows
names(covrs)
covrs <- covrs[!duplicated(covrs[, c(3,4,6,7)]), ]

# Check for duplicates in all but value
any(duplicated(covrs[,c(3,4,6,7)])) # station, year, month, variable

#covars <- covars%>% filter(duplicated(covars[,1:4]) == F) ### CHECK LATER
head(covrs)

# Make data wide and calculate rain

covs.mo <- covrs %>% unique() %>% 
  mutate(variable2 = case_when(variable %in% c("ppt", "pt", "ppt-add.csv", "ppt-add2.csv", "ppt.nc-ad.csv") ~ "ppt",
                                     variable %in% c("snowmelt", "snowmelt-add.csv", "snowmelt.nc-ad.csv") ~ "snowmelt",
                                     variable %in% c("temm", "temm-add.csv", "temm.nc-ad.csv") ~ "temm",
                                     variable %in% c("snowfall", "snowfall-add.csv", "snowfall.nc-ad.csv") ~ "snowfall")) %>% ungroup() %>% dplyr::select(-c(variable)) %>% 
  group_by(station, year, month, variable2) %>% filter(values == max(values, na.rm = T)) %>% unique() %>%
  pivot_wider(names_from = variable2, values_from = values)  %>%
     mutate(rain = ppt - snowfall,
             temp = temm-273.15)
head(covs.mo)

covs.mo$date <- as.Date(paste(covs.mo$year, covs.mo$month, "15", sep = "-"), format = '%Y-%m-%d')
co.mo <- covs.mo %>% dplyr::select(station, year, month, date, ppt, snowfall, snowmelt, rain, temp)

write.csv(co.mo, "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/All_ERA5_monthly_AK_CAN.csv", row.names = F)

co.mo <- read.csv("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/All_ERA5_monthly_AK_CAN.csv")

### Calculate yearly average
covs.yr <- covs.mo %>% group_by(station, year) %>%
  summarise(n = length(ppt), ppt = sum(ppt)*1000, snowfall = sum(snowfall)*1000, rain = sum(rain)*1000, temp = mean(temp), snowmelt = sum(snowmelt)*1000) %>% # Convert to mm
  filter(n == 12) %>% dplyr::select(-c(n))
head(covs.yr)

range <- covs.yr %>% pivot_longer(cols = 3:7) %>% group_by(station, name) %>% summarise(start = min(year), end = max(year), n = length(name))

# Export yearly
write.csv(covs.yr, "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/All_covariates_yearly.csv", row.names = F)


#### Check data completeness

# AK 60 years
ak.60 <- covs.mo %>% filter(year > 1966, station %in% c(15258000, 15276000, 15290000, 15292700,15484000,
                                                        15072000, 15514000, 15511000, 15515500, 15356000))
ak.60$date <- as.Date(paste(ak.60$year, ak.60$month, "15", sep = "-"))

ak.60 %>% ggplot(aes(date, station, col = is.na(ppt))) + geom_point()

#ak.60 %>% filter(station == 15258000) %>% ggplot(aes(date, temp)) + geom_point()

# Canada
covs.mo %>% filter(station == "10RA001") %>% ggplot(aes(date, station, col = is.na(rain))) + geom_point()
covs.mo %>% filter(station == "07FD001") %>% ggplot(aes(date, station, col = is.na(ppt))) + geom_point()
covs.mo %>% filter(station == "04DC001") %>% ggplot(aes(date, station, col = is.na(snowmelt))) + geom_point()
covs.mo %>% filter(station == "15024800") %>% ggplot(aes(date, station, col = is.na(snowfall))) + geom_point()

covs.mo %>% filter(station %in% unique(Qm$station.x),) %>% ggplot(aes(date,station, col = is.na(temp))) + geom_point()

covs.mo %>% filter(station == "09EA004") %>% ggplot(aes(year, station, col = is.na(values))) + geom_jitter()

covs.mo %>% filter(station== "09EA004") %>% ggplot(aes(as.Date(paste(year,month,'15',sep='-')),values)) + geom_point()

