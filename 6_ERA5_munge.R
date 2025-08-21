#### Covariate organization
# Compile data from .csv file output from clipping ERA5 data to watersheds
# Some parts of code must be adjusted for different climate variables
library(tidyverse)
library(tidyhydat)
library(ggplot2)



########### Combine all hourly data for each variable and calculate daily and monthly averages
# Combine files for each year 
filelist = list.files(path = "Output_data/", pattern = "*snowmelt*", # Adjust variable name
                      full.names = TRUE, recursive = T)
df_list <- lapply(filelist, read.csv) 
names(df_list) <- basename( filelist )

df_list_named <- map2(df_list, names(df_list), ~mutate(.x, source = .y))

# Make all columns characters
df_list_named <- lapply(df_list_named, function(x) {x[] <- lapply(x, as.character); x})
cov.raw <- df_list_named %>% reduce(full_join)

# Extract variable from file name
cov.raw$variable <- substr(cov.raw$source, start = 6, stop = 200)
cov.raw$variable <- gsub(".([0-9]+)|.nc.csv","",cov.raw$variable)

# Remove duplicated time steps
covs.hrly <- cov.raw %>% dplyr::select(date_UTC, time, station, variable, value) 

# Check for duplicates in all but value
any(duplicated(covs.hrly[,1:4]))

names(covs.hrly) <- c("date_UTC", "hour_UTC", "station", "variable", "value")
covs.hr <- covs.hrly

### Take monthly averages
# Based on timezone of discharge data

# Match stations and UTC-X
meta.ak <- read.csv("Data/river_metadata.csv") %>% filter(NameNom == "") %>% 
  rename(station = StationNum) %>%
  dplyr::select(station) %>% unique()
meta.ak$UTC_minus <- 9 # All Alaska stations in same time zone
meta.ak$station <- as.character(meta.ak$station)
meta.can <- hy_stations() %>% filter(STATION_NUMBER %in% covs.hr$station) %>%
  rename(station = STATION_NUMBER) %>%
  mutate(UTC_minus = case_when(PROV_TERR_STATE_LOC %in% c("YT", "BC") ~ 8, # Adjust time zones for Canadian provinces
                               PROV_TERR_STATE_LOC %in% c("NT", "AB", "NU") ~ 7,
                               PROV_TERR_STATE_LOC %in% c("MB", "SK") ~ 6,
                               PROV_TERR_STATE_LOC %in% c("ON", "QC") ~ 5,
                               PROV_TERR_STATE_LOC %in% c("NL") ~ 4
  )) %>% dplyr::select(station, PROV_TERR_STATE_LOC, UTC_minus)

meta <- full_join(meta.ak, meta.can)

# Adjust datetime from UTC to individual time zones
covs.hr$station <- as.character(covs.hr$station)
covs.hr <- left_join(covs.hr, meta)
covs.hr$time <- ifelse(nchar(covs.hr$hour_UTC) == 1, paste(0,covs.hr$hour_UTC, sep = ""), covs.hr$hour_UTC)

covs.hr$datetime_UTC <- as.POSIXct(paste(covs.hr$date_UTC, " ", covs.hr$time, ":00:00", sep = ""), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
covs.hr$datetime_adj <- covs.hr$datetime_UTC-covs.hr$UTC_minus* 3600
covs.hr$date_adj <- as.Date(covs.hr$datetime_adj)
covs.hr$year <- year(covs.hr$datetime_adj)
covs.hr$month <- month(covs.hr$datetime_adj)

### Calculate monthly average
covs.hr$value <- as.numeric(covs.hr$value)
covs.mo <- covs.hr %>% group_by(station, year, month, variable) %>% 
  summarise(n = length(variable), value = if_else(first(variable) == 'temp', mean(value), sum(value))) %>% 
  filter(n > 600) %>% dplyr::select(-c(n)) 

# Export monthly data
write.csv(covs.mo, "Output_data/Snowmelt_monthly.csv", row.names = F)

# Daily
covs.day <- covs.hr %>% mutate(jday = yday(datetime_adj)) %>% group_by(station, year, month, jday, variable) %>%
  summarise(n = length(variable), value = if_else(first(variable) %in% c("temp", "temm", "temm-add.csv", "temm.nc-ad.csv"), mean(value), sum(value))) %>% 
  filter(n == 24) %>% dplyr::select(-c(n))

write.csv(covs.day, "Output_data/Snowmelt_daily.csv", row.names = F)




########### Combine monthly data for all climate variables
filelist = list.files(path = "Output_data/", pattern = "*monthly*", full.names = TRUE)
filelist

covars0 <- lapply(filelist, read.csv)
covars <- lapply(covars0, function(x) {x[] <- lapply(x, as.character); x}) %>% reduce(full_join) %>% unique() 
covars$value <- as.numeric(covars$value)
covrs <- covars %>% filter(is.na(value) == F,is.na(station) == F,is.na(year) == F,is.na(month) == F)

# Remove duplicated rows
covrs <- covrs[!duplicated(covrs[,1:4]), ]

# Check for duplicates in all but value
any(duplicated(covrs[,1:4])) # station, year, month, variable

# Make data wide and calculate rain
covs.mo <- covrs %>% unique() %>%
  pivot_wider(names_from = variable, values_from = value)  %>%
  mutate(ppt2 = ifelse(ppt < 0, 0, ppt), # Replace negative precipitation with zero
         snow2 = ifelse(snowfall < 0, 0, snowfall), # Replace negative snowfall with zero
         snow3 = ifelse(ppt2 <= 0, 0, snow2), # If total precipitation is zero, make snowfall also zero
         snow4 = ifelse(snow3 > ppt2, ppt2, snow3), # If snowfall is greater than total precipitation, make equal to total precipitation
         rain = ppt2 - snow4, # Calculate rain
         temp = temp-273.15) %>% select(station,year,month,ppt2,snow4,snowmelt,rain,temp) %>%
  rename(snowfall = snow4, ppt = ppt2)

# Add date
covs.mo$date <- as.Date(paste(covs.mo$year, covs.mo$month, "15", sep = "-"), format = '%Y-%m-%d')

# Export all combined monthly climate data
write.csv(co.mo, "Output_data/All_ERA5_monthly_AK_CAN.csv", row.names = F)
