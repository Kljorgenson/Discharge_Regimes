# Packages
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyhydat)
library(tidyverse)


#### Metadata and catchment characteristics
### Ecoregions
ecoreg <- read.csv("Input_data/Ecoregions.csv") %>% rename(station = StationNum) %>% rename(ecoregion1 = NA_L1NAME) %>%
  mutate(ecoregion1 = paste(substr(ecoregion1, 1, 1), tolower(substr(ecoregion1, 2, nchar(ecoregion1))),  sep=""),
         ecoregion1 = case_when(ecoregion1 == "Hudson plain" ~ "Hudson Plain", # Capitalize ecoregion names
                                ecoregion1== "Northern forests" ~ "Northern Forests",
                                ecoregion1== "Northwestern forested mountains" ~ "Northwestern Forested Mountains",
                                ecoregion1== "Marine west coast forest" ~ "Marine West Coast Forest",
                                ecoregion1 %in% c("Tundra", "Taiga") ~ ecoregion1))

head(ecoreg)
tail(ecoreg)

# Assign the ecoregion that covers the largest area within the catchment
eco.largest <- ecoreg %>% group_by(station) %>% rename(ecoreg1_percent = PERCENTAGE) %>% 
  filter(ecoreg1_percent == max(ecoreg1_percent)) %>%
  # Add site in Hudson Plain ecoregion to Northern Forests ecoregion
  mutate(across(ecoregion1, ~ifelse(ecoregion1 == "Hudson Plain", "Northern Forests", .)))



### Discharge metadata
# List of all Alaskan stations
stations <- read.csv("Input_data/USGS_stations.csv")
stations$station <- as.character(stations$station)

# Catchment area and latitude and longitude of catchment centroid from spatial analysis in ArcGIS
met <- read.csv("Input_data/river_metadata.csv") %>% rename(station = StationNum, area_km2 = Area_km2, latitude.c = Latitude_c, longitude.c = Longitude_c) %>%
  left_join(stations) %>% mutate(name = ifelse(is.na(name), NameNom,name)) %>%
  select(station, name, area_km2, latitude.c, longitude.c) %>%
  mutate(area_km2_log = log(area_km2))
tail(met)


### Broad permafrost groupings
# Percent of catchment area in each permafrost category from spatial analysis in ArcGIS
perm <- read.csv("Input_data/Permafrost_percent.csv") %>% dplyr::select(StationNum, EXTENT, AREA, PERCENTAGE)
names(perm) <- c("station", "permafrost_class", "permafrost_area", "permafrost_percent")

# Add an 'Unfrozen' category for catchments without any unfrozen area, and assign the area to zero
non.perm <- perm %>% group_by(station) %>% filter(sum(permafrost_percent) < 100) %>% 
  summarise(permafrost_percent = 100-sum(permafrost_percent))

non.perm$permafrost_class <- "Unfrozen"
non.perm$permafrost_area <- 0

# Assign the permafrost category that covers the largest area within the catchment
perm <- full_join(perm, non.perm) %>% group_by(station) %>% filter(permafrost_percent == max(permafrost_percent)) %>%
  mutate(across(permafrost_class, ~ifelse(permafrost_class == "I", "S", .))) # Combine I and S classes


### Glacier area from spatial analysis in ArcGIS
glaciers <- read.csv("Input_data/Glacier_area.csv") %>% dplyr::select(StationNum, AREA, PERCENTAGE) %>%
  rename(station = StationNum, glacier_area_km2 = AREA, glacier_percent = PERCENTAGE)

head(glaciers)


### Climate variables from ERA5 data
## Select start and end dates to clip the climate data data to the same date range and the discharge anomaly data
# Historical period
range.60 <- read.csv("Output_data/DFFT_metrics_1970-1992_ak_can.csv")%>% 
  rename(name = site) %>% left_join(met) %>% select(name, station, start.date, end.date)

# Recent period
range.20 <- read.csv("Output_data/DFFT_metrics_20yrs_ak_can.csv") %>% 
  rename(name = site) %>% left_join(met) %>% select(name, station, start.date, end.date)

## Read in the ERA5 data, clip it to the date range, take the sum (precipitation, snowfall, rain, snowmelt) or mean (temperature), 
# convert from meters to millimeters for precipitation variables, and take the average during the period for each catchment.
covs.an.60 <- read.csv("Output_data/All_ERA5_monthly_AK_CAN.csv") %>%  
  left_join(range.60) %>% mutate(dat = as.Date(date)) %>% filter(date >= start.date & date <= end.date) %>% 
  filter(year >= 1970 & year <= 1992, min(year, na.rm = T) < 1980) %>% group_by(station, year) %>%
  summarise(n = length(ppt), ppt = sum(ppt)*1000, snowfall = sum(snowfall)*1000, rain = sum(rain)*1000, temp = mean(temp), snowmelt = sum(snowmelt)*1000) %>% # Convert to mm
  filter(n == 12) %>% dplyr::select(-c(n)) %>% summarise(ppt.60 = mean(ppt, na.rm = T),rain.60 = mean(rain, na.rm = T),snowfall.60 = mean(snowfall, na.rm = T),temp.60 = mean(temp, na.rm = T),snowmelt.60 = mean(snowmelt, na.rm = T))

covs.an.20 <- read.csv("Output_data/All_ERA5_monthly_AK_CAN.csv") %>% 
  left_join(range.20) %>% mutate(dat = as.Date(date)) %>% filter(date >= start.date & date <= end.date) %>% 
  filter(year >= 2000) %>% group_by(station, year) %>%
  summarise(n = length(ppt), ppt = sum(ppt)*1000, snowfall = sum(snowfall)*1000, rain = sum(rain)*1000, temp = mean(temp), snowmelt = sum(snowmelt)*1000) %>% # Convert to mm
  filter(n == 12) %>% dplyr::select(-c(n)) %>% summarise(ppt.20 = mean(ppt, na.rm = T),rain.20 = mean(rain, na.rm = T),snowfall.20 = mean(snowfall, na.rm = T),temp.20 = mean(temp, na.rm = T),snowmelt.20 = mean(snowmelt, na.rm = T))


any(is.na(covs.an.60$ppt.60))

### Calculate snowmelt timing and precipitation seasonality index
## Precipitation seasonality index for June-October
# Clip ends of date ranges to exclude years without full June-October season
range20 <- range.20 %>% mutate(jday.start = yday(start.date), jday.end = yday(end.date)) %>%
  mutate(start.20 = ifelse(jday.start <= 152, year(start.date), year(start.date) + 1),
         end.20 = ifelse(jday.end >= 304, year(end.date), year(end.date) - 1))%>% select(station,start.20,end.20)
range60 <- range.60 %>% mutate(jday.start = yday(start.date), jday.end = yday(end.date)) %>%
  mutate(start.60 = ifelse(jday.start <= 152, year(start.date), year(start.date) + 1),
         end.60 = ifelse(jday.end >= 304, year(end.date), year(end.date) - 1))%>% select(station,start.60,end.60)

range <- full_join(range20,range60) %>% select(station,start.20,end.20,start.60,end.60)

# Read in precipitation data
rain <- read.csv("Output_data/Precipitation_daily.csv") %>%
  filter(is.na(value) == F) 

# Clip data to date range and calculate seasonality index
si <- rain %>% filter(year >= '1970' & year <= '2020') %>% filter(year <= '1992' | year >= '2000') %>% 
  filter(jday >= 152 & jday <= 304) %>% # June-Oct
  mutate(period = ifelse(year > '1992', 'SI.20', 'SI.60')) %>% # Assign period
  group_by(station,year) %>% mutate(sum = sum(value)) %>% mutate(si.part = 1/sum*abs(value-(sum/(304-151)))) %>%
  full_join(range) %>%
  filter(ifelse(period == 'SI.20', year >= start.20 & year <= end.20, year >= start.60 & year <= end.60)) %>%
  group_by(station, period) %>% summarise(SI = sum(si.part)) %>%
  pivot_wider(values_from = SI, names_from = period)
head(si)

## Snowmelt timing, calculated as the centroid of the snowmelt peak occurs (mean day weighted by the snowmelt value) 
# Clip ends of date ranges to exclude years without March-June
range_20 <- range.20 %>% mutate(jday.start = yday(start.date), jday.end = yday(end.date)) %>%
  mutate(start.20 = ifelse(jday.start <= 60, year(start.date), year(start.date) + 1),
         end.20 = ifelse(jday.end >= 182, year(end.date), year(end.date) - 1))%>% select(station,start.20,end.20)
range_60 <- range.60 %>% mutate(jday.start = yday(start.date), jday.end = yday(end.date)) %>%
  mutate(start.60 = ifelse(jday.start <= 60, year(start.date), year(start.date) + 1),
         end.60 = ifelse(jday.end >= 182, year(end.date), year(end.date) - 1))%>% select(station,start.60,end.60)

range2 <- full_join(range_20,range_60) %>% select(station,start.20,end.20,start.60,end.60)

# Read in snowmelt data and calculate mean centroid as measure of snowmelt timing
melt <- read.csv("Output_data/Snowmelt_daily.csv") %>% 
  mutate(date = as.Date(jday - 1, origin = paste0(year, "-01-01"))) %>%
  filter(jday >= 60 & jday < 182, value > 0) %>% # Filter for March-June to exclude snowmelt following early fall storms
  group_by(station,year) %>% 
  summarise(center = weighted.mean(jday,value)) %>% # Calculate weighted mean
  filter(year >= '1970' & year <= '2020') %>% filter(year <= '1992' | year >= '2000') %>% 
  mutate(period = ifelse(year > '1992', 'snowmelt.t.20', 'snowmelt.t.60')) %>% # Assign period
  full_join(range2) %>%
  filter(ifelse(period == 'snowmelt.t.20', year >= start.20 & year <= end.20, year >= start.60 & year <= end.60)) %>% # Filter for years within date range
  group_by(period, station) %>%
  summarise(snowmelt.t = mean(center)) %>% # Take the average timing of snowmelt for each period
  pivot_wider(values_from = snowmelt.t, names_from = period)
head(melt)


### Join all metadata/catchment characteristics together
meta1 <- full_join(glaciers, perm) %>% full_join(met) %>% full_join(covs.an.60) %>% 
  full_join(covs.an.20) %>% full_join(si) %>% full_join(melt)
meta_full <- full_join(meta1, eco.largest) %>% filter(is.na(ecoregion1) == F)

# Replace NA with zero for glacier area and add unfrozen values for permafrost class
meta_full$glacier_area_km2 <- ifelse(is.na(meta_full$glacier_area_km2), 0, meta_full$glacier_area_km2)
meta_full$glacier_percent <- ifelse(is.na(meta_full$glacier_percent), 0, meta_full$glacier_percent)
meta_full$permafrost_class <- ifelse(is.na(meta_full$permafrost_class), "Unfrozen", meta_full$permafrost_class)
meta_full$permafrost_percent <- ifelse(is.na(meta_full$permafrost_percent), 100, meta_full$permafrost_percent)

# Classify as a glacial river if glacier % cover >= 1%
meta_full$glacial <- ifelse(meta_full$glacier_percent < 1, "non-glacial", "glacial")


# Assign short names to all rivers in dataset
nicknames <- read.csv("Input_data/Nicknames.csv")
meta_full <- left_join(meta_full, nicknames) 

# Only keep metadata for rivers selected for final analysis
meta_final <- meta_full %>% 
  filter(name %in% range.20$name) 

### Save river metadata
write.csv(meta_final, "Output_data/Metadata_AK_CAN.csv", row.names = F)


