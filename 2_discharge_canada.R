### Download and sort all Water Survey of Canada discharge (m3/s) data
library(tidyhydat)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(lubridate)
library(tidyr)
library(readxl)
library(tidyverse)

### Download HYDAT database only once
#download_hydat()

## All station metadata
stations <- hy_stations()



######### Discharge data available over long-term time period
## Starting list of all stations within study region
stn_range <- hy_stn_data_range() %>% filter(DATA_TYPE == "Q") %>% filter(Year_to >= "2015", Year_from <= "1976")
stations.60.all <- hy_stations() %>% filter(STATION_NUMBER %in% stn_range$STATION_NUMBER, LATITUDE > 51)

# Clip to study region
northern_region <- st_read("Input_data/Northern_region_shape/Northern_region.shp")
stations.60.all <- st_as_sf(stations.60.all, coords = c("LONGITUDE","LATITUDE"), crs = "NAD83")
stations.60.all <- st_transform(stations.60.all, crs = st_crs(northern_region))
stations.60.region <- st_filter(stations.60.all, northern_region)
stations.60.region <- left_join(stations.60.region, stn_range) %>% filter(HYD_STATUS == "ACTIVE" | Year_to >= "2015")

# Load discharge data from HYDAT
Q.ca <- hy_daily(station_number = stations.60.region$STATION_NUMBER)

# Remove stage data and reformat
q.ca.60.all <- Q.ca %>% filter(Parameter == "Flow")
q.ca.60.all$Date <- as.Date(q.ca.60.all$Date)
q.ca.60.all$year <- format(q.ca.60.all$Date, "%Y")
names(q.ca.60.all) <- c('station', 'date', 'parameter', 'discharge', 'symbol', 'year')

## Add preliminary data from some stations - acquired by request from Water Survey of Canada
file.list <- list.files(path = "Input_data/Preliminary_Q_CAN/", pattern = ".xlsx", full.names = T)
q.prelim <- file.list %>%
  map_dfr(~read_excel(.x, skip = 14))%>% select('ISO 8601 UTC', Value, station)
names(q.prelim) <- c('date', 'discharge', 'station')
q.prelim$date <- as.Date(q.prelim$date)
q.prelim.daily <- q.prelim %>% group_by(date, station) %>% summarise(discharge = mean(discharge, rm.na = T)) %>% filter(station != "10MA003")
head(q.prelim.daily)

# Plot to check quality
q.prelim.daily %>% group_by(station) %>% ggplot(aes(date, discharge)) + geom_point() + facet_wrap(~station, scales = "free_y")

# Join prelim data with all data
q.ca.60.all <- full_join(q.ca.60.all, q.prelim.daily)


## Fill all missing daily time steps: remove all NA to remove trailing and leading NAs and then fill in all missing daily timesteps
q.ca.60.all <- q.ca.60.all %>% filter(is.na(discharge) == F) %>% 
  group_by(station) %>% complete(date = seq(min(date),max(date), by = "day"))
                                                           
# Join station metadata
stations2 <- stations[,c(1,2,3,7,8,9)]
names(stations2) <- c('station', 'name_full', 'location', 'latitude', 'longitude', 'area')
q.60.ca <- left_join(q.ca.60.all, stations2)
length(unique(q.60.ca$station))
unique(q.60.ca$name_full)

# Range
Q.range <- q.60.ca %>% group_by(station, name_full) %>% summarise(start = min(date), end = max(date))

# Select potential data
q.60 <- q.60.ca %>% filter(date >= "1967-01-01", ) %>%
  filter(date >= '1970-03-15' | !station %in% c('07BF002', '07HA003'), # Trim data for stations with bad data just at the beginning
         date >= '1972-01-01' | !station == '06BD001',
         date >= '1972-09-01' | !station == '09BC004',
         date >= '1970-09-17' | !station == '08DA005',
         date > '1972-08-29' | !station == '06BA002',
         date > '1968-09-29' | !station == '05UF004',
         date > '1976-05-26' | !station == '06DA005') 


# Add missing river names
q.60[q.60$station == "10NC001", 7] <- "ANDERSON RIVER BELOW CARNWATH RIVER"
q.60[q.60$station == "09AB001", 7] <- 'YUKON RIVER AT WHITEHORSE'
sort(unique(q.60$name_full))


### Assess candidate time series for inclusion in the analysis using the following criteria: 1) <30% of days missing over the entire record, 2) no gaps > 3 months, 3) not missing more than four consecutive years for any month, and 4) nested stations comprised < 50% of the catchment area of downstream stations. 
## Total NA < 30%
miss <- q.60 %>% group_by(name_full, station) %>% summarise(n = length(discharge),
                                            n.na = sum(is.na(discharge)),
                                            per.na = n.na/n) %>% filter(per.na > .30)
miss

## By month NA > 30%
q.60$month <- month(q.60$date)
q.60$year <- year(q.60$date)
q.60.m <- q.60 %>% group_by(name_full, station, year, month) %>% mutate(miss = sum(is.na(discharge))) %>%
  summarise(discharge = ifelse(max(miss)<=10, mean(discharge, na.rm = T), NA)) # <= 10 days missing each month
miss.m <- q.60.m %>% group_by(name_full, station, month) %>% summarise(n = length(discharge),
                                           n.na = sum(is.na(discharge)),
                                           per.na = n.na/n) %>% filter(per.na > .30)
sort(unique(miss.m$name_full))
sort(unique(miss$station))

## Not missing more than 4 years in a row for any month
gaps.m <- q.60.m %>% filter(discharge %in% NA) %>% group_by(name_full, station, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name_full, station, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)

# Pardon key rivers
gaps.m <- gaps.m %>% filter(!station %in% c('10LA002', '10RC001', '10RA001', '10KA001', '06LC001', '10MC002'))

sort(unique(gaps.m$station))

## Identify data gaps greater than 4 months
q.60.m$date <- as.Date(paste(q.60.m$year, q.60.m$month, "15", sep = "-"))
gaps <- q.60.m %>% filter(discharge %in% NA) %>% group_by(name_full, station) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name_full, station, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4)

# Pardon key rivers
gaps <- gaps %>% filter(!station %in% c('10LA002', '10KA001', '10LC014', '10MC002', "06LC001", "10PB001", "10LC007", "09EA004", "10LC003"),
                        !name_full %in% "MACKENZIE RIVER AT FORT SIMPSON")
sort(unique(gaps$name_full))

## Select final data by remove stations not meeting criteria: >50% catchment outside study region, damed rivers, and catchments >50% of downstream catchments
q_60 <- q.60 %>% filter(!station %in% c(gaps$station, gaps.m$station, miss.m$station, miss$station)) %>%
            filter(!station %in% c("05UE005", "07BJ001", "07GH002", "07GE001", "07GJ001", # >50% outside region
            "05UE005", "07BJ001", "07GH002", "07GE001", "07GJ001", "04FC001", "05UF006", "05UF007", "07FB006", "07FB009","07GD004", "08EC013", # >50% outside region
        "07EF001", "07FA004", "07FD010", "07HA001", "07FD002", # Contain large dams: Peace River
        "07QD007", # Taltson River Dam
        "06EA006", "06EB004", "06FD001", "06EA002", "06FB001", # Remove all Churchhill except farthest upstream (above dams)
        "06DD002", # Whitesands Dam on Reindeer River
        "05TG001", # Wuskwatim Dam on Burntwood River: Except Durocher et al. left it in.                   
          "09CD001",  '10ED001', "10KA001", # Remove sites on same rivers
            "10BE001", "10LC002", "09AB001")) # Liard at lower crossing, MACKENZIE RIVER (EAST CHANNEL) AT INUVIK, YUKON RIVER AT WHITEHORSE
            

q_60 <- q.60 %>% filter(station %in% q.60_m$station)                                

Q.range <- q_60 %>% group_by(station, name_full) %>% summarise(start = min(date), end = max(date))

# Save daily discharge data
write.csv(q_60, "Output_data/Q60_daily_can.csv", row.names = F)

                



                
############ All sites for recent period (2000-2022)
## Stations with recent ~20 years
stn_range2 <- hy_stn_data_range() %>% filter(DATA_TYPE == "Q") %>% filter(Year_to >= "2010", Year_from <= "2007", RECORD_LENGTH >= 10)
stn_range2
stations.20.all <- hy_stations() %>% filter(STATION_NUMBER %in% stn_range2$STATION_NUMBER)

# Clip to study region
stations.20.all <- st_as_sf(stations.20.all, coords = c("LONGITUDE","LATITUDE"), crs = "NAD83")
stations.20.all <- st_transform(stations.20.all, crs = st_crs(northern_region))
stations.20.region <- st_filter(stations.20.all, northern_region)

# Load discharge data from HYDAT
Q.ca.20 <- hy_daily(station_number = stations.20.region$STATION_NUMBER)
Q.ca.20$Date <- as.Date(Q.ca.20$Date)
Q.ca.20$year <- format(Q.ca.20$Date, "%Y")
head(Q.ca.20)

# Remove stage data
Q.ca.20.all <- Q.ca.20 %>% filter(Parameter == "Flow", Date >= "2000-01-01")
names(Q.ca.20.all) <- c('station', 'date', 'parameter', 'discharge', 'symbol', 'year')

## Add preliminary data from some stations
file.list <- list.files(path = "Input_data/Preliminary_Q_CAN/", pattern = ".xlsx", full.names = T)
q.prelim <- file.list %>% map_dfr(~read_excel(.x, skip = 14)) %>% select('ISO 8601 UTC', Value, station)
names(q.prelim) <- c('datetime_UTC', 'discharge', 'station')
q.prelim$date <- as.Date(q.prelim$datetime_UTC)
q.prelim.daily <- q.prelim %>% group_by(date, station) %>% summarise(discharge = mean(discharge, rm.na = T)) %>% 
  filter(date >= "2000-01-01")
q.prelim.daily$version <- 'prelim'
head(q.prelim.daily)

# Plot to check quality
q.prelim.daily %>% group_by(station) %>% ggplot(aes(date, discharge)) + geom_point() + facet_wrap(~station, scales = "free_y")

# Join prelim data with all data
Q.ca.20.all <- full_join(Q.ca.20.all, q.prelim.daily) %>% group_by(station, date) %>% 
  filter(length(discharge) == 1 | !version %in% 'prelim')

# Add additional data from eastern Canada
adds <- read.csv("Input_data/Additional_Q_CAN/daily_03NE011_03NE012.csv", skip = 1) %>% filter(PARAM == 1) %>% select(ID, Date, Value)
names(adds) <- c("station", "date", "discharge")
adds$date <- as.Date(adds$date)
Q.ca.20.all <- Q.ca.20.all %>% filter(!station %in% c('03NE011', '03NE012')) %>%
               full_join(adds)

# Remove all NA to remove trailing and leading NAs, then fill all missing daily time steps
Q.ca.20.all <- Q.ca.20.all %>% filter(is.na(discharge) == F, date >= "2000-01-01") %>% group_by(station) %>% 
               complete(date = seq(min(date), max(date), by = "day"))

# Join station metadata
stations20 <- stations[,c(1,2,3,7,8,9)]
names(stations20) <- c('station', 'name_full', 'location', 'latitude', 'longitude', 'area')
q.ca.20 <- left_join(Q.ca.20.all, stations20) %>% filter(date >= "2000-01-01")
length(unique(q.ca.20$station))
unique(q.ca.20$name_full)

# Filter for start and end dates
Q.range20 <- q.ca.20 %>% group_by(station, name_full) %>% summarise(start = min(date), end = max(date))
keep <- Q.range20 %>% filter(start < "2006-01-01" & end > '2015-01-01')
q.ca.20 <- q.ca.20 %>% filter(station %in% keep$station)

### Assess candidate time series for inclusion in the analysis using the following criteria: 1) <30% of days missing over the entire record, 2) no gaps > 3 months, 3) not missing more than four consecutive years for any month, and 4) nested stations comprised < 50% of the catchment area of downstream stations. 
# Overall missing < 30%
miss.20 <- q.ca.20 %>% group_by(name_full, station) %>% summarise(n = length(discharge),
                                           n.na = sum(is.na(discharge)),
                                           per.na = n.na/n) %>% filter(per.na > .30)
miss.20

# By month missing < 30%
q.ca.20$month <- month(q.ca.20$date)
q.ca.20$year <- year(q.ca.20$date)
q.20.m <- q.ca.20 %>% group_by(name_full, station, year, month) %>% mutate(miss = sum(is.na(discharge))) %>%
  summarise(discharge = ifelse(max(miss)<=10, mean(discharge, na.rm = T), NA)) # <= 10 days missing per month
miss.20.m <- q.20.m %>% group_by(name_full, station, month) %>% summarise(n = length(discharge),
                                                    n.na = sum(is.na(discharge)),
                                                    per.na = n.na/n) %>% filter(per.na > .30)

# Not missing more than 4 years in a row for any month
gaps.20.m <- q.20.m %>% filter(discharge %in% NA) %>% group_by(name_full, station, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name_full, station, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)
gaps.20.m %>% filter(station %in% q_60$station) %>% select(station) %>% unique()

# Identify data gaps greater than 4 months
q.20.m$date <- as.Date(paste(q.20.m$year, q.20.m$month, "15", sep = "-"))
gaps <- q.20.m %>% filter(discharge %in% NA) %>% group_by(name_full, station) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name_full, station, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4)
# Pardon key stations
gaps <- gaps %>% filter(!station %in% c("10LA002", "06LC001", "10GC001", "10LC003", '10RA001'))

# Additional stations to remove because of dams, <50% within study region, or >50% overlap with other stations
rm <- c("05UE005", "07BJ001", "07GH002", "07GE001", "07GJ001", "04FC001", "05UF006", "05UF007", "07FB006", "07FB009","07GD004", "08EC013", "05UB009", "04CA002", "04DA002", "05UB009", # These are <50% inside region
        "07EF001", "07FA004", "07FD010", "07HA001", "07FD002", # Contain large dams: Peace River
          "07QD007", # Taltson River Dam
        "06EA006", "06EB004", "06FD001", "06EA002", "06FB001", # Remove all Churchhill except farthest upstream (above dams)
        "05UE005", "07BJ001", "07GH002", "07GE001", "07GJ001", "04FC001", "05UF006", "05UF007", "07FB006", "07FB009","07GD004", "08EC013", # These are <50% inside region
        "07EF001", "07FA004", "07FD010", "07HA001", "07FD002", # Contain large dams: Peace River
        "07QD007", # Taltson River Dam
        "06EA006", "06EB004", "06FD001", "06EA002", "06FB001", # Remove all Churchhill except farthest upstream (above dams)
        "06DD002", # Whitesands Dam on Reindeer River
        "05TG001",
        "06DD002", # Whitesands Dam on Reindeer River
        "05TG001", # Wuskwatim Dam on Burntwood River
        "09EB001", # Yukon at Dawson
        "07ED001","10FB006","10LC002","08AA010", "03NE001","07AA008","08AB002","10ED001", "10FB006", "10LC002", "07SA008","09AB001","09CD001", # >50% of downstream watershed
        "08DB014", # Only goes through 2017
        "10BE001", "10ED001", # Liard
        "08BB005") # AK

## Finalize data selection
q_20 <- q.ca.20 %>% filter(!station %in% c(miss.20$station, gaps.20.m$station, gaps$station, miss.20.m$station, rm))

# Make sure no missing any stations with historical data
q_60 %>% ungroup() %>% select(station, name_full) %>% filter(!station %in% q_20$station) %>% unique()

# Join missing names
stats <- hy_stations() %>% select(STATION_NUMBER, STATION_NAME) %>% rename(station = STATION_NUMBER, name_full = STATION_NAME)
q_20.1 <- left_join(q_20, stats)

# Check range
q_20.1 %>% group_by(station) %>% summarise(start = min(date), end = max(date)) %>% View()

write.csv(q_20.1, "Output_data/Q20_daily_can.csv")

