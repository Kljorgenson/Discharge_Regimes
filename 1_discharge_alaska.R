### Download and sort all USGS discharge (cfs) data
library(dataRetrieval)
library(dplyr)
library(ggplot2)
library(data.table)
library(zoo)


### Download and format data for all USGS river gauges
# Select all Alaska river gauges
stations <- read.csv("Input_data/USGS_stations.csv")
stations$station <- as.character(stations$station)

# Download discharge data from USGS for all stations
USGS.Q <- readNWISdv(siteNumbers = stations$station,
                      parameterCd = "00060") # USGS code for discharge (cfs)

# Rename columns from USGS names
USGS.Q <- renameNWISColumns(USGS.Q)

# Join river names
names(USGS.Q) <- c("agency", "station", "date", "cfs_0", "flow_cd")
Q.usgs0 <- left_join(USGS.Q, stations) %>% mutate(cfs = ifelse(flow_cd %in% c("A", "A e", "A R", "A [4]"), cfs_0, NA)) # Filter for quality data

# Format dates
Q.usgs0$year <- format(Q.usgs0$date, "%Y")
Q.usgs0$date <- as.Date(Q.usgs0$date)


## Fill in NA rows for all missing daily time steps
# First remove all NA to remove trailing and leading NAs
Q.usgs0 <- Q.usgs0 %>% filter(is.na(cfs) == F)

# Fill in all missing daily timesteps
Q.usgs <- Q.usgs0 %>% filter(is.na(date) == F) %>%
  group_by(name) %>% tidyr::complete(date = seq(min(date),
                               max(date), by = "day"))

# Start and end dates by site
Q.range <- Q.usgs %>% group_by(name) %>% summarise(start = min(date), end = max(date))





####### Discharge data available over long-term time period
## Select stations with discharge data that starts before 1980 and ends before 2015
Q.long <- Q.usgs %>% group_by(name) %>% filter(date > "1966-01-01") %>%
  summarise(start = min(date), end = max(date)) %>% 
  filter(start < "1980-01-01" & end > "2015-01-01")

# Filter for the stations selected above
Q.usgs.l <- Q.usgs %>% filter(name %in% unique(Q.long$name),
                              date > "1966-01-01")

### Remove rivers not fitting selection criteria
## Total missing daily values is < 30%
miss <- Q.usgs.l %>% group_by(name) %>% summarise(n = length(cfs),
                                                            n.na = sum(is.na(cfs)),
                                                            per.na = n.na/n) %>% filter(per.na > .30)

## Total missing monthly values is < 30%
# Take monthly averages for months missing data for <= 10 days
Q.usgs.l$month <- month(Q.usgs.l$date)
Q.usgs.l$year <- year(Q.usgs.l$date)
Q.usgs.l.m <- Q.usgs.l %>% group_by(name, year, month) %>% mutate(miss = sum(is.na(cfs))) %>%
  summarise(cfs = ifelse(max(miss)<=10, mean(cfs, na.rm = T), NA))

# Select sites missing data for > 30% of all months
miss.m <- Q.usgs.l.m %>% group_by(name, month) %>% summarise(n = length(cfs),
                                                                       n.na = sum(is.na(cfs)),
                                                                       per.na = n.na/n) %>% filter(per.na > .30)

## Not missing more than 4 years in a row for any individual month
gaps.m <- Q.usgs.l.m %>% filter(cfs %in% NA) %>% group_by(name, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)

## Not containing consecutive gaps > 4 months
Q.usgs.l.m$date <- as.Date(paste(Q.usgs.l.m$year, Q.usgs.l.m$month, "15", sep = "-"))
gaps <- Q.usgs.l.m %>% filter(cfs %in% NA) %>% group_by(name) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4)


## View all stations that meet the criteria
stats <- Q.usgs.l.m %>% filter(!name %in% c(miss$name, miss.m$name, gaps.m$name, gaps$name)) %>% ungroup() %>% select(name) %>% unique()

## Select stations to include. Stikine not included because of later start date
Q.usgs.l <- Q.usgs.l %>% filter(date >= "1967-01-01", name %in% c('KUPARUK', 'CHENA_1', 'FISH_SE', 'KENAI_1', 'L_CHENA', 'SALCHA', 'SHIP', 'TALKEETNA', 'TANANA_3', 'YUKON_1', 'YUKON_4', 'L_SUSITNA'),
                                date >= '1974-05-28' | name != 'KUPARUK') # Remove poor early data for Kuparuk

# Export daily discharge data
write.csv(Q.usgs.l, "Data/Q60_daily.csv", row.names = F)




####### Discharge data over recent time period (2000-2022)
## Select stations with discharge data that starts before 2000 and ends before 2015
# Select stations
Q.20 <- Q.usgs %>% group_by(name) %>%
  summarise(start = min(date), end = max(date)) %>% 
  filter(start < "2005-01-01" & end > "2015-01-01")

# Filter discharge data for the stations selected above
Q.usgs.20 <- Q.usgs %>% filter(name %in% unique(Q.20$name),date >= "2000-01-01") %>%
  filter(date > '2001-03-31' | !name == 'YUKON_3', # Trim data for stations with bad data just at the beginning of the timeseries
         date > '2001-05-16' | !name == 'SAWMILL_2',
         date > '2001-04-30' | !name == 'WILLOW',
         date > '2003-09-30' | !name == 'TAIYA',
         date > '2000-05-19' | !name == 'MATANUSKA') 


### Remove rivers not fitting selection criteria
## Total missing daily values is < 30%
miss <- Q.usgs.20 %>% group_by(name) %>% summarise(n = length(cfs),
                                                  n.na = sum(is.na(cfs)),
                                                  per.na = n.na/n) %>% filter(per.na > .30)

## Total missing monthly values is < 30%
# Take monthly averages for months missing data for <= 10 days
Q.usgs.20$month <- month(Q.usgs.20$date)
Q.usgs.20$year <- year(Q.usgs.20$date)
Q.usgs.20.m <- Q.usgs.20 %>% group_by(name, year, month) %>% mutate(miss = sum(is.na(cfs))) %>%
  summarise(cfs = ifelse(max(miss)<=10, mean(cfs, na.rm = T), NA))
miss.m <- Q.usgs.20.m %>% group_by(name, month) %>% summarise(n = length(cfs),
                                                             n.na = sum(is.na(cfs)),
                                                             per.na = n.na/n) %>% filter(per.na > .30)


## Not missing more than 4 years in a row for any individual month
gaps.m <- Q.usgs.20.m %>% filter(cfs %in% NA) %>% group_by(name, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)

## Not containing consecutive gaps > 4 months
Q.usgs.20.m$date <- as.Date(paste(Q.usgs.20.m$year, Q.usgs.20.m$month, "15", sep = "-"))
gaps <- Q.usgs.20.m %>% filter(cfs %in% NA) %>% group_by(name) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4) %>% filter(gap_start < "2020-01-01")
gaps <- gaps %>% filter(!name %in% c("KUSKOKWIM_1", "LEMON", "SUSITNA"))


## Select discharge data for all stations that meet the criteria
Q.usgs_20 <- Q.usgs.20 %>% filter(date >= "2000-01-01", !name %in% c(miss.m$name, miss$name, gaps.m$name, gaps$name),
                                      !name %in% c("TANANA_2", "KENAI_2", "KENAI_3", "CHENA_4")) # Remove stations with catchment area overlap > 50%
Q.usgs_20 <- merge(Q.usgs_20, stations)


# Export daily discharge data
write.csv(Q.usgs_20, "Output_data/Q20_daily.csv", row.names = F)
