### Download and sort all USGS cfs data
library(dataRetrieval)
library(dplyr)
library(ggplot2)
library(data.table)
library(zoo)



## USGS code:
# cfs: 00060
# Mean: 00003
# Median: 00008 

## Select AK gauges
stations <- read.csv("Input_data/USGS_stations.csv")
stations$station <- as.character(stations$station)

## Download data from USGS
USGS.Q <- readNWISdv(siteNumbers = stations$station,
                      parameterCd = "00060")

# Rename columns from USGS names
USGS.Q <- renameNWISColumns(USGS.Q)
head(USGS.Q)

# Join river names
names(USGS.Q) <- c("agency", "station", "date", "cfs_0", "flow_cd")
Q.usgs <- left_join(USGS.Q, stations) %>% mutate(cfs = ifelse(flow_cd %in% c("A", "A e", "A R", "A [4]"), cfs_0, NA))
Q.usgs$year <- format(Q.usgs$date, "%Y")
head(Q.usgs)

## Export
write.csv(Q.usgs, "Output_data/All_USGS_Q.csv")


### Read back in all USGS cfs data
Q.usgs0 <- read.csv("Output_data/All_USGS_Q.csv")

Q.usgs0$date <- as.Date(Q.usgs0$date)


## Fill in NA rows for all missing daily time steps
# First remove all NA to remove trailing and leading NAs
Q.usgs0 <- Q.usgs0 %>% filter(is.na(cfs) == F)

# Fill in all missing daily timesteps
Q.usgs <- Q.usgs0 %>% filter(is.na(date) == F) %>%
  group_by(name) %>% tidyr::complete(date = seq(min(date),
                               max(date), by = "day"))


## Create dataframe with data ranges and gaps for each river
# Start and end dates by site
Q.range <- Q.usgs %>% group_by(name) %>% summarise(start = min(date), end = max(date))


####### 55 year window
Q.long <- Q.usgs %>% group_by(name) %>% filter(date > "1966-01-01") %>%
  summarise(start = min(date), end = max(date)) %>% 
  filter(start < "1980-01-01" & end > "2015-01-01")
Q.usgs.l <- Q.usgs %>% filter(name %in% unique(Q.long$name),
                              date > "1966-01-01")

## Check percent missing
# Total NA < 30%
miss <- Q.usgs.l %>% group_by(name) %>% summarise(n = length(cfs),
                                                            n.na = sum(is.na(cfs)),
                                                            per.na = n.na/n) %>% filter(per.na > .30)
miss

# By month NA > 30%
Q.usgs.l$month <- month(Q.usgs.l$date)
Q.usgs.l$year <- year(Q.usgs.l$date)
Q.usgs.l.m <- Q.usgs.l %>% group_by(name, year, month) %>% mutate(miss = sum(is.na(cfs))) %>%
  summarise(cfs = ifelse(max(miss)<=10, mean(cfs, na.rm = T), NA))
miss.m <- Q.usgs.l.m %>% group_by(name, month) %>% summarise(n = length(cfs),
                                                                       n.na = sum(is.na(cfs)),
                                                                       per.na = n.na/n) %>% filter(per.na > .30)

miss.m

# Not missing more than 4 years in a row for any month
gaps.m <- Q.usgs.l.m %>% filter(cfs %in% NA) %>% group_by(name, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)
gaps.m

# Identify data gaps greater than 3 months
Q.usgs.l.m$date <- as.Date(paste(Q.usgs.l.m$year, Q.usgs.l.m$month, "15", sep = "-"))
gaps <- Q.usgs.l.m %>% filter(cfs %in% NA) %>% group_by(name) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4)
gaps 

unique(Q.usgs.l.m$name)
Q.range <- Q.usgs.l.m %>% group_by(name) %>% summarise(start = min(date), end = max(date))
Q.range

stats <- Q.usgs.l.m %>% filter(!name %in% c(miss$name, miss.m$name, gaps.m$name, gaps$name)) %>% ungroup() %>% select(name) %>% unique()

# Select stations to include. Stikine not included because of later start date
Q.usgs.l_m <- Q.usgs.l.m %>% filter(date >= "1967-01-15", name %in% c('KUPARUK', 'CHENA_1', 'FISH_SE', 'KENAI_1', 'L_CHENA', 'SALCHA', 'SHIP', 'TALKEETNA', 'TANANA_3', 'YUKON_1', 'YUKON_4', 'L_SUSITNA'),
                                    date > '1974-05-15' | name != 'KUPARUK')
  
Q.usgs.l <- Q.usgs.l %>% filter(date >= "1967-01-01", name %in% c('KUPARUK', 'CHENA_1', 'FISH_SE', 'KENAI_1', 'L_CHENA', 'SALCHA', 'SHIP', 'TALKEETNA', 'TANANA_3', 'YUKON_1', 'YUKON_4', 'L_SUSITNA'),
                                date >= '1974-05-28' | name != 'KUPARUK')

unique(Q.usgs.l_m$name)

write.csv(Q.usgs.l_m, "Data/Q60_monthly.csv", row.names = F)
write.csv(Q.usgs.l, "Data/Q60_daily.csv", row.names = F)




####### 20 year window
Q.20 <- Q.usgs %>% group_by(name) %>%
  summarise(start = min(date), end = max(date)) %>% 
  filter(start < "2005-01-01" & end > "2015-01-01")
Q.usgs.20 <- Q.usgs %>% filter(name %in% unique(Q.20$name),date > "2000-01-01") %>%
  filter(date > '2001-03-31' | !name == 'YUKON_3',
         date > '2001-05-16' | !name == 'SAWMILL_2',
         date > '2001-04-30' | !name == 'WILLOW',
         date > '2003-09-30' | !name == 'TAIYA',
         date > '2000-05-19' | !name == 'MATANUSKA') # Trim data for stations with bad data just at the beginning

Q.range.20 <- Q.usgs.20 %>% group_by(name) %>% summarise(start = min(date), end = max(date))

## Check percent missing
# Total NA < 30%
miss <- Q.usgs.20 %>% group_by(name) %>% summarise(n = length(cfs),
                                                  n.na = sum(is.na(cfs)),
                                                  per.na = n.na/n) %>% filter(per.na > .30)
miss

# By month NA > 10%
Q.usgs.20$month <- month(Q.usgs.20$date)
Q.usgs.20$year <- year(Q.usgs.20$date)
Q.usgs.20.m <- Q.usgs.20 %>% group_by(name, year, month) %>% mutate(miss = sum(is.na(cfs))) %>%
  summarise(cfs = ifelse(max(miss)<=10, mean(cfs, na.rm = T), NA))
miss.m <- Q.usgs.20.m %>% group_by(name, month) %>% summarise(n = length(cfs),
                                                             n.na = sum(is.na(cfs)),
                                                             per.na = n.na/n) %>% filter(per.na > .30)

miss.m

# Not missing more than 4 years in a row for any month
gaps.m <- Q.usgs.20.m %>% filter(cfs %in% NA) %>% group_by(name, month) %>% mutate(group = cumsum(c(TRUE, diff(year) > 4))) %>%
  group_by(name, month, group) %>% summarise(n = length(group), gap_start = min(year), gap_end = max(year)) %>% filter(n >4)
gaps.m 

# Identify data gaps greater than 4 months
Q.usgs.20.m$date <- as.Date(paste(Q.usgs.20.m$year, Q.usgs.20.m$month, "15", sep = "-"))
gaps <- Q.usgs.20.m %>% filter(cfs %in% NA) %>% group_by(name) %>% mutate(group = cumsum(c(TRUE, diff(date) >70))) %>%
  group_by(name, group) %>% summarise(n = length(group), gap_start = min(date), gap_end = max(date)) %>% filter(n >4) %>% filter(gap_start < "2020-01-01")
gaps <- gaps %>% filter(!name %in% c("KUSKOKWIM_1", "LEMON", "SUSITNA"))
gaps



Q.range <- Q.usgs.20 %>% group_by(name) %>% summarise(start = min(date), end = max(date))
Q.range

Q.usgs_20 <- Q.usgs.20 %>% filter(date >= "2000-01-01", !name %in% c(miss.m$name, miss$name, gaps.m$name, gaps$name),
                                      !name %in% c("TANANA_2", "KENAI_2", "KENAI_3", "CHENA_4"))
Q.usgs_20 <- merge(Q.usgs_20, stations)

Q.usgs_20 %>% select(name, name_full, station) %>% unique()

write.csv(Q.usgs_20, "Output_data/Q20_daily.csv", row.names = F)
