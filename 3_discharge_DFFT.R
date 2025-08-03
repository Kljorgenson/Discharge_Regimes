library(ggplot2)
library(dplyr)
library(lubridate)
library(discharge)
library(tidyverse)
library(data.table)
library(zoo)

###### Use the discrete fast Fourier transform (DFFT) to quantify seasonal and stochastic variation of river discharge
#### Recent period
### Combine discharge data from Alaska and Canada
# Alaska river discharge
Qm <- read.csv("Output_data/Q20_daily.csv") %>%
  dplyr::select(date, name, cfs) %>% 
  mutate(discharge = cfs*0.028316846592) # Covert from CFS to m^3/s
head(Qm)

# Canada river discharge
Qm_can <- read.csv("Output_data/Q20_daily_can.csv")%>% 
  dplyr::select(name_full, station, date, discharge) %>% rename(name = name_full)
head(Qm_can) 

# Join AK and CAN data
Q.m <- full_join(Qm, Qm_can) %>% filter(is.na(name) == F)
Q.m <- Q.m %>% mutate(date = as.Date(date,format = "%Y-%m-%d")) 

# Trim date ranges
Q.m <- Q.m %>% group_by(name) %>% filter(date >= "2000-01-01" & date <= "2022-12-31") %>%
  # Extra filtering to remove the start of the record for problem sites
  filter(date >= '2000-05-01' | name != 'MATANUSKA', date >= '2004-04-01' | name != 'MARGUERITE_2', date >= '2005-01-01' | name != 'GRASS RIVER ABOVE STANDING STONE FALLS')
head(Q.m)

# Date ranges
range <- Q.m %>% group_by(name, station) %>% summarise(start = min(date), end = max(date), length = (end-start)/365)
max(range$start)
min(range$end)

# Format data for 'discharge' package: Date (YYY-MM-DD) in column1, discharge in column 2
Q.dat <- Q.m %>% dplyr::select(site = name, date, discharge)
Q.dat <- as.data.frame(Q.dat)


#### Run DFFT ####
## Function for running DFFT
flow.func <- function(Site){
dat <- Q.dat %>% filter(site == Site) %>% select(date, discharge) %>% arrange(date)
dat.flow <- asStreamflow(dat, river.name = Site, max.na = 5000) # Remove cutoff for maximum NA values allowed
dat.seas<-fourierAnalysis(dat.flow)

data <- dat.seas$signal
data$site <- paste(Site)
data$seasonal <- dat.seas$seasonal 
data <- data %>% mutate(across(seasonal, ~ifelse(seasonal == 1, "yes", "no")))
return(data)
}

# Run DFFT for all sites
sites <- unique(Q.dat$site)
df.list <- lapply(sites, flow.func)
flow.dat <- rbindlist(df.list)
head(flow.dat)

# Check date range
# fourierAnalysis trims a year off the end of Alces and Trail, the earliest ending sites
flow.dat %>% filter(is.na(resid.sig)==F) %>% group_by(site) %>% summarise(start = min(date), end = max(date)) %>% View()



#### Calculate discharge regime summary metrics
### Temporal variation in high flows (σHF)
# Loop through all rivers
lis <- list()
for( i in sites){
  dat <- Q.dat %>% filter(site == i) %>% dplyr::select(date, discharge) %>% arrange(date)
  dat.flow <- asStreamflow(dat, river.name = unique(dat$site), max.na = 900) 
  df <- data.frame(n.floods = sigmaHighFlows(dat.flow)$n.floods, # Functions from the 'discharge' package
                   sigma.hfb = sigmaHighFlows(dat.flow)$sigma.hfb)
  
  df$start.date <- dat.flow$start
  df$end.date <- dat.flow$end
  df$site <- i
  lis[[i]] <- df
}

sigmaF <- rbindlist(lis)


### Seasonality DFFT metrics (SNR, Arms)
# Function for saving summary metrics from DFFT
metric.func <- function(Site){
  dat <- Q.dat %>% filter(site == Site) %>% dplyr::select(date, discharge) %>% arrange(date)
  dat.flow <- asStreamflow(dat, river.name = unique(dat$site), max.na = 5000) # Remove cutoff for maximum NA values allowed
  dat.seas<-fourierAnalysis(dat.flow)

  
  data <- dat.seas$signal 
  HSAF <- getHSAF(data$resid.sig, data$year) %>% summarise(HSAF = mean(HSAF, na.omit = T)) # Function from 'discharge' package for
  
  data <- data.frame(site = Site, seasonal = dat.seas$seasonal, rms.signal = ifelse(is.null(dat.seas$rms$rms.signal) == TRUE, "NA", dat.seas$rms$rms.signal),
                     rms.noise = dat.seas$rms$rms.noise, snr = ifelse(is.null(dat.seas$rms$snr) == TRUE, "NA", dat.seas$rms$snr))
  
  return(data)
}


df.list <- lapply(sites, metric.func)
metric.dat <- rbindlist(df.list)

# Join these metrics with σHF from above
metric.dat <- full_join(metric.dat, sigmaF)
head(metric.dat)

# Save summary of metric data
write.csv(metric.dat, "Output_data/DFFT_metrics_20yrs_ak_can.csv", row.names = F)



### Metrics calculated yearly
## High spectral anomaly magnitude (HSAM) and timing of HSAM during open water period (March-October)
# Select seasonal peak for time reference
peak <- flow.dat %>% group_by(site) %>% filter(pred2==max(pred2)) %>% select(site,jday) %>%
  rename(peak = jday)

# Select maximum annual anomaly (HSAM) and its timing
SAM.sum <- full_join(flow.dat,peak) %>% filter(jday >= 60 & jday <= 304) %>% group_by(site, year)  %>%
  filter(length(resid.sig) > 200) %>% # Remove years without full summer period
  filter(resid.sig == max(resid.sig)) %>% 
  mutate(timing = jday-peak) %>% dplyr::select(date,year,jday,timing,resid.sig,peak)

head(SAM.sum)
any(is.na(SAM.sum))

# Export HSAM and HSAM timing data
write.csv(SAM.sum, "Output_data/SAM_20_ak_can.csv", row.names = F)



### Flood duration
# Function
event.func <- function(Site){
  dat <- flow.dat %>% filter(site == Site)
  events <- independentEvents(cutoff.val = 0, data = dat, data.column =10, below.cutoff=F) # Function from 'discharge' package
  events$site <- unique(dat$site)
  return(events)
}

df.list <- lapply(sites, event.func)
event.dat <- rbindlist(df.list)
head(event.dat)

write.csv(event.dat, "Output_data/events.csv", row.names = F)




#### Plot of seasonal signal and discharge anomalies calculated from DFFT

# Make plot of all DFFT data

# Make list of seasonal signal plots
flow.plot <- function(Site){
  dat <- Q.dat %>% filter(site == Site) %>% select(date, discharge) %>% arrange(date)
  dat.flow <- asStreamflow(dat, river.name = Site, max.na = 5000)
  dat.seas<-fourierAnalysis(dat.flow)
  
  png(file=paste("Figures/DFFT_fits/DFFT_seasonal_signal_2000-2022_", Site, ".png", sep = ""), width = 1800, height = 1800, pointsize = 40)
  plot(dat.seas)
  title(paste(Site), line = 0.6)
  dev.off()
}

lapply(sites, flow.plot)






################################################################################################################

###### Historical period
# Alaska discharge data
Qm <- read.csv("Output_data/Q60_daily.csv") %>%
  dplyr::select(date, name, cfs) %>% 
  mutate(discharge = cfs*0.028316846592) # Covert from CFS to m^3/s
head(Qm)

# Canada discharge data
Qm_can <- read.csv("Output_data/Q60_daily_can.csv")%>% 
  dplyr::select(name_full, station, date, discharge) %>% rename(name = name_full)
head(Qm_can) 

# Join AK and CAN data
Q.m <- full_join(Qm, Qm_can)
Q.m <- Q.m %>% mutate(date = as.Date(date,format = "%Y-%m-%d"))

# Select date range
Q.m <- Q.m %>% group_by(name) %>% filter(date >= "1970-01-01" & date <= "1992-12-31")
head(Q.m)

# Range
range <- Q.m %>% group_by(name) %>% summarise(start = min(date), end = max(date))
max(range$start)
min(range$end)

# Format for discharge package
Q.dat <- Q.m %>% dplyr::select(site = name, date, discharge)
Q.dat <- as.data.frame(Q.dat)

# Save list of historical period stations
meta1 <- read.csv("Input_data/river_metadata.csv") %>% rename(station = StationNum, name = NameNom)%>% mutate(station = as.character(station))
meta2 <- read.csv("Input_data/USGS_stations.csv") %>% mutate(station = as.character(station))
stations <- full_join(meta1, meta2) %>% select(station, name) %>% filter(name %in% Q.m$name)%>% 
  select(station, name) %>% unique()
write.csv(stations, 'Output_data/stations_historical.csv')


#### Run DFFT ####
sites <- unique(Q.dat$site)
## Function for running DFFT
flow.func <- function(Site){
  dat <- Q.dat %>% filter(site == Site) %>% dplyr::select(date, discharge)
  dat.flow <- asStreamflow(dat, river.name = unique(dat$site), max.na = 5000) # Remove cutoff for maximum NA values allowed
  dat.seas<-fourierAnalysis(dat.flow)
  
  data <- dat.seas$signal
  data$site <- paste(Site)
  data$seasonal <- dat.seas$seasonal 
  data <- data %>% mutate(across(seasonal, ~ifelse(seasonal == 1, "yes", "no")))
  return(data)
}

df.list <- lapply(sites, flow.func)
flow.dat <- rbindlist(df.list)
head(flow.dat)

# Range
range <- flow.dat %>% group_by(site) %>% summarise(start = min(date), end = max(date))
max(range$start)
min(range$end)


#### Calculate discharge regime summary metrics
### Temporal variation in high flows (σHF)
# Loop through all rivers
lis <- list()
for( i in sites){
  dat <- Q.dat %>% filter(site == i) %>% dplyr::select(date, discharge)
  dat.flow <- asStreamflow(dat, river.name = unique(dat$site), max.na = 5000) 
  df <- data.frame(sigma.hfb = sigmaHighFlows(dat.flow)$sigma.hfb)
  
  df$start.date <- dat.flow$start
  df$end.date <- dat.flow$end
  df$site <- i
  lis[[i]] <- df
}

sigmaF <- rbindlist(lis)


### Seasonality DFFT metrics (SNR, Arms)
# Function for saving summary metrics from DFFT
metric.func <- function(Site){
  dat <- Q.dat %>% filter(site == Site) %>% dplyr::select(date, discharge)
  dat.flow <- asStreamflow(dat, river.name = unique(dat$site), max.na = 2000)
  dat.seas<-fourierAnalysis(dat.flow)
  
  data <- dat.seas$signal 
  HSAF <- getHSAF(data$resid.sig, data$year) %>% summarise(HSAF = mean(HSAF, na.omit = T)) # Function from 'discharge' package for
  
  data <- data.frame(site = Site, seasonal = dat.seas$seasonal, rms.signal = ifelse(is.null(dat.seas$rms$rms.signal) == TRUE, "NA", dat.seas$rms$rms.signal),
                     rms.noise = dat.seas$rms$rms.noise, snr = ifelse(is.null(dat.seas$rms$snr) == TRUE, "NA", dat.seas$rms$snr))
  
  return(data)
}

df.list <- lapply(sites, metric.func)
metric.dat <- rbindlist(df.list)

# Join these metrics with σHF from above
metric.dat <- full_join(metric.dat, sigmaF)
head(metric.dat)

# Save summary of metric data
write.csv(metric.dat, "Output_data/DFFT_metrics_1970-1992_ak_can.csv", row.names = F)




### Metrics calculated yearly
## High spectral anomaly magnitude (HSAM) and timing of HSAM during open water period (March-October)
# Use peaks from recent period for reference for HSAM timing
SAM.sum <- left_join(flow.dat,peak) %>% filter(jday >= 60 & jday <= 304) %>% group_by(site, year)  %>% 
  filter(length(resid.sig) > 200) %>% # Remove years without full summer period
  filter(resid.sig == max(resid.sig, na.rm = T)) %>% 
  mutate(timing = jday-peak) %>% dplyr::select(site, date,year,jday,timing,resid.sig,peak)

head(SAM.sum)
any(is.na(SAM.sum))

# Export HSAM and HSAM timing data
write.csv(SAM.sum, "Output_data/SAM_1970-1992_ak_can.csv", row.names = F)

## Event duration
# Function for event duration
event.func <- function(Site){
  dat <- flow.dat %>% filter(site == Site)
  
  events <- independentEvents(cutoff.val = 0, data = dat, data.column =10, below.cutoff=F)
  events$site <- unique(dat$site)
  
  return(events)
}

df.list <- lapply(sites, event.func)
event.dat <- rbindlist(df.list)
head(event.dat)

write.csv(event.dat, "Output_data/events_1970-1992.csv", row.names = F)
