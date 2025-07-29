## Compile and clip ERA5 climate data

library(ecmwfr)
library(ncdf4)
library(raster) # package for raster manipulation
library(sf) # package for geospatial analysis
library(terra)
library(tidyr)
library(dplyr)
library(data.table)
library(stats)
library(stars)
library(foreach)

# list files
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Raw/Region/", pattern="Era5-temp-*", full.names = TRUE)
listfile 

#for(i in 33:length(listfile)){
for(i in 1:32){
## Covert to raster stack
d0 <- rast(listfile[i])

## Read in watershed boundaries from shapefile
meta <- read.csv("Data/Metadata_AK_CAN.csv")
stations60 <- read.csv("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/Discharge_data/Q_60_monthly_AK_CAN.csv")%>% 
  left_join(meta) %>% select(station) %>% unique()

shape.watershed <- st_read("Data/Watershed_export/All_watersheds_final_9.shp") %>% 
                  filter(StationNum %in% stations60$station) # For 1967-1998, only long-term stations

# Fix invalid geometries
shape.watershed <- st_make_valid(shape.watershed)

# Reproject to match raster CRS
shape.watershed <- st_transform(shape.watershed, crs(d0))

# Convert to terra SpatVector
shape_vect <- vect(st_zm(shape.watershed)) # Drop z dimension

Sys.time()
temp = raster::extract(d0, shape.watershed, method="simple", fun=mean, 
                       exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)
Sys.time()

# Convert SpatialPolygonDataframe into data frame
temp_df<- as.data.frame(temp)

# Make data long
temp_df <- gather(temp_df, dateUTC, value, 2:(ncol(temp_df)), factor_key = T) 
temp_df$date0 <- substr(temp_df$dateUTC, 17, 40) # 15 for ppt and snowfall, 17 for snowmelt
temp_df$date0 <- gsub("(.*)e\\.(.*)", "\\1\\e+\\2", temp_df$date0)
temp_df$datetime_UTC <- as.POSIXct(as.numeric(temp_df$date0), origin = "1970-01-01", tz = "UTC")
temp_df$date_UTC <- as.Date(temp_df$datetime_UTC)
temp_df$time <- hour(temp_df$datetime_UTC)


df <- data.frame(ID = unique(temp_df$ID), station = unique(shape_vect$StationNum))
temp.dat <- merge(temp_df, df)

write.csv(temp.dat, paste("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/Processed/Region/",substr(start = 72, stop = 200, listfile[i]),".csv", sep = ""), row.names = F)
}
