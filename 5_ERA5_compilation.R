## Compile and clip ERA5 climate data
# Some parts of code must be adjusted for each different variable
library(ncdf4)
library(raster) # package for raster manipulation
library(sf) # package for geospatial analysis
library(terra)
library(tidyr)
library(dplyr)
library(data.table)

# list files containing ERA5 raster data (.nc)
listfile <- list.files(path="Output_data/", pattern="Era5-temp-*", full.names = TRUE)
listfile 

for(i in 1:length(listfile)){
## Covert to raster stack
d0 <- rast(listfile[i])

## Read in watershed boundaries from shapefile
shape.watershed1 <- st_read("Input_data/Watershed_export/All_watersheds_1.shp")
shape.watershed2 <- st_read("Input_data/Watershed_export/All_watersheds_2.shp")
shape.watershed <- rbind(shape.watershed1, shape.watershed2)

# Fix invalid geometries
shape.watershed <- st_make_valid(shape.watershed)

# Reproject to match raster CRS
shape.watershed <- st_transform(shape.watershed, crs(d0))

# Extract average value for each hour and watershed
temp = raster::extract(d0, shape.watershed, method="simple", fun=mean, 
                       exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)


# Convert SpatialPolygonDataframe into data frame
temp_df<- as.data.frame(temp)

# Make data long
temp_df <- gather(temp_df, dateUTC, value, 2:(ncol(temp_df)), factor_key = T)

# Extract date
temp_df$date0 <- substr(temp_df$dateUTC, 17, 40) # Clip date out of string: start is 15 for ppt and snowfall, 16 for temperature, 17 for snowmelt

# Convert date from scientific notation if needed
temp_df$date0 <- gsub("(.*)e\\.(.*)", "\\1\\e+\\2", temp_df$date0)

# Format date
temp_df$datetime_UTC <- as.POSIXct(as.numeric(temp_df$date0), origin = "1970-01-01", tz = "UTC")
temp_df$date_UTC <- as.Date(temp_df$datetime_UTC)
temp_df$time <- hour(temp_df$datetime_UTC)

# Add station IDs to dataframe
df <- data.frame(ID = unique(temp_df$ID), station = unique(shape_vect$StationNum))
temp.dat <- merge(temp_df, df)

# Export data
write.csv(temp.dat, paste("Output_data/",substr(start = 72, stop = 200, listfile[i]),".csv", sep = ""), row.names = F)
}
