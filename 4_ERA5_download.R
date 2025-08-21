## Compile and clip ERA5 climate data
library(ecmwfr)


#### Import ERA5 monthly data from CDS
# Region except NL: c(73, -169, 48, -79)
# NL: c(57, -64, 54, -61)

wf_set_key(key = ) # Enter key obtained from cds.climate.copernicus.eu

years <- c(as.character(1970:2022)) # Select date range

for(i in years){
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels",
    "product_type"   = "reanalysis",
    "variable"       = c('2m_temperature'), # Select climate variable
    "year"           = as.character(i),
    "month"          = c('01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12'),
    "day"            = c('01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30', '31'),
    "time"           = c('00:00', '01:00', '02:00','03:00', '04:00', '05:00',  '06:00', '07:00', '08:00',  '09:00', '10:00', '11:00',  '12:00', '13:00', '14:00',  '15:00', '16:00', '17:00','18:00', '19:00', '20:00','21:00', '22:00', '23:00'),
    "area"           = c(73, -169, 48, -79), # Select region
    "format"         = "netcdf",
    "target"         = paste("Era5-temp-", i, ".nc", sep = "") # Name file
  )
  
  ncfile <- wf_request(
    request = request,   
    transfer = TRUE,  
    path = "Output_data/",
    verbose = FALSE
  )
}
