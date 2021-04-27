# This code example demonstrates how to run the runauto function 
# in microclima to predict microclimate conditions 'anywhere and 
# at anytime' using a digital elevation model (downloaded from Amazon Web 
# Services). The example is developed for Macquarie Island in the 
# sub-Antarctic, but the model could be run for any location.

#~# Set to UTC
Sys.setenv(TZ='UTC')

#~# Install microclima and NicheMapR from github
# Need to install these two packages from github
#install.packages("remotes")
#remotes::install_github("ilyamaclean/microclima")
#remotes::install_github("mrke/NicheMapR")

#~# Load required packages
library(microclima)
library(NicheMapR)
library(elevatr)
library(tidyverse)

#~# Source a dem for Macquarie Island using the elevatr package and AWS

# Grid coordinates defining Macquarie Island
macca_coords <- expand.grid(x = c(485447.5, 496727.5),
                            y = c(3929828, 3962628))

# Projection in well-know text (WKT) format
macca_wkt <- 
  'PROJCS["WGS 84 / UTM zone 57S",
       GEOGCS["WGS 84",
              DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                             AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
              PRIMEM["Greenwich",0,
                     AUTHORITY["EPSG","8901"]],
              UNIT["degree",0.01745329251994328,
                   AUTHORITY["EPSG","9122"]],
              AUTHORITY["EPSG","4326"]],
       UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
       PROJECTION["Transverse_Mercator"],
       PARAMETER["latitude_of_origin",0],
       PARAMETER["central_meridian",159],
       PARAMETER["scale_factor",0.9996],
       PARAMETER["false_easting",500000],
       PARAMETER["false_northing",10000000],
       AUTHORITY["EPSG","32757"],
       AXIS["Easting",EAST],
       AXIS["Northing",NORTH]]'

# Download raster for Macquarie Island
macca_dem <- get_elev_raster(locations = macca_coords, prj = macca_wkt, 
                             z = 12, src = 'aws', clip = "bbox")

# Set areas below 0 zero elevation as the sea (ideally, needs clipping 
# to island boundary)
macca_dem[macca_dem < 0] <- NA

# Aggregate to approx 100 x 100 m
macca_dem <- aggregate(macca_dem, 9, mean) 
plot(macca_dem)

#~# Run microclimate model using microclima's runauto function 

# (see package documentation for further details, e.g. ?runauto)
# For this example, run models at 2 m and 4 cm to illustrate the 
# difference in predicted temperature caused by the height at which the
# predictions are made. Habitat is set to "Barren or sparsely vegetated".

# Predictions at 2 m above the surface
macca_2m <- runauto(r = macca_dem, 
                    dstart = "01/10/2017", 
                    dfinish = "02/10/2017", 
                    hgt = 2, # 2 m
                    l = NA,
                    x = NA,
                    habitat = 16, # Barren or sparsely vegetated
                    r.is.dem = TRUE,
                    coastal = FALSE)

# Predictions at 4 cm above the surface
macca_4cm <- runauto(r = macca_dem, 
                    dstart = "01/10/2017", 
                    dfinish = "02/10/2017", 
                    hgt = 0.04, # 4 cm
                    l = NA,
                    x = NA,
                    habitat = 16, # Barren or sparsely vegetated
                    r.is.dem = TRUE,
                    coastal = FALSE)

# Extract values for a point and plot as time series
site <- data.frame(x = 490472.0, y = 3939745)
coordinates(site) <- ~ x + y

#~# Extract both time series 
# 2 m = blue and 4 cm = red 
temp_pred <- as.data.frame(matrix(nrow = 24, ncol = 2)) 
for(j in 11:35) {
  
  # 2 m predictions
  hrTmp2m <- if_raster(macca_2m$temp[,,j], macca_dem)
  temp_pred[j-11,1] <- raster::extract(hrTmp2m, site)
  
  # 4 cm predictions
  hrTmp4cm <- if_raster(macca_4cm$temp[,,j], macca_dem)
  temp_pred[j-11,2] <- raster::extract(hrTmp4cm, site)
  
}

#~# Plot for a 24 hr period (UTC+11)
ggplot(temp_pred) + 
  geom_line(aes(x = 0:23, y = V1), colour = 'blue') + 
  geom_line(aes(x = 0:23, y = V2), colour = 'Red') +
  geom_text(aes(x = 0, y = temp_pred[1,2], label = "4 cm"), vjust = 2) +
  geom_text(aes(x = 0, y = temp_pred[1,1], label = "2 m"), vjust = 2) +
  theme_minimal(base_size = 14) + 
  xlab("Time (24 hr clock; UTC + 11)") +
  ylab(bquote("Temperature ("*degree*"C)"))



