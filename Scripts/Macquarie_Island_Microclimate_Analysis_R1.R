# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     Identifying plausible microrefugia across a threatened sub-Antarctic
#       ecosystem (Macquarie Island's fell-fields)
#
#       By David Baker (April 2020)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear environment
rm(list = ls())

# Set to UTC
Sys.setenv(TZ='UTC')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     SECTION 1) Set up for analysis                                          
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Set file path
function_path <- "Scripts/Macquarie_Island_Microclimate_functions_R1.R"

#~# Source for functions
source(function_path)

#~# Install microclima
#install.packages("remotes")
#remotes::install_github("ilyamaclean/microclima")
#remotes::install_github("mrke/NicheMapR")

#~# Load packages
packageList <- c(
  # Data handling
  "tidyverse", "readxl", "gtools", "plyr",
  # Spatial data package
  "rgdal", "raster", "rgeos", "sf",
  # Parsing date and time
  "lubridate",
  # plotting
  "gtable", "cowplot", "see", "ggthemes","grid","gridExtra",
  # Microclimate modelling
  "microclima", "NicheMapR",#"elevatr","rts","RCurl",
  # Statistical modelling
  "brms", "DHARMa", "gstat", "Metrics"
)
packageList[which(loadLibrary(packageList)==FALSE)]

#~# Macquarie island projection
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

#~# Elevation layer 
# Elevation: output from the Digital Elevation Model (DEM) at 5 x 5m 
dem <- raster("Data/mi_dem.grd")
# Currently issues a warning because of changes between PROJ 4 and PROJ 6 
# (i.e. dropping support for +towgs). 
list_coordOps(paste0(proj4string(dem), " +type=crs"), "EPSG:32757")
# i.e. Best instantiable operation has only ballpark accuracy
# The correct projection is given is wkt format above and we'll update the 
# projection accordingly.
crs(dem) <- macca_wkt
list_coordOps(paste0(proj4string(dem), " +type=crs"), "EPSG:32757")
# i.e. Best instantiable operation has accuracy: 0 m

# Elevation layer aggregate to 30 x 30 m
dem_l <- aggregate(dem, 20, mean)
dem_l <- trim(dem_l)

# Convert to shapefile of island boundary
dem_l2 <- dem_l
dem_l2[!is.na(dem_l2)] <- 0
macca_boundary <- rasterToPolygons(dem_l2, dissolve = TRUE) 
macca_boundary <- gSimplify(macca_boundary, 200, topologyPreserve=TRUE)

#~# Elevation
elev <- na.omit(data.frame(coordinates(dem_l), elev=getValues(dem_l)))

#~# Site coordinates
# ibuttons to use (sorted by Cath)
ibutton_select_correct <- c("BS-C_59A_H", "BS-E_60A_H", "BS1B_3A_H",
                            "BS2A_17B_H", "BS2B_24A_H", "BS3A_53A_H", 
                            "BS3B_39A_H", "BS4B_6A_H",  "BS4C_12A_H",
                            "BS5A_46B_H", "BS5B_48A_H", "BS6A_42A_H", 
                            "BS6B_44A_H", "BS7A_11A_H", "BS7C_14A_H",
                            "BS8A_34A_H", "GS1A_30A_H", "GS2C_27A_H",
                            "GS2D_25A_H", "GS3A_38A_H", "GS3B_45A_H",
                            "GS4B_29A_H", "GS4D_4A_H",  "GS5A_50A_H", 
                            "GS5D_52A_H", "GS6C_58A_H" ,"GS6D_31A_H",
                            "GS7B_8A_H",  "GS7C_1A_H",  "GS8A_26B_H", 
                            "GS8B_20A_H", "PS-D_61A_H", "PS1D_35A_H",
                            "PS2A_23A_H", "PS2X_41A_H", "PS3A_55A_H",
                            "PS4A_9A_H",  "PS4D_10A_H", "PS5A_40A_H",
                            "PS5B_51A_H", "PS7A_21A_H", "PS7B_2A_H",  
                            "PS7D_22A_H", "PS8B_5A_H" , "PS8C_15A_H",
                            "PS8D_16A_H", "RS-B_57A_H", "RS1A_43A_H", 
                            "RS2A_19A_H", "RS2D_18A_H", "RS3A_54A_H",
                            "RS3B_37A_H", "RS3X_49A_H", "RS4A_7A_H" ,
                            "RS5B_47A_H", "RS5D_56A_H", "RS6A_32A_H",
                            "RS7B_62A_H", "RS7C_13A_H", "RS8A_36A_H", 
                            "RS8B_28A_H", "RS8C_33A_H")

#~# ibuttons
ibut <- na.omit(get(load("Data/ibutton_extracted_all.Rdata")))
#' There are multiple sensors at each site. Pick "H" in preference to "T" at
#' sites where both co-occur
ibut <- ibut[ibut$SITE_CODE %in% ibutton_select_correct,]

#~# Coordinates of ibutton sites
ibut_coord <- 
  read.csv("Data/Macquarie_ibutton_and_topographic_Variables.csv")[,c('x','y','SITE_CODE')]
ibut_coord <- ibut_coord[!duplicated(ibut_coord),]
ibut_coord <- ibut_coord[ibut_coord$SITE_CODE %in% ibutton_select_correct,]

#~# Spatial object
coords <- ibut_coord
coordinates(coords) <- ~x+y
proj4string(coords) <- 
  CRS("+proj=utm +zone=57 +south +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Vector of months for the snow free season
mths <- c("Oct","Nov","Dec","Jan","Feb","Mar")
#~#~#~#  Dates
tabDates <- data.frame(mth = mths, 
                       start = c("30/09/2017","31/10/2017","30/11/2017",
                                 "01/01/2017","31/01/2017","28/02/2017"),
                       finish = c("31/10/2017","30/11/2017","31/12/2017",
                                  "31/01/2017","28/02/2017","30/03/2017"),
                       start2 = c("01/10/2017","01/11/2017","01/12/2017",
                                  "01/01/2017","01/02/2017","01/03/2017"))

#~# Things to always keep in the global environment
keep_files <- c("keep_files","shareDir","packList","function_path",
                "macca_boundary","dem","dem_l","elev","ibut","ibut_coord",
                "coords","tabDates","mths")

#~# Tidy workplace and just keep files that are needed later
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 2) Create observation datasets
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Load back in functions
source(function_path)

#~# Monthly summary statistics
ibut_oct <- observed_summary_stats(ibut, "01/10/2017", "31/10/2017", 
                                   coords, "Oct")
ibut_nov <- observed_summary_stats(ibut, "01/11/2017", "30/11/2017", 
                                   coords, "Nov")
ibut_dec <- observed_summary_stats(ibut, "01/12/2017", "31/12/2017", 
                                   coords, "Dec")
ibut_jan <- observed_summary_stats(ibut, "01/01/2017", "31/01/2017", 
                                   coords, "Jan")
ibut_feb <- observed_summary_stats(ibut, "01/02/2017", "28/02/2017", 
                                   coords, "Feb")
ibut_mar <- observed_summary_stats(ibut, "01/03/2017", "30/03/2017", 
                                   coords, "Mar")
ibut_mth_sum <- do.call(rbind, list(ibut_oct,ibut_nov,ibut_dec,
                                    ibut_jan,ibut_feb,ibut_mar))
save(ibut_mth_sum, file = "Outputs/ibutton_monthly_summary_statistics.Rdata")

#~# Tidy workplace and just keep files that are needed later
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 3) Run models (AWS, CRA, LR)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Here, we predict microclimate using models driven by 
#' 1) automated weather station (AWS)
#' 2) climate reanalysis (CRA) data, and a 
#' 3) lapse rate model across
#' 
#' Macquarie Island and evaluate the 
#' 1) predictions using in situ (iButton) observation data.
#' 2) Predictions are made hourly at 4 cm above the ground - iButtons record 4
#'    hourly.
#'  
#' Assessments are made for monthly and daily time periods.

#~# Load back in functions
source(function_path)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~# Run microclimate model for Months using AWS or CRA
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#~# Run at 4cm with AWS+CRA data
oct_hrly_04cm_aws <- run_micro_model('30/09/2017', '31/10/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(oct_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/oct_hrly_04cm_aws.Rdata", 
     compress = 'xz')
nov_hrly_04cm_aws <- run_micro_model('31/10/2017', '30/11/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(nov_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/nov_hrly_04cm_aws.Rdata", 
     compress = 'xz')
dec_hrly_04cm_aws <- run_micro_model('30/11/2017', '31/12/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(dec_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/dec_hrly_04cm_aws.Rdata", 
     compress = 'xz')
jan_hrly_04cm_aws <- run_micro_model('01/01/2017', '31/01/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(jan_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/jan_hrly_04cm_aws.Rdata", 
     compress = 'xz')
feb_hrly_04cm_aws <- run_micro_model('31/01/2017', '28/02/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(feb_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/feb_hrly_04cm_aws.Rdata", 
     compress = 'xz')
mar_hrly_04cm_aws <- run_micro_model('28/02/2017', '30/03/2017', 0.04, 
                                     useObs = TRUE, dem_l, h_NCEP = NULL)
save(mar_hrly_04cm_aws, 
     file = "Outputs/Monthly_predictions_04cm/mar_hrly_04cm_aws.Rdata", 
     compress = 'xz')

#~# Run at 4cm with CRA data
oct_hrly_04cm_cra <- run_micro_model('30/09/2017', '31/10/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(oct_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/oct_hrly_04cm_cra.Rdata", 
     compress = 'xz')
nov_hrly_04cm_cra <- run_micro_model('31/10/2017', '30/11/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(nov_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/nov_hrly_04cm_cra.Rdata", 
     compress = 'xz')
dec_hrly_04cm_cra <- run_micro_model('30/11/2017', '31/12/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(dec_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/dec_hrly_04cm_cra.Rdata", 
     compress = 'xz')
jan_hrly_04cm_cra <- run_micro_model('01/01/2017', '31/01/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(jan_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/jan_hrly_04cm_cra.Rdata", 
     compress = 'xz')
feb_hrly_04cm_cra <- run_micro_model('31/01/2017', '28/02/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(feb_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/feb_hrly_04cm_cra.Rdata", 
     compress = 'xz')
mar_hrly_04cm_cra <- run_micro_model('28/02/2017', '30/03/2017', 0.04, 
                                     useObs = FALSE, dem_l, h_NCEP = NULL)
save(mar_hrly_04cm_cra, 
     file = "Outputs/Monthly_predictions_04cm/mar_hrly_04cm_cra.Rdata", 
     compress = 'xz')

#~# Run lapse rate model
stationData <- read.csv("Data/BoM_Data_MI/2020-04-02-DS016_AADC_AWS/HM01X_Data_300004_5559781380.txt")
oct_lr <- run_lapse_rate(stationData, '30/09/2017', '31/10/2017', 'Oct')
nov_lr <- run_lapse_rate(stationData, '31/10/2017', '30/11/2017', 'Nov')
dec_lr <- run_lapse_rate(stationData, '30/11/2017', '31/12/2017', 'Dec')
jan_lr <- run_lapse_rate(stationData, '31/12/2016', '31/01/2017', 'Jan')
feb_lr <- run_lapse_rate(stationData, '31/01/2017', '28/02/2017', 'Feb')
mar_lr <- run_lapse_rate(stationData, '28/02/2017', '30/03/2017', 'Mar')
lr_dat <- do.call(rbind, list(oct_lr, nov_lr, dec_lr, jan_lr, feb_lr, mar_lr))
save(lr_dat, 
     file = 'Outputs/lapse_rate_daily_estimates_Oct_Mar.Rdata', 
     compress = 'xz')

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 4) Evaluate models
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Load back in functions
source(function_path)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~# Monthly statistics
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

ibut_mth_sum <- get(load("Outputs/ibutton_monthly_summary_statistics.Rdata"))
lr_mth_all <- get(load('Outputs/lapse_rate_daily_estimates_Oct_Mar.Rdata'))

#~#~#~# Evaluate AWS
mth_eval_aws <- lapply(mths, 
                       function(mth)  
                         model_pred_err("AWS", mth, ibut_mth_sum, dem_l)) 
mth_eval_aws <- do.call(rbind, mth_eval_aws)
mth_eval_aws$mth <- factor(mth_eval_aws$mth, levels = mths)

#~#~#~# Evaluate CRA
mth_eval_cra <- lapply(mths, 
                       function(mth)
                         model_pred_err("CRA", mth, ibut_mth_sum, dem_l))
mth_eval_cra <- do.call(rbind, mth_eval_cra)
mth_eval_cra$mth <- factor(mth_eval_cra$mth, levels = mths)

#~#~#~# Evaluate lapserate
mth_eval_lr <- lapply(mths, 
                      function(mth)
                        lr_pred_err("LR", lr_mth_all, mth, ibut_mth_sum))
mth_eval_lr <- do.call(rbind, mth_eval_lr)
mth_eval_lr$mth <- factor(mth_eval_lr$mth, levels = mths)

#~#~#~# Combine data
cs <- c("SITE_CODE", "x", "y", "n_days", 
        "mean_err", "max_err",  "min_err",  
        "mth", "weights", "type")
mth_eval <- rbind(mth_eval_aws[ ,cs], mth_eval_cra[ ,cs], mth_eval_lr[ ,cs])
mth_eval <- mth_eval %>% distinct() 
save(mth_eval, file = "Outputs/monthly_error_rates_Oct_Mar.Rdata")

# Plot error distributions
mth_eval$type[mth_eval$type %in% "AWS"] <- "AWS+CRA" # Label as AWS+CRA now
t_mean_plot <- density_panel_plot(mth_eval, "mean_err", "mean")
t_max_plot <- density_panel_plot(mth_eval, "max_err", "max")
t_min_plot <- density_panel_plot(mth_eval, "min_err", "min")
t_err_panel <- grid.arrange(t_min_plot, t_mean_plot, t_max_plot, ncol = 3, 
                            bottom = textGrob("Error (predicted - observed)"))
ggsave(filename = "Outputs/Figures/Figure_1_terror_panel.png", 
       plot = t_err_panel, width = 8, height = 8)

# Monthly RMSE
rmse_aws <- lapply(mths, function(mth) { 
  tmp <- mth_eval_aws[mth_eval_aws$mth == mth,]
  rmse_aws <- error_evaluation_mth(tmp) 
})
rmse_aws <- do.call(rbind, rmse_aws)
rmse_cra <- lapply(mths, function(mth) { 
  tmp <- mth_eval_cra[mth_eval_cra$mth == mth,]
  rmse_aws <- error_evaluation_mth(tmp) 
})
rmse_cra <- do.call(rbind, rmse_cra)
rmse_lr <- lapply(mths, function(mth) { 
  tmp <- mth_eval_lr[mth_eval_lr$mth == mth,]
  rmse_lr <- error_evaluation_mth(tmp) 
})
rmse_lr <- do.call(rbind, rmse_lr)
rmse_mth_all <- do.call(rbind,list(rmse_aws,rmse_cra,rmse_lr))

#~#~#~# Tables 2 (formatted)
rmse_mth_tab <- rmse_mth_all %>% 
  pivot_wider(id_cols = c("mth"), names_from = c(type), 
              values_from = c(rmse_min, rmse_mean, rmse_max))
rmse_mth_tab <- 
  rmse_mth_tab[,c('mth', 'rmse_min_AWS', 'rmse_mean_AWS', 'rmse_max_AWS', 
                         'rmse_min_CRA', 'rmse_mean_CRA', 'rmse_max_CRA', 
                         'rmse_min_LR', 'rmse_mean_LR', 'rmse_max_LR')]
write.csv(rmse_mth_tab, file = "Outputs/Table_2_monthly_rmse.csv", 
          row.names = FALSE)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~# Daily statistics
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
          
site_codes <- unique(ibut_coord$SITE_CODE)

#~#~#~# Daily quantiles
lapply(1:nrow(tabDates), function(x) {
  
  res <- site_time_series(cra_dat = NULL, aws_dat = NULL, 
                          obsDat = ibut, 
                          lr_dat = lr_mth_all, 
                          obsCoord = ibut_coord, 
                          dstart = tabDates$start2[x], 
                          dfinish = tabDates$finish[x], 
                          mth =  tabDates$mth[x], 
                          grain = 100,
                          dem_rast = dem_l,
                          outPathFile = "Outputs/",
                          outPathFig = "Outputs/Figures/TimeSeries/")
    
} )
ts_all <- lapply(tabDates$mth, function(mth) { 
  tmp <- get(load(paste0("Outputs/Daily_statistics_",mth,'.Rdata')))
})
ts_all <- do.call(rbind, ts_all)
save(ts_all, file = "Outputs/Timeseries_temperatures_daily_all.Rdata")

#~# RMSE daily mean weighted
# AWS
rmse_tmean_aws <- lapply(mths, 
                         function(mth) 
                           model_evaluation_day(ts_all, "AWS", "tmean", mth) )
rmse_tmean_aws <- do.call(rbind, rmse_tmean_aws)
rmse_tmin_aws <- lapply(mths, 
                        function(mth) 
                          model_evaluation_day(ts_all, "AWS", "tmin", mth) )
rmse_tmin_aws <- do.call(rbind, rmse_tmin_aws)
rmse_tmax_aws <- lapply(mths, 
                        function(mth) 
                          model_evaluation_day(ts_all, "AWS", "tmax", mth) )
rmse_tmax_aws <- do.call(rbind, rmse_tmax_aws)
rmse_aws <- do.call(rbind,list(rmse_tmean_aws,rmse_tmin_aws,rmse_tmax_aws))
# CRA
rmse_tmean_cra <- lapply(mths, 
                         function(mth) 
                           model_evaluation_day(ts_all, "CRA", "tmean", mth) )
rmse_tmean_cra <- do.call(rbind, rmse_tmean_cra)
rmse_tmin_cra <- lapply(mths, 
                        function(mth) 
                          model_evaluation_day(ts_all, "CRA", "tmin", mth) )
rmse_tmin_cra <- do.call(rbind, rmse_tmin_cra)
rmse_tmax_cra <- lapply(mths, 
                        function(mth) 
                          model_evaluation_day(ts_all, "CRA", "tmax", mth) )
rmse_tmax_cra <- do.call(rbind, rmse_tmax_cra)
rmse_cra <- do.call(rbind,list(rmse_tmean_cra,rmse_tmin_cra,rmse_tmax_cra))
# LR
rmse_tmean_lr <- lapply(mths, 
                        function(mth) 
                          model_evaluation_day(ts_all, "LR", "tmean", mth) )
rmse_tmean_lr <- do.call(rbind, rmse_tmean_lr)
rmse_tmin_lr <- lapply(mths, 
                       function(mth) 
                         model_evaluation_day(ts_all, "LR", "tmin", mth) )
rmse_tmin_lr <- do.call(rbind, rmse_tmin_lr)
rmse_tmax_lr <- lapply(mths, 
                       function(mth) 
                         model_evaluation_day(ts_all, "LR", "tmax", mth) )
rmse_tmax_lr <- do.call(rbind, rmse_tmax_lr)
rmse_lr <- do.call(rbind,list(rmse_tmean_lr,rmse_tmin_lr,rmse_tmax_lr))
# Combine
rmse_all <- do.call(rbind,list(rmse_aws, rmse_cra, rmse_lr))
write.csv(rmse_all, "Outputs/rmse_r_daily.csv")

#~# RMSE daily statistics plot
rmse_all$lab[rmse_all$quant == "tmean"] <- "'T'[mean]"
rmse_all$lab[rmse_all$quant == "tmin"] <- "'T'[min]"
rmse_all$lab[rmse_all$quant == "tmax"] <- "'T'[max]"
rmse_all$lab <- factor(rmse_all$lab, 
                       levels = c("'T'[min]","'T'[mean]","'T'[max]"))
rmse_all$mth <- factor(rmse_all$mth, levels = mths)
rmse_all$type[rmse_all$type %in% "AWS"] <- "AWS+CRA" # Label as AWS+CRA now

#~# RMSE plots
rmse_daily_p <- 
  ggplot(rmse_all, aes(x = mth, y = rmse_m, group = type)) + 
  geom_point(aes(colour = type, shape = type), 
             position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(x=mth, ymax = (rmse_m + rmse_sd), ymin = (rmse_m - rmse_sd), 
                    colour = type, width = 0.6), 
                position = position_dodge(width=0.75)) + 
  scale_color_colorblind() +
  scale_y_continuous(limits = c(0,NA)) +
  theme_bw(base_size = 14)  %+replace% 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  facet_wrap(~lab, nrow = 1, scale = "free", labeller = "label_parsed") +
  xlab("Month") + ylab(bquote("Across site mean\n RMSE (\u00B1 1 SD)"))

#~# Correlation coefficient plot
rmse_all$u_sd <- rmse_all$r_m + rmse_all$r_sd
rmse_all$u_sd[rmse_all$u_sd > 1] <- 1
r_daily_p <- ggplot(rmse_all, aes(x = mth, y = r_m, group = type)) + 
  geom_point(aes(colour = type, shape = type), 
             position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(x=mth, ymax = u_sd, ymin = (r_m - r_sd), 
                    colour = type, width = 0.6), 
                position = position_dodge(width = 0.75)) + 
  scale_color_colorblind() +
  scale_y_continuous(limits = c(0,1)) + 
  theme_bw(base_size = 14)  %+replace% 
  theme(legend.position = c(0.56,0.25),  legend.title = element_blank()) +
  facet_wrap(~lab, nrow = 1, scale = "free", labeller = "label_parsed") +
  xlab("Month") + ylab(bquote("Across site mean\n r (\u00B1 1 SD)"))

#~# Panel
gA <- ggplotGrob(rmse_daily_p)
gB <- ggplotGrob(r_daily_p)
grid::grid.newpage()
daily_err_pan <- grid.arrange(rbind(gA, gB))
ggsave(filename = "Outputs/Figures/Figure_2_rmse_r_daily_panel.png", 
       plot = daily_err_pan, width = 10, height = 6)

#~# Daily quantiles plots
daily_quantile_plots(ts_all, 'tmean', 'AWS')
daily_quantile_plots(ts_all, 'tmax', 'AWS')
daily_quantile_plots(ts_all, 'tmin', 'AWS')
daily_quantile_plots(ts_all, 'tmean', 'CRA')
daily_quantile_plots(ts_all, 'tmax', 'CRA')
daily_quantile_plots(ts_all, 'tmin', 'CRA')
daily_quantile_plots(ts_all, 'tmean', 'LR')
daily_quantile_plots(ts_all, 'tmax', 'LR')
daily_quantile_plots(ts_all, 'tmin', 'LR')

#~# Tables
rmse_tab <- read.csv("Outputs/rmse_r_daily.csv") %>% 
  mutate(value_m = round(rmse_m,2), value_sd = round(rmse_sd,2)) %>%
  mutate(value = paste0(value_m, " (", value_sd, ")")) %>%
  pivot_wider(id_cols = "mth", 
              names_from = c(type, quant),  
              values_from = value)
rmse_tab$Metric <- "RMSE"
rmse_tab <- rmse_tab[,c('Metric','mth', 'AWS_tmin', 'AWS_tmean', 'AWS_tmax', 
                        'CRA_tmin', 'CRA_tmean', 'CRA_tmax', 'LR_tmin', 
                        'LR_tmean', 'LR_tmax' )]
r_tab <- read.csv("Outputs/rmse_r_daily.csv") %>% 
  mutate(value_m = round(r_m,2), value_sd = round(r_sd,2)) %>%
  mutate(value = paste0(value_m, " (", value_sd, ")")) %>%
  pivot_wider(id_cols = "mth", 
              names_from = c(type, quant), 
              values_from = value)
r_tab$Metric <- "r"
r_tab <- r_tab[,c('Metric','mth', 'AWS_tmin', 'AWS_tmean','AWS_tmax','CRA_tmin', 
                  'CRA_tmean', 'CRA_tmax', 'LR_tmin', 'LR_tmean', 'LR_tmax' )]
rmse_r_tab <- rbind(rmse_tab, r_tab)
write.csv(rmse_r_tab, "Outputs/Table_S1_daily_rmse.csv", row.names = FALSE)

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 5) Spatial prediction errors using iButton data
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Load back in functions
source(function_path)

#~# Time series data
ts_all <- get(load("Outputs/Timeseries_temperatures_daily_all.Rdata"))

#~# Errors over entire period for three approaches
ts_04cm_p <- ts_all %>% 
  right_join(ibut_coord, ., by = "SITE_CODE") %>%
  dplyr::select(x, y, Day, type, err_mean, err_min, err_max, SITE_CODE, mth) %>%
  pivot_longer(cols = c( err_mean, err_min, err_max)) 
ts_04cm_p$id[grep(ts_04cm_p$name, pattern = "_max")] <- "'T'[max]"
ts_04cm_p$id[grep(ts_04cm_p$name, pattern = "_mean")] <- "'T'[mean]"
ts_04cm_p$id[grep(ts_04cm_p$name, pattern = "_min")] <- "'T'[min]"
ts_04cm_p_sum <- ts_04cm_p %>% 
  group_by(SITE_CODE, id, type, x, y) %>%
  dplyr::summarise(m_err = mean(value)) %>%
  drop_na() 
ts_04cm_p_sum$err_sig <- factor(sign(ts_04cm_p_sum$m_err), levels = c(-1,1))
ts_04cm_p_sum$err <- abs(ts_04cm_p_sum$m_err)
ts_04cm_p_sum$type[ts_04cm_p_sum$type %in% "AWS"] <- "AWS+CRA" # Label as AWS+CRA now

#~# Plot
mb_sf <- st_as_sf(macca_boundary)
mc_err_plot <- ggplot() + 
  geom_sf(data = mb_sf, fill='Grey85') +
  geom_point(aes(x=x, y=y, size=err, colour=err_sig), 
             data = ts_04cm_p_sum, shape=19,alpha=0.7)  + 
  scale_colour_manual(values=c('#7570b3','#d95f02'), breaks=c('-1','1'), 
                      labels=c('-','+'), drop=FALSE, name='Err. sign') + 
  scale_size_area(breaks=seq(0,3,1), max_size = 3, limits=c(0,3), 
                  name='|Err.| (\u00B0C)') + 
  theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank(),
                                                 axis.text.x = element_text(angle = 90)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  facet_grid(id ~ type, labeller = 'label_parsed') 
ggsave(filename = 'Outputs/Figures/Figure_3_spatial_error_panel.png', 
       plot = mc_err_plot, width = 6, height = 10)

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 6) Daily statistics models
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Load back in functions
source(function_path)

#~# Time series data
ts_all <- get(load("Outputs/Timeseries_temperatures_daily_all.Rdata"))

#~# Daily observation drivers of variations

# Build model data frame with environmental variables
ts_enviro <- merge(ibut_coord, ts_all, by = "SITE_CODE", all.y = TRUE)

# Calculate aspect from 5m dem, aggregate to 100m and adjust to centre north 
# on zero
aspect_5m <- terrain(dem, opt = "aspect", unit = "degrees", neighbors = 8)
aspect_100m <- aggregate(aspect_5m, 20, mean)
aspect_100m[aspect_100m > 180] <- aspect_100m[aspect_100m > 180] * -1
ts_enviro$Aspect <- raster::extract(aspect_100m, ts_enviro[, c("x", "y")])
ts_enviro$Aspect <- (ts_enviro$Aspect - mean(ts_enviro$Aspect)) /  sd(ts_enviro$Aspect)

# Calculate slope from 5m dem
slope_5m <- terrain(dem, opt = "slope", unit = "degrees", neighbors = 8)
slope_100m <- aggregate(slope_5m, 20, mean)
ts_enviro$Slope <- raster::extract(slope_100m, ts_enviro[, c("x", "y")])
ts_enviro$Slope <- (ts_enviro$Slope - mean(ts_enviro$Slope)) /  sd(ts_enviro$Slope)

# Distance from AWS (convert to km)
aws_posit <- data.frame(lon = 495913.49, lat = 3960910.70)
distFromAWS <- distanceFromPoints(dem, aws_posit)
distFromAWS <- distFromAWS / 1000
ts_enviro$DistToAWS <- raster::extract(distFromAWS, ts_enviro[, c("x", "y")])
ts_enviro$DistToAWS <- (ts_enviro$DistToAWS - mean(ts_enviro$DistToAWS)) /  sd(ts_enviro$DistToAWS)

# Elevation
ts_enviro$Elevation <-  raster::extract(dem_l, ts_enviro[, c("x", "y")])
ts_enviro$Elevation <- (ts_enviro$Elevation - mean(ts_enviro$Elevation)) /  sd(ts_enviro$Elevation)

# Rename type to model
names(ts_enviro)[which(names(ts_enviro) == 'type')] <- "Mod"

#~# Tmean error
# Gaussian
tmean_gaus <- brm(err_mean ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                    DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE), 
               data = ts_enviro, family = gaussian(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
               chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmean_gaus <- add_criterion(tmean_gaus, "loo")
save(tmean_gaus, file = "Outputs/brms_tmean_hs_model_gaussian.Rdata")
# Student
tmean_stud <- brm(err_mean ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                    DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE),
               data = ts_enviro, family = student(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
               chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmean_stud <- add_criterion(tmean_stud, "loo")
save(tmean_stud, file = "Outputs/brms_tmean_hs_model_student.Rdata")
# Evaluate
loo_compare(tmean_gaus,tmean_stud)
pp_check(tmean_stud)
bayes_R2(tmean_stud)

#~# Tmean error
# Gaussian
tmin_gaus <- brm(err_min ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                   DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE), 
              data = ts_enviro,  family = gaussian(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
              chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmin_gaus <- add_criterion(tmin_gaus, "loo")
save(tmin_gaus, file = "Outputs/brms_tmin_hs_model_gaussian.Rdata")
# Student
tmin_stud <- brm(err_min ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                   DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE),
                 data = ts_enviro,  family = student(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
                 chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmin_stud <- add_criterion(tmin_stud, "loo")
save(tmin_stud, file = "Outputs/brms_tmin_hs_model_student.Rdata")
# Evaluate
loo_compare(tmin_gaus,tmin_stud)
pp_check(tmin_stud)
bayes_R2(tmin_stud)

#~# Tmean error
# Gaussian
tmax_gaus <- brm(err_max ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                   DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE),
              data = ts_enviro, family = gaussian(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
              chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmax_gaus <- add_criterion(tmax_gaus, "loo")
save(tmax_gaus, file = "Outputs/brms_tmax_hs_model_gaussian.Rdata")
# Student
tmax_stud <- brm(err_max ~ DistToAWS + Elevation + Aspect + Slope + Mod + 
                   DistToAWS:Mod + Elevation:Mod + Aspect:Mod + Slope:Mod + (1 | SITE_CODE),
                 data = ts_enviro, family = student(), prior = set_prior(horseshoe(df = 3, par_ratio = 0.1)),
                 chains = 4, cores = 4, control = list(adapt_delta = .99, max_treedepth = 12) )
tmax_stud <- add_criterion(tmax_stud, "loo")
save(tmax_stud, file = "Outputs/brms_tmax_hs_model_student.Rdata")
# Evaluate
loo_compare(tmax_gaus,tmax_stud)
pp_check(tmax_stud)
bayes_R2(tmax_stud)

#~#~# Extract regression coefficients  
load("Outputs/brms_tmin_hs_model_student.Rdata")
load("Outputs/brms_tmean_hs_model_student.Rdata")
load("Outputs/brms_tmax_hs_model_student.Rdata")
tmean_brms <- brms_est(tmean_stud, "'T'[mean]")
tmin_brms <- brms_est(tmin_stud, "'T'[min]")
tmax_brms <- brms_est(tmax_stud, "'T'[max]")
coeff_all <- do.call(rbind, list(tmean_brms, tmin_brms, tmax_brms))

#~#~# Format into table
coeff_table <- coeff_all %>% 
  dplyr::select(Var, Model, fname) %>%
  pivot_wider(id_cols = Var, names_from = Model, values_from = fname)
write.csv(coeff_table, "Outputs/Table_S2_Model_coefficients_table.csv", row.names = FALSE)

#~#~# Create plot
coeff_all$Model <- factor(coeff_all$Model, 
                          levels = c("'T'[min]", "'T'[mean]", "'T'[max]"))
coeff_plot <-  ggplot(coeff_all) + 
  geom_point2(aes(x = Var, y = Estimate)) +
  geom_errorbar(aes(x = Var, ymin = l95, ymax = u95), width = 0.3) +
  geom_hline(yintercept = 0) +
  theme_bw(base_size = 14) %+replace% theme(axis.title.y = element_blank()) +
    ylab("Standardised regression coefficients") +
  coord_flip() + facet_wrap(~ Model,nrow = 1, labeller = label_parsed)
ggsave(filename = "Outputs/Figures/Figure_4_Regression_coefficients.png",  
       width = 9, height = 5)

#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 7) Scale evaluation
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Test the effect of scale on model predictions over limited domain

#~# Load back in functions
source(function_path)

#~# Pre-download and process hourly NCEP data
h_NCEP_oct <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("30/09/2017", "31/10/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_oct, file = "Outputs/Scale/h_NCEP_oct.Rdata")
h_NCEP_nov <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("31/10/2017", "30/11/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_nov, file = "Outputs/Scale/h_NCEP_nov.Rdata")
h_NCEP_dec <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("30/11/2017", "31/12/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_dec, file = "Outputs/Scale/h_NCEP_dec.Rdata")
h_NCEP_jan <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("01/01/2017", "31/01/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_jan, file = "Outputs/Scale/h_NCEP_jan.Rdata")
h_NCEP_feb <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("31/01/2017", "28/02/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_feb, file = "Outputs/Scale/h_NCEP_feb.Rdata")
h_NCEP_mar <- hourlyNCEP(lat = -54.3, long = 158.57, 
                         tme = as.POSIXct(c("28/02/2017", "30/03/2017"), 
                                          format = "%d/%m/%Y", tz = 'UTC'))
save(h_NCEP_mar, file = "Outputs/Scale/h_NCEP_mar.Rdata")
f <- list.files("Outputs/Scale/", pattern = 'h_NCEP', full.names = TRUE)
for(i in 1:length(f)) load(f[[i]])

#~# Run models over limited domain at several scales and evaluate quantiles
for(i in 2:nrow(tabDates)) {
  
  scale_comp_run(as.character(tabDates$start[i]), 
                 as.character(tabDates$finish[i]), 
                 as.character(tabDates$start2[i]),
                 grains = c(100, 50),  
                 mth = as.character(tabDates$mth[i]),  
                 dem_in = dem, 
                 h_NCEP = get(paste0("h_NCEP_", tolower(tabDates$mth[i]))),
                 Xmin = 483000, Xmax = 500000, Ymin = 3945000, Ymax = 4000000, "N")
}

#~# 100 m vs 50 m 
f <- list.files("Outputs/Scale/", pattern = 'summary', full.names = TRUE)
scale_sum <- lapply(f, function(x) {
  tmp <- get(load(x))
})
scale_sum <- do.call(rbind, scale_sum)
scale_sum$grain <- factor(scale_sum$grain, levels = c("100", "50", "25"))
scale_sum_l <- pivot_longer(scale_sum, cols = c("tmean", "tmin", "tmax"), 
                            names_to = 'Quantile')
scale_sum_l$id[grep(scale_sum_l$Quantile, pattern = "max")] <- "'T'[max]"
scale_sum_l$id[grep(scale_sum_l$Quantile, pattern = "mean")] <- "'T'[mean]"
scale_sum_l$id[grep(scale_sum_l$Quantile, pattern = "min")] <- "'T'[min]"
scale_sum_l$type[scale_sum_l$type %in% "AWS"] <- "AWS+CRA" # Label as AWS+CRA now

#~# Create plot
scale_plot <- scale_sum_l %>% 
  pivot_wider(id_cols = c(SITE_CODE, Day, mth, id, type), 
              values_from = value, names_from = grain, names_prefix = 'x') %>%
  ggplot() + geom_point2(aes(x=x100, y=x50, colour = mth), alpha = 0.3) + 
  scale_color_tableau() +
  facet_grid(id ~ type, scales = "fixed", labeller = 'label_parsed') +
  theme_bw(base_size = 14) %+replace% theme(legend.position = 'bottom', 
                                            legend.title = element_blank()) +
  xlab("100 m spatial grain") + ylab("50 m spatial grain")
ggsave(filename = 'Outputs/Figures/scale_error_panel.png', 
       plot = scale_plot, width = 8, height = 8)

#~# Correlation
r_scale <- scale_sum_l %>% 
  pivot_wider(id_cols = c(SITE_CODE, Day, mth, id, type), 
              values_from = value, names_from = grain, names_prefix = 'x') %>%
  mutate(error = x100 - x50) %>%
  group_by(id, type) %>%
  dplyr::summarise(r = cor(x100, x50), error_m = mean(error))
write.csv(r_scale, "Outputs/Scale_test_correlation_coef.csv", row.names = FALSE)


#~# Tidy workplace
rm(list=ls()[!ls() %in% c(keep_files)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   SECTION 8) Ibutton deployment plots
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Summarise deployment dates
firstLastDate <- ibut %>%
  group_by(SITE_CODE) %>%
  dplyr::summarise(firstDate = min(Date_iButton, na.rm = TRUE),
            lastDate = max(DATE_2, na.rm = TRUE)) %>%
  arrange(firstDate) 

# First and last deployment date
range(firstLastDate$firstDate)
range(firstLastDate$lastDate)

# Plot
deploy <- ggplot(firstLastDate) + 
  geom_segment(aes(x = firstDate, xend = lastDate, y = 1:62, yend = 1:62), size = 3, alpha = 0.2) +
  theme_minimal(base_size = 14) %+replace% theme(axis.text.y = element_blank()) +
  xlab("Deployment period") + ylab("iButton site")
ggsave(filename = "Outputs/Figures/Figure_S1_Deployment_Dates.png", plot = deploy,
       width = 6, height = 12)


#~#~#~#~#~#~#~#~#~#~#~#~#~# END OF SCRIPT #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
