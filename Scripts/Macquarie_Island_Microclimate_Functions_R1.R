# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     Identifying plausible microrefugia across a threatened sub-Antarctic ecosystem 
#     (Macquarie Island's fell-fields) By David Baker (April 2020)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Installing / load packages 
#'
#' @param packList A list of required packages.
#'
#' @export
loadLibrary <- function(packList){
  
  packLoaded <- list.files(.libPaths())
  lapply(packList,function(x){
    if(!(x %in% packLoaded)){install.packages(x)}
    require(x,character.only = TRUE)
  })
  
}

#' Extracting seasonal observed temperature statistics  
#' 
#' @param obsTmp In situ iButton temperature time-series data.
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' @param coordI Coordinates of iButtons.
#' @param mth Month label.
#' 
#' @import plyr
#' 
#' @export Dataframe with monthly summary statistics for temperature.
observed_summary_stats <- function(obsTmp, dstart, dfinish, coordI, mth) {
  
  # Remove start/end
  obsTmp_trim <- remove_start_end(obsTmp, dstart, dfinish)
  
  # Mean seasonal/monthly stats
  tmp_daily <- ddply(obsTmp_trim, .(SITE_CODE, Day), 
                     summarise, 
                     tmean = mean(Temp), 
                     tmin = min(Temp), 
                     tmax = max(Temp),
                     diurnRange = max(Temp) - min(Temp))
  tmp_sum <- ddply(tmp_daily, .(SITE_CODE), 
                   summarise, 
                   n_days = length(tmean),
                   mtmean_obs = mean(tmean, na.rm = TRUE),
                   mtmin_obs = mean(tmin, na.rm = TRUE),
                   mtmax_obs = mean(tmax, na.rm = TRUE),
                   diurnRange = mean(diurnRange, na.rm = TRUE))
  tmp_sum <- na.omit(merge(coordI, tmp_sum, by = 'SITE_CODE', all.y = TRUE))
  tmp_sum$mth <- mth
  tmp_sum

}

#' Remove start and end days that fall with month 
#' 
#' @param obsTmp In situ iButton temperature time-series data.
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' 
#' @export Trimmed timeseries data (where necessary).
remove_start_end <- function(dat, dstart, dfinish) {
  
  dstart <- as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC")
  dfinish <- as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC")
  
  #~# Extract observed temps for given dates and site
  dat <- dat[dat$Date >= as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC") & 
               dat$Date <= as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC"),]
  dat$Day <- as.POSIXct(format(dat$Date, format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC") 
  
  #~# Loop through sites 
  datAll <- lapply(unique(dat$SITE_CODE), function(x) {
    sDat <- dat[dat$SITE_CODE == x,]
    fdd <- as.POSIXct(format(min(sDat$Date, na.rm = TRUE), 
                             format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC")
    ldd <- as.POSIXct(format(max(sDat$Date, na.rm = TRUE), 
                             format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC") 
    if(fdd > dstart & fdd < dfinish) {
      sDat <- sDat[!(1:nrow(sDat) %in% which(sDat$Day <= (fdd + (24*3600)))),]
    }
    if(ldd > dstart & ldd < dfinish) {
      sDat <- sDat[!(1:nrow(sDat) %in% which(sDat$Day >= (ldd - (24*3600)))),]
    }
  })
  datAll <- do.call(rbind, datAll)
  
}

#' Extracting seasonal observed temperature statistics  
#' 
#' @param climDat Temperature predictions from microclima (as dat$temps).
#' @param sites A spatial object of the ibutton sites, used to extract predictions from
#' model outputs.
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#' 
#' @export Summary of monthly statistics predicted by microclima.
extract_season_stats <- function(climDat, sites, dem_in = dem_l) {
  
  ind <- c(rep(1,13),sort(rep(seq(2, (dim(climDat)[3]/24)), 24)))
  day_tmean <- as.data.frame(matrix(nrow = length(sites), ncol = (max(ind)-2))) 
  day_tmin <- as.data.frame(matrix(nrow = length(sites), ncol = (max(ind)-2)))  
  day_tmax <- as.data.frame(matrix(nrow = length(sites), ncol = (max(ind)-2))) 
  
  for(i in 2:(max(ind) - 1)) {
    dayClim <- climDat[,,which(ind == i)]
    temp_pred_a <- as.data.frame(matrix(nrow = length(sites), ncol = 24)) 
    for(j in 1:24) {
      dayClim_i <- if_raster(dayClim[,,j], dem_in)
      temp_pred_a[,j] <- raster::extract(dayClim_i, sites)
    }
    day_tmean[,(i-1)] <- apply(temp_pred_a, 1, function(x) mean(x, na.rm = TRUE) )
    day_tmin[,(i-1)] <- apply(temp_pred_a, 1, function(x) min(x, na.rm = TRUE) )
    day_tmax[,(i-1)] <- apply(temp_pred_a, 1, function(x) max(x, na.rm = TRUE) )
  }
  
  #Summarise
  tmean <- rowMeans(day_tmean, na.rm = TRUE)
  tmin <- rowMeans(day_tmin, na.rm = TRUE)
  tmax <- rowMeans(day_tmax, na.rm = TRUE)
  diurnRange <- rowMeans(day_tmax - day_tmin, na.rm = TRUE)
  
  # Return data.frame
  datOut <- cbind.data.frame(SITE_CODE = sites$SITE_CODE, tmean, tmin, tmax, diurnRange) 

}

#' Calculate relative humidity 
#' 
#' Calculation based on the formula given by BOM (accompanying the AWS data)
#' that takes dry bulb air temperature and dew point temperature was inputs.
#' 
#' @param TDry Dry bulb air temperature (degC).
#' @param TDew Dew point temperature (degC).
#' 
#' @return Relative humidity values.
#' @export 
rel_hum <- function(TDry, TDew) { 
  rh <- (100*exp(1.8096+((17.2694*TDew) / (237.3+TDew))))/exp(1.8096+((17.2694*TDry) / (237.3+TDry))) 
}

#' Convert relation humidity to specific humidity
#' 
#' This function is from https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R
#' but altered so the temperature is input in degC rather the in K and then converted to degC.
#' 
#' @param rh Relative humidity (proportion)
#' @param Tc Temperature (degees C)
#' @param press Air pressure (pascals)
#' 
#' @return Specific humidity values.
#' @export
rh2sh <- function(rh, Tc, press) {
  es <- 6.112 * exp((17.67 * Tc) / (Tc + 243.5))
  e <- rh * es
  p_mb <- press / 100
  qair <- (0.622 * e) / (p_mb - (0.378 * e))
} 

#' Prepare input data to model using Automatic Weather Station (AWS) data
#' 
#' Downloads and interpolates NCEP data and then processes AWS values and substitute 
#' into the NCEP dataframe the AWS derived variables for temperature, pressure, winddir, 
#' windspeed, and humidity.
#' 
#' @param stationData The BOM AWS hourly data.
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#'
#' @return Dataframe of hourly NCEP + AWS data.
#' @export 
station_data_prepare <- function(stationData, dstart, dfinish){
  
  #~# Convert start and finish dates
  dstart <- as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC')
  dfinish <- as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')

  # Select required columns
  stationData <- stationData[,c(
    'Year.Month.Day.Hour.Minutes.in.YYYY.2', 'MM.2', 'DD.2', 'HH24.2', 'MI.format.in.Universal.coordinated.time',
    'Mean.sea.level.pressure.in.hPa',
    'Air.Temperature.in.degrees.C',
    'Dew.point.temperature.in.degrees.C',
    'Wind.direction.in.degrees.true',
    'Wind.speed.in.km.h',
    'Speed.of.maximum.windgust.in.last.10.minutes.in..km.h',
    'Precipitation.since.9am.local.time.in.mm')]
  
  # NOTE: the date time is in UTC
  stationData$DateTime <- apply(stationData, 1, function(x) paste(x[c(1:3)], collapse = '/') )
  stationData$DateTime <- paste(stationData$DateTime,stationData$HH24,
                                stationData$MI.format.in.Universal.coordinated.time, sep = ":")
  stationData$DateTime <- as.POSIXct(stationData$DateTime, format = "%Y/%m/%d:%H:%M")
  
  # Select required columns
  stationData <- stationData[,!(names(stationData) %in% 
                                  c('Year.Month.Day.Hour.Minutes.in.YYYY.2', 'MM.2', 'DD.2', 'HH24.2', 'MI.format.in.Universal.coordinated.time'))]
  
  # Rename for consistency with NCEP data
  colnames(stationData) <- c('Pressure','AirTemp','DewPoint','WindDir','WindSpeed','WindGust','Rainfall','Date')
  
  #~# Subset to required period
  stationData <- stationData[stationData$Date >= dstart & 
                               stationData$Date <= dfinish, ]
  
  #~# Pressure from hPa (hectopascals) to Pa (pascals)
  stationData$Pressure <- stationData$Pressure * 100
  
  #~# Windspeed from km/h to m/sec (1 m per sec = 3.6 km per hour)
  stationData$WindSpeed <- stationData$WindSpeed * (1/3.6)
  
  #~# Convert to humidity
  stationData$humidity_rel <- apply(stationData, 1, function(x) {
    rel_hum(as.numeric(x['AirTemp']),as.numeric(x['DewPoint']))
  })
  stationData$humidity_spec <- apply(stationData, 1, function(x) {
    rh2sh(as.numeric(x['humidity_rel']),as.numeric(x['AirTemp']),as.numeric(x['Pressure']))
  })

  #~# Edit to include specific columns and rename
  stationDataHr <- stationData[,c('Date', 'AirTemp', 'Pressure', 'WindDir', 'WindSpeed', 'humidity_spec')]
  names(stationDataHr) <- c('obs_time','temperature','pressure', 'winddir', 'windspeed','humidity')
  stationDataHr
  
}

#' Run microclima 
#' 
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' @param ht The height of predictions.
#' @param useObs Logical (TRUE/FALSE) value indicating whether NCEP only (FALSE) or NCEP+AWS (TRUE) 
#' should be used.
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#' @param h_NCEP Optionally can provide the pre-obtained NCEP data for the time period.
#' 
#' @return microclima output.
#' @export 
run_micro_model <- function(dstart, dfinish, ht, useObs = FALSE, dem_in, h_NCEP = NULL) {
  
  #~# Downloading NCEP 
  if(is.null(h_NCEP)) {
    h_NCEP <- hourlyNCEP(lat = -54.3, long = 158.57, tme = as.POSIXct(c(dstart, dfinish), format = "%d/%m/%Y", tz = 'UTC'))
  }
  
  if(useObs == TRUE) {
    
    #~# Macca station data
    stationData <- read.csv("Data/BoM_Data_MI/2020-04-02-DS016_AADC_AWS/HM01X_Data_300004_5559781380.txt")
    stationDataHr <- station_data_prepare(stationData, dstart, dfinish)
  
    #~# Merge in with NCEP data
    h_NCEPwtStat <- h_NCEP[, !(names(h_NCEP) %in% c('temperature','pressure', 'winddir', 'windspeed','humidity'))]
    h_NCEPwtStat <- merge(h_NCEPwtStat, stationDataHr, by = 'obs_time', all.x = TRUE)
    h_NCEPwtStat <- h_NCEPwtStat[,names(h_NCEP)]
    h_NCEPwtStat$winddir <- as.numeric(h_NCEPwtStat$winddir)

    #~# Correcting for missing values by inserting NCEP values instead
    for(x in 2:6) h_NCEPwtStat[,x][is.na(h_NCEPwtStat[,x])] <- h_NCEP[,x][is.na(h_NCEPwtStat[,x])] 

    #~# Set data 
    h_NCEP <- h_NCEPwtStat
      
  }  
  
  # #~# Run microclimate model - Downloading NCEP @ 4cm ht = 0.04
  hrly_temp <- runauto(r = dem_in, dstart = dstart, dfinish = dfinish,  
                       hgt = ht, l = NA, x = NA,
                       habitat = 16, hourlydata = h_NCEP, 
                       dailyprecip = NA,
                       use.raster = TRUE, r.is.dem = TRUE, coastal = FALSE,
                       continuous = FALSE, steps = 16, plot.progress = FALSE, summarydata = FALSE)
  
  #plot(if_raster(hrly_temp$temps[,,1],dem_l))
  
  return(hrly_temp)
  
}

#' Run lapserate model
#' 
#' Use the lapserate function from microclima to calculate lapse rate adjusted temperature 
#' based on the supplied dem.
#' 
#' @param stationData The AWS hourly data as an input.
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' @param mth Month label.
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#'  
#' @return lapse rate prediction across the landscape.
#' @export 
run_lapse_rate <- function(stationData, dstart, dfinish, mth, dem_in = dem_l) {
  
  #~# Convert start and finish dates
  dstart <- as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC")
  dfinish <- as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC")
  
  #~# Station data
  stationDataHr <- na.omit(station_data_prepare(stationData, dstart, dfinish))

  # Calculate lapse rate
  stationDataHr$lapseRate <- sapply(1: nrow(stationDataHr), function(x) {
    lapseRate <- microclima::lapserate(stationDataHr$temperature[x], stationDataHr$humidity[x], stationDataHr$pressure[x])
    lapseRate
  })

  # Extract lapse rate corrected temperature for each site / time / date
  temp_lr <- lapply(1:nrow(stationDataHr), function(x){
    hr_lr <-  dem_in * stationDataHr$lapseRate[x] + stationDataHr$temperature[x]
    lr_ibut <- cbind(ibut_coord, Temp = raster::extract(hr_lr, ibut_coord[,1:2]))
    lr_ibut <- lr_ibut[,c('SITE_CODE', 'Temp')]
    lr_ibut$obs_time <- stationDataHr$obs_time[x]
    lr_ibut
  })
  temp_lr <- do.call(rbind,temp_lr)
  
  # Add UTC + 11 time
  temp_lr$Date <- temp_lr$obs_time + (60*60*11)
  temp_lr$DayCont <- as.numeric(as.POSIXct(temp_lr$obs_time , format="%d/%m/%Y", tz = "UTC")) 
  temp_lr$Day <- as.numeric(as.POSIXct(format(temp_lr$Date, format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC")) / (24*3600)
  temp_lr <- temp_lr[temp_lr$Date >= (dstart + 24*60*60) & temp_lr$Date <= dfinish,]
  temp_lr$mth <- mth
  temp_lr

}

#' Weighted Root mean squared error (RMSE)
#' 
#' @param actual Observed values from the ibuttons.
#' @param predicted Prediced values from models.
#' @param weight Number of days of observations.
#' 
#' @return Weighted RMSE for each ibutton site.
#' @export 
weighted.rmse <- function(actual, predicted, weight) {
  
  rmse_w <- sqrt(sum((predicted-actual)^2*weight)/sum(weight))
  
}

#' Calculate monthly error model predictions
#' 
#' Calculate monthly prediction errors between in situ observations and 
#' microclima predictions.
#' 
#' @param type Model type label (AWS or CRA).
#' @param mth Month label.
#' @param obsStats Observed monthly statistics from observed_summary_stats.
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#' 
#' @return Monthly prediction errors
#' @export 
model_pred_err <- function(type, mth, obsStats, dem_in) {
  
  # Climate data
  climDat <- get(load(paste0("Outputs/Monthly_predictions_04cm/",tolower(mth),"_hrly_04cm_",tolower(type),".Rdata")))
  climDat <- climDat$temps
  
  # Calculate seasonal statistics from predictions
  pre <- extract_season_stats(climDat, obsStats, dem_in)
  
  # Obs dataframe
  obs <- na.omit(as.data.frame(obsStats[obsStats$mth == mth, ]))
    
  # Combine and remove na columns
  dat_all <- merge(pre, obs, by = "SITE_CODE")

  # Error rates
  dat_all$mean_err <-  dat_all$tmean - dat_all$mtmean_obs 
  dat_all$max_err <- dat_all$tmax - dat_all$mtmax_obs 
  dat_all$min_err <- dat_all$tmin - dat_all$mtmin_obs
  dat_all$mth <- mth
  dat_all$type <- type
  dat_all$weights <-  dat_all$n_days / sum(dat_all$n_days, na.rm = TRUE)
  dat_all
  
}

#' Calculate monthly error from lapserate model
#' 
#' Calculate monthly prediction errors between in situ observations and 
#' LR model predictions.
#' 
#' @param type Model type label (LR).
#' @param lrDat
#' @param mth Month label.
#' @param obsStats Observed monthly statistics from observed_summary_stats.
#' 
#' @return Monthly prediction errors
#' @export 
lr_pred_err <- function(type, lrDat, mth, obsStats) {
  
  # Obs dataframe
  obs <- na.omit(as.data.frame(obsStats[obsStats$mth == mth, ]))
  
  # Lapse rate data
  lrDat <- lrDat[lrDat$mth == mth,]
  lrs <- ddply(lrDat, .(SITE_CODE, Day), summarise, 
               tmean = mean(Temp), 
               tmin = min(Temp), 
               tmax = max(Temp))
  lrs <- ddply(lrs, .(SITE_CODE), summarise,
               tmean = mean(tmean, na.rm = TRUE),
               tmin = mean(tmin, na.rm = TRUE),
               tmax = mean(tmax, na.rm = TRUE))
  
  # Combine and remove na columns
  dat_all <- merge(lrs, obs, by = "SITE_CODE")
  
  # Error rates
  dat_all$mean_err <- dat_all$tmean - dat_all$mtmean_obs 
  dat_all$max_err <- dat_all$tmax - dat_all$mtmax_obs
  dat_all$min_err <- dat_all$tmin - dat_all$mtmin_obs
  dat_all$mth <- mth
  dat_all$type <- type
  dat_all$weights <-  dat_all$n_days / sum(dat_all$n_days, na.rm = TRUE)
  dat_all
  
}

#' Create density panels
#' 
#' @param dat Monthly temperature prediction errors data.
#' @param q_err Quantile required (Tmin, Tmean, Tmax).
#' @param title Plot title.
#' 
#' @return ggplot object 
#' @export 
density_panel_plot <- function(dat, q_err, title){

  # Clip Tmin because of large errors on two sites in Jan and Dec
  if(q_err == "min_err") dat <- dat[dat[[q_err]] > -.95 & dat[[q_err]] < 3.5, ]# dat[dat[[q_err]] <= -.95 | dat[[q_err]] >= 3.5, ]
  err_p <- ggplot(dat) + 
    geom_density(aes_string(x = q_err, weight = "weights", fill = "type"), trim = TRUE, alpha = 0.5) +
    geom_point2(aes_string(x = q_err, y = -0.2, alpha = "n_days", colour = "type"), shape = 108, size = 4, alpha = 1) +
    geom_vline(xintercept = 0, linetype = 1, size = 1, colour = 'red') +
    geom_vline(xintercept = c(-1,1), linetype = 5) +
    geom_label(aes(x = -Inf, y = Inf, label = mth), fill = 'white', vjust = 1.5, hjust = -0.2, size = 4) +
    scale_color_colorblind() + 
    scale_fill_colorblind() +
    facet_wrap(~mth, scales = 'fixed', nrow = 6) +
    theme_minimal(base_size = 14) %+replace% theme(
                                     axis.title = element_blank(),
                                     axis.text.y = element_blank(),
                                     strip.background = element_blank(),
                                     strip.text.x = element_blank(),
                                     plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
                                     legend.position = c(0.2,0.05), 
                                     legend.title = element_blank()) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
    ggtitle(bquote(T[.(title)]*' error')) 
  if(q_err != 'max_err') err_p <- err_p + theme(legend.position = "none")
  err_p
  
}

#' Monthly weighted root mean squared error (RMSE)
#' 
#' @param dat Monthly temperature prediction and observation data.
#' 
#' @return dataframe of weighted rmse 
#' @export 
error_evaluation_mth <- function(dat){

  # RMSE 
  rmse_min <- weighted.rmse(dat$mtmin_obs, dat$tmin, dat$n_days)
  rmse_mean <- weighted.rmse(dat$mtmean_obs, dat$tmean, dat$n_days)
  rmse_max <- weighted.rmse(dat$mtmax_obs, dat$tmax, dat$n_days)
  datOut <- data.frame(mth = unique(dat$mth), 
                    type = unique(dat$type),
                    rmse_min = rmse_min,
                    rmse_mean = rmse_mean,
                    rmse_max = rmse_max)

}

#' Daily across site evaluation statistics
#' 
#' Weighted mean across sites of RMSE and Pearson's correlation coefficient 
#' for daily statistic 
#' 
#' @param dat Daily temperature prediction and observation data.
#' @param type Model type label (CRA, AWS, LR).
#' @param quant Tmin, Tmean, or Tmax.
#' @param mth Month label.
#'  
#' @return dataframe of weighted rmse 
#' @export 
model_evaluation_day <- function(dat, type, quant, mth) {
  
  site_rmse <- lapply(unique(dat$SITE_CODE), function(x) {
    
    tmp <- dat[dat$SITE_CODE == x & dat$type == type & dat$mth %in% mth, ]
    rmse_d <- Metrics::rmse(tmp[[quant]], tmp[[paste0(quant,"_obs")]])
    corr_d <- cor(tmp[[quant]], tmp[[paste0(quant,"_obs")]], method = "pearson")
    res <- data.frame(SITE_CODE = x, N = nrow(tmp),
                      mth = mth, rmse = rmse_d, r = corr_d)
    
  })
  site_rmse <- do.call(rbind, site_rmse)

  # Weighted mean
  rmse_mean_w <- weighted.mean(site_rmse$rmse, site_rmse$N)
  rmse_sd_w <- sqrt(weighted.mean( (site_rmse$rmse - rmse_mean_w)^2, site_rmse$N ))

  # Pearson's correlation coefficient
  r_mean_w <- weighted.mean(site_rmse$r, site_rmse$N)
  r_sd_w <- sqrt(weighted.mean( (site_rmse$r - r_mean_w)^2, site_rmse$N ))
  
  # Res
  res <- data.frame(type = type,
                    quant = quant,
                    mth = unique(site_rmse$mth),
                    rmse_m = rmse_mean_w,
                    rmse_sd = rmse_sd_w,
                    r_m = r_mean_w,
                    r_sd = r_sd_w)

}

#' Daily time-series data and plots
#' 
#' Plots the time series for each prediction and the observed
#' (mainly for checking purposes rather than inference)
#' 
#' @param cra_dat CRA predictions from microclima.
#' @param aws_dat AWS predictions from microclima.
#' @param obsDat iButton observation data.
#' @param lr_dat Lapse rate predictions.
#' @param obsCoord Coordinates of observation sites.
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' @param mth Month label.
#' @param grain Horizontal grain (e.g. 100 x 100 m).
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#' @param outPathFile 
#' @param outPathFig 
#'  
#' @return dataframe of weighted rmse 
#' @export 
site_time_series <- function(cra_dat = NULL, aws_dat = NULL, obsDat, lr_dat, 
                             obsCoord = ibut_coord, dstart, dfinish, mth, 
                             grain = 100, dem_in = dem_l,
                             outPathFile = "Outputs/",
                             outPathFig = "Outputs/Figures/TimeSeries/") {
  
  # Climate data
  if(is.null(cra_dat)) {
    cra_dat <- get(load(paste0("Outputs/Monthly_predictions_04cm/",tolower(mth),"_hrly_04cm_cra.Rdata")))
    cra_dat <- cra_dat$temps
    aws_dat <- get(load(paste0("Outputs/Monthly_predictions_04cm/",tolower(mth),"_hrly_04cm_aws.Rdata")))
    aws_dat <- aws_dat$temps
  }
  
  # Sites that intersect extent
  obsCoord_sub <- obsCoord[obsCoord$y >= extent(dem_in)[3,] & obsCoord$y <= extent(dem_in)[4,],]
  
  # Site daily 
  dailySite <- lapply(unique(obsCoord_sub$SITE_CODE), function(siteCode) {
  
    print(siteCode)
    
    #~# Coordinates
    x <- obsCoord[which(obsCoord$SITE_CODE==siteCode),'x']
    y <- obsCoord[which(obsCoord$SITE_CODE==siteCode),'y']
    
    #~# Convert to datetime 
    dstart <- as.POSIXct(dstart,format="%d/%m/%Y",tz='UTC') 
    dfinish <- as.POSIXct(dfinish,format="%d/%m/%Y",tz='UTC') 
    if(mth == "Jan") dstart <- dstart + 24*60*60 
  
    #~# Time day start
    d1 <- as.numeric(as.POSIXct(dstart,format="%d/%m/%Y",tz='UTC'))
    d2 <- as.numeric(as.POSIXct(dfinish,format="%d/%m/%Y",tz='UTC'))
    
    #~# Extract site specific information
    site <- data.frame(siteCode, x, y)
    coordinates(site) <- ~x + y
    nHours <- dim(aws_dat)[3]
    proj4string(site) <- CRS('+proj=utm +zone=57 +south +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
    
    #~# Observed data (already at UTC+11)
    obs <- na.omit(obsDat[obsDat$SITE_CODE == siteCode & 
                            obsDat$Date >= dstart & obsDat$Date <= dfinish,])
    
    if(nrow(obs) > 0 ) {
      
      obs$DayCont <- as.numeric(as.POSIXct(obs$Date, format="%d/%m/%Y", tz = "UTC")) 
      obs$Date <- as.POSIXct(obs$DayCont,tz = "UTC", origin = "1970/01/01")
      obs$Day <- as.numeric(as.POSIXct(format(obs$Date,format="%d/%m/%Y"), format="%d/%m/%Y", tz = "UTC")) / (24*3600)
      obs$DayCont <- obs$DayCont / d1 
      obs <- obs[,c('DayCont','Day','Temp')]
      obs$Mod <- 'Obs'

      #~# AWS i <- 1
      temp_aws <- rep(NA, nHours)
      for(i in 1:nHours){
        climDat_i <- if_raster(aws_dat[,,i], dem_in)
        temp_aws[i] <- raster::extract(climDat_i, site)
      }
      # Set time and convert to UTC+11
      aws <- data.frame(hr = 1:nHours, Temp = temp_aws)
      aws$DayCont <- apply(aws , 1, function(x) d1 + ((x[1]-1)*60*60) ) - 13*60*60 # UTC+11 time
      aws$Date <- as.POSIXct(aws$DayCont,tz = "UTC", origin = "1970/01/01")
      aws$Day <- as.numeric(as.POSIXct(format(aws$Date, format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC")) / (24*3600)
      aws <- aws[aws$Date >= dstart & aws$Date <= dfinish,]
      aws$DayCont <- aws$DayCont / d1
      aws <- aws[,c('DayCont','Day','Temp')]
      aws$Mod <- 'AWS'

      #~# CRA
      temp_cra <- rep(NA,nHours)
      for(i in 1:nHours){
        climDat_i <- if_raster(cra_dat[,,i], dem_in)
        temp_cra[i] <- raster::extract(climDat_i,site)
      }
      # Set time and convert to UTC+11
      cra <- data.frame(hr = 1:nHours, Temp = temp_cra)
      cra$DayCont <- apply(cra , 1, function(x) d1 + ((x[1]-1)*60*60) ) - 13*60*60 # UTC+11 time
      cra$Date <- as.POSIXct(cra$DayCont,tz = "UTC", origin = "1970/01/01")
      cra$Day <- as.numeric(as.POSIXct(format(cra$Date, format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC")) / (24*3600)
      cra <- cra[cra$Date >= dstart & cra$Date <= dfinish,]
      cra$DayCont <- cra$DayCont / d1
      cra <- cra[,c('DayCont','Day','Temp')]
      cra$Mod <- 'CRA'
      
      #~# Lapse rate model (datetime is recorded in UTC but reflects UTC+11; convert to UTC+11)
      lr <- lr_dat[lr_dat$SITE_CODE == siteCode & lr_dat$mth == mth,]
      lr$obs_time <- lr$obs_time + (60*60*11)
      lr$DayCont <- as.numeric(as.POSIXct(lr$obs_time , format="%d/%m/%Y", tz = "UTC")) 
      lr$Date <- as.POSIXct(lr$DayCont,tz = "UTC", origin = "1970/01/01")
      lr$Day <- as.numeric(as.POSIXct(format(lr$Date, format = "%d/%m/%Y"), format = "%d/%m/%Y", tz = "UTC")) / (24*3600)
      lr <- lr[lr$Date >= dstart & lr$Date <= dfinish,]
      lr$DayCont <- lr$DayCont / d1
      lr <- lr[,c('DayCont','Day','Temp')]
      lr$Mod <- 'LR'
     
      ### Plot
      dat_p <- do.call(rbind, list(aws, cra, lr, obs))
      dat_p$Mod <- factor(dat_p$Mod, levels = c('Obs', 'AWS', 'CRA', 'LR'))
      lmin <- min(dat_p$DayCont)
      lmax <- max(dat_p$DayCont)
      trend_p <- ggplot(dat_p) +
        geom_line(aes(x = DayCont, y = Temp, colour = Mod)) + 
        scale_x_continuous(breaks = c(lmin, lmax), labels = c(dstart,dfinish)) +
        scale_colour_colorblind(name = "Model") + 
        theme_bw(base_size = 16) + 
        xlab("Day in month") + ylab(bquote("Temperature ("* degree*C*")")) +
        ggtitle(paste0('Site code: ', siteCode, ' (', mth,', ', grain, 'm)'))
      sc <- sub("-","_", siteCode)
      ggsave(filename = paste0(outPathFig,sc,"_",mth,"_",grain,".png"), plot = trend_p, width = 14, height = 7)
      
      #~#~# Calculate daily mean, min, max evaluation statistics
      #~# Obs vs Pred
      # Summarise prediction data by day
      sum_day_aws <- aws %>%
        group_by(Day) %>%
        summarise(tmean = mean(Temp, na.rm = TRUE),
                  tmin = min(Temp,na.rm = TRUE),
                  tmax = max(Temp, na.rm = TRUE)) 
      sum_day_aws$type <- 'AWS'
      # Summarise prediction data by day
      sum_day_cra <- cra %>%
        group_by(Day) %>%
        summarise(tmean = mean(Temp, na.rm = TRUE),
                  tmin = min(Temp,na.rm = TRUE),
                  tmax = max(Temp, na.rm = TRUE)) 
      sum_day_cra$type <- 'CRA'
      # Summarise observation data by day
      sum_day_lr <- lr %>%
        group_by(Day) %>%
        summarise(tmean = mean(Temp, na.rm = TRUE),
                  tmin = min(Temp,na.rm = TRUE),
                  tmax = max(Temp, na.rm = TRUE)) 
      sum_day_lr$type <- 'LR'
      # Sum prediction
      sum_pred <- do.call(rbind, list(sum_day_aws, sum_day_cra, sum_day_lr))
  
      # Summarise observation data by day
      sum_day_obs <- obs %>%
        group_by(Day) %>%
        summarise(tmean_obs = mean(Temp, na.rm = TRUE),
                  tmin_obs = min(Temp,na.rm = TRUE),
                  tmax_obs = max(Temp, na.rm = TRUE)) 
  
      # Merge by day 
      sum_day <- merge(sum_pred, sum_day_obs, by = 'Day', all.x = TRUE)
      sum_day <- na.omit(sum_day)
  
      # Error
      sum_day$err_mean <- sum_day$tmean - sum_day$tmean_obs 
      sum_day$err_min <- sum_day$tmin - sum_day$tmin_obs 
      sum_day$err_max <- sum_day$tmax - sum_day$tmax_obs 
  
      # Add site code
      sum_day$SITE_CODE <- siteCode
      sum_day$mth <- mth
      sum_day$grain <- grain
      sum_day
      
    }
    
  })
    
  dailySite <- do.call(rbind, dailySite)
  save(dailySite, file = paste0(outPathFile,"Daily_statistics_",mth,"_", grain ,'.Rdata'), compress = "xz")
  return(dailySite)
}

#' Daily quantile scatterplots (for supplementary material)
#' 
#' Site by site daily scatterplots for each quantile.
#' 
#' @param dat
#' @param quant
#' @param type
#'  
#' @return dataframe of weighted rmse 
#' @export 
daily_quantile_plots <- function(dat, quant, type) {
  
  tmp <- dat[dat$type == type, ]
  
  pp <- ggplot(tmp) + 
    geom_point(aes_string(x = paste0(quant,'_obs'), y = quant, colour = 'mth'), alpha = 0.5, size = 0.95) + 
    scale_color_tableau() +
    geom_abline(slope = 1, intercept = 0, colour = 'grey70') +
    facet_wrap(~SITE_CODE, ncol = 5) +   
    theme_lucid(base_size = 14) %+replace% theme(strip.background = element_blank(),strip.text.x = element_blank(),
          legend.position = 'bottom',legend.direction = "horizontal", legend.title = element_blank()) +
    xlab(bquote("Observed "*T[.(quant)])) + ylab(bquote("Predicted "*T[.(quant)])) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  ggsave(filename = paste0("Outputs/Figures/Daily_",quant,"_v_obs_",type,".png"), plot = pp, width = 8, height = 12)

}

#' Run limited domain microclimate model and multiple scales
#' 
#' 
#' @param dstart Start date of monthly time period in dd/mm/yyyy format.
#' @param dfinish End date of monthly time period in dd/mm/yyyy format.
#' @param date_start2 Start date of monthly time period in dd/mm/yyyy format for plotting (i.e
#' should start on first of the month, whereas dstart usually starts the last day of the previous
#' month).
#' @param grain Horizontal grain (e.g. 100 x 100 m).
#' @param mth Month label.
#' @param dem_in The digital elevation model (dem) used to run the microclimate models.
#' @param h_NCEP Optionally can provide the pre-obtained NCEP data for the time period.
#' @param Xmin Minimum longitude of domain.
#' @param Xmax Maximum longitude of domain.
#' @param Ymin Minimum latitude of domain.
#' @param Ymax Maximum latitude of domain.
#' @param location Name of location (e.g. N or S).
#' 
#' @return microclima predictions
#' @export 
scale_comp_run <- function(date_start, date_end, date_start2, grains = c(100, 50), 
                           mth, dem_in = dem, h_NCEP = NULL,
                           Xmin, Xmax, Ymin, Ymax, location) {
  
  #~# Downloading NCEP 
  if(is.null(h_NCEP)) {
    h_NCEP <- hourlyNCEP(lat = -54.3, long = 158.57, tme = as.POSIXct(c(date_start, date_end), format = "%d/%m/%Y", tz = 'UTC'))
  }
  
  for(i in grains) {
    
    #~# grain m
    dem_ag <- aggregate(dem_in, (i/5), mean)
    dem_gr_N <-  trim(crop(dem_ag, extent(Xmin, Xmax, Ymin, Ymax)))
    # AWS
    gr_N_aws <- run_micro_model(date_start, date_end, ht = 0.04, useObs = TRUE, dem_gr_N, h_NCEP = h_NCEP)
    save(gr_N_aws, file = paste0("Outputs/Scale/",mth,"_",i,"m_",location,"_aws.Rdata"), compress = 'xz')
    # CRA
    gr_N_cra <- run_micro_model(date_start, date_end, ht = 0.04, useObs = FALSE, dem_gr_N, h_NCEP = h_NCEP)
    save(gr_N_cra, file = paste0("Outputs/Scale/",mth,"_",i,"m_",location,"_cra.Rdata"), compress = 'xz')
    # Nov LR
    gr_N_lr <- run_lapse_rate(date_start, date_end, mth, dem_gr_N)
    save(gr_N_lr, file = paste0("Outputs/Scale/",mth,"_",i,"m_",location,"_lr.Rdata"), compress = 'xz')
    # Summarise
    gr_N_sum <- site_time_series(cra_dat = gr_N_cra$temps, aws_dat = gr_N_aws$temps,
                                 ibut, gr_N_lr, ibut_coord, 
                                 date_start2, date_end, mth, i,
                                 dem_in = dem_gr_N,  
                                 outPathFile = "Outputs/Scale/", 
                                 outPathFig = "Outputs/Scale/Figures/")
    save(gr_N_sum, file = paste0("Outputs/Scale/",mth,"_",i,"m_summary_",location,".Rdata"), compress = 'xz')
    
  }
 
  return(NULL)
  
}

#' Extract brms regression coefficients
#' 
#' @param mod brms model object.
#' @param quant Temperature qunatile (Tmin, Tmean, Tmax)
#' 
#' @return dataframe
#' @export 
brms_est <- function(mod, quant) {
  
  fixed <- summary(mod)$fixed
  random <- summary(mod)$random[[1]]
  df <- rbind.data.frame(fixed,random) %>%
    rownames_to_column(var = "Var") %>%
    dplyr::select('Var','Estimate', 'l-95% CI', 'u-95% CI') %>%
    dplyr::rename('l95' = 'l-95% CI', 'u95' = 'u-95% CI') %>%
    mutate(Var = factor(Var, levels = rev(Var))) %>%
    mutate(Model = quant) %>%
    mutate(fname = paste0(round(Estimate,2), 
                          " (", round(l95,2), ", ", round(u95,2), ")"))
  
}

