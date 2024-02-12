# SCRIPT TO REPRODUCE THE WORK IN MENGOLI ET. AL. 2024 GBC PAPER SUBMISSION

rm(list = ls())
cat('\014')
# library(netcdf)----
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(lubridate)


# DATA WITH SIMULATED SOIL MOISTURE FROM SPLASHv1 UPLOADING----
# SIMULATED DATA FROM: https://github.com/LEMONTREE-Land-Ecosystem-Models/global_subdaily_gpp
# SET DIRECTORY
setwd('C:/Desktop/third_peper/DavidO_GPP_global/daily_soilMoist_HPC')

# OPEN DATA
  nc_dataSM <- nc_open('soil_moisture_2010.nc')
  # nc_dataSM <- nc_open('soil_moisture_2011.nc')
  # nc_dataSM <- nc_open('soil_moisture_2015.nc')
  # nc_dataSM <- nc_open('soil_moisture_2016.nc')

# extract one variable
  thetaWrel <- ncvar_get(nc_dataSM, "soilmstress_mengoli")
  
  # conversion (relative moisture --> water content)
  thetaWp = 0.151   #m3/m3
  theta <- (thetaWp + 0.15*thetaWrel)
  
  nc_close(nc_dataSM)
  rm(nc_dataSM, thetaWrel)
  
# CLAY CC DATA UPLOADING ----
  setwd('C:/Desktop/third_peper/DavidO_GPP_global/')
  nc_cc <- nc_open('collo_fac.nc')

# lon <- ncvar_get(nc_cc, "longitude")
# lat <- ncvar_get(nc_cc, "latitude")
  crs <- ncvar_get(nc_cc, "layer")
  
  nc_close(nc_cc)
  rm(nc_cc)

# NEW VARIABLE CREATION IN NETCDF----
# DATA WITH SIMULATED DAILY GPP FROM SUB-DAILY PMODEL UPLOADING----
# SIMULATED DATA FROM: https://github.com/LEMONTREE-Land-Ecosystem-Models/global_subdaily_gpp
# SET DIRECTORY  
setwd('C:/Desktop/third_peper/DavidO_GPP_global/daily_GPP_HPC/2010_2019')

# open data
  nc_data <- nc_open('daily_gpp_2010.nc')
  # nc_data <- nc_open('daily_gpp_2011.nc')
  # nc_data <- nc_open('daily_gpp_2015.nc')
  # nc_data <- nc_open('daily_gpp_2016.nc')

# extract one variable
  var <- ncvar_get(nc_data, "subdaily_gpp_with_mengoli_beta")
  
  nc_close(nc_data)
  rm(nc_data)

# new var1 (BE) each day
  var1 <- 0.47 * var
  
# new var2 (BE) over the year
  var2 <- rowSums(var1[,,1:365], dims=2, na.rm=T)             
  var2_na = is.na(var1)
  var2_na <- rowSums(var2_na, dims=2)
  pos_na = which(var2_na == 365)
  if(length(pos_na)>0) var2[pos_na] = NA
  rm(var2_na, pos_na)
  
  
# MEAN DAILY TEMPERATURE CALCULATION TO COMPUTE var3 (ks)----
# METEO DATA FROM: https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.20d54e34?tab=overview
# SET DIRECTORY    
setwd('C:/Desktop/data_for_Carlo_global/Tair_K/2010')
  
# importo i file con le temperature e calcolo la media
# per ogni pixel per tutte le ore del giorno
# MEAN DAILY TEMP CALCULATION
Tair_mean = var 
conta_strati = 1
for (ncT in sort(list.files(pattern = 'Tair_WFDE5_CRU',full.names = F,recursive = T))) {
  nc_dataT = nc_open(ncT)
  Tair <- ncvar_get(nc_dataT, "Tair")
  
  for (d in seq(1,dim(Tair)[3], by = 24)){
    Tair_mean2 = rowMeans(Tair[,,d:(d+23)], dims=2, na.rm =T) - 273.15              
    var1b = Tair[,,d:(d+23)]
    var2_na = is.na(var1b)
    var2_na <- rowSums(var2_na, dims=2)
    pos_na = which(var2_na == dim(var1b)[3])
    if(length(pos_na)>0) Tair_mean2[pos_na] = NA
    rm(var2_na)
    Tair_mean[,,conta_strati] = Tair_mean2
    conta_strati = conta_strati + 1 
  }
  
  nc_close(nc_dataT)
  rm(nc_dataT,Tair)
  rm(var1b,Tair_mean2)
}
rm(ncT, d)
rm(conta_strati)


# TEMPERATURE FUNCTION CALCULATION (TO THEN COMPUTE ks)-----

for (Q10 in c (1.00)){
  
  # TEST SENSITIVITY Q10 TERM----
  # temperature function
  f_t = Tair_mean 
  
  for (n in seq(1,dim(Tair_mean)[3])){
    f_t[,,n] = (Q10)^((Tair_mean[,,n] - 15)/10)
  }
  rm(n)
  
}
  
# MOISTURE FUNCTION COMPUTATION TO THEN COMPUTE ks-----
    # moisture function
    # with a constant soil porosity
    phi = 0.463
    thetaOpt = 0.65 *phi
    a1 = crs
  
    fm1 <- var
    for (n in seq(1,dim(var)[3])){
  
      theta1 = theta[,,n]
  
      fm1d = fm1[,,n]
  
      pos = which(theta1 < thetaOpt)
      if (length(pos) > 0)
        fm1d[pos] = ((0.1 + thetaOpt) / (0.1 + theta1[pos])) * (theta1[pos]/thetaOpt)^(1+ 2*a1[pos])
  
      pos = which(theta1 >= thetaOpt)
      if (length(pos) > 0)
        fm1d[pos] = ((phi - theta1[pos])/(phi - thetaOpt))
  
      pos = which(fm1d < 0)
      if (length(pos) > 0)
        fm1d[pos] = 0
  
      pos = which(theta1 >= thetaOpt)
      if (length(pos) > 0)
        fm1d[pos] = fm1d[pos]^(0.75)
  
      fm1[,,n] = fm1d
  
    }
  
  rm(n, pos, a1, fm1d, theta1)
  
 
# new var3 (ks) each day
    var3 = f_t
    for (n in seq(1,dim(var3)[3]))
  	var3[,,n] <- f_t[,,n] * fm1[,,n]
    
    rm(fm1, f_t, n)
    
# new parameter value over the year
# var_gt <- rowSums(var[,,1:365], dims=2, na.rm=T)
    var_gt <- var2
    var_ks <- rowSums(var3[,,1:365], dims=2, na.rm=T)          
    
    var1b = var3
    var2_na = is.na(var1b)
    var2_na <- rowSums(var2_na, dims=2)
    pos_na = which(var2_na == dim(var1b)[3])
    if(length(pos_na)>0) var_ks[pos_na] = NA
    rm(var2_na,var1b)
    
    Csoil <- var_gt /var_ks
    
  rm(var_gt, var_ks, pos_na)
    
# new var4 (RH) each day
    var4 = var3
    for (n in seq(1, dim(var3)[3]))
      var4[,,n] <- Csoil *var3[,,n]
    rm(n)
    
# new var5 (NEE) each day
    var5 <- var4 - var1            
    
  # UPSCALE NEE (here called GPP_agg) from 0.5° to 1°-----
    GPP_agg = c()
    for (t in seq(1, dim(var5)[3])){
      var_x = var5[,,t]
      row2 = c()
      for (n2 in seq (1, ncol(var_x), by = 2)){
        row1 = c()
        for (n1 in seq(1, nrow(var_x), by =2)){
          row1 = c(row1, mean(var_x[n1:(n1+1),n2:(n2+1)], na.rm=T) ) 
        }
        row2 = c(row2, row1)
      }
     GPP_agg = c(GPP_agg,row2)
    }
    # print('lunghezza vettore GPP_agg')
    # print(length(GPP_agg))
    GPP_agg = array(GPP_agg, dim = c(nrow(var5)/2,ncol(var5)/2,dim(var5)[3]))
    rm(row1, row2, t,var_x, n1, n2)
    
    # fixing one issue
    GPP_agg[,, 365] = GPP_agg[,, 364]
    # GPP_agg[,, 366] = GPP_agg[,, 365]  # year 2016
  
  # MONTHLY AGGREGATION-----
    giorni_del_mese = c(31, 28, 31, 30, 31, 30, 31, 31,30,31,30,31)
    if (dim(GPP_agg)[3] == 366)
      giorni_del_mese = c(31, 29, 31, 30, 31, 30, 31, 31,30,31,30,31)
   
    Var_mean = GPP_agg
    
    inizio_mese = 1
    cont_mese = 1
      for (d in giorni_del_mese){
        fine_mese = inizio_mese +d -1
        print(inizio_mese)
        print(fine_mese)
        
        Tair_mean2 = rowSums (GPP_agg[,,inizio_mese:fine_mese], dims =2, na.rm=T)
        var1b = GPP_agg[,,inizio_mese:fine_mese]
        var2_na = is.na(var1b)
        var2_na <- rowSums(var2_na, dims=2)
        pos_na = which(var2_na == dim(var1b)[3])
        if(length(pos_na)>0) Tair_mean2[pos_na] = NA
        rm(var2_na)
        
        Var_mean[,,cont_mese] = Tair_mean2  
        inizio_mese = fine_mese +1
        cont_mese = cont_mese +1
        
      }
    rm(cont_mese, inizio_mese, fine_mese, giorni_del_mese, d)
    rm(GPP_agg, Tair_mean2, pos_na, var1b)
    
    Var_mean = Var_mean[,,1:12]
    
# END OF THE SCRIPT; FINAL OUTPUT: MONTHLY NEE at 1°x 1°
    
# SIMPLE GRAPHS----- 
  library(ncdf4) # package for netcdf manipulation
  library(raster) # package for raster manipulation
  library(rgdal) # package for geospatial analysis
  library(ggplot2) # package for plotting
    
# INVERSION PRODUCT (OBSERVATIONS) LOADING
# DATA DOWNLOADED FROM: https://zenodo.org/records/5829774; Jiang et al., 2022
  
  setwd('C:/Desktop/third_peper/NEE_DATA/GCAS2021_gridded_fluxes')
  
  nc_dataNEE <- nc_open('GCAS2021_monthlyflux_1x1.nc')

  lon <- ncvar_get(nc_dataNEE, "lon")
  lon[lon > 180] <- lon[lon > 180] - 360 
  
  lat <- ncvar_get(nc_dataNEE, "lat", verbose = F)
  t <- ncvar_get(nc_dataNEE, "datetime")
  
  ndvi.array <- ncvar_get(nc_dataNEE, "bio_opt") # store the data in a 3-dimensional array
 
  dim(ndvi.array) 

  # packages
  library(lattice)
  library(RColorBrewer)
  
  # SELECT 12 MONTHS OF 2010
  ndvi.array = ndvi.array[,,1:12]

# MASK IMPOSITION
pos_na = which(is.na(Var_mean[,,10])==1)   # MONTH: 10 YEAR: 2010
ndvi.slice = ndvi.array[,,10]
ndvi.slice[pos_na] = NA

# VALUES RANGE
# varM = Var_mean[,,10]
# range(varM, na.rm = T)
# ndvi.slice = ndvi.array[,,10]
# range(ndvi.slice, na.rm = T)

# # SIMULATIONS
r <- raster(t(Var_mean[,,10]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')

plot(r, zlim = c(-200, 300), col=(rev(brewer.pal(11,"RdBu"))), box = FALSE) #Dec
grid(nx = NULL, ny = NULL,
     lty = 2, col = "gray", lwd = 1)

title(main = "Dec 2015 NEE SIM", cex.main = 1,   font.main= 2, col.main= "black", line = 0.25)
title(ylab="Latitude", line=2.5, cex.lab=1, font.lab  = 2)
title(xlab="Longitude", line=2.5, cex.lab=1, font.lab  = 2)

# # OBSERVATIONS
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')

plot(r, zlim = c(-200, 300), col=(rev(brewer.pal(11,"RdBu"))), box = FALSE) #Dec
grid(nx = NULL, ny = NULL,
     lty = 2, col = "gray", lwd = 1)

title(main = "Dec 2015 NEE SIM", cex.main = 1,   font.main= 2, col.main= "black", line = 0.25)
title(ylab="Latitude", line=2.5, cex.lab=1, font.lab  = 2)
title(xlab="Longitude", line=2.5, cex.lab=1, font.lab  = 2)




  
