# Created on Oct 28, 2018
# Based on the script for CCF Prec-GPP Zhang in the cluster created by MRU on May 3, 2018
# to run CCF for Prec and Maiac EVI
# Modified on Feb 2, 2019 to use arrays, take only negative timelags, change the NA thershold and remove the log() in prec
# Modified on Feb 14, 2019 to run for SIF

library("ncdf4") 
library(zoo)
library(raster)
library(rasterVis)
library(gpclib)
library(maptools)
library(maps)
library(forecast)
library(astsa)
library(TSA)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(itertools)
library(lubridate)


# netcdf:
# ncname <- "/Volumes/MRUDATA/Data/SIFTER/SIF_Amazone_Koren2018"
# ncfname <- paste(ncname, ".nc", sep = "")
# ncin <- nc_open(ncfname)
# print(ncin) 
# 
# # SIFTER:
# sifter.files <- list.files(path="/Volumes/Maria/sif/orig_data/", pattern = "*.nc$", full.names = TRUE)
# sifter.stack <- stack(sifter.files)
# plot(sifter.stack[[2]]*1000)
# sifter.mask <- mask(sifter.stack, continents)
# plot(sifter.mask*1000)
# # the data with the Amazon mask:
# sifter.amazon.file <- stack("/Volumes/MRUDATA/Data/SIFTER/SIF_Amazone_Koren2018.nc")
# plot(sifter.amazon.file[[1]] )
# 
# # GOME-2 GFZ (daily, 2007-2016, 0.5x0.5):
# sif.gfz <- stack("/Volumes/MRUDATA/Data/SIF/gome2/GFZ/GOME-2_SIF_2007.nc")
# sif.gfz[sif.gfz < -999] <- NA
# plot(sif.gfz[[28]])
# sif.gfz.month2 <- mean(sif.gfz[[32:60]], na.rm=TRUE)
# 
# # SCIAMACHY
# sif.sci <- stack("/Volumes/MRUDATA/Data/SIF/SCIA/SCIA_SIF_08_2002-01_2012.nc")
# sif.sci[sif.sci < -999] <- NA
# plot(sif.sci[[1]])

# RADIATION data for resampling to 1 degree
radnc <- "/depot/jsdukes/data/Maria/radiation/ForSIF/CERES_EBAF-Surface_Ed4.0_Subset_200701-201712.nc" # clm gpp netcdf file name
radin <- brick(radnc)  # read as a brick of rasters because it has a third dimension time
#radin # my clm gpp file is called clmgppin
#print(clmgppin)
radin <- rotate(radin[[1]])


# GOME-2 v27 (Mary)
sif.files <- list.files(path="/depot/jsdukes/data/Maria/sif/SIF_Mary_OneDrive_1_10-31-2018/gome2/MetopA_Level3/", pattern = "*.nc$", full.names = TRUE)
sif.mary <- stack(sif.files)  # Jan 2007 to Dec 2017
sif.mary[sif.mary == -999.00] <- NA
sif.mary1deg <- resample(sif.mary, radin)
sif.arrays <- as.array(sif.mary1deg)
# look at stuff:
# plot(sif.mary)
# continents <- readOGR("/Volumes/MRUDATA/Data/LandCover/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp")
# sif.mary.crop <- crop(sif.mary, extent(continents))
# sif.mary.mask <- mask(sif.mary.crop, continents)
# plot(sif.mary.mask)


#precipitation:
# precip.files <- list.files(path="/Volumes/MRUDATA/Time series work/precipitation_trmm3b43/02072018/trmm_tmpa_download/downloads/", pattern = "*.nc4$", full.names = FALSE)
# preciptest <- raster("/Volumes/MRUDATA/Time series work/precipitation_trmm3b43/02072018/trmm_tmpa_download/downloads/3B43.20000101.7A.HDF.nc4")
# for(f in precip.files) {
#   trmm_infile <- sprintf("/Volumes/MRUDATA/Time series work/precipitation_trmm3b43/02072018/trmm_tmpa_download/downloads/%s", f)  #path of the precipitation file
#   trmmin_raster <- raster(trmm_infile) # read precipiation netcdf file as a raster
#   trmmin.t = t(trmmin_raster)  # t() for transpose
#   trmmin.flipy = flip(trmmin.t, direction = 2) # flip(,direction=2) to flip my data in the y direction
#   trmmin.t.flipxy = flip(trmmin.flipy, direction = 1)  # flip(,direction=1) to flip my data in the x direction
#   crs(trmmin.t.flipxy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   #set projection
#   basicfilename<-gsub(".nc4.*","",f)   #string using the "f" to use later
#   resampledprecname <- sprintf("/Volumes/MRUDATA/Time series work/precipitation_trmm3b43/02072018/trmm_tmpa_download/downloads/2000_2017_forsif/%s_resampledsif.nc4", basicfilename) #path of the dierctory and filename to use to save the new files
#   new.rasterA = projectRaster(trmmin.t.flipxy, sif.mary, filename=resampledprecname) #define the precipiation rasters projection and extent as it is in the sif raster and save the new raster in a new file using the name and path above
# }
#Create the stack using the resampled rasters:
trmmres <- list.files(path="/depot/jsdukes/data/Maria/trmm_tmpa_download/prec_trmmresampled_forsif_1deg/", pattern = "*resampledsif1deg.gri$", full.names = TRUE)  #path of the directory where the new files are saved
trmmres_st <- stack(trmmres[85:216])   #create a stack of rasters with the new files from 2007 to 2017
prec.arrays <- as.array(trmmres_st)

# location stuff:
#cells <- cellFromRowCol(trmmres_st, rownr = 1:180, colnr=1:360)
#xytrops <- xyFromCell(trmmres_st,cells)
#xytrops
#tropics <- which(xytrops[,2] > -20 & xytrops[,2] < 20, arr.ind=TRUE)
#tropicslatitudes <- xytrops[194,2]  #tropics at 0.5 are [128:207,2] and at 1deg are [71:110,2]
#tropicslongitudes <- xytrops[250,1]  # 1:720

#=== For loop to run the CCF in all the tropical pixels ===
ci <- 0.95
n.used <- 132   #calculated from ccfpregpp$n.used
significance_level <- qnorm((1 + ci)/2)/sqrt(n.used)

listofdates <- seq(from=as.Date("2007-01-01"), to=as.Date("2017-12-01"), by="month")

#=== parallel! ===
detectCores()
nCores <- 4 #manually for non-cluster machines
#nCores <- as.numeric(Sys.getenv('NSLOTS')) #for use on cluster
cl <- makeCluster(4)
registerDoParallel(cl)
#cl
#print('im fine til here')
getDoParWorkers()


precsifpeaks.df <- foreach(  #seq_len(dim(prec.arrays)[2]), 
  chunk = isplitVector(seq_len(dim(prec.arrays)[2]), chunks=getDoParWorkers()),.combine='rbind',.packages=c("foreach")
                       ) %dopar% {
  foreach(b=chunk, .combine='rbind', .packages=c('sp','raster','forecast','TSA')) %:%  #first or outer "vector" to iterate on: longitudes vector
    foreach(a=c(71:110), .combine='rbind') %dopar% {   # second or inner "vector" to iterate on: lats
    lonlat <- cbind(b,a)
    
    trmmppttrop <- prec.arrays[a, b, ]    #before trmmppttrop <- extract(trmmres_st, lonlat) 
    gpptrop <- sif.arrays[a, b, ]    #before gpptrop <- extract(evi.resampled.ras, lonlat)
    
    if (length(gpptrop[is.na(gpptrop)==TRUE]) > (length(gpptrop)/3) | round(sum(gpptrop, na.rm=TRUE)) == 0){
      print("lots of nas")
      peakcorr <- NA
      lineofpeak <- NA
      lagofpeak <- NA
      sig <- NA
    }
    else if (sum(gpptrop==0,na.rm=TRUE)> (length(gpptrop)/3)){
      print("lots of zeros")
      peakcorr <- NA
      lineofpeak <- NA
      lagofpeak <- NA
      sig <- NA
    }
    else{
      print("not many nas or zeros")
      trmmppttrop.ts <- ts(trmmppttrop, start=2007, frequency = 12)
      gppfinal.ts <- ts(gpptrop, start=2007, frequency = 12)
      
      pptspline <- zoo::na.spline(trmmppttrop.ts, maxgap=3)
      gppspline <- zoo::na.spline(gppfinal.ts, maxgap=3)
      pptspline_ts <- ts(pptspline, start=2007, frequency = 12)
      gppspline_ts <- ts(gppspline, start=2007, frequency = 12)
      
      if(length(pptspline[is.na(pptspline)==TRUE]) > length(pptspline)/3 | length(gppspline[is.na(gppspline)==TRUE]) > length(gppspline)/3){
        peakcorr <- NA
        lineofpeak <- NA
        lagofpeak <- NA
        sig <- NA
      }
      
      else{
        
        arimamodel_ppt <- auto.arima(log(pptspline_ts+0.000001))
        prewhitened_pptgpp <- tryCatch(prewhiten(x=log(pptspline_ts+0.000001), y=gppspline_ts, x.model=arimamodel_ppt, lag.max=6, plot=FALSE), error=function(s) NA)
        
        if(is.na(prewhitened_pptgpp)==TRUE){
          peakcorr <- NA
          lineofpeak <- NA
          lagofpeak <- NA
          sig <- NA
        }
        else{
          #plot(autoprewhiten2$ccf, main="auto2", ci.type="white")
          ppt.gpp.ccf.df1 <- data.frame(lag=prewhitened_pptgpp$ccf$lag, CCF=prewhitened_pptgpp$ccf$acf)
          significants <- which(abs(ppt.gpp.ccf.df1$CCF[3:7]) > abs(significance_level))
          
          if (length(significants) > 0){
            #print("significant")
            peakcorr <- ppt.gpp.ccf.df1$CCF[2 + which.max( abs(ppt.gpp.ccf.df1$CCF[3:7]) )]
            lineofpeak <- which(ppt.gpp.ccf.df1$CCF == peakcorr, arr.ind=TRUE)
            lagofpeak <- ppt.gpp.ccf.df1[lineofpeak,1]
            sig <- 1
          } else{
            #print("not significant")
            peakcorr <- ppt.gpp.ccf.df1$CCF[2 + which.max( abs(ppt.gpp.ccf.df1$CCF[3:7]) )]
            lineofpeak <- which(ppt.gpp.ccf.df1$CCF == peakcorr, arr.ind=TRUE)
            lagofpeak <- ppt.gpp.ccf.df1[lineofpeak,1]
            sig <- 0
          }
        }
      }
    }
    corename <- paste(Sys.info()[['nodename']],Sys.getpid(),sep=',')
    data.frame(lonlat,peakcorr=peakcorr, lagofpeak=lagofpeak, sig=sig)
  }
}

write.table(precsifpeaks.df, 'ccf_precsif_1360_2_nodiff.csv', sep=",", row.names=FALSE)

