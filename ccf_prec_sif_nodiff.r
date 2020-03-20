# Created on Oct 28, 2018
# to run CCF for Prec and SIF

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

# RADIATION data for resampling to 1 degree
radnc <- "/depot/jsdukes/data/Maria/radiation/ForSIF/CERES_EBAF-Surface_Ed4.0_Subset_200701-201712.nc" # clm gpp netcdf file name
radin <- brick(radnc)  # read as a brick of rasters because it has a third dimension time
#radin # my clm gpp file is called clmgppin
#print(clmgppin)
radin <- rotate(radin[[1]])

# GOME-2 v27
sif.files <- list.files(path="/depot/jsdukes/data/Maria/sif/SIF_Mary_OneDrive_1_10-31-2018/gome2/MetopA_Level3/", pattern = "*.nc$", full.names = TRUE)
sif.mary <- stack(sif.files)  # Jan 2007 to Dec 2017
sif.mary[sif.mary == -999.00] <- NA
sif.mary1deg <- resample(sif.mary, radin)
sif.arrays <- as.array(sif.mary1deg)

#Create the stack using the resampled rasters:
trmmres <- list.files(path="/depot/jsdukes/data/Maria/trmm_tmpa_download/prec_trmmresampled_forsif_1deg/", pattern = "*resampledsif1deg.gri$", full.names = TRUE)  #path of the directory where the new files are saved
trmmres_st <- stack(trmmres[85:216])   #create a stack of rasters with the new files from 2007 to 2017
prec.arrays <- as.array(trmmres_st)

ci <- 0.95
n.used <- 132   #calculated from ccfpregpp$n.used
significance_level <- qnorm((1 + ci)/2)/sqrt(n.used)
listofdates <- seq(from=as.Date("2007-01-01"), to=as.Date("2017-12-01"), by="month")

#=== Parallel loop ===
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

write.table(precsifpeaks.df, 'ccf_precsif.csv', sep=",", row.names=FALSE)

