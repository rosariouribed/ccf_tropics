# Created by MRU on Oct 28, 2018
# to run CCF for Prec and SIF

#======== Load libraries ========

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

#======== Load and prepare data ========

## Contact MRU at uribem@purdue.edu about data availability
# Precipitation data:
trmmres <- list.files(path="/depot/jsdukes/data/Maria/trmm_tmpa_download/prec_trmmresampled_forsif_1deg/", pattern = "*resampledsif1deg.gri$", full.names = TRUE)  #path of the directory where the new files are saved
trmmres_st <- stack(trmmres[85:216])   # Raster stack Jan 2007 to Dec 2017
prec.arrays <- as.array(trmmres_st)   # transform to 2d matrix for faster computations

# SIF from GOME-2 v27:
sif.files <- list.files(path="/depot/jsdukes/data/Maria/sif/SIF_Mary_OneDrive_1_10-31-2018/gome2/MetopA_Level3/", pattern = "*.nc$", full.names = TRUE)
sif.gome <- stack(sif.files)  # Raster stack Jan 2007 to Dec 2017
sif.gome[sif.gome == -999.00] <- NA
sif.gome1deg <- resample(sif.gome, trmmres_st)
sif.arrays <- as.array(sif.gome1deg)  # transform to 2d matrix for faster computations


#======== For loop with CCF analysis ========

# Explore and define some stuff for parallel for loop:
detectCores()
nCores <- 4 #manually for non-cluster machines
cl <- makeCluster(nCores)
registerDoParallel(cl)
getDoParWorkers()

# Info for confidence intervals in CCF:
ci <- 0.95
n.used <- 132   #calculated from ccfpregpp$n.used
significance_level <- qnorm((1 + ci)/2)/sqrt(n.used)
listofdates <- seq(from=as.Date("2007-01-01"), to=as.Date("2017-12-01"), by="month")

# Loop for CCF analysis:
precsifpeaks.df <- foreach(  #seq_len(dim(prec.arrays)[2]), 
  chunk = isplitVector(seq_len(dim(prec.arrays)[2]), chunks=getDoParWorkers()),.combine='rbind',.packages=c("foreach")
                       ) %dopar% {
  foreach(b=chunk, .combine='rbind', .packages=c('sp','raster','forecast','TSA')) %:%  #first or outer "vector" to iterate on longitudes vector (array columns)
    foreach(a=c(71:110), .combine='rbind') %dopar% {   # second or inner "vector" to iterate on latitudes. #71:110 are the tropics latitudes rows
    lonlat <- cbind(b,a)
    
    trmmppttrop <- prec.arrays[a, b, ]    # similar to using extract(trmmres_st, lonlat) in a raster, but faster!
    gpptrop <- sif.arrays[a, b, ]
    
    # When lot's of NAs or zeros in the data:
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
    # When not too many NAs or zeros in the data:
    else{
      print("not many nas or zeros")
      trmmppttrop.ts <- ts(trmmppttrop, start=2007, frequency = 12)   # Read as time series object
      gppfinal.ts <- ts(gpptrop, start=2007, frequency = 12)
      pptspline <- zoo::na.spline(trmmppttrop.ts, maxgap=3)   # Fill short gaps in the time series using spline
      gppspline <- zoo::na.spline(gppfinal.ts, maxgap=3)
      pptspline_ts <- ts(pptspline, start=2007, frequency = 12)   # Read as time series object
      gppspline_ts <- ts(gppspline, start=2007, frequency = 12)
      
      # When still NAs after spline:
      if(length(pptspline[is.na(pptspline)==TRUE]) > length(pptspline)/3 | length(gppspline[is.na(gppspline)==TRUE]) > length(gppspline)/3){
        peakcorr <- NA
        lineofpeak <- NA
        lagofpeak <- NA
        sig <- NA
      }
      
      # When not too many NAs after spline:
      else{ 
        arimamodel_ppt <- auto.arima(log(pptspline_ts+0.000001))   # Fit the ARIMA model to the climate variable
        prewhitened_pptgpp <- tryCatch(prewhiten(x=log(pptspline_ts+0.000001), y=gppspline_ts, x.model=arimamodel_ppt, lag.max=6, plot=FALSE), error=function(s) NA)   # Prewhiten and run CCF using the identified ARIMA model
        # If prewhitening didn't work and CCF couldn't be run for pixel:
        if(is.na(prewhitened_pptgpp)==TRUE){
          peakcorr <- NA
          lineofpeak <- NA
          lagofpeak <- NA
          sig <- NA
        }
        # If prewhitening worked and the CCF was run:
        else{
          ppt.gpp.ccf.df1 <- data.frame(lag=prewhitened_pptgpp$ccf$lag, CCF=prewhitened_pptgpp$ccf$acf) # Dataframe with lags and correlation coefficients
          significants <- which(abs(ppt.gpp.ccf.df1$CCF[3:7]) > abs(significance_level))  # Check if there were significant correlations
          if (length(significants) > 0){
            print("significant")
            peakcorr <- ppt.gpp.ccf.df1$CCF[2 + which.max( abs(ppt.gpp.ccf.df1$CCF[3:7]) )]   # Strongest correlation
            lineofpeak <- which(ppt.gpp.ccf.df1$CCF == peakcorr, arr.ind=TRUE)
            lagofpeak <- ppt.gpp.ccf.df1[lineofpeak,1]  # Lag of strongest correlation
            sig <- 1  # Flag as significant
          } else{
            print("not significant")
            peakcorr <- ppt.gpp.ccf.df1$CCF[2 + which.max( abs(ppt.gpp.ccf.df1$CCF[3:7]) )]   # Strongest correlation
            lineofpeak <- which(ppt.gpp.ccf.df1$CCF == peakcorr, arr.ind=TRUE)
            lagofpeak <- ppt.gpp.ccf.df1[lineofpeak,1]  # Lag of strongest correlation
            sig <- 0  # Flag as not significant
          }
        }
      }
    }
    corename <- paste(Sys.info()[['nodename']],Sys.getpid(),sep=',')
    data.frame(lonlat,peakcorr=peakcorr, lagofpeak=lagofpeak, sig=sig)  # Add line to results data frame
  }
}

write.table(precsifpeaks.df, 'ccf_precsif.csv', sep=",", row.names=FALSE)

