# compute attenuation rates based on Lagrangian comparisons
lagrangianComparison <- function(reachInfoPath,
                                 sampleTimesPath,
                                 sampleConcentrationsPath,
                                 pInterfPath,
                                 outTable,
                                 dateTimeFormat,
                                 stationUp,
                                 stationDown){
  # check order of stations
  if(stationUp >= stationDown)
    stop('stations must be numbered in downstream direction, stationUp must have smaller index as stationDown')

  # check that stations are consecutive
  if((stationDown - stationUp) > 1)
    stop('stations must be consecutive: stationUp = i, stationDown = i+1')
    
  # reach information (reach lengths, flow velocities and travel times)
  reachInfo <- dplyr::tbl_df(read.table(reachInfoPath, 
                                        header=TRUE, 
                                        sep=';',
                                        colClasses='numeric'))

  # read and format interference probabilities computed with probInterf.R
  probInterf <- dplyr::tbl_df(read.table(pInterfPath,
                                         header=TRUE,
                                         sep=';',
                                         stringsAsFactors=FALSE))
  probInterf <- probInterf %>% dplyr::mutate_if(is.integer, .funs=as.double)
  
  # read and format concentrations in samples
  sampleConcentrations <- dplyr::tbl_df(read.table(sampleConcentrationsPath,
                                                   header=TRUE,
                                                   sep=';',
                                                   colClasses=c('numeric',
                                                                'numeric',
                                                                'character',
                                                                'numeric'),
                                                   stringsAsFactors=FALSE))
  
  if(stationUp >= max(sampleConcentrations$station))
    stop('stationUp cannot be at downstream end of reach')
  
  # read and format sample information
  sampleTimes <- dplyr::tbl_df(read.table(sampleTimesPath, sep=';', 
                                          header=TRUE,
                                          colClasses=c('numeric',
                                                       'character', 
                                                       'character',
                                                       'numeric',
                                                       'numeric'),
                                          stringsAsFactors=FALSE))
  sampleTimes$tBeg <- as.POSIXct(sampleTimes$tBeg, format=dateTimeFormat)
  sampleTimes$tEnd <- as.POSIXct(sampleTimes$tEnd, format=dateTimeFormat)
  
  # build sample IDs based on station and sample number
  sampleTimes$sampleID <- paste(sampleTimes$station, sampleTimes$sample, sep='_')
  sampleConcentrations$sampleID <- paste(sampleConcentrations$station, 
                                         sampleConcentrations$sample, sep='_')
  
  # join sampleConcentrations, sampleTimes and probInterf
  sampleConcTime <- sampleTimes %>%
    dplyr::inner_join(sampleConcentrations, 
                      by='sampleID') %>%
    dplyr::inner_join(probInterf, by='sampleID') %>%
    dplyr::mutate(tMid=0.5*(as.numeric(tBeg) + as.numeric(tEnd)),
                  station=station.x,
                  sample=sample.x,
                  pInterfMax=apply(X=probInterf[, -grep(pattern='sampleID', 
                                                       x=names(probInterf))], 
                                   MARGIN=1, 
                                   FUN=max)) %>%
    dplyr::select(station, sample, sampleID, tBeg, tEnd, tMid,
                  dx0, variable, value, pInterfMax)

  # grab data for each station separately
  waveUp <- sampleConcTime %>% 
    dplyr::filter(station==stationUp) %>% 
    dplyr::arrange(tMid)
  waveDown <- sampleConcTime %>% 
    dplyr::filter(station==stationDown) %>% 
    dplyr::arrange(tMid)

  # travel-time conform comparison and rate computation based on
  # k = ((c0 - cInterp)/c0)/t:

  # empty data.frame's to store results
  kUpDown <- as.data.frame(matrix(nrow=nrow(waveUp), 
                             ncol=8, 
                             dimnames=list(NULL, 
                                           c("k", "shiftedSample", "left", "right", 
                                             "interfSample", "interfLeft", "interfRight",
                                             "tStart"))))
  kDownUp <- as.data.frame(matrix(nrow=nrow(waveDown), 
                                  ncol=8, 
                                  dimnames=list(NULL, 
                                                c("k", "shiftedSample", "left", "right", 
                                                  "interfSample", "interfLeft", "interfRight",
                                                  "tStart"))))
  
  # travel time in reach [hours]
  dt <- reachInfo$dt[reachInfo$reach==stationUp] 
  
  # downstream direction:
  
  # loop over samples in waveUp
  for(i in 1:nrow(waveUp)){
    
    # shift samples downstream using corresponding travel time
    B1B2 <- waveUp$tMid[i] + dt*3600
    
    # is sample shifted from waveUp surrounded by samples in waveDown?
    # if so, interpolate concentration between left and right neighbors in waveDown 
    # and compute k
    loc <- B1B2 - waveDown$tMid
    if(length(which(loc>0))>0 & length(which(loc<0))>0){
      left <- waveDown[max(which(loc>0)), ]
      right <- waveDown[max(which(loc>0))+1, ]
      if(!is.na(left$value) & !is.na(right$value)){
        B2 <- rbind(left, right)
        B1B2 <- as.data.frame(approx(x=B2$tMid, y=B2$value, xout=B1B2),
                              col.names=c("tMid", "value"))
        kUpDown$k[i] <- ((B1B2$value - waveUp$value[i])/waveUp$value[i])/dt
      }
      kUpDown$shiftedSample[i] <- waveUp$sampleID[i]
      kUpDown$interfSample[i] <- waveUp$pInterfMax[i]
      kUpDown$left[i] <- left$sampleID
      kUpDown$right[i] <- right$sampleID
      kUpDown$interfLeft[i] <- left$pInterfMax
      kUpDown$interfRight[i] <- right$pInterfMax
      kUpDown$tStart[i] <- waveUp$tMid[i]
    }else{
      kUpDown$shiftedSample[i] <- waveUp$sampleID[i]
    }
  }
  kUpDown$shiftDirection <- 'downstream'
  
  # repeat in upstream direction:
  for(i in 1:nrow(waveDown)){
    B2B1 <- waveDown$tMid[i] - dt*3600
    loc <- B2B1 - waveUp$tMid
    if(length(which(loc>0))>0 & length(which(loc<0))>0){
      left <- waveUp[max(which(loc>0)), ]
      right <- waveUp[max(which(loc>0))+1, ]
      if(!is.na(left$value) & !is.na(right$value)){
        B1 <- rbind(left, right)
        B2B1 <- as.data.frame(approx(x=B1$tMid, y=B1$value, xout=B2B1),
                              col.names=c("tMid", "value"))
        kDownUp$k[i] <- ((waveDown$value[i] - B2B1$value)/B2B1$value)/dt
      }
      kDownUp$shiftedSample[i] <- waveDown$sampleID[i]
      kDownUp$interfSample[i] <- waveDown$pInterfMax[i]
      kDownUp$left[i] <- left$sampleID
      kDownUp$right[i] <- right$sampleID
      kDownUp$interfLeft[i] <- left$pInterfMax
      kDownUp$interfRight[i] <- right$pInterfMax
      kDownUp$tStart[i] <- waveDown$tMid[i] - dt*3600
    }else{
      kDownUp$shiftedSample[i] <- waveDown$sampleID[i]
    }
  }
  kDownUp$shiftDirection <- 'upstream'
  
  # bind all k data frames and re-format tStart
  kk <- kUpDown %>%
    rbind(kDownUp) %>%
    dplyr::mutate(tStart=as.POSIXct(tStart, origin='1970-01-01 00:00:00'))

  # write output table
  if(!is.na(outTable))
    write.table(kk, file=outTable, 
                quote=FALSE, 
                sep=";",
                col.names=TRUE, 
                row.names=FALSE)
  
  
}
