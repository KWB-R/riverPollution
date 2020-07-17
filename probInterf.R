# interference analysis
probInterf <- function(latQpath, latQxPath,
                       latQerror,
                       reachInfoPath,
                       sampleTimesPath,
                       outTable,
                       dateTimeFormat)
{
  library(dplyr)
  
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
  
  # read and format lateral input discharges
  latQ <- dplyr::tbl_df(read.table(latQpath, 
                                   header=TRUE,
                                   sep=';', 
                                   check.names=FALSE,
                                   colClasses='character'))
  latQ$dateTime <- as.POSIXct(latQ$dateTime, format=dateTimeFormat)
  for(i in 2:ncol(latQ)){latQ[, i] <- as.numeric(dplyr::pull(latQ, i))}

  # read lateral input downstream locations
  latQx <- dplyr::tbl_df(read.table(latQxPath, 
                                    header=TRUE, 
                                    sep=';', 
                                    colClasses='character'))
  latQx$dx0 <- as.numeric(latQx$dx0)
  latQx$reach <- as.numeric(latQx$reach)
  
  # reach information (reach lengths, flow velocities and travel times)
  reachInfo <- dplyr::tbl_df(read.table(reachInfoPath, 
                                        header=TRUE, 
                                        sep=';',
                                        colClasses='numeric'))

  # add reach information to latQx
  latQx <- latQx %>% dplyr::inner_join(reachInfo, by='reach')
  
  # make data frame to store arrival times of all samples at all latQ
  arrivalTimes <- matrix(nrow=nrow(sampleTimes), 
                         ncol=2*nrow(latQx),
                         dimnames=list(NULL,
                                       c(rbind(paste0('tBeg_', latQx$id), 
                                               paste0('tEnd_', latQx$id)))))
  
  # loop through all samples and compute arrival times at all latQ
  for(i in 1:nrow(sampleTimes)){
    
    # grab reach of sampling station where current composite sample was taken
    # reach i is downstream of station i
    stationi <- sampleTimes$station[i]
    
    # grab coordinate [m] of station where current composite sample was taken
    xStationi <- sampleTimes$dx0[i]
    
    # accumulated time [hours] in full reaches upstream or downstream of the sampling
    # point needed since total travel time = tFullReaches + tCurrentReach. this is how 
    # we account for the fact that each reach has its own flow velocity
    tFullReaches <- 0
    
    # distance of full reaches upstream or downstream of the sampling point
    dxFullReaches <- 0
    
    # initialize counter for no. full reaches btw. CSO and sampling station
    nr <- 0
    
    # initialize counter for signaling reach change
    nrOld <- 0
    
    # vectors to store arrival times for current composite sample
    tBegj <- rep(NA, times=nrow(latQx))
    tEndj <- rep(NA, times=nrow(latQx))
    
    # loop through all latQx
    for(j in 1:nrow(latQx)){
      
      # grab downstream distance of current latQ [m]
      latQxj <- latQx$dx0[j]
      
      # grab reach of current latQ
      latQreachj <- latQx$reach[j]
      
      # compute distance latQ - sampling station
      dxlatQstation <- latQxj - xStationi
      
      # where is latQ relative to station?
      if(dxlatQstation > 0){ # positive -> latQ is downstream
        
        # no. of full reaches btw. the reach where the current latQ is located
        # and the sampling station
        nr <- abs(stationi - latQreachj)
        
        # are there full reaches btw. station and latQ? -> update tFullReaches and 
        # dxFullReaches
        if(nr != nrOld){
          
          # get travel time and flow distance for intervening full reaches
          tFullReaches <- sum(reachInfo$dt[match(seq(from=stationi, 
                                                     to=(latQreachj-1),
                                                     by=1),
                                                 reachInfo$reach)])
          dxFullReaches <- sum(reachInfo$reachLength[match(seq(from=stationi, 
                                                               to=(latQreachj-1),
                                                               by=1),
                                                           reachInfo$reach)])
          
          # update nrOld
          nrOld <- nr
        }
        
      }else{ # negative -> latQ upstream
        
        # no. of full reaches btw. the reach where the current CSO is located and the 
        # sampling station
        nr <- stationi - latQreachj - 1
        
        # are there full reaches btw. station and latQ? -> update tFullReaches and 
        # dxFullReaches
        if(nr != nrOld){
          
          # get travel time and flow distance for intervening full reaches.
          # for the case of latQ in reach upstream of station (nr=0), the ifelse
          # makes tFullReaches, dxFullReaches and nrOld = 0. then, when the loop gets to 
          # the 1st latQ in the reach downstream of station (nr also = 0 but on the 
          # positive side), no update is made for tFullReaches, dxFullReaches and nrOld,
          # so computations on the positive side continue with tFullReaches, dxFullReaches 
          # and nrOld = 0, until a new downstream reach is reached
          tFullReaches <- ifelse(nr==0, 
                                 0,
                                 -sum(reachInfo$dt[match(seq(from=(latQreachj+1),
                                                             length=nr), 
                                                         reachInfo$reach)]))
          dxFullReaches <- ifelse(nr==0,
                                  0,
                                  -sum(reachInfo$reachLength[match(seq(from=(latQreachj+1), 
                                                                       length=nr), 
                                                                   reachInfo$reach)]))
          
          # update nrOld
          nrOld <- nr
        }
      }
      
      # compute arrival time of sample i at CSO[j]
      tBegj[j] <- sampleTimes$tBeg[i] + tFullReaches*3600 + 
        (dxlatQstation-dxFullReaches)/latQx$v[j]*3600
      tEndj[j] <- sampleTimes$tEnd[i] + tFullReaches*3600 + 
        (dxlatQstation-dxFullReaches)/latQx$v[j]*3600
    }
    
    # store arrival times of sample i at all latQ
    arrivalTimes[i, ] <- c(rbind(tBegj, tEndj))
  }
  
  # convert table of arrival times to data frame
  arrivalTimes <- cbind(sampleID=sampleTimes$sampleID, 
                        data.frame(arrivalTimes), 
                        stringsAsFactors=FALSE)
  
  # join table of composite samples and arrival times
  sampleTimes <- sampleTimes %>%
    dplyr::inner_join(arrivalTimes, by='sampleID')

  # probabilistic approach to interference using density kernels on predicted latQ activity
  # times:
  
  # make time axis extending from beginning to end of sampling period (5-min. resolution)
  allt <- c(sampleTimes$tBeg, latQ$dateTime)
  tt <- seq(from=min(allt), to=max(allt), by=5*60)
  
  # make matrix to store interference likelihoods for all samples at all latQ
  interfProb <- matrix(0, nrow=nrow(sampleTimes), 
                       ncol=nrow(latQx),
                       dimnames=list(sampleTimes$sampleID, 
                                     latQx$id))
  
  # Approach: identify latQ activity times (latQ>0) predicted by hydraulic model and set
  # kernels on this points. Probability density is zero at all remaining times.
  # Using gaussian kernel, "diffuminate" times of latQ activity using model errors 
  # in latQ start and end from comparison with measurements at one of the latQ.
  # Kernel parameter bandwidth ("bw") = prediction error in laQ activity times for
  # latQ start and end added in quadrature = sdError,
  # sdError = sqrt( 1/n * sum(xi - xbar)^2 ); with xi - xbar = hydraulic model errors
  sdError <- sqrt(1/2 * (latQerror[1]^2 + latQerror[2]^2))*3600
  
  # function to integrate probability density curve
  activityProb <- function(data, integrLimits){
    
    # filter data between infiltration limits
    data <- data[data[, 1] >= integrLimits[[1]] &
                   data[, 1] <= integrLimits[[2]], ]
    
    # integrate with trapezoidal rule
    l1 <- data[1:(nrow(data)-1), 2]
    l2 <- data[2:nrow(data), 2]
    hh <- as.numeric(data[2:nrow(data), 1]) - as.numeric(data[1:(nrow(data)-1), 1])
    pp <- sum((l1+l2)/2*hh)
    
    return(pp)
  }
  
  # loop through all samples and all latQ
  for(i in 1:nrow(sampleTimes)){
    
    for(j in 1:nrow(latQx)){
      
      # grab data for current latQ
      latQjName <- latQx$id[j]
      latQj <- latQ[, c("dateTime", latQjName)]
      
      # grab times of latQ activity
      tActivity <- latQj$dateTime[dplyr::pull(latQj, 2)>0]
      
      # if latQ was active at all during sampling, compute probability of interference
      # by integrating kernel-based density:
      
      # latQ at all active?
      if((length(tActivity)>0)){
        
        # latQ active during sampling?
          
          # compute probability density for activity times predicted by model. use
          # kernel as mentioned above
          tActivityDiffuse <- density(x=as.numeric(tActivity), 
                                      bw=sdError, 
                                      adjust=1, 
                                      kernel="gaussian",
                                      cut=3)
          
          # make data frame with output from density calculation
          tActivityDiffuse <- data.frame(dateTime=tActivityDiffuse$x,
                                         probDens=tActivityDiffuse$y)
          
          # extend to cover full samling period
          tActivityDiffuse <- as.data.frame(approx(x=tActivityDiffuse$dateTime,
                                                   y=tActivityDiffuse$probDens, 
                                                   xout=tt),
                                            col.names=c("dateTime", "probDens"))
          tActivityDiffuse[is.na(tActivityDiffuse$probDens), 2] <- 0 
          
          # grab arrival time of samplei at latQj
          arrTimeij <- arrivalTimes[i, grepl(pattern=latQjName, 
                                             x=names(arrivalTimes))]
          
          # interference probability is the total probability between tBeg and 
          # tEnd of arrival of sample i
          interfProb[i, j] <- activityProb(data=tActivityDiffuse, 
                                           integrLimits=arrTimeij)
      }
    }
  }
  
  # write output table
  if(!is.na(outTable)){
    interfProb <- as.data.frame(interfProb)
    interfProb$sampleID <- rownames(interfProb)
    write.table(interfProb, file=outTable, 
                quote=FALSE, 
                sep=";",
                col.names=TRUE, 
                row.names=FALSE)
  }
}
