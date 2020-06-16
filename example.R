source('......../probInterf.R')
probInterf(latQpath='...../latQ.txt',
           latQxPath='...../latQx.txt',
           reachInfoPath='...../reachInfo.txt',
           latQerror=c(1.3, 0.32),
           sampleTimesPath='....../sampleTimes.txt',
           outTable='...../probInterf.txt',
           dateTimeFormat='%m/%d/%Y %H:%M')

source('........./lagrComp.R')
lagrangianComparison(reachInfoPath='..../reachInfo.txt',
                     sampleTimesPath='...../sampleTimes.txt',
                     sampleConcentrationsPath='....../sampleConcentrations.txt',
                     pInterfPath='...../probInterf.txt',
                     outTable='...../kCond_1-2.txt',
                     dateTimeFormat='%m/%d/%Y %H:%M',
                     stationUp=1,
                     stationDown=2)
