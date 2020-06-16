source('c:/gitRepositories/riverPollution/probInterf.R')
probInterf(latQpath='c:/gitRepositories/riverPollution/latQ.txt',
           latQxPath='c:/gitRepositories/riverPollution/latQx.txt',
           reachInfoPath='c:/gitRepositories/riverPollution/reachInfo.txt',
           latQerror=c(1.3, 0.32),
           sampleTimesPath='c:/gitRepositories/riverPollution/sampleTimes.txt',
           outTable='c:/gitRepositories/riverPollution/probInterf.txt',
           dateTimeFormat='%m/%d/%Y %H:%M')

source('c:/gitRepositories/riverPollution/lagrComp.R')
lagrangianComparison(reachInfoPath='c:/gitRepositories/riverPollution/reachInfo.txt',
                     sampleTimesPath='c:/gitRepositories/riverPollution/sampleTimes.txt',
                     sampleConcentrationsPath='c:/gitRepositories/riverPollution/sampleConcentrations.txt',
                     pInterfPath='c:/gitRepositories/riverPollution/probInterf.txt',
                     outTable='c:/gitRepositories/riverPollution/kCond_1-2.txt',
                     dateTimeFormat='%m/%d/%Y %H:%M',
                     stationUp=1,
                     stationDown=2)
