source('c:/gitRepositories/riverPollution/probInterf.R')
probInterf(latQpath='c:/gitRepositories/riverPollution/latQ.txt',
           latQxPath='c:/gitRepositories/riverPollution/latQx.txt',
           reachInfoPath='c:/gitRepositories/riverPollution/reachInfo.txt',
           latQerror=c(1.3, 0.32),
           sampleTimesPath='c:/gitRepositories/riverPollution/sampleTimes.txt',
           outTable='c:/gitRepositories/riverPollution/probInterf.txt',
           dateTimeFormat='%m/%d/%Y %H:%M')


