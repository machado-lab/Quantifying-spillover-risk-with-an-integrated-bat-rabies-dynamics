#### clear Console = ctrl + L

#library(base)
library(rstudioapi)
library(tidyverse)

# Set file location as working directory
setwd(dirname(getActiveDocumentContext()$path)) # => library(rstudioapi)

# Erase workspace
rm(list=ls())

########################################################################
TimeUnit <- "days"
BeginTime <- as.difftime(format(Sys.time(), '%H:%M:%S'), units=TimeUnit);

########################################################################
# To identify for which parameter settings we obtained the results #TODO: change before new execution
IDexec <- "E12"

################################### #TODO: function were changed => change according to it this code
##### LOAD FUNCTIONS
source("FunctionsFit.R")

###################################
##### LOAD DATA
# ### Load sample data
# load("samplecattledata_elevation.RData")
# farms <- data.frame(samplecattledata_elevation)
# load("samplebatdata_SDMelevation.RData")
# roosts<- data.frame(samplebatdata_SDMelevation)

### Load full data
load("allbatdata_SP.RData")
roosts<- data.frame(allbatdata_SP)
load("allcattledata.RData")
farms <- data.frame(allcattledata)

###################################
##### PREPARE DATA FOR SIMS
roosts<- roosts%>%
  filter(abrigo!="Empty") # => library(tidyverse)
roosts<- roosts%>%
  filter(abrigo!="Overnight") # We do not expect transmission from or to Overnight roosts
farms <- subset(farms,farms$SALDO>0) # since we do not include movement of Cattle, we remove from data Farms with no Cattle
# roosts<- roosts%>%
#   filter(UF=="SP") ## done for 'allbatdata' and saved as 'allbatdata_SP'

### Create data frames
dtR <- data.frame(ID=roosts$ID,x=roosts$LON_ORIGIN,y=roosts$LAT_ORIGIN,sdm=roosts$sdm,elev=roosts$elevation, group_type=roosts$abrigo)
dtR <- dtR %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates
dtF <- data.frame(ID=farms$ID,x=farms$LON_ORIGIN,y=farms$LAT_ORIGIN,elev=farms$elevation,totalAnimals = farms$SALDO)
dtF <- dtF %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates
rm(farms, roosts)

### Remove data with missing values
## roosts
missR <- which(dtR$sdm>1|dtR$sdm<0|dtR$elev<0)
if(length(missR)>0){
  dtR <- dtR[-missR,]
}
## farms
missF <- which(dtF$elev<0)
if(length(missF)>0){
  dtF <- dtF[-missF,]
}
rm(missR, missF)

### Approximation of the ranges for lon/lat coordinates
coordRanges <- data.frame(lonL = -53.108752, lonU = -44.161761, latL = -25.303938, latU = -19.781184)
toler <- 6

### Correct typos in coordinates
correctedTmp <- correctCoords(dtR, dtF, coordRanges, toler)
dtR$x <- correctedTmp[[1]]
dtR$y <- correctedTmp[[2]]
dtF$x <- correctedTmp[[3]]
dtF$y <- correctedTmp[[4]]
rm(correctedTmp)

### Remove data with typos in coordinates that we are not able to correct
dtF <- dtF[which(dtF$y > coordRanges$latL - toler & dtF$y < coordRanges$latU + toler &
                   dtF$x > coordRanges$lonL - toler & dtF$x < coordRanges$lonU + toler),]
dtR <- dtR[which(dtR$y > coordRanges$latL - toler & dtR$y < coordRanges$latU + toler &
                       dtR$x > coordRanges$lonL - toler & dtR$x < coordRanges$lonU + toler),]

### Cut data:
dtR <- dtR[which((summary(dtF$x)[2])<dtR$x & dtR$x<(summary(dtF$x)[5]) & (summary(dtF$y)[2])<dtR$y & dtR$y<(summary(dtF$y)[5])),]
dtF <- dtF[which((summary(dtF$x)[2])<dtF$x & dtF$x<(summary(dtF$x)[5]) & (summary(dtF$y)[2])<dtF$y & dtF$y<(summary(dtF$y)[5])),]

dtR <- dtR[which((summary(dtF$x)[2])<dtR$x & dtR$x<(summary(dtF$x)[5]) & (summary(dtF$y)[2])<dtR$y & dtR$y<(summary(dtF$y)[5])),]
dtF <- dtF[which((summary(dtF$x)[2])<dtF$x & dtF$x<(summary(dtF$x)[5]) & (summary(dtF$y)[2])<dtF$y & dtF$y<(summary(dtF$y)[5])),]

## Amounts of bats in roosts
bachTotal <-20
harTotal <- 100
dtR$totalAnimals <- ifelse(dtR$group_type == 'Bachelor', bachTotal,
                        ifelse(dtR$group_type == 'Harem', harTotal,
                               NaN))
if (length( which(is.nan(dtR$totalAnimals)))>0) {
  nameRep<-paste(IDexec, "_Report_NumbBats.txt",sep="")
  sink(nameRep, append=TRUE) #open sink file and add output
  cat("There are NaN in dtR$totalAnimals:", which(is.nan(dtR$totalAnimals)))
}
rm(bachTotal, harTotal)

summR <- summary(dtR)
nR <- nrow(dtR)
summF <- summary(dtF)
nF <- nrow(dtF)

###################################
## BUILD NETWORKS
# transmission between roosts
roostsM <- buildRRnetwork(dtR)
# transmission from roosts to farms
RoostFarmM <- buildRFnetwork(dtR,dtF)

###################################
## PARAMETER SETTINGS FOR ALL SIMS
# maximal time to run the code
MaxTime <- 2 # in TimeUnit (set with BeginTime)
# number of runs for Fitting Phase 1, each for random First Infected Roost
n1Runs <- 2 #TODO: change
# number of different betaRR priors generated
n1Priors <- 20 #TODO: 30000
# number of runs for Fitting Phase 2, each for random posterior betaRR and equilibrium of infected roosts
n2Runs <- 5 #TODO: change
# number of different betaRF priors generated
n2Priors <- 300 #TODO: 30000
#number of simulations post fitting (for each posteriors) // TODO: or numb of simulation each for random posteriors?
nrepeats_postFitting <- 20 #TODO: 30000
# BETAs RANGES TO GENERATE PRIORs:
interval_priorsRR_L <- 1
interval_priorsRR_U <- 600
interval_priorsRF_L <- 1e-10
interval_priorsRF_U <- 450
# we do not need betaFR, betaFF, since there is no transmission from farms (transmission only via bat bites)
#
# recovery rate of a roost
gamR <- 1/(365/2) # half a year
#
settings <- data.frame(MaxTime, TimeUnit, n1Runs, n1Priors, n2Runs, n2Priors, interval_priorsRR_L, interval_priorsRR_U, interval_priorsRF_L, interval_priorsRF_U, gamR, nrepeats_postFitting)
rm(MaxTime, TimeUnit,n1Runs, n1Priors, n2Runs, n2Priors, interval_priorsRR_L, interval_priorsRR_U, interval_priorsRF_L, interval_priorsRF_U, gamR, nrepeats_postFitting)
# prevalence R=0.01 (1% sampled bats infected),
# incidence F~6/(1000000/mean(dtF$totalAnimals)) (6 outbreaks on 1 million of cattle head)
# since we start with no infected farm, and not suppose recovery => incidence F=prevalence F
# we will fit the number of infected roosts and farms after a year:
# we need it name as 'niR', 'niF', respectively for predict function to know the fitted value of the same named variables
fitData <- data.frame(niR=0.01*nrow(dtR), niF=6/(1000000/mean(dtF$totalAnimals))*nrow(dtF)) 
######################################################################################################################################
save(IDexec, BeginTime, summR, nR, summF, nF, roostsM, RoostFarmM, settings, fitData, file=paste(IDexec, "_DataSettings.RData", sep="")) #
######################################################################################################################################

# GENERATE BETA RR PRIORS
betasRXprior <- data.frame(betaRR=rep(0, settings$n1Priors))
betasRXprior$betaRR <- runif(settings$n1Priors, settings$interval_priorsRR_L, settings$interval_priorsRR_U)

# # REPORT
# nameRep<-paste(IDexec, "_Report_RRbetaRange.txt",sep="")
# sink(nameRep, append=TRUE) #open sink file and add output
#   #
#   cat(paste("range to generate betaRR priors: (", settings$interval_priorsRR_L, "; ", settings$interval_priorsRR_U, ") \n", sep=""))
#   cat("min, 1st Qu., Median, Mean, 3rd Qu., Max \n")
#   cat("betaRR priors: \n")
#   cat(summary(betasRXprior$betaRR), "\n")
#   #
# closeAllConnections()

# TO SAVE RESULTS
SimData <- data.frame(betaRR=c(), niR=c(), FirstInfRoost=c(), statusR=c())
save(SimData, file = paste(IDexec, "_RRtransmResults.RData", sep=""))

#### FITTING PHASE 1 ####
#RUN 'n1Runs' SIMULATIONS OF 10 YEAR RR TRANSMISSION FROM RANDOMLY SELECTED FIRST INF ROOST (FOR EACH PRIOR)
for (simRun in 1:settings$n1Runs) {
  ## Set 1st infected roost (=FirstInfRoost)
  FirstInfRoost <- sample(1:nrow(dtR), 1) # random first infected roost (which row)
  # print(c('First infected roost: ', FirstInfRoost))
  ## Fit the model
  simRRtransm(BeginTime, dtR, roostsM, betasRXprior, settings, fitData)
}

# LOAD SIMULATED DATA
SimData <- load(paste(IDexec, "_RRtransmResults.RData", sep=""), RRdata <- new.env())

# FIT THE RR TRANSMISSION RATE
betasRXpost <- fitRR(RRdata$SimData, fitData)

#TODO: before phase 2 add betaRF# ...$betaRF=rep(0, settings_pars$n2Priors)
### END ###

# # REPORT
# nameRep<-paste(IDexec, "_Report_RRbetaRange.txt",sep="")
# sink(nameRep, append=TRUE) #open sink file and add output
# #smt came up, I have to cancel the meeting, sorry for so late to giving you know. Pleas, reschedule the meeting for tomorrow afternoon.
# cat("\n__________________________________________________________\n__________________________________________________________\n")
# #
# closeAllConnections()
