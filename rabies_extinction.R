# rabies_extinction.R
#------------------------------------------------------------------------------#
# How often does rabies go extinct from the metapopulation in the case of      #
# particular interventon strategy (eg. in the absence of interventon, ie       #
# interv = "No intervention")?                                                  #
#                                                                              #
# How often there is no spillover, how often with no detected spillover to     #
# farms there are some undetected infections in farms and/or roosts, how often #
# the only infection is the initial one?                                                         #
#------------------------------------------------------------------------------#

# Clear the console in RStudio
cat("\f")
# Erase workspace
rm(list=ls())
# Set file location as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
## TO SET
SdmInit = "High"
nSimsPerPost <- 50
# For which intervention we want to 
interv = "No intervention"
# Where you have saved the results data
folder_withdata = paste0("../sims/sims_", nSimsPerPost)
################################################################################
if(nSimsPerPost==10){
  WhichSims <- "1.1to2_2.1to2"
} else if(nSimsPerPost==20){
  WhichSims <- "1.1to3_2.1to3"
} else if(nSimsPerPost==30){
  WhichSims <- "1.1to4_2.1to4"
} else if(nSimsPerPost==50){
  WhichSims <- "1.1to5_2.1to5"
} else {
  stop('Invalid number of simulations per posterior (nSimsPerPost).')
}
load(paste0(folder_withdata, "/AllsimsAllIntervFinal", SdmInit,
            "_", WhichSims, ".RData"))

R1_i = which(colnames(ResultsAllAtoD)=="R1")
F1_i = which(colnames(ResultsAllAtoD)=="F1")
Flast_i = ncol(ResultsAllAtoD)-1
roosts_indexes = R1_i:(F1_i-1)
farms_indexes = F1_i:(length(ResultsAllAtoD)-1)

Res_interv = subset(ResultsAllAtoD, ResultsAllAtoD$intervention == interv)
n_sims = nrow(Res_interv)

### EASY ONE
# SHOW
# 1) How often there is no spillover to farms (which can include sims where
# rabies is present but only in roosts)
Res_interv_nospillover = subset(Res_interv, Res_interv$outbreaksF == 0)
n_NoSpillover = nrow(Res_interv_nospillover)

# 2) How often there are some roosts infected and/or farms exposed but no
# detected outbreak in farm yet?
# calc. numb of undetected infections in farms and/or roosts
Statuses_nointerv_nospillover = Res_interv_nospillover[,R1_i:Flast_i]
Res_interv_nospillover$n_infections =
  rowSums(Statuses_nointerv_nospillover > 0)
#
n_SomeInfections = length(which(Res_interv_nospillover$n_infections>0))

# 3) How often only the FirstInfRoost is still infected (which means there is no
# other infection since the induction of rabies virus, but this one is still
# infected after one year simulated).
# save status of FirstInfRoost
Res_interv_nospillover$status_FirstInfRoost =
  rep(-2,nrow(Res_interv_nospillover))
for(row_res in 1:nrow(Statuses_nointerv_nospillover)){
  Res_interv_nospillover$status_FirstInfRoost[row_res] =
    Statuses_nointerv_nospillover[row_res,
                                  Res_interv_nospillover$FirstInfRoost[
                                    row_res
                                  ]
    ]
}
#
n_FirstInfOnly = length(which(
  Res_interv_nospillover$n_infections == 1 &
    Res_interv_nospillover$status_FirstInfRoost == 1
))

# 4) How often there is full extinction, ie no farm and no roost infected after
# one year simulated.        
n_Extinction = n_NoSpillover - n_SomeInfections
