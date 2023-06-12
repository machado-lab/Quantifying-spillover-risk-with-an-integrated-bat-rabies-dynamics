# load sims with results of 1st or 2nd phase of fitting
# keep sims with fitting parameter in +/- range from fitting value
# do ABC fitting
# return posteriors (1st phase: betaRR and R statuses; 2nd phase: betaRF)

library(rstudioapi) # getActiveDocumentContext()
library(data.table) # fread()
# Set file location as working directory
setwd(dirname(getActiveDocumentContext()$path))
# Erase workspace
rm(list=ls())

# TO SET
################################################################################
# ID
IDexecAll <- "BD2a"
# For how many files do you want to process?
nfiles <- 2
# What is the order of the first one?
IDexecPart <- 1
# How many priors do you have?
#G5.1=1000; G5.2=780; G5.3=608; G5.4=112 // 1to4=2500
#testG5.1=10; testG5.2=171; testG5.3=12
#BD2a.1= 5*29; BD2a.2=30*29
nsims <- 35*29
# What distance from the fitting value do you want to accept to keep it as
#posteriors?
range <- 5
################################################################################

### LOAD FUNCTIONS
source("FncsABC.R")

### # LOAD SETTINGS of simulations (fitData, ...)
load(paste0(IDexecAll, "_DataSettings.RData"))

#####
# CREATE ONE FILE FOR ALL RESULTS IF DOES NOT EXIST YET 
#####
if(settings$fitted==0){
  fileName <- "RR"
} else if(settings$fitted==1){
  fileName <- "RF"
} else if(settings$fitted==2){
  fileName <- "FINAL"
} else {
  print("ERROR: The acceptable values of 'fitted' are 0, 1, or 2.")
  stop()
}

if(!file.exists(paste0(IDexecAll, ".", IDexecPart, "to",
                       IDexecPart+nfiles-1, "_", fileName, "resultsALL.csv"))){
  source("loadsims.R")
}

### CHANGE FITTING VALUES IF THIS IS TEST USING SUBSET OF DATA
if(nR < 4167){
  fitData <- data.frame(niR=41.7, outbreaksF=226)
}
### LOAD SIMULATED DATA
sims <- fread(paste0(IDexecAll, ".", IDexecPart, "to",
                     IDexecPart+nfiles-1, "_", fileName, "resultsALL.csv"))
# check if the data really include all the simulations you require
if(nrow(sims)!=nsims){
  print(paste("ERROR: You required to get posteriors of", nsims,
  "simulations, but the loaded data include results of", nrow(sims),
  "simulations."))
  stop()
} else {
  # for 1st phase of fitting: add niR
  if(settings$fitted==0){
    sims <- results(sims)
  }
  ########
  # GET CORRECTED POSTERIORS
  ########
  posts <- fit(sims, fitData, range, settings$fitted)
  # for 2nd phase of fitting: posts is a vector, but I need data.frame
  if(settings$fitted==1){
    posts <- data.frame(betaRF = posts)
  }
  fwrite(posts, file = paste0(IDexecAll, ".", IDexecPart, "to",
                                IDexecPart+nfiles-1, "_beta", fileName,
                                "posteriors_range", range, ".csv"),
         sep = ";", row.names = FALSE, col.names = TRUE
  )
  if(settings$fitted==0){
    postsRR <- posts
    save(postsRR, file = paste0(IDexecAll, ".", IDexecPart, "to",
                               IDexecPart+5nfiles-1, "_beta", fileName,
                               "posteriors_range", range, "_v2.RData"),
         version = 2)
  } else { #TODO: do it general, but same time have different name
                # postsRR X postsRF
    postsRF <- posts 
    save(postsRF, file = paste0(IDexecAll, ".", IDexecPart, "to",
                                IDexecPart+nfiles-1, "_beta", fileName,
                                "posteriors_range", range, "_v2.RData"),
         version = 2)
  }
}

#TODO:
# PLOT POSTERIORS

# TABLE POSTERIORS
  
