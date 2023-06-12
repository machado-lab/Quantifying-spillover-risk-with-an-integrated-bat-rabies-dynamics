# creates 1 file including 1 or more requiring files with results
# general (RR & RF)

library(tidyverse) # ggplot()
require(data.table) # setDT()
library(gridExtra) # gred.arrange()

### LOAD FUNCTIONS
source("FncsPlotReport.R")

### LOAD RESULTS & create a csv with all results in one
if(file.exists(paste0(IDexecAll, ".", IDexecPart, "to",
                           IDexecPart+nfiles-1, fileName, "resultsALL.csv"))){
  file.remove(paste0(IDexecAll, ".", IDexecPart, "to",
                           IDexecPart+nfiles-1, fileName, "resultsALL.csv"))
}
sims <- fread(paste0(IDexecAll, ".", IDexecPart, "_", fileName, "results.csv"))
col_names <- colnames(sims)
fwrite(sims, file = paste0(IDexecAll, ".", IDexecPart, "to",
                           IDexecPart+nfiles-1, "_", fileName, "resultsALL.csv"),
       append=T)
if(nfiles > 1){
  for(i in (IDexecPart+1):(IDexecPart+nfiles-1)){
    addsims = fread(paste0(IDexecAll, ".", i, "_", fileName, "results.csv"))
    fwrite(addsims, file = paste0(IDexecAll, ".", IDexecPart, "to",
                                  IDexecPart+nfiles-1, "_", fileName, "resultsALL.csv"),
           append=T)
    sims <- Map(c, sims, addsims)
  }
} else {
  sims <- Map(c, sims)
}
sims <- as.data.frame(matrix(unlist(sims), nrow=length(unlist(sims[1]))))
colnames(sims) <- headerSims(settings$fitted)
################################################################################
