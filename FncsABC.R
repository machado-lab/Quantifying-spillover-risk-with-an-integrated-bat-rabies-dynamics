# functions for ABC_CorrectBetas.R (and related files)

# HEADER OF SIMS
################
headerSims <- function(fitted_n) {
  # fitted_n <- fitted
  # rm(fitted)
  if(fitted_n==0) {
    X <- c("R")
    nX <- c(nR)
    inits <- c("FirstInfRoost", "betaRR")
  } else {
    X=c("R", "F")
    nX <- c(nR, nF)
    inits <- c("betaRR", "betaRF")
  }
  StatusesList <- list()
  j <- 0
  for(pop in 1:length(X)){
    for(i in 1:(nX[pop])){
      j <- j + 1
      name <- paste0(X[pop], i)
      StatusesList[j] <- name
    }
  }
  if(fitted_n==0){
    names_res <- c(inits, StatusesList)
  } else{
    names_res <- c(inits, "outbreaksF", StatusesList)
  }
  return(names_res)
}
################

# STATUSES OF ROOSTS TO NUMBS OF INFECTIOUS (niR)
#################################################
results <- function(sims_res){
  sims_res$niR <- rep(0, nrow(sims_res))
  for(i in 1:nrow(sims_res)){
    sims_res$niR[i] <- length(which(sims_res[i,3:(nR+2)]==1))
  }
  return(sims_res)
}
#################################################
# # TODO:
# results <- function(sims_res, n_setVar){
#   sims_res$niR <- rep(0, nrow(sims_res))
#   for(i in 1:nrow(sims_res)){
#     sims_res$niR[i] <- length(which(sims_res[i,(n_setVar+1):(nR+n_setVar)]==1))
#   }
#   # return(sims_res)
#   return(sims_res$niR)
# }
# # MORE DETAILS:
# results <- function(sims_matrix, n_setVar){
#   sims_res <- sims[,1:n_setVar]
#   sims_res$niR <- rep(0, nrow(sims_matrix))
#   sims_res$npastiR <- rep(0, nrow(sims_matrix))
#   sims_res$nneveriR <- rep(0, nrow(sims_matrix))
#   sims_res$niF <- rep(0, nrow(sims_matrix))
#   sims_res$neF <- rep(0, nrow(sims_matrix))
#   sims_res$npastiF <- rep(0, nrow(sims_matrix))
#   sims_res$nneveriF <- rep(0, nrow(sims_matrix))
#   for(i in 1:nrow(sims_res)){
#     # Roosts
#     sims_res$niR[i] <- length(which(sims_matrix[i,(n_setVar+1):(nR+n_setVar)]==1))
#     sims_res$npastiR[i] <- length(which(sims_matrix[i,(n_setVar+1):(nR+n_setVar)]==-1))
#     sims_res$nneveriR[i] <- length(which(sims_matrix[i,(n_setVar+1):(nR+n_setVar)]==0))
#     # Farms
#     sims_res$niF[i] <- length(which(sims_matrix[i,(nR+n_setVar+1):(nR+nF+n_setVar)]==2))
#     sims_res$neF[i] <- length(which(sims_matrix[i,(nR+n_setVar+1):(nR+nF+n_setVar)]==1))
#     sims_res$npastiF[i] <- length(which(sims_matrix[i,(nR+n_setVar+1):(nR+nF+n_setVar)]==-1))
#     sims_res$nneveriF[i] <- length(which(sims_matrix[i,(nR+n_setVar+1):(nR+nF+n_setVar)]==0))
#   }
#   # return(sims_res)
#   return(sims_res)
# }


# ABC FITTING
#############
fit <- function(Sim_Data, fit_Data, range_fit, n_fitted) {
  #######################
  # Sim_Data <- sims
  # fit_Data <- fitData
  # range_fit <- range
  # n_fitted <- fitted
  #
  # rm(sims, fitData, range, fitted)
  ###########################################################
  ## fit the fitted parameter so that the fitting value at the end of sim is as in
  #fit_Data
  #
  if(n_fitted==0){
    fittingPar <- "niR"
  } else if(n_fitted==1){
    fittingPar <- "outbreaksF"
  }
  # distance to fitting value = abs(Z--score) /using empirical sd/:
  # double brackets: calling a column by substitutional variable; ie evaluate
  #what is inside; ie use column saved in substitutional variable
  Sim_Data$rho <- sqrt(((Sim_Data[[fittingPar]]-fit_Data[[fittingPar]])
                        / sd(Sim_Data[[fittingPar]]))^2)
  # # We keep the best fitting ones so that all are in required 'range' +/- from
  #the fitting value
  SimKeptData <- subset(Sim_Data,
                        abs(Sim_Data[[fittingPar]]-fit_Data[[fittingPar]])
                        <= range_fit)
  error <- max(SimKeptData$rho)
  ## This regression bit does not change results
  #(as the fitted.values=predicted)
  # Apply the weighted linear regression for the fitted parameter
  #E(theta|S(x))
  # Weight the values kept:
  weight <- (1/error)*(1-(SimKeptData$rho/error)^2)
  if(n_fitted==0){
    lrbeta <- lm(betaRR~niR,SimKeptData,weights=weight)
    beta <- "betaRR"
  } else {
    lrbeta <- lm(betaRF~outbreaksF,SimKeptData,weights=weight)
    beta <- "betaRF"
  }
  rm(error, weight)
  ## Calculate E(theta|S(y))
  predbeta <- predict(lrbeta,fit_Data)
  ## Correct the betas
  cbeta <- SimKeptData[[beta]]-lrbeta$fitted.values+predbeta
  ## Posteriors
  post <- SimKeptData
  post[[beta]] <- cbeta
  # keep only positive betas
  post <- subset(post, post[[beta]]>0)
  ## RETURN beta & if fitted==0 statuses of R
  if(n_fitted==0){
    return(post[,2:(length(post[1])-2)])
  } else if(n_fitted==1){
    return(post[[beta]])
  }
}
#############

#-------------------------------------------------------------------------------
# 
# # PLOT
# ######
# plot <- function(propsToKeep_Vec, sims_res, n_cols, plot_bool,
#                                   fillTable_bool, bestFitted_tab=NULL,
#                                   range_R=NULL){
#   ###############################
#   # propsToKeep_Vec <- proportionsToKeep
#   # sims_res <- sims
#   # n_cols <- ncols
#   # fillTable_bool <- T F
#   ##
#   # bestFitted_tab <- bestFitted
#   # range_R <- rangeR
#   # rm(proportionsToKeep, sims, ncols, bestFitted, rangeR)
#   ###############################
#   propsToKeep_Vec <- sort(unique(propsToKeep_Vec))
#   # TODO: change it to vectors and only always prolong the vectors (ie will not
#   #need unlist when calling it)
#   Ls <- list()
#   Us <- list()
#   linenames <- list()
#   for(p in propsToKeep_Vec){
#     error <- quantile(sims_res$rho,probs=p)
#     # if sd(sims_res[[nI]]) != 0, i.e. if simulated nI are not all the same
#     if(error != Inf){
#       #TODO: BELOW BACK IF we want to Keep sims with rho less than error (when
#       #less or equal to error, we don't need this)
#       # if (error==min(sims_res$rho)) {
#       #   SimKeptData <- subset(sims_res,sims_res$rho==error)
#       # } else {
#       SimKeptData <- subset(sims_res,sims_res$rho<=error)#TODO: or less than error?
#       pPos <- which(propsToKeep_Vec==p)
#       Ls[pPos] <- min(SimKeptData[[nI]])
#       Us[pPos] <- max(SimKeptData[[nI]])
#       linenames[pPos] <- paste0(p*100, "% interval")
#       if(fillTable_bool){
#         # fill in the table
#         for(c in 1:n_cols){
#           bestFitted_tab[pPos,] <- sorttointerval(SimKeptData, i_int=c,
#                                                   last=(c==n_cols),
#                                                   bestFitted_tab[pPos,]
#           )
#         }
#       }
#       rm(pPos)
#     }
#   }
#   if(fillTable_bool){
#     # add column with rowsums to the table:
#     bestFitted_tab <- cbind(bestFitted_tab, rowSums(bestFitted_tab))
#     colnames(bestFitted_tab)[ncol(bestFitted_tab)] <- "Sum"
#     print(bestFitted_tab)
#     ### Start writing to an output file
#     sink(paste0(IDexecAll, "_ResReport", length(sims_res[[1]]), "sims.txt"))
#     # report the results tabled
#     print(bestFitted_tab)
#     cat("\n")
#     # report if it is or is not in range
#     for(p in propsToKeep_Vec[propsToKeep_Vec != 1]){
#       pPos <- which(propsToKeep_Vec==p)
#       reportInRange(unlist(Ls[pPos]), unlist(Us[pPos]), p*100, range_R)
#     }
#     ### Stop writing to the file
#     sink()
#   }
#   if(plot_bool){
#     # plot
#     pngName <- paste0(IDexecAll, "_plotRRtransm",
#                       ifelse(max(propsToKeep_Vec) < 1, "PosteriorsNotCorrected",
#                              "All"), length(sims_res[[1]]), "sims.png"
#     )
#     setDT(SimKeptData)[, freq := .N, by = .(betaRR, nI)]
#     png(pngName)
#     plotRRtransm_allPost<- ggplot(SimKeptData, aes(betaRR, nI, col=freq)) +
#       # fitted value
#       geom_point() + 
#       scale_color_gradient(low="yellow", high = "red") +
#       geom_hline(aes(yintercept=fitData[[nI]], linetype='fitted value'), colour='red') #+
#     for (p in propsToKeep_Vec[propsToKeep_Vec != 1]) {
#       pPos <- which(propsToKeep_Vec==p)
#       plotRRtransm_allPost <- plotRRtransm_allPost +
#         # Low and Upper limits of the kept data
#         geom_hline(aes(yintercept=!!unlist(Ls[pPos]),
#                        linetype=!!unlist(linenames[pPos])), colour = linecols[pPos]
#         ) +
#         geom_hline(aes(yintercept=!!unlist(Us[pPos]),
#                        linetype=!!unlist(linenames[pPos])), colour = linecols[pPos]
#         )
#     }
#     nprops <- length(propsToKeep_Vec[propsToKeep_Vec != 1])
#     plotRRtransm_allPost <- plotRRtransm_allPost +
#       scale_linetype_manual(name = "kept data", values = c("dotted",
#                                                            rep("dashed",nprops-1),
#                                                            "solid"),
#                             breaks = c("fitted value", linenames[1:nprops]),
#                             guide = guide_legend(
#                               override.aes = list(color =c("red",
#                                                            linecols[1:nprops]
#                               )
#                               )
#                             )
#       )
#     # theme(
#     #   # center, thicken and enlarge the title
#     #   plot.title = element_text(hjust = 0.5, face = "bold", size = (30)),
#     #   # setup the legend
#     #   legend.position = c(0.82, 0.8), legend.key.height = unit(1.5, "cm"),
#     #   legend.key.width = unit(1, "cm"),
#     #   legend.title = element_text(face = "bold", size = (20)),
#     #   legend.text = element_text(size = (15)), 
#     #   legend.background = element_rect(fill="lightblue"),
#     #   # axis title size
#     #   axis.title=element_text(size=23)
#     # )
#     print(plotRRtransm_allPost)
#     dev.off()
#   }
# }
# ######
# 
# # TABLE
# #######
# 
# #######
# 
