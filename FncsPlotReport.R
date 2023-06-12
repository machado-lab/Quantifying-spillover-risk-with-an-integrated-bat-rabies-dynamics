# functions for loadsims.R (ready for RR & RF)

# PLOT THE SPECIFIED PROPORTIONS OF BEST FITTED DATA
####################################################
plot_fillTable_report <- function(propsToKeep_Vec, sims_res, n_cols, plot_bool,
                                     fillTable_bool, bestFitted_tab=NULL,
                                     range_accept=NULL){
  # ##############################
  # propsToKeep_Vec <- proportionsToKeep
  # sims_res <- sims
  # n_cols <- ncols
  # plot_bool <- T # or F
  # fillTable_bool <- T # or F
  # #
  # bestFitted_tab <- bestFitted
  # range_accept <- range
  # rm(proportionsToKeep, sims, ncols, bestFitted, range)
  ###############################
  propsToKeep_Vec <- sort(unique(propsToKeep_Vec))
  # TODO: change it to vectors and only always prolong the vectors (ie will not
  #need unlist when calling it)
  Ls <- list()
  Us <- list()
  linenames <- list()
  #### z-score #TODO: generalise
  if(settings$fitted==0){
    sims_res$rho <- sqrt(((sims_res$niR-fitData$niR)/sd(sims_res$niR))^2)
  } else if(settings$fitted==1){
    sims_res$rho <- sqrt(((sims_res$outbreaksF-fitData$outbreaksF)/sd(sims_res$outbreaksF))^2)
  }
  for(p in propsToKeep_Vec){
    error <- quantile(sims_res$rho,probs=p)
    # if sd(sims_res$niR) != 0, i.e. if simulated niR are not all the same
    if(error != Inf){
      #TODO: BELOW BACK IF we want to Keep sims with rho less than error (when
      #less or equal to error, we don't need this)
      # if (error==min(sims_res$rho)) {
      #   SimKeptData <- subset(sims_res,sims_res$rho==error)
      # } else {
      #TODO: like this or less than error?
      SimKeptData <- subset(sims_res,sims_res$rho<=error)
      pPos <- which(propsToKeep_Vec==p)
      # known parameter 'fittingPar' is used to fit the unknown parameter
      #'fittedPar'
      if(settings$fitted==0){
        fittingPar <- "niR"
        fittedPar <- "betaRR"
        fileName <- "RR"
      } else if(settings$fitted==1){
        fittingPar <- "outbreaksF"
        fittedPar <- "betaRF"
        fileName <- "RF"
      } else {
        fileName <- "FINAL"
      }
      ResToSort <- SimKeptData[which(names(SimKeptData)==fittingPar)]
      ValsToSort <- SimKeptData[which(names(SimKeptData)==fittedPar)]
      Ls[pPos] <- min(ResToSort)
      Us[pPos] <- max(ResToSort)
      linenames[pPos] <- paste0(p*100, "% interval")
      if(fillTable_bool){
        # fill in the table
        for(c in 1:n_cols){
          bestFitted_tab[pPos,c] <- sorttointerval(ValsToSort, i_int=c,
                                              n_cols
          )
        }
      }
      rm(pPos)
    }
  }
  if(fillTable_bool){
    # add column with rowsums to the table:
    bestFitted_tab <- cbind(bestFitted_tab, rowSums(bestFitted_tab))
    colnames(bestFitted_tab)[ncol(bestFitted_tab)] <- "Sum"
    print(bestFitted_tab)
    ### Start writing to an output file
    sink(paste0(IDexecAll, "_ResReport", length(sims_res[[1]]), "sims.txt"))
    # report the results tabled
    print(bestFitted_tab)
    cat("\n")
    # report if it is or is not in range
    for(p in propsToKeep_Vec[propsToKeep_Vec != 1]){
      pPos <- which(propsToKeep_Vec==p)
      reportInRange(unlist(Ls[pPos]), unlist(Us[pPos]), p*100, range_accept)
    }
    ### Stop writing to the file
    sink()
  }
  if(plot_bool){
    # plot
    pngName <- paste0(IDexecAll, "_plot", fileName, "transm",
                      ifelse(max(propsToKeep_Vec) < 1, "PosteriorsNotCorrected",
                             "All"), length(sims_res[[1]]), "sims.png"
                      )
    # there needs to be by = .(get(substituion(A)), get)substitution(B)), or
    #.(A,B), or c("A", "B") - but with get() it is not possible to have there ""
    #therefore the dot allows us to make it vector of strings this way
    setDT(SimKeptData)[, freq := .N, by = .(get(fittedPar), get(fittingPar))]
    png(pngName)
    plotRRtransm_allPost<- ggplot(SimKeptData, aes(betaRR, niR, col=freq)) +
      # fitted value
      geom_point() + 
      scale_color_gradient(low="yellow", high = "red") +
      geom_hline(aes(yintercept=fitData$niR, linetype='fitted value'),
                 colour='red')
    for (p in propsToKeep_Vec[propsToKeep_Vec != 1]) {
      pPos <- which(propsToKeep_Vec==p)
      plotRRtransm_allPost <- plotRRtransm_allPost +
        # Low and Upper limits of the kept data
        geom_hline(aes(yintercept=!!unlist(Ls[pPos]),
                   linetype=!!unlist(linenames[pPos])), colour = linecols[pPos]
        ) +
        geom_hline(aes(yintercept=!!unlist(Us[pPos]),
                   linetype=!!unlist(linenames[pPos])), colour = linecols[pPos]
        )
    }
    nprops <- length(propsToKeep_Vec[propsToKeep_Vec != 1])
    plotRRtransm_allPost <- plotRRtransm_allPost +
      scale_linetype_manual(name = "kept data",
                            values = c("dotted", rep("dashed",nprops-1),
                                       "solid"),
                            breaks = c("fitted value", linenames[1:nprops]),
                            guide = guide_legend(
                              override.aes = list(color =c("red",
                                                           linecols[1:nprops]
                                                           )
                                                  )
                              )
                            )
      # theme(
      #   # center, thicken and enlarge the title
      #   plot.title = element_text(hjust = 0.5, face = "bold", size = (30)),
      #   # setup the legend
      #   legend.position = c(0.82, 0.8), legend.key.height = unit(1.5, "cm"),
      #   legend.key.width = unit(1, "cm"),
      #   legend.title = element_text(face = "bold", size = (20)),
      #   legend.text = element_text(size = (15)), 
      #   legend.background = element_rect(fill="lightblue"),
      #   # axis title size
      #   axis.title=element_text(size=23)
      # )
    print(plotRRtransm_allPost)
    dev.off()
  }
}
####################################################

# SORT TO INTERVAL
##################
sorttointerval <- function(Vals_ToSort, i_int, n_cols, L_limit=settings$beta_L){
                           #TODO: add , intLength = 1){
  # Vals_ToSort <- ValsToSort
  # i_int <- the order of interval
  # last <- F #(otherwise) T if it is the last interval
  # rowtofill <- bestFitted_tab[pPos,]
  if(i_int!=n_cols){
    return(length(which(
      (Vals_ToSort >= floor(L_limit)+i_int-1 &
        (Vals_ToSort < floor(L_limit)+(i_int)))
    )))
  } else {
    return(length(which(
      (Vals_ToSort >= floor(L_limit)+i_int-1) &
        (Vals_ToSort <= floor(L_limit)+(i_int))
    )))
  }
}
##################

# WRITE DOWN IF IT IS OR IS NOT IN RANGE
########################################
reportInRange <- function(L, U, percentkept, range){
  inrange <- L >= (fitData$niR - range) & U <= (fitData$niR + range)
  cat(paste0("When keeping ", percentkept,
             "%,\n the deviation from fitted value is",
             ifelse(inrange, "", " not"), " in +/-", range, " range.\n",
             ifelse(inrange, "", paste0(" It is in range: ",
                                        ifelse(U > fitData$niR, "+", ""),
                                        U - fitData$niR,
                                        ifelse(L > fitData$niR, "/+", "/"),
                                        L - fitData, ".\n"))))
}
########################################
