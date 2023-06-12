
# Set file location as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#TODO: change the order of interventions in legend to V&C, V, C, no

########## TO SET ###
# What is different for the simulations compared in the histograms
#"interventions" or "sdm" values?
Diff <- "interventions"

if(Diff=="interventions"){
  #TODO: create a figure with all 3 images in it as a for loop
  SdmInit <- "Low" # "High", "Middle", or "Low"
} else if(Diff=="sdm"){
  #TODO: create a figure with all 4 images in it as a for loop
  Interventions <- "" # "CullVacc", "Vacc", "Cull", or "NoInterv"
} else {
  stop('Invalid argument for Diff, only "interventions" or "sdm" are valid.')
}
# How many number of simulations per posterior you want to use for histograms?
nSimsPerPost <- 10
# What format to use to save histogram, "pdf" or "png"?
Format <- "png"
#####################
if(Diff=="interventions"){
  WhatResA <- paste0("FinalCullVacc", SdmInit)
  WhatResB <- paste0("FinalVacc", SdmInit)
  WhatResC <- paste0("FinalCull", SdmInit)
  WhatResD <- paste0("FinalNoInterv", SdmInit)
} else { # Diff=="sdm"
  WhatResA <- paste0("Final", Interventions, "High")
  WhatResB <- paste0("Final", Interventions, "Middle")
  WhatResC <- paste0("Final", Interventions, "Low")
}
#
if(nSimsPerPost==10){
  WhichSims <- "1.1to2_2.1to2"
} else if(nSimsPerPost==20){
  WhichSims <- "1.1to3_2.1to3"
} else {
  stop('Invalid number of simulations per posterior (nSimsPerPost).')
}
# load All results (TODO: and settings?)
# A
load(paste0("../sims/sims_", nSimsPerPost, "/Allsims", WhatResA, "_BD",
            WhatResA, WhichSims, ".RData"))
resultsAllA <- resultsAll
# load(paste0(WhatResA, "Results/BD", WhatResA, "1_DataSettings.RData"))
# settingsA <- settings
# B
load(paste0("../sims/sims_", nSimsPerPost, "/Allsims", WhatResB, "_BD",
            WhatResB, WhichSims, ".RData"))
resultsAllB <- resultsAll
# load(paste0(WhatResB, "Results/BD", WhatResB, "1_DataSettings.RData"))
# settingsB <- settings
# C
load(paste0("../sims/sims_", nSimsPerPost, "/Allsims", WhatResC, "_BD",
            WhatResC, WhichSims, ".RData"))
resultsAllC <- resultsAll
# load(paste0(WhatResC, "Results/BD", WhatResC, "1_DataSettings.RData"))
# settingsC <- settings
if(Diff=="interventions"){
  load(paste0("../sims/sims_", nSimsPerPost, "/Allsims", WhatResD, "_BD",
              WhatResD, WhichSims, ".RData"))
  resultsAllD <- resultsAll
  # load(paste0(WhatResD, "Results/BD", WhatResD, "1_DataSettings.RData"))
  # settingsD <- settings
}
rm(resultsAll)
# rm(settings)


if(Diff=="interventions"){
  ########
  resultsAllA$intervention <- "Vaccination & Culling"
  resultsAllB$intervention <- "Vaccination"
  resultsAllC$intervention <- "Culling"
  resultsAllD$intervention <- "No intervention"
  
  ResultsAllAtoD <- rbind(resultsAllA, resultsAllB, resultsAllC, resultsAllD)
  ResultsAllAtoD$intervention <-
    factor(ResultsAllAtoD$intervention,
           levels=c("Vaccination & Culling", "Vaccination", "Culling",
                    "No intervention"))

  # HISTOGRAM
  library(ggplot2) # ggplot()
  if(Format=="pdf"){
    pdf(paste0("hists_sdm", SdmInit, ".pdf"),width=7,height=4.64)
  } else if(Format=="png"){
    png(paste0("hists_sdm", SdmInit, ".png"),width=7,height=4.64, unit="in",
        res=200)
  } else {
    stop('Unacceptable format, set "pdf" or "png"')
  }
  par(font=2, cex = .75, cex.axis=0.5, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,
      mgp=c(3,0.4,0))
  # histograms
  #TODO: max of y automatic (max of all frequencies)
  x_max <- 250
  nIntervals <- 30
  p <- ggplot(ResultsAllAtoD, aes(x=outbreaksF, y=..density.., color=intervention),
         #TODO: fix the following, this obviously doesn't work
         show.legend = (SdmInit=="Low")) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, 0.15)) +
    # scale_x_continuous(breaks=seq(0, x_max, x_max/nIntervals)) +
    # xlim(-0.1, x_max) +
    geom_histogram(fill="white", alpha=0.2, position="identity",
                   bins=(nIntervals+1), boundary = 0) + 
    labs(title=paste(SdmInit, "sdm"), x="# Outbreaks in farms", y = "Density") +
    theme(plot.title = element_text(hjust = 0.5, face="bold"),
          legend.title = element_text(face="bold"))
  p
} else {
  stop("TODO: Comparing hists for different sdm values is not finished yet.")
}

dev.off()

tmp <- ggplot_build(p)
tmp$data[[1]]$xmin
tmp$data[[1]]$xmax
tmp$data[[1]]$density
length(which(ResultsAllAtoD$intervention=="Vaccination & Culling" & ResultsAllAtoD$outbreaksF >=0 & ResultsAllAtoD$outbreaksF <7))/1070

# if(Diff=="interventions"){
#   #TODO TABLE WITH DATA OF HISTOGRAM
#   if(Format=="pdf"){
#     pdf(paste0("HistsTab_sdm", SdmInit, ".pdf"),width=7,height=4.64)
#   } else if(Format=="png"){
#     png(paste0("HistsTab_sdm", SdmInit, ".png"),width=7,height=4.64, unit="in",
#         res=200)
#     library(grid)
#     library(gridExtra)
#     # Function to extract legend
#     # https://stackoverflow.com/a/13650878/496488
#     g_legend <- function(plot_name){
#       tmp <- ggplot_gtable(ggplot_build(plot_name))
#       leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#       legend <- tmp$grobs[[leg]]
#       return(legend)}
#     # Extract the legend as a separate grob
#     leg = g_legend(p)
#     # Create a table grob
#     tab = t(ResultsAllAtoD)
#     tab = tableGrob(tab, rows=NULL)
#     tab$widths <- unit(rep(1/ncol(tab), ncol(tab)), "npc")
#     # Lay out plot, legend, and table grob
#     grid.arrange(arrangeGrob(nullGrob(), 
#                              p + guides(fill=FALSE) + 
#                                theme(axis.text.x=element_blank(),
#                                      axis.title.x=element_blank(),
#                                      axis.ticks.x=element_blank()),
#                              widths=c(1,8)), 
#                  arrangeGrob(arrangeGrob(nullGrob(),leg,heights=c(1,10)),
#                              tab, nullGrob(), widths=c(6,20,1)),
#                  heights=c(4,1))
#     
#   } else {
#     stop('Unacceptable format, set "pdf" or "png"')
#   }
# 
# } else {
#   stop("TODO: Comparing hists for different sdm values is not finished yet.")
# }
# 
# dev.off()

