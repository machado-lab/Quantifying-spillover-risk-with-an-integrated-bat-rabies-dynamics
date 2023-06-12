library(base)
library(rstudioapi)
library(grImport)
library(lattice)
library(tidyverse)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
rm(list=ls())  #Erase workspace

source("FunctionsPrelim.R") #load functions
load("allcattledata_v2.RData")
cattle <- data.frame(allcattledata)
load("allbatdata_SP_v2.RData")
bats<- data.frame(allbatdata_SP)


##For sample data
# dtc <- data.frame(ID=cattle$ID,x=cattle$LON_ORIGINAL,y=cattle$LAT_ORIGINAL,elev=cattle$elevation,totalAnimals = cattle$SALDO)
# dtc <- dtc %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates
# dtb <- data.frame(ID=bats$ID,x=bats$LON_ORIGINAL,y=bats$LAT_ORIGINAL,sdm=bats$sdm,elev=bats$elevation, group_type=bats$abrigo)
# dtb <- dtb %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates

##For full data
dtc <- data.frame(ID=cattle$ID,x=cattle$LON_ORIGIN,y=cattle$LAT_ORIGIN,elev=cattle$elevation,totalAnimals = cattle$SALDO)
dtc <- dtc %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates
dtb <- data.frame(ID=bats$ID,x=bats$LON_ORIGIN,y=bats$LAT_ORIGIN,sdm=bats$sdm,elev=bats$elevation, group_type=bats$abrigo)
dtb <- dtb %>% distinct(x,y, .keep_all = TRUE) # rmv same coordinates

##Cut data:
dtb <- dtb[which((-48.43)<dtb$x & dtb$x<(-47.15) & (-22.57)<dtb$y & dtb$y<(-21.87)),]
dtc <- dtc[which((-48.43)<dtc$x & dtc$x<(-47.15) & (-22.57)<dtc$y & dtc$y<(-21.87)),]

bach_amount <-20
har_amount <- 100
overn_amount <- 0
dtb$numb_bats <- ifelse(dtb$group_type == 'Bachelor', bach_amount,
                        ifelse(dtb$group_type == 'Harem', har_amount,
                               overn_amount))
#Using simulated dataset
set.seed(1)
dt <- generatePops(20)
roostM  <- buildRoadNetwork(dt)
RoostFarmM <- buildRoostNetwork(dt)

# roadM <- buildRoadNetworkR(dtc) ## no transmission between farms
#roostM <- buildRoostNetworkR(dtb)
#RoostFarmM <- buildRoostToFarmNetworkR(dtb,dtc)

#statusC <- rep(0,nrow(dtc)) # Infected status (0 uninf, 1 inf)
#statusB <- rep(0,nrow(dtb)) # Infected status (0 uninf, 1 inf)

statusC <- rep(0,nrow(dt)) # Infected status (0 uninf, 1 inf)
statusB <- rep(0,nrow(dt)) # Infected status (0 uninf, 1 inf)



gamB <- 1/25 #two weeks to recover
gamC <- 0  #We run the code until the first C is infected, so there is no recovery ##before: 1/6 #ie 8 weeks to recover

###################################
## seed 1 roost, run x timesteps, until first farm infected, check how many roosts got infected
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
#cl <- makeCluster(cores[1]-1, type = "FORK") #not to overload your computer
cl <- makeCluster(cores[1]-1)

nrepeats <- 20 #TODO: bigger numb in future (eg 30000)
betaBB <- runif(nrepeats,10^-6,10^-4)
betaBC <- runif(nrepeats,10^-6,10^-4)
betaCB <- rep(0,nrepeats) # no CB transmission
betaCC <- rep(0,nrepeats) # no CC transmission


maxstep <- rep(0,nrepeats)

registerDoParallel(cl)

maxstep <- foreach(i=1:nrepeats, .combine=cbind) %dopar% {
  #print(c("Start repeat: ",i))
  statusC <- rep(0,nrow(dt)) # Infected status (0 uninf, 1 inf)
  statusB <- rep(0,nrow(dt)) # Infected status (0 uninf, 1 inf)

  nfi <- sum(statusC) #number of farm infected
  statusB[4] <- 1 #insert random infection in bats

  step <- 0
  i=1:nrepeats
  while(nfi < 1){
    #print(c("Step: ",step))
    tmp <-stepFunction(statusB,statusC,betaBB[i],betaBC[i],betaCB[i],betaCC[i],gamB,gamC)
    statusB <- tmp[(length(statusC)+1):length(tmp)]
    statusC <- tmp[1:length(statusC)]
    step <- step+1
    maxstep[i] <- step
    nfi <- sum(statusC)
    if (step == 100) break
  }
  return(maxstep[i])
}

stopCluster(cl)

# Number of maxstep (V1) in dependence of transmission rates BB (V2) & BC (V3)
#save.image("Results.RData")
info <- matrix(c(maxstep,betaBB,betaBC), ncol=3)
## Using ggplot
ggplot(as.data.frame(info), aes(V3, V2, col=V1)) +
  geom_point(alpha=0.5) +
  scale_color_gradient(low="red", high="yellow")
#####

########################################
## Do quick ABC fitting
## fit to 1 infected farm per 5 timesteps aprox #TODO: Later, fit to reasonable summary statistics.
Sim <- data.frame(betaBB,betaBC,maxstep=c(maxstep))
data <- data.frame(maxstep=5)
rho<-sqrt(((Sim$maxstep-data$maxstep)/sd(Sim$maxstep))^2) # distance to "target", 5, (using empirical sd) = abs(Z-score)
Sim$rho <- rho

# We keep the lowest rho # -- 10% quantile, ie probs=0.1 #TODO: Later, to decide which quantile (5-10%) we want. Lower quantile, more precise value, but less data keeping.
## i.e. best fitting ones
error <- quantile(Sim$rho,probs=0.1)
SimKept <- subset(Sim,Sim$rho<error)
# Weight the values saved #
weight<-(1/error)*(1-(SimKept$rho/error)^2)

## This regression bit does not change results (as the fitted.values=predicted)
# Apply the weighted linear regression for the 2 fitted parameters E(theta|S(x))
lrbetaBB <- lm(betaBB~maxstep,SimKept,weights=weight)
lrbetaBC <- lm(betaBC~maxstep,SimKept,weights=weight)

# Calculate E(theta|S(y))
predbetaBB<-predict(lrbetaBB,data)
predbetaBC<-predict(lrbetaBC,data)

# Correct the parameters
cbetaBB<-SimKept$betaBB-lrbetaBB$fitted.values+predbetaBB
cbetaBC<-SimKept$betaBC-lrbetaBC$fitted.values+predbetaBC
###########################################################

################# Lets pick the highest density point and run a few sims forward 
######## to look for the likely spread in bats
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1, type = "FORK") #not to overload your computer

nrepeats_postFitting <- 100 #3000
betaBB <- median(cbetaBB)
betaBC <- median(cbetaBC)
betaCB <- 0
#betaCC <- rep(0,nrepeats_postFitting) #TODO: later setted up as constant = 0, so it is pointless in here?
betaCC <- 0

registerDoParallel(cl)

BBinf <- foreach(i=1:nrepeats_postFitting, .combine=cbind) %dopar% {
  #print(c("Start repeat: ",i))
  statusC <- rep(0,nrow(dtc)) # Infected status (0 uninf, 1 inf)
  statusB <- rep(0,nrow(dtb)) # Infected status (0 uninf, 1 inf)

  nfi <- sum(statusC) #number of farm infected
  statusB[1] <- 1 #insert random infection in bats
  #TODO: Shouldn't it really be random and not always the same roost?

  step <- 0
  while(nfi < 1){
    #print(c("Step: ",step))
    tmp <-stepFunction(statusB,statusC,betaBB,betaBC,betaCB,betaCC,gamB,gamC)
    statusB <- tmp[(length(statusC)+1):length(tmp)]
    statusC <- tmp[1:length(statusC)]
    step <- step+1
    nfi <- sum(statusC)
    if (step == 100) break
  }
  return(sum(statusB))
}

stopCluster(cl)

#save.image("Results2.RData")

###################################
### Plots for panel ###
pdf("batCommunityInfections.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

barplot(table(BBinf)/sum(table(BBinf)),xlim=c(0,13),ylim=c(0,.5), ylab = "",xlab="", main="", col="grey",
        space = 0, axes = F, axisnames = F)
axis(1,at=seq(0,12,by=2)+.5,labels=seq(0,12,by=2))
axis(2)

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("Additional number of roosts infected",side=1,cex=1,line=1.2)

dev.off()

######## Do a small run
statusB <- rep(0,nrow(dtb)) # Statuses B setted to 0 (not infected) -- not necessary, but useful if I want to do a small run again
Reqired_nInfB <- 1 #4 TODO: Choose the number reasonably.
#dtb <- dtb %>% filter(y< -21.7) # reduce data
#dtc <- dtc %>% filter(y< -21.7) # reduce data
while(sum(statusB)<Reqired_nInfB){
  statusC <- rep(0,nrow(dtc)) # Infected status (0 uninf, 1 inf)
  statusB <- rep(0,nrow(dtb)) # Infected status (0 uninf, 1 inf)

  nfi <- sum(statusC) #number of farm infected
  statusB[1] <- 1 #insert random infection in bats
  
  step <- 0
  while(nfi < 1){
    #print(c("Step: ",step))
    tmp <-stepFunction(statusB,statusC,betaBB,betaBC,betaCB,betaCC,gamB,gamC)
    statusB <- tmp[(length(statusC)+1):length(tmp)]
    statusC <- tmp[1:length(statusC)]
    step <- step+1
    nfi <- sum(statusC)
    if (step == 100) break
  }
  print(c("Numb of infected B: ", sum(statusB)))
}

### Plot for subset of data
# IDb <- which(dtb$y< -21.7 & dtb$x >-46.9)
# IDc <- which(dtc$y< -21.7 & dtc$x >-46.9)
# 
# pdf("SpatialDistribution.pdf",width=4.64,height=4.64)
# par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
# plot(NA, ylim = c(-22,-21.68), xlim=c(-46.95,-46.6), pch=19, axes = F,ylab = "",xlab="")
# points(dtb$x[IDb],dtb$y[IDb],col=statusB[IDb]+1,pch=19)
# points(dtc$x[IDc],dtc$y[IDc],col=statusC[IDc]+1,pch=17)
# points(dtb$x[1],dtb$y[1],col=3,pch=19) #the first infected root
# 
# axis(1)
# axis(2)
# 
# mtext("Latitude",side=2,cex=1,line=1.2)
# mtext("Longitude",side=1,cex=1,line=1.2)
# 
# legend("bottomleft", c("Roosts","Farms"),
#        col=1, 
#        pch=c(19,17), bty='n',cex=0.75)
# 
# dev.off()


### Plot that we want
ymin <- min(min(dtb$y), min(dtc$y))
ymax <- max(max(dtb$y), max(dtc$y))
xmin <- min(min(dtb$x), min(dtc$x))
xmax <- max(max(dtb$x), max(dtc$x))
empty_margin <- 0.05 # to see the legend properly
pdf("SpatialDistribution.pdf")
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
plot(dtb$x,dtb$y,pch=19, ylim = c(ymin-empty_margin,ymax+empty_margin), xlim=c(xmin-empty_margin, xmax+empty_margin))
points(dtb$x,dtb$y,col=statusB+1,pch=19)
points(dtc$x,dtc$y,col=statusC+1,pch=17)
points(dtb$x[1],dtb$y[1],col=3,pch=19) #the first infected root - TODO: Do we want it red or green?

axis(1)
axis(2)

mtext("Latitude",side=2,cex=1,line=1.2)
mtext("Longitude",side=1,cex=1,line=1.2)

legend("bottomleft", c("Roosts","Farms"),
       col=1,
       pch=c(19,17), bty='n',cex=0.75)

dev.off()

### THE END ###


###################################################################################################################
# ############################### the rest is not finished, we will not do this now
# # another plot  
# plot(dtb$x,dtb$y,pch=19)
# points(dtb$x,dtb$y,col=nStatusB+1,pch=19)
# points(dtc$x,dtc$y,col=nStatusC+1,pch=17)
# # TODO: What was the idea? What should be the nStatuses?
# 
# ###################################################################
# ### load farm image #TODO: to farms look like cows and roots as bats, or smt like that - All the remaining code is related to this?
# PostScriptTrace("farm3.eps") # creates .xml in the working directory
# farmIcon <- readPicture("farm3.eps.xml")
# farmIcon <- readPicture("icons.eps.xml")
# 
# plot(dt$x,dt$y, pty = 's', type = 'n', xlim=c(0,100),ylim=c(0,100))
# 
# xx = grconvertX(x = dt$x, from = 'user', to = 'ndc')
# yy = grconvertY(y = dt$y, from = 'user', to = 'ndc')
# grid.symbols(farmIcon, x=xx,y=yy, size=.05)
# x = dt$x
# y = dt$y
# xyplot(y~x, xlab = "x", ylab = "y",
#        panel = function(x, y) {
#          grid.symbols(farmIcon, x, y, units = "native", size = unit(10, "mm"))
#        })
# 
# points(dt$x,dt$y,col=dt$type, pch=19)
# for(i in 1:nrow(dt)){
#   for(j in 1:nrow(dt)){
#     if (roadM[i,j] != 0){
#       segments(dt$x[i],dt$y[i],dt$x[j],dt$y[j],col=2)
#     }
#     if (roostM[i,j] != 0){
#       segments(dt$x[i],dt$y[i],dt$x[j],dt$y[j],col=1,lty=2)
#     }
#   }
# }
# 
# Status <- rep(0,nrow(dt)) # Infected status (0 uninf, 1 inf)
# 
