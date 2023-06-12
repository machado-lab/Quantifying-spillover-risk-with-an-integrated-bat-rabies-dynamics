## Functions for PrelimSims

###############################################################
# FOR RANDOM POP

 generatePops <- function(nPops=10){
   x <- runif(nPops,0,100)
   y <- runif(nPops,0,100)
   name <- LETTERS[1:nPops] #Only nPops < 27
   pop <- round(rnorm(nPops,100,sd=10))
   type <- rbinom(nPops,1,.5)+1 #1 roosts, 2 farms
   
   return(data.frame(name,pop,x,y,type))
 }

# ## For the pops random above
 buildRoadNetwork <- function (dt){
#   # random weight (0 to 1) only for type 2 (farms); when x is type 2 (the farm Im interested in)
   roadNetwork <- sapply(1:nrow(dt), function(x){(dt$type-1)*runif(nrow(dt),0,1)*(dt$type[x]-1)})
   diag(roadNetwork)<-0 # set diagonal to zero
   return(roadNetwork)
 }

 buildRoostNetwork <- function(dt){
#   # weight is pure euclidean distance, only for type 1 (roosts); when x is a roost 
   roostNetwork <- sapply(1:nrow(dt),function(i){(-dt$type+2)*sqrt((dt$x-dt$x[i])^2+(dt$y-dt$y[1])^2)/100*(-dt$type[i]+2)})
   return(roostNetwork)
 }

###############################################################
## For the read data loaded

## We don't need the network between farm, since we do not suppose the transmission between farms due to the movement
# buildRoadNetworkR <- function (dt){
#   # random weight (0 to 1) accross the network, on average each farm is connected to 70% of the network
#   roadNetwork <- sapply(1:nrow(dt), function(x){rbinom(nrow(dt),1,.7)*runif(nrow(dt),0,1)})
#   diag(roadNetwork)<-0 # set diagonal to zero
#   return(roadNetwork)
# }

### Build network between roosts according to distance & number of bats & sdm 
buildRoostNetworkR <- function(dt){
  # weight is 1/euclidean distance
  roostNetwork_dist <- sapply(1:nrow(dt),function(i){1/sqrt((dt$x-dt$x[i])^2+(dt$y-dt$y[i])^2)})
  number_bats_M <- matrix(dt$numb_bats, length(dt$numb_bats), 1)
  amount_weight <- number_bats_M %*% c(number_bats_M)
  sdm_vec <- matrix(dt$sdm, length(dt$sdm), 1)
  average_sdm_weight <- (sdm_vec %*% c(sdm_vec))/2
  roostNetwork_DistNumb <- roostNetwork_dist * amount_weight * average_sdm_weight
  roostNetwork_DistNumb[roostNetwork_dist<11.119] <- 0 # limit transmission to 10km distance
  diag(roostNetwork_DistNumb)<-0 # set diagonal to zero (no transmission to itself)
  return(roostNetwork_DistNumb)
}

### Build network from Roosts to Farms according to distance & altitude
# (but limited by 10km distance independent of the altitude)
buildRoostToFarmNetworkR <- function(dt_B, dt_C){
  # weight is 1/euclidean distance (without altitudes) if less then ~10km; rows=Roosts, columns=Farm
  roostToFarmNetwork_dist <- sapply(1:nrow(dt_C),function(i){1/sqrt((dt_B$x-dt_C$x[i])^2+(dt_B$y-dt_C$y[i])^2)})

  # change weight to 1/euclidean distance INCLUDING ALTITUDES and it is 0 if C is higher than B
  roostToFarmNetwork_DistAltit <- sapply(1:nrow(dt_C),function(i){ifelse(dt_B$elev<dt_C$elev[i],0,1/sqrt((dt_B$x-dt_C$x[i])^2+(dt_B$elev-dt_C$elev[i])^2))})
  roostToFarmNetwork_DistAltit[roostToFarmNetwork_dist<11.119] <- 0 # limit transmission to 10km of original distance
  return(roostToFarmNetwork_DistAltit)
  #
  
  ## If we want to OMIT THE ALTITUDES
  # roostToFarmNetwork_dist[roostToFarmNetwork_dist<11.119] <- 0 # limit transmission to 10km of original distance
  # return(roostToFarmNetwork_dist)
  ##
}

stepFunction <- function(statusB,statusC,betaBB,betaBC,betaCB,betaCC,gamB,gamC){
 
 # updating status of farms (cattle)
  nStatusC <- sapply(1:length(statusC),function(x){return(ifelse(statusC[x]==1,pRecov(x,gamC),pInfectC(x,betaBC,betaCC)))})

 # updating status of bats
  nStatusB <- sapply(1:length(statusB),function(x){return(ifelse(statusB[x]==1,pRecov(x,gamB),pInfectB(x,betaCB,betaBB)))})

 return(c(nStatusC,nStatusB))
}

## Probability of recovery (same function for cattle & bats)
pRecov <- function(ID,gam){
  return(ifelse(runif(1) < 1-exp(-gam),0,1))
}

## Probability of infection cattle (due to bats) # We now do not allow cattle to be infected by cattle
pInfectC <- function(ID,betaBC,betaCC){
# only BC transmission:
  beta <- betaBC * sum(RoostFarmM[,ID] * statusB)
  return(ifelse(runif(1) < 1-exp(-beta),1,0))
}

## Probability of infection bats (due to bats) # We now do not allow bats to be infected by cattle
pInfectB <- function(ID,betaBB,betaCB){
# CB & BB transmission:  
  #infC <- which(statusC==1)
  #dist <- sapply(1:length(infC), function(i){1/sqrt((dtb$x[ID]-dtc$x[infC[i]])^2+(dtb$y[ID]-dtc$y[infC[i]])^2)})
  #dist[dist<11.119] <- 0 # limit transmission to 10km
  #beta <- ifelse(length(infC)>0,betaCB * sum(dist), 0) + betaBB * sum(roostM[,ID] * statusB)
# only BB transmission:
	beta <- betaBB * sum(roostM[,ID] * statusB)
  return(ifelse(runif(1) < 1-exp(-beta),1,0))
}

##########Functions to generate random df for bats and cattle:
#####Create function for simulated df- one for cattle one for bats
#For bats
bat_generatePops <- function(nPops=10) {
  x <- runif(nPops, 0, 100)
  y <- runif(nPops, 0, 100)
  sdm <- runif(nPops, min=0, max=1)
  ID <- LETTERS[1:nPops]
  elev <- round(rnorm(nPops, 100, sd=10))
  group_type <- sample(c("Bachelor", "Harem", "Overnight"), 10, replace = TRUE)
  
  return(data.frame(ID, x, y, sdm, elev, group_type))
}

#For cattle:
cattle_generatePops <- function(nPops=10) {
  x <- runif(nPops, 0, 100)
  y <- runif(nPops, 0, 100)
  #sdm <- runif(nPops, min=0, max=1)
  ID <- LETTERS[1:nPops]
  elev <- round(rnorm(nPops, 100, sd=10))
  #group_type <- sample(c("Bachelor", "Harem", "Overnight"), 10, replace = TRUE)
  totalAnimals <- round(rnorm(nPops, 100, sd=10))
  
  return(data.frame(ID, x, y, elev, totalAnimals))
}
