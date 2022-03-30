# title: Simulation code for evaluating forest assisted migration under climate change
# author: "Yibiao Zou"
# last updated date: "2022/1/14"


rm(list=ls()) # clear everything
# install.packages()
require(pracma)
require(OpenImageR)
require(rgdal)
require(irr)
library(cvAUC)
library(ROCR)

###############################################################
#-------------------loading data set---------------------------
#-----dispersal kernel------
load(".../GDK_4d.RData")
# GDK, two dimensional dispersal kernel

#-------Climate data during clmate change-------
load(".../DEGD_AM_G_SSP585_Lin.RData")
# DEGD_AM_G, growing degree days for each site and each year during 100 years clmate change

load(".../Dindex_G_SSP585_Lin.RData")
# Dindex_G, drought index for each site and each year during 100 years clmate change

load(".../FireFreq_G_SSP585_Lin.RData")
# FireFreq_G, 3D array of annual fire frequency for each site and each year during 100 years clmate change, for the matrix of each year, row1:no fire; row2:light fire; row3: severe fire

load(".../minWinTemp_G_SSP585_Lin.RData")
# minWinTemp_G, minimum winter temperature for each site and each year during 100 years clmate change

# record number of sites in the study region
nlandscape <- length(Dindex_G[100,])
#-----------assised migration destination-----------
load(".../CSDestination_SSP585_N.RData") 
AMDestination <- CSDestination
# assisted migration destination for each species, which are calculated based on species climatic tolerance parameters (5% or 95% quantile)

#-----------steady state forest before climate change-------
load(".../GFinOutputMicro-1.RData")
# GFinOutputMicro: a young forest used to test the model (40 years)

#-------------Initialize AM related parameters-------
spy_Strat <- 3
spy_Seed <- F

#######################################################################
#--------------------model simulation function-------------------

modelSimulateGW <- function(Forbiome,overTime,P,X,nyear,management="NULL"){ # Global warming scenario
  ecoSymOutput <- ecoSymOutputSetup(P,X)
  nyear1 <- nyear + 100
  for(y in (nyear+1):nyear1){
    dy <- y - nyear
    CCOutput <- timeStep(Forbiome,P,X,DEGD=DEGD_AM_G[dy,],
                         Dindex=Dindex_G[dy,],FireFreq=FireFreq_G[,,dy],
                         minWinTemp=minWinTemp_G[dy,],GDK=GDK,management=management,y,nyear=nyear,Biomass=overTime$Biomass[,,(y-1)])
    Forbiome <- CCOutput$Forbiome
    P <- CCOutput$P
    X <- CCOutput$X
    overTime <- snapShot(Forbiome,overTime,y,P,X)
    P$AMpara$CurPopSize <- rowSums(overTime$Abundance[,,y])
    P$AMpara$CurBasArea <- rowSums(overTime$BasArea[,,y])
    P$AMpara$CurRange <- rowSums(overTime$Range[,,y])
    
    ecoSymOutput <- ecoSymUpd(ecoSymOutput,overTime$Abundance[,,y],overTime$Biomass[,,y],dy)
    
    P$AMpara$Abundance <- overTime$Abundance[,,y]
    newSPool <- overTime$Range[,,y] 
    P$SPool <- updSPool(P$SPool,newSPool)
    P$kDMH <- updDMH(P,Dindex_G[dy,])
  }
  return(list(Forbiome=Forbiome,overTime=overTime,ecoSymOutput=ecoSymOutput,P=P,X=X,nyear=nyear))
}

# ecosystem level output update
ecoSymUpd <- function(ecoSymOutput,Abundance,Biomass,dy){
  nspec <- nrow(Abundance)
  nlandscape <- ncol(Abundance)
  totalSiteAbun <- colSums(Abundance)
  totalSiteBiom <- colSums(Biomass)
  totalSpecAbun <- rowSums(Abundance)
  totalSpecBiom <- rowSums(Biomass)
  totalAbun <- sum(Abundance)
  totalBiom <- sum(Biomass)
  
  RatioAbun <- sapply(1:nlandscape, function(i) Abundance[,i]/totalSiteAbun[i])
  RatioBiom <- sapply(1:nlandscape, function(i) Biomass[,i]/totalSiteBiom[i])
  SIAlpha_A <- sapply(1:nlandscape, function(i) RatioAbun[,i]^2)
  SIAlpha_B <- sapply(1:nlandscape, function(i) RatioBiom[,i]^2)
  SHAlpha_A <- sapply(1:nlandscape, function(i) (-RatioAbun[,i]*log(RatioAbun[,i])))
  SHAlpha_B <- sapply(1:nlandscape, function(i) (-RatioBiom[,i]*log(RatioBiom[,i])))
  
  ecoSymOutput$Alpha_A[dy,] <- 1-colSums(SIAlpha_A)
  ecoSymOutput$Alpha_B[dy,] <- 1-colSums(SIAlpha_B)
  ecoSymOutput$Gamma_A[dy] <- 1-sum((totalSpecAbun/totalAbun)^2)
  ecoSymOutput$Gamma_B[dy] <- 1-sum((totalSpecBiom/totalBiom)^2)
  ecoSymOutput$ShaAlpha_A[dy,] <- colSums(SHAlpha_A)
  ecoSymOutput$ShaAlpha_B[dy,] <- colSums(SHAlpha_B)
  ecoSymOutput$ShaGamma_A[dy] <- sum(-(totalSpecAbun/totalAbun)*log(totalSpecAbun/totalAbun))
  ecoSymOutput$ShaGamma_B[dy] <- sum(-(totalSpecBiom/totalBiom)*log(totalSpecBiom/totalBiom))
  ecoSymOutput$TBiomass[dy] <- sum(Biomass) 
  return(ecoSymOutput)
}

ecoSymOutputSetup <- function(P,X){
  nspec <- P$nspec
  nlandscape <- X$nlandscape
  Alpha_A <- matrix(0,nrow=100,ncol=nlandscape) # calculate alpha diversity based on abundance
  Alpha_B <- matrix(0,nrow=100,ncol=nlandscape)
  Gamma_A <- rep(0,100)
  Gamma_B <- rep(0,100)
  ShaAlpha_A <- matrix(0,nrow=100,ncol=nlandscape)
  ShaAlpha_B <- matrix(0,nrow=100,ncol=nlandscape)
  ShaGamma_A <- rep(0,100)
  ShaGamma_B <- rep(0.100)
  TBiomass <- rep(0,100)
  ecoSymOutput <- list(Alpha_A=Alpha_A,Alpha_B=Alpha_B,
                       Gamma_A=Gamma_A,Gamma_B=Gamma_B,
                       ShaAlpha_A=ShaAlpha_A,ShaAlpha_B=ShaAlpha_B,
                       ShaGamma_A=ShaGamma_A,ShaGamma_B=ShaGamma_B,
                       TBiomass=TBiomass)
  return(ecoSymOutput)
}

snapShot <- function(Forbiome,overTime,y,P,X){
  nid <- P$nid
  nspec <- P$nspec
  preThreshold <- X$preThreshold
  ha <- X$ha
  nlandscape <- X$nlandscape
  preDensity <- X$preDensity
  for(l in 1:nlandscape){
    if(isempty(Forbiome[[l]])){
      overTime$BasArea[,l,y] <- 0
      overTime$Biomass[,l,y] <- 0
      overTime$Abundance[,l,y] <- 0
    }else{
      Forbiome[[l]] <- as.matrix(Forbiome[[l]])
      bas <- 0.25*pi*Forbiome[[l]][nid$dbh,]^2/10000
      bas[bas < preThreshold] <- 0
      biomass <- 0.2*Forbiome[[l]][nid$dbh,]^2.4 + 
        Forbiome[[l]][nid$kC1,]*Forbiome[[l]][nid$kA1,]*Forbiome[[l]][nid$dbh,]^Forbiome[[l]][nid$kA2,] # dry stemwood biomass + foliage biomass
      overTime$BasArea[,l,y] <- sapply(1:nspec, function(i) sum(bas[Forbiome[[l]][nid$id,]==i]))
      overTime$Biomass[,l,y] <- sapply(1:nspec, function(i) sum(biomass[Forbiome[[l]][nid$id,]==i]))
      overTime$Abundance[,l,y] <- sapply(1:nspec, function(i) sum((Forbiome[[l]][nid$id,]==i)*1))
    }
  }
  overTime$Range[,,y] <- (overTime$BasArea[,,y] > preDensity*ha)*1
  return(overTime)
}

timeStep <- function(Forbiome,P,X,DEGD,Dindex,FireFreq,minWinTemp,GDK,management="NULL",y,nyear,validation=F,Biomass=NULL){
  SPool <- P$SPool
  nlandscape <- X$nlandscape
  nspec <- P$nspec
  dkSPool <- t(sapply(1:nspec, function(i) DimTrans_LDD(SPool[i,],nlandscape,GDK[[i]])))
  
  # check whether there is fire (no fire at the first stage)
  if(!isempty(FireFreq)){ 
    FireRegime <- sapply(1:nlandscape,function(i) sample(0:2,size=1,prob=FireFreq[,i]))
  }else{FireRegime <- NULL}
  X$FireRegime <- FireRegime
  kB <- P$kDMH$kB
  kHmax <- P$kDMH$kHmax
  
  if(management=="AM"){
    dy=y-nyear
    P$LandFireLY[FireRegime > 0] <- dy
    AMOut <- AM(Forbiome,P,X,dy=dy,Biomass=Biomass,dkSPool,Strategy=spy_Strat,Seed=spy_Seed)
    Forbiome <- AMOut$Forbiome
    P <- AMOut$P
    dkSPool <- AMOut$dkSPool
  }
  
  Forbiome <- lapply(1:nlandscape, function(i) 
    lifeCycle(Forbiome[[i]],P,X,DEGD[i],Dindex[i],FireRegime[i],minWinTemp[i],dkSPool[,i],kB[,i],kHmax[,i]))
  
  return(list(Forbiome=Forbiome,P=P,X=X))
}

# LifeCycle of each tree: recruitment -> growth -> mortality
lifeCycle <- function(Forest,P,X,DEGD,Dindex,FireRegime,minWinTemp,dkSPool,kB,kHmax){
  if((sum(dkSPool)==0)&&(isempty(Forest))){return(Forest)} # dkSPool == 0 only if Forest == numeric(0). And because dkSPool == 0, there will be no recruitment
  if(!isempty(Forest)){Forest <- as.matrix(Forest)}
  Forest <- Recruitment(Forest,P,X,minWinTemp,dkSPool,kB)
  if(!isempty(Forest)){
    Forest <- Growth(Forest,P,X,DEGD,Dindex,kB,kHmax)
    Forest <- Mortality(Forest,P,X,FireRegime)
    if(!isempty(Forest)){Forest <- as.matrix(Forest)}
  }
  return(Forest)
}

DimTrans_LDD <- function(ODVector,nlandscape,GDK){ # function input: one dimension vector
  TDMatrix <- matrix(NA,nrow=nlandscape/12,ncol=12)
  # NSite <- nlandscape
  FinVector <- rep(NA,nlandscape)
  for(i in 1:nlandscape){
    TDMatrix[i] <- ODVector[i]
  }
  BDMatrix <- matrix(0, nrow=(nlandscape/12+2), ncol=14)
  BDMatrix[2:(nlandscape/12+1), 2:13] <- TDMatrix
  NewMatrix <- convolution(image=BDMatrix, kernel=GDK, mode = "same")
  TDMatrix <- NewMatrix[2:(nlandscape/12+1), 2:13]
  for(j in 1:nlandscape){
    FinVector[j] <- TDMatrix[j]
  }
  return(FinVector)
}

#--------------Recruitment-------------------------- 
Recruitment <- function(Forest,P,X,minWinTemp,dkSPool,kB){
  nid <- P$nid
  nspec <- P$nspec
  plotVector <- X$plotVector
  nplot <- X$nplot
  phi <- P$phi
  maxNewtrees <- X$maxNewtrees
  params <- P$params
  
  if(isempty(Forest)){AL_birth <- rep(1,nplot)
  }else{
    la <- leafArea(Forest,nid)
    sumla_birth <- sapply(1:nplot, function(i) sum(la[Forest[nid$pLocation,]==i]))
    AL_birth <- sapply(1:nplot, function(i) phi*exp(-0.25*sumla_birth[i]/833))
  }
  
  # filter
  FilterFlag <- matrix(NA,nspec,nplot)
  flag2 <- minWinTemp < P$kWiTN 
  flag3 <- minWinTemp > P$kWiTX # winter temperature should not be too low or too high
  for(i in 1:nplot){
    flag1 <- P$kLy > AL_birth[i]
    FilterFlag[,i] <- flag1|flag2|flag3 # light availability on the forest floor should not be too low
  }
  Filterindex <- sapply(1:nplot, function(i) sum(dkSPool[!FilterFlag[,i]]) > 0)
  Filter_plotVector <- plotVector[Filterindex]
  
  if(length(Filter_plotVector)>0){
    newTrees_plot_total <- sapply(1:length(Filter_plotVector), function(i) sample(2:maxNewtrees,size=1))
    sum_newTrees = sum(newTrees_plot_total)
    newForest = matrix(0,nrow=19,ncol=sum_newTrees)
    # Determine the species_id, initial size, and maximum age etc.
    newForest[2,] <- rep(Filter_plotVector,newTrees_plot_total)
    newForest[1,] <- sapply(1:sum_newTrees, function(i) newForest[1,i]=sample(P$spec_id*(!FilterFlag[,newForest[nid$pLocation,i]]),size = 1,
                                                                              prob = dkSPool*(!FilterFlag[,newForest[nid$pLocation,i]])/sum(dkSPool[!FilterFlag[,newForest[nid$pLocation,i]]])))
    newForest[3:15,] <- params[,newForest[1,]]
    newForest[16,] <- P$birthSize 
    gs <- newForest[nid$kSmin,] + newForest[nid$kE1,]*(1-AL_birth[newForest[nid$pLocation,]])
    kC <- -gs/kB[newForest[1,]]
    newForest[17,] <- 137 + kB[newForest[1,]]*(1-exp(kC*newForest[nid$dbh,]))
    # Combine the new matrix with the old one
    Forest <- cbind(Forest,newForest)
  }
  return(Forest)
}

#-----------------Growth-------------------------
Growth <- function(Forest,P,X,DEGD,Dindex,kB,kHmax){
  # S <- P$S
  # DEGD <- DEGD + S * rnorm(1,0,1) # growing degree days
  a <- P$a
  ntrees <- ncol(Forest)
  nid <- P$nid
  phi <- P$phi
  
  if(ntrees > 0){
    dinc = rep(NA,ntrees) # Actual growth increment
    gDDGF = rep(NA,ntrees) # growing degree days limiting factor
    gALGF = rep(NA,ntrees) # light availability limiting factor
    gSMGF = rep(NA,ntrees) # soil moisture limiting factor
    
    # light availability of the trees
    SLAR <- SLARC(Forest,nid=nid,ntrees,X$nplot)
    AL = phi*exp(-SLAR/833*0.25) # available light
    AL[AL > 1] <- 1
    
    # calculate amount of growth for each tree
    hMax <- kHmax[Forest[nid$id,]]
    ht <- Forest[nid$ht,]
    gs <- Forest[nid$kSmin,] + Forest[nid$kE1,]*(1-AL)
    kC <- -gs/kB[Forest[nid$id,]]
    fh <- gs*(1-(ht-137)/(hMax-137))
    fh[fh < 0] <- 0
    
    # potential increment
    dnc <- Forest[nid$G,]*Forest[nid$dbh,]*(1-ht/hMax)/(2*ht+fh*Forest[nid$dbh,])
    dnc[dnc < 0] <- 0
    
    # effect of temperature
    gDDGF <- 1-exp((Forest[nid$dMin,]-DEGD)*a)
    gDDGF[gDDGF<0] <- 0
    # effect of water defficiency/drought
    gSMGF <- sapply(1:ntrees, function(i) max(c(1-Dindex/Forest[nid$kDrT,i],0))^(1/2))
    # effect of light availability
    gALGF <- 1-exp(-4.64*(AL-0.05)) + (Forest[nid$kLa,] - 1)*((2.24*(1-exp(-1.136*(AL-0.08))))-(1-exp(-4.64*(AL-0.05))))/8
    gALGF[gALGF < 0] <- 0
    
    # Actual increment
    dinc <- dnc*(gDDGF*gALGF*gSMGF)^(1/3)
    # individuals that don't grow much are set to 0 and considered positive for "no growth"
    dinc[dinc<=0.1*dnc] <- 0 # use a new threshold: 0.1*dnc, rather than use 1mm 
    Forest[nid$nogrowth,dinc<=0.1*dnc] <- 1
    # Finally, add the growth onto the DBH of the tree
    Forest[nid$dbh,] = Forest[nid$dbh,] + dinc
    Forest[nid$ht,] <- ht + fh*dinc
  }
  return(Forest)
}

#-----------------Mortality-----------------------
Mortality <- function(Forest,P,X,FireRegime){
  ntrees <- ncol(Forest)
  unifs <- runif(3*ntrees,0,1)  
  randDeathCheck <- unifs[1:ntrees] # These are for age-related death
  noGrowDeathCheck <- unifs[(ntrees+1):(2*ntrees)] # These are for no-growth/stress-induced death
  randDeath <- P$randDeath
  noGrowDeath <- P$noGrowDeath
  nid <- P$nid
  fsev <- X$fsev
  
  if(!isempty(FireRegime)){
    FireDeathCheck <- unifs[(2*ntrees+1):(3*ntrees)] # These are for fire-induced death  
    
    # Fire index: to check the firestate of each site and the fire tolerance class of each species
    Findex1 <- Forest[nid$kFi,]==1 & (FireRegime==1|FireRegime==2)
    Findex21 <- Forest[nid$kFi,]==2 & FireRegime==1
    Findex22 <- Forest[nid$kFi,]==2 & FireRegime==2
    Findex31 <- Forest[nid$kFi,]==3 & FireRegime==1
    Findex32 <- Forest[nid$kFi,]==3 & FireRegime==2
    Findex41 <- Forest[nid$kFi,]==4 & FireRegime==1
    Findex42 <- Forest[nid$kFi,]==4 & FireRegime==2
    
    # Fire induced mortality
    FireDP <- rep(0,ntrees) # Fire induced death probability
    if(any(Findex1)){FireDP[Findex1] <- 1}
    if(any(Findex21)){FireDP[Findex21] <- exp(((-(1-fsev[1])*0.00202)-0.00053)*Forest[nid$dbh,Findex21])}
    if(any(Findex22)){FireDP[Findex22] <- exp(((-(1-fsev[2])*0.00202)-0.00053)*Forest[nid$dbh,Findex22])}
    if(any(Findex31)){FireDP[Findex31] <- exp(((-(1-fsev[1])*0.02745)-0.00255)*Forest[nid$dbh,Findex31])}
    if(any(Findex32)){FireDP[Findex32] <- exp(((-(1-fsev[2])*0.02745)-0.00255)*Forest[nid$dbh,Findex32])}
    if(any(Findex41)){FireDP[Findex41] <- exp(-0.00053*Forest[nid$dbh,Findex41]) - 0.5 - (1-fsev[1])*0.5}
    if(any(Findex42)){FireDP[Findex42] <- exp(-0.00053*Forest[nid$dbh,Findex42]) - 0.5 - (1-fsev[2])*0.5}
    FireDP[FireDP<0] <- 0
    # Identify which individuals die
    deadTrees <- which(randDeathCheck<randDeath/Forest[nid$ageMax,] |  # because of age-related death
                         (noGrowDeathCheck<noGrowDeath & Forest[nid$nogrowth,]==1) | # or because of no-growth/stress-induced death
                         FireDeathCheck<FireDP) # because of fire
  }else{
    deadTrees <- which(randDeathCheck<randDeath/Forest[nid$ageMax,] |  # because of age-related death
                         (noGrowDeathCheck<noGrowDeath & Forest[nid$nogrowth,]==1))
  }
  if(!isempty(deadTrees)){
    Forest <- Forest[,-deadTrees]
  }
  
  # all nogrowth is reset
  if(!isempty(Forest)&is.matrix(Forest)){
    Forest[nid$nogrowth,] = 0  # 1 -> not enough growth; 0 -> enough growth
  }
  
  return(Forest)
}

# Calculate light availability for each tree 
SLARC <- function(Forest,nid,ntrees,nplot){
  la <- leafArea(Forest,nid=nid)
  iht = floor(Forest[nid$ht,]/10) # the unit here is decimeter, used to reduce the overestimate competition among trees that have similar height
  # calculate light availability for each tree; trees are only shaded by higher trees in the same site
  SLAR = rep(0,ntrees)
  for(i in 1:nplot){
    SI <- Forest[nid$pLocation,] == i # plot index
    if(length(iht[SI]) > 0){
      ihtMax <- max(iht[SI])
      sumla <- sapply(1:ihtMax, function(h) sum(la[SI & iht>h])) # sum leaf area/shading above each tree
      SLAR[SI] <- sumla[iht[SI]]
    }
  }
  return(SLAR)
}

# calculate leaf area based on dbh and allometric parameter
leafArea <- function(Forest,nid){
  # leaf area of the trees
  la <- Forest[nid$kC2,]*Forest[nid$kA1,]*Forest[nid$dbh,]^Forest[nid$kA2,]
  return(la)
}

# function to calculate AM distance
AMDistance <- function(TargetSite,j,LLT=130){
  x1 <- j%/%LLT
  y1 <- j%%LLT
  x2 <- TargetSite%/%LLT
  y2 <- TargetSite%%LLT
  AMDST <- ((x1-x2)^2+(y1-y2)^2)^0.5
  return(AMDST)
}

#----------------Assisted migration----------------
AM <- function(Forbiome,P,X,dy,Biomass,dkSPool,Strategy=1,Seed=F){
  CurPopSize <- P$AMpara$CurPopSize
  OriPopSize <- P$AMpara$OriPopSize
  CurBasArea <- P$AMpara$CurBasArea
  OriBasArea <- P$AMpara$OriBasArea
  CurRange <- P$AMpara$CurRange
  OriRange <- P$AMpara$OriRange
  SiteBiom <- colSums(Biomass)
  # AMpropo <- P$AMpropo
  AMpropo <- 0.4
  nid <- P$nid
  kB <- P$kDMH$kB
  landscape <- seq(1,X$nlandscape,by=1)
  nplot <- X$nplot
  # Abundance <- overTime$Abundance # new
  # minAMdbh <- X$minAMdbh 
  # maxAMdbh <- X$maxAMdbh
  minAMdbh <- 1.27
  maxAMdbh <- 25.4
  Abundance <- P$AMpara$Abundance
  r1 <- 0.3
  r2 <- 0.3
  r3 <- 0.3
  LandFireLY <- P$LandFireLY
  
  # AMspec <- as.vector(P$spec_id[CurPopSize < r1*OriPopSize|CurBasArea < r2*OriBasArea|CurRange < r3*OriRange|CurRange<=2])
  AMspec <- as.vector(P$spec_id[CurPopSize < r1*OriPopSize|CurBasArea < r2*OriBasArea|CurRange < r3*OriRange])
  AMits <- 4 # AM interval (Sensitivity analysis variable)
  AMcoy <- 4  # AM continuous years (Sensitivity analysis variable)
  AMSiteNum <- 10
  # I might assume each year for each species, there is maximum 1000 individual seedlings available for AM
  AMseedlingNum <- round(100*nplot/12) # mean number of indiviual seedling per hectare that would be assisted migrated
  SeedAdvantage <- 30
  
  AMlog <- P$AMpara$AMlog
  AMSiteLog <- P$AMSiteLog
  FireRegime <- X$FireRegime
  AMState <- P$AMState
  
  # if(length(AMspec)==0){return(list(Forbiome=Forbiome,P=P,dkSPool=dkSPool))}
  if(length(AMspec)!=0){
    AMspecUp <- AMspec[AMState[AMspec,dy]==0]
    if((dy+AMcoy+AMits-1)<=100){
      AMState[AMspecUp,dy:(dy+AMcoy-1)] <- 1
      AMState[AMspecUp,(dy+AMcoy):(dy+AMcoy+AMits-1)] <- -1 # AM frozen period
    }else if((dy+AMcoy+AMits-1)>100 && (dy+AMcoy)<=100){
      AMState[AMspecUp,dy:(dy+AMcoy-1)] <- 1
      AMState[AMspecUp,(dy+AMcoy):100] <- -1
    }else if((dy+AMcoy)>100){
      AMState[AMspecUp,dy:100] <- 1
    }
  }
  
  P$AMState <- AMState
  AMspecA <- P$spec_id[AMState[,dy]==1]
  for(ts in AMspecA){
    tick <- 0
    # AMsite <- landscape[Abundance[ts,]>0 & AMDestination[ts,]==0]
    TargetSite <- landscape[AMDestination[ts,]==1]
    FireTargetSite <- landscape[AMDestination[ts,]==1 & FireRegime>0]
    AMsite <- landscape[Abundance[ts,]>0 & AMDestination[ts,]==0]
    LandFireLY[AMDestination[ts,]==0] <- -1 
    if(Seed==F){ # seedling assisted migration
      for(i in 1:AMSiteNum){
        j <- sample(AMsite,size=1,prob=rep(1/length(AMsite),length(AMsite)))
        AMForest <- matrix(0,nrow=19,ncol=AMseedlingNum)
        AMForest[1,] <- rep(ts,AMseedlingNum)
        AMForest[2,] <- sample(1:nplot,size=AMseedlingNum,prob=rep(1/nplot,nplot),replace=T)
        AMForest[3:15,] <- P$params[,AMForest[1,]]
        AMForest[16,] <- 2.5 
        gs <- AMForest[nid$kSmin,] + 0.75*AMForest[nid$kE1,]
        kC <- -gs/kB[AMForest[1,]]
        AMForest[17,] <- 137 + kB[AMForest[1,]]*(1-exp(kC*AMForest[nid$dbh,]))
        # post-fire would be added here
        if(Strategy==1){ # least competition/open canopy assisted migration
          RsLocation<- order(SiteBiom[TargetSite])[i]
          Rs <- TargetSite[RsLocation]
          # browser()
          if(length(Rs)==0||is.na(Rs)){next}
          Forbiome[[Rs]] <- cbind(Forbiome[[Rs]],AMForest)
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==2){ # minimum distance assisted migration
          RsLocation <- order(abs(AMDistance(TargetSite,j)))[i]
          Rs <- TargetSite[RsLocation]
          if(length(Rs)==0||is.na(Rs)){next}
          # browser()
          Forbiome[[Rs]] <- cbind(Forbiome[[Rs]],AMForest)
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==3){ # Post-fire assisted migration strategy
          Rs <- order(LandFireLY,decreasing=T)[i]
          if(length(Rs)==0||is.na(Rs)){next}
          Forbiome[[Rs]] <- cbind(Forbiome[[Rs]],AMForest)
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==4){ # Least-fire assisted migration strategy
          LandFireLY[LandFireLY < 0] <- 101 
          Rs <- order(LandFireLY,decreasing=F)[i]
          if(length(Rs)==0||is.na(Rs)){next}
          Forbiome[[Rs]] <- cbind(Forbiome[[Rs]],AMForest)
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }
        tick <- tick + 1
      }
    }
    if(Seed==T){ # seed assisted migration
      for(i in 1:AMSiteNum){
        j <- sample(AMsite,size=1,prob=rep(1/length(AMsite),length(AMsite)))
        if(Strategy==1){
          RsLocation<- order(SiteBiom[TargetSite])[i]
          Rs <- TargetSite[RsLocation]
          if(length(Rs)==0||is.na(Rs)){next}
          dkSPool[ts,Rs] <- SeedAdvantage
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==2){
          RsLocation <- order(abs(AMDistance(TargetSite,j)))[i]
          Rs <- TargetSite[RsLocation]
          if(length(Rs)==0||is.na(Rs)){next}
          dkSPool[ts,Rs] <- SeedAdvantage
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==3){
          Rs <- order(LandFireLY,decreasing=T)[i]
          if(length(Rs)==0||is.na(Rs)){next}
          dkSPool[ts,Rs] <- SeedAdvantage
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }else if(Strategy==4){
          LandFireLY[LandFireLY < 0] <- 101 
          Rs <- order(LandFireLY,decreasing=F)[i]
          if(length(Rs)==0||is.na(Rs)){next}
          dkSPool[ts,Rs] <- SeedAdvantage
          AMSiteLog[dy,Rs,ts] <- AMSiteLog[dy,Rs,ts] + 1
        }
        tick <- tick + 1
      }
    }
    if(tick>0){P$AMpara$AMlog[ts,dy] <- 1}
  }
  P$AMpara$AMrecord[dy] <- 1
  P$AMSiteLog <- AMSiteLog
  return(list(Forbiome=Forbiome,P=P,dkSPool=dkSPool))
}

# update species pool at the end of each year; might be edited*
updSPool <- function(SPool,newSPool){
  SPool <- ((SPool + newSPool) >=1)*1
  return(SPool)
}

updDMH <- function(P,Dindex){
  kB <- P$kDMH$kB
  kHmax <- P$kDMH$kHmax
  kRedMax <- P$kRedMax
  kDrT <- P$kDrT
  hMax <- P$params[3,]
  nlandscape <- length(Dindex)
  for(i in 1:nlandscape){
    d = (kDrT-Dindex[i])/kDrT
    d[d<0] <- 0
    kHmax[,i] <- d*hMax*(1-kRedMax)+hMax*kRedMax
    kB[,i] <- kHmax[,i] - 137 
  }
  return(list(kHmax=kHmax,kB=kB))
}

#################################################################
#---------------------Simulation experiment and output---------
#------------climate change and AM scenario------
BCCOutput <- GFinOutputMicro$validationTest # extract steady state forest and species-level output recording data-frame
BCCOutput$P$AMState <- matrix(0,nrow=23,ncol=100)
BCCOutput$P$AMSiteLog <- array(0,dim=c(100,1560,23))
BCCOutput$P$LandFireLY <- rep(0,1560)

# run AM simulation
CCAMTest1 <- modelSimulateGW(BCCOutput$Forbiome,BCCOutput$overTime,BCCOutput$P,
                            BCCOutput$X,BCCOutput$nyear)

# store output
nyear <- CCAMTest1$nyear
CCAMTestSO <- list(Biomass=CCAMTest1$overTime$Biomass[,,(nyear+1):(nyear+100)],
                   Abundance=CCAMTest1$overTime$Abundance[,,(nyear+1):(nyear+100)],
                   Range=CCAMTest1$overTime$Range[,,(nyear+1):(nyear+100)],
                   AMlog=CCAMTest1$P$AMpara$AMlog,
                   ecoSymOutput=CCAMTest1$ecoSymOutput)

save(CCAMTestSO, file=".../CCAMTestSO.RData")
