# Gwaai Haanas sea otter meta-population model
# Summary:
# 
rm(list = ls())
# Set User Parameters  ---------------------------------------------------------
reps = 100           # Number replications for population sims (should use at least 100)
Nyrs = 25            # Number of years to project population dynamics
ImigratOpt = 1       # Immigration option: 0 = none, 1 = low, 2 = high
MnImLo = 1           # Mean annual immigrants with low immigration option
MnImHi = 3           # Mean annual immigrants with high immigration option
V_mn = 2.5           # Population front asymptotic wavespeed, km/yr, min  
V_mx = 5.5           # Population front asymptotic wavespeed, km/yr, max
K_mean = 3           # Baseline mean K (can modify as fxn of habitat variables)
K_sig = 1.5          # Standard deviation in local K (variation over space)
sig = 0.05           # Environmental stochasticity (std dev in log-lambda)
rmax = log(1.22)     # Maximin rate of growth = 22% per year
theta = 1            # theta parameter for theta-logistic (1 = Ricker model, >1 = delayed DD) 
Yr1 = 2018           # Calendar Year to begin simulations at
Initpop = 10         # Number of animals in initial population (at least 2 adult females)
Initblk = c(1)       # List of initially occupied bloaks: e.g. c(1) = Block 1 only
# ~~~~~~END User parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Load necessary libraries -------------------------------------------------
library(gtools)
library(BMS)
library(boot)
library(ggplot2)
library(ggmap)
library(reshape2)
library(rgdal)
library(raster)
#
# Load files --------------------------------------------------------------
data = read.csv("GHBlockdata.csv", header = TRUE)
Cdata = read.csv("GHCelldata.csv", header = TRUE)
Demdat = read.csv("RandDem.csv", header = TRUE)
Distmat =  read.csv("Distmat.csv", header = FALSE)
# Probabilities of dispersal from each block for each age/sex class
DispP = read.csv("GHDispProb.csv", header = TRUE) 
# Inter-pop movement matrices: pairwise prob of dispersal based on 
#  pairwise distances and dispersal kernels for each age/sex class
destJF = read.csv("GHDispMatJF.csv", header = FALSE) 
destAF = read.csv("GHDispMatAF.csv", header = FALSE);
destJM = read.csv("GHDispMatJM.csv", header = FALSE);
destAM = read.csv("GHDispMatAM.csv", header = FALSE);
Habavg = read.csv("HabAvg.csv", header = TRUE);
# HAIDA <- raster("HGland.grd") # Raster of Gwaii Haanas, NAD_1983_Albers (BC Environment)
#
# Process data ----------------------------------------------------------------
Dispers = 2.5  # Over-Dispersion param for Neg Binomial # immigrants per year
pparLo = Dispers/(Dispers+MnImLo/length(Initblk))
pparHi = Dispers/(Dispers+MnImHi/length(Initblk))
Years = c(Yr1:(Yr1+Nyrs-1))  
Yrs = seq(1:Nyrs)
Years = Yrs-1+Yr1
P = dim(Distmat)[1]  # number blocks (or sub-populations)
N0 = numeric(length = P)
for (i in 1:length(Initblk)){
  N0[Initblk[i]] = round(Initpop)/length(Initblk)
}
#  Create Data variables for calculating relative density of Hab cells:
PUID = Cdata$PU_ID
Blk = Cdata$BlockID
area = Cdata$Area
dep = Cdata$DEPTH
fetch = Cdata$Fetch
botm = as.character(Cdata$BT_Code)
# Load params for habitat density at K fxn
params = read.csv("Hab_params.csv")
parms = params$Parms
# Define Kcalc function
Kcalc <- function(PUID,Blk,area,dep,botm,fetch,parms,Kmn){
  b <- parms  
  # Hab dens multiplier fxn: exp(b1*X1 + b2*X2... + bn*Xn)
  # Depth part of fxn: b1*(-1*dep) - b2*dep^2
  # where sum(area*exp(b1*(-1*dep)-b2*(dep^2)))/sum(area) =~ 1
  Depfxn = b[1]*(-1*dep) - b[2]*dep^2
  # Fetch part of fxn
  fch = fetch - mean(fetch)
  Fchfxn = b[3]*(fch) 
  btfxn = numeric(length = length(botm))
  btfxn[which(botm=="1")] = b[4]; btfxn[which(botm=="1a")] = b[5]; btfxn[which(botm=="1b")] = b[6]
  btfxn[which(botm=="2")] = b[7]; btfxn[which(botm=="2a")] = b[8]; btfxn[which(botm=="2b")] = b[9]
  btfxn[which(botm=="3")] = b[10]; btfxn[which(botm=="3a")] = b[11]; btfxn[which(botm=="3b")] = b[12]
  btfxn[which(botm=="0")] = b[13]
  mult = exp(Depfxn + Fchfxn + btfxn) 
  # NOTE: sum((mult*area))/sum(area) should equal approx 1
  Kdns = numeric(length = max(Blk))
  Ktot = numeric(length = max(Blk))
  AreaB = numeric(length = max(Blk))
  for (i in 1:max(Blk)){
    ii = which(Blk == i)
    Kdns[i] = mean(Kmn*mult[ii])
    Ktot[i] = sum(Kmn*area[ii]*mult[ii])
    AreaB[i] = sum(area[ii])
    Kdns[i] = Ktot[i]/sum(area[ii])
  }
  Habdns = data.frame(PUID = PUID, Reldens = mult)
  Ktab = data.frame(Block = seq(1:max(Blk)), Area = AreaB, 
                    Kdns = Kdns, Ktot = Ktot)
  result <- list(Ktab=Ktab,Habdns=Habdns)
  return(result)
}
# 
# Dispersal probabilities
disp = matrix(data = NA,nrow = 4, ncol = P)
disp[1,] = DispP$Jf
disp[2,] = DispP$Af
disp[3,] = DispP$Jm
disp[4,] = DispP$Am
#
# Estimate K for each Block: 
# NOTE: this will eventually be a function of habitat parameters
tmp = Kcalc(PUID,Blk,area,dep,botm,fetch,parms,K_mean); Ktab=tmp$Ktab; Habdns=tmp$Habdns
Kmn = Ktab$Kdns; KV = K_sig^2
muK = log(Kmn/sqrt(1+KV/Kmn^2)); sigK = sqrt(log(1+KV/Kmn^2))
Areahab = Ktab$Area
# *Create a default matrix with lambda ~1.22 (rmax)
A = matrix(c(
  0.7176,    0.4608,    0,         0,
  0.2524,    0.9800,    0,         0,
  0,         0.4608,    0.7132,    0,
  0,         0,         0.2468,    0.9700),byrow=T,ncol=4)    
W=eigen(A)$vector;          # W=matrix of right eigenvectors 
lambdas=eigen(A)$values;    # lambdas=vector of eigenvalues
lambda=max(Re(lambdas))		# lambda1=dominant eigenvalue, real part only
w=abs(W[,1])					      # w=stable distribution, unscaled
sad = w/sum(w)                # w=stable distribution, scaled
#
#~~~~~~~~~~~~~~ End data processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Run Simulations -----------------------------------------------------------
# Initialize arrays to store results:
N = array(data = 0, c(P,Nyrs,reps))
Nmn = matrix(0,nrow=P,ncol=Nyrs)
Nmn[,1] = N0
zvec = matrix(data = 0,nrow = 4, ncol = 1)

# Cycle through reps
for (r in 1:reps){
  # Number of years of population establishment (before range expansion begins)
  YrsInit = round(runif(1,2,10))
  # Determine wavespeed
  V = runif(1,V_mn,V_mx)
  alpha = 1/V
  # Reinitialize Block occupation matrix (keeps track of which blocks occupied each year)
  BlokOcc = matrix(data = 0,nrow = P, ncol = Nyrs)
  # Re-initialize K vector
  K = numeric(length = P)
  # Re-Initialize population vector at year 1 (using N0) for each rep
  n = array(data = 0, c(4,P,Nyrs))
  nt = matrix(data = 0,nrow = 4, ncol = 1) 
  nd = matrix(data = 0,nrow = 4, ncol = 1)
  # If uncertainty in params, include here:
  # tmp = Kcalc(PUID,Blk,area,dep,botm,fetch,parms,KDmean); Ktab=tmp$Ktab; Kmn = Ktab$Kdns
  # muK = log(Kmn/sqrt(1+KV/Kmn^2)); sigK = sqrt(log(1+KV/Kmn^2))
  K = numeric(length = P)
  for (i in 1:P){
    K[i] = rlnorm(1,muK[i],sigK[i])*Areahab[i]
    # n[1:4, i, 1] = rmultinom(1, N0[i], sad) 
    if(N0[i]>0){
      n[1:4, i, 1] =c(0,2,0,0)+ rmultinom(1, N0[i]-2, c(.1,0,.4,.5))
      BlokOcc[i,1] = 1
    }else{
      n[1:4, i, 1] = c(0,0,0,0)
    }
  }
  N[,1, r] = colSums(n[, , 1])
  for (y in 2:Nyrs) {
    #  First, determine if any new blocks occupied (allows for range expansion)
    oc = which(BlokOcc[,y-1]==1)
    BlokOcc[oc,y] = 1      
    # Find all unoccupied cells
    noc = which(BlokOcc[,y-1]==0)
    # S_noc = data$Block[noc]
    # Loop through unoccupied cells, see if they could be colonized
    if(sum(noc)>0){
      for(k in 1:length(noc)){
        j = noc[k]
        nAd = data$NAdj[j]
        # for each adjacent cell, see if its duration of occupation
        # is greater than alpha*Distance between centroids (ie would expect 
        #  range expansion into this new habitat)
        for (a in 1:nAd){
          Adblk = data[j,5+a] # Assumes NAdj is in collumn 5 
          if (sum(BlokOcc[Adblk,1:(y-1)])>0){
            dst = Distmat[Adblk,j]  
            # Probability of colonization: logit fxn of years occupied relative to 
            #   movement of population front given assymptotic wave speed
            # NOTE: assymtotic wave speed only obtained once exp growth is occuring
            probexp = inv.logit(2*(sum(BlokOcc[Adblk,1:(y-1)])-alpha*dst))*max(0,sign(y-YrsInit))
            BlokOcc[j,y] = max(BlokOcc[j,y],rbinom(1,1,probexp))
          }
        }
      }
    }
    # Next, step through blocks and compute dynamics for all occupied blocks
    for (i in 1:P){
      if (BlokOcc[i, y] == 1){
        # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
        lamstoch = max(.95,min(1.22, round(exp(rmax*(1-(N[i,y-1,r]/K[i])^theta)+rnorm(1,0,sig)),2)))
        # NOTE: early surveys indicate slow growth over first few years (allee effect):
        idxs = which(round(Demdat$Lam,2) == lamstoch)
        j = sample(idxs,1)
        br = Demdat$br2[j]; wr = Demdat$wr2[j];
        fsj = Demdat$fs1[j]; fsa = Demdat$fs2[j]; 
        msj = Demdat$ms1[j]; msa = Demdat$ms2[j];
        gf = Demdat$Gf[j]; gm = Demdat$Gm[j];
        nt[1:4,] = n[1:4,i,y-1]  
        # Account for demographic stochasticity
        if(sum(nt)<25 & nt[2]>0){ 
          juvs = round(nt[2]*br*wr*fsa)
          Fjuvs = rbinom(1,juvs,.5)
          Mjuvs = juvs-Fjuvs
          FF= Fjuvs/nt[2]
          FM= Mjuvs/nt[2]  
        }else{
          FF=br*wr*fsa
          FM=br*wr*fsa 
        }
        # Construct matrix               
        AP = matrix(c(
          fsj*(1-gf),  FF,      0,           0,    
          fsj*gf,      fsa,     0,           0,
          0,           FM,      msj*(1-gm),  0,
          0,           0,       msj*gm,      msa),nrow=4,ncol=4,byrow = T)
        #
        # Next lines do matrix multiplication (within-block demog transitions)
        nt1 = round(AP%*%nt)
        # ***NOTE: NEXT LINES ACCOUNT FOR IMMIGRATION 
        #  (Assumes outside immigrants arrive at initially occupied block(s) only)
        if (ImigratOpt > 0) {
          if (ImigratOpt == 1 & BlokOcc[i, 1] == 1) {
            NImm = round(rnbinom(1,Dispers,pparLo))
            ni = rmultinom(1, NImm, c(.1,.05,.4,.45))
          }else if (ImigratOpt == 2 & BlokOcc[i, 1] == 1){
            NImm = round(rnbinom(1,Dispers,pparHi))
            ni = rmultinom(1, NImm, c(.1,.05,.4,.45))
          }else{
            ni = c(0, 0, 0, 0)
          }
        } else {
          ni = c(0, 0, 0, 0)
        }
        nt1 = pmax(zvec,nt1 + ni)
        # Next, Calculate number of dispersers (with stochasticity)
        # ****NOTE: no dispersal happens until pop established, and no more than
        # 1/2 of block residsents can disperse in one year
        if (y > YrsInit & sum(BlokOcc[,y])>1){
          nd[1] = min(floor(nt1[1]/2),rpois(1,nt1[1]*disp[1,i]))
          nd[2] = min(floor(nt1[2]/2),rpois(1,nt1[2]*disp[2,i]))
          nd[3] = min(floor(nt1[3]/2),rpois(1,nt1[3]*disp[3,i]))
          nd[4] = min(floor(nt1[4]/2),rpois(1,nt1[4]*disp[4,i]))
        } else {
          nd[1:4] = zvec
        }
        n[1:4,i,y] = n[1:4,i,y] + nt1-nd
        #
        # Now distribute dispersers randomly among "currently occupied" blocks 
        # with probabilities determined appropriately for each age/sex class
        if (sum(nd) > 0){  
          # JF dispersal
          JFprob = BlokOcc[,y]*destJF[,i] 
          JFprob = JFprob/sum(JFprob)
          ndJF = rmultinom(1, nd[1], JFprob)
          n[1,,y] = n[1,,y] + t(ndJF)
          # AF dispersal
          AFprob = BlokOcc[,y]*destAF[,i] 
          AFprob = AFprob/sum(AFprob)
          ndAF = rmultinom(1, nd[2], AFprob)
          n[2,,y] = n[2,,y] + t(ndAF)     
          # JM dispersal
          JMprob = BlokOcc[,y]*destJM[,i] 
          JMprob = JMprob/sum(JMprob)
          ndJM = rmultinom(1, nd[3], JMprob)
          n[3,,y] = n[3,,y] + t(ndJM)
          # AM dispersal
          AMprob = BlokOcc[,y]*destAM[,i] 
          AMprob = AMprob/sum(AMprob)
          ndAM = rmultinom(1, nd[4], AMprob)
          n[4,,y] = n[4,,y] + t(ndAM)            
        }
      }
    }
    # Tabulate sub-population abundance in each block 
    N[,y,r] = colSums(n[,,y])
    Nmn[,y] = Nmn[,y] + N[,y,r]/reps
  }
}
Dmn = Nmn
for (y in 1:Nyrs){
  Dmn[,y] = Nmn[,y]/Areahab
}
#  Do some plots ---------------------------------------------------------
# Heatmap of Density vs Coastal Block
dfDens = data.frame(Blocks = (data$BlockID),Dmn)
colnames(dfDens) = c('Block',as.character(Years))
df_Dens <- melt(dfDens, id.vars = "Block")
names(df_Dens)[2:3] <- c("Year", "Density")
dfDens = df_Dens[with(df_Dens,order(Block,Year)),]; rm(df_Dens)
dfDens$Block = as.factor(dfDens$Block)

plt1 = ggplot(dfDens, aes(Year, Block)) +
  geom_tile(aes(fill = Density), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  xlab("Year in Future") +
  ylab("Coastal Block #") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Mean Expected Density",title="Projected Density by Block (otters/km2)")
print(plt1)

# Trend plot of abundance over time
# Calculate mean trend and CI (use bootstrap CI, 1000 reps)
# First, entire Metapopulation (all of SEAK):
mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)
CIL.fun <- function(dat, idx) quantile(dat[idx], 0.05, na.rm = TRUE)
CIH.fun <- function(dat, idx) quantile(dat[idx], 0.95, na.rm = TRUE)
means = numeric(length=Nyrs)
SEmean = numeric(length=Nyrs)
Lo = numeric(length=Nyrs)
Hi = numeric(length=Nyrs)
CImn = matrix(0,nrow = Nyrs,ncol=2)
means[1] = Nmn[1,1]
Lo[1] = Nmn[1,1]
Hi[1] = Nmn[1,1]
CImn[1,1:2] = Nmn[1,1]
for(y in 2:Nyrs){
  Nsum = colSums(N[,y,])
  bootobj = boot(Nsum, mean.fun, R=1000, sim="ordinary")
  means[y] = median(bootobj$t)
  SEmean[y] = sd(bootobj$t)
  tmp = boot.ci(bootobj, type="bca"); CImn[y,] = tmp$bca[4:5]
  bootobj = boot(Nsum, CIL.fun, R=100, sim="ordinary")
  Lo[y] = median(bootobj$t)
  bootobj = boot(Nsum, CIH.fun, R=100, sim="ordinary")
  Hi[y] = median(bootobj$t)  
}
Pop_Overall <- data.frame(Year=Years,Mean=means,lower=Lo,upper=Hi,
                          SEmean=SEmean,CImeanLo=CImn[,1],CImeanHi=CImn[,2])
write.csv(Pop_Overall,'Results_GHtot.csv',row.names = FALSE)

titletxt = paste0("Sea Otter Population Projection, ", Nyrs," Years")
plt2 = (ggplot(Pop_Overall, aes(Year, Mean))+
         geom_line(data=Pop_Overall)+
         geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
         geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
         xlab("Year") +
         ylab("Expected Abundance") +
         ggtitle(titletxt, subtitle="All of Gwaii Haanas"))
print(plt2)

# Summary For each block 
randsmp = sample(seq(1,reps),1000,replace = TRUE)
meansALL = numeric(length=P*Nyrs)
LoALL = numeric(length=P*Nyrs)
HiALL = numeric(length=P*Nyrs)
YearsALL = numeric(length=P*Nyrs)
BlockIDs = character(length=P*Nyrs)
for (p in 1:P){
  tmp = matrix(nrow=1000,ncol=Nyrs)
  for(r in 1:1000){
    tmp[r,] <- (N[p,,randsmp[r]])
  }
  Lo = numeric(length=Nyrs)
  Hi = numeric(length=Nyrs)
  means = numeric(length=Nyrs)
  means[1] = mean(tmp[,1])
  Lo[1] = mean(tmp[,1])
  Hi[1] = mean(tmp[,1])
  for(y in 2:Nyrs){
    dens = density(tmp[,y],adjust = 3)
    CI = quantile(dens, probs=c(.05, .95))
    # Lo[y] <- max(0,as.numeric(CI[1]))
    Lo[y] <- quantile(tmp[,y], probs=c(.05))
    Hi[y] <- round(as.numeric(CI[2]))
    means[y] <- mean(tmp[,y])
  }
  meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
  LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
  HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
  YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
  BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Block[p],Nyrs))
}
Hab_Blocks <- data.frame(Block = BlockIDs, 
                         Year=YearsALL,Mean=meansALL, 
                         lower=LoALL,upper=HiALL,Density=dfDens$Density )
write.csv(Hab_Blocks,'Results_GHblocks.csv',row.names = FALSE)

Habdns$DensT = 0
BlkDnsT = Hab_Blocks$Density[Hab_Blocks$Year==max(Years)]
for (i in 1:length(Habdns$PUID)){
  Habdns$DensT[i] = Habdns$Reldens[i]*sum(BlkDnsT*as.numeric(Habavg[i,2:(P+1)]))
}
write.csv(Habdns,'Results_GHcells.csv',row.names = FALSE)
Cdata$DensT = Habdns$DensT

# plot map of Results ---------------------------------------------------------
GHland<-readOGR("GH_Land_poly.shp", layer="GH_Land_poly")
GHland_df <- fortify(GHland)
map <- ggplot() +
  geom_polygon(data = GHland_df, 
            aes(x = long, y = lat, group = group)) +
  
  # geom_point(aes(x=Xcoord, y=Ycoord, color=DensT), data=Cdata, alpha=1, size=1, color="grey20") +
  geom_point(aes(x=Xcoord, y=Ycoord, color=DensT), data=Cdata, alpha=1, size=1.5, shape=15)+
  scale_colour_gradientn("Mean Density, Final Year", 
                         colours=c( "#f9f3c2","#660000"))+ # change color scale
  coord_equal(ratio=1) # square plot to avoid the distortion
print(map) 


