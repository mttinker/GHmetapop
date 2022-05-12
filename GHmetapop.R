# Gwaai Haanas sea otter meta-population model
# EXPLANATION
# Model simulates growth of a recently colonized sea otter population at Gwaii Haanas.
# Dynamics modeled for P coastal blocks (sub-populations) linked by dispersal, with
# results of simulations then downscaled to 1km grid of coastal habitat cells. 
# Stage-structed matrix model includes density-dependence, environmental stochasticity,
# demographic stochasticity, immigration and dispersal, range expansion via diffusion,
# and habitat-based variation in local equilibrium densities.
# User provides model parameters, which are informed by analyses of data from
# sea otter populations elsewhere in BC and SE Alaska
# 
rm(list = ls())
# Set User Parameters  ---------------------------------------------------------
reps = 500           # Number replications for population sims (should use at least 100)
Nyrs = 25           # Number of years to project population dynamics
ImmRt = 0          # Immigration rate (avg. net new immigtants per year): 0 = none
V_sp = 3             # Population front asymptotic wavespeed, km/yr, minimum  
Emax = 5             # Maximum years before pop "established" (and before range expansion begins)
K_mean = 3.5         # Overall mean K density (modified as fxn of habitat variables)
K_sig = 1            # Stochasticity in K (std. deviation in K density)
sig = 0.05           # Environmental stochasticity (std. deviation in log-lambda)
rmax = 0.2           # Maximum rate of growth: default = log(1.22), or 22% per year
theta = 1            # theta parameter for theta-logistic (1 = Ricker model, >1 = delayed DD) 
Initpop = 20         # Number of animals in initial population (at least 2 adult females)
Initblk = c(1)# List of initially occupied bloaks: e.g. c(1) = Block 1 only
# ~~~~~~END User parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Load necessary libraries -------------------------------------------------
# NOTE: Ensure following packages are installed
library(gtools)
library(mvtnorm)
library(boot)
library(ggplot2)
library(rgdal)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(reshape2)
library(doParallel)
library(abind)
forloop_packages<-c('gtools','mvtnorm')
#
# Load files ----------------------------------------------------------------
data = read.csv("./data/GHBlockdata.csv", header = TRUE)  # Data for coastal blocks
Cdata = read.csv("./data/GHCelldata.csv", header = TRUE)  # Data for habitat cells
Demdat = read.csv("./data/RandDem.csv", header = TRUE)    # Stochastic vital rates
Distmat =  as.matrix(read.csv("./data/Distmat.csv", header = FALSE)) # Inter-blk LCP distances
# Probabilities of dispersal from each block for each age/sex class
DispP = read.csv("./data/GHDispProb.csv", header = TRUE) 
# Inter-pop movement matrices: pairwise prob of dispersal based on 
#  pairwise LCP distances and dispersal kernels for each age/sex class
destJF = read.csv("./data/GHDispMatJF.csv", header = FALSE) 
destAF = read.csv("./data/GHDispMatAF.csv", header = FALSE);
destJM = read.csv("./data/GHDispMatJM.csv", header = FALSE);
destAM = read.csv("./data/GHDispMatAM.csv", header = FALSE);
# Load parameters for habitat density at K function
params = read.csv("./data/Hab_params.csv")
# Load GIS data for plotting results (NAD_1983_Albers)
load("./data/GISdata.rdata")
#
# Process data ---------------------------------------------------------------
Yr1 = as.numeric(format(Sys.Date(), "%Y"))
V_mn = max(.5,V_sp-1)
V_mx = min(8,V_sp+1)
# Inverse-scale param: over-disperseion for Neg Binomial # immigrants per year
InvScale = 1
ppar = InvScale/(InvScale+ImmRt)
KV = K_sig^2  # Variance in K density
rmax = min(log(1.22),rmax)
Years = c(Yr1:(Yr1+Nyrs-1))  
Yrs = seq(1:Nyrs)
Years = Yrs-1+Yr1
P = nrow(Distmat)  # number blocks (or sub-populations)
# Initialize population vector
initblks = rep(0,P)
N0 = rep(0,P)
initblks[Initblk] = 1
N0[Initblk] = 2
N0 = N0 + rmultinom(1,max(1,Initpop-2*length(Initblk)), initblks)
#  Create Data variables for calculating relative density of Hab cells:
Cdata$PU_ID = Cdata$Cell_ID
PUID = Cdata$PU_ID
Blk = Cdata$BlockID
area = Cdata$Area
dep = Cdata$Depth
fetch = Cdata$Fetch
c1 = which(colnames(Cdata)=="BoP_1")
cF = which(colnames(Cdata)=="BoP_3b")
botm = as.matrix(Cdata[,c1:cF])
rm(c1,cF)
egrass = Cdata$Eelgrass
parms = params$Parms
Nparms = length(parms)
# Define the "Kcalc" function, estimates local K density based on habitat variables
Kcalc <- function(PUID,Blk,area,dep,botm,fetch,egrass,parms,Kmn){
  # NOTE: Hab dens multiplier fxn: exp(b1*X1 + b2*X2... + bn*Xn)
  #   Multiplier is used to adjust local K density for each cell 
  b <- parms # Create parameter vector  
  # Depth part of fxn: b1*(-1*dep) - b2*dep^2   (depth is in negative values)
  #   where sum(area*exp(b1*(-1*dep)-b2*(dep^2)))/sum(area) =~ 1
  Depfxn = b[1] + b[2]*(-1*dep) - b[3]*dep^2
  # Bottom patch part of function
  Btfxn = b[4]*botm[,1] + b[5]*botm[,2] + b[6]*botm[,3] +
          b[7]*botm[,4] + b[8]*botm[,5] + b[9]*botm[,6] +
          b[10]*botm[,7] + b[11]*botm[,8] + b[12]*botm[,9]
  #  Btfxn = numeric(length = nrow(botm))
  #  Btfxn[which(botm=="1")] = b[4]; Btfxn[which(botm=="1a")] = b[5]; Btfxn[which(botm=="1b")] = b[6]
  #  Btfxn[which(botm=="2")] = b[7]; Btfxn[which(botm=="2a")] = b[8]; Btfxn[which(botm=="2b")] = b[9]
  #  Btfxn[which(botm=="3")] = b[10]; Btfxn[which(botm=="3a")] = b[11]; Btfxn[which(botm=="3b")] = b[12]
  # Eelgrass part of function
  EGfxn =  b[13]*(egrass) 
  # Fetch part of fxn
  Fchfxn = b[14]*(fetch - 25) 
  # Combine to terms to create multiplier for each hab cell
  mult = exp(Depfxn + Btfxn + EGfxn + Fchfxn) 
  # NOTE: sum((mult*area))/sum(area) should equal approx 1 (to maintain overall K_mean)
  Kdns = numeric(length = max(Blk))
  Ktot = numeric(length = max(Blk))
  AreaB = numeric(length = max(Blk))
  # Loop thru coastal blocls to estimate total K and K_density
  for (i in 1:max(Blk)){
    ii = which(Blk == i) # Select all hab cells in this Block
    # Total K for block, Ktot: 
    #  summed product of cell density (adjusted by multiplier) and cell area
    Ktot[i] = sum(Kmn*mult[ii]*area[ii]) 
    AreaB[i] = sum(area[ii]) # Total area of all hab cells in Block
    Kdns[i] = Ktot[i]/sum(area[ii]) # Mean density for block
  }
  Habdns = data.frame(PUID = PUID, Reldens = mult)
  Ktab = data.frame(Block = seq(1:max(Blk)), Area = AreaB, 
                    Kdns = Kdns, Ktot = Ktot)
  result <- list(Ktab=Ktab,Habdns=Habdns)
  return(result)
}
# Estimate mean K for each Block, and mu/sig for log-normal sampling of K: 
tmp = Kcalc(PUID,Blk,area,dep,botm,fetch,egrass,parms,K_mean); 
Ktab=tmp$Ktab; Habdns=tmp$Habdns; 
Areahab = Ktab$Area
cellscale = tmp$Habdns$Reldens
# 
# Environmental stochasticity: calc SIGMA for correlated random effects
SIGMA = as.matrix((sig^2)*exp(-.005*Distmat))
MU_Z = rep(0,P)
zvec = matrix(data = 0,nrow = 4, ncol = 1)
#
# Dispersal probabilities
disp = matrix(data = NA,nrow = 4, ncol = P)
disp[1,] = DispP$Jf
disp[2,] = DispP$Af
disp[3,] = DispP$Jm
disp[4,] = DispP$Am
#
# *Create a default matrix with lambda ~1.22 (rmax)
A = matrix(c(
  0.5307,    0.4183,         0,         0,
  0.4093,    0.9700,         0,         0,
  0,         0.4183,    0.5110,         0,
  0,         0,         0.3690,    0.9200),byrow=T,ncol=4)    
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
# N = array(data = 0, c(P,Nyrs,reps))
# Nmn = matrix(0,nrow=P,ncol=Nyrs)
# Nmn[,1] = N0
# Set up for parallel processing
acomb <- function(...) abind(..., along=3)
ncores = min(20,detectCores()-4)
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)
# Cycle through reps
Nreps <- foreach(r=1:reps, .combine='acomb', .multicombine=TRUE,
                 .packages=forloop_packages) %dopar% {
  N = matrix(0,nrow = P,ncol = Nyrs)
  # Number of years of population establishment (before range expansion begins)
  E = round(runif(1,round(Emax/2),Emax))
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
  parms = runif(Nparms,params$ParmLo,params$ParmHi)
  # *** NOTE: If including uncertainty in params, uncomment next 2 lines:
  tmp = Kcalc(PUID,Blk,area,dep,botm,fetch,egrass,parms,K_mean); Kmn =tmp$Ktab$Kdns
  muK = log(Kmn/sqrt(1+KV/Kmn^2)); sigK = sqrt(log(1+KV/Kmn^2))
  K = numeric(length = P)
  for (i in 1:P){
    K[i] = rlnorm(1,muK[i],sigK[i])*Areahab[i]
    # n[1:4, i, 1] = rmultinom(1, N0[i], sad) 
    if(N0[i]>0 & length(Initblk)==1){
      n[1:4, i, 1] =c(0,2,0,0)+ rmultinom(1, max(0,N0[i]-2), sad)
      BlokOcc[i,1] = 1
    }else if(N0[i]>0 & length(Initblk)>1) {  
      n[1:4, i, 1] =c(0,1,0,0)+ rmultinom(1, max(0,N0[i]-1), sad)
      BlokOcc[i,1] = 1
    }else{
      n[1:4, i, 1] = c(0,0,0,0)
    }
  }
  N[,1] = colSums(n[, , 1])
  # env. stoch., correlated random annual deviations 
  eps = rmvnorm(Nyrs, MU_Z, SIGMA)
  for (y in 2:Nyrs) {
    #  First, determine if any new blocks occupied (allows for range expansion)
    oc = which(BlokOcc[,y-1]==1)
    BlokOcc[oc,y] = 1      
    # Find all unoccupied cells
    noc = which(BlokOcc[,y-1]==0)
    # S_noc = data$Block[noc]
    # Loop through unoccupied cells, see if they could be colonized
    if(sum(noc) > 0 & y > E){
      for(k in 1:length(noc)){
        j = noc[k]
        iib = which(Distmat[,j]<50 & Distmat[,j]>0 & BlokOcc[,(y-1)]==1)
        nAd = length(iib)
        # for each neighbouring block, see if its duration of occupation
        # is greater than alpha*Distance between centroids (ie would expect 
        #  range expansion into this new habitat)
        if(nAd==0){
          BlokOcc[j,y] = 0
        }else{
          for (a in 1:nAd){
            Adblk = iib[a]
            if (sum(BlokOcc[Adblk,1:(y-1)])>0){
              dst = Distmat[Adblk,j]  
              # Probability of colonization: logit fxn of years occupied relative to 
              #   movement of population front given assymptotic wave speed
              # NOTE: assymtotic wave speed only obtained after establishment phase
              probexp = inv.logit(2*(sum(BlokOcc[Adblk,1:(y-1)])-alpha*dst)+1.5)
              BlokOcc[j,y] = max(BlokOcc[j,y],rbinom(1,1,probexp))
            }
          }
        }
      }
    }
    # Calculate annual immigration and distribute randomly among occupied blocks
    if(ImmRt==0){
      NImm = rep(0,P)
    }else{
      NImm = rnbinom(1,InvScale,ppar)
      NImm = rmultinom(1,NImm, BlokOcc[,y])
    }
    # Next, step through blocks and compute dynamics for all occupied blocks
    for (i in 1:P){
      if (BlokOcc[i, y] == 1){
        # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
        lamstoch = max(.95,min(1.22, round(exp(rmax*(1-(N[i,y-1]/K[i])^theta)+eps[y,i]),2)))
        # NOTE: account for "population establishment" phase (allee effect):
        if (y <= E){
          lamstoch = max(.95,min(1.22, round(lamstoch^0.2,2)))
        } 
        idxs = which(round(Demdat$Lam,2)==lamstoch)
        j = sample(idxs,1)
        br = Demdat$br2[j]; wr = Demdat$wr2[j];
        fsj = Demdat$fs1[j]; fsa = Demdat$fs2[j]; 
        msj = Demdat$ms1[j]; msa = Demdat$ms2[j];
        gf = Demdat$Gf[j]; gm = Demdat$Gm[j];
        nt[1:4,1] = n[1:4,i,y-1]  
        FF=(br/2)*wr*fsa
        FM=FF         
        # Construct matrix               
        AP = matrix(c(
          fsj*(1-gf),  FF,      0,           0,    
          fsj*gf,      fsa,     0,           0,
          0,           FM,      msj*(1-gm),  0,
          0,           0,       msj*gm,      msa),nrow=4,ncol=4,byrow = T)
        #
        # Next lines do matrix multiplication (within-block demog transitions)
        if (sum(nt)>0){
          nt1 = AP%*%nt
          if (sum(nt)>50){
            nt1 = round(nt1)
          }else{
            # randomly round up or down for small pops (demographic stochasticity)
            nt1[1] = sample(c(ceiling,floor),1)[[1]](nt1[1])
            nt1[2] = sample(c(ceiling,floor),1)[[1]](nt1[2])
            nt1[3] = sample(c(ceiling,floor),1)[[1]](nt1[3])
            nt1[4] = sample(c(ceiling,floor),1)[[1]](nt1[4])
          }
        }else{
          nt1 = zvec
        }
        # NEXT LINES ACCOUNT FOR IMMIGRATION or RE-INTRODUCTION 
        #  (randomly assign age/sex class to immigrants)
        if (NImm[i]>0) {
          ni = rmultinom(1, NImm[i], sad)              
        } else {
          ni = zvec
        }
        # Next, Calculate number of dispersers (with stochasticity)
        # ****NOTE: no dispersal happens until pop established, 
        # must be other occupied blocks nearby,
        #  some individuals of each age class must remain resident
        recip_poss = length(which(BlokOcc[,y]*Distmat[,i]>0 & BlokOcc[,y]*Distmat[,i]<150))
        if (y > E & recip_poss>0 & sum(nt1)>0){
          nd[1] = min(floor(.5*nt1[1]),rpois(1,nt1[1]*disp[1,i]))
          nd[2] = min(floor(.33*nt1[2]),rpois(1,nt1[2]*disp[2,i]))
          nd[3] = min(floor(.9*nt1[3]),rpois(1,nt1[3]*disp[3,i]))
          nd[4] = min(floor(.5*nt1[4]),rpois(1,nt1[4]*disp[4,i]))
        } else {
          nd[1:4] = zvec
        }
        n[1:4,i,y] = n[1:4,i,y] + nt1 + ni - nd
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
    N[,y] = colSums(n[,,y])
    # Nmn[,y] = Nmn[,y] + N[,y,r]/reps
  }
  return(N)
}
N = Nreps
Nmn = apply(Nreps,c(1,2),mean)
Dmn = Nmn
for (y in 1:Nyrs){
  Dmn[,y] = Nmn[,y]/Areahab
}
rm(Nreps)
#  Do some plots ---------------------------------------------------------
# Heatmap of Density vs Coastal Block
if(length(Initblk)==1){
  plotlab = paste(c(paste0("Initial population of ",Initpop, " located in section"),Initblk),
                  collapse=" ")
}else{
  plotlab = paste(c(paste0("Initial population of ",Initpop, " divided among sections"),Initblk),
                collapse=" ")
}
dfDens = data.frame(Blocks = (data$BlockID),Dmn)
colnames(dfDens) = c('Block',as.character(Years))
df_Dens <- melt(dfDens, id.vars = "Block")
names(df_Dens)[2:3] <- c("Year", "Density")
dfDens = df_Dens[with(df_Dens,order(Block,Year)),]; rm(df_Dens)
dfDens$Block = as.factor(dfDens$Block)
maxD <- ceiling(100*max(dfDens$Density))/100
plt1 = ggplot(dfDens, aes(Year, Block)) +
  geom_tile(aes(fill = Density), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue",limits=c(0, maxD)) +
  xlab("Year in Future") +
  ylab("Coastal Section #") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14,face="bold"),
        axis.title=element_text(size=12),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Mean Expected Density",
       title=paste0("Projected Density by Section (otters/km2) after ", Nyrs," Years"),
       subtitle=plotlab) +
  ggtitle(paste0("Projected Density by Section (otters/km2) after ", Nyrs," Years"),
          subtitle=plotlab)
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
means[1] = mean(colSums(N[,1,]))
Lo[1] = mean(colSums(N[,1,]))
Hi[1] = mean(colSums(N[,1,]))
CImn[1,1:2] = sum(Nmn[,1])
for(y in 2:Nyrs){
  Nsum = colSums(N[,y,])
  bootobj = boot(Nsum, mean.fun, R=500, sim="ordinary")
  means[y] = median(bootobj$t)
  SEmean[y] = sd(bootobj$t)
  CImn[y,1] = means[y] - 1.96*SEmean[y]
  CImn[y,2] = means[y] + 1.96*SEmean[y]
  # tmp = boot.ci(bootobj, type="bca", conf = 0.90); CImn[y,] = tmp$bca[2:3]
  bootobj = boot(Nsum, CIL.fun, R=500, sim="ordinary")
  Lo[y] = median(bootobj$t)
  bootobj = boot(Nsum, CIH.fun, R=500, sim="ordinary")
  Hi[y] = median(bootobj$t)  
}
Pop_Overall <- data.frame(Year=Years,Mean=means,lower=Lo,upper=Hi,
                        SEmean=SEmean,CImeanLo=CImn[,1],CImeanHi=CImn[,2])
titletxt = paste0("Projected Sea Otter Population, ", Nyrs," Years")
maxN = ceiling(Pop_Overall$upper[Nyrs]/100)*100
plt2 = (ggplot(Pop_Overall, aes(Year, Mean))+
         geom_line(data=Pop_Overall)+
         geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
         geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
         ylim(0,maxN) +  
         xlab("Year") +
         ylab("Expected Abundance") +
         ggtitle(titletxt, subtitle=plotlab)) + theme_classic(base_size = 12)
print(plt2)

# Compute observed rate of range spread:
# NOTES: 
# - if starting from one focal area (e.g. south end of island), then there are 
#  two "mostly independent" coastlines (east and west coast) that range front
#  is moving along, so observed rate of range spread should be ~ 2x V_sp
# - block considered "occupied" when >2 otters present, on average
# - range extent corrected for coastline complexity, to approximate 1-D coast
#
CoastStr = 2*mean(sort(Distmat[,1],decreasing=T)[2:3])
Coastblk = numeric(length = P)
for (i in 1:P){
  Coastblk[i] = mean(sort(Distmat[,i],decreasing=F)[2:3])
}
CoastAll = sum(Coastblk)
CoastCrct = CoastStr/CoastAll
Range = numeric(length = Nyrs) 
for (i in 1:Nyrs){
  ii = which(Nmn[,i]>2)
  tmp = numeric()
  for(j in 1:length(ii)){
    tmp[j] = mean(sort(Distmat[,ii[j]],decreasing=F)[2:3])
  }
  Range[i] = sum(tmp)*CoastCrct
}
df_Rngspr = data.frame(Year = Years[(Emax+1):min(80,Nyrs)],
                       Range_ext = Range[(Emax+1):min(80,Nyrs)])
ftV = lm(Range_ext ~ Year, data=df_Rngspr)
summary(ftV)
plt3 = ggplot(data = df_Rngspr, aes(x=Year,y=Range_ext)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  labs(x="Year",y="Coastal range extent (km)",
       title = paste("Observed range spread over", Nyrs,"Years,", "V_sp = ",V_sp, "km/yr"),
       subtitle = paste("Adj R2 = ",signif(summary(ftV)$adj.r.squared, 3),
                     ", Slope =",signif(ftV$coef[[2]], 3),
                      ", Approx Corrected Wave Spd (if 1 starting area) =", 
                      signif(ftV$coef[[2]]*.5,2))) +
  theme_classic()
print(plt3)

# Output summary tables ---------------------------------------------------
# Simmulation summary (Overall):
# write.csv(Pop_Overall,paste0('Results_GHtot_',savename,'.csv'),row.names = FALSE)
# Calculate summary For each block 
randsmp = sample(seq(1,reps),1000,replace = TRUE)
meansALL = numeric(length=P*Nyrs)
LoALL = numeric(length=P*Nyrs)
HiALL = numeric(length=P*Nyrs)
YearsALL = numeric(length=P*Nyrs)
BlockIDs = character(length=P*Nyrs)
BlockArea = numeric(length=P*Nyrs)
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
    # Optional: load library BMS, and fit density curve
    # dens = as.numeric(density(tmp[,y],adjust = 3))
    # CI =  quantile.density(dens, probs=c(.05, .95))
    # Lo[y] <- max(0,as.numeric(CI[1]))
    # Lo[y] <- max(0,as.numeric(CI[1]))
    Lo[y] <- quantile(tmp[,y], probs=c(.05))
    Hi[y] <- quantile(tmp[,y], probs=c(.95))
  }
  meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
  LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
  HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
  YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
  BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Block[p],Nyrs))
  BlockArea[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = rep(Areahab[p],Nyrs)
}
Hab_Blocks <- data.frame(Block = BlockIDs, Area = BlockArea, 
                         Year=YearsALL, Mean=meansALL, 
                         lower=LoALL,upper=HiALL,Density=dfDens$Density,
                         DensityLO = LoALL/BlockArea, DensityHI = HiALL/BlockArea)
Hab_Blocks_Fin = Hab_Blocks[which(Hab_Blocks$Year==max(Hab_Blocks$Year)),]
# Output block summaries of simulation results 
write.csv(Hab_Blocks,'./results/Results_HGblocksproject.csv',row.names = FALSE)
#
Celldens = numeric(length = nrow(Cdata))
for (i in 1:P){
  ii = which(Cdata$BlockID==i)
  Celldens[ii] = Hab_Blocks_Fin$Density[i]*cellscale[ii]
}
CellDensProject = data.frame(CellID = Cdata$Cell_ID,BlockID=Cdata$BlockID,
                             Cellscale=cellscale, Celldens = Celldens)
write.csv(CellDensProject,'./results/Results_CellDensProject.csv',row.names = FALSE)
#
# MAP OUTPUT:  ------------------------------------------------------
endvals = Hab_Blocks[Hab_Blocks$Year==max(YearsALL),c(1,7)]
HGblk = merge(HGblk, endvals, by.x='BlockID', by.y = 'Block')
ggplot() + 
  geom_polygon(data=HGblk, aes(x = long, y = lat, fill = Density, color=Density, group = group),
               alpha = 1,size = 1) +
  scale_fill_continuous(low = "#fff7ec", high = "#7F0000") + 
  scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_x_continuous(name = "Longitude (degrees)") + # , breaks = NULL, labels = NULL
  scale_y_continuous(name = "Latitude (degrees)") + # , breaks = NULL, labels = NULL
  geom_polygon(data = HGlnd, aes(x=long,y=lat,fill=piece,group=group),
               color="wheat4", fill="cornsilk1",size = 0.1) +   
  # north(HGlnd,location = "topright") +
  # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
  #         transform = FALSE, location = "bottomleft") +
  ggtitle(titletxt) +
  # coord_equal(ratio=1) + 
  coord_map("conic", lat0 = 18, xlim = c(-133.3, -130.8), ylim=c(51.8, 54.3)) +
  theme_minimal()

