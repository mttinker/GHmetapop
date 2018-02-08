# Gwaai Haanas sea otter meta-population model
library(gtools)
#
# Load files --------------------------------------------------------------
data = read.csv("Blockdata.csv", header = TRUE)
Demdat = read.csv("RandDem.csv", header = TRUE)
SEmoves = read.csv("SE_Ottermoves.csv", header = TRUE)
Distmat =  read.csv("Distmat.csv", header = FALSE)
#
# Inter-pop movement matrives: pairwise prob of dispersal based on 
#  pairwise distances and dispersal kernels for each age/sex class
# destJF = read.csv("destJF.csv", header = FALSE);
# destAF = read.csv("destAF.csv", header = FALSE);
# destJM = read.csv("destJM.csv", header = FALSE);
# destAM = read.csv("destAM.csv", header = FALSE);
#
# Set User Parameters  ---------------------------------------------------------
reps = 25            # Number replications for population sims
Nyrs = 10            # Number of years to project population dynamics
ImigratOpt = 0       # Immigration option: 0 = none, 1 = low, 2 = high
V_mn = 3.5           # Population front asymptotic wavespeed, mean  
V_sd = .75           # Population front asymptotic wavespeed, mean
Nstg = 4             # Number of age/sex classes 
sig = 0.05           # Environmental stochasticity (std dev in log-lambda)
vals = dim(Distmat)  # determine number of blocks and years
P = dim(Distmat)[1]  # number blocks (or sub-populations)
rmax = log(1.22)     # Maximin rate of growth = 23% per year
theta = 0.75         # theta parameter for theta-logistic (value of 1 = basic Ricker model) 
Yr1 = 2018           # Calenadar Year to begin simulations at
# ~~~~~~END User parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Process data ----------------------------------------------------------------
Years = c(Yr1:(Yr1+Nyrs-1))  
Yrs = seq(1:Nyrs)
Years = Yrs-1+Yr1
N0 = data$Startpop
S = 4 # Number of age/sex classes 
# Dispersal probabilities
disp = matrix(data = NA,nrow = 4, ncol = P)
disp[1,] = data$Jfdisp
disp[2,] = data$Afdisp
disp[3,] = data$Jmdisp
disp[4,] = data$Amdisp
#
# Estimate K for each Block: 
# NOTE: this will eventually be a function of habitat parameters
Kmn = numeric(length=P)
Ksd = numeric(length=P)
KV = numeric(length=P)
Areahab = numeric(length=P)
for(i in 1:P){
  Areahab[i] = data$Area[i] 
  Kmn[i] = data$Kdens[i] 
  Ksd[i] = 1
  KV[i] = Ksd[i]^2
}
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
  YrsInit = round(runif(1,2,5))
  # Determine wavespeed
  V = max(0.5,rnorm(1,V_mn,V_sd))
  alpha = 1/V
  # Reinitialize Block occupation matrix (keeps track of which blocks occupied each year)
  BlokOcc = matrix(data = 0,nrow = P, ncol = Nyrs)
  # Re-initialize K vector
  K = numeric(length = P)
  # Re-Initialize population vector at year 1 (using N0) for each rep
  n = array(data = 0, c(4,P,Nyrs))
  nt = matrix(data = 0,nrow = 4, ncol = 1) 
  nd = matrix(data = 0,nrow = 4, ncol = 1)
  for (i in 1:P){
    # n[1:4, i, 1] = rmultinom(1, N0[i], sad) 
    if(N0[i]>0){
      n[1:4, i, 1] =c(0,2,0,0)+ rmultinom(1, N0[i]-2, c(.1,0,.4,.5))
      BlokOcc[i,1] = 1
    }else{
      n[1:4, i, 1] = c(0,0,0,0)
    }
    KDmean = rlnorm(1,log(Kmn[i]/sqrt(1+KV[i]/Kmn[i]^2)),sqrt(log(1+KV[i]/Kmn[i]^2)))
    K[i] = KDmean*Areahab[i]
  }
  rnd = rbinom(Nyrs, 1, .5) # random coin-toss for rounding (demog stochasticity)
  N[,1, r] = colSums(n[, , 1])
  for (y in 2:Nyrs) {
    #  First, determine if any new blocks occupied (allows for range expansion)
    oc = which(BlokOcc[,y-1]==1)
    BlokOcc[oc,y] = 1      
    # Find all unoccupied cells
    noc = which(BlokOcc[,y-1]==0)
    # S_noc = data$Block[noc]
    # Loop through unoccupied cells, see if they could be colonized
    for(k in 1:length(noc)){
       j = noc[k]
       nAd = data$Nadj[j]
       # for each adjacent cell, see if its duration of occupation
       # is greater than alpha*Distance between centroids (ie would expect 
       #  range expansion into this new habitat)
      for (a in 1:nAd){
        Adblk = data[j,10+a]
        dst = Distmat[Adblk,j]  
        # Probability of colonization: logit fxn of years occupied relative to 
        #   movement of population front given assymptotic wave speed
        probexp = inv.logit(2*(sum(BlokOcc[Adblk,1:y])-alpha*dst))*max(0,sign(y-YrsInit))
        BlokOcc[j,y] = max(BlokOcc[j,y],rbinom(1,1,probexp))
      }
    }
    # Next, step through blocks and compute dynamics for all occupied blocks
    for (i in 1:P){
      if (BlokOcc[i, y] == 1){
        # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
        lamstoch = max(.95,min(1.22, round(exp(rmax*(1-(N[i,y-1,r]/K[i])^theta)+rnorm(1,0,sig)),2)))
        # NOTE: early surveys indicate slow growth over first few years (allee effect):
        if(y<=YrsInit){
          lamstoch = max(1,round(lamstoch*0.95,2))
        }
        idxs = which(round(Demdat$Lam,2) == lamstoch)
        j = sample(idxs,1)
        br = Demdat$br2[j]; wr = Demdat$wr2[j];
        fsj = Demdat$fs1[j]; fsa = Demdat$fs2[j]; 
        msj = Demdat$ms1[j]; msa = Demdat$ms2[j];
        gf = Demdat$Gf[j]; gm = Demdat$Gm[j];
        AP = matrix(c(
          fsj*(1-gf),  br*wr*fsa,      0,           0,    
          fsj*gf,      fsa,            0,           0,
          0,           br*wr*fsa,      msj*(1-gm),  0,
          0,           0,              msj*gm,      msa),nrow=4,ncol=4,byrow = T)
        #
        # Next lines do matrix multiplication (within-block demog transitions)
        nt[1:4,] = n[1:4,i,y-1]
        if (rnd[y]==1){
          nt1 = floor(AP%*%nt)   # randomly round up or down to integers
        } else  {
          nt1 = round(AP%*%nt)
        }
        # ***NOTE: NEXT LINES ACCOUNT FOR IMMIGRATION 
        #  (Assumes outside immigrants arrive at initially occupied block(s) only)
        if (ImigratOpt > 0) {
          if (ImigratOpt == 1 & BlokOcc[i, 1] == 1) {
            NImm = round(rgamma(1,.75,.5))
            ni = rmultinom(1, NImm, c(.1,.05,.4,.45))
          }else if (ImigratOpt == 2 & BlokOcc[i, 1] == 1){
            NImm = round(rgamma(1,1.5,.5))
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
        # 1/3 of block residsents can disperse in one year
        if (y > YrsInit){
          nd[1] = min(floor(nt1[1]/3),rpois(1,nt1[1]*disp[1,i]))
          nd[2] = min(floor(nt1[2]/3),rpois(1,nt1[2]*disp[2,i]))
          nd[3] = min(floor(nt1[3]/3),rpois(1,nt1[3]*disp[3,i]))
          nd[4] = min(floor(nt1[4]/3),rpois(1,nt1[4]*disp[4,i]))
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
# Calculate mean trend and CI (use bootstrap CI, 1000 reps)
# First, entire Metapopulation (all of SEAK):
randsmp = sample(seq(1:reps),1000, replace=TRUE)
tmp = matrix(nrow=1000,ncol=Nyrs)
for(r in 1:1000){
  tmp[r,] <- colSums(N[,,randsmp[r]])
}
RepSums <- data.frame(tmp)
Lo = numeric(length=Nyrs)
Hi = numeric(length=Nyrs)
means = numeric(length=Nyrs)
means[1] <- mean(RepSums[,1])
for(y in 2:Nyrs){
  dens = density(RepSums[,y],adjust = 7)
  CI = quantile(dens, probs=c(.025, .975))
  Lo[y] <- max(0,as.numeric(CI[1]))
  Hi[y] <- as.numeric(CI[2])
  means[y] <- mean(RepSums[,y])
}
Pop_Overall <- data.frame(Year=Years,Mean =means, lower=Lo,upper=Hi)

# For each block (the BIG block table)
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
  means[1] <- mean(tmp[,1])
  for(y in 2:Nyrs){
    dens = density(tmp[,y],adjust = 7)
    CI = quantile(dens, probs=c(.025, .975))
    Lo[y] <- max(0,as.numeric(CI[1]))
    Hi[y] <- as.numeric(CI[2])
    means[y] <- mean(tmp[,y])
  }
  meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
  LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
  HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
  YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
  BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Block[p],Nyrs))
}
Hab_Blocks <- data.frame(Block = BlockIDs, 
                         Year=YearsALL,Mean =meansALL, 
                         lower=LoALL,upper=HiALL)

library(ggplot2)
plt = (ggplot(Pop_Overall, aes(Year, Mean))+
         geom_line(data=Pop_Overall)+
         geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.3)+
         xlab("Year") +
         ylab("Expected Population Size") +
         ggtitle(titltxt,subtitle="All of Gwaii Haanas"))
print(plt)