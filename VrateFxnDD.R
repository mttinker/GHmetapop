# VrateFxnDD: Script to determine fxns for calculating baseline vital rates with density dependence
#
# EXPLANATION
# This script generates random sets of vital rates for a 2-stage, 2-sex sea otter pop matrix
# by interpolating between a set of rates for a HIGH growth pop (low density pop) and a set
# of rates for a LOW growth pop (high density pop at K), rates from Kodiak I and Amchitka I, 
# respectively (taken from from Monson et al 2000 (Oikos) and Monson 1995 Masters thesis)
# A funcitonal relationship between the interpolation parameter and Lambda is then estimated,
# as well as the relationship between Lambda and the growth transition probability params
# needed to parameterize 2-stage matrix (G = prob of transitioning from sub-ault to adult).
# These functions can then be used to facilitate simulations of population dynamics using 
# Density-dependent matrix projections with environmental stochasticity
# 
# Start by defining function to build and evaluate a sea otter population matrix 
somatrix2st <- function(vrates){
  #SOMATRIX function to generate 4x4 projection matrix for sea otters
  #  from user-supplied vital rates (2 stages for 2 sexes)
  # EXPLANATION:
  # Produces 2-sex matrix for 2 stages (age classes),
  #  Sub-adults: 0.5-3.5 yr olds (3 year stage duration)
  #  Adults: 3.5-19.5 year olds (16 year stage duration)
  # INPUT: vector of 8 elements containing vital rates:
  #  br = birth rates for stage 1 & 2
  #  wr = weaning success rates for stage 1 & 2
  #  fs = female survival rates for stage 1 & 2
  #  ms = male survival rates for stage 1 & 2
  # OUTPUT:
  # lambda = algebraic value of lambda
  # G = growth rate (sub-adult to adult transition prob)  
  # sad = associated stable age distribution  
  # M = 4 by 4 projection matrix
  br = vrates[1:2]
  wr = vrates[3:4]
  fs = vrates[5:6]
  ms = vrates[7:8]
  SD = 3; # stage duration for sub-adult
  lm = 1; # Initial lambda estimate (use iterations to stabalize)
  for (i in 1:4){
    gf = ((fs[1]/lm)^SD - (fs[1]/lm)^(SD-1))/((fs[1]/lm)^SD -1)
    gm = ((ms[1]/lm)^SD - (ms[1]/lm)^(SD-1))/((ms[1]/lm)^SD -1)
    GF = fs[1]*gf
    PF = fs[1]*(1-gf)
    GM = ms[1]*gm
    PM = ms[1]*(1-gm)
    R = (br/2)*wr*fs
    M = matrix(c(PF,    R[2],   0,      0,    
                 GF,    fs[2],  0,      0,
                 0,     R[2],   PF,     0,
                 0,     0,      GF,     ms[2]),byrow=T,ncol=4) 
    lm=eigen(M)$values[1]    # lambdas=vector of eigenvalues
  }
  lambda = lm 
  W=eigen(M)$vectors          # W=matrix of right eigenvectors 
  w=abs(W[,1])					      # w=stable distribution, unscaled
  ssd = w/sum(w)                # w=stable distribution, scaled
  Gr = c(gf, gm)
  result <- list(lam=lambda,G = Gr,SSD = ssd,M=M)
  return(result)
}
# Generate random matrices ------------------------------
reps = 1000
vrates = matrix(0,nrow=reps,ncol = 8)  # br, wr, fs, ms]
Gvals = matrix(0,nrow=reps,ncol = 2)
lams = numeric()
btar = numeric(); 
L = c(0, .98, 0, .50, .74, .85, .73, .83) 
H = c(0, .98, 0, .90, .96, .98, .93, .95) 
# R-max [Kodiak]
vrates[1,] = H 
rslt = somatrix2st(vrates[1,])   
lams[1] = rslt$lam; Gvals[1,] = rslt$G;  btar[1]= 0
# Stable/Declining pop [Amchitka]
vrates[2,] = L 
rslt = somatrix2st(vrates[2,]) 
lams[2] = rslt$lam; Gvals[2,] = rslt$G;  btar[2]= 0
##
rpar = cbind(runif(reps,0,18),runif(reps,0,18))
btar[2] = 1
for (i in 3:reps){
  btar[i] = rbeta(1,rpar[i,1]+2,rpar[i,2]+2)
  vrates[i,] = L*btar[i] + H*(1-btar[i])
  rslt = somatrix2st(vrates[i,]) 
  lams[i] = rslt$lam; Gvals[i,] = rslt$G;  
}
lams2 = lams^2
Demdat = data.frame(Lam=lams,br1 = vrates[,1],br2 = vrates[,2],
                    wr1 = vrates[,3],wr2 = vrates[,4],
                    fs1 = vrates[,5],fs2 = vrates[,6],
                    ms1 = vrates[,7],ms2 = vrates[,8],
                    Gf = Gvals[,1], Gm = Gvals[,2])
write.csv(Demdat,'RandDem.csv',row.names = F)
Grf = Gvals[,1]; Grm = Gvals[,2]
plot(lams,btar)
plot(lams,Gvals[,1])
plot(lams,Gvals[,2])
fit1 = lm(btar~lams)
summary(fit1)
fit2 = lm(Grf~lams+lams2)
summary(fit2)
fit3 = lm(Grm~lams+lams2)
summary(fit3)



