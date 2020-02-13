runsims <- function(){
  # Load files ---------------------------------------------------------
  data = read.csv("./data/GHBlockdata.csv", header = TRUE)  # Data for coastal blocks
  Cdata = read.csv("./data/GHCelldata.csv", header = TRUE)  # Data for habitat cells
  Demdat = read.csv("./data/RandDem.csv", header = TRUE)    # Stochastic vital rates
  Distmat =  read.csv("./data/Distmat.csv", header = FALSE) # Inter-blk LCP distances
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
  source("Kcalc.r")
  # User Parameters  ---------------------------------------------------------
  # (set via UI inputs)
  reps = input$Reps        # Number replications for population sims (should use at least 100)
  Nyrs = input$Nyrs        # Number of years to project population dynamics
  Emax = input$Ephase      # Maximum years before pop "established" (before range expansion begins)
  ImmRt = input$ImmRt      # Immigration rate (avg. net new immigtants per year): 0 = none
  V_sp = input$V_sp        # Population front asymptotic wavespeed, km/yr, minimum  
  K_mean = input$Kmean     # Overall mean K density (modified as fxn of habitat variables)
  K_sig = input$Ksig       # Stochasticity in K (std. deviation in K density)
  sig = input$Estoch       # Environmental stochasticity (std. deviation in log-lambda)
  rmax = input$Rmax        # Maximum rate of growth: default = log(1.22), or 22% per year
  theta = input$Theta      # theta parameter for theta-logistic (1 = Ricker model, >1 = delayed DD) 
  Initpop = input$InitN    # Number of animals in initial population (at least 2 adult females)
  Initblk = input$InitSect # List of initially occupied bloaks: e.g. c(1) = Block 1 only
  #
  # Process data and set up sims --------------------------------------
  Initblk = as.numeric(Initblk)
  if(length(Initblk)==0){
    Initblk = 1
  }
  Yr1 = as.numeric(format(Sys.Date(), "%Y"))
  V_mn = max(.5,V_sp-1)
  V_mx = min(8,V_sp+1)
  Dispers = 2.5  # Over-Dispersion param for Neg Binomial # immigrants per year
  KV = K_sig^2  # Variance in K density
  # rmax = min(log(1.22),rmax)
  ppar = Dispers/(Dispers+ImmRt/length(Initblk))
  Years = c(Yr1:(Yr1+Nyrs-1))  
  Yrs = seq(1:Nyrs)
  Years = Yrs-1+Yr1
  P = nrow(Distmat)  # number blocks (or sub-populations)
  # Initialize population vector
  N0 = numeric(length = P)
  for (i in 1:length(Initblk)){
    N0[Initblk[i]] = round(Initpop/length(Initblk))
  }
  #  Create Data variables for calculating relative density of Hab cells:
  PUID = Cdata$Cell_ID
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
  # Estimate mean K for each Block, and mu/sig for log-normal sampling of K: 
  tmp = Kcalc(PUID,Blk,area,dep,botm,fetch,egrass,parms,K_mean); 
  Ktab=tmp$Ktab; Habdns=tmp$Habdns; Areahab = Ktab$Area
  # NOTE: sum((Habdns$Reldens*area))/sum(area) should =~ 1
  #  sum((Habdns$Reldens*area))/sum(area)
  # *** NOTE: if including uncertainty in hab params, comment out next 2 lines:
  # Kmn = Ktab$Kdns; 
  # muK = log(Kmn/sqrt(1+KV/Kmn^2)); sigK = sqrt(log(1+KV/Kmn^2))
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
  N = array(data = 0, c(P,Nyrs,reps))
  Nmn = matrix(0,nrow=P,ncol=Nyrs)
  Nmn[,1] = N0
  zvec = matrix(data = 0,nrow = 4, ncol = 1)
  #
  withProgress(message = 'Running simulations, please be patient...', 
               session=session, 
  { # for loop for simulations
    for (r in 1:reps){
      incProgress(0.8*(1/reps),session=session)
      # Number of years of population establishment (before range expansion begins)
      E = round(runif(1,round(Emax/4),Emax))
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
                  # NOTE: assymtotic wave speed only obtained once exp growth is occuring
                  probexp = inv.logit(2*(sum(BlokOcc[Adblk,1:(y-1)])-alpha*dst))*max(0,sign(y-E))
                  BlokOcc[j,y] = max(BlokOcc[j,y],rbinom(1,1,probexp))
                }
              }
            }
          }
        }
        # Next, step through blocks and compute dynamics for all occupied blocks
        for (i in 1:P){
          if (BlokOcc[i, y] == 1){
            # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
            lamstoch = max(.95,min(1.22, round(exp(rmax*(1-(N[i,y-1,r]/K[i])^theta)+rnorm(1,0,sig)),2)))
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
            nt[1:4,] = n[1:4,i,y-1]  
            # Account for demographic stochasticity in R
            if(sum(nt)<50 & nt[2]>0){ 
              juvs = rbinom(1,nt[2],br*wr*fsa)
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
            # NEXT LINES ACCOUNT FOR IMMIGRATION 
            #  (Assumes outside immigrants arrive at initially occupied block(s) only)
            if (ImmRt > 0) {
              if (BlokOcc[i, 1] == 1) {
                NImm = round(rnbinom(1,Dispers,ppar))
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
            if (y > E & sum(BlokOcc[,y])>1){
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
    # Process results --------------------------------------------------
    #
    # Density over time
    Dmn = Nmn
    for (y in 1:Nyrs){
      Dmn[,y] = Nmn[,y]/Areahab
    }
    dfDens = data.frame(Blocks = (data$BlockID),Dmn)
    colnames(dfDens) = c('Block',as.character(Years))
    df_Dens <- melt(dfDens, id.vars = "Block")
    names(df_Dens)[2:3] <- c("Year", "Density")
    dfDens = df_Dens[with(df_Dens,order(Block,Year)),]; rm(df_Dens)
    dfDens$Block = as.factor(dfDens$Block)
    
    # Trend plot of abundance over time
    # Calculate mean trend and CI (use bootstrap CI, 1000 reps)
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
      incProgress(0.15*(1/Nyrs),session=session)
      Nsum = colSums(N[,y,])
      bootobj = boot(Nsum, mean.fun, R=500, sim="ordinary")
      means[y] = median(bootobj$t)
      SEmean[y] = sd(bootobj$t)
      tmp = boot.ci(bootobj, type="bca", conf = 0.90); CImn[y,] = tmp$bca[4:5]
      bootobj = boot(Nsum, CIL.fun, R=500, sim="ordinary")
      Lo[y] = median(bootobj$t)
      bootobj = boot(Nsum, CIH.fun, R=500, sim="ordinary")
      Hi[y] = median(bootobj$t)  
    }
    Pop_Overall <- data.frame(Year=as.integer(Years),Mean=means,lower=Lo,upper=Hi,
                              SEmean=SEmean,CImeanLo=CImn[,1],CImeanHi=CImn[,2])
    #
    # Calculate summary of abundance and density by block/year, with CI
    randsmp = sample(seq(1,reps),1000,replace = TRUE)
    meansALL = numeric(length=P*Nyrs)
    LoALL = numeric(length=P*Nyrs)
    HiALL = numeric(length=P*Nyrs)
    YearsALL = numeric(length=P*Nyrs)
    BlockIDs = character(length=P*Nyrs)
    BlockArea = numeric(length=P*Nyrs)
    for (p in 1:P){
      incProgress(.05*(1/P),session=session)
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
        # dens = as.numeric(density(tmp[,y],adjust = 3))
        # CI =  quantile.density(dens, probs=c(.05, .95))
        # Lo[y] <- max(0,as.numeric(CI[1]))
        # Lo[y] <- max(0,as.numeric(CI[1]))
        Lo[y] <- quantile(tmp[,y], probs=c(.05))
        Hi[y] <- quantile(tmp[,y], probs=c(.95))
        means[y] <- mean(tmp[,y])
      }
      meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
      LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
      HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
      YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
      BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Block[p],Nyrs))
      BlockArea[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = rep(Areahab[p],Nyrs)
    }
    Hab_Blocks <- data.frame(Block = BlockIDs, Area = BlockArea, 
                             Year=as.integer(YearsALL), Mean=meansALL, 
                             lower=LoALL,upper=HiALL,Density=dfDens$Density,
                             DensityLO = LoALL/BlockArea, DensityHI = HiALL/BlockArea)
    Hab_Blocks_Fin = Hab_Blocks[which(Hab_Blocks$Year==max(Hab_Blocks$Year)),]     
    datasim = Hab_Blocks_Fin[,c(1,7)]
    HGblk1 = merge(HGblk, datasim, by.x='BlockID', by.y = 'Block')
    titletxt <- paste0("Sea Otter Population Projection, ", Nyrs," Years")
    mapplot2 = ggplot() + 
      geom_polygon(data=HGblk1, aes(x = long, y = lat, fill = Density, color=Density, group = group),
                   alpha = 1,size = 1) +
      scale_fill_continuous(low = "#fff7ec", high = "#7F0000") + 
      scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
      scale_x_continuous(name = "East-west (m)") + # , breaks = NULL, labels = NULL
      scale_y_continuous(name = "North-south (m)") + # , breaks = NULL, labels = NULL
      geom_polygon(data = HGlnd, aes(x=long,y=lat,fill=piece,group=group),
                   color="wheat4", fill="cornsilk1",size = 0.1) +   
      north(HGlnd,location = "topright") +
      scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
               transform = FALSE, location = "bottomleft") +
      ggtitle(titletxt) +
      coord_equal(ratio=1) + theme_minimal()  
  })
  # Update "values" structure
  values$Pop_Overall <- Pop_Overall
  values$dfDens <- dfDens
  values$Hab_Blocks <- Hab_Blocks
  values$Hab_Blocks_Fin <- Hab_Blocks_Fin
  # 
  ggsave("./www/mapdens.png",plot=mapplot2,width = 8,height = 8,dpi=600)
  return(datasim)
}
