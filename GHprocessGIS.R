# Script to process spatial data for Gwaii Haanas sea otter populaiton model
#
# Load necessary libraries
library(gdistance)
library(raster)
library(rgdal) 
library(sp) 
library(fields)
library(fitdistrplus)
#
# Load Data ----------------------------------------------------------------------
# Data for coastl blocks and habitat cells around Haida Gwaii
Bdata = read.csv("./data/GHBlockdata.csv", header = TRUE)
Cdata = read.csv("./data/GHCelldata.csv", header = TRUE)
Cdata$PU_ID = Cdata$Cell_ID
# Import movement data from radio tagged sea otters in southern SE Alaska
# to use for fitting sea otter dispersal kernels:
SEmoves = read.csv("./data/SE_Ottermoves.csv", header = TRUE)
CAmoves = read.csv("./data/ottmovesBSMB.csv", header = TRUE)
ir1 = sample.int(nrow(SEmoves),2000,replace=T)
ir2 = sample.int(nrow(CAmoves),10000,replace=F)
ottmoves = data.frame(dist = SEmoves$LCD_km[ir1],Sex=SEmoves$SexN[ir1],
                Age=SEmoves$Ageglass[ir1])
ottmoves =rbind(ottmoves, 
                data.frame(dist = CAmoves$move[ir2],Sex=CAmoves$Sex[ir2]+1, 
                Age=CAmoves$Ageclass[ir2]+1))

# matrix of x-y coordinates of coastl block centroids
BLK = as.matrix(cbind(Bdata$Xcoord,Bdata$Ycoord))
# matrix of x-y coordinates of grid cell centroids
CELL = as.matrix(cbind(Cdata$Xcoord,Cdata$Ycoord))

# Load Raster Map of Haida Gwaii and process for LCP distances ------------------
HAIDA <- raster("./data/HGland.grd") # Raster of Gwaii Haanas, NAD_1983_Albers (BC Environment)
plot(HAIDA) #SEAK Raster
points(BLK) #plots Block Centroid points
# Process data for estimating least cost path (LCP) distances between points
trHAIDA<- transition(1/HAIDA, mean, directions=8) ##8 directions 
#correct distances (map distortion) for diagonal movements and large area of study site
trHAIDA=geoCorrection(trHAIDA) 
# Compute euclidean distance and LCP distance matrices 
EUC=rdist(BLK) #euclidian distances
LCD=costDistance(trHAIDA,BLK,BLK) #calculates LCD matrix, pairwise block-block distances 
Distmat = LCD/1000   # Convert pairwise LCP distances to km
# Save matrix as csv file to load in metapopulaiton model
write.table(Distmat, file = "./data/Distmat.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")

# Compute Dispersal parameters--------------------------------------------------
Distmat =  as.matrix(read.csv("./data/Distmat.csv", header = FALSE)) # Inter-blk LCP distances
# Fit sea otter dispersal kernels for 4 age/sex classes (M/F juv and adult)
#  and use these to compute a) probability of emmigration from each block,
#  and b) dispersal probability matrix, cell i,j = prob that otter emmigrating from 
#  block i ends up in block j (0's along diagonal)
#  ****CODE THIS AND SAVE RESULTS
Wpar = matrix(0,nrow = 4,ncol = 2)
Xpar = matrix(0,nrow = 4,ncol = 1)
cntr = 0
for (i in 1:2){
  for (j in 1:2){
    cntr = cntr+1
    # ii = which(SEmoves$SexN==i & SEmoves$Ageglass==j)
    # xi = pmax(0.1,SEmoves$LCD_km[ii])
    ii = which(ottmoves$Sex==i &ottmoves$Age==j)
    xi = pmax(0.1,ottmoves$dist[ii])
    fitd = fitdist(xi,"weibull")
    fitx = fitdist(xi,"exp")
    plot(fitd)
    # cdfcomp(list(fitd,fitx), addlegend=TRUE)
    Wpar[cntr,1] = fitd$estimate[1]
    Wpar[cntr,2] = fitd$estimate[2]  
    Xpar[cntr,1] = fitx$estimate  
  }
}
# Save DispP to GHDispProb.csv
NBlk = length(Bdata$BlockID)
DispP = data.frame(Block = Bdata$BlockID)
DispP$Jf = numeric(length =NBlk)
DispP$Af = numeric(length =NBlk) 
DispP$Jm = numeric(length =NBlk)
DispP$Am = numeric(length =NBlk)
DispPx = DispP
for (i in 1:NBlk){
  mndst = 0.75*mean(sort(Distmat[,i],decreasing=F)[2:3])
  DispP$Jf[i] = 1-pweibull(mndst,Wpar[1,1],Wpar[1,2])
  DispP$Af[i] = 1-pweibull(mndst,Wpar[2,1],Wpar[2,2])
  DispP$Jm[i] = 1-pweibull(mndst,Wpar[3,1],Wpar[3,2])
  DispP$Am[i] = 1-pweibull(mndst,Wpar[4,1],Wpar[4,2])
  DispPx$Jf[i] = 1-pexp(mndst,Xpar[1,1])
  DispPx$Af[i] = 1-pexp(mndst,Xpar[2,1])
  DispPx$Jm[i] = 1-pexp(mndst,Xpar[3,1])
  DispPx$Am[i] = 1-pexp(mndst,Xpar[4,1])
}
Disp = DispP
write.csv(Disp,"./data/GHDispProb.csv",row.names = FALSE)
# Inter-pop dispersal matrices for each age/sex class
#  (save each as matrix using write.table)
GHDispMatJF = 1-pweibull(Distmat,Wpar[1,1],Wpar[1,2])
diag(GHDispMatJF) = 0
GHDispMatAF = 1-pweibull(Distmat,Wpar[2,1],Wpar[2,2])
diag(GHDispMatAF) = 0
GHDispMatJM = 1-pweibull(Distmat,Wpar[3,1],Wpar[3,2])
diag(GHDispMatJM) = 0
GHDispMatAM = 1-pweibull(Distmat,Wpar[4,1],Wpar[4,2])
diag(GHDispMatAM) = 0
SumJF = colSums(GHDispMatJF)
SumAF = colSums(GHDispMatAF)
SumJM = colSums(GHDispMatJM)
SumAM = colSums(GHDispMatAM)
for (i in 1:NBlk){
  GHDispMatJF[,i] = GHDispMatJF[,i]/SumJF[i]
  GHDispMatAF[,i] = GHDispMatAF[,i]/SumAF[i]
  GHDispMatJM[,i] = GHDispMatJM[,i]/SumJM[i]
  GHDispMatAM[,i] = GHDispMatAM[,i]/SumAM[i]  
}
write.table(GHDispMatJF, file = "./data/GHDispMatJF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatAF, file = "./data/GHDispMatAF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatJM, file = "./data/GHDispMatJM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatAM, file = "./data/GHDispMatAM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
# Calculate pairwise LCP distances from each hab cell to each block centroid -----
#  (for use in interpolating densities by weighted averaging) 
# NOTE: NO LONGER USING THIS
# LCD=costDistance(trHAIDA,CELL,BLK) #calculates LCD matrix, distances cell to block centroid 
# DistmatHab = as.data.frame(LCD/1000)
# #
# # Create "HabAvg" matrix: uses inverse distance^2 weighting,
# # for row i, col j, value represents proportional contribution for block j 
# # on habitat call i
# tmp = DistmatHab^2
# HabAvg = 1/(tmp+0.000001)
# Divisor = rowSums(HabAvg)
# for (i in 1:NBlk){
#   HabAvg[,i] = HabAvg[,i]/Divisor
# }
# HabAvg = data.frame(cbind(Cdata$PU_ID,HabAvg))
# colnames(HabAvg) = c('PU_ID',as.character(Bdata$BlockID))  
# write.csv(HabAvg,"HabAvg.csv",row.names = FALSE)  
detach(package:gdistance)
detach(package:raster)