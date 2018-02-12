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
Bdata = read.csv("GHBlockdata.csv", header = TRUE)
Cdata = read.csv("GHCelldata.csv", header = TRUE)
# Import movement data from radio tagged sea otters in southern SE Alaska
# to use for fitting sea otter dispersal kernels:
SEmoves = read.csv("SE_Ottermoves.csv", header = TRUE)
# matrix of x-y coordinates of coastl block centroids
BLK = as.matrix(cbind(Bdata$Xcoord,Bdata$Ycoord))
# matrix of x-y coordinates of coastl block centroids
CELL = as.matrix(cbind(Cdata$Xcoord,Cdata$Ycoord))

# Load Raster Map of Haida Gwaii and process for LCP distances ------------------
HAIDA <- raster("HGland.grd") # Raster of Gwaii Haanas, NAD_1983_Albers (BC Environment)
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
write.table(Distmat, file = "Distmat.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")

# Compute Dispersal parameters--------------------------------------------------
# Fit sea otter dispersal kernels for 4 age/sex classes (M/F juv and adult)
#  and use these to compute a) probability of emmigration from each block,
#  and b) dispersal probability matrix, cell i,j = prob that otter emmigrating from 
#  block i ends up in block j (0's along diagonal)
#  ****CODE THIS AND SAVE RESULTS
Wpar = matrix(0,nrow = 4,ncol = 2)
cntr = 0
for (i in 1:2){
  for (j in 1:2){
    cntr = cntr+1
    ii = which(SEmoves$SexN==i & SEmoves$Ageglass==j)
    xi = pmax(0.1,SEmoves$LCD_km[ii])
    fitd = fitdist(xi,"weibull")
    plot(fitd)
    cdfcomp(fitd, addlegend=FALSE)
    Wpar[cntr,1] = fitd$estimate[1]
    Wpar[cntr,2] = fitd$estimate[2]  
  }
}

# Save DispP to GHDispProb.csv
NBlk = length(Bdata$BlockID)
DispP = data.frame(Block = Bdata$BlockID)
DispP$Jf = numeric(length =NBlk)
DispP$Af = numeric(length =NBlk) 
DispP$Jm = numeric(length =NBlk)
DispP$Am = numeric(length =NBlk)
for (i in 1:NBlk){
  dst = numeric(length = Bdata$NAdj[i])
  for (j in 1:Bdata$NAdj[i]){
    dst[j] = Distmat[Bdata[i,5+j],i]
  }
  mndst = median(dst)
  DispP$Jf[i] = 1-pweibull(mndst,Wpar[1,1],Wpar[1,2])
  DispP$Af[i] = 1-pweibull(mndst,Wpar[2,1],Wpar[2,2])
  DispP$Jm[i] = 1-pweibull(mndst,Wpar[3,1],Wpar[3,2])
  DispP$Am[i] = 1-pweibull(mndst,Wpar[4,1],Wpar[4,2])
}
Disp = DispP
Disp[,2:5] = Disp[,2:5]^1.5
write.csv(Disp,"GHDispProb.csv",row.names = FALSE)
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
write.table(GHDispMatJF, file = "GHDispMatJF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatAF, file = "GHDispMatAF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatJM, file = "GHDispMatJM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(GHDispMatAM, file = "GHDispMatAM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")

# Calculate pairwise LCP distances from each hab cell to each block centroid -----
#  (for use in interpolating densities by weighted averaging) 
#
NCell = length(Cdata$PU_ID)
LCD=costDistance(trHAIDA,CELL,BLK) #calculates LCD matrix, distances cell to block centroid 
DistmatHab = as.data.frame(LCD/1000)
#
# Create "HabAvg" matrix: uses inverse distance^2 weighting,
# for row i, col j, value represents proportional contribution for block j 
# on habitat call i
tmp = DistmatHab^2
HabAvg = 1/(tmp+0.000001)
Divisor = rowSums(HabAvg)
for (i in 1:NBlk){
  HabAvg[,i] = HabAvg[,i]/Divisor
}
# Back-up code - if LCD calculaiton explodes with large number of cell values
# N100 = floor(Ncell/100)
# Nremain = Ncell-N100*100
# LCD=costDistance(trHAIDA,CELL,BLK) #calculates LCD matrix, distances cell to block centroid 
# 
# for (i in 2:N100){
#   LCD=costDistance(trHAIDA,CELL[1:100],BLK)  
# 
#   
# }
  
  
