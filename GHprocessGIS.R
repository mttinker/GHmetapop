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
# Data for coastl blocks around Haida Gwaii
data = read.csv("GHBlockdata.csv", header = TRUE)
# Import movement data from radio tagged sea otters in southern SE Alaska
# to use for fitting sea otter dispersal kernels:
SEmoves = read.csv("SE_Ottermoves.csv", header = TRUE)
# matrix of x-y coordinates of coastl block centroids
BLK = as.matrix(cbind(data$Xcoord,data$Ycoord))

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
LCD=costDistance(trHAIDA,BLK,BLK) #calculates LCD matrix
Distmat = LCD/1000   # COnvert ditances to km
# Save matrix as csv file to load in metapopulaiton model
write.table(Distmat, file = "Distmat.csv.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
# Compute Dispersal parameters--------------------------------------------------
# Fit sea otter dispersal kernels for 4 age/sex classes (M/F juv and adult)
#  and use these to compute a) probability of emmigration from each block,
#  and b) dispersal probability matrix, cell i,j = prob that otter emmigrating from 
#  block i ends up in block j (0's along diagonal)
#  ****CODE THIS AND SAVE RESULTS

# Save DispP to GHDispProb.csv
NBlk = length(data$BlockID)
DispP = data.frame(Block = data$BlockID)
DispP$Jf = numeric(length =NBlk)
DispP$Af = numeric(length =NBlk) 
DispP$Jm = numeric(length =NBlk)
DispP$Am = numeric(length =NBlk)


# Inter-pop dispersal matrices for each age/sex class
#  (save each as matrix using wrtie.table)
GHDispMatJF = matrix(data = 0,nrow = NBlk,ncol = NBlk)
GHDispMatAF = matrix(data = 0,nrow = NBlk,ncol = NBlk)
GHDispMatJM = matrix(data = 0,nrow = NBlk,ncol = NBlk)
GHDispMatAM = matrix(data = 0,nrow = NBlk,ncol = NBlk)

