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