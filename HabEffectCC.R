# Analysis of sea otter habitat associations, central coast BC
# Depth dxn fit -----------------------------------------------
library(ggplot2)
library(gtools)
library(stats)
data = read.csv("CC_depth_obs.csv", header = TRUE)  # Data on coastal otters
#logistic Regression Model
dat = data[data$Depth<=40,]
dat$Depth2 = dat$Depth^2
# data$Depth3 = data$Depth^3
fit <- glm(Sighting~Depth+Depth2,data=dat,family=binomial(link="logit"),weights = Count)
summary(fit)
xi = seq(0,40)
newdat = data.frame(Depth=xi,Depth2 = xi^2)
preddat = predict(fit, newdata=newdat, se.fit=TRUE)
yi = predict(fit,newdata=newdat,type="response")
y_lo = exp(preddat$fit-1.96*preddat$se.fit)/(1+exp(preddat$fit-1.96*preddat$se.fit))
y_hi = exp(preddat$fit+1.96*preddat$se.fit)/(1+exp(preddat$fit+1.96*preddat$se.fit))
plot(dat$Depth,dat$Sighting,main="Otter use vs Depth, BC Central Coast",
     xlab="Depth",ylab="Probability of Sighting")
lines(xi,yi,col="blue")
lines(xi,y_lo,lty=2)
lines(xi,y_hi,lty=2)

# BoP dxn fit -----------------------------------------------
dat2 = read.csv("CC_BoP.csv", header = TRUE)  # Data on coastal otters
dat2$DepthCode = factor(dat2$DepthCode, levels = c("0to5", "5to10", "10to20", "20to50"))
dat2$count = dat2$count_tota
# dat2$count[dat2$Obs==0] = 2
fit2 <- glm(Obs~HabCode+DepthCode,data=dat2,family=binomial(link="logit"),weights = sqrt(count)) #,weights = count
# fit2 <- glm(Obs~HabCode+DepthCode+HabCode*DepthCode,data=dat2,family=binomial(link="logit"),weights = count)
summary(fit2)
newdat2 = unique(dat2[,3:4])
newdat2a = cbind(HabCode=unique(dat2[,3]),data.frame(DepthCode=rep("0to5",length(unique(dat2[,3])))))
newdat2b = cbind(data.frame(HabCode=rep("1a",length(unique(dat2[,4])))),DepthCode=unique(dat2[,4]))
preddat2a = predict(fit2, newdata=newdat2a, se.fit=TRUE)
preddat2b = predict(fit2, newdata=newdat2b, se.fit=TRUE)
# dfHb = data.frame(BoP = newdat2a$HabCode, Mean = 1.3*exp(preddat2a$fit)/(1+exp(preddat2a$fit)),
#                   Lo = 1.3*exp(preddat2a$fit-1.96*preddat2a$se.fit)/(1+exp(preddat2a$fit-1.96*preddat2a$se.fit)),
#                   Hi = 1.3*exp(preddat2a$fit+1.96*preddat2a$se.fit)/(1+exp(preddat2a$fit+1.96*preddat2a$se.fit)))

dfHb = data.frame(BoP = newdat2a$HabCode, Mean = preddat2a$fit,
                  SE = preddat2a$se.fit ) 
dfHb$Parest = numeric(length=9)
dfHb$Par_Lo = numeric(length=9)
dfHb$Par_Hi = numeric(length=9)
tmp = inv.logit(mean(dfHb$Mean))
for (j in 1:9){
  tmp2 = inv.logit(rnorm(10000,dfHb$Mean[j],dfHb$SE[j]))/tmp
  dfHb$Parest[j] = log(mean(tmp2))
  dfHb$Par_Lo[j] = log(quantile(tmp2,0.025))
  dfHb$Par_Hi[j] = log(quantile(tmp2,0.975))
}

plt1 = ggplot(dfHb, aes(BoP, Parest)) + 
  geom_col() +  
  geom_errorbar(aes(ymin = Par_Lo, ymax = Par_Hi), width=0.2) +
  xlab("Bottom Type") +
  ylab("Log(Relative Density)") +
  ggtitle("Habitat Effects on Density, Central Coast")
print(plt1)

dfDp = data.frame(Depth = newdat2b$DepthCode, Mean = preddat2b$fit,
                  SE = preddat2b$se.fit)
dfDp$Parest = numeric(length=4)
dfDp$Par_Lo = numeric(length=4)
dfDp$Par_Hi = numeric(length=4)
tmp = inv.logit(mean(dfDp$Mean))
for (j in 1:4){
  tmp2 = inv.logit(rnorm(10000,dfDp$Mean[j],dfDp$SE[j]))/tmp
  dfDp$Parest[j] = log(mean(tmp2))
  dfDp$Par_Lo[j] = log(quantile(tmp2,0.025))
  dfDp$Par_Hi[j] = log(quantile(tmp2,0.975))
}

plt2 = ggplot(dfDp, aes(Depth, Parest)) + 
  geom_col() +  
  geom_errorbar(aes(ymin = Par_Lo, ymax = Par_Hi), width=0.2) +
  xlab("Depth Class") +
  ylab("Log(Relative Density)") +
  ggtitle("Depth Effects on Density, Central Coast")  
print(plt2)




