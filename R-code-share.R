##read in data
rm(list=ls())
library(prospectr)
library(pls)

rootdir <- "C:/Users/yge2/Desktop/"
setwd(paste0(rootdir, "VisNIR leaf property experiments/Workshop"))
labdata <- read.table("labdata.csv", header=T, sep=",")
specR <- read.table("specR.csv", header=T, sep=",")

wl <- 350:2500 # wavelength 350 to 2500 nm
wl10 <- seq(410, 2490, 10) # wavelength every 10 nm from 410 to 2500 nm, 
# center wavelength 410, 420, 430,...


specR.avg <- binning(specR[, 56:2145], bin.size=10) # start from 405 nm, end at 2494 nm, average every 10 nm
matplot(wl10, t(specR.avg), type="l") # examine all reflectance readings
# average the reflectance spectra to every 10 nm 

N1 <- 300
set.seed(10)
val.set <- seq(3, N1, by=3)
cal.set <- (1:N1)[-val.set]

################### vegetation indices############

# GNDVI (Green NDVI) 550 and 800 nm
# MCARI chlorophyll index 

b550 <- which(wl10 == 550)
b800 <- which(wl10 == 800)

GNDVI <- (specR.avg[, b800]-specR.avg[, b550])/(specR.avg[, b800]+specR.avg[, b550])
plot(GNDVI, labdata$N)

b850 <- which(wl10 == 850)
b730 <- which(wl10 == 730)
b570 <- which(wl10 == 570)

MCARI <- ((specR.avg[, b850]-specR.avg[,b730])-0.2*(specR.avg[, b850]-specR.avg[,b570]))/specR.avg[, b730]
plot(MCARI, labdata$chl)

### results of using MCARI index for chlorophyll estimation
lm.chl <- lm(labdata$chl[cal.set] ~ MCARI[cal.set]+I(MCARI[cal.set]^2))
pred.lm.chl <- lm.chl$coefficients[1]+lm.chl$coefficients[2]*MCARI[val.set]+lm.chl$coefficients[3]*MCARI[val.set]^2
plot(labdata$chl[val.set], pred.lm.chl, pch=19, cex=1.2,
     xlab = "measured chl", ylab="predicted chl")
abline(0,1,lwd=2,col=2)
summary(lm(pred.lm.chl~labdata$chl[val.set])) #R2 = 0.91
sqrt(mean((pred.lm.chl-labdata$chl[val.set])^2)) #RMSE = 31.7 umol/m2

lm.N <- lm(labdata$N[cal.set] ~ GNDVI[cal.set]+I(GNDVI[cal.set]^2))
pred.lm.N <- lm.N$coefficients[1]+lm.N$coefficients[2]*GNDVI[val.set]+lm.N$coefficients[3]*GNDVI[val.set]^2
plot(labdata$N[val.set], pred.lm.N, pch=19, cex=1.2,
     xlab="measured N", ylab="redicted N")
abline(0,1,lwd=2, col=2)
summary(lm(pred.lm.N ~ labdata$N[val.set])) #R2 = 0.67
sqrt(mean((pred.lm.N - labdata$N[val.set])^2)) #RMSE = 0.376%

# develop plsr model for chlorophyll use calibration set
df.chl <- data.frame(Y=labdata$chl, X=I(as.matrix(specR.avg)))
plsr.chl <- plsr(Y~X, data=df.chl[cal.set,], ncomp=20, val="LOO")
validationplot(plsr.chl, val.type="RMSEP", estimate = c("train","CV"), type="b")
# lowest RMSEP is at ncomp=17
pred.plsr.chl <- predict(plsr.chl, newdata=df.chl[val.set,], ncomp=12)
plot(labdata$chl[val.set], pred.plsr.chl, pch=19, cex=1.2, 
     xlab="measured chl", ylab="predicted chl")
abline(0, 1, lwd=2, lty=1, col=2)
RMSEP(plsr.chl, newdata=df.chl[val.set,], ncomp=12) # RMSE = 24.02 umol/m2
R2(plsr.chl, newdata=df.chl[val.set,], ncomp=12) #R2 = 0.95

# develop plsr model for N - using calibration set
df.N <- data.frame(Y=labdata$N, X=I(as.matrix(specR.avg)))
plsr.N <- plsr(Y~X, data=df.N[cal.set, ], ncomp=20, val="LOO")
validationplot(plsr.N, val.type="RMSEP", estimate = c("train","CV"), type="b")
RMSEP(plsr.N, newdata=df.N[val.set, ], ncomp=11) #RMSE=0.289
R2(plsr.N, newdata=df.N[val.set,], ncomp=11)  #R2=0.80

pred.plsr.N <- predict(plsr.N, newdata=df.N[val.set, ], ncomp=11)
plot(labdata$N[val.set], pred.plsr.N, pch=19, cex=1.2,
     xlab="measured N", ylab="predicted N")
abline(0, 1, lwd=2, lty=1, col=2)

# LMA modeling and estimation with PLSR
df.LMA <- data.frame(Y=labdata$LMA, X=I(as.matrix(specR.avg)))
plsr.LMA <- plsr(Y~X, data=df.LMA[cal.set,], ncomp=20, val="LOO")
validationplot(plsr.LMA, val.type="RMSEP", estimate = c("train","CV"), type="b")
RMSEP(plsr.LMA, newdata=df.LMA[val.set,], ncomp=13)
R2(plsr.LMA, newdata=df.LMA[val.set,], ncomp=13)
plot(plsr.LMA, newdata=df.LMA[val.set,], ncomp=13, asp=1, line=T)


# Leaf S concentration modeling and estimation with PLSR
df.S <- data.frame(Y=labdata$S, X=I(as.matrix(specR.avg)))
plsr.S <- plsr(Y~X, data=df.S[cal.set, ], ncomp=20, val="LOO")
validationplot(plsr.S, val.type="RMSEP", estimate = c("train","CV"), type="b")
RMSEP(plsr.S, newdata=df.S[val.set,], ncomp=11)
R2(plsr.S, newdata=df.S[val.set,], ncomp=11)

Pred.plsr.S <- predict(plsr.S, newdata=df.S[val.set, ], ncomp=11)
plot(labdata$S[val.set], Pred.plsr.S, pch=19, cex=1.2,
     xlab="measured S", ylab="predicted S")
abline(0, 1, lwd=2, lty=1, col=2)