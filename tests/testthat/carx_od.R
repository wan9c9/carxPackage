trueX = c(0.2,0.4)
nAR = 2
trueAR = c(-0.28,0.25)
trueSigma = 0.60
lcl = -1
ucl = 1
dat = simulateCarx(100, trueAR, trueX, trueSigma, lcl, ucl, seed=0)
std= trueSigma
dat$y[100]=dat$y[100]+5*std
dat$y[40]=dat$y[40]+3*std
dat$y[80]=dat$y[80]+4*std
dat$y[60]=dat$y[60]+2.9*std
plot(dat$y,type="l")
model0 <- carx(dat$y, dat$x,censorIndicator=dat$censorIndicator,
               lowerCensorLimit = dat$lowerCensorLimit,
               upperCensorLimit = dat$upperCensorLimit,
               nAR = nAR,
               getCI = FALSE)
rslt <- outlier(model0)

data(phosphorusWestForkCedarRiverAtFinchford)
data1 <- phosphorusWestForkCedarRiverAtFinchford
std=sd(data1$logP)
data1$logP[100]=data1$logP[100]+3*std
data1$logP[40]=data1$logP[40]+3*std
data1$logP[80]=data1$logP[80]+1.5*std
tmp <- carx(logP~logQ, data=data1, censorIndicator=data1$censorIndicator,
            lowerCensorLimit=data1$lowerCensorLimit,upperCensorLimit=NULL,
            nAR=2, getCI=FALSE)

rslt <- outlier(tmp,data1)
plot(rslt$updatedModel)
