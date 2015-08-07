
nAR = 2
nObs = 200
trueX = c(0.2,0.4)
trueAR = c(-0.28,0.25)
trueSigma = 0.60
lcl = -1
ucl = 1
dat = simulateCarx(nObs, trueAR, trueX, trueSigma, lcl, ucl, seed=0)

dframe <- as.data.frame(list(y = dat$y,  = dat$x))

carx(y ~ x1 + x2, data=as.data.frame(dat$x),
     censorIndicator = dat$censorIndicator,
     lowerCensorLimit = dat$lowerCensorLimit,
     upperCensorLimit = dat$upperCensorLimit,
     nAR = nAR)

