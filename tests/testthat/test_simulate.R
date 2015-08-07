load_all()
rm(list=ls())

nObs = 100
trueX = c(0.2,0.4)
trueAR = c(-0.28,0.25)
trueSigma = 0.60
lcl = -1
ucl = 1

#debug(carx.sim.cenTS)
ts = carx.sim.cenTS(100, trueAR, trueX, trueSigma, lcl, ucl, seed=0)

m <- carx(value~X1+X2,data=ts,CI.compute=FALSE)

n.ahead=1
newdata <- matrix(rnorm(
p <- predict(m,n.ahead=3,)

