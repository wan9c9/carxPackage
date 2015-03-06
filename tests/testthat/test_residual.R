#require(carx)
#in this section, we test the independence of the residuals calculated from our method
load_all()
options(verbose=TRUE)
message("testing calculation of residuals")
nObs <- 500
rslt <- simulation.carx(nRep=1,nObs=nObs,trueSigma=0.2, lcl=-0.7,ucl=0.7, fullEstimation=F)
obj <- rslt$object

r <- residuals(obj)
svg("testRes_resAcf.svg")
acf(r[-(1:obj$nAR)])
dev.off()

#f <- fitted(obj)
#fr <- f + r*obj$sigma
