load_all()
options(verbose=FALSE)
message("testing diagnostics")
nObs <- 100
rslt <- simulation.carx(nRep=1,nObs=nObs,trueSigma=0.2, lcl=-0.7,ucl=0.7, fullEstimation=F)
obj <- rslt$object
#plot(obj)

nRep=500
nLag=20
dfChisq = nLag - 6
gof <- goodnessOfFit.carx(obj,nLag=nLag,nRep=nRep,nCores=6)
qqplot(qchisq(ppoints(nRep),df=dfChisq ),gof$bVal)
qqline(gof$bVal,distribution= function(p) qchisq(p,df=dfChisq),prob=c(0.1,0.6),col=2)


