load_all()
options(verbose=TRUE)
message("testing plot")
nObs <- 100
rslt <- simulation.carx(nRep=1,nObs=nObs,trueSigma=0.2, lcl=-0.7,ucl=0.7, fullEstimation=F)
obj <- rslt$object
plotData <- TRUE
if(plotData)
{
	svg("testRes_comPlot.svg")
	plot(obj)
	lines(obj$y)
	dev.off()
}


