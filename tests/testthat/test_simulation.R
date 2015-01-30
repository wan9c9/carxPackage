require(carxPackage)

#! test double censoring
message("test double censoring")
rslt <- simulation.carx(nRep=1,fullEstimation=F)
summary(rslt)
plot(rslt$object,main="CARX simulation test plot")
plot(rslt$object,main="CARX simulation test plot",saveFig="testPlot.eps")
#plot(rslt$object,transformFun=exp, main="simulation CARX plot")

#! test lower censoring
message("test lower censoring")
rslt <- simulation.carx(nRep=1,ucl=Inf,fullEstimation=F)
#summary(rslt)
#print(rslt)
plot(rslt$object,main="CARX simulation test plot")
plot(rslt$object,main="CARX simulation test plot",saveFig="testPlot.svg")
#plot(rslt$object,transformFun=exp, main="simulation CARX plot")

