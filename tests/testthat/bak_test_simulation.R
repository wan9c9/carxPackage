#require(carx)

#! test double censoring
message("test double censoring")
rslt <- simulation.carx(nRep=1,nObs=1000,fullEstimation=F)
#summary(rslt)
#print(rslt)
obj <- rslt$object
print(obj)
f <- fitted(obj)
r <- residuals(obj)
fr <- f + r*obj$sigma

plot(obj)
lines(fr,col="red") #impose fitted + residuals
points(fr,col="red",pch='*') #impose fitted + residuals
lines(obj$y,col="green",pch=20) #impose original series
points(obj$y,col="green",pch=20) #impose original series

ts.plot(obj$y - fr)
acf(r[-(1:obj$nAR)])

#find uncensored index
idx <- NULL
for(i in (obj$nAR+1):obj$nObs)
{
	if( all(obj$censorIndicator[i:(i-obj$nAR)] == 0) )
		idx <- c(idx,i)
}
plot(obj$y[idx], fr[idx], main="The plot should be a straight line")

plot(obj)
#test predict.carx
n.ahead <- 10
newdata <- matrix(rnorm(n.ahead*rslt$object$nX),nrow=n.ahead,ncol=rslt$object$nX)
print(newdata)
pred <- predict(rslt$object,newdata,n.ahead)
print(pred)
t <- 1:n.ahead
plot(t,pred$fit)
lines(t,pred$ci[,1])
lines(t,pred$ci[,2])
dev.copy2eps(file="compPlot.eps")
dev.off()

#test plot
#plot(rslt$object,main="CARX simulation test plot")
#plot(rslt$object,main="CARX simulation test plot",saveFig="testPlot.svg")

#plot(rslt$object,transformFun=exp, main="simulation CARX plot")

#! test lower censoring
#message("test lower censoring")
#rslt <- simulation.carx(nRep=1,ucl=Inf,fullEstimation=F)
#summary(rslt)
#print(rslt)
#plot(rslt$object,main="CARX simulation test plot")
#plot(rslt$object,main="CARX simulation test plot",saveFig="testPlot.svg")
#plot(rslt$object,transformFun=exp, main="simulation CARX plot")

