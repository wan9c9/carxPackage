#require(carx)

#! test double censoring
message("test double censoring")
rslt <- simulation.carx(nRep=1,fullEstimation=F)
summary(rslt)
print(rslt)
obj <- rslt$object
f <- fitted(obj)
r <- residuals(obj)
print(obj)
resid(obj)

fr <- f + r*obj$sigma
#find uncensored index
idx <- NULL
for(i in (obj$nAR+1):obj$nObs)
{
	if( all(obj$censorIndicator[i:(i-obj$nAR)] == 0) )
		idx <- c(idx,i)
}
plot(obj$y[idx], fr[idx], main="The plot should be a straight line")

plot(obj)
lines(fr,col="red") #impose fitted + residuals
lines(obj$y,col="green") #impose original series

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

