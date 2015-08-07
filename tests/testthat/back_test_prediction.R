
load_all()


a <- arima.sim(n=10,list(ar=c(0.8,-0.5)), sd=0.1)

#test for predictARX
#pure AR
p <- predictARX(y=c(1,-1),x=c(0,0),prmtrX=c(0),prmtrAR=c(0.5),sigma=1,n.ahead=1,newxreg=c(0),CI.level=0.95)
expect_equal(p$fit, -0.5)
expect_equal(p$se.fit,1)

p <- predictARX(y=c(1,-1),x=c(0,0),prmtrX=c(0),prmtrAR=c(0.5,-0.2),sigma=1,n.ahead=1,newxreg=c(0),CI.level=0.95)
expect_equal(p$fit, -0.7)
expect_equal(p$se.fit,1)

p <- predictARX(y=c(1,-1),x=matrix(c(1,-1,1,-1),nrow=2,ncol=2),prmtrX=c(1,2),
                prmtrAR=c(0.5,-0.2),sigma=1,n.ahead=2,newxreg=matrix(c(1,-1,1,-1),nrow=2,ncol=2),
                CI.level=0.95)
# x =    
#       [,1] [,2]
#[1,]    1    1
#[2,]   -1   -1

expect_equal(p$fit, c(4.4,-2.7))
expect_equal(p$se.fit,c(1.0,1.25))


nObs <- 200
n.ahead <- 10

cts <- carx.sim.cenTS(nObs=nObs,sigmaEps=0.5)
plot(cts)

currentTS <- cts[1:(nObs-n.ahead)]
futureTS <- cts[-(1:(nObs-n.ahead))]

m <- carx(value~X1+X2,data=currentTS,p=2,CI.compute=FALSE)
#debug(predict.carx)
p <- predict(m,n.ahead=n.ahead,newxreg=coredata(futureTS))

plot(coredata(futureTS)[,"value"],ylim=range(data.frame(p,coredata(futureTS))))
lines(p$fit)
lines(p$ci[,1],col="red")
lines(p$ci[,2],col="red")

###############################################
cts <- carx.sim.cenTS(nObs=nObs,prmtrX=0.2,x=as.matrix(rep(1,nObs),nrow=nObs,ncol=1),sigmaEps=0.5)
plot(cts)

currentTS <- cts[1:(nObs-n.ahead)]
futureTS <- cts[-(1:(nObs-n.ahead))]

m <- carx(value~1,data=currentTS,p=2,CI.compute=FALSE)
#debug(predict.carx)
p <- predict(m,n.ahead=n.ahead,newxreg=coredata(futureTS))

plot(coredata(futureTS)[,"value"],ylim=range(data.frame(p,coredata(futureTS))))
lines(p$fit)
lines(p$ci[,1],col="red")
lines(p$ci[,2],col="red")


