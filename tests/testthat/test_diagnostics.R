load_all()
options(verbose=TRUE)
message("testing diagnostics")
#nObs <- 500
nLag <- 10
#rslt <- simulation.carx(nRep=1,nObs=nObs,trueSigma=0.2, lcl=-1.0,ucl=1.0, fullEstimation=F)
rslt <- simulation.carx(nRep=1)
obj <- rslt$object
res <- residuals(obj)
ac <- acf(res)
lb <- lbStat(res,nLag)


#res <- residuals(obj,method="rr")

nRep=1000
nLag=20
dfChisq = nLag - 6
gof <- goodnessOfFit.carx(obj,nLag=nLag,nRep=nRep,nCores=6)
qqplot(qchisq(ppoints(nRep),df=dfChisq ),gof$bVal)
qqline(gof$bVal,distribution= function(p) qchisq(p,df=dfChisq),prob=c(0.1,0.9),col="red")
qqline(gof$bVal,distribution= function(p) qchisq(p,df=dfChisq-1),prob=c(0.1,0.9),col="blue")
qqline(gof$bVal,distribution= function(p) qchisq(p,df=dfChisq+1),prob=c(0.1,0.9),col="blue")


#A simulation study to assess the gof proposed.
nRep= 1000
require(foreach)
require(doMC)
registerDoMC(7)
foreach(idx=1:nRep,.combine='c')
	

simulation2.carx <-function(trueEV=c(0.2,0.4), 
													 trueAR=c(-0.28,0.25), 
													 trueSigma=0.60,
													 lcl=-1,
													 ucl=1,
													 nObs=100,
													 nRep=1000,
													 alpha=0.95,
													 nLag = 10,
													 fullEstimation = FALSE,
													 testPower=0)
{
	message(c("Simulation study begins at ",date()))
	t0 <- proc.time()
	require(matrixStats)
	args <- commandArgs(TRUE)
	if(length(args) > 0)
	{
		for(i in 1:length(args))
		{
			eval(parse(text = args[i]))
		}
	}
	lowercl <- rep(lcl, nObs)
	uppercl <- rep(ucl, nObs)

  
	#test if AR is stationary
	roots <- polyroot(c(1,-trueAR))
	message("roots of polynomials of trueAR, see if it's stationary")
	print(roots)
	if(any(abs(roots)<=1))
	{
		message("supplied trueAR is not statoinary")
	}
	
	nAR <- length(trueAR)
  trueAR1 <- numeric(nAR+1)
	if(testPower!=0)
	{
		trueAR1[1] <- trueAR[1] + testPower
		for(i in 2:nAR)
			trueAR1[i] <- trueAR[i] - testPower*trueAR[i-1]
		trueAR1[nAR+1] <- -testPower*trueAR[nAR]
	}
	else
		trueAR1 <- trueAR
	message("original arPrmtr")
	print(trueAR)
	message("disturbed arPrmtr")
	print(trueAR1)
	roots <- polyroot(c(1,-trueAR1))
	message("roots of polynomials of trueAR, see if it's stationary")
	print(roots)


	truePrmtr <- c(trueEV,trueAR,trueSigma)
  nX <- length(trueEV)
	nPrmtr <- length(truePrmtr)

	prmtrEstd <- matrix(0,nRep,nPrmtr)
	prmtrInit <- matrix(0,nRep,nPrmtr)
	estCI <- matrix(0,nRep,2*nPrmtr)
	coverage <- numeric(nPrmtr)

	#' return numberOfReplication, censorRate, prmrtrInit, prmtrEstd 
	#' ( lower_ci, upper_ci, coverageIndicator)
	simEst <-function(iRep)
	{
		message(sprintf("Rep:%i",iRep))

		dat <- carx.simulate(nObs, trueAR1, trueEV, trueSigma, lowercl, uppercl, seed=37513*iRep)
		rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$lowerCensorLimit,dat$upperCensorLimit, nAR, getCI=fullEstimation,skipIndex=seq(1,nAR))
		ret <- c(iRep,rslt$censorRate,rslt$prmtrInit,rslt$prmtrEstd)
		if(fullEstimation)
		{
			ci <- rslt$CI
			ret <- c(ret, ci[,1])
			ret <- c(ret, ci[,2])
			coverage <- (truePrmtr >= ci[,1])*(truePrmtr <= ci[,2])
			ret <- c(ret,coverage)
		}
		
		rslt0 <- carx(dat$y,dat$x,rep(0,nObs),NULL,NULL, nAR, getCI=fullEstimation,skipIndex=seq(1,nAR))
		eta0 <- dat$y - dat$x%*%rslt0$prmtrX
		rsdl0 <- numeric(nObs)
		for(i in (nAR+1):nObs)
			rsdl0[i] <- eta0[i] - eta0[(i-1):(i-nAR)]%*%rslt0$prmtrAR
		rtest0 <- lbStat(rsdl0,nLag)

		rsdl <- residuals(rslt)
		#print(rsdl)
		rtest <- lbStat(rsdl,nLag)
		#print(rtest)
		ret <- c(ret,rtest0,rtest)
		#print(ret)
		list(object=rslt,summary=ret)
	}

	if(nRep == 1)
	{
		rslt <- simEst(nRep)
		message(sprintf("Time used:"))
		print(proc.time()-t0)
		message(c("Simulation study ends at ",date())) 
		return(rslt)
	}

	iter <- 1:nRep
	rslt <- mclapply(iter,simEst,mc.cores=detectCores())
	#print(rslt)

	if(fullEstimation) 
		summaryList <- matrix(nrow=nRep,ncol=(4+5*nPrmtr))
	else
		summaryList <- matrix(nrow=nRep,ncol=(4+2*nPrmtr))
	i <- 1
	objects <- NULL
	for(r in rslt) 
	{
		objects <-c(objects,r$object)
		#print(r$summary)
		summaryList[i,] <- r$summary
		i <- i+1
	}
	#print(summaryList)

	rslt <- summaryList
	colm <- colMeans(rslt)
	avgCR <- colm[2]
	avgEstd <- colm[(2+nPrmtr+1):(2+2*nPrmtr)]
	std <- colSds(rslt)
	std <- std[(2+nPrmtr+1):(2+2*nPrmtr)]

	if(fullEstimation)
		coverageRate <- colm[(2+4*nPrmtr+1):(2+5*nPrmtr)]
	else
		coverageRate <- rep(NaN, nPrmtr)

	rejectionRate0 <- sum(rslt[,dim(rslt)[2]-1] > qchisq(alpha, nLag-nAR))/nRep
	rejectionRate <- sum(rslt[,dim(rslt)[2]] > qchisq(alpha, nLag-nAR))/nRep

	summary <- cbind(truePrmtr,avgEstd,std,coverageRate)
	colnames(summary) <- c('true', 'avg.estd','std.err','coverage.rate')
	xnames <- paste0('X',1:length(trueEV))
	rnames <- c(xnames,paste0('AR',1:nAR),"sigma")
	rownames(summary) <- rnames

	rslt <- list(nRep=nRep,
							 #objects = objects,
							 averageCensorRate=avgCR,
							 alpha = alpha,
							 allResult = summaryList,
							 rejectionRate0 = rejectionRate0,
							 rejectionRate = rejectionRate,
							 summary = summary)

	message(sprintf("\nReplication: %i\n", nRep))
	message(sprintf("\nAverage censor rate: %f\n", avgCR))
	message(sprintf("\n rejection rate0: %f\n", rejectionRate0))
	message(sprintf("\n rejection rate: %f\n", rejectionRate))
	message("trueAR1")
	print(trueAR1)

	message("                     Simulation Summary        ")
	if(testPower != 0)
	{
		message("true arPrmtr are")
		print(trueAR1)
		message("WARNING: the true value in the following table are wrong")
	}
	message("     true          meanEstd       stdError     coverageRate ")
	for(i in 1:nPrmtr)
	{
		m <- c(sprintf("%10.3f ",truePrmtr[i]),
					 sprintf("\t%10.3f ",avgEstd[i]),
					 sprintf("\t%10.3f ",std[i]), 
					 sprintf("\t%10.3f ",coverageRate[i])
					 )
		message(m)
	}
	message(sprintf("Time used:"))
	print(proc.time()-t0)
	message(c("Simulation study ends at ",date()))

	return(rslt)
}


#debug(simulation2.carx)
rslt <- simulation2.carx(testPower = 0.1,nRep=1000)

options(verbose=F)
#delta <- seq(0.0,0.4,0.1)
delta <- seq(0.8,0.9,0.1)
power0 <- numeric(length(delta))
avgCR <- numeric(length(delta))
power <- numeric(length(delta))
#allRslt <- NULL
for(i in 1:length(delta) )
{
	rslt <- simulation2.carx(testPower = delta[i],nRep=1000)
	#allRslt <- c(allRslt, rslt)
	avgCR[i] <- rslt$averageCensorRate
	power0[i] <- rslt$rejectionRate0
	power[i] <- rslt$rejectionRate
}
print(delta)
print(avgCR)
print(power0)
print(power)

plot(delta,power)
lines(delta,power0)
