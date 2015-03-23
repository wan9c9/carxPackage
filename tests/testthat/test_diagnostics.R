.First()
load_all()
options(verbose=TRUE)
message("testing diagnostics")

simulation2.carx <-function(trueEV=c(0.2,0.4), 
													 trueAR=c(-0.28,0.25), 
													 trueSigma=0.60,
													 lcl=-1,
													 ucl=Inf,
													 nObs=c(100,200),
													 nRep=1000,
													 alpha=0.95,
													 nLag = c(10,20),
													 fullEstimation = FALSE,
													 targetCensorRate = c(0.15,0.3),
													 testPower=0)
{
	message(c("Power study begins at ",date()))
	t0 <- proc.time()
	require(matrixStats)
	
	uppercl <- ucl
  
	#test if AR is stationary
	roots <- polyroot(c(1,-trueAR))
	if(any(abs(roots)<=1))
	{
	  message("roots of polynomials of trueAR, see if it's stationary")
		message("supplied trueAR is not statoinary")
	  print(roots)
		return
	}
	
	nAR <- length(trueAR)
  trueAR1 <- numeric(nAR+1)

	if(testPower!=0) #perturb trueAR according to delta
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

	#roots <- polyroot(c(1,-trueAR1))
	#message("roots of polynomials of trueAR, see if it's stationary")
	#print(roots)
	#truePrmtr <- c(trueEV,trueAR,trueSigma)

	#find censorLimit
	dat <- carx.simulate(10000, trueAR1, trueEV, trueSigma, -1, 1, seed=0)
	lowercl <- quantile(dat$y,targetCensorRate)

	#' ( lower_ci, upper_ci, coverageIndicator)
	simEst <-function(iRep)
	{
		ret <- NULL
		message(sprintf("Rep:%i",iRep))
		for(n in nObs)
		{
			for(lcl in lowercl)
			{ 
				ret <- c(n,lcl)
				#message(sprintf("n:%f, cl:%f",n,lcl))
				dat <- carx.simulate(n, trueAR1, trueEV, trueSigma, lcl, uppercl, seed=37513*iRep)
				rslt0 <- carx(dat$y,dat$x,rep(0,n),NULL,NULL, nAR, getCI=fullEstimation)
				rsdl0 <- residuals(rslt0)
				for(l in nLag) ret <- c(ret, lbStat(rsdl0,l)>qchisq(alpha, l-nAR))

			  rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$lowerCensorLimit,dat$upperCensorLimit, nAR, getCI=fullEstimation)
				rsdl <- residuals(rslt)
				for(l in nLag) ret <- c(ret, lbStat(rsdl,l)>qchisq(alpha, l-nAR))
			}
		}
		ret
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
	#return(rslt)
	#message("rslt")
	#print(rslt)


	summaryList <- matrix(nrow=nRep,ncol=(length(nObs)*length(lowercl)*length(nLag)*2))
	i <- 1
	for(r in rslt) 
	{
		summaryList[i,] <- as.vector(r)
		i <- i+1
	}
	#print(summaryList)
	rslt <- summaryList
	colm <- colMeans(summaryList)
	message(sprintf("Time used:"))
	print(proc.time()-t0)
	message(c("Simulation study ends at ",date())) 
	return(colm)
}


#debug(simulation2.carx)
#rslt <- simulation2.carx(testPower = 0.0,nRep=1000)

delta <- seq(0.0,0.8,0.1)
powermat <- matrix(nrow=length(delta),ncol = 1+16)
powermat[,1] <- delta
for(i in 1:length(delta) )
{
	rslt <- simulation2.carx(testPower = delta[i],nRep=1000)
	powermat[i,-1] <- rslt
}

##########################################

	
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

