.First()
library(devtools)
load_all()
options(verbose=F)
message("testing diagnostics\n")

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
													 testPower=0,
													 debug=F)
{
	if(debug) message(c("Power study begins at ",date()))
	t0 <- proc.time()
	require(matrixStats)
	
	uppercl <- ucl
  
	#test if AR is stationary
	roots <- polyroot(c(1,-trueAR))
	if(debug)
	{
		message("roots of original AR")
		print(roots)
	}
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
  
	if(debug)
	{
		message("original arPrmtr")
		print(trueAR)
		message("disturbed arPrmtr")
		print(trueAR1)
		roots <- polyroot(c(1,-trueAR1))
		message("roots of polynomials of trueAR, see if it's stationary")
		print(roots)
	}

	#find censorLimit
	dat <- carx.simulate(10000, trueAR1, trueEV, trueSigma, -1, 1, seed=0)
	lowercl <- quantile(dat$y,targetCensorRate)

	simEst <-function(iRep)
	{
		ret <- NULL
		if(debug) message(sprintf("Rep:%i",iRep))
		for(n in nObs)
		{
			i = 1
			for(lcl in lowercl)
			{ 
				#message(sprintf("n:%f, cl:%f",n,lcl))
				dat <- carx.simulate(n, trueAR1, trueEV, trueSigma, lcl, uppercl, seed=37513*iRep)
				if(i == 1)
				{
					rslt0 <- carx(dat$y,dat$x,rep(0,n),NULL,NULL, nAR, getCI=fullEstimation)
					rsdl0 <- residuals(rslt0)
					for(l in nLag) ret <- c(ret, lbStat(rsdl0,l)>qchisq(alpha, l-nAR))
				  i = 2
				}
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
	if(debug)
	{
		message("rslt")
		print(rslt)
	}


	summaryList <- matrix(nrow=nRep,ncol=length(nObs)*6)
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
	if(debug) print(colm)
	if(debug) message(c("Simulation study ends at ",date())) 
	return(colm)
}


#debug(simulation2.carx)

singleTest <- F
if(singleTest)
{
	rslt <- simulation2.carx(testPower = 0.0,nRep=2)
	print(rslt)
	write.table(rslt,file="testRslt",row.names=F,col.names=F)
}

batch <- T
if(batch)
{
    delta <- seq(0.9,0.9,0.1)
    powermat <- matrix(nrow=length(delta),ncol = 1+12)
    powermat[,1] <- delta
    for(i in 1:length(delta) )
    { 
			message(sprintf("delta: %f",delta[i]))
			rslt <- simulation2.carx(testPower = delta[i],nRep=1000)
			print(rslt)
			powermat[i,-1] <- rslt
    }
    write.table(powermat,file="powerMat",row.names=F,col.names=F)
}
print(powermat)
