lbStat <- function(res,nLag)
{
	n <- length(res)
	ac <- acf(res,lag.max=nLag,type="correlation",plot=FALSE)
	#print(ac)
	#print(ac[1:nLag])
	ac <- ac$acf[2:(nLag+1)]
	#print(ac)
	v <- seq.int(n-1,n-nLag)
	ret <- n*(n+2)*sum(ac^2/v)
	ret
}

lbStat2 <- function(res,nLag)
{
	n <- length(res)
	ac <- sapply(nLag, function(l) acf(res,lag.max=l,type="correlation",plot=FALSE)$acf[1:l])
	v <- seq(n-1,n-nLag,by=-1)
	ret <- n*(n+2)*sum(ac^2/v)
	ret
}




goodnessOfFit.carx  <- function(object,nLag=10, seed = NULL, bootstrap=FALSE, nRep=10, nCores=NULL,...)
{

	if(!is.null(seed))
		set.seed(seed)

	message("Calling goodnessOfFit.carx")
	res <- residuals(object)[-(1:object$nAR)]
	tStat <- lbStat(res,nLag)

	if(bootstrap)
	{
		require(foreach)
		require(doMC)
		if(is.null(nCores)) 
			registerDoMC()
		else
			registerDoMC(nCores)

		message(sprintf("forks: %i",getDoParWorkers()))

		simTStat <- foreach(i=1:nRep,.combine='c') %dopar%
		{
			message(sprintf("Diagnostics: nRep = %i",i))
			sim <- carx.simulate(object$nObs, object$prmtrAR, object$prmtrX,object$sigma,
													object$lowerCensorLimit,object$upperCensorLimit,x = object$x)
			obj <- carx(sim$y, sim$x, sim$censorIndicator, 
									object$lowerCensorLimit, object$upperCensorLimit, object$nAR, 
									prmtrX=object$prmtrX,prmtrAR=object$prmtrAR,sigmaEps=object$sigma,
									getCI=FALSE)
			tmpRes <- residuals(obj)[-(1:object$nAR)]
			lbStat(tmpRes,nLag)
		}
		pVal <- sum(simTStat > tStat)/nRep
		return(list(tVal = tStat, pVal=pVal, bVal = simTStat))
	}
	else
	{
		return(list(tVal = tStat, pVal = pchisq(tStat,df = nLag - object$nAR )))
	}
		
}

#debug(goodnessOfFit.carx)
