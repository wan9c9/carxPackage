#' returns the Ljung-Box test statistic for a series with a given lag
#'
#' \code{lbStat} returns the Ljung-Box test statistic for a series with a given lag.
#' @param res the series.
#' @param nLag number of lags used.
#' @return the Ljung-Box test statistic.
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




#' determine the goodness of fit of a fitted \code{carx} object
#'
#' \code{goodnessOfFit.carx} determines the goodness of fit of a fitted \code{carx} model.
#' @param object a fitted \code{carx} object.
#' @param nLag number of lags of the residuals used.
#' @param seed seed for random number generator.
#' @param bootstrap indicates whether to use the bootstrap to get the distribution of the test statistic.
#' @param nRep number of replications in bootstrap procedure.
#' @param nCores number of cores of CPU to be used in bootstrap, nCores>1 can speed up the computation. Default = \code{NULL}. See also \code{registerDoMC} in \pkg{doMC}.
#' @param ... not used.
#' @return a list of the test statistic \code{tVal} and its p-value \code{pVal}.
#' @import doParallel
#' @export
goodnessOfFit.carx  <- function(object,nLag=10, seed = NULL, bootstrap=FALSE, nRep=1000, nCores=NULL,...)
{

	if(!is.null(seed))
		set.seed(seed)
	#message("Calling goodnessOfFit.carx")
	res <- residuals(object)[-(1:object$p)]
	tStat <- lbStat(res,nLag)

	if(bootstrap)
	{
		#if(is.null(nCores))
			#doMC::registerDoMC()
		#else
		doMC::registerDoMC(nCores)
		#message(sprintf("forks: %i",foreach::getDoParWorkers()))
		simTStat <- foreach::foreach(i=1:nRep,.combine='c') %dopar%
		{
			message(sprintf("Diagnostics: nRep = %i",i))
			sim <- carx.sim(object$nObs, object$prmtrAR, object$prmtrX,object$sigma, object$lcl,object$ucl,x = object$x)
			obj <- carx(sim$y, sim$x, sim$ci, object$lcl, object$ucl, object$p, prmtrX=object$prmtrX,prmtrAR=object$prmtrAR,sigmaEps=object$sigma, getCI=FALSE)
			tmpRes <- residuals(obj)[-(1:object$p)]
			lbStat(tmpRes,nLag)
		}
		pVal <- sum(simTStat > tStat)/nRep
		return(list(tVal = tStat, pVal=pVal, bVal = simTStat))
	}
	else
	{
		return(list(tVal = tStat, pVal = pchisq(tStat,df = nLag - object$p )))
	}

}

#debug(goodnessOfFit.carx)
