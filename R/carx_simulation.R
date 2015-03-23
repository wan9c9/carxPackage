#' provides a simulation study for \code{carx}.
#'
#' \code{simulation.carx} uses provided parameters and other settings to perform a simulation study.
#' Note that the exogenous variables are simulated from independent standard normal.
#' @param trueEV the true parameter (coefficients) for exogenous variable.
#' @param trueAR the true AR parameter.
#' @param lcl the lower censor limit.
#' @param ucl the upper censor limit.
#' @param nObs the sample size.
#' @param nRep the number of replications in the simulation study.
#' @param fullEstimation a bool value indicating if confidence interval by bootstrap 
#'        is performed.
#' @return If nRep is 1, then the returned value is a list of the \code{carx} object and
#' a summary of the simulated result. If nRep > 1, then it returns a list of (nRep, objects,
#' averageCensorRate, alpha, summary), where
#' nRep is the number of replication as is supplied;
#' objects is a vector of estimated \code{carx} objects;
#' averageCensorRate is the mean of censor rates over all replications;
#' alpha is the alpha value for the confidence interval;
#' summary is a summary matrix of the simulated results.
#'
#' @examples
#' simulation.carx(nRep=1) # nRep == 1
#' simulation.carx() #using default settings
#' simulation.carx(ucl=Inf) # no right censoring
#' simulation.carx(lcl=-Inf) # no left censoring

simulation.carx <-function(trueEV=c(0.2,0.4), 
													 trueAR=c(-0.5,0.3,-0.1), 
													 trueSigma=0.2,
													 lcl=-1,
													 ucl=1,
													 nObs=100,
													 nRep=2,
													 alpha=0.95,
													 nLag = 10,
													 fullEstimation=F)
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


	nAR <- length(trueAR)
	truePrmtr <- c(trueEV,trueAR,trueSigma)
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
		dat <- carx.simulate(nObs, trueAR, trueEV, trueSigma, lowercl, uppercl, seed=37513*iRep)
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
		#print("HH")
		rsdl <- residuals(rslt)
		#print(rsdl)
		rtest <- lbStat(rsdl,nLag)
		#print(rtest)
		ret <- c(ret,rtest)
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

	if(fullEstimation) 
		summaryList <- matrix(nrow=nRep,ncol=(3+5*nPrmtr))
	else
		summaryList <- matrix(nrow=nRep,ncol=(3+2*nPrmtr))
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

	rejectionRate <- sum(rslt[,dim(rslt)[2]] > qchisq(alpha, nLag-nAR))/nRep

	summary <- cbind(truePrmtr,avgEstd,std,coverageRate)
	colnames(summary) <- c('true', 'avg.estd','std.err','coverage.rate')
	xnames <- paste0('X',1:length(trueEV))
	rnames <- c(xnames,paste0('AR',1:nAR),"sigma")
	rownames(summary) <- rnames

	rslt <- list(nRep=nRep,
							 objects = objects,
							 averageCensorRate=avgCR,
							 alpha = alpha,
							 allResult = summaryList,
							 summary = summary)

	message(sprintf("\nReplication: %i\n", nRep))
	message(sprintf("\nAverage censor rate: %f\n", avgCR))
	message(sprintf("\n rejection rate: %f\n", rejectionRate))

	message("                     Simulation Summary        ")
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


