
#' @export
carx <- function(y,...) UseMethod("carx")


#' The default estimation method for CARX
#'
#' Use given data and other settings to estimate the parameters of CARX, and optionally compute the standard error or confidence intervals of parameter estimates by parametric bootstrap.
#' @param y a vector of possibly censored responses, only uncensored observation
#' marked by zeros of \code{ci} are used.
#' @param x a matrix of covariates, or some object which can be coerced to matrix.
#' @param ci a vector of -1,0,1's indicating that the corresponding y is
#' left-censored, not censored, and right-censored respectively. Default = NULL, indicating no censoring.
#' @param lcl a vector of lower censoring limits, or a number assuming constant limit, default = NULL, implying no lower censoring limit.
#' @param ucl a vector of upper censoring limits, or a number assuming constant limit, default = NULL, implying no upper censoring limit.
#' @param p the order of AR model for the regression errors, default = 1.
#' @param prmtrX the initial value for the parameter of \code{x}, default = \code{NULL}.
#' @param prmtrAR the initial value for the AR coefficients, default = \code{NULL}.
#' @param sigmaEps the initial value for the standard deviation of the errors of the AR process, default = \code{NULL}.
#' @param tol the tolerance in estimating the parameters, default = 1.0e-4.
#' @param max.iter maximum number of iterations allowed when estimating parameters, default = 500.
#' @param CI.compute bool value to indicate if the confidence interval for the
#' parameter is to be constructed, default = \code{FALSE}.
#' @param CI.level numeric value in (0,1) to get the \code{CI.level} confidence interval for the parameter, default = 0.95.
#' @param b number of bootstrap samples when estimating confidence interval for the parameter, default = 1000.
#' @param cenTS an optional argument to store \code{cenTS} object, which might be relevant to data used, used in \code{carx.formula}, default = NULL.
#' @param verbose bool value indicates whether to print intermediate information, default = FALSE.
#' @param ... not used.
#' @return a CARX object of the estimated model.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' model0 <- carx(y=dat$y, x=dat[,c("X1","X2")], ci=dat$ci, lcl=dat$lcl, ucl=dat$ucl, p=2)
#' #or simply call
#' model0 <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' plot(model0)
#' tsdiag(model0)

carx.default <- function(y,x,ci=NULL,lcl=NULL,ucl=NULL,
       p=1,prmtrX=NULL,prmtrAR=NULL,sigmaEps=NULL,
			 tol=1e-4,max.iter=500,CI.compute=FALSE,CI.level=0.95,b=1000,
			 cenTS=NULL,verbose=FALSE,...)
{
  #message("Enter carx.default.")
	verbose <- verbose || options()$verbose
	nObs <- length(y)

  finiteY <- which(is.finite(y))

	#standardize censoreIndicator
  if(is.null(ci))
    ci <- rep(0,length(y))
	ci[finiteY][ ci[finiteY]>0 ] <- 1; ci[finiteY][ ci[finiteY]<0 ] <- -1
	noCensor  <- all(ci[finiteY] == 0)

	if(is.null(lcl)) # no lower censoring
	{
		if(any(ci[finiteY]<0))
			stop("Error in data: lcl is null but there exist left-censored data.")
		else
			lcl = rep(-Inf,nObs)
	} else{
		if(length(lcl) == 1)
		{
			if(!is.nan(lcl))
				lcl <- rep(lcl,nObs)
			else
			{
				warning("lcl is NaN, I will set it to be -Inf.")
				lcl <- rep(-Inf,nObs)
			}
		}
		else #length > 1
		{
			if(length(lcl) != nObs)
				stop("Error: The dimesion of lower censor limit doesn't match that of y.")
		}
	}

	if(is.null(ucl)) #no upper censoring
	{
		if(any(ci[finiteY]>0)) #check
			stop("Error in data: ucl is null but there exist right-censored data.")
		else
			ucl = rep(Inf,nObs)
	} else{
		if(length(ucl) == 1)
		{
			if(!is.nan(ucl))
				ucl <- rep(ucl,nObs)
			else
			{
				warning("ucl is NaN, I will set it to be Inf.")
				ucl <- rep(Inf,nObs)
			}
		}
		else
		{
			if(length(ucl) != nObs)
				stop("Error: The dimesion of upper censor limit doesn't match that of y.")
		}
	}

  if(!is.null(x))
	{
		if(!is.matrix(x)) # x may be a vector
			x <- as.matrix(x)
		if(dim(x)[1] != nObs){
			stop(" The dimension of x doesn't match that of y.")
			return(NULL)
		}
		nX <- dim(x)[2]
		externalVariable <- x
		if(all(x==1))
			xIsOne <- TRUE
		else
			xIsOne <- FALSE
	}else
	{
		warning("x is null, I will set x = 1s, i.e., representing the intercept")
		nX <- 1
		externalVariable <- rep(1,nObs)
		xIsOne <- TRUE
	}

  #check for finite rows in data and construct skipIndex
  finiteRows <- is.finite(y) & (apply(externalVariable, 1, function(x){all(is.finite(x))}))
  skipIndex <- rep(0,length(y))
  skipIndex[1:p] <- 1
  for(i in 1:length(y))
  {
    if(!finiteRows[i])
      skipIndex[i:(i+p)] <- 1
  }
  skipIndex <- which(skipIndex > 0)
  if(verbose)
  {
    message("Skip index is constructed as ")
    message(paste(skipIndex,sep=', '))
  }
	nSkip <- length(skipIndex)

	if(any(lcl[-skipIndex] >= ucl[-skipIndex]))
		stop("Error in censor limit: some lower censor limits are bigger than upper censor limits.")

	ret = list(y = y,
		   x = externalVariable,
		   xIsOne = xIsOne,
		   ci = ci,
		   lcl = lcl,
		   ucl = ucl,
       CI.compute = CI.compute,
       tol = tol,
       max.iter = max.iter,
       CI.level = CI.level,
       b = b,
		   skipIndex = skipIndex,
       finiteRows = finiteRows,
       cenTS=cenTS,
       verbose = verbose
		   )
	#print(ret)

	#parameters
	#generic
	if(is.null(prmtrX) || is.null(prmtrAR) || is.null(sigmaEps) )
	{
		prmtrX <- numeric(nX)
		prmtrAR <- numeric(p)
		sigmaEps <- numeric(1)
		needInit <- TRUE
	}
	else
		needInit <- FALSE

	#special to store estimated values
	prmtrXEstd <- numeric(nX)
	prmtrAREstd <- numeric(p)
	sigmaEpsEstd <- numeric(1)

	trend <- numeric(nObs)
	covEta <- matrix(nrow=p+1, ncol = p+1)
	wkMean <- matrix(nrow=nObs, ncol = p+1)
	wkCov <- array(0, dim=c(nObs,p+1,p+1))
	res <- numeric(nObs)

  getNPrmtr <- function(){
		return (nX + p +  1)
	}

  pVal <- numeric(getNPrmtr())


	# working space
	resetWK <- function()
	{
		trend <<- numeric(nObs)
		covEta <<- matrix(nrow=p+1, ncol = p+1)
		wkMean <<- matrix(nrow=nObs, ncol = p+1)
		wkCov <<- array(0, dim=c(nObs,p+1,p+1))
		res <<- numeric(nObs)
    pVal <<- numeric(getNPrmtr())

		#initialize wkMean
		for(idx in seq(1,nObs)[-skipIndex])
			wkMean[idx,] <<- y[idx:(idx-p)]
	}
	#-------------------------------


	getPrmtr <- function(){
		ret <- c(prmtrX, prmtrAR, sigmaEps)
		return(ret)
	}

	getCensorRate <- function(){
		return(sum(abs(ci[finiteRows]))*1.0/sum(finiteRows))
	}


	setInitPrmtrForBootstrap <- function()
	{ # set current parameter to be the estimated
		prmtrX   <<-   prmtrXEstd
		prmtrAR   <<-	prmtrAREstd
		sigmaEps  <<-	sigmaEpsEstd
	}

	setEstdPrmtr <- function()
	{ # store the estimated parameter
		prmtrXEstd   <<- prmtrX
		prmtrAREstd   <<- prmtrAR
		sigmaEpsEstd  <<- sigmaEps
	}

	getEstdPrmtr <- function(){
		ret <- c(prmtrXEstd, prmtrAREstd, sigmaEpsEstd)
		return(ret)
	}

	updateTrend <- function(){
		trend <<- externalVariable%*%prmtrX
	}

	updateCovEta <- function(){
		covEta <<- computeCovAR(prmtrAR, sigmaEps)
	}

	#calculate the expectations
	eStep <- function()
	{
		#if(verbose) message("E-step begin.")
		# expection step
		# need to 1) calculate the joint distriubtion of AR process
		#         2) update the trend function/observations
		#         3) calculate the conditional expection & covariances of the censored observations
		updateCovEta()
		updateTrend()
		for(idx in seq(1,nObs)[-skipIndex])
		{
			tmpCensorIndicator <- ci[idx:(idx-p)]
			nCensored <- sum(abs(tmpCensorIndicator))
			if(nCensored>0)
			{# at least one is censored
				if( nCensored < (p+1) )
				{
					conditionalIndex <- which(tmpCensorIndicator==0) #indices for observed data
					tmpY <- y[idx:(idx-p)][conditionalIndex]
					tmpM <- trend[idx:(idx-p)]
					cdist <- conditionalDistMvnorm(tmpY,conditionalIndex,tmpM,covEta)
					tmpMean <- cdist$'mean'
					tmpVar <- cdist$'var'
				}else{
					tmpMean <- trend[idx:(idx-p)]
					tmpVar <- covEta
				}
				tmpLower <- rep(-Inf,length = nCensored)
				tmpUpper <- rep(Inf,length = nCensored)
				censored <- tmpCensorIndicator[tmpCensorIndicator!=0]

				# lower limit is upper censor limit
				tmpLower[censored>0] <- ucl[idx:(idx-p)][tmpCensorIndicator>0]
				# upper limit is lower censor limit
				tmpUpper[censored<0] <- lcl[idx:(idx-p)][tmpCensorIndicator<0]
				#if(abs(det(tmpVar)) < 0.001) warning("covariance matrix is nearly singluar.")
				ret <- tmvtnorm::mtmvnorm(tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)
				wkMean[idx,tmpCensorIndicator!=0] <<- ret$'tmean'
				wkCov[idx,tmpCensorIndicator!=0,tmpCensorIndicator!=0] <<- ret$'tvar'
			}
			#if(verbose) message("E-step done.")
		}
	}

	eStepNaive <- function()
	{
		if(verbose) message("entering eStepNaive")
		if(!noCensor)
		{
			for(idx in seq(1,nObs)[-skipIndex])
			{
				tmpCensorIndicator <- ci[idx:(idx-p)]
				wkMean[idx,tmpCensorIndicator>0] <<- ucl[idx:(idx-p)][tmpCensorIndicator>0]
				wkMean[idx,tmpCensorIndicator<0] <<- lcl[idx:(idx-p)][tmpCensorIndicator<0]
				wkCov[idx,,] <<- 0
			}
		}
		if(verbose) message("E-step Naive done.")
	}

	updatePrmtrAR <- function(){
		wkMean2 <- wkMean
		for(i in 1:(p+1)){
			wkMean2[(p+1):nObs,i] <- wkMean2[(p+1):nObs,i] - trend[(p+2-i):(nObs+1-i)]
		}
		mat <- matrix(0,p,p)
		vec <- numeric(p)
		for(idx in seq(1,nObs)[-skipIndex]){
			vec <- vec + wkMean2[idx,1]*wkMean2[idx,-1] + wkCov[idx,-1,1]
			mat <- mat + outer(wkMean2[idx,-1],wkMean2[idx,-1]) + wkCov[idx,-1,-1]
		}

		tmp <- solve(mat,vec)
		if(any(abs(polyroot(c(1,-tmp)))<=1.0)) #tmp is not feasible
		{
			grad <- mat%*%prmtrAR - vec
			tmp <- prmtrAR - grad
			i <- 1
			while(any(abs(polyroot(c(1,-tmp)))<=1.0) && i <= 10)
			{
				tmp <- prmtrAR - grad/2^i
				i <- i+1
			}
			tmp <- as.vector(tmp)
		}
		return(tmp)
	}

	updatePrmtrEV <- function()
	{
		if(p > 1){
			tmpY <- wkMean[,1] - wkMean[,-1]%*%prmtrAR
		}
		else{
			tmpY <- wkMean[,1] - wkMean[,-1]*prmtrAR
		}
		tmpY <- tmpY[-skipIndex]

		tmpX <- matrix(0,nObs-nSkip,nX)
		index <- seq(1,nObs)[-skipIndex]
		for(idx in 1:(nObs-nSkip)){
			if (p >1){
				tmpX[idx,] <- externalVariable[(index[idx]),] - prmtrAR %*% externalVariable[(index[idx]-1):(index[idx]-p),]
			}
			else{
				tmpX[idx,] <- externalVariable[(index[idx]),] - prmtrAR * externalVariable[(index[idx]-1):(index[idx]-p),]
			}
		}
		ret <- solve(t(tmpX) %*% tmpX, t(tmpX)%*%tmpY )

		return(ret)
	}

	updateSigmaEps <- function(){
		idx <- seq(1,nObs)[-skipIndex]
		vec <- wkMean[idx,1] - trend[idx]
		for(i in 1:p){
			vec <- vec - prmtrAR[i]*(wkMean[idx,i+1] - trend[(idx-i)])
		}
		ret <- sum(vec^2)
		tmpVec <- c(1,-prmtrAR)
		for(i in idx){
			ret <- ret + as.numeric( tmpVec%*%wkCov[i,,]%*%tmpVec )
		}
		return( sqrt(ret/(nObs-nSkip)) )
	}

	mStep <- function(){
    prevPrmtr <- c(prmtrX,prmtrAR,sigmaEps)
		delta <- 0

		newPrmtrAR <- updatePrmtrAR()
		delta <- delta + sum(abs(newPrmtrAR - prmtrAR))
		prmtrAR <<- newPrmtrAR

		newPrmtrEV <- updatePrmtrEV()
		delta <- delta + sum(abs(newPrmtrEV - prmtrX))
		prmtrX <<- newPrmtrEV

		newSigmaEps <- updateSigmaEps()
		delta <- delta + abs(newSigmaEps - sigmaEps)
		sigmaEps <<- newSigmaEps
		
    return(delta/sqrt(sum(prevPrmtr)^2))
	}

	estimatePrmtr <- function(tol,max.iter)
	{
		if(verbose) message("estimating parameters")
		delta <- 1.0
		nIter <- 1
		while(delta > tol && nIter < max.iter)
		{
			eStep()
			delta <- mStep()
			if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s", nIter, delta, toString(round(getPrmtr(),digits = 4)) ) )
			nIter <- nIter + 1
		}
		if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s \nParameter estimated\n", nIter-1, delta, toString(round(getPrmtr(),digits=4)) ) )
	}

	initPrmtr <- function(tol,max.iter)
	{
		if(verbose) message("Initializing parameters")
		delta <- 1.0
		nIter <- 1

		prmtrX <<- numeric(nX)
		prmtrAR <<- numeric(p)
		sigmaEps <<- 1

		eStepNaive()
		while(delta > tol && nIter < max.iter)
		{
			updateTrend() # needed since it is done to update the trend which is done in the eStep for our method
			delta <- mStep()
			if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s", nIter, delta, toString(round(getPrmtr(),digits=4) ) ))
			nIter <- nIter + 1
		}
		if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s \nParameter initialized\n", nIter-1, delta, toString(round(getPrmtr(),digits=4) ) ))
	}

	exptdLogLik <- function()
	{
		val <- -(nObs - nSkip)*(1+log(sigmaEpsEstd^2))/2
		val
	}

	bootstrapSample <- function()
	{
		eps <- stats::rnorm(nObs,0, sigmaEps)

		updateTrend()

		eta <- numeric(nObs)
		eta[skipIndex] <- eps[skipIndex]
		for(i in seq(1,nObs)[-skipIndex]){
			eta[i] <- eta[(i-1):(i-p)] %*% prmtrAR + eps[i]
			y[i] <<- trend[i] + eta[i]
		}
		ci <<- rep(0,nObs)
		ci[y<lcl] <<- -1
		ci[y>ucl] <<- 1
		if(verbose) message(paste0("\ncensor rate: ", sum(abs(ci))/nObs))
	}

	bootstrapCI <- function(CI.level,b)
	{
		#message(sprintf('Bootstraping CI'))
		yOriginal <- y #copy y as it will be overwritten
		tmpResult <- matrix(0,b,getNPrmtr())
		#pb <- txtProgressBar(1,b,style=3)
		for(i in 1:b){
			#message(sprintf('Bootstraping CI %i/%i',i,b))
			#setTxtProgressBar(pb,i)
			setInitPrmtrForBootstrap()
			bootstrapSample()
			resetWK()
			estimatePrmtr(tol,max.iter)
			tmpResult[i,] <- getPrmtr()
		}
		#close(pb)
		qv <- c((1-CI.level)/2,1-(1-CI.level)/2)
		for(j in 1:(getNPrmtr())){
			CI[j,] <<- stats::quantile(tmpResult[,j],qv)
      pval <- mean(tmpResult[,j]>0)
      pVal[j] <<- 2*min(pval,1-pval)
		}
		covMat <<- stats::cov(tmpResult)
		y <<- yOriginal #set back observed data
	}

	# begin execution
	resetWK()
	#message("initializing parameters")
	if(needInit || noCensor )  # If no censoring, it is estimating an AR-X
		initPrmtr(tol, max.iter)
	prmtrInit <- getPrmtr()
	if(!noCensor)
		estimatePrmtr(tol,max.iter)
	setEstdPrmtr()

	coeff <- c(prmtrXEstd,prmtrAREstd)
	xnames <- colnames(x)
	if(is.null(xnames))
		xnames <- paste0('X',1:dim(x)[2])
	names(coeff) <- c(xnames,paste0('AR',1:p))

	ret$coefficients = coeff
	ret$prmtrX = prmtrXEstd
	ret$prmtrAR = prmtrAREstd
	ret$sigma = sigmaEpsEstd

	ret$censorRate = getCensorRate()
	rnames <- c(names(coeff),"sigma")
	ret$prmtrInit = prmtrInit
	ret$prmtrEstd = getEstdPrmtr()
	names(ret$prmtrInit) <- rnames
	names(ret$prmtrEstd) <- rnames
	ret$nObs = nObs
	ret$logLik = exptdLogLik()
	ret$p = p
	ret$nX = nX
	ret$npar = getNPrmtr()
	ret$aic = 2*(-ret$logLik + ret$npar)
	ret$call = match.call()

	if(CI.compute){
		#rnames <- c(names(coeff),"sigma")
		CI <- matrix(nrow=getNPrmtr(),ncol=2)
		covMat <- matrix(nrow=getNPrmtr(),ncol=getNPrmtr())
		bootstrapCI(CI.level,b)
		rownames(CI) <- rnames
		colnames(CI) <- c(sprintf("%4.2f%%",100*(1-CI.level)/2), sprintf("%4.2f%%",100*(1+CI.level)/2))
		rownames(covMat) <- rnames
		colnames(covMat) <- rnames
		ret$CI <- CI
		ret$vcov <- covMat
    ret$pVal <- pVal
	}else{
		ret$CI <- NULL
		ret$vcov <- NULL
    ret$pVal <- NULL
	}
	class(ret) <- "carx"
  #message("Exit carx.default.")
	ret
}

#' Provide a simple formula interface to the \code{carx} method
#'
#' This interface will use the supplied \code{formula} and data provided by \code{data} and other arguments in \code{...} to invoke the \code{carx.default} method.
#' @param formula the formula.
#' @param data the data, can be a \code{list}, \code{data.frame}, or a \code{cenTS} object.
#' @param ... other arguments supplied to \code{\link{carx.default}}.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' model0 <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)

carx.formula <- function(formula, data=list(),...)
{
  #message("Enter carx.formula")
  vars <- list(...)
  nvars <- names(vars)

  if('cenTS' %in% class(data))
  {
    data2 <- data.frame(zoo::coredata(data))
    vars$cenTS <- data
  }
  else
    data2 <- data

	mf <- stats::model.frame(formula=formula,data=data2,na.action=NULL)
	y <- stats::model.response(mf)
	x <- stats::model.matrix(attr(mf,"terms"),data=mf)

  if( "lcl" %in% names(data2) )
    vars$lcl <- data2$lcl

  if( "ucl" %in% names(data2) )
    vars$ucl <- data2$ucl

  if( "ci" %in% names(data2) )
    vars$ci <- data2$ci

  toPass <- c(list(y=y,x=x),vars)
	#est <- carx.default(y,x,...)
  est <- do.call(carx.default,toPass)
	est$call <- match.call()
  #print(est$call)
	est$formula <- formula
  est$data <- data
  #message("Exit carx.formula")
	est
}

#' Return the quasi-log-likelihood of a \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param ... not used.
#' @return the quasi-log-likelihood.
#' @export
logLik.carx <- function(object,...)
{
	ret <- object$logLik
	class(ret) <- 'logLik.carx'
	ret
}

#' The AIC of a fitted  \code{carx} object
#'
#' Return the AIC of a fitted \code{carx} object based on the quasi-log-likelihood.
#' @param object a fitted  \code{carx} object.
#' @param ... not used.
#' @param k penalty for the number of parameters, default = 2.
#' @return the AIC value
#' @export
AIC.carx <- function(object,...,k=2)
{
	val <- -2*object$logLik + k*object$npar
	#class(val) <- "AIC.carx"
	val
}



#' Print a short description of the fitted model
#' @param x a fitted model object.
#' @param ... not used.
#' @return none.
#' @export
print.carx <- function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(x$coefficients)

	cat("\nResidual (innovation) standard deviation:\n")
	print(x$sigma)

	cat("\nCensoring rate:\n")
	print(x$censorRate)
	cat("\nSample size:\n")
	print(x$nObs)

	cat("\nNumber of parameters:\n")
	print(x$npar)

	cat("\nQuasi-log-likelihood:\n")
	print(x$logLik)

	cat("\nAIC:\n")
	print(x$aic)


	#cat("\nInitial estimates:\n")
	#print(x$prmtrInit)
	#cat("\nFinal estimates:\n")
	#print(x$prmtrEstd)
	if(x$CI.compute){
		cat(paste0("\nConfidence interval:\n"))
		print(x$CI)
		cat("\nVariance-covariance matrix:\n")
		print(x$vcov)
    cat(paste0("N.B.: Confidence intervals and variance-covariance matrix are based on ", x$b, " bootstrap samples.\n"))
	}
}

#' Summarize the fitted \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param ... not used.
#' @return a summary.
#' @export
summary.carx <- function(object,...)
{
	numDig <- function(x,d1,d2){
		y <- x
		y[abs(x) > 0.1^d1] <- round(x[abs(x) > 0.1^d1],d1)
		y[abs(x) < 0.1^d1] <- round(x[abs(x) < 0.1^d1],d2)
		y
	}

	if(is.null(object$CI))
  {
		tab <- cbind(Estimate = object$prmtrEstd)
	}else
  {
		se <- sqrt(diag(object$vcov))
		tVal <- c(coef(object),object$sigma)/se

		est <- coef(object)
		tab <- cbind(Estimate = est,
			     StdErr =  se[1:(object$npar-1)],
			     lowerCI = object$CI[1:(object$npar-1),1],
			     upperCI = object$CI[1:(object$npar-1),2],
           p.value = object$pVal[1:(object$npar-1)]
			     )
	}
	res <- list(call=object$call,
		    coefficients=tab,
		    aic = AIC.carx(object))
	class(res) <- "summary.carx"
	res
}

#' Print the summary of an \code{carx} object
#' @param x a summary of an \code{carx} object.
#' @param ... not used.
#' @return none.
#' @export
print.summary.carx <- function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")

	cat("\nCoefficients:\n")
	#print(x$coefficients)
  pval <- 'p.value' %in% colnames(x$coefficients)
  stats::printCoefmat(x$coefficients, P.values = pval, has.Pvalue = pval)
 cat(sprintf("Note: confidence intervals are based on %i bootstraps.  P.values are one-sided.\n", x$b))

	cat("\nAIC:\n")
	print(x$aic)
}


