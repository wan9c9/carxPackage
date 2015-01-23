#' \code{carx}: A package for estimating the parameters of Censored Auto-Regressive model with eXogenous covariates (CARX).
#'
#' \code{carx}: A package for estimating the parameters of Censored Auto-Regressive model with eXogenous covariates (CARX), 
#' , or it can be thought as regression models with possibly censored responses and autoregressive residuals
#' \code{carx} allows left-censoring, right-censoring or double-censoring of the response variable \code{y},
#' where whether the corresponding \code{y_t} is censored or not is recorded in the variable \code{censorIndicator}, 
#' which takes value -1, 0, and 1 if it is left-censored, not censored, and right-censored, respectively.
#' The censoring limit is recorded in the variable \code{censorLimit}.
#'
#' @docType package
#' @name carx
NULL


library(stats)
library(tmvtnorm)
#library(debug)
library(ggplot2)
#library(R2HTML)
#library(zoo)
#library(matrixStats)


#' \code{conditionalDistMvnorm} calculates the conditional mean & variance of a random vector following multivariate normal distribution.
#'
#' \code{conditionalDistMvnorm} calculates the conditional mean & variance of a sub-vector
#' of a vector given the rest of known elements, following multivariate normal distribution.
#' This function calculates the conditional mean & variance of a multivariate normal distribution.
#' with (meanVec, varMat), conditional on y at indices conditionalIndex.
#'
#' @param y a vector to be conditioned on
#' @param conditionalIndex the index to be conditioned on
#' @param meanVec the mean vector the joint multivariate normal distribution
#' @param varMat the variance-covariance matrix of the joint multivariate normal distribution
#' @return a list consisting of 'mean' and 'var' representing the conditional mean and variance respectively.
#' @keywords conditional distribution multivariate normal
#' @export
#' @examples
#' conditionalDistMvnorm(c(-0.5,0.5), c(2,4), c(1,2,3,4),matrix(
#' c(1,0.3,0.2,0.1, 0.3,1,-0.1,0.3,0.2,-0.1,1,0.1,0.1,0.3,0.1,1),
#' nrow=4,ncol=4,byrow=TRUE))

conditionalDistMvnorm <- function(y, conditionalIndex, meanVec, varMat)
{
	  sigma11 <- varMat[-conditionalIndex,-conditionalIndex]
	  sigma22 <- varMat[ conditionalIndex, conditionalIndex]
	  sigma12 <- varMat[-conditionalIndex, conditionalIndex]
	  sigma21 <- varMat[ conditionalIndex,-conditionalIndex]
	  invSigma22 <- solve(sigma22)
	  mNew <- meanVec[-conditionalIndex] + as.vector(sigma12 %*% invSigma22 %*%(y - meanVec[conditionalIndex]))
	  vNew <- sigma11 - sigma12 %*% invSigma22 %*% sigma21

	  return (list('mean'=mNew,'var' = vNew))
}

#' compute the covariance matrix of {\eta_t,\cdots,\eta_{t-p}} of the AR model
#' @param order the order of the AR model
#' @param sigmaEps the standard deviation of the residuals of the AR model
#' @return the covariance matrix needed
computeCovAR <- function(order, sigmaEps)
{
	  n <- length(order)+1
	  val <- ARMAacf(ar=order, lag.max=n)
	  val <- as.vector(val)

	  mat <- matrix(nrow=n,ncol=n)
	  for(i in 1:n)
	  {
			mat[i,i] <- 1
			if(i > 1){
				  for(j in 1:(i-1))
				  {
						mat[i,j] <- val[abs(i-j)+1]
						mat[j,i] <- mat[i,j]
				  }
			}
	  }
	  v <- sigmaEps^2/(1-order%*%val[2:n])
	  mat <- v[1,1]*mat
	  return(mat)
}


carx <- function(x,...) UseMethod("carx")

#' \code{carx.default} is the default method for CARX.

#' \code{carx.default} uses given data and other settings to compute the estimated parameters of CARX,
#' and optionally compute the standard error or confidence intervals of parameter estimates
#' by parametric bootstrap.

#' @param y a vector of regressors
#' @param x a matrix of covariances
#' @param censorIndicator a vector of -1,0,1's indicating that the corresponding y is
#' left-censored, not censored, or right-censored.
#' @param censorLimit a vector of censor limits for each y
#' @param lowerCensorLimit a vector of lower censor limits for each y or a number assuming uniform limit
#' @param upperCensorLimit a vector of upper censor limits for each y or a number assuming uniform limit
#' @param nAR the order of AR model for the regression errors
#' @param tol the tolerance in estimating the parameters, default = 1.0e-4
#' @param max.iter maximum number of iterations allowed when estimating parameters, default = 500
#' @param getCI bool value to indicate if the confidence interval for the parameter is needed.
#' @param alpha numeric value (0,1) to get the confidence interval for the parameter
#' @param nBootstrapSample number of bootstrap samples when estimating confidence interval for the parameter, default = 1000
#' @param useGoodRes bool value to indicate if use estimated residuals to bootstrap the confidence interval, default = FALSE
#' @param skipIndex a vector of indices indicating indices to be skipped, as calculating the conditional log-likelihood need some initial values to start, also is the initial values to calculate the conditional log-likelihood, useful if there are multiple segment of series in the whole series.
#' @return a CARX object of the estimated model


carx.default <- function(y,x,censorIndicator,lowerCensorLimit,upperCensorLimit,nAR,
						 tol=1e-4,max.iter=100,getCI=FALSE,alpha=0.95,nBootstrapSample=1000,
						 useGoodRes=FALSE,skipIndex=NULL,verbose=FALSE)
{
	  verbose <- verbose || options()$verbose
	  nObs <- length(y)
	  
	  censorIndicator[ censorIndicator>0 ] <- 1
	  censorIndicator[ censorIndicator<0 ] <- -1

	  if(is.null(lowerCensorLimit))
	  {
			if(any(censorIndicator<0))
				  stop("Error in data: lowerCensorLimit is null but there exist left-censored data.")
			else
				  lowerCensorLimit = rep(-Inf,nObs)
	  } else{
			if(length(lowerCensorLimit) == 1)
			{ 
				  if(!is.nan(lowerCensorLimit)) 
						lowerCensorLimit <- rep(lowerCensorLimit,nObs)
				  else
				  {
						warning("lowerCensorLimit is NaN, I will set it to be -Inf.")
						lowerCensorLimit <- rep(-Inf,nObs)
				  }
			}
			else
			{
				  if(length(lowerCensorLimit) != nObs)
						stop("Error: The dimesion of lower censor limit doesn't match that of y.")
			}
	 }

	  if(is.null(upperCensorLimit))
	  {
			if(any(censorIndicator>0))
				  stop("Error in data: upperCensorLimit is null but there exist right-censored data.")
			else
				  upperCensorLimit = rep(Inf,nObs)
	  } else{
			if(length(upperCensorLimit) == 1)
			{
				  if(!is.nan(upperCensorLimit))
						upperCensorLimit <- rep(upperCensorLimit,nObs) 
				  else 
				  {
						warning("upperCensorLimit is NaN, I will set it to be Inf.")
						upperCensorLimit <- rep(Inf,nObs)
				  }
			}
			else
			{
				  if(length(upperCensorLimit) != nObs)
						stop("Error: The dimesion of upper censor limit doesn't match that of y.")
			}
	  }
	  if(any(lowerCensorLimit > upperCensorLimit))
			stop("Error in censor limit: some lower censor limit is bigger than upper censor limit.")
      
	  if(!is.null(x))
	  {
			nEV <- dim(x)[2]
			if(dim(x)[1] != nObs){
				  stop(" The dimension of x doesn't match that of y.")
				  return(NULL)
			}
	        externalVariable <- x
	  }else
	  {
			nEV <- 1
			externalVariable <- rep(1,nObs)
	  }

      if(is.null(skipIndex))
	  {
			skipIndex <- seq(1,nAR)
			if(verbose) 
			{
				  message("skip index is constructed as ") 
				  message(skipIndex)
			}
	  }
	  nSkip <- length(skipIndex)

	  ret = list(y = y,
				 x = externalVariable,
				 censorIndicator = censorIndicator,
				 lowerCensorLimit = lowerCensorLimit,
				 upperCensorLimit = upperCensorLimit,
				 skipIndex = skipIndex
				 )

	  #parameters
	  #generic
	  prmtrAR <- numeric(nAR)
	  sigmaEps <- numeric(1)
	  prmtrEV <- numeric(nEV)

	  #special to store estimated values
	  prmtrAREstd <- numeric(nAR)
	  sigmaEpsEstd <- numeric(1)
	  prmtrEVEstd <- numeric(nEV)



	  # working data
	  resetWK <- function()
	  {
			trend <<- numeric(nObs)
			covEta <<- matrix(nrow=nAR+1, ncol = nAR+1)
			wkMean <<- matrix(nrow=nObs, ncol = nAR+1)
			wkCov <<- array(0, dim=c(nObs,nAR+1,nAR+1))
			res <<- numeric(nObs)

			#initialize wkMean
			for(idx in seq(1,nObs)[-skipIndex])
				  wkMean[idx,] <<- y[idx:(idx-nAR)] 
	  }
	  #-------------------------------

	  getNPrmtr <- function(){
			return (nEV + nAR +  1)
	  }

	  getPrmtr <- function(){
			ret <- c(prmtrEV, prmtrAR, sigmaEps)
			return(ret)
	  }

	  getCensorRate <- function(){
			return(sum(abs(censorIndicator))*1.0/nObs)
	  }


	  setInitPrmtrForBootstrap <- function()
	  {
			prmtrEV   <<-   prmtrEVEstd
			prmtrAR   <<-	prmtrAREstd
			sigmaEps  <<-	sigmaEpsEstd
	  }

	  setEstdPrmtr <- function()
	  {
			prmtrEVEstd   <<- prmtrEV
			prmtrAREstd   <<- prmtrAR
			sigmaEpsEstd  <<- sigmaEps
	  }

	  getEstdPrmtr <- function(){
			ret <- c(prmtrEVEstd, prmtrAREstd, sigmaEpsEstd)
			return(ret)
	  }


	  updateTrend <- function(){
			trend <<- externalVariable%*%prmtrEV
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
			#         3) calculat the expection & covariances of the censored observations
			updateCovEta()
			updateTrend()
			for(idx in seq(1,nObs)[-skipIndex])
			{
				  tmpCensorIndicator <- censorIndicator[idx:(idx-nAR)]
				  nCensored <- sum(abs(tmpCensorIndicator))
				  if(nCensored>0)
				  {# at least one is censored
						if( nCensored < (nAR+1) )
						{
							  conditionalIndex <- which(tmpCensorIndicator==0) #indices for observed data
							  tmpY <- y[idx:(idx-nAR)][conditionalIndex]
							  tmpM <- trend[idx:(idx-nAR)]
							  cdist <- conditionalDistMvnorm(tmpY,conditionalIndex,tmpM,covEta)
							  tmpMean <- cdist$'mean'
							  tmpVar <- cdist$'var'
						}else{
							  tmpMean <- trend[idx:(idx-nAR)]
							  tmpVar <- covEta
						}
						tmpLower <- rep(-Inf,length = nCensored)
						tmpUpper <- rep(Inf,length = nCensored)
						censored <- tmpCensorIndicator[tmpCensorIndicator!=0]

						tmpLower[censored>0] <- upperCensorLimit[idx:(idx-nAR)][tmpCensorIndicator>0]
						tmpUpper[censored<0] <- lowerCensorLimit[idx:(idx-nAR)][tmpCensorIndicator<0]
						ret <- mtmvnorm(tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)

						wkMean[idx,tmpCensorIndicator!=0] <<- ret$'tmean'
						wkCov[idx,tmpCensorIndicator!=0,tmpCensorIndicator!=0] <<- ret$'tvar'
				  }
				  #if(verbose) message("E-step done.")
			}
	  }

	  eStepNaive <- function()
	  {
			if(verbose) message("entering eStepNaive")
			for(idx in seq(1,nObs)[-skipIndex])
			{
				  tmpCensorIndicator <- censorIndicator[idx:(idx-nAR)]
				  wkMean[idx,tmpCensorIndicator>0] <<- upperCensorLimit[idx:(idx-nAR)][tmpCensorIndicator>0]
				  wkMean[idx,tmpCensorIndicator<0] <<- lowerCensorLimit[idx:(idx-nAR)][tmpCensorIndicator<0]
				  wkCov[idx,,] <<- 0
			}
			#if(verbose) message("E-step Naive done.")
	  }

	  setResiduals <- function()
	  {
			res <<- numeric(nObs)
			res[skipIndex] <<- NaN
			res[-skipIndex] <<- y[-skipIndex] - fittedValues[-skipIndex]
			res[censorIndicator>0] <<- upperCensorLimit[censorIndicator>0] - fittedValues[censorIndicator>0]
			res[censorIndicator<0] <<- lowerCensorLimit[censorIndicator<0] - fittedValues[censorIndicator<0]
	  }

	  
	  getResiduals <- function(){
			return(res)
	  }

	  updatePrmtrAR <- function(){
			wkMean2 <- wkMean
			for(i in 1:(nAR+1)){
				  wkMean2[(nAR+1):nObs,i] <- wkMean2[(nAR+1):nObs,i] - trend[(nAR+2-i):(nObs+1-i)]
			}
			mat <- matrix(0,nAR,nAR)
			vec <- numeric(nAR)
			for(idx in seq(1,nObs)[-skipIndex]){
				  vec <- vec + wkMean2[idx,1]*wkMean2[idx,-1] + wkCov[idx,-1,1]
				  mat <- mat + outer(wkMean2[idx,-1],wkMean2[idx,-1]) + wkCov[idx,-1,-1]
			}
			return( solve(mat,vec) )
	  }

	  updatePrmtrEV <- function(){
			if(nAR > 1){
				  tmpY <- wkMean[,1] - wkMean[,-1]%*%prmtrAR
			}
			else{
				  tmpY <- wkMean[,1] - wkMean[,-1]*prmtrAR
			}
			tmpY <- tmpY[-skipIndex]

			tmpX <- matrix(0,nObs-nSkip,nEV)
			index <- seq(1,nObs)[-skipIndex]
			for(idx in 1:(nObs-nSkip)){
				  if (nAR >1){
						tmpX[idx,] <- externalVariable[(index[idx]),] - prmtrAR %*% externalVariable[(index[idx]-1):(index[idx]-nAR),]
				  }
				  else{
						tmpX[idx,] <- externalVariable[(index[idx]),] - prmtrAR * externalVariable[(index[idx]-1):(index[idx]-nAR),]
				  }
			}
			ret <- solve(t(tmpX) %*% tmpX, t(tmpX)%*%tmpY )

			return(ret)
	  }

	  updateSigmaEps <- function(){
			idx <- seq(1,nObs)[-skipIndex]
			vec <- wkMean[idx,1] - trend[idx]
			for(i in 1:nAR){
				  vec <- vec - prmtrAR[i]*(wkMean[idx,i+1] - trend[(idx-i)])
			}
			ret <- sum(vec^2)
			tmpVec <- c(1,-prmtrAR)
			for(i in idx){
				  ret <- ret + as.numeric( tmpVec%*%wkCov[i,,]%*%tmpVec )
			}
			return( sqrt(ret/(nObs-length(skipIndex))) )
	  }

	  mStep <- function(){
			delta <- 0

			newPrmtrAR <- updatePrmtrAR()
			delta <- delta + sum(abs(newPrmtrAR - prmtrAR))
			prmtrAR <<- newPrmtrAR

			newPrmtrEV <- updatePrmtrEV()
			delta <- delta + sum(abs(newPrmtrEV - prmtrEV))
			prmtrEV <<- newPrmtrEV

			newSigmaEps <- updateSigmaEps()
			delta <- delta + abs(newSigmaEps - sigmaEps)
			sigmaEps <<- newSigmaEps

			#return( delta/sqrt(sum(prmtrAR^2)+sum(prmtrEV^2)+sum(sigmaEps^2)) )
			return( delta )
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
				  if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s", nIter, delta, toString(getPrmtr()) ) )
				  nIter <- nIter + 1
			}
			if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s, \nParameter estimated\n", nIter-1, delta, toString(getPrmtr()) ) )
	  }

	  initPrmtr <- function(tol,max.iter)
	  {
			if(verbose) message("initializing parameters")
			delta <- 1.0
			nIter <- 1

			prmtrEV <<- numeric(nEV)
			prmtrAR <<- numeric(nAR)
			sigmaEps <<- 1

			eStepNaive()
			while(delta > tol && nIter < max.iter)
			{
				  updateTrend() # needed since it is done to update the trend which is done in the eStep for our method
				  delta <- mStep()
				  if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s", nIter, delta, toString(getPrmtr()) ) )
				  nIter <- nIter + 1
			}
			if(verbose) message(sprintf("Iter:%i, delta: %f, prmtr: %s \nParameter initialized\n", nIter-1, delta, toString(getPrmtr()) ) )
	  }

	  exptdLogLik <- function()
	  {
			#wkMean2 <- wkMean
			#for(i in 1:(nAR+1)){
			#	wkMean2[(nAR+1):nObs,i] <- wkMean2[(nAR+1):nObs,i] - trend[(nAR+2-i):(nObs+1-i)]
			#}
			#tmpVec <- c(1,-prmtrAR)

			#ret <- sum( (wkMean2 %*% tmpVec)^2)
			#for(i in (nAR+1):nObs)
			#	ret <- ret + t(tmpVec)%*%wkCov%*%tmpVec

			val <- -(nObs - nSkip)*(1+log(sigmaEpsEstd^2))/2
			val
	  }



	  setGoodResiduals <- function()
	  {
			#select residuals where no censored data is involved in calculation,
			# good residuals may be of little amount if censor rate is high and AR order is large
			tmp <- rep(FALSE,nObs)
			for(i in seq(1,nObs)[-skipIndex])
				  tmp[i] <- !sum(censorIndicator[i:(i-nAR)])
			goodRes <<- getResiduals()[tmp]
			nG <- length(goodRes)
			pct <- nG/nObs
			#message(sprintf('\nNumber of good residuals:%d, %4.2f%% of %d\n',nG,pct*100,nObs))
	  }

	  bootstrapSample <- function(useGoodRes)
	  {
			if(useGoodRes)
				  eps <- sample(goodRes, nObs, replace=TRUE)
			else
				  eps <- rnorm(nObs,0, sigmaEps)

			updateTrend()

			eta <- numeric(nObs)
			eta[skipIndex] <- eps[skipIndex]
			#y[1:nAR] <<- eps[1:nAR] # this step is not necessary
			for(i in seq(1,nObs)[-skipIndex]){
				  eta[i] <- eta[(i-1):(i-nAR)] %*% prmtrAR + eps[i]
				  y[i] <<- trend[i] + eta[i]
			}
			censorIndicator <<- rep(0,nObs)
			censorIndicator[y<lowerCensorLimit] <<- -1
			censorIndicator[y>upperCensorLimit] <<- 1
			if(verbose) message(paste0("\ncensor rate: ", sum(abs(censorIndicator))/nObs))
	  }

	  bootstrapCI <- function(alpha,nBootstrapSample)
	  {
			#message(sprintf('Bootstraping CI'))
			#pb <- txtProgressBar(1,nBootstrapSample,style=3)
			yOriginal <- y #copy y as it will be overwritten
			tmpResult <- matrix(0,nBootstrapSample,getNPrmtr())
			for(i in 1:nBootstrapSample){
			      message(sprintf('Bootstraping CI %i/%i',i,nBootstrapSample))
				  #setTxtProgressBar(pb,i)
				  setInitPrmtrForBootstrap()
				  bootstrapSample(useGoodRes)
				  resetWK()
				  estimatePrmtr(tol,max.iter)
				  tmpResult[i,] <- getPrmtr()
			}
			#close(pb)
			qv <- c((1-alpha)/2,1-(1-alpha)/2)
			for(j in 1:getNPrmtr()){
				  ci[j,] <<- quantile(tmpResult[,j],qv)
			}
			covMat <<- cov(tmpResult)
			y <<- yOriginal #set back observed data
	  }

	  setFitted <- function(){
			updateCovEta()
			updateTrend()
			fittedValues <<- numeric(nObs)
			for(idx in seq(1,nObs)[-skipIndex])
			{
				  wkm <- y[(idx-1):(idx-nAR)]
				  tmpCensorIndicator <- censorIndicator[(idx-1):(idx-nAR)]
				  nCensored <- sum(abs(tmpCensorIndicator))
				  if(nCensored)
				  {# at least one is censored
						if( nCensored < nAR )
						{
							  conditionalIndex <- which(tmpCensorIndicator==0)
							  tmpY <- y[(idx-1):(idx-nAR)][conditionalIndex]
							  tmpM <- trend[(idx-1):(idx-nAR)]
							  cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta[-1,-1])
							  tmpMean <- cdist$'mean'
							  tmpVar <- cdist$'var'
						}else{
							  tmpMean <- trend[(idx-1):(idx-nAR)]
							  tmpVar <- covEta[-1,-1]
						}

						tmpLower <- rep(-Inf,length = nCensored)
						tmpUpper <- rep(Inf,length = nCensored)
						censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
						tmpLower[censored>0] <- upperCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator>0]
						tmpUpper[censored<0] <- lowerCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator<0]
						ret <- mtmvnorm(tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper,doComputeVariance=FALSE)

						wkm[tmpCensorIndicator!=0] <- ret$'tmean'
				  }
				  fittedValues[idx] <<- trend[idx] + prmtrAR%*%(wkm-trend[(idx-1):(idx-nAR)])
			}
			fittedValues[skipIndex] <<- NaN
	  }
	  getFitted <- function(){
			return(fittedValues)
	  }
      
	  # begin execution
	  resetWK()
	  #message("initializing parameters")
	  initPrmtr(tol, max.iter)
	  prmtrInit <- getPrmtr()
	  estimatePrmtr(tol,max.iter)
	  setEstdPrmtr()
	  setFitted() #must be called before setResiduals, because setResiduals use fittedValues returned by it
	  setResiduals()

	  coeff <- c(prmtrEVEstd,prmtrAREstd)
	  xnames <- colnames(x)
	  if(is.null(xnames))
			xnames <- paste0('X',1:dim(x)[2])
	  names(coeff) <- c(xnames,paste0('AR',1:nAR))

	  ret$coefficients = coeff
	  ret$prmtrEV = prmtrEVEstd
	  ret$prmtrAR = prmtrAREstd
	  ret$sigma = sigmaEpsEstd

	  ret$fitted.values = getFitted()
	  ret$residuals = getResiduals()
	  ret$censorRate = getCensorRate()
	  rnames <- c(names(coeff),"sigma")
	  ret$prmtrInit = prmtrInit
	  ret$prmtrEstd = getEstdPrmtr()
      names(ret$prmtrInit) <- rnames
	  names(ret$prmtrEstd) <- rnames
	  ret$nObs = nObs
	  ret$logLik = exptdLogLik()
	  ret$nAR = nAR
	  ret$npar = getNPrmtr()
	  ret$aic = 2*(-ret$logLik + ret$npar)
	  ret$call = match.call()

	  if(getCI){
			if(useGoodRes) setGoodResiduals()
			rnames <- c(names(coeff),"sigma")
			ci <- matrix(nrow=getNPrmtr(),ncol=2)
			covMat <- matrix(nrow=getNPrmtr(),ncol=getNPrmtr())
			bootstrapCI(alpha,nBootstrapSample)
			rownames(ci) <- rnames
			colnames(ci) <- c(sprintf("%4.2f%%",100*(1-alpha)/2), sprintf("%4.2f%%",100*(1+alpha)/2))
			rownames(covMat) <- rnames
			colnames(covMat) <- rnames
			ret$CI <- ci
			ret$vcov <- covMat
	  }else{
			ret$CI <- NULL
			ret$vcov <- NULL
	  }
	  class(ret) <- "carx"
	  ret
}


carx.formula <- function(formula, data=list(),...)
{
	  mf <- model.frame(formula=formula,data=data)
	  x <- model.matrix(attr(mf,"terms"),data=mf)
	  y <- model.response(mf)

	  est <- carx.default(y,x,...)
	  est$call <- match.call()
	  est$formula <- formula
	  est
}

carx.simulate <- function(nObs, prmtrAR, prmtrEV, sigmaEps, lowerCensorLimit, upperCensorLimit, seed=0){
	  nAR <- length(prmtrAR)
	  nEV <- length(prmtrEV)

	  set.seed(seed)
	  eps <- rnorm(nObs,0, sigmaEps)

	  externalVariable <- matrix(rnorm(nObs*nEV), nrow= nObs, ncol = nEV)
	  trend <- externalVariable%*%prmtrEV

	  eta <- numeric(nObs)
	  y <- numeric(nObs)
	  eta[1:nAR] <- eps[1:nAR]

	  y[1:nAR] <- eps[1:nAR]
	  for(i in (nAR+1):nObs){
			eta[i] <- eta[(i-1):(i-nAR)] %*% prmtrAR + eps[i]
			y[i] <- trend[i] + eta[i]
	  }
	  censorIndicator <- rep(0,nObs)
	  censorIndicator[y<lowerCensorLimit] <- -1
	  censorIndicator[y>upperCensorLimit] <- 1

	  if(options()$verbose) message(paste0("simulated series: censor rate: ", sum(abs(censorIndicator))/nObs))
	  ret <- list(y = y,
				  x = externalVariable,
				  censorIndicator=censorIndicator,
				  lowerCensorLimit=lowerCensorLimit,
				  upperCensorLimit=upperCensorLimit
				  )
	  ret
}



logLik.carx <- function(x,...)
{
	  ret <- x$logLik
	  class(ret) <- 'logLik.carx'
	  ret
}

carx.aic <- function(x,...)
{
	  val <- 2*(-x$logLik + x$npar)
	  val
}


#' print a short description of the fitted model (only a few lines)
#' @param x a fitted model object
print.carx <- function(x,...)
{
	  cat("Call:\n")
	  print(x$call)
	  cat("\nCoefficients:\n")
	  print(x$coefficients)

	  cat("\nSigma of errors of AR process:\n")
	  print(x$sigma)
	  cat("\nVariance of errors of AR process:\n")
	  print(x$sigma^2)

	  cat("\nCensor rate:\n")
	  print(x$censorRate)
	  cat("\n#obs:\n")
	  print(x$nObs)

	  cat("\nnpar:\n")
	  print(x$npar)

	  cat("\nlogLik:\n")
	  print(x$logLik)

	  cat("\nAIC:\n")
	  print(x$aic)


	  cat("\nprmtrInit:\n")
	  print(x$prmtrInit)
	  cat("\nprmtrEstd:\n")
	  print(x$prmtrEstd)
	  if(!is.null(x$CI)){
			cat("\nconfidence interval\n")
			print(x$CI)
			cat("\nvariance-covariance matrix\n")
			print(x$vcov)

	  }

}

#' summarize the fitted model on parameters, residuals and model fit
summary.carx <- function(object,...){
	  numDig <- function(x,d1,d2){
			y <- x
			y[abs(x) > 0.1^d1] <- round(x[abs(x) > 0.1^d1],d1)
			y[abs(x) < 0.1^d1] <- round(x[abs(x) < 0.1^d1],d2)
			y
	  }

	  if(is.null(object$CI)){
			tab <- cbind(Estimate = c(coef(object),object$sigma))
			#tab <- cbind(Estimate = numDig(c(coef(object),object$sigma),2,4))
	  }else{
			se <- sqrt(diag(object$vcov))
			tVal <- c(coef(object),object$sigma)/se

			est <- coef(object)
			est["sigma"] <- object$sigma
			tab <- cbind(Estimate = est,
						 StdErr =  se,
						 lowerCI = object$CI[,1],
						 upperCI = object$CI[,2]
						 )
			#tab <- cbind(Estimate = numDig(c(coef(object),object$sigma),2,4),
			#	 StdErr = numDig(se,2,4),
			#	 lowerCI = numDig(object$CI[,1],2,4),
			#	 upperCI = numDig(object$CI[,2],2,4)
			#	 )
	  }
	  res <- list(call=object$call,
				  coefficients=tab,
				  aic = carx.aic(object))
	  class(res) <- "summary.carx"
	  res
}

carx.prmtrStr <- function(object,pC=NULL){
	  prmtrStr <- ""
	  prmtrs <- summary(object)$coefficients

	  npar <- object$npar
	  if(!is.null(object$CI))
	  {
			tmp <- t(rbind(object$prmtrEstd,t(object$CI)))
			if(is.null(pC)){
				  for( i in 1:npar ){
						prmtrStr <- paste( prmtrStr, ,sprintf("%s (%s, %s)", signif(tmp[i,1],3),signif(tmp[i,2],2),signif(tmp[i,3],2)) )
						if( i < npar ){
							  prmtrStr <- paste( prmtrStr ," & ")
						}
				  }
			}
			else{
				  #get AR prmtrs
				  prmtrStr <- paste(prmtrStr, "\\pbox{4cm}{")
				  for( i in 1:pC[1]){
						vals <- prmtrs[paste0("AR",i),]
						prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
						if( i < pC[1])
							  prmtrStr <- paste(prmtrStr,"\\\\")
				  }
				  prmtrStr <- paste(prmtrStr,"} &")
				  prmtrStr <- paste(prmtrStr, "\\pbox{4cm}{")
				  if(pC[2] == 1){
						#get intercepts

						vals <- prmtrs["(Intercept)",]
						prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
				  }
				  if(pC[2]==4){
						for( i in 1:pC[2]){
							  vals <- prmtrs[paste0("as.factor(season)",i),]
							  prmtrStr <- paste( prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
							  if( i < pC[2] )
									prmtrStr <- paste(prmtrStr,"\\\\")
						}
				  }
				  prmtrStr <- paste(prmtrStr,"} &")

				  # get trend
				  prmtrStr <- paste(prmtrStr, "\\pbox{5.5cm}{")
				  vals <- prmtrs["tInMonth",]
				  prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(12*vals["Estimate"],3)),toString(signif(12*vals["lowerCI"],2)),toString(signif(12*vals["upperCI"],2))) )
				  prmtrStr <- paste(prmtrStr,"} &")


				  #get coefficients for logQ
				  prmtrStr <- paste(prmtrStr, "\\pbox{4cm}{")
				  if(pC[4] == 1)
				  {
						#no seasonal effect
						vals <- prmtrs["logQ",]
						prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
				  }
				  else{

						for( i in 1: pC[4]){
							  #get names
							  name <- sprintf("logQ:as.factor(season)%i",i)
							  if(is.na( match(name, names(prmtrs[,"Estimate"]))))
									name <- sprintf("as.factor(season)%i:logQ",i)

							  vals <- prmtrs[name, ]
							  prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
							  if( i < pC[4])
									prmtrStr <- paste(prmtrStr,"\\\\")

						}
				  }
				  prmtrStr <- paste(prmtrStr,"} &")

				  # sigma
				  prmtrStr <- paste(prmtrStr, "\\pbox{4cm}{")
				  vals <- prmtrs["sigma",]
				  prmtrStr <- paste(prmtrStr, sprintf("%s (%s, %s) ", toString(signif(vals["Estimate"],3)),toString(signif(vals["lowerCI"],2)),toString(signif(vals["upperCI"],2))) )
				  prmtrStr <- paste(prmtrStr,"} ")
			}

	  }else{
			tmp <- object$prmtrEstd
			for( i in 1:npar ){
				  if(abs(tmp[i]) > 0.01)
						prmtrStr <- paste(prmtrStr , sprintf("%.2f ",tmp[i]) )
				  else
						prmtrStr <- paste(prmtrStr , sprintf("%.4f ",tmp[i]) )
				  if( i < npar ){
						prmtrStr <- paste(prmtrStr , " & ")
				  }
			}
	  }
	  prmtrStr
}


predict.carx <- function(object,newdata=NULL,...)
{
	  if(is.null(newdata))
			y <- fitted(object)
	  else{
			if(!is.null(object$formula)){
				  x <- model.matrix(object$formula,newdata)
			} else
			{
				  x <- newdata
			}
			message("not implemented")
	  }
	  y
}



plotData <- function( carxData,timeAxis=NULL)
{
	  y <- object$y
	  lcl <- object$lowerCensorLimit
	  ucl <- object$upperCensorLimit

	  t <- seq(1,length(y))
	  ts.plot(t,y,col="black")
	  lines(t,lcl,col="red")
	  lines(t,ucl,col="red")

}





plot.carx <- function(object,transformFun=NULL,xAxisVar=NULL,xlab="",ylab="",saveFig="", outliers=NULL,...)
{
	  setEPS()
	  postscript(saveFig)

	  yh <- predict(object)
	  y <- object$y
	  lcl <- object$lowerCensorLimit
	  ucl <- object$upperCensorLimit
	  ylim <- range(c(y,yh,censorLimit),na.rm=TRUE)


	  #y[object$censorIndicator] <- NA
	  y[object$censorIndicator] <- object$censorLimit[object$censorIndicator]

	  if(is.null(xAxisVar))
	  {
			xAxisVar <- 1:length(yh)
	  }

	  xrange <- range(xAxisVar)
	  if(is.null(transformFun))
	  {
			yrange <- ylim
			plot(as.zoo(as.ts(zoo(y, xAxisVar))), lty=1,xlab=xlab,ylab=ylab,ylim=ylim)
			#plot(xAxisVar,y,'l',lty=1,xlab=xlab,ylab=ylab,ylim=ylim)
			lines(as.zoo(as.ts(zoo(yh, xAxisVar))), lty=2,col='blue') #xlab=xlab,ylab=ylab,ylim=ylim)
			#lines(xAxisVar,yh,'l',lty=2,col='blue')
			lines(as.zoo(as.ts(zoo(object$censorLimit, xAxisVar))),lty=1,col="red")
			#lines(xAxisVar[object$censorIndicator],object$censorLimit[object$censorIndicator],col="red",pch='*')
	  }else{
			yrange <- transformFun(ylim)
			plot(as.zoo(as.ts(zoo(transformFun(y), xAxisVar))), lty=1,xlab=xlab,ylab=ylab,ylim=transformFun(ylim))
			#plot(as.zoo(as.ts(zoo(transformFun(y),xAxisVar))),'l',lty=1,xlab=xlab,ylab=ylab)
			#plot(xAxisVar,transformFun(y),'l',lty=1,xlab=xlab,ylab=ylab,ylim=transformFun(ylim))
			lines(as.zoo(as.ts(zoo(transformFun(yh), xAxisVar))), lty=2,col='blue') #xlab=xlab,ylab=ylab,ylim=ylim)
			#lines(xAxisVar,transformFun(yh),lty=2,col='blue')
			lines(as.zoo(as.ts(zoo(transformFun(object$censorLimit), xAxisVar))),lty=1,col="red")
			#lines(xAxisVar,transformFun(object$censorLimit),lty=1,col="red")
			#lines(xAxisVar[object$censorIndicator],transform(object$censorLimit[object$censorIndicator]),col="red",pch='*')
	  }
	  legend(xrange[2]*(0.7),yrange[2]*0.7,legend=c(ylab,'Fitted value','Censor limit'),lty=c(1,2,1),col=c('black','blue','red'))
	  if(!is.null(outliers)) abline(v=xAxisVar[outliers],col="red",lty=2)
	  dev.off()
	  #if(saveFig != "")
	  #{
	  #dev.copy2eps(file=saveFig)
	  #dev.off()
	  #}

}

plot.residuals <- function(object,x=NULL,saveFig="",xlab="",ylab="",classify.by=NULL,type="l",lty=1)
{
	  #if(saveFig != "") jpeg(saveFig)
	  y <- residuals(object)/object$sigma
	  if(is.null(x))
	  {
			x <- 1:length(y)
			xlab <- "index"
	  }else{
			if(typeof(x) == "character")
			{
				  if(x == 'fitted'){
						x <- predict(object)
						xlab <- "Fitted value"
				  }
			}
	  }
	  plot(x,y,xlab=xlab,ylab=ylab,type=type,lty=lty)
	  abline(h=0,col="black",lty=3)
	  #abline(h=1.96*object$sigma,col="red",lty=3)
	  #abline(h=-1.96*object$sigma,col="red",lty=3)

	  if(saveFig != "")
	  {
			dev.copy2eps(file=saveFig)
			dev.off()
	  }
}

plotResAcf <- function(object,saveFig="")
{
	  setEPS()
	  postscript(saveFig)
	  acf(residuals(object),na.action=na.pass,main="")
	  dev.off()
	  #if(saveFig != "")
	  #{
	  #dev.copy2eps(file=saveFig)
	  #dev.off()
	  #}
}


print.summary.carx <- function(x,...)
{
	  cat("Call:\n")
	  print(x$call)
	  cat("\n")

	  cat("\nCoefficients:\n")
	  print(x$coefficients)

	  cat("\nAIC:\n")
	  print(x$aic)
}



summarizeResult <- function(sampleSize,rslt, prefix){
	  nRep <- dim(rslt)[1]
	  nPrmtr <- dim(rslt)[2]
	  #par(mfrow=c(nPrmtr,3))
	  for( i in 1:nPrmtr)
	  {
			jpeg(paste0(prefix,'_nSample',sampleSize,'_nRep',nRep,'_plot',i,'_hist.jpg'))
			hist(rslt[,i])
			dev.off()
			jpeg(paste0(prefix,'_nSample',sampleSize,'_nRep',nRep,'_plot',i,'_density.jpg'))
			plot( density(rslt[,i]) )
			dev.off()
			jpeg(paste0(prefix,'_nSample',sampleSize,'_nRep',nRep,'_plot',i,'_qqnorm.jpg'))
			qqnorm(rslt[,i])
			dev.off()
	  }
}


carx.singleSimulation <- function(iRep=1)
{
	  trueEV <- c(0.2,0.4)
	  trueAR <- c(0.1,0.3,-0.2)
	  trueSigma <- sqrt(0.5)
	  sampleSize <- 100
	  lcl <- -1.0
	  ucl <- 1.5
	  #iRep <- 1
	  fullEstimation <- T

	  args <- commandArgs(TRUE)
	  if(length(args) > 0)
	  {
			for(i in 1:length(args))
			{
				  eval(parse(text = args[i]))
			}
	  }


	  nAR <- length(trueAR)
	  truePrmtr <- c(trueAR,trueEV,trueSigma)
	  nPrmtr <- length(truePrmtr)
	  lowercl <- rep(lcl, sampleSize)
	  uppercl <- rep(ucl, sampleSize)


	  dat <- carx.simulate(sampleSize, trueAR, trueEV, trueSigma, lowercl, uppercl, seed=37513*iRep)
	  #print(dat)

	  rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$lowerCensorLimit,dat$upperCensorLimit, nAR, getCI=fullEstimation,skipIndex=seq(1,nAR))

	  ret <- c(rslt$prmtrEstd,rslt$prmtrInit)
	  if(fullEstimation)
	  {
			ci <- rslt$CI
			ret <- c(ret, ci[,1])
			ret <- c(ret, ci[,2])
			coverage <- (truePrmtr >= ci[,1])*(truePrmtr <= ci[,2])
			ret <- c(ret,coverage)
	  }

	  d <- paste0('./sim_n',sampleSize,'_l_',lcl,'_u_',ucl)
	  dir.create(d,showWarnings=FALSE)
	  cdir <- getwd()
	  setwd(d)
	  f <- file(paste0(toString(iRep),".txt"))
	  writeLines(paste(ret,collapse=' '),f)
	  close(f)
	  setwd(cdir)
	  rslt
}

carx.simulation <-function(trueEV=c(0.2,0.4), 
						   trueAR=c(0.1,0.3,-0.2), 
						   trueSigma=sqrt(0.5),
						   lcl=-1,
						   ucl=1,
						   sampleSize=100,
						   nRep=100,
						   fullEstimation=T)
{
	  args <- commandArgs(TRUE)
	  if(length(args) > 0)
	  {
			for(i in 1:length(args))
			{
				  eval(parse(text = args[i]))
			}
	  }
      lowercl <- rep(lcl, sampleSize)
	  uppercl <- rep(ucl, sampleSize)


	  nAR <- length(trueAR)
	  truePrmtr <- c(trueAR,trueEV,trueSigma)
	  nPrmtr <- length(truePrmtr)

	  prmtrEstd <- matrix(0,nRep,nPrmtr)
	  prmtrInit <- matrix(0,nRep,nPrmtr)
	  estCI <- matrix(0,nRep,2*nPrmtr)
	  coverage <- numeric(nPrmtr)

	  simEst <-function(iRep)
	  {
			message(iRep)
			dat <- carx.simulate(sampleSize, trueAR, trueEV, trueSigma, lowercl, uppercl, seed=37513*iRep)
			rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$lowerCensorLimit,dat$upperCensorLimit, nAR, getCI=fullEstimation ,skipIndex=seq(1,nAR))
            ret <- c(iRep,rslt$censorRate,rslt$prmtrInit,rslt$prmtrEstd)
			if(fullEstimation)
			{
				  ci <- rslt$CI
			  ret <- c(ret, ci[,1])
				  ret <- c(ret, ci[,2])
				  coverage <- (truePrmtr >= ci[,1])*(truePrmtr <= ci[,2])
				  ret <- c(ret,coverage)
			}
			ret
	  }

	  iter <- 1:nRep
	  rslt <- mclapply(iter,simEst,mc.cores=detectCores())
	  rslt
}



# tests
#rslt <- carx.singleSimulation()
#rslt
#summary(rslt)

rslt <- carx.simulation()
print(rslt)
#for(i in 1:10) 
	  #carx.singleSimulation(i)
