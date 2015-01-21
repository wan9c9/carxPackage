library(stats)
library(tmvtnorm)
library(debug)
library(ggplot2)
library(R2HTML)
library(zoo)
library(matrixStats)


#'calculate the conditional mean & variance of a multivariate normal distribution 
#' with (mean, varMat), conditional on y
#' @y a vector to be conditioned on
#' @conditionalIndex the index to be conditioned on
#' @meanVec the mean vector the joint multivariate normal distribution
#1 @varMat the variance-covariance matrix of the joint multivariate normal distribution
conditionalDistMvnorm <- function(y,conditionalIndex,meanVec,varMat)
{
	sigma11 <- varMat[-conditionalIndex,-conditionalIndex]
	sigma22 <- varMat[ conditionalIndex, conditionalIndex]
	sigma12 <- varMat[-conditionalIndex, conditionalIndex]
	sigma21 <- varMat[ conditionalIndex,-conditionalIndex]
	invSigma22 <- solve(sigma22)
	mNew <- meanVec[-conditionalIndex] + as.vector(sigma12 %*% invSigma22 %*%(y - mnVec[conditionalIndex]))
	vNew <- sigma11 - sigma12 %*% invSigma22 %*% sigma21

	return (list('mean'=mNew,'var' = vNew))
}

#' compute the covariance matrix of {\eta_t,\cdots,\eta_{t-p}} of the AR process
#' pAR the order of the AR model
#' sigmaEps the standard deviation of the residuals of the AR model
computeCovAR <- function(pAR, sigmaEps)
{
	n <- length(pAR)+1
	val <- ARMAacf(ar=pAR, lag.max=n)
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
	v <- sigmaEps^2/(1-pAR%*%val[2:n])
	mat <- v[1,1]*mat
	return(mat)
}


carx <- function(x,...) UseMethod("carx")

carx.default <-
	function(y,x,censorIndicator,censorLimit,nAR,tol=1e-4,max.iter=500,getCI=FALSE,alpha=0.95,nBootstrapSample=1000, useGoodRes=FALSE,skipIndex=NaN){

	nObs <- dim(x)[1]
	nEV <- dim(x)[2]
	if(length(y) != nObs){
		message(" dimension doesn't match!")
		return(NULL)
	}

	if( !is.vector(skipIndex))
	{
		skipIndex <- seq(1,nAR)
		message("skip index is constructed as ")
		message(skipIndex)
	}
	#if(!is.nan(skipIndex))
	#{
	#	for(i in restartIndex)
	#		skipIndex <- c(skipIndex,seq(i,i+nAR-1))
	#}
	nSkip <- length(skipIndex)


	censorIndicator <- as.logical(censorIndicator)

	ret = list(y = y,
		   x = x,
		   censorIndicator = censorIndicator,
		   censorLimit = censorLimit,
		   skipIndex = skipIndex
		   )

	#parameters
	prmtrAR <- numeric(nAR)
	sigmaEps <- numeric(1)
	prmtrEV <- numeric(nEV)

	prmtrAREstd <- numeric(nAR)
	sigmaEpsEstd <- numeric(1)
	prmtrEVEstd <- numeric(nEV)

	externalVariable <- x


	# working data
	resetWK <- function(){
		trend <<- numeric(nObs)
		covEta <<- matrix(nrow=nAR+1, ncol = nAR+1)
		wkMean <<- matrix(nrow=nObs, ncol = nAR+1)
		wkCov <<- array(0, dim=c(nObs,nAR+1,nAR+1))
		res <<- numeric(nObs)
	}
	#-------------------------------

	getNPrmtr <- function(){
		return (nAR + nEV + 1)
	}

	getPrmtr <- function(){
		ret <- c(prmtrAR, prmtrEV, sigmaEps)
		return(ret)
	}

	getCensorRate <- function(){
		return(sum(censorIndicator)*1.0/length(censorIndicator))
	}


	setInitPrmtrForBootstrap <- function()
	{
		prmtrAR   <<-		prmtrAREstd
		prmtrEV   <<-     prmtrEVEstd
		sigmaEps  <<-	sigmaEpsEstd
	}

	setEstdPrmtr <- function()
	{
		prmtrAREstd   <<- prmtrAR 
		prmtrEVEstd   <<- prmtrEV
		sigmaEpsEstd <<- sigmaEps
	}

	getEstdPrmtr <- function(){
		ret <- c(prmtrAREstd, prmtrEVEstd, sigmaEpsEstd)
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
		#print("E-step begin.")
		# expection step
		# need to 1) calculate the joint distriubtion of AR process
		#         2) update the trend function/observations
		#         3) calculat the expection & covariances of the censored observations
		updateCovEta()
		updateTrend()
		for(idx in seq(1,nObs)[-skipIndex])
		{
			wkMean[idx,] <<- y[idx:(idx-nAR)] #be careful that the wk*'s are persistent so keep them                     
			#evaluated every time for safety
			wkCov[idx,,] <<- 0                # so here
			tmpCensorIndicator <- censorIndicator[idx:(idx-nAR)]
			#print(c(idx, sum(tmpCensorIndicator), nAR+1))
			nCensored <- sum(tmpCensorIndicator)
			if(nCensored)
			{# at least one is censored
				if( nCensored < (nAR+1) )
				{
					conditionalIndex <- which(!tmpCensorIndicator)
					tmpY <- y[idx:(idx-nAR)][conditionalIndex]
					tmpM <- trend[idx:(idx-nAR)]
					cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta)
					tmpMean <- cdist$'mean'
					tmpVar <- cdist$'var'
				}else{
					tmpMean <- trend[idx:(idx-nAR)]
					tmpVar <- covEta
				}
				tmpCensorLimit <- censorLimit[idx:(idx-nAR)][tmpCensorIndicator]
				ret <- mtmvnorm(tmpMean,tmpVar,upper=tmpCensorLimit)

				wkMean[idx,tmpCensorIndicator] <<- ret$'tmean'      
				wkCov[idx,tmpCensorIndicator,tmpCensorIndicator] <<- ret$'tvar' 
			} 
			#message("E-step done.") 
		}
	}

	eStepNaive <- function(){
		#message("entering eStepNaive")
		#updateTrend()
		for(idx in seq(1,nObs)[-skipIndex]){
			tmpCensorIndicator <- censorIndicator[idx:(idx-nAR)]
			#print(c(idx, sum(tmpCensorIndicator), nAR+1))
			wkMean[idx,] <<- y[idx:(idx-nAR)]
			wkMean[idx,tmpCensorIndicator] <<- censorLimit[idx:(idx-nAR)][tmpCensorIndicator]
			wkCov[idx,,] <<- 0                
		}
		#print("E-step Naive done.")
	}

	setResiduals <- function(){
		res <<- numeric(nObs)
		res[skipIndex] <<- NaN
		res[-skipIndex] <<- y[-skipIndex] - fittedValues[-skipIndex]
		res[censorIndicator] <<- censorLimit[censorIndicator] - fittedValues[censorIndicator]
	}

	setResiduals_old <- function(){
		wkMean2 <- wkMean
		index <- seq(1,nObs)[-skipIndex]
		for(i in 1:(nAR+1)){
			wkMean2[index,i] <- wkMean2[index,i] - trend[(index+1-i)]
		}
		ret <- as.vector(wkMean2 %*% c(1,-prmtrAREstd))
		#ret[is.na(ret)] <- 0  # the residuals on skipIndex are not defined and assigned to be zero
		ret[skipIndex] <-  0 # the residuals on skipIndex are not defined and assigned to be zero
		#message(sprintf("number of nan's in residuals: %d\n",sum(is.na(ret))))
		res <<- ret
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
		#tmpY <- tmpY[(nAR+1):nObs]
		tmpY <- tmpY[-skipIndex]

		#tmpX <- matrix(0,nObs-nAR,nEV)
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

	updateSigmaEps0 <- function(){
		vec <- wkMean[(nAR+1):nObs,1] - trend[(nAR+1):nObs]
		for(i in 1:nAR){
			vec <- vec - prmtrAR[i]*(wkMean[(nAR+1):nObs,i+1] - trend[(nAR+1-i):(nObs-i)])
		}
		ret <- sum(vec^2)
		tmpVec <- c(1,-prmtrAR)
		for(i in (nAR+1):nObs){
			ret <- ret + as.numeric( tmpVec%*%wkCov[i,,]%*%tmpVec )
		}
		return( sqrt(ret/(nObs-nAR)) )
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


	estimatePrmtr <- function(tol,max.iter){
		#message("estimating parameters")
		delta <- 1.0  
		nIter <- 1
		while(delta > tol && nIter < max.iter){
			#print(paste0("nIter: ",nIter))
			eStep()
			delta <- mStep()
			#message(sprintf("Iter:%i , delta: %f, prmtr: %s \n", nIter, delta, toString(getPrmtr()) ) )
			nIter <- nIter + 1
		}
		message(sprintf("Iter:%i , delta: %f, prmtr: %s, \nParameter estimated\n", nIter-1, delta, toString(getPrmtr()) ) )
	}

	initPrmtr <- function(tol,max.iter){
		message("initializing parameters")
		delta <- 1.0  
		nIter <- 1

		prmtrAR <<- numeric(nAR)
		prmtrEV <<- numeric(nEV)
		sigmaEps <<- 1

		eStepNaive()
		while(delta > tol && nIter < max.iter){
			updateTrend() # needed since it is done to update the trend which is done in the eStep for our method
			#print(paste0("nIter: ",nIter))
			delta <- mStep()
			#message(sprintf("Iter:%i , delta: %f, prmtr: %s \n", nIter, delta, toString(getPrmtr()) ) )
			nIter <- nIter +1
		}
		message(sprintf("Iter:%i , delta: %f, prmtr: %s \nParameter initialized\n", nIter-1, delta, toString(getPrmtr()) ) )
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



	setGoodResiduals <- function(){
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
			eps <- rnorm(nObs,0, sigmaEps )

		updateTrend()

		eta <- numeric(nObs)
		eta[skipIndex] <- eps[skipIndex]
		#y[1:nAR] <<- eps[1:nAR] # this step is not necessary
		for(i in seq(1,nObs)[-skipIndex]){
			eta[i] <- eta[(i-1):(i-nAR)] %*% prmtrAR + eps[i]
			y[i] <<- trend[i] + eta[i]
		}
		censorIndicator <<- (y < censorLimit)
		#message(paste0("censor rate: ", sum(censorIndicator)/nObs))
	}

	bootstrapCI <- function(alpha,nBootstrapSample)
	{
		yOriginal <- y #copy y as it will be overwritten
		tmpResult <- matrix(0,nBootstrapSample,getNPrmtr())
		for( i in 1:nBootstrapSample ){
			message(sprintf('bootstraping CI ... %i/%i',i,nBootstrapSample))
			setInitPrmtrForBootstrap()
			bootstrapSample(useGoodRes)
			estimatePrmtr(tol,max.iter)
			tmpResult[i,] <- getPrmtr()
		}
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
			nCensored <- sum(tmpCensorIndicator)
			if(nCensored)
			{
				# at least one is censored 
				if( nCensored < nAR )
				{ 
					conditionalIndex <- which(!tmpCensorIndicator)
					tmpY <- y[(idx-1):(idx-nAR)][conditionalIndex] 
					tmpM <- trend[(idx-1):(idx-nAR)] 
					cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta[-1,-1]) 
					tmpMean <- cdist$'mean' 
					tmpVar <- cdist$'var' 
				}else{ 
					tmpMean <- trend[(idx-1):(idx-nAR)] 
					tmpVar <- covEta[-1,-1]
				}
				tmpCensorLimit <- censorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator]
				ret <- mtmvnorm(tmpMean,tmpVar,upper=tmpCensorLimit) 
				wkm[tmpCensorIndicator] <- ret$'tmean'
			}
			fittedValues[idx] <<- trend[idx] + prmtrAR%*%(wkm-trend[(idx-1):(idx-nAR)])
		}
		fittedValues[skipIndex] <<- NaN
	}
	getFitted <- function(){
		return(fittedValues)
	}

	#message("initializing parameters")
	resetWK()
	initPrmtr(tol, max.iter)
	prmtr0 <- getPrmtr()
	estimatePrmtr(tol,max.iter)
	setEstdPrmtr()
	setFitted() #must be called before setResiduals, because setResiduals use fittedValues returned by it
	setResiduals()

	coeff <- c(prmtrAREstd,prmtrEVEstd)
	names(coeff) <- c(paste0('AR',1:nAR), colnames(x))

	ret$coefficients = coeff
	ret$sigma = sigmaEpsEstd
	ret$prmtrAR = prmtrAREstd
	ret$prmtrEV = prmtrEVEstd
	ret$fitted.values = getFitted()
	ret$residuals = getResiduals()
	ret$censorRate = getCensorRate()
	ret$prmtr0 = prmtr0
	ret$prmtr1 = getEstdPrmtr()
	ret$nObs = nObs
	ret$logLik = exptdLogLik()
	ret$nAR = nAR
	ret$npar = getNPrmtr()
	ret$aic = 2*(-ret$logLik + ret$npar)
	ret$call = match.call()

	if(getCI){
		if(useGoodRes) setGoodResiduals()
		ci <- matrix(nrow=getNPrmtr(),ncol=2)
		covMat <- matrix(nrow=getNPrmtr(),ncol=getNPrmtr())
		bootstrapCI(alpha,nBootstrapSample)
		colnames(covMat) <- rownames(covMat) <- c(names(coeff),"sigma")
		ret$CI <- ci
		ret$vcov <- covMat
	}else{
		ret$CI <- NULL
		ret$vcov <- NULL
	}

	class(ret) <- "carx"
	ret
}



outlierDetection <- function(model){
	message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/model$nObs
	eps <- rnorm(nSample,0,model$sigma)
	trend <- model$x%*%model$prmtrEV
	covEta <- computeCovAR(model$prmtrAR, model$sigma)
	nObs <- model$nObs
	nAR <- model$nAR
	prmtrAR <- model$prmtrAR
	skipIndex <- model$skipIndex
	y <- model$y
	censorIndicator <- model$censorIndicator
	censorLimit <- model$censorLimit
	y[censorIndicator] <- censorLimit[censorIndicator]


	pValues <- numeric(nObs)
	pValues[skipIndex] <- 1

	for(idx in seq(1,nObs)[-skipIndex])
	{
		#message(sprintf("checking %i",idx))
		wkm <- y[(idx-1):(idx-nAR)]
		tmpCensorIndicator <- censorIndicator[(idx-1):(idx-nAR)]
		nCensored <- sum(tmpCensorIndicator)
		if(nCensored)   #at least one is censored
		{
			if( nCensored < nAR ) #not all are censored
			{ 
				conditionalIndex <- which(!tmpCensorIndicator)
				tmpY <- y[(idx-1):(idx-nAR)][conditionalIndex] 
				tmpM <- trend[(idx-1):(idx-nAR)]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta[-1,-1])
				tmpMean <- cdist$'mean' 
				tmpVar <- cdist$'var' 
			}else{ 
				tmpMean <- trend[(idx-1):(idx-nAR)]
				tmpVar <- covEta[-1,-1]
			}
			tmpCensorLimit <- censorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator]
			#print(tmpMean)
			#print(tmpVar)
			smpl <- rtmvnorm(nSample,tmpMean,tmpVar,upper=tmpCensorLimit,algorithm="gibbs")
			smpl <- as.matrix(smpl)
			#print(smpl)
			ySmpl <- numeric(nSample)
			for(i in 1:nSample){
				wkm[tmpCensorIndicator] <- smpl[i,]
				ySmpl[i] <- trend[idx] + (wkm - trend[(idx-1):(idx-nAR)])%*%prmtrAR + eps[i]
			}
			pU <- sum(ySmpl > y[idx])/nSample
			pL <- sum(ySmpl < y[idx])/nSample
		}
		else{
			r <- y[idx]-trend[idx] - (wkm-trend[(idx-1):(idx-nAR)])%*%prmtrAR
			r <- r/model$sigma
			pU <- pnorm(r,lower.tail=FALSE)
			pL <- pnorm(r,lower.tail=TRUE)
		}
		pValues[idx] <- min(pU,pL)
	}
	minP <- min(pValues)
	if( minP <= threshold ){
		i <- which(pValues == minP)
		i
	} else
		-1
}

setCensorLimit <- function(cl,nObs) { 
	if(length(cl) == 1){
		censorLimit <- rep(cl,nObs)
	}
	else {
		if(length(cl) == nObs){
			censorLimit <- cl
		}
		else{
			print("error censor limit")
		}
	}
	censorLimit
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

carx.simulate <- function(nObs, prmtrAR, prmtrEV, sigmaEps, censorLimit, seed=0){
	nAR <- length(prmtrAR)
	nEV <- length(prmtrEV)

	cl <- setCensorLimit(censorLimit,nObs)

	set.seed(seed)
	eps <- rnorm(nObs,0, sigmaEps )

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
	censorIndicator <- (y < cl)

	message(paste0("censor rate: ", sum(censorIndicator)/nObs))
	ret <- list( y = y, x = externalVariable, censorIndicator=censorIndicator,censorLimit=cl)
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

	cat("\n#npar:\n")
	print(x$npar)

	cat("\n#logLik:\n")
	print(x$logLik)

	cat("\n#AIC:\n")
	print(x$aic)


	cat("\n#prmtr0:\n")
	print(x$prmtr0)
	cat("\n#prmtr1:\n")
	print(x$prmtr1)
	if(!is.null(x$CI)){
		cat("\n# confidence interval\n")
		print(x$CI)
		cat("\n# variance-covariance matrix\n")
		print(x$vcov)

	}

}

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
		tmp <- t(rbind(object$prmtr1,t(object$CI)))
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
		tmp <- object$prmtr1
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

plot.carx <- function(object,transformFun=NULL,xAxisVar=NULL,xlab="",ylab="",saveFig="", outliers=NULL,...)
{
	setEPS()
	postscript(saveFig)

	yh <- predict(object)
	censorLimit <- object$censorLimit
	y <- object$y
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


carx.singleSimulation <- function()
{
    trueAR <- c(0.1,0.3,-0.2)
    trueEV <- c(0.2,0.4)
    trueSigma <- sqrt(0.5)
    sampleSize <- 100
    climit <- -100
    iRep <- 1
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

	dat <- carx.simulate(sampleSize, trueAR, trueEV, trueSigma, climit, seed=37513*iRep)

	censorRate <- sum(dat$y < climit)
	rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$censorLimit,nAR, getCI=fullEstimation,skipIndex=seq(1,nAR))
  
  ret <- c(rslt$prmtr1,rslt$prmtr0)
	if(fullEstimation)
	{
		ci <- rslt$CI
		ret <- c(ret, ci[,1])
		ret <- c(ret, ci[,2])
		coverage <- (truePrmtr >= ci[,1])*(truePrmtr <= ci[,2]) 
		ret <- c(ret,coverage)
	}
 
	d <- paste0('./sim_n',sampleSize,'_c',climit)
	dir.create(d,showWarnings=FALSE)
	cdir <- getwd()
	setwd(d)
	f <- file(paste0(toString(iRep),".txt"))
	writeLines(paste(ret,collapse=' '),f)
	close(f)
	setwd(cdir)
}
	



carx.simulation <-function(sampleSize=100,climit=0,nRep=1,fullEstimation=FALSE){
	args <- commandArgs(TRUE)
  if(length(args) > 0)
	{
		for(i in 1:length(args))
		{
			eval(parse(text = args[i]))
		}
	}
	dir.create(paste0('./sim_n',sampleSize,'_nRep', nRep,'_c',climit))
	cdir <- getwd()
	setwd(paste0('./sim_n',sampleSize,'_nRep', nRep,'_c',climit))
	

	trueAR <- c(0.1,0.3,-0.2)
	trueEV <- c(0.2,0.4)
	trueSigma <- sqrt(0.5)

	nAR <- length(trueAR)
	truePrmtr <- c(trueAR,trueEV,trueSigma)
	nPrmtr <- length(truePrmtr)

	estPrmtr <- matrix(0,nRep,nPrmtr)
	estPrmtrNaive <- matrix(0,nRep,nPrmtr)
	estCI <- matrix(0,nRep,2*nPrmtr)
	coverage <- numeric(nPrmtr)
	#t0 <- proc.time()
	censorRate <- 0
	for(i in 1:nRep)
	{
		message(sprintf("nRep: %i",i))
		dat <- carx.simulate(sampleSize, trueAR, trueEV, trueSigma, climit, seed=37513*i)
		censorRate <- censorRate + sum(dat$y < climit)
		#break
		#write.csv(dat$y,"y2.dat") #write.csv(dat$x,"y2.dat") #write.csv(dat$y,"y2.dat")
		rslt <- carx.default(dat$y,dat$x,dat$censorIndicator,dat$censorLimit,nAR, getCI=fullEstimation ,skipIndex=seq(1,nAR))
		#print(rslt)
		estPrmtr[i,] <- c(rslt$coefficients,rslt$sigma)
		estPrmtrNaive[i,] <- rslt$prmtr0

		if(fullEstimation)
		{
			ci <- rslt$CI
			estCI[i,] <- as.vector(ci)
			coverage <- coverage + (truePrmtr >= ci[,1])*(truePrmtr <= ci[,2]) 
		}
	}
	censorRate <- censorRate/(nRep*sampleSize)
	message(sprintf("climit %f, censor rate %f \n", climit, censorRate))

	write.csv(estPrmtr, paste0('n_',sampleSize,'_nRep', nRep,'_c',climit,'_estdPrmtr.csv'), row.names=FALSE,col.names=FALSE)
	write.csv(estPrmtrNaive, paste0('n_',sampleSize,'_nRep', nRep,'_c',climit,'_estdPrmtrNaive.csv'), row.names=FALSE,col.names=FALSE)
	if(fullEstimation)
	{
		coverage <- coverage*1.0/nRep
		write.csv(coverage, paste0('n_',sampleSize,'_nRep', nRep,'_c',climit,'_coverageRate.csv'), row.names=FALSE,col.names=FALSE)
		write.csv(estCI, paste0('n_',sampleSize,'_nRep', nRep,'_c',climit,'_estdCI.csv'), row.names=FALSE,col.names=FALSE)
	}

  
	m0 <- colMeans(estPrmtrNaive)
	s0 <- colSds(estPrmtrNaive)
	m <- colMeans(estPrmtr)
	s <- colSds(estPrmtr)

	cat(censorRate, file="resultStr.txt",append=TRUE,sep="\n")
  
	rsltStr0 <- ""
	rsltStr <- ""

	for( i in seq(1,nPrmtr))
	{
		rsltStr0 <- paste0(rsltStr0, sprintf( "%6.4f (%6.4f)",m0[i],s0[i]))
		rsltStr <- paste0(rsltStr, sprintf( "%6.4f (%6.4f) [%6.4f%%]",m[i],s[i],coverage[i]))
		if(i < nPrmtr)
		{ 
			rsltStr0 <- paste0(rsltStr0, " & " )
			rsltStr <- paste0(rsltStr, " & " )
		}
	}
	print(rsltStr)
	cat(rsltStr,file="resultStr.txt",append=TRUE,sep="\\\\\n")
	cat(rsltStr0,file="resultStr.txt",append=TRUE,sep="\\\\\n")

	summarizeResult(sampleSize,estPrmtr,paste0('n',sampleSize,'_nRep', nRep,'_c',climit))
	summarizeResult(sampleSize,estPrmtrNaive,paste0('n',sampleSize,'_nRep', nRep,'_c',climit,'_naive'))
	setwd(cdir)
	return(rsltStr)
}

#carx.singleSimulation()
