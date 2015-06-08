#' compute the residuals of a fitted \code{carx} object
#'
#' \code{residuals.carx} computes the simulated residuals from a fitted \code{carx} object.
#' @param object a fitted \code{carx} object.
#' @param type a string indicates which type of residual is to be returned.
#' "raw" returns the simulated residuals;
#' "pearson" returns the raw residuals divided by estimated standard error of the residuals.
#' @param ... not used.
#' @return the simulated residuals.
#' @export
residuals.carx <- function(object,type="pearson",...)
{
	nObs <- object$nObs
	nAR <- object$nAR
	y <- object$y

	trend <- object$x%*%object$prmtrX
	eta <- object$y - trend
	for(idx in 1:nAR)
	{
		if(object$censorIndicator[idx] > 0)
			y[idx] = object$upperCensorLimit[idx]
		else if(object$censorIndicator[idx] < 0 )
				y[idx] = object$lowerCensorLimit[idx]
	}

	for(idx in (nAR+1):nObs)
	{
		if(object$censorIndicator[idx] == 0)
		{
			#message(sprintf("idx %i, not censored",idx))
			next
		}

		if(all(object$censorIndicator[(idx-1):(idx-nAR)]==0))
    {
		  #message(sprintf("idx %i, fast ",idx))
			tmpMean <- trend[idx] + eta[(idx-1):(idx-nAR)]%*%object$prmtrAR
			if(object$censorIndicator[idx] > 0)
				y[idx] <- tmvtnorm::rtmvnorm(1, mean=c(tmpMean), sigma=c(object$sigma),lower=object$upperCensorLimit[idx],upper=Inf,algorithm="gibbs")
			else
				y[idx] <- tmvtnorm::rtmvnorm(1, mean=c(tmpMean), sigma=c(object$sigma),lower=-Inf, upper = object$lowerCensorLimit[idx],algorithm="gibbs")
			next
		}

		#message(sprintf("idx %i, slow",idx))
		iStart <- 1
		for(i in (idx-nAR):1)
		{
			if(all(object$censorIndicator[i:(i+nAR-1)]==0))
			{
				iStart <- i
				break
			}
		}

		nStart <- idx - iStart + 1
		#message(sprintf("idx: %i, iStart: %i, nStart: %i",idx,iStart,nStart))
		tmpCensorIndicator <- object$censorIndicator[idx:iStart] #reverse order
		nCensored <- sum(tmpCensorIndicator!=0)
		covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart)
		if( nCensored < nStart )
		{
			conditionalIndex <- which(tmpCensorIndicator==0)
			tmpY <- object$y[idx:iStart][conditionalIndex]
			cdist <- conditionalDistMvnorm(tmpY, conditionalIndex, trend[idx:iStart], covEta)
			tmpMean <- cdist$'mean'
			tmpVar <- cdist$'var'
		}else
		{
			tmpMean <- trend[idx:iStart]
			tmpVar <- covEta
		}
		tmpLower <- rep(-Inf,length = nCensored) #( y[idx], censored obs)
		tmpUpper <- rep(Inf,length = nCensored)
		censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
		tmpLower[censored>0] <- object$upperCensorLimit[idx:iStart][tmpCensorIndicator>0]
		tmpUpper[censored<0] <- object$lowerCensorLimit[idx:iStart][tmpCensorIndicator<0]
		ysim <- tmvtnorm::rtmvnorm(1, mean=tmpMean, sigma=tmpVar, lower=tmpLower, upper=tmpUpper,algorithm="gibbs")
		y[idx] <- ysim[1]
  }
	#print(y)
	#print(object)
	#m1 <- arimax(y,order=c(object$nAR,0,0),xreg = data.frame(object$x),include.mean=FALSE,transform.pars=FALSE,init=c(object$prmtrAR,object$prmtrX),method="ML")
	#print(m1)
	m2 <- carx(y, object$x, rep(0,nObs), NULL, NULL, object$nAR, getCI=FALSE)
	#print(m2)
	rsdl <- numeric(nObs)
	eta <- y - object$x%*%m2$prmtrX
	for(idx in (nAR+1):nObs)
	{
		rsdl[idx] <- eta[idx] - eta[(idx-1):(idx-nAR)]%*%m2$prmtrAR
	}
	if(type=="raw")
		return(rsdl)
	if(type=="pearson")
	        return(rsdl/m2$sigma)
}

