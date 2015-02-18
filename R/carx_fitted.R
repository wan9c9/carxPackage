#' get the fitted values of a \code{carx} object
fitted.carx <- function(object,...)
{
	message("Calling fitted.carx")
	nObs  <- object$nObs
	nAR  <- object$nAR

	ret <- rep(NA,nObs)
	trend <- object$x%*%object$prmtrX
	eta <- object$y - trend

	for(idx in (nAR+1):nObs)
	{
		#message(sprintf("calculating %i",idx))
		if(all(object$censorIndicator[(idx-1):(idx-nAR)]==0))
		{
			ret[idx] <- trend[idx] + eta[(idx-1):(idx-nAR)]%*%object$prmtrAR
		}
		else
		{
			#message(sprintf("Index %i is censored",idx))
			#find the beginning index of nAR consecutive observations
			iStart <- 1
			for(i in (idx-nAR):1)
			{
				if(all(object$censorIndicator[i:(i+nAR-1)]==0))
				{
					iStart <- i
					break
				}
			} 

			nStart <- idx - iStart
			#message(sprintf("idx: %i, iStart: %i, nStart: %i",idx,iStart,nStart))
			tmpCensorIndicator <- object$censorIndicator[(idx-1):iStart] #looking back
			nCensored <- sum(tmpCensorIndicator!=0)
			covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart+1)
			if( nCensored < nStart )
			{
				conditionalIndex <- which(tmpCensorIndicator==0) + 1
				tmpY <- object$y[idx:iStart][conditionalIndex]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex, trend[idx:iStart], covEta)
				tmpMean <- cdist$'mean'
				tmpVar <- cdist$'var'
			}else 
			{
				tmpMean <- trend[idx:iStart]
				tmpVar <- covEta
			}

			tmpLower <- rep(-Inf,length = nCensored+1) #( y[idx], censored obs)
			tmpUpper <- rep(Inf,length = nCensored+1)
			censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
			tmpLower[-1][censored>0] <- object$upperCensorLimit[(idx-1):iStart][tmpCensorIndicator>0]
			tmpUpper[-1][censored<0] <- object$lowerCensorLimit[(idx-1):iStart][tmpCensorIndicator<0]
			mn <- mtmvnorm(mean=tmpMean, sigma=tmpVar,lower=tmpLower,upper=tmpUpper,doComputeVariance=FALSE)
			ret[idx] <- mn$tmean[1]
		}
	}
	ret
}
