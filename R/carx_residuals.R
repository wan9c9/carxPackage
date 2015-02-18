#' get the residuals of a fitted \code{carx} object
residuals.carx <- function(object,...)
{
	message("Calling residuals.carx")
	nObs  <- object$nObs
	nAR  <- object$nAR

	rsdl <- rep(NA,nObs)
	trend <- object$x%*%object$prmtrX
	eta <- object$y - trend

	for(idx in (nAR+1):nObs)
	{
		#message(sprintf("calculating %i",idx))
		if(all(object$censorIndicator[idx:(idx-nAR)]==0))
		{
			rsdl[idx] <- (eta[idx] - eta[(idx-1):(idx-nAR)]%*%object$prmtrAR)/object$sigma
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

			currentCensorIndicator <- object$censorIndicator[idx]
			nStart <- idx - iStart
			#message(sprintf("idx: %i, iStart: %i, nStart: %i",idx,iStart,nStart))
			if(nStart == nAR && all(object$censorIndicator[(idx-1):iStart]==0)) #y[idx] is censored, but y[(idx-1):(idx-nAR)] are observed
			{ 
				mu <- eta[(idx-1):(idx-nAR)]%*%object$prmtrAR
				if(currentCensorIndicator > 0)
				{
					climit <- object$upperCensorLimit[idx] - trend[idx]
					val <- pnorm(climit, mu, object$sigma, lower.tail=FALSE)
					rsdl[idx] <- qnorm(runif(1,0,val),lower.tail=FALSE)
				}
				else
				{
					climit <- object$lowerCensorLimit[idx]- trend[idx]
					val <- pnorm(climit, mu, object$sigma, lower.tail=TRUE)
					rsdl[idx] <- qnorm(runif(1,0,val), lower.tail=TRUE)
				}
			} else 
			{
				#message("censoring occurs in latest nAR obs")
				tmpCensorIndicator <- object$censorIndicator[(idx-1):iStart] #reverse order
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

				if(currentCensorIndicator == 0)
				{ 
					pval <- ptmvnorm.marginal(object$y[idx], 1, mean=tmpMean, sigma=tmpVar,lower=tmpLower,upper=tmpUpper)
					#print(pval)
					rsdl[idx] <- qnorm(pval,lower.tail=TRUE)
					#print(rsdl[idx])
				}
				else
				{
					if(currentCensorIndicator > 0)
					{ 
						pval <- ptmvnorm.marginal(object$upperCensorLimit[idx], 1, mean=tmpMean, sigma=tmpVar,lower=tmpLower,upper=tmpUpper)
						rsdl[idx]  <- qnorm(runif(1, 0, 1-pval), lower.tail=FALSE)
					}
					else
					{
						pval <- ptmvnorm.marginal(object$lowerCensorLimit[idx], 1, mean=tmpMean, sigma=tmpVar,lower=tmpLower,upper=tmpUpper)
						rsdl[idx]  <- qnorm(runif(1, 0, pval), lower.tail=TRUE)
					}
				}
			}
		}
	}
	rsdl
}

#debug(residuals.carx)
