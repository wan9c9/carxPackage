
#' detecting the outlier of the response data modelled by \code{carx}

ot.carx <- function(model){
	#message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/model$nObs
	eps <- rnorm(nSample,0,model$sigma)
	trend <- model$x%*%model$prmtrX
	covEta <- computeCovAR(model$prmtrAR, model$sigma,lag=model$nAR)
	nObs <- model$nObs
	nAR <- model$nAR
	prmtrAR <- model$prmtrAR
	skipIndex <- model$skipIndex
	y <- model$y
	censorIndicator <- model$censorIndicator
	lowerCensorLimit <- model$lowerCensorLimit
	upperCensorLimit <- model$upperCensorLimit

	y[censorIndicator>0] <- upperCensorLimit[censorIndicator>0]
	y[censorIndicator<0] <- lowerCensorLimit[censorIndicator<0]

	pValues <- numeric(nObs)
	pValues[skipIndex] <- 1

	for(idx in seq(1,nObs)[-skipIndex])
	{
		#message(sprintf("checking %i",idx))
		wkm <- y[(idx-1):(idx-nAR)]
		tmpCensorIndicator <- censorIndicator[(idx-1):(idx-nAR)]
		nCensored <- sum(tmpCensorIndicator!=0)
    censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
		if(nCensored)   #at least one is censored
		{
			if( nCensored < nAR ) #not all are censored
			{ 
				conditionalIndex <- which(tmpCensorIndicator==0) #indices of known values
				tmpY <- y[(idx-1):(idx-nAR)][conditionalIndex] #known values
				tmpM <- trend[(idx-1):(idx-nAR)]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta)
				tmpMean <- cdist$'mean'
				tmpVar <- cdist$'var'
			}else{ 
				tmpMean <- trend[(idx-1):(idx-nAR)]
				tmpVar <- covEta
			}
      tmpLower <- rep(-Inf,length = nCensored)
      tmpUpper <- rep(Inf,length = nCensored)
      tmpLower[censored>0] <- model$upperCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator>0]
      tmpUpper[censored<0] <- model$lowerCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator<0]
			smpl <- tmvtnorm::rtmvnorm(nSample,tmpMean,tmpVar,lower=tmpLower,upper=tmpUpper,algorithm="gibbs")
			smpl <- as.matrix(smpl)
			ySmpl <- numeric(nSample)
			for(i in 1:nSample){
				wkm[tmpCensorIndicator!=0] <- smpl[i,]
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



outlierDetection.carx <- function(object)
{
  ot <- 1
  nOL <- 0
  olv <- NULL

  while(ot>0)
  {
    ot <- ot.carx(object)
    if(ot != -1) # find an outlier
    { 
      if(nOL == 0) # no outlier yet
        olv <- ot
      else
      {
        if(is.vector(olv) && ot %in% olv)
          break
        olv <- c(olv,ot)
      }
      nOL <- nOL + 1
      oi <- numeric(object$nObs)
      oi[ot] <- 1
      newVar <-paste0("OutlierIndicator_",ot)
      data1 <- object$x
      data1[,newVar] <- oi
			newFormula <- update(formula(object), as.formula(paste("~.+",newVar)))
		  newObj <- carx(newFormula, data = data1, 
                     censorIndicator=object$censorIndicator,
                     lowerCensorLimit=object$lowerCensorLimit,
                     upperCensorLimit=object$upperCensorLimit,
                     nAR=object$nAR,
                     prmtrX = object$prmtrX,
                     prmtrAR = object$prmtrAR,
                     sigmaEps = object$sigma,
                     tol = object$tol,
                     max.iter = object$max.iter,
                     getCI=FALSE,
                     alpha = object$alpha,
                     nBootstrapSample = object$nBootstrapSample,
                     skipIndex=object$skipIndex,
                     verbose = FALSE
                     )
    }
  }
  newObj
}
