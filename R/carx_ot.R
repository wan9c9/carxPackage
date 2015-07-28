
#' detecting the outlier of the response data modelled by \code{carx}

ot.carx <- function(object)
{
	#message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/object$nObs
	eps <- rnorm(nSample,0,object$sigma)
	trend <- object$x%*%object$prmtrX
	covEta <- computeCovAR(object$prmtrAR, object$sigma,lag=object$nAR)
	nObs <- object$nObs
	nAR <- object$nAR
	prmtrAR <- object$prmtrAR
	skipIndex <- object$skipIndex
	y <- object$y
	censorIndicator <- object$censorIndicator
	lowerCensorLimit <- object$lowerCensorLimit
	upperCensorLimit <- object$upperCensorLimit

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
      tmpLower[censored>0] <- object$upperCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator>0]
      tmpUpper[censored<0] <- object$lowerCensorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator<0]
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
			r <- r/object$sigma
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

# We need to run the carx function and have the original object first. 
detect.outlier.carx <- function(object,data,oi.prefix="OI_")
{
  newObj=object
  newFormula = formula(object)
  oiVec=NULL
  ot=ot.carx(object)
  ot = 1
  # If data is not provided:
  if(missing(data))
  {
    while(ot !=-1) 
    {
      ot=ot.carx(object)
      if(!is.null(oiVec) & ot %in% oiVec)
        break
      oi=numeric(newObj$nObs)
      oi[ot]=1
      newVar=paste0(oi.prefix,ot)
      if(is.null(newFormula) | missing(data))
      {
      newx=data.matrix(data.frame(newObj$x,newVar=oi))
      newObj=carx(newObj$y,newx,
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
      else
      {
        newData = data
        newData[,newVar] <- oi
        newFormula = update(newFormula, as.formula(paste("~.+",newVar)))
		    newObj <- carx(newFormula, data = newData,
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
      oiVec=c(oiVec,ot)
    }
    list("updatedModel"= tmp,"outlierIndices"=oiVec)
 }
}



