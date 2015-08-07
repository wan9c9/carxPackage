
#' detecting the outlier of the response data modelled by \code{carx}

ot.carx <- function(object)
{
	#message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/object$nObs
	eps <- rnorm(nSample,0,object$sigma)
	trend <- object$x%*%object$prmtrX
	covEta <- computeCovAR(object$prmtrAR, object$sigma,lag=object$p)
	nObs <- object$nObs
	p <- object$p
	prmtrAR <- object$prmtrAR
	skipIndex <- object$skipIndex
	y <- object$y
	ci <- object$ci
	lcl <- object$lcl
	ucl <- object$ucl

	y[ci>0] <- ucl[ci>0]
	y[ci<0] <- lcl[ci<0]

	pValues <- numeric(nObs)
	pValues[skipIndex] <- 1

	for(idx in seq(1,nObs)[-skipIndex])
	{
		#message(sprintf("checking %i",idx))
		wkm <- y[(idx-1):(idx-p)]
		tmpCensorIndicator <- ci[(idx-1):(idx-p)]
		nCensored <- sum(tmpCensorIndicator!=0)
    censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
		if(nCensored)   #at least one is censored
		{
			if( nCensored < p ) #not all are censored
			{ 
				conditionalIndex <- which(tmpCensorIndicator==0) #indices of known values
				tmpY <- y[(idx-1):(idx-p)][conditionalIndex] #known values
				tmpM <- trend[(idx-1):(idx-p)]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta)
				tmpMean <- cdist$'mean'
				tmpVar <- cdist$'var'
			}else{ 
				tmpMean <- trend[(idx-1):(idx-p)]
				tmpVar <- covEta
			}
      tmpLower <- rep(-Inf,length = nCensored)
      tmpUpper <- rep(Inf,length = nCensored)
      tmpLower[censored>0] <- object$ucl[(idx-1):(idx-p)][tmpCensorIndicator>0]
      tmpUpper[censored<0] <- object$lcl[(idx-1):(idx-p)][tmpCensorIndicator<0]
			smpl <- tmvtnorm::rtmvnorm(nSample,tmpMean,tmpVar,lower=tmpLower,upper=tmpUpper,algorithm="gibbs")
			smpl <- as.matrix(smpl)
			ySmpl <- numeric(nSample)
			for(i in 1:nSample){
				wkm[tmpCensorIndicator!=0] <- smpl[i,]
				ySmpl[i] <- trend[idx] + (wkm - trend[(idx-1):(idx-p)])%*%prmtrAR + eps[i]
			}
			pU <- sum(ySmpl > y[idx])/nSample
			pL <- sum(ySmpl < y[idx])/nSample
		}
		else{
			r <- y[idx]-trend[idx] - (wkm-trend[(idx-1):(idx-p)])%*%prmtrAR
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
		NULL
}

# We need to run the carx function and have the original object first. 
outlier <- function(object,data=NULL,oi.prefix="OI_") UseMethod("outlier")

outlier.carx <- function(object,data=NULL,oi.prefix="OI_")
{
  newObj = object
  newFormula = NULL
  tryCatch({newFormula=formula(object)},error=function(e){newFormula=NULL})

  oiVec = NULL
  ot = -1
  while(TRUE)
  {
    ot=ot.carx(object)
    if(is.null(ot) || (!is.null(oiVec) && ot %in% oiVec))
      break
    oi = numeric(newObj$nObs)
    oi[ot] = 1
    newVar = paste0(oi.prefix,ot)
    if(is.null(newFormula) | is.null(data))
    {
      newx=data.frame(newObj$x)
      newx[,newVar] <- oi
      newObj=carx(object$y,newx,
                 ci=object$ci,
                 lcl=object$lcl,
                 ucl=object$ucl,
                 p=object$p,
                 prmtrX = c(object$prmtrX,0),
                 prmtrAR = object$prmtrAR,
                 sigmaEps = object$sigma,
                 tol = object$tol,
                 max.iter = object$max.iter,
                 getCI=FALSE,
                 CI.level = object$CI.level,
                 b = object$b,
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
                   ci=object$ci,
                   lcl=object$lcl,
                   ucl=object$ucl,
                   p=object$p,
                   prmtrX = c(object$prmtrX,0),
                   prmtrAR = object$prmtrAR,
                   sigmaEps = object$sigma,
                   tol = object$tol,
                   max.iter = object$max.iter,
                   getCI=FALSE,
                   CI.level = object$CI.level,
                   b = object$b,
                   skipIndex=object$skipIndex,
                   verbose = FALSE
                   )
    }
    oiVec=c(oiVec,ot)
  }
  list("updatedModel"= newObj,"outlierIndices"=oiVec)
}
