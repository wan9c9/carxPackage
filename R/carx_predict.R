#' Provide predictions with fitted \code{carx} object
#'
#' \code{predict.carx} provides an method to predict the future values of an fitted
#' \code{carx} object with given new observations in \code{x}.
#' @param object A fitted \code{carx} object.
#' @param newdata The new observations for the coverates \code{x}.
#' If there is no covariates, the value can be assigned to be \code{NULL}.
#' Otherwise, an matrix of new observations is required for give predictions.
#' @param n.ahead The number of steps ahead the user wants to predict, default = 1.
#' @param level The level used to construct the Confidence interval, default = 0.95.
#' @param nRep The number of replications to be performed when censoring exists in the last \code{nAR}
#' observations, default = 1000.
#' @param ... not used.
#' @return A list consisting of \code{fit}, \code{se.fit}, and \code{ci} representing the predictions,
#' standard errors of predictions, and confidence intervals respectively.
#' @export
predict.carx <- function(object,newdata=NULL,n.ahead=1,level=0.95,nRep=1000,...)
{
  #tt <- terms(object)
  #if(!inherits(object,"carx")){
  #warning("calling predict.lm(<fake-carx-object>)...")
  #}
  #if(missing(newdata) || is.null(newdata)){
  #mm <- X <- model.matrix(object)
  #mmDone <- TRUE
  #}
  #else{
  #Terms <- delete.reponse(tt)
  #m <- model.frame(Terms,newdata,na.action=na.action

  if( is.null(newdata) && !object$xIsOne )
    stop("ERROR: newdata supplied is NULL, but the x data in model is not ones.")

  if( n.ahead < 1)
	  stop("ERROR: n.ahead must be greater than or equal to 1.")

  nAR <- object$nAR
  nObs <- object$nObs
  probs <- c((1-level)/2,(1+level)/2)
  qntl <- matrix(nrow=n.ahead,ncol=2)

  if(object$xIsOne) #no x, add intercept
  {
    newdata <- as.matrix(rep(1,n.ahead))
  }

  if(dim(newdata)[1] != n.ahead)
    stop("ERROR: number of rows in x doesn't equal to n.ahead.")

	#find the beginning index of nAR consecutive observations
	iStart <- 1
	for(i in (nObs-nAR+1):1)
	{
		if(all(object$censorIndicator[i:(i+nAR-1)]==0))
		{
			iStart <- i
			break
		}
	}

	nStart <- nObs - iStart + 1
	newdata <- rbind(object$x[iStart:nObs,],newdata)
	eta <- object$y[iStart:nObs] - object$x[iStart:nObs,]%*%object$prmtrX
	eta <- c(eta, rep(0,n.ahead))
  yPred <- c(object$y[iStart:nObs], rep(0,n.ahead))
	if(nStart == nAR)
	{
		#no censoring in latest nAR observations
	  message("no censoring in latest nAR obs")
		coefs <- matrix(rep(n.ahead*n.ahead),nrow=n.ahead,ncol=n.ahead)
		predSE <- rep(0,n.ahead)
		for(i in 1:n.ahead)
		{
			eta[iStart+i] <- eta[(iStart+i-1):i]%*%object$prmtrAR
			yPred[iStart+i] <- newdata[iStart+i,]%*%object$prmtrX + eta[iStart+i]
			coefs[i,i] <- 1
			if(i>1)
			{
				if( min(i-1,n.ahead-1,nAR) > 0)
				{
					for( j in 1:min(i-1,n.ahead-1,nAR) )
						coef[i,] <- coef[i,] + object$prmtrAR[j]*coef[i-j,]
				}
			}
			predSE[i] <- object$sigma^2*sum(coef[i,]^2)
		}
		yPred <- yPred[-(1:iStart)]
		q <- qnorm(probs)
		qntl[,1] <- yPred - predSE*q[1]
		qntl[,2] <- yPred - predSE*q[2]
			## prediction error?
		}
	else
	{
		message("censoring occurs in latest nAR obs")
		tmpCensorIndicator <- object$censorIndicator[nObs:iStart] #reverse order
		nCensored <- sum(tmpCensorIndicator!=0)
    covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart)
    trend <- as.vector(newdata[nStart:1,]%*%object$prmtrX)
    if( nCensored < nStart )
    {
      conditionalIndex <- which(tmpCensorIndicator==0)
      tmpY <- object$y[nObs:iStart][conditionalIndex]
      cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,trend,covEta)
      tmpMean <- cdist$'mean'
      tmpVar <- cdist$'var'
    }else{
      tmpMean <- trend
      tmpVar <- covEta
    }

    tmpLower <- rep(-Inf,length = nCensored)
    tmpUpper <- rep(Inf,length = nCensored)
    censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
    tmpLower[censored>0] <- object$upperCensorLimit[nObs:iStart][tmpCensorIndicator>0]
    tmpUpper[censored<0] <- object$lowerCensorLimit[nObs:iStart][tmpCensorIndicator<0]

    yCensored <- tmvtnorm::rtmvnorm(nRep,tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)
    eps <- matrix(rnorm(nRep*n.ahead,0,object$sigma),nrow=nRep,ncol=n.ahead)
    etaFuture <- matrix(nrow=nRep,ncol=nStart+n.ahead)

    for(iRep in 1:nRep)
    {
      etaFuture[iRep,nStart:1] <- object$y[nObs:iStart]
      etaFuture[iRep,nStart:1][tmpCensorIndicator!=0] <- yCensored[iRep,]
      etaFuture[iRep,nStart:1] <- etaFuture[iRep,nStart:1] - trend

      for(i in 1:n.ahead)
        etaFuture[iRep,nStart+i] <- etaFuture[iRep,(nStart+i-1):(nStart+i-nAR)]%*%object$prmtrAR + eps[iRep,i]
    }
    center <- newdata[-(1:nStart),]%*%object$prmtrX
    for(i in 1:n.ahead)
	    qntl[i,] <- quantile(etaFuture[,nStart+i],probs=probs)
    qntl[,1] <- qntl[,1] + center
    qntl[,2] <- qntl[,2] + center
    yPred <-  colMeans(etaFuture[,-(1:nStart)])
    yPred <-  yPred + center
    predSE <- matrixStats::colSds(etaFuture[,-(1:nStart)])
  }
  list("fit"=as.vector(yPred),"se.fit"=predSE,"ci"=qntl)
}

#debug(predict.carx)
