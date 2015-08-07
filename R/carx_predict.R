
#' compute the coefficients of MA representation of an AR model
#'

predictARX <- function(y,x,prmtrX,prmtrAR,sigma,n.ahead,newxreg,CI.level=0.95)
{
  nObs <- length(y)
  p <- length(prmtrAR)
  stopifnot( nObs >= p )
  iStart <- nObs-p+1 #use only latest p observations
  nStart <- p 
  x <- as.matrix(x)
  newxreg <- as.matrix(newxreg)
	newxreg <- rbind(as.matrix(x[iStart:nObs,]),newxreg)
	eta <- y[iStart:nObs] - as.matrix(x[iStart:nObs,])%*%prmtrX
	eta <- c(eta, rep(0,n.ahead))
  yPred <- c(y[iStart:nObs], rep(0,n.ahead))
  coefs <- matrix(rep(0,n.ahead*n.ahead),nrow=n.ahead,ncol=n.ahead)
  predSE <- rep(0,n.ahead)
  for(i in 1:n.ahead)
  {
    eta[nStart+i] <- sum(eta[(nStart+i-1):(nStart+i-p)]*prmtrAR)
    yPred[nStart+i] <- newxreg[nStart+i,]%*%prmtrX + eta[nStart+i]
    coefs[i,i] <- 1
    if(i>1)
    {
      if( min(i-1,n.ahead-1,p) > 0)
      {
        for( j in 1:min(i-1,n.ahead-1,p) )
          coefs[i,] <- coefs[i,] + prmtrAR[j]*coefs[i-j,]
      }
    }
    predSE[i] <- sigma^2*sum(coefs[i,]^2)
  }
  yPred <- yPred[-(1:p)]
  probs <- c((1-CI.level)/2,(1+CI.level)/2)
  q <- qnorm(probs)
  qntl <- matrix(nrow=n.ahead,ncol=2)
  qntl[,1] <- yPred + predSE*q[1]
  qntl[,2] <- yPred + predSE*q[2]
  list("fit"=as.vector(yPred),"se.fit"=predSE,"ci"=qntl)
}


#' Provide predictions with fitted \code{carx} object
#'
#' \code{predict.carx} provides an method to predict the future values of an fitted
#' \code{carx} object with given new observations in \code{x}.
#' @param object A fitted \code{carx} object.
#' @param newxreg The new observations for the coverates \code{x}.
#' If there is no covariates, the value can be assigned to be \code{NULL}.
#' Otherwise, an matrix of new observations is required for give predictions.
#' @param n.ahead The number of steps ahead the user wants to predict, default = 1.
#' @param CI.level The CI.level used to construct the Confidence interval, default = 0.95.
#' @param nRep The number of replications to be performed when censoring exists in the last \code{p}
#' observations, default = 1000.
#' @param ... not used.
#' @return A list consisting of \code{fit}, \code{se.fit}, and \code{ci} representing the predictions,
#' standard errors of predictions, and confidence intervals respectively.
#' @export
predict.carx <- function(object,newxreg=NULL,n.ahead=1,CI.level=0.95,nRep=1000,na.action=NULL,...)
{
  
  if( n.ahead < 1)
	  stop("ERROR: n.ahead must be greater than or equal to 1.")

  if(object$xIsOne) #no x, add intercept
  {
    newxreg <- as.matrix(rep(1,n.ahead))
  }
  else
  {
    if(is.null(newxreg))
      stop("ERROR: newxreg supplied is NULL, but the x data in model is not ones.")
    
    newxreg <- as.matrix(newxreg)
    if(dim(newxreg)[1] != n.ahead)
        stop("New data doesn't have the same row of data as n.ahead.")

    #if formula is in the model, the user is more likely to supply the raw data,
    tryCatch({fml=formula(object)},error=function(e){fml=NULL})
    if(is.null(fml))
    {
      #no formula, the new data must have the same dimension as x
      #if(is.vector(newxreg))
      if(dim(newxreg)[2] != dim(object$x)[2])
        stop("New data doesn't have the same number of variables as object$x.")
    }
    else
    {
      if(dim(newxreg)[2] != dim(object$x)[2])
      {
        message("I am trying to combine the formula and your supplied data.")
        mf <- model.frame(formula=getCovariateFormula(fml),data=as.data.frame(newxreg),na.action=NULL)
        newxreg <- model.matrix(attr(mf,"terms"),data=mf)

      }
    }
  }
 
  p <- object$p
  nObs <- object$nObs
  probs <- c((1-CI.level)/2,(1+CI.level)/2)
  qntl <- matrix(nrow=n.ahead,ncol=2)

  if(dim(newxreg)[1] != n.ahead)
    stop("ERROR: number of rows in x doesn't equal to n.ahead.")

	#find the beginning index of p consecutive observations
	iStart <- 1
	for(i in (nObs-p+1):1)
	{
		if(all(object$ci[i:(i+p-1)]==0))
		{
			iStart <- i
			break
		}
	}

  stopifnot(all(is.finite(object$y[iStart:nObs]))) #this means there exists na values from iStart which we cannot handle right now.

	nStart <- nObs - iStart + 1
	
	if(nStart == p)
	{
	  message("latest p observations are not censored")
    predictARX(object$y,object$x,object$prmtrX,object$prmtrAR,object$sigma,n.ahead,newxreg,CI.level=CI.level)
  }
	else
	{
	  message("latest p observations are censored")
	  eta <- object$y[iStart:nObs] - object$x[iStart:nObs,]%*%object$prmtrX
  	eta <- c(eta, rep(0,n.ahead))
    yPred <- c(object$y[iStart:nObs], rep(0,n.ahead))
    predSE <- rep(0,n.ahead)
	  newxreg <- rbind(object$x[iStart:nObs,],newxreg)
		tmpCensorIndicator <- object$ci[nObs:iStart] #reverse order
		nCensored <- sum(tmpCensorIndicator!=0)
    covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart)
    trend <- as.vector(newxreg[nStart:1,]%*%object$prmtrX)
    if( nCensored < nStart )
    {
      conditionalIndex <- which(tmpCensorIndicator==0)
      tmpY <- object$y[nObs:iStart][conditionalIndex]
      cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,trend,covEta)
      tmpMean <- cdist$'mean'
      tmpVar <- cdist$'var'
    }else
    {
      tmpMean <- trend
      tmpVar <- covEta
    }

    tmpLower <- rep(-Inf,length = nCensored)
    tmpUpper <- rep(Inf,length = nCensored)
    censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
    tmpLower[censored>0] <- object$ucl[nObs:iStart][tmpCensorIndicator>0]
    tmpUpper[censored<0] <- object$lcl[nObs:iStart][tmpCensorIndicator<0]

    yCensored <- tmvtnorm::rtmvnorm(nRep,tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)
    eps <- matrix(rnorm(nRep*n.ahead,0,object$sigma),nrow=nRep,ncol=n.ahead)
    etaFuture <- matrix(nrow=nRep,ncol=nStart+n.ahead)

    for(iRep in 1:nRep)
    {
      etaFuture[iRep,nStart:1] <- eta[nStart:1]
      etaFuture[iRep,nStart:1][tmpCensorIndicator!=0] <- yCensored[iRep,] - trend[tmpCensorIndicator!=0]

      for(i in 1:n.ahead)
        etaFuture[iRep,nStart+i] <- etaFuture[iRep,(nStart+i-1):(nStart+i-p)]%*%object$prmtrAR + eps[iRep,i]
    }
    center <- newxreg[-(1:nStart),]%*%object$prmtrX
    for(i in 1:n.ahead)
	    qntl[i,] <- quantile(etaFuture[,nStart+i],probs=probs)
    qntl[,1] <- qntl[,1] + center
    qntl[,2] <- qntl[,2] + center
    yPred <-  colMeans(etaFuture[,-(1:nStart)])
    yPred <-  yPred + center
    predSE <- matrixStats::colSds(etaFuture[,-(1:nStart)])
    list("fit"=as.vector(yPred),"se.fit"=predSE,"ci"=qntl)
  }
}

#debug(predict.carx)
