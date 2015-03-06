carx.simulate <- function(nObs, prmtrAR, prmtrX, sigmaEps, lowerCensorLimit, upperCensorLimit, x = NULL, seed=NULL)
{
	nAR <- length(prmtrAR)
	nX <- length(prmtrX)
  
	if(!is.null(seed)) 
		set.seed(seed)
	eps <- rnorm(nObs,0, sigmaEps)
  
	if(is.null(x))
		x <- matrix(rnorm(nObs*nX), nrow= nObs, ncol = nX)
	trend <- x%*%prmtrX

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
		    x = x,
		    censorIndicator=censorIndicator,
		    lowerCensorLimit=lowerCensorLimit,
		    upperCensorLimit=upperCensorLimit
		    )
	ret
}


