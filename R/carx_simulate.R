#' simulates a sample data for \code{carx}.
#'
#' \code{simulate.carx} uses provided parameters and other settings to simulate a series of data.
#'
#' @param nObs number of observations to be simulated.
#' @param prmtrAR the AR parameter.
#' @param prmtrX the parameter for X.
#' @param sigmaEps the standard deviation for the white noises of the AR process.
#' @param lowerCensorLimit the lower censor limit.
#' @param upperCensorLimit the upper censor limit.
#' @param x optional matrix for X.
#' @param seed optional to set the seed of random number generator used by \code{R}.
#' @return a list of simulated \code{y}, \code{x}, \code{censorIndicator}, \code{lowerCensorLimit} and \code{upperCensorLimit}.
#' @export
#' @examples
#' nObs = 100
#' trueX = c(0.2,0.4)
#' trueAR = c(-0.28,0.25)
#' trueSigma = 0.60
#' lcl = -1
#' ucl = 1
#' dat = simulateCarx(100, trueAR, trueX, trueSigma, lcl, ucl, seed=0)


simulateCarx <- function(nObs, prmtrAR, prmtrX, sigmaEps, lowerCensorLimit, upperCensorLimit, x = NULL, seed=NULL)
{
	nAR <- length(prmtrAR)
	nX <- length(prmtrX)

	if(!is.null(seed))
		set.seed(seed)
	eps <- rnorm(nObs,0, sigmaEps)

	if(is.null(x))
  {
		x <- matrix(rnorm(nObs*nX), nrow= nObs, ncol = nX)
    #colnames(x) <- paste0("X",seq(1,length(nX)))
  }

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


