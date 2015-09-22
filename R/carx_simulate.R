#' Simulate a sample data for \code{carx}
#'
#' Use provided parameters and other settings to simulate a series of data.
#'
#' @param nObs number of observations to be simulated.
#' @param prmtrAR the AR parameter.
#' @param prmtrX the parameter for X.
#' @param sigmaEps the standard deviation for the white noises of the AR process.
#' @param lcl the lower censor limit.
#' @param ucl the upper censor limit.
#' @param x optional matrix for X, default = \code{NULL}, in which case X will be simulated from standard normal distribution with dimensions determined by \code{nObs} and \code{prmtrX}.
#' @param seed optional to set the seed of random number generator used by \code{R}.
#' @return a data frame of simulated \code{y}, \code{x}, \code{ci}, \code{lcl} and \code{ucl}.
#' @export
#' @examples
#' dat = carxSim()

carxSim <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigmaEps=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL)
{
	p <- length(prmtrAR)
	nX <- length(prmtrX)

	if(!is.null(seed))
		set.seed(seed)
	eps <- stats::rnorm(nObs,0, sigmaEps)

	if(is.null(x))
  {
		x <- matrix(stats::rnorm(nObs*nX), nrow= nObs, ncol = nX)
  }

  if(is.null(colnames(x)))
    colnames(x) <- paste0("X",seq(1,nX))

	trend <- x%*%prmtrX

	eta <- numeric(nObs)
	y <- numeric(nObs)
	eta[1:p] <- eps[1:p]

	y[1:p] <- eps[1:p]
	for(i in (p+1):nObs){
		eta[i] <- eta[(i-1):(i-p)] %*% prmtrAR + eps[i]
		y[i] <- trend[i] + eta[i]
	}
	ci <- rep(0,nObs)
	ci[y<lcl] <- -1
	ci[y>ucl] <- 1

  if(options()$verbose) message(paste0("simulated series: censor rate: ", sum(abs(ci))/nObs))
	ret <- list(y = y,
		    ci=ci,
		    lcl=lcl,
		    ucl=ucl,
		    x = x
		    )
  ret <- data.frame(ret)
  colnames(ret) <- c("y","ci","lcl","ucl",colnames(x))
	ret
}


#' simulate a sample \code{\link{cenTS}} data for \code{carx}
#'
#' Use provided parameters and other settings to simulate a series of data as a \code{cenTS} object.
#'
#' @inheritParams carxSim
#' @param value.name the name of the response series
#' @param end.date the date of the last observation, default = \code{Sys.date()}.
#' @return a \code{cenTS} object with regressors.
#' @seealso \code{\link{carxSim}}.
#' @export
#' @examples
#' cts = carxSimCenTS()
carxSimCenTS <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigmaEps=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL, value.name = 'y', end.date=Sys.Date())
{
  ret <- carxSim(nObs,prmtrAR, prmtrX, sigmaEps, lcl, ucl, x, seed)
  #ret is a data.frame
  names(ret) <- c("value",names(ret)[-1])
  ret <- as.list(ret)
  ret$order.by <- end.date+seq(-nObs,-1,by=1)
  #listx$value <- ret$y
  ret$value.name <- value.name

  val <- do.call(cenTS,ret)
  val
}
