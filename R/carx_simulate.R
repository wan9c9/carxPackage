#' simulates a sample data for \code{carx}.
#'
#' \code{simulate.carx} uses provided parameters and other settings to simulate a series of data.
#'
#' @param nObs number of observations to be simulated.
#' @param prmtrAR the AR parameter.
#' @param prmtrX the parameter for X.
#' @param sigmaEps the standard deviation for the white noises of the AR process.
#' @param lcl the lower censor limit.
#' @param ucl the upper censor limit.
#' @param x optional matrix for X.
#' @param seed optional to set the seed of random number generator used by \code{R}.
#' @return a list of simulated \code{y}, \code{x}, \code{ci}, \code{lcl} and \code{ucl}.
#' @export
#' @examples
#' nObs = 100
#' trueX = c(0.2,0.4)
#' trueAR = c(-0.28,0.25)
#' trueSigma = 0.60
#' lcl = -1
#' ucl = 1
#' dat = carx.sim(100, trueAR, trueX, trueSigma, lcl, ucl, seed=0)
#' cts = carx.sim.cenTS(100, trueAR, trueX, trueSigma, lcl, ucl, seed=0)

#simulate <- function(nObs, prmtrAR, prmtrX, sigmaEps, lcl, ucl, x = NULL, seed=NULL) UseMethod("simulate")
carx.sim <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigmaEps=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL)
{
	p <- length(prmtrAR)
	nX <- length(prmtrX)

	if(!is.null(seed))
		set.seed(seed)
	eps <- rnorm(nObs,0, sigmaEps)

	if(is.null(x))
  {
		x <- matrix(rnorm(nObs*nX), nrow= nObs, ncol = nX)
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
		    x = x,
		    ci=ci,
		    lcl=lcl,
		    ucl=ucl
		    )
	ret
}



carx.sim.cenTS <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigmaEps=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL, end.date=Sys.Date())
{
  ret <- carx.sim(nObs,prmtrAR, prmtrX, sigmaEps, lcl, ucl, x, seed)
  
  listx <- lapply(seq_len(ncol(ret$x)), function(i) ret$x[,i])
  names(listx) <- colnames(ret$x)
  
  listx$order.by <- end.date+seq(-nObs,-1,by=1)
  listx$value <- ret$y
  listx$lcl <- ret$lcl
  listx$ucl <- ret$ucl
  
  val <- do.call(cenTS,listx)
  val 
}
