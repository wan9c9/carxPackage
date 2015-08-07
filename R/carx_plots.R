#' plot a fitted \code{carx} object
#'
#' \code{plot.carx} plots a fitted \code{carx} object with other settings
#' the y axis will be the response, the x axis can be supplied by the user.
#' @param x a fitted \code{carx} object.
#' @param FUN an optional function to be applied to the data related to the
#' responses of the object, e.g., if a log transformation has been applied to the data, an exp function can be supplied so that the data plotted are at the same scale of the orginal data, default = \code{NULL}.
#' @param xAxisVar an optional vector to be plotted as the x variable, default =
#' \code{NULL} corresponds time series plot, other vectors are assumed to be
#' sequential.
#' @param xlab the label for x axis, default = "Index".
#' @param ylab the label for y axis, default = "Response".
#' @param outliers a \code{NULL} value or a vector, each element of which indicates a location at which the response is treated as an outlier, default = \code{NULL}, i.e., no outlier
#' @param ... other parameters supplied to the generic function \code{plot}.
#' @export
plot.carx <- function(x,FUN=identity, xAxisVar=NULL, xlab="Index", ylab="Response", outliers=NULL, ...)
{
	object <- x
	if(is.null(xAxisVar)) xAxisVar <- 1:length(object$y)
	xrange <- range(xAxisVar)
	lgd <- c(ylab,'Fitted value')
	plty <-c(1,2)
	pcol <-c('black','blue')

	yh <- sapply(fitted(object), FUN)
	y <- sapply(object$y,  FUN)

	validLcl <- ifelse(any(is.finite(object$lcl)), TRUE, FALSE)
	validUcl <- ifelse(any(is.finite(object$ucl)), TRUE, FALSE)
  
	if(validLcl)
	{
		lcl <- sapply(object$lcl, FUN)
		lgd <- c(lgd,"Lower censor limit")
		plty <- c(plty,3)
		pcol <- c(pcol,'red')
	}
	else
		lcl <- -Inf
	
	if(validUcl)
	{
	  ucl <- sapply(object$ucl, FUN)
		lgd <- c(lgd,"Upper censor limit")
		plty <- c(plty,5)
		pcol <- c(pcol,'red')
	}
	else
		ucl <- Inf

	ylim <- range(c(y,yh,lcl,ucl),na.rm=TRUE,finite=T)

	if( any(object$ci>0) )
		y[object$ci>0] <- ucl[object$ci>0]
	if( any(object$ci<0) )
		y[object$ci<0] <- lcl[object$ci<0]

	
	#plot(xAxisVar, y, type="n", xaxt="n", yaxt="n")
	#this.legend.size <- legend("topright",legend=lgd,lty=plty,col=pcol,plot=FALSE)
	#ylim[2] <- 1.04*(ylim[2] + this.legend.size$rect$h)
	ylim[2] <- 1.3*ylim[2]

	plot(xAxisVar, y, type="l", lty=1, xlab=xlab, ylab=ylab, ylim=ylim, col='black',...)
	lines(xAxisVar, yh, lty=2, col='blue')

	if(validLcl)
	{
		lines(xAxisVar, lcl, lty=3, col="red")
		points(xAxisVar[object$ci<0],lcl[object$ci<0],pch=2)
	}
	if(validUcl)
	{
		lines(xAxisVar, ucl, lty=5, col="red")
		points(xAxisVar[object$ci>0], ucl[object$ci>0],pch=6)
	}

	legend("topright",legend=lgd,lty=plty,col=pcol)
	if(!is.null(outliers)) abline(v=xAxisVar[outliers],col="red",lty=2)
}

#' Plotting the residual of a fitted \code{carx} object
#'
#' Plot a fitted \code{carx} object with other settings
#' the y axis will be the response, the x axis can be supplied by the user
#' @param object a fitted \code{carx} object.
#' @param residualType the type of the residual.
#' @param x indicates x variable, 
#' default = \code{NULL}, the time series plot of residuals will be plotted; 
#' if x is "fitted", the residuals against the fitted values; 
#' if x is a vector, the residuals against the supplied vector.
#' @param ... other parameters to be supplied to \code{plot}.
#' @export

plotResiduals.carx <- function(object, residualType="pearson", x=NULL,xlab="Index",ylab="Residual",...)
{
	y <- residuals(object,residualType)
	if(is.null(x))
	{
		x <- 1:length(y)
		xlab <- "Index"
	}else{
		if(typeof(x) == "character")
		{
			if(x == 'fitted'){
				x <- fitted(object)
				xlab <- "Fitted value"
			}
		}
	}
	plot(x,y,xlab=xlab,ylab=ylab,...)
	abline(h=0,col="black",lty=3)
}

#' Plotting the ACF of the residuals of a fitted \code{carx} object
#'
#' Plot the ACF of the residuals of a fitted \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param ... other parameters to be supplied to \code{acf}.
#' @return no value is returned, a figure will either be sent to display or be saved.
#' @export
plotResAcf.carx <- function(object,...)
{
	acf(residuals(object),na.action=na.pass,...)
}


tsdiag <- function(object,...) UseMethod("tsdiag")

tsdiag.carx <- function(object,gof.lag,tol=0.1,col="red",omit.initial=TRUE,...)
{
  opar = par(mfrow = c(3, 1), mar = c(4, 4, 4, 3) + 0.1, oma = c(1, 0, 2, 0))
  n = object$nObs
  if (missing(gof.lag)) 
      lag.max = 10 * log10(n)
  else 
    lag.max = gof.lag

  residuals = residuals(object,type="raw")
  if (omit.initial) 
      residuals = window(residuals, start = time(residuals)[object$p + 1])
  std.res = residuals/object$sigma
  n = length(std.res)
  h1 = qnorm(0.025/n)
  plot(std.res, ylab = "Standardized Residuals", type = "p", ...)
  abline(h = h1, lty = 2, col = col)
  abline(h = -h1, lty = 2, col = col)
  abline(h = 0)
  acf(as.numeric(residuals), lag.max = lag.max, ylab = "ACF of Residuals", 
      ci.col = col, main = "", ...)
  lbv = rep(NA, lag.max)
  k = object$p+1
  for (i in k:lag.max) {
      lbv[i] = Box.test(std.res, lag = i, type="Ljung-Box",fitdf = object$p)$p.value
  }
  plot(y = lbv, x = 1:lag.max, ylim = c(0, 1), pch = 21, ylab = "P-values", 
      xlab = "Number of lags", ...)
  abline(h = 0.05, lty = 2, col = col)
  par(opar)
  invisible()
}
