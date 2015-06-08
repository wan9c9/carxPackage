#' plot a fitted \code{carx} object
#'
#' \code{plot.carx} plots a fitted \code{carx} object with other settings
#' the y axis will be the response, the x axis can be supplied by the user.
#' @param x a fitted \code{carx} object.
#' @param Fun an optional function to be applied to the data related to the
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

	validLcl <- ifelse(any(is.finite(object$lowerCensorLimit)), TRUE, FALSE)
	validUcl <- ifelse(any(is.finite(object$upperCensorLimit)), TRUE, FALSE)
  
	if(validLcl)
	{
		lcl <- sapply(object$lowerCensorLimit, FUN)
		lgd <- c(lgd,"Lower censor limit")
		plty <- c(plty,3)
		pcol <- c(pcol,'red')
	}
	else
		lcl <- -Inf
	
	if(validUcl)
	{
	  ucl <- sapply(object$upperCensorLimit, FUN)
		lgd <- c(lgd,"Upper censor limit")
		plty <- c(plty,5)
		pcol <- c(pcol,'red')
	}
	else
		ucl <- Inf

	ylim <- range(c(y,yh,lcl,ucl),na.rm=TRUE,finite=T)

	if( any(object$censorIndicator>0) )
		y[object$censorIndicator>0] <- ucl[object$censorIndicator>0]
	if( any(object$censorIndicator<0) )
		y[object$censorIndicator<0] <- lcl[object$censorIndicator<0]

	
	plot(xAxisVar, y, type="n", xaxt="n", yaxt="n")
	this.legend.size <- legend("topright",legend=lgd,lty=plty,col=pcol,plot=FALSE)
	ylim[2] <- 1.04*(ylim[2] + this.legend.size$rect$h)

	plot(xAxisVar, y, type="l", lty=1, xlab=xlab, ylab=ylab, ylim=ylim, col='black',...)
	lines(xAxisVar, yh, lty=2, col='blue')

	if(validLcl)
	{
		lines(xAxisVar, lcl, lty=3, col="red")
		points(xAxisVar[object$censorIndicator<0],lcl[object$censorIndicator<0],pch=2)
	}
	if(validUcl)
	{
		lines(xAxisVar, ucl, lty=5, col="red")
		points(xAxisVar[object$censorIndicator>0], ucl[object$censorIndicator>0],pch=6)
	}

	legend("topright",legend=lgd,lty=plty,col=pcol)
	if(!is.null(outliers)) abline(v=xAxisVar[outliers],col="red",lty=2)
}

#' Plotting the residual of a fitted \code{carx} object
#'
#' Plot a fitted \code{carx} object with other settings
#' the y axis will be the response, the x axis can be supplied by the user
#' @param object a fitted \code{carx} object.
#' @param x indicates x variable, default = \code{NULL}, it will plot the time series plot the residuals; if x is "fitted", it will plot the residuals against the fitted values; if x is a vector, it will plot the residual against the vector.
#' @param xlab the label for x axis, default = "Index".
#' @param ylab the label for y axis, default = "Residuals".
#' @param saveFig indicates whether to save the figure, default = \code{NULL} indicates the plot will be sent to the display, otherwise, saveFig need to be a file name to save the file, with postfix ".eps", ".jpg" or ".svg".
#' @param width the width of the figure, if to be saved.
#' @param height the height of the figure, if to be saved.
#' @param units the units width and height of the figure.
#' @param res the resolution of the figure, only needed if the figure is to be saved in "jpg" format.
#' @param ... other parameters to be supplied to \code{residuals.carx} or \code{plot}.
#' @return no value is returned, a figure will either be sent to display or be saved.
#' @export

plotResiduals.carx <- function(object,x=NULL,xlab="",ylab="",saveFig=NULL,width=6,height=5,units="in",res=300,...)
{
	if(!is.null(saveFig))
	{
		if(substring(saveFig,nchar(saveFig)-3) == ".eps")
		{
			setEPS()
			postscript(saveFig)
		}
		else if (substring(saveFig,nchar(saveFig)-3) == ".svg")
			svg(saveFig)
    else if (substring(saveFig,nchar(saveFig)-3) == ".jpg")
			jpeg(saveFig,width=width,height=height,units=units,res=res,quality=100)
	}

	y <- residuals(object,...)
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
  if(!is.null(saveFig))
		dev.off()
}

#' Plotting the ACF of the residuals of a fitted \code{carx} object
#'
#' Plot the ACF of the residuals of a fitted \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param saveFig indicates whether to save the figure, default = \code{NULL} indicates the plot will be sent to the display, otherwise, saveFig need to be a file name to save the file, with postfix ".eps", ".jpg" or ".svg".
#' @param width the width of the figure, if to be saved.
#' @param height the height of the figure, if to be saved.
#' @param units the units width and height of the figure.
#' @param res the resolution of the figure, only needed if the figure is to be saved in "jpg" format.
#' @param main the main title in the plot.
#' @param ... other parameters to be supplied to \code{residuals.carx} or \code{plot}.
#' @return no value is returned, a figure will either be sent to display or be saved.
#' @export
plotResAcf <- function(object,saveFig=NULL,width=6,height=5,units="in",res=300,main="",...)
{
	if(!is.null(saveFig))
	{
		if(substring(saveFig,nchar(saveFig)-3) == ".eps")
		{
			setEPS()
			postscript(saveFig)
		}
		else if (substring(saveFig,nchar(saveFig)-3) == ".svg")
			svg(saveFig)
    else if (substring(saveFig,nchar(saveFig)-3) == ".jpg")
			jpeg(saveFig,width=width,height=height,units=units,res=res,quality=100)
	}
	acf(residuals(object,...),na.action=na.pass,main=main,...)
  if(!is.null(saveFig))
		dev.off()
}

