
plot.carx <- function(object,transformFun=NULL,xAxisVar=NULL,xlab="Index",ylab="Observations",saveFig=NULL, outliers=NULL,...)
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
	}

	yh <- fitted(object)
	y <- object$y
	lcl <- object$lowerCensorLimit
	ucl <- object$upperCensorLimit
	ylim <- range(c(y,yh),na.rm=TRUE)

	#y[object$censorIndicator] <- NA
	if( any(object$censorIndicator>0) ) 
		y[object$censorIndicator>0] <- object$upperCensorLimit[object$censorIndicator>0]
	if( any(object$censorIndicator<0) ) 
		y[object$censorIndicator<0] <- object$lowerCensorLimit[object$censorIndicator<0]

	if(is.null(xAxisVar))
	{
		xAxisVar <- 1:length(yh)
	}

	xrange <- range(xAxisVar)
	if(is.null(transformFun))
	{
		plot(as.zoo(as.ts(zoo(y, xAxisVar))), type="n", xaxt="n", yaxt="n")
		this.legend.size <- legend("topright",legend=c(ylab,'Fitted value','Lower censor limit','Upper censor limit'),lty=c(1,2,3,5),col=c('black','blue','red','red'),plot=FALSE)
		ylim <- 1.04*(ylim + this.legend.size$rect$h)

		plot(as.zoo(as.ts(zoo(y, xAxisVar))), lty=1,xlab=xlab,ylab=ylab,ylim=ylim,col='black',...)
		lines(as.zoo(as.ts(zoo(yh, xAxisVar))), lty=2,col='blue')

		lines(as.zoo(as.ts(zoo(object$lowerCensorLimit, xAxisVar))),lty=3,col="red")
		lines(as.zoo(as.ts(zoo(object$upperCensorLimit, xAxisVar))),lty=5,col="red")
		points(xAxisVar[object$censorIndicator<0],object$lowerCensorLimit[object$censorIndicator<0],pch=2)
		points(xAxisVar[object$censorIndicator>0],object$upperCensorLimit[object$censorIndicator>0],pch=6)
	}else{
		ylim <- transformFun(ylim)
		plot(as.zoo(as.ts(zoo(transformFun(y), xAxisVar))), type="n",xaxt="n",yaxt="n")
		this.legend.size <- legend("topright",legend=c(ylab,'Fitted value','Lower censor limit','Upper censor limit'),lty=c(1,2,3,5),col=c('black','blue','red','red'),plot=FALSE)
		ylim <- 1.04*(ylim + this.legend.size$rect$h)
		plot(as.zoo(as.ts(zoo(transformFun(y), xAxisVar))), lty=1,xlab=xlab,ylab=ylab,ylim=ylim,...)
		lines(as.zoo(as.ts(zoo(transformFun(yh), xAxisVar))), lty=2,col='blue')
		lines(as.zoo(as.ts(zoo(transformFun(object$lowerCensorLimit), xAxisVar))),lty=3,col="red")
		lines(as.zoo(as.ts(zoo(transformFun(object$upperCensorLimit), xAxisVar))),lty=5,col="red")
		points(xAxisVar[object$censorIndicator<0],transformFun(object$lowerCensorLimit[object$censorIndicator<0]),pch=2)
		points(xAxisVar[object$censorIndicator>0],transformFun(object$upperCensorLimit[object$censorIndicator>0]),pch=6)
	}
	legend("topright",legend=c(ylab,'Fitted value','Lower censor limit','Upper censor limit'),lty=c(1,2,3,5),col=c('black','blue','red','red'))
	if(!is.null(outliers)) abline(v=xAxisVar[outliers],col="red",lty=2)
	if(!is.null(saveFig))
		dev.off()
}

plot.residuals <- function(object,x=NULL,saveFig="",xlab="",ylab="",classify.by=NULL,type="l",lty=1)
{
	#if(saveFig != "") jpeg(saveFig)
	y <- residuals(object)/object$sigma
	if(is.null(x))
	{
		x <- 1:length(y)
		xlab <- "index"
	}else{
		if(typeof(x) == "character")
		{
			if(x == 'fitted'){
				x <- predict(object)
				xlab <- "Fitted value"
			}
		}
	}
	plot(x,y,xlab=xlab,ylab=ylab,type=type,lty=lty)
	abline(h=0,col="black",lty=3)
	#abline(h=1.96*object$sigma,col="red",lty=3)
	#abline(h=-1.96*object$sigma,col="red",lty=3)

	if(saveFig != "")
	{
		dev.copy2eps(file=saveFig)
		dev.off()
	}
}

plotResAcf <- function(object,saveFig="")
{
	setEPS()
	postscript(saveFig)
	acf(residuals(object),na.action=na.pass,main="")
	dev.off()
	#if(saveFig != "")
	#{
	#dev.copy2eps(file=saveFig)
	#dev.off()
	#}
}

plotData <- function( carxData,timeAxis=NULL)
{
	y <- object$y
	lcl <- object$lowerCensorLimit
	ucl <- object$upperCensorLimit

	t <- seq(1,length(y))
	ts.plot(t,y,col="black")
	lines(t,lcl,col="red")
	lines(t,ucl,col="red")

}
