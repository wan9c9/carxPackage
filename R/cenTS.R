require(xts)

#' create censored time series class
#' create censored time series class \code{cenTS}, for each series, it consists of the values,
#' the name of which can be specified by the user, and by default is "value',
#'  and lower/upper censoring limit denoted by \code{lcl} and \code{ucl} respectively. 
#'  It can also store related variables in the \code{xreg} which is a \code{list}, right now all variables values are assumed to be of the same lengh of and thus aligned with the main censored time series. 
#' @param order.by the index vector
#' @param value the value vector
#' @param value.name the name of the value, default = "value"
#' @param lcl the vector of lower censoring limit 
#' @param ucl the vector of upper censoring limit
#' @return a \code{cenTS} object

cenTS <- function(value, order.by, value.name = "value", lcl=NULL,ucl=NULL,...)
{
  #step0: check ... variables
  xreg <- list(...)
  xregNames <- names(xreg)

  if(value.name %in% xregNames | "lcl" %in% xregNames | "ucl" %in% xregNames)
    stop(paste("Variable names '",value.name, "', 'lcl', and 'ucl' are reserved for cenTS, but there is at least one of them has(ve) appearing in the list of variables"))
  xreg <- data.frame(xreg)
  #for(x in xreg)
  #{
  #  if(length(x) != length(order.by))
  #    stop("variables in xreg must be of the same length as the time series!")
  #}
  ci <- rep(0,length(value))
  ci[!is.finite(value)] <- NA

  if(!is.null(lcl))
  {
    if(length(lcl)==1)
      lcl <- rep(lcl,length(value))
    idx <- is.finite(value) & is.finite(lcl)
    idx2 <- value[idx]<=lcl[idx]
    ci[idx][idx2] <- -1
    value[idx][idx2] <- lcl[idx][idx2]
  }
  
  if(!is.null(ucl))
  {
    if(length(ucl)==1)
      ucl <- rep(ucl,length(value))
    idx <- is.finite(value) & is.finite(ucl)
    idx2 <- value[idx] >= ucl[idx]
    ci[idx][idx2] <- 1
    value[idx][idx2] <- ucl[idx][idx2]
  }

  val <- data.frame(value.name=value)
  names(val) <- c(value.name)

  val$lcl = lcl
  val$ucl = ucl
  val$ci = ci
  #ret <- xts(val,order.by=order.by)
  ret <- xts(data.frame(val,xreg),order.by=order.by)


  attr(ret,"xreg") <- colnames(xreg)
  attr(ret,"value.name") <- value.name
  attr(ret,"censoring.rate") <- mean(abs(ci[is.finite(ci)]))

  class(ret) <- c('cenTS','xts','zoo')
  invisible(ret)
}

print.cenTS <- function (object)
{
  print(as.xts(object))
  cat(paste("censoring rate:",round(attributes(object)$censoring.rate,3),"\n"))
}


#' returns the \code{xreg} part of the \code{cenTS} object
#' @param object an \code{cenTS} object
#' @return the list in \code{xreg}
#' 
xreg <- function(object) UseMethod("xreg")

xreg.cenTS <- function(object)
{
  as.xts(object[,attributes(object)$xreg])
}


#' plot a \code{cenTS} object
#' 
plot.cenTS <- function(x, type = "l", auto.grid = TRUE, major.ticks = "auto", 
    minor.ticks = TRUE, major.format = TRUE, bar.col = "grey", 
    candle.col = "white", ann = TRUE, axes = TRUE, ...) 
{
    series.title <- deparse(substitute(x))
    value.name <- attributes(x)$value.name
    if(value.name != "value")
      series.title <- value.name
    
    ep <- axTicksByTime(x, major.ticks, format.labels = major.format)
    otype <- type
    
    xycoords <- xy.coords(.index(x), x[,value.name])
    ylim <- range(coredata(x),na.rm=TRUE)
    plot(xycoords$x, xycoords$y, type = type, axes = FALSE, ann = FALSE, ylim=ylim,...)

   
    if(!is.null(x$lcl)) 
    {
      lines(xycoords$x,x$lcl,col="red")
      idx <- is.finite(x[,value.name]) & is.finite(x$lcl)
      idx2 <- x[,value.name][idx] <= x$lcl[idx]
      points(xycoords$x[idx][idx2],x$lcl[idx][idx2],pch=2)
    }
    
    if(!is.null(x$ucl)) 
    {
      lines(xycoords$x,x$ucl,col="red")
      idx <- is.finite(x[,value.name]) & is.finite(x$ucl)
      idx2 <- x[,value.name][idx] >= x$ucl[idx]
      points(xycoords$x[idx][idx2],x$ucl[idx][idx2],pch=6)
    }

    if (auto.grid) {
        abline(v = xycoords$x[ep], col = "grey", lty = 4)
        grid(NA, NULL)
    }
    dots <- list(...)
    if (axes) {
        if (minor.ticks) 
            axis(1, at = xycoords$x, labels = FALSE, col = "#BBBBBB", 
                ...)
        axis(1, at = xycoords$x[ep], labels = names(ep), las = 1, 
            lwd = 1, mgp = c(3, 2, 0), ...)
        axis(2, ...)
    }
    box()
    if (!"main" %in% names(dots)) 
        title(main = series.title)
    do.call("title", list(...))
}



