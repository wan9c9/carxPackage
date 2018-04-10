#' Create a censored time series object of \code{cenTS} class
#'
#' Create a censored time series response object of \code{cenTS} class. Default name of the response is "value",
#'  with the vectors of lower/upper censoring limits denoted by \code{lcl} and \code{ucl} respectively.
#'  The vector of censoring indicators, i.e., \code{ci}, is part of the \code{cenTS} object.
#'  Additional related variables can be stored and provided in the construction function, whose names
#'  are stored in \code{xreg}. All variable values are assumed to be of the same length of and thus
#'  aligned with the censored response time series. \code{cenTS} inherits from [xts::xts].
#' @param order.by the index vector, must be a vector of time/date.
#' @param value the value vector.
#' @param lcl the vector of lower censoring limits, or a single numeric representing the constant limit.
#'  Default = \code{NULL} indicating no lower limit.
#' @param ucl the vector of upper censoring limits, or a single numeric representing the constant limit.
#'  Default = \code{NULL} indicating no upper limit.
#' @param ci the vector of censoring indicators whose value is "N"("L","R","I")
#' if the corresponding response is not censored  (left, right, interval censored).
#' Default = \code{NULL}, in which case, the function will compute \code{ci} by \code{value}, \code{lcl}
#' and \code{ucl}. It is the user's responsibility to check the consistency of the data.
# If \code{ci} is not \code{NULL}, the function will check the consistency of the data,
# assuming the observed values less (greater) than or equal to left (right) censoring limits are censored,
# and are observed otherwise. The function will stop if inconsistent results are found.
#' @param value.name the name of the value, default = "value".
#' @param ... additional variables, must be able to be coerced to a \code{data.frame}.
#' @return a \code{cenTS} object, any censored observation will be replaced by its corresponding censoring limit.
#' @export

#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  print(ts)
#'  print(xreg(ts))
#'  plot(ts)
#'
# \dontrun{
# #wrong call, case 1
#ts <- cenTS(value=c(1,-2,1,NA,0),
#            order.by=as.Date(strDates,"%Y-%m-%d"),
#            lcl=c(-3,-2,-1,-1,0),
#            ucl=c(3,2,1,1,1),
#           ci =c(1,-1,1,NA,-1)
#            ci =c("R","L","R","N","L")
#)
##wrong call, case 2
#ts <- cenTS(value=c(1,-2,1,NA,0),
#            order.by=as.Date(strDates,"%Y-%m-%d"),
#            lcl=c(-3,-2,-1,-1,0),
#            ucl=c(3,2,1,1,1),
#           ci =c(1,-1,1,NA,-1)
#            ci =c("R","L","R","N","L")
#)
#
#
##wrong call, case 3
#ts <- cenTS(value=c(1,-2,1,NA,0),
#            order.by=as.Date(strDates,"%Y-%m-%d"),
#            lcl=c(-3,-2,-1,-1,0),
#            ucl=c(3,2,1,1,1),
#           ci =c(0,-1,0,NA,-1)
#            ci =c("N","L","R","N","L")
#)
# }
#
cenTS <- function(value, order.by,
                  lcl = NULL,ucl = NULL,
                  ci = NULL,
                  value.name = "value",
                  ...)
{
  #step0: check ... variables
  xreg <- list(...)
  if(length(xreg)>0)
  {
    xregNames <- names(xreg)

    if(value.name %in% xregNames | "lcl" %in% xregNames | "ucl" %in% xregNames)
      stop(paste("Variable names '",value.name, "', 'lcl', and 'ucl' are reserved for cenTS, but there is at least one of them has(ve) appeared in the list of ... variables."))
    xreg <- data.frame(xreg)
  }
  else
    xreg <- NULL
  #for(x in xreg)
  #{
  #  if(length(x) != length(order.by))
  #    stop("variables in xreg must be of the same length as the time series!")
  #}
  hasCI <- TRUE
  if(is.null(ci))
  {
    ci <- rep("N",length(value))
    hasCI <- FALSE
  }
  ci[!is.finite(value)] <- "N"

  #if(!is.null(lcl))
  #{
  #  if(length(lcl)==1)
  #    lcl <- rep(lcl,length(value))
  #  #check left censored data
  #  idx <- is.finite(value) & ci=="L"
  #  if(length(idx)>1)
  #  {
  #    idx2 <- value[idx] <= lcl[idx] #find sub-index
  #    if(hasCI)
  #    {
  #      test <- idx2 != (ci[idx]=="L")
  #      if(any(test))
  #        stop("Inconsistency found in data, at index ",seq(1,length(value))[idx][test])
  #    }
  #    else
  #      ci[idx][idx2] <- "L"
  #    value[idx][idx2] <- lcl[idx][idx2]
  #  }
  #}

  #if(!is.null(ucl))
  #{
  #  if(length(ucl)==1)
  #    ucl <- rep(ucl,length(value))
  #  idx <- is.finite(value) & is.finite(ucl)
  #  idx2 <- value[idx] >= ucl[idx]
  #  if(hasCI)
  #  {
  #    test <- idx2 != (ci[idx]=="R")
  #    if(any(test))
  #      stop("Inconsistency found in data, at index ",seq(1,length(value))[idx][test])
  #  }
  #  else
  #    ci[idx][idx2] <- "R"
  #  value[idx][idx2] <- ucl[idx][idx2]
  #}

  val <- data.frame(value.name=value)
  names(val) <- c(value.name)

  val$lcl = lcl
  val$ucl = ucl
  #val$ci = ci
  if(is.null(xreg))
  {
    ret <- xts::xts(data.frame(val),order.by=order.by)
    attr(ret,"xreg") <- NULL
  }
  else
  {
    ret <- xts::xts(data.frame(val,xreg),order.by=order.by)
    attr(ret,"xreg") <- colnames(xreg)
  }
  #ret = merge(ret,ci)
  #colnames(ret) = c(colnames(ret),"ci")
  attr(ret,"ci") = as.character(ci)
  
  attr(ret,"value.name") <- value.name
  attr(ret,"censoring.rate") <- mean((ci!="N")[is.finite(ret[,value.name])])
  #attr(ret,"data") = cbind(val,ci)
  class(ret) <- c('cenTS','xts','zoo')
  invisible(ret)
}

#' Print a \code{cenTS} object
#' @param x a \code{cenTS} object.
#' @param ... not used.
#' @return none.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  print(ts)
#'

print.cenTS <- function (x,...)
{
  #print(rbind(xts::as.xts(x),attributes(x)$ci))
  if("lcl" %in% names(x) | "ucl" %in% names(x))
  {
    ci = attributes(x)$ci 
    print(cbind(data.frame(x),ci))
  }else
  {
    print(x)
  }
  cat(paste("\nCensoring rate:",round(attributes(x)$censoring.rate,4),"\n"))
}


#' Return the \code{xreg} part of the \code{cenTS} object
#' @param object a \code{cenTS} object.
#' @return the list in \code{xreg}.
#' @seealso \code{\link{cenTS}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  xreg(ts)
#'

xreg <- function(object) UseMethod("xreg")

#' Return the \code{xreg} part of the \code{cenTS} object
#' @param object a \code{cenTS} object.
#' @return the list in \code{xreg}.
#' @seealso \code{\link{cenTS}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  xreg(ts)
xreg.cenTS <- function(object)
{
  if(!is.null(attributes(object)$xreg))
    xts::as.xts(object[,attributes(object)$xreg])
  else
    NULL
}


#' Plot a \code{cenTS} object
#' @param x a \code{cenTS} object.
#' @param type,auto.grid,major.ticks,minor.ticks,major.format,bar.col,candle.col,ann,axes,ylim,main,...
#' standard parameters to control the plot.
#' @seealso \code{\link[xts]{plot.xts}}.
#' @export
#' @examples
#' strDates <- c("2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05")
#' ts <- cenTS(value=c(1,-2,1,NA,0),
#'             order.by=as.Date(strDates,"%Y-%m-%d"),
#'             lcl=c(-3,-2,-1,-1,0),
#'             ucl=c(3,2,1,1,1),
#'             x=c(1,1,1,1,1),
#'             y=c(2,2,2,2,2))
#'  plot(ts)
plot.cenTS <- function(x, type = "l", auto.grid = TRUE, major.ticks = "auto",
    minor.ticks = TRUE, major.format = TRUE, bar.col = "grey",
    candle.col = "white", ann = TRUE, axes = TRUE,ylim=NULL,main=NULL, ...)
{
    ci = attributes(x)$ci
    value.name <- attributes(x)$value.name
    series.title <- series.title <- deparse(substitute(x))
    if(value.name != "value")
      series.title <- value.name
    if(!is.null(main))
      series.title <- main

    ep <- xts::axTicksByTime(x, major.ticks, format.labels = major.format)
    otype <- type

    xy <- grDevices::xy.coords(xts::.index(x), x[,value.name])
    ylim0 <- range(zoo::coredata(x)[,intersect(colnames(zoo::coredata(x)),c(value.name,"lcl","ucl"))],na.rm=TRUE,finite=TRUE)
    if(is.null(ylim))
      ylim <- ylim0
    else
    {
      ylim[1] <- min(ylim0[1],ylim[1])
      ylim[2] <- max(ylim0[2],ylim[2])
    }

    graphics::plot(xy[["x"]], xy[["y"]], type = type, axes = FALSE, ann = FALSE, ylim=ylim,...)

    if(!is.null(x$lcl) && any(is.finite(x$lcl)))
      graphics::lines(xy$x,x$lcl,col="red",lty=3)
    if(!is.null(x$ucl) && any(is.finite(x$ucl)))
      graphics::lines(xy$x,x$ucl,col="red",lty=5)
    
    upTriangle=2
    downTriangle=4
    idxL = which(ci=="L")
    if(length(idxL)>0)
      graphics::segments(xy$x[idxL],x$lcl[idxL],xy$x[idxL],ylim[1],lty=3)
    idxR = which(ci=="R")
    if(length(idxR)>0)
      graphics::segments(xy$x[idxR],x$ucl[idxR],xy$x[idxR],ylim[2],lty=5)
    idxI = which(ci=="I")
    if(length(idxI)>0)
      graphics::segments(xy$x[idxI],x$lcl[idxI],xy$x[idxI], x$ucl[idxI],lty=5)

    if (auto.grid) {
        graphics::abline(v = xy$x[ep], col = "grey", lty = 4)
        graphics::grid(NA, NULL)
    }
    dots <- list(...)
    if (axes) {
        if (minor.ticks)
            graphics::axis(1, at = xy$x, labels = FALSE, col = "#BBBBBB",
                ...)
        graphics::axis(1, at = xy$x[ep], labels = names(ep), las = 1,
            lwd = 1, mgp = c(3, 2, 0), ...)
        graphics::axis(2, ...)
    }
    graphics::box()
    if (!"main" %in% names(dots))
        graphics::title(main = series.title)
    do.call("title", list(...))
}
