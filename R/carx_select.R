#' Select the best CARX model from a variety of formulas with maximal AR order by AIC
#'
#' Select the best CARX model from a formula or a list of formulas with maximal AR order by AIC.
#' @param formulas a list of model formulas, or a single formula.
#' @param data the data for CARX.
#' @param max.ar maximal AR order.
#' @param detect.outlier logical to determine if outlier detection is performed before the AIC for a particular model formula is computed, default = \code{FALSE}.
#' @param ... other arguments to be supplied, if not null, it will be call with the selected model and data. Examples include {CI.compute=TRUE}, which will cause the function to re-estimate the model with confidence intervals computed, as in the selection part, the confidence interval is not computed.
#' @export
#' @return a list consisting of:
#' \itemize{
#' \item{\code{fitted}}{ the fitted CARX object of the model with the smallest AIC.}
#' \item{\code{aicMat}}{ the matrix of AIC where rows correspond to the model formulas and columns correspond to AR order.}
#' }
#' @examples
#' dataSim <- carxSimCenTS(nObs=100)
#' fmls <- list(M1=y~X1,M2=y~X1+X2,M3=y~X1+X2-1)
#' \dontrun{carxSelect(y~X1,max.ar=3,data=dataSim)}
#' \dontrun{carxSelect(formulas=fmls,max.ar=3,data=dataSim)}
#'
carxSelect <- function(formulas, max.ar, data=list(), detect.outlier=F,verbose=FALSE,...)
{
  if(typeof(formulas) == 'language')
    formulas <- list(M1=formulas)
  if(is.null(formulas))
    names(formulas) <- paste0("M",1:length(formulas))

  nModel <- length(formulas)
  aics <- matrix(nrow=nModel,ncol=max.ar)
  rownames(aics) <- names(formulas)
  colnames(aics) <- paste0("AR",seq(1,max.ar))
  saic0 <- NULL
  m0 <- NULL

  i  <- 0
  for(f in formulas)
  {
    i <- i + 1
    for( p in 1:max.ar)
    {
      #message(paste("Trying model fomula:",deparse(f), ", with AR order", p))
      tmp <- carx(f,data=data, p=p, CI.compute=FALSE,...)
      if(detect.outlier)
        tmp <- outlier(tmp)
      a <- AIC(tmp)
      if(verbose) message(paste0("Model formula:", deparse(formula(tmp)), ", AR order:",tmp$p, ", AIC: ", round(a,digits=4)))
      aics[i,p] <- a
      if(is.null(saic0))
      {
        saic0 <- a
        m0 <- tmp
      }
      else
      {
        if(a < saic0)
        {
          saic0 <- a
          m0 <- tmp
        }
      }
    }
  }
  if(length(list(...))>0)
  {
    m1 <- carx(formula(m0), data=m0$data,p=m0$p,...)
    #outlier indices will be destroyed by the above command)
    m1$outlier.indices <- m0$outlier.indices
    m1$outlier.prefix <- m0$outlier.prefix
    m0 <- m1
  }
  list(fitted = m0, aicMat = aics)
}
