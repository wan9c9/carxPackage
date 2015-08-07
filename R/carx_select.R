#' Select the best CARX model from a variety of specifications with maximal AR order by AIC
#'
#' Select the best CARX model from a variety of specifications with maximal AR order by AIC
#' @param formulas a list of model formulas
#' @param data the data for CARX
#' @param max.ar maximal AR order
#' @param detect.outlier logical to determine if outlier detection is performed before the AIC for a particular model formula is computed.
#' @return a list consisting of:
#' \itemize{
#' \item{\code{fitted}}{ the fitted CARX object of the model with the smallest AIC}
#' \item{\code{detect.outlier}}{the argument \code{detect.outlier}}
#' \item{\code{aicMat}}{ the matrix of AIC where rows correspond to the model formulas and columns correspond to AR order}
#' }
carx.select <- function(formulas, max.ar, data=list(), detect.outlier=F)
{
  nModel <- length(formulas)
  aics <- matrix(nrow=nModel,ncol=max.ar)
  rownames(aics) <- names(formulas)
  saic0 <- NULL
  m0 <- NULL
 
  i  <- 0
  for(f  in formulas)
  { 
    i  <- i + 1
    for( p in 1:max.ar)
    {
      #print(f)
      #print(p)
      tmp <- carx(f,data=data, p=p, CI.compute=FALSE)
      if(detect.outlier)
        tmp <- outlier(tmp)$updatedModel
      a <- AIC(tmp)
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
  list(fitted = m0, detect.outlier=detect.outlier, aicMat = aics)
}
