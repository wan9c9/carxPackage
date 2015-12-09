
#' @return 0: null hypothesis that no trend is present in the model is not rejected, 1: null is rejected.
testTrendCarx <- function(object,level=0.05,nb=19,k=5,verbose=TRUE)
{
  require(mgcv)
  #browser()
  ptm0 <- proc.time()
  computeTestStat <- function(res,k)
  {
    #ss <- smooth.spline(res)
    ss <- mgcv::gam(res~s(seq(1,length(res)),bs='cr',k=k))
    f <- fitted(ss)
    c <- cor(res,f)
    c^2
  }
  res <- residuals(object)
  ts0 <- computeTestStat(res,k)
  bts <- numeric(nb)
  for(iRep in 1:nb)
  {
    simData <- carxSim(nObs=object$nObs,
                       prmtrAR=object$prmtrAR,
                       prmtrX=object$prmtrX,
                       sigma=object$sigma,
                       lcl=object$lcl,
                       ucl=object$ucl,
                       x=object$x)
    bObj <- carx(y=simData$y,
                 x=object$x,
                 ci=simData$ci,
                 lcl=simData$lcl,
                 ucl=simData$ucl,
                 p = object$p,
                 prmtrX=object$prmtrX,
                 prmtrAR=object$prmtrAR,
                 sigma=object$sigma
                 )
    resp <- residuals(bObj)
    bts[iRep] <- computeTestStat(resp,k)
  }
  #val <- ifelse((sum(bts>ts0)+1)/(nb+1) <= level, 1,0)
  val <- (sum(bts>ts0)+1)/(nb+1)
	ptm1 <- proc.time() - ptm0
  if(verbose){
    hist(bts,breaks=as.integer(nb/10),main=sprintf("nb:%4d,p-value: %5.2f ",nb, val))
		abline(v=ts0,col="red")
		message(sprintf("Test done in %5.1f secs, with %d bootstraps. Result: %5.3f",ptm1[3],nb,val))
	}
  val
}
