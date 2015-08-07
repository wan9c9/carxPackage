

load_all()

#test select the AR order


singleTestSelectAROrder <- function(iter)
{
  seed <- 1375911
  cts <- carx.sim.cenTS(seed=iter*seed)
  m <- carx.select(list(f1=as.formula(value~X1+X2)),max.ar=4, data=cts)
  m$fitted$p
}


nRep <- 100
orders <- numeric(nRep)
orders <- parallel::mclapply(1:nRep, singleTestSelectAROrder,mc.cores=parallel::detectCores()-1)
hist(unlist(orders))
