dat = carxSimCenTS(nObs=2000,lcl=-0.5,ucl=0.5,intervalCensoring=TRUE)
mdl <- carx(y=dat$y, x=dat[,c("X1","X2")], ci=dat$ci, lcl=dat$lcl, ucl=dat$ucl, p=2)
