library(carx)
context("test cenTS")

load_all()
strDates <- c("2000-01-01",
              "2000-01-02",
              "2000-01-03",
              "2000-01-04",
              "2000-01-05")

#debug(cenTS)

ts <- cenTS(value=c(1,-2,1,NA,0),order.by=as.Date(strDates,"%Y-%m-%d"),
            lcl=c(-3,-2,-1,-1,0),ucl=c(3,2,1,1,1),
            x=c(1,1,1,1,1),y=c(2,2,2,2,2))
print(ts)
print(xreg(ts))
plot(ts)


ts2 <- cenTS(value=c(1,-2,1,NA,0),order.by=as.Date(strDates,"%Y-%m-%d"),
            lcl=c(-3,-2,-1,-1,0),ucl=c(3,2,1,1,1),
            x=data.frame( x=c(1,1,1,1,1),y=c(2,2,2,2,2))
            )
print(ts2)
print(xreg(ts2))
plot(ts2)

ts3 <- cenTS(value=c(1,-2,1,NA,0),order.by=as.Date(strDates,"%Y-%m-%d"),
            lcl=c(-3,-2,-1,-1,0),ucl=c(3,2,1,1,1),
            x=data.frame( x=c(1,1,1,1,1),y=c(2,2,2,2,2)),
            z=rep(0,5)
            )
print(ts3)
print(xreg(ts3))
plot(ts3)

