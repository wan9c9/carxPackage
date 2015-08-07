sfo8903 <- read.csv('/home/chao/workspace/carxPackage/data/sfo8903.csv',header=TRUE)
print(sfo8903)
sfo8903.dat = sfo8903[,3]
sfo8903.dat[516]=log(120)#NaN
sfo8903.dat[540]=2.55 # mean of adjacent values, previous NaN
sfo8903.dat[694]=log(120)#NaN

so8903.dat = sfo8903.dat * (sfo8903.dat>0) + 0.1 * (sfo8903.dat<=0)
sfo8903.dat.log = log(sfo8903.dat)

#prepare model data
y <- sfo8903.dat.log
nObs <- length(y)
lowerCensorLimit <- rep(-Inf,nObs)
upperCensorLimit <- rep(log(120),nObs)
x <- rep(1,nObs) #intercept
censorIndicator <- rep(0,nObs)
censorIndicator[y>=log(120)]  <- 1
#no left censoring
nAR <- 2
model <- carx(y, x, censorIndicator, lowerCensorLimit, upperCensorLimit, nAR, verbose=T)
plot(model)

