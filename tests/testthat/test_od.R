set.seed(1)
data_len=c(100,200,300)
ot_mag=c(1,2,4,6,8)


for (i in data_len){
  nullNumber=0
  correctPredict=0
  falsePredict=0
  wrongNumber=0
  for (k in 1:10){
    trueX = c(0.2,0.4)
    nAR = 2
    trueAR = c(-0.28,0.25)
    trueSigma = 0.60
    lcl = -1
    ucl = 1
    for (j in ot_mag){
      dat = simulateCarx(i, trueAR, trueX, trueSigma, lcl, ucl)
      dat$y[i/2]=dat$y[i/2]+j*sd(dat$y)
      model <- carx(dat$y, dat$x,censorIndicator=dat$censorIndicator,
                    lowerCensorLimit = dat$lowerCensorLimit,
                    upperCensorLimit = dat$upperCensorLimit,
                    nAR = nAR,
                    getCI = FALSE)
      tmp=od(model)
      if (length(tmp$outlierIndices==1)){
        if(tmp$outlierIndices==i/2){
          correctPredict=correctPredict+1
        }
        else{
          falsePredict=falsePredict+1
        }
      }
      else if(length(tmp$outlierIndices>1)){
        wrongNumber=wrongNumber+1
      }
      else{
        nullNumber=nullNumber+1
      }
    } 
  }
  report=NULL
  report=list("nullNumber"=nullNumber,"correctPredict"=correctPredict,
              "falsePredict"=falsePredict,"wrongNumber"=wrongNumber)
  print(report)
}









