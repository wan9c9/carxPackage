
latexFigureStr <- function(figureFile, option, caption, label){
	s <- "\\begin{figure}[!htbp]\n"
	s <- paste(s , "\\includegraphics[",option, "]{",figureFile,"}\n")
	s <- paste(s, "\\caption{",caption,"\\label{",label,"}}\n")
	s <- paste(s,"\\end{figure}\n")
	s
}

parseData <- function(raw){
	raw <- as.matrix(raw)
	d <- dim(raw)
	nObs <<- d[1]
	y <<- raw[,1]
	censorIndicator <<- (raw[,2] != 0)
	censorLimit <<- raw[,3]
	externalVariable <<- raw[,4:d[2]]
	externalVariable <<- cbind(externalVariable, externalVariable[,3]^2)
	#externalVariable <<- cbind(externalVariable, sign(externalVariable[,3])*(sqrt(abs(externalVariable[,3]))))
	#externalVariable[,3] <<- log( 1 + exp( externalVariable[,3]) )
	nEV <<- dim(externalVariable)[2]
}

loadData <- function(fileDest){
	raw <- as.matrix(read.csv(fileDest, header=TRUE))
	d <- dim(raw)
	nObs <<- d[1]
	y <<- raw[,2]
	#y <<- log( 1 + exp(raw[,1]) )
	censorIndicator <<- (raw[,2] != 0)
	censorLimit <<- raw[,3]
	externalVariable <<- raw[,4:d[2]]
	externalVariable <<- cbind(externalVariable, externalVariable[,3]^2)
	#externalVariable <<- cbind(externalVariable, sign(externalVariable[,3])*(sqrt(abs(externalVariable[,3]))))
	#externalVariable[,3] <<- log( 1 + exp( externalVariable[,3]) )
	nEV <<- dim(externalVariable)[2]
}


generateLatexFile <- function(fileDes){
	#for generating latex figure file
	ids <- unlist(strsplit(fileDes,'_'))
	storetID <- ids[2]
	stationName <- unlist( strsplit(ids[3],'[.]'))[1]

	tmpFile <- paste0("../report4v5/",storetID,".tex")
	system(paste0("cp ../report4v5/genericLatexFigure.tex ",tmpFile))
	system(paste0("sed -i -e 's/wideSiteNumber/",storetID,"/g' ",tmpFile))
	system(paste0("sed -i -e 's/wideSiteName/",stationName,"/g' ",tmpFile))
}


realDataAnalysis <- function(fileDes, arOrder,models = NULL, outputPrefix="",modelNames=NULL,outPutFile=NULL,...){ 
	rawData <- read.csv(fileDes, header=TRUE, stringsAsFactors=F,as.is=T)
	rawData$logP <- as.numeric(rawData$logP)
	rawData$censorLimit <- as.numeric(rawData$censorLimit)
	rawData$logQ <- as.numeric(rawData$logQ)
	nObs <- dim(rawData)[1]
	ids <- unlist(strsplit(fileDes,'_'))
	storetID <- ids[2]
	stationName <- unlist( strsplit(ids[3],'[.]'))[1]
	rslt= list




	if(outputPrefix == "") outputPrefix <- storetID 

	cdir <- getwd()

	wd <- paste0('./',storetID)
	dir.create(wd)
	setwd(wd)
	if(!is.null(models)){
		wd <- paste0('./',paste0(c(models)))
		dir.create(wd)
		setwd(wd)
	}
	print(getwd())
	if(is.null(outPutFile))
		htmlFile <- HTMLInitFile(getwd(),filename="index", BackGroundColor="#BBBBEE")
	else
		htmlFile <- outPutFile
	HTML(paste("Analysing ", storetID, stationName))
	#HTMLStart(getwd(),filename="summary")

	# prepare scatter plot
	rawData$season = 1
	rawData$season[rawData$month > 3] = 2
	rawData$season[rawData$month > 6] = 3
	rawData$season[rawData$month > 9] = 4

	#plot the scatter plot marked by season
	scatterBySeason = T
	if(scatterBySeason){
		setEPS()
		postscript(paste0(storetID,"_scatterPlotMonth.eps"))
		plot(rawData$logQ,rawData$logP,xlab="log(Q)",ylab="log(P)",type="n",main="")
		for(i in 1:4)
		{
			lP <- rawData$logP[rawData$season==i]
			lQ <- rawData$logQ[rawData$season==i]

			if (length(lP) != 0)
			{
				m <- lm(lP~lQ)
				lines(lQ,fitted(m),col=i)
				text(lQ,lP,label=i,col=i)
			}
		}
		message("scatterBySeason done!")
		dev.off()
	}

	#return()
	plotScatter = T
	if(plotScatter){
		setEPS()
		postscript(paste0(storetID,"_scatterPlot.eps"))
		plot(rawData$logQ,rawData$logP,xlab="log(Q)",ylab="log(P)",type='n')
		title( paste0('Scatter: ', storetID, '( ', stationName, ') ') )
		for(i in 1:12)
		{
			lP <- rawData$logP[rawData$month==i]
			lQ <- rawData$logQ[rawData$month==i]
			if (length(lP) != 0)
				text(lQ,lP,label=i,col=i)
		}
		dev.off()
	}


	if(is.null(models))
	{
		s1 <- logP ~ tInMonth + logQ
		s2 <- logP ~ tInMonth + logQ:as.factor(season)
		s3 <- logP ~ as.factor(season) + tInMonth + logQ - 1
		s4 <- logP ~ as.factor(season) + tInMonth + logQ:as.factor(season) -1
		s <- c(s1,s2,s3,s4)
	}else
		s <- models
	aic <- NULL
	nAR <- NULL
	spec <- NULL
	HTML("----------- Selecting Model -------------")

	#construct skipIndex so that the criteria based AIC are consistent
	skipIndex <- seq(1,arOrder)
	t <- rawData$tInMonth
	for( i in 2:length(t))
	{
		if( t[i] - t[i-1] > 1)
			skipIndex <- c(skipIndex, seq(i,i+arOrder-1))
	}
	nGap <- sum(diff(t)!=1)
	#HTML(paste("number of gaps ",nGap))

	for( i in 1:arOrder)
	{
		j = 0
		for(m in s)
		{
			j = j+1
			HTML(" *****************")
			HTML(paste(c("spec: ", m)))
			HTML(paste("AR order", i))
			data1 <- rawData
			message(sprintf("AR order: %d, Model: %s",i, format(m)))
			tmp <- carx(m, data=data1, censorIndicator=data1$censorIndicator,censorLimit=data1$censorLimit,nAR=i,getCI=FALSE,skipIndex=skipIndex,...)
			HTML(summary(tmp))
			censorRate <- tmp$censorRate
			fn <- paste0(outputPrefix,"_AR",i,"_M",j)
			plot(tmp,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="log(P)",saveFig=paste0(fn,"_obs1.eps"))
			HTML(paste("AIC:", tmp$aic))
			if(is.null(aic))
			{
				aic <- tmp$aic
				nAR <- i
				spec <- m
				fData <- data1
				fOlv <- NULL
				if(j == 1) 
					pC <- c(i,1,1,1,1)
				else if( j == 2)
					pC <- c(i,1,1,4,1)
				else if( j == 3)
					pC <- c(i,4,1,1,1)
				else
					pC <- c(i,4,1,4,1)

			} else{
				if(tmp$aic < aic)
				{ 
					aic <- tmp$aic
					nAR <- i
					spec <- m
					fData <- data1
					fOlv <- NULL
					if(j == 1) 
						pC <- c(i,1,1,1,1)
					else if( j == 2)
						pC <- c(i,1,1,4,1)
					else if( j == 3)
						pC <- c(i,4,1,1,1)
					else
						pC <- c(i,4,1,4,1)
				}
			}
			#detecting outlier
			ot <- 1 #initialize ot flag
			nOL <- 0
			while(ot > 0){
				ot <- outlierDetection(tmp)
				print(ot)
				if(ot != -1){ 
					if(nOL == 0)
						olv <- ot
					else{
						if(is.vector(olv) && ot %in% olv)
							break
				olv <- c(olv,ot)
					}
					nOL <- nOL +1
					oi <- numeric(tmp$nObs)
					oi[ot] <- 1
					newVar <-paste0("OutlierIndicator",ot)
					data1[,newVar] <- oi
					HTML(paste(newVar, "is added."))
					m <- update(m, as.formula(paste("~.+",newVar)))
					message(sprintf("AR order: %d, Model: %s",i, format(m)))
					tmp <- carx(m, data=data1, censorIndicator=data1$censorIndicator,censorLimit=data1$censorLimit,nAR=i,getCI=FALSE,skipIndex=skipIndex,...)
					HTML(summary(tmp))
					fn <- paste0(fn,"_OT_",ot)
					plot(tmp,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="log(P)",saveFig=paste0(fn,"_obs1.eps"),outliers=ot)
					#readline("Press <return to continue")
					#plot(tmp,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="log(P)") #,saveFig=paste0("../../report4_v4/",outputPrefix,"_obs1.eps"))
					HTML(paste(c("spec: ", m)))
					HTML(paste("AIC:", tmp$aic))
					if(is.null(aic))
					{
						aic <- tmp$aic
						nAR <- i
						spec <- m
						fData <- data1
						fOlv <- olv
						if(j == 1) 
							pC <- c(i,1,1,1,1)
						else if( j == 2)
							pC <- c(i,1,1,4,1)
						else if( j == 3)
							pC <- c(i,4,1,1,1)
						else
							pC <- c(i,4,1,4,1)

					} else{
						if(tmp$aic < aic)
						{ 
							aic <- tmp$aic
							nAR <- i
							spec <- m
							fData <- data1
							fOlv <- olv
							if(j == 1) 
								pC <- c(i,1,1,1,1)
							else if( j == 2)
								pC <- c(i,1,1,4,1)
							else if( j == 3)
								pC <- c(i,4,1,1,1)
							else
								pC <- c(i,4,1,4,1)
						}
					}
				}
			}
		} 
	}

	#re-construct skipIndex
	skipIndex <- seq(1,nAR)
	t <- rawData$tInMonth
	for( i in 2:length(t))
	{
		if( t[i] - t[i-1] > 1)
			skipIndex <- c(skipIndex, seq(i,i+nAR-1))
	}


	tmpModel <- carx( spec, data=fData, censorIndicator=fData$censorIndicator,censorLimit=fData$censorLimit,nAR=nAR,getCI=F,skipIndex=skipIndex,...)

	plotPartialResidual = T
	if(plotPartialResidual){
		co <- coef(tmpModel)
		resi <- residuals(tmpModel)
		#pr <- rawData$logP - co["tInMonth"]*rawData$tInMonth #get partial residuals by excluding the time trend
		#rawData$pr <- pr #assign to rawData dataframe

		#inter <- co["(Intercept)"]

		setEPS()
		postscript(paste0("../../report4v5/",storetID,"_partialResidualPlot.eps"))
		plot(rawData$logQ,rawData$logP,xlab="log(Q)",ylab="Partial Residual",type='n')

		for(i in 1:4)
		{
			lQ <- rawData$logQ[rawData$season==i]
			res <- resi[rawData$season==i]

			if (length(lQ) != 0)
			{
				if(pC[2] == 4 )
					interS <- co[sprintf("as.factor(season)%d",i)]
				else 
					interS <- co["(Intercept)"]
				if(pC[4] == 4)
				{
					multi <- co[sprintf("logQ:as.factor(season)%d",i)]
					if(is.na(multi))
						multi <- co[sprintf("as.factor(season)%d:logQ",i)]
				}

				else 
					multi <- co["logQ"]
					f <- interS + multi*lQ
					lines(lQ,f,col=i)
					text(lQ,f+res,label=i,col=i)
			}
		}
		dev.off()
	}
	#return()


	HTML("----------- Selected Model -------------")
	HTML(paste(c("storetID: ", storetID)))
	HTML(paste(c("stationName: ", stationName)))
	HTML(paste0(c("date range: ", rawData$year[1],"-",rawData$month[1], " to  ", rawData$year[nObs],"-",rawData$month[nObs])))
	HTML(paste("number of obs: ",nObs))
	HTML(paste("censor rate: ",censorRate))
	HTML(paste("number of gaps: ",nGap))
	HTML(paste( c("\nmodel spec: ",spec)))
	if(!is.null(fOlv))
		HTML(paste(c("Outliers (tInMonth): ", rawData$tInMonth[fOlv])))
	HTML(paste("model AR order:",nAR))
	HTML(paste("model npar:",tmpModel$npar))
	prmtrStr <- carx.prmtrStr(tmpModel,pC)
	HTML(prmtrStr)


	HTML(summary(tmpModel))
	#message(summary(tmpModel)$coefficients)
	#message(prmtrStr)
	summaryStr <- paste(sprintf("\\pbox{6cm}{%s } & %s \\\\ \\hline",stationName,prmtrStr))

	message("printing plots for fitted model...")
	#plot(tmpModel,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="log(P)",saveFig=paste0(outputPrefix,"_obs1.eps"))
	plot(tmpModel,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="log(P)",saveFig=paste0("../../report4v5/",outputPrefix,"_obs1.eps"), outliers=fOlv)
	#plot(tmpModel,transformFun=exp,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="P",saveFig=paste0(outputPrefix,"_obs2.eps"))
	plot(tmpModel,transformFun=exp,xAxisVar=rawData$tInMonth,xlab="Time (in month)", ylab="P",saveFig=paste0("../../report4v5/",outputPrefix,"_obs2.eps"),outliers=fOlv)
	plot.residuals(tmpModel,"fitted",paste0("../../report4v5/",storetID, "_resVsFitted.eps"),ylab="Standardized Residual",type="p")
	#plot.residuals(tmpModel,"fitted",paste0(storetID, "_resVsFitted.eps"),ylab="Residual",type="p")
	plot.residuals(tmpModel,rawData$tInMonth,paste0("../../report4v5/",storetID,"_tsRes.eps"),xlab="Time (in month) ",ylab="Standardized Residual")
	#plot.residuals(tmpModel,rawData$tInMonth,paste0(storetID,"_tsRes.eps"),xlab="Time (in month) ",ylab="Residual")
	plot.residuals(tmpModel,rawData$logQ,paste0("../../report4v5/",storetID,"_resVsQ.eps"),xlab="log(Q)",ylab="Standardized Residual",type="p")
	#plot.residuals(tmpModel,rawData$logQ,paste0(storetID,"_resVsQ.eps"),xlab="log(Q)",ylab="Residual",type="p")
	plotResAcf(tmpModel, paste0("../../report4v5/",outputPrefix,"_resAcf.eps"))
	#plotResAcf(tmpModel, paste0(outputPrefix,"_resAcf.eps"))
	ret <- list(storetID=storetID,
		    stationName=stationName,
		    nAR = nAR,
		    model=spec,
		    prmtrStr=prmtrStr,
		    summaryStr=summaryStr,
		    est = tmpModel,

		    summary=summary(tmpModel)
		    )
	save(tmpModel,file="model.RData")
	setwd(cdir)
	ret
}

outlierDetection <- function(model){
	message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/model$nObs
	eps <- rnorm(nSample,0,model$sigma)
	trend <- model$x%*%model$prmtrEV
	covEta <- computeCovAR(model$prmtrAR, model$sigma)
	nObs <- model$nObs
	nAR <- model$nAR
	prmtrAR <- model$prmtrAR
	skipIndex <- model$skipIndex
	y <- model$y
	censorIndicator <- model$censorIndicator
	censorLimit <- model$censorLimit
	y[censorIndicator] <- censorLimit[censorIndicator]


	pValues <- numeric(nObs)
	pValues[skipIndex] <- 1

	for(idx in seq(1,nObs)[-skipIndex])
	{
		#message(sprintf("checking %i",idx))
		wkm <- y[(idx-1):(idx-nAR)]
		tmpCensorIndicator <- censorIndicator[(idx-1):(idx-nAR)]
		nCensored <- sum(tmpCensorIndicator)
		if(nCensored)   #at least one is censored
		{
			if( nCensored < nAR ) #not all are censored
			{ 
				conditionalIndex <- which(!tmpCensorIndicator)
				tmpY <- y[(idx-1):(idx-nAR)][conditionalIndex] 
				tmpM <- trend[(idx-1):(idx-nAR)]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta[-1,-1])
				tmpMean <- cdist$'mean' 
				tmpVar <- cdist$'var' 
			}else{ 
				tmpMean <- trend[(idx-1):(idx-nAR)]
				tmpVar <- covEta[-1,-1]
			}
			tmpCensorLimit <- censorLimit[(idx-1):(idx-nAR)][tmpCensorIndicator]
			#print(tmpMean)
			#print(tmpVar)
			smpl <- rtmvnorm(nSample,tmpMean,tmpVar,upper=tmpCensorLimit,algorithm="gibbs")
			smpl <- as.matrix(smpl)
			#print(smpl)
			ySmpl <- numeric(nSample)
			for(i in 1:nSample){
				wkm[tmpCensorIndicator] <- smpl[i,]
				ySmpl[i] <- trend[idx] + (wkm - trend[(idx-1):(idx-nAR)])%*%prmtrAR + eps[i]
			}
			pU <- sum(ySmpl > y[idx])/nSample
			pL <- sum(ySmpl < y[idx])/nSample
		}
		else{
			r <- y[idx]-trend[idx] - (wkm-trend[(idx-1):(idx-nAR)])%*%prmtrAR
			r <- r/model$sigma
			pU <- pnorm(r,lower.tail=FALSE)
			pL <- pnorm(r,lower.tail=TRUE)
		}
		pValues[idx] <- min(pU,pL)
	}
	minP <- min(pValues)
	if( minP <= threshold ){
		i <- which(pValues == minP)
		i
	} else
		-1
}

