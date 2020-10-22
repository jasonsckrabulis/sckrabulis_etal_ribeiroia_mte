# Sensitivity analysis and To optimization figure

#_
# Figure S1: Sensitivity analysis

# Constants
K<-273.15
k<-8.62*10^-5

# Plotting functions

# Simulated resistance equation
resistPlot <- function(x, rTo, Ear, Thr){
	Tor <- 19+K
	Ehr <- 3.25
	y <- rTo*exp(-(Ear/k*(1/(x+K)-1/(Tor))))*(1+exp(Ehr/k*(1/Thr-1/(x+K))))^-1
	ifelse(y>1,1,y)
}

# Infectivity curve, using estimates from the best fit model in 'Cercaria swimming speed.txt'
infectPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- (-1.382)
	infect <-iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1+C
	y <- ifelse(infect>=1, 1, infect)
	ifelse(y<0, 0, y)
}

# Interaction model (infectivity-resistance)
fullPlot <- function(x, rTo, Ear, Thr){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <-(-1.382)
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1+C
	thresh <- ifelse(infect>=1, 1, infect)
	Ehr <- 3.25
	Tor <- 19+K
	resist <- rTo*exp(-(Ear/k*(1/(x+K)-1/(Tor))))*(1+exp(Ehr/k*(1/Thr-1/(x+K))))^-1
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Plot parameters
par(mfrow=c(3,2),omi=c(.4,.4,0,.4))

#A: Linear change in rTo for resistance

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot curves at incremental inreases to rTo
	curve(resistPlot(x,rTo=.1,Ear=.65,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, xlab="", ylab="", xaxt="n", lty=6, col="grey65")
	curve(resistPlot(x,rTo=.15,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=5, col="grey55")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=4, col="grey45")
	curve(resistPlot(x,rTo=.25,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=3, col="grey35")
	curve(resistPlot(x,rTo=.3,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=2, col="grey25")
	curve(resistPlot(x,rTo=.35,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=1)
	axis(side=1, labels=FALSE)
	legend(32.5, 1.05, c('.1','.15','.2','.25','.3','.35'), lty=c(6,5,4,3,2,1), lwd=c(2,2,2,2,2,2),	col=c('grey65','grey55','grey45','grey35','grey25','black'), bty='n')
	mtext(LETTERS[1], side=3, adj=.05, line=-2, cex=1)

#B: Inserting resistances from Panel A into interaction model

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot constant infectivty curve
	curve(infectPlot, from=10, to=40, ylim=c(0,1), lwd=3, xlab="", ylab="", xaxt="n", yaxt="n")
	
	# Plot resultant interaction curves (infect-resist)
	curve(fullPlot(x,rTo=.1,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=6, col="grey65")
	curve(fullPlot(x,rTo=.15,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=5, col="grey55")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=4, col="grey45")
	curve(fullPlot(x,rTo=.25,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=3, col="grey35")
	curve(fullPlot(x,rTo=.3,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=2, col="grey25")
	curve(fullPlot(x,rTo=.35,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=1)
	axis(side=1, labels=FALSE)
	axis(side=4)
	mtext(LETTERS[2], side=3, adj=.05, line=-2, cex=1)

#C: Linear change in Thr for resistance
	
	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot curves and incremental increases to Thrr
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=299.15), from=10, to=40, ylim=c(0,1), lwd=2, xlab="", ylab="", xaxt="n", lty=6, col="grey65")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=301.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=5, col="grey55")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=4, col="grey45")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=305.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=3, col="grey35")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=307.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=2, col="grey25")
	curve(resistPlot(x,rTo=.2,Ear=.65,Thr=309.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=1)
	axis(side=1, label=FALSE)
	mtext(LETTERS[3], side=3, adj=.05, line=-2, cex=1)
	legend(32.5, 1.05, c('26','28','30','32','34','36'), lty=c(6,5,4,3,2,1), lwd=c(2,2,2,2,2,2), col=c('grey65','grey55','grey45','grey35','grey25','black'), bty='n')

#D: Inserting resistances from Panel C into interaction model

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot constant infectivity curve
	curve(infectPlot, from=10, to=40, ylim=c(0,1), lwd=3, xlab="", ylab="", yaxt="n", xaxt="n")
	
	# Plot resultant interaction curves (infect-resist)
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=299.15), from=10, to=40, lwd=2, add=TRUE, lty=6, col="grey65")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=301.15), from=10, to=40, lwd=2, add=TRUE, lty=5, col="grey55")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=303.15), from=10, to=40, lwd=2, add=TRUE, lty=4, col="grey45")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=305.15), from=10, to=40, lwd=2, add=TRUE, lty=3, col="grey35")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=307.15), from=10, to=40, lwd=2, add=TRUE, lty=2, col="grey25")
	curve(fullPlot(x,rTo=.2,Ear=.65,Thr=309.15), from=10, to=40, lwd=2, add=TRUE, lty=1)
	axis(side=1, label=FALSE)
	axis(side=4)
	mtext(LETTERS[4], side=3, adj=.05, line=-2, cex=1)

#E: Linear change in Ear for resistance

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot curves with incremental increases to Earr
	curve(resistPlot(x,rTo=.2,Ear=.2,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, xlab="", ylab="", lty=6, col="grey65")
	curve(resistPlot(x,rTo=.2,Ear=.4,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=5, col="grey55")
	curve(resistPlot(x,rTo=.2,Ear=.6,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=4, col="grey45")
	curve(resistPlot(x,rTo=.2,Ear=.8,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=3, col="grey35")
	curve(resistPlot(x,rTo=.2,Ear=1,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=2, col="grey25")
	curve(resistPlot(x,rTo=.2,Ear=1.2,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=1)
	axis(side=1, labels=FALSE)
	legend(32.5, 1.05, c('.2','.3','.4','.5','.6','.7'), lty=c(6,5,4,3,2,1), lwd=c(2,2,2,2,2,2), col=c('grey65','grey55','grey45','grey35','grey25','black'), bty='n')
	mtext(LETTERS[5], side=3, adj=.05, line=-2, cex=1)

#F: Inserting resistances from Panel E into interaction model

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	# Plot constant infectivity curve
	curve(infectPlot, from=10, to=40, ylim=c(0,1), lwd=3, xlab="", ylab="", yaxt="n")
	
	# Plot resultant interaction curves (infect-resist)
	curve(fullPlot(x,rTo=.2,Ear=.2,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=6, col="grey65")
	curve(fullPlot(x,rTo=.2,Ear=.4,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=5, col="grey55")
	curve(fullPlot(x,rTo=.2,Ear=.6,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=4, col="grey45")
	curve(fullPlot(x,rTo=.2,Ear=.8,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=3, col="grey35")
	curve(fullPlot(x,rTo=.2,Ear=1,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=2, col="grey25")
	curve(fullPlot(x,rTo=.2,Ear=1.2,Thr=303.15), from=10, to=40, ylim=c(0,1), lwd=2, add=TRUE, lty=1)
	axis(side=4)
	mtext(LETTERS[6], side=3, adj=.05, line=-2, cex=1)

	mtext(expression(paste("Performance temperature (",degree,"C)")), side=1, line=1.5, outer=TRUE)
	mtext("Simulated resistance", side=2, line=1.5, outer=TRUE)
	mtext("Predicted proportion of encysted metacercariae", side=4, line=1.5, outer=TRUE)

#_
# Figure S2: Toc, Tor, and ToR optimization for best fit models

# Use Altman et al. 2016 data on Dryad
cystData <- read.csv(file.choose(),header=TRUE)
cystData <- cystData[order(cystData$PerfTemp),]
cystData$logAdjPropEncysted <- log(cystData$AdjPropEncysted)
cystData$AccCenter <- cystData$AccTemp-mean(cystData$AccTemp)

# Use Altman et al. 2016 data on Dryad
propClearData <- read.csv(file.choose(),header=TRUE)
propClearData$PerfTempK <- propClearData$PerfTemp+K
propClearData$AccTempK <- propClearData$AccTemp+K
propClearData$AccCenter <- propClearData$AccTemp-mean(propClearData$AccTemp)
propClearData <- subset(propClearData, PropCleared1!="NA")
propClearData$logPropClear <- log(propClearData$PropCleared1+1)

# Use 'Uninfected tadpole respiration.csv'
respData <- read.csv(file.choose(),header=TRUE)
respData$PerfTempK <- respData$PerfTemp+K
respData$AccTempK <- respData$AccTemp+K
respData <- subset(respData, AccTemp!="NA"&corO2.Time.Mass!="NA")
respData$logcorO2.Time.Mass <- log(respData$corO2.Time.Mass+1)
respData$AccCenter <- respData$AccTemp-mean(respData$AccTemp)

# Constants
k <- 8.62*10^-5
K <- 273.15

# Respiration best fit model ToR, taken from 'Tadpole respiration.csv'
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EdR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EdR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EdR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), start=c(EdR=3, bRTo=.4, bER=.6, mER=0, qER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
A <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])

# Metacercaria encystment best fit model Tor, taken from 'Metacercaria encystment.txt'
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infect <- function(x, iTo, C){
		Ei <- .40734
		Edi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Edi/k*(1/Thi-1/(x+K)))))^-1+C
		ifelse(y>=1, 1, y)
	}
	resistMarrAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))
	}
	model <- nls(formula=AdjPropEncysted ~ infect(x=PerfTemp, iTo=2.58397, C=(-1.28218))-resistMarrAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo=0, bEr, mEr), start=c(brTo=.4, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))
	vectAIC[j] <- AIC(model)
}
B <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])


# Metacercaria clearance best fit model Toc, taken from 'Metacercaria clearance.txt'
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Edc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Edc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Edc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
C <- data.frame(temp=range/10,AIC=vectAIC[vectAIC!=0])

# Plot parameters
par(mfrow=c(3,1), omi=c(.4,.4,0,0))

#A: Respiration best model (quad ea acclimation effect) to=14.55
	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	plot(A$temp, A$AIC, type="l", lwd=2, xlab="", xaxt="n", ylab="")
	axis(side=1, labels=FALSE)
	mtext(LETTERS[1], side=3, adj=.025, line=-1.75, cex=1)

#B: Encystment best model (line ea acclimation effect) to=19
	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	plot(B$temp, B$AIC, type="l", lwd=2, xlab="", xaxt="n", ylab="")
	axis(side=1, labels=FALSE)
	mtext(LETTERS[2], side=3, adj=.025, line=-1.75, cex=1)

#C: Clearance best model (quad ea acclimation effect) to=26.9
	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	plot(C$temp, C$AIC, type="l", lwd=2, xlab="", ylab="")
	mtext(LETTERS[3], side=3, adj=.025, line=-1.75, cex=1)

	mtext("AIC values", side=2, line=1.5, outer=TRUE)
	mtext(expression(paste("Standardization temperature T"[0]*" (",degree,"C)")), side=1, line=1.5, outer=TRUE)