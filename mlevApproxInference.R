#######################################################################################
## Name:         mlevApproxInference.R
## Article :     Streamlined mean field variational Bayes for longitudinal and multilevel
##               data analysis
## Authors :     Cathy Y. Y. Lee and Matt P. Wand
## Last updated: 29 JULY 2015
## Purpose:      As described in http://matt-wand.utsacademics.info/LeeWand.pdf
## R Version:    3.2.0 (2015-04-16)                                                             
## Input data files:  ---                                                      
## Output data files: 2LevModAccRes_Gaussian.txt;  2LevModCovRes_Gaussian.txt
##                    2LevModTimeRes_Gaussian.txt; 2LevModMCMC_Gaussian.pdf
##                    2LevModlogML_Gaussian.pdf;   2LevModfits_Gaussian.pdf
##                    accuBoxplots_Gaussian.pdf
# Required R packages: MASS, magic, Matrix, rstan, gdata, lattice, splines, KernSmooth
########################################################################################

## Clear R memory:

rm(list=ls()) 

## Set working and output directories:

cat("Enter your own working directory and output path.\n")

setwd ("/Users/cathylee/Documents/UTS/PhD/Thesis/Publications/mlevpap/mlevpap.revision/Programs")

output.path <- paste(getwd(),"/Results/",sep="")

## Load required R packages:

libname <- c("MASS","magic","Matrix","rstan","gdata","lattice","splines","KernSmooth")
lapply(libname, require, character.only=T)

## Load user-written functions:

source(paste(getwd(),"/Code/Functions/ZOSull.r",sep=""))
source(paste(getwd(),"/Code/Functions/mlevMCMC.r",sep=""))
source(paste(getwd(),"/Code/Functions/mlevMFVB.r",sep=""))
source(paste(getwd(),"/Code/Functions/summMCMC.r",sep=""))
source(paste(getwd(),"/Code/Functions/accVarApp.r",sep=""))

wait <- function()
{
   cat("Hit return to continue\n")
   ans <- readline()
   invisible()
}

########################################################
###           Require user specification             ###
########################################################

# Set start and end points of the simulation study:

isimSTT <- 1 # Start number
isimEND <- 1 # End number 
responseTypeVal <- "Gaussian" # Gaussian or Bernoulli response

## Set flags for type of study (TRUE/FALSE):

doCoverSim <- T               # Coverage results  
doAccurSim <- T               # Accuracy results
doRecordTime <- T             # Computation time results
doStreamlinedVal <- T         # Naive or Streamlined approach          

## Set flags for simulation (TRUE/FALSE):

doMCMC <- T            # Markov chain Monte Carlo (MCMC) sampling via Stan
doMFVB <- T            # Mean field variational Bayes (MFVB) approximation                  
plotFunctionFits <- T  # Plot penalised spline function fits
plotLowerBound <- T    # Plot lower bound on the marginal log-likelihood
doMFVBvsMCMC <- T      # Comparison between MFVB and MCMC 
maxIterVal <- 200      # Maximum iteration for variational algorithms

## Set flags for Rstan (TRUE/FALSE):

useSavedFit <- T       # Use saved data for MCMC sampling after compilation once
compileCode <- !useSavedFit

## Set flags for creating pdf plots (TRUE/FALSE):

createPDF <- T         # Create pdf outputs?
numPlots <- isimEND    # Number of pdf outputs
userwait <- !createPDF # Allow pausing in R scripts?
      
## Set up plotting parameters:

ng <- 101     # Grid size
lwdVal <- 2   # Line width
cexVal <- 1.5 # Font size

## Set MCMC parameters:

nBurnin <- 100 # Number of burn-ins
nIter <- 100   # Number of iterations
nThin <- 1     # Number of thinnings

########################################################
###         End of user specification                ### 
########################################################

MCMC.plots <- vector(numPlots, mode='list')
accu.plots <- vector(numPlots, mode='list')
fits.plots <- vector(numPlots, mode='list')
logML.plots <- vector(numPlots, mode='list')

## Set up file for coverage results:

if (doCoverSim)
{
   cat("ARE YOU SURE THAT YOU WANT TO OVERWRITE 2LevModCovRes.txt (y/n)????\n")
   ans <- readline()
   filenameCov <- paste(output.path,"2LevModCovRes_",responseTypeVal,".txt",sep="")
   if (ans=="y")
   {
      if (responseTypeVal == "Gaussian") 
         titleLine <- c("isim","Quin1","Quin2","Quin3","Quin4","betax",
                        "SigmaR11","SigmaR22","SigmaR12","sigsqeps")
      if (responseTypeVal == "Bernoulli") 
         titleLine <- c("isim","Quin1","Quin2","Quin3","Quin4",
                        "betax","SigmaR11","SigmaR22","SigmaR12")
      write(titleLine,filenameCov,n=length(titleLine),append=FALSE)
   }
}

## Set up file for accuracy results:

if (doAccurSim)
{
   cat("ARE YOU SURE THAT YOU WANT TO OVERWRITE 2LevModAccRes.txt (y/n)????\n")
   ans <- readline()
   filenameAcc <- paste(output.path,"2LevModAccRes_",responseTypeVal,".txt",sep="")
   if (ans=="y")
   {
      if (responseTypeVal == "Gaussian")
         titleLine <- c("isim","Quin1","Quin2","Quin3","Quin4","betax","SigmaR11",
                        "SigmaR22","SigmaR12","sigsqeps")
      if (responseTypeVal == "Bernoulli")
         titleLine <- c("isim","Quin1","Quin2","Quin3","Quin4","betax","SigmaR11",
                        "SigmaR22","SigmaR12")
      write(titleLine,filenameAcc,n=length(titleLine),append=FALSE)
   }
}

## Set up file for MFVB and MCMC computation times:

if (doRecordTime)
{
   cat("ARE YOU SURE THAT YOU WANT TO OVERWRITE 2LevModTimeRes.txt (y/n)????\n")
   ans <- readline()
   filenameTime <- paste(output.path,"2LevModTimeRes_",responseTypeVal,".txt",sep="")
   if (ans=="y")
   {   
      titleLine <- c("isim","TimeMCMC","TimeMFVB")
      write(titleLine,filenameTime,n=length(titleLine),append=FALSE)
   }
   timeMCMC <- rep(NA,isimEND-isimSTT+1)
   timeMFVB <- rep(NA,isimEND-isimSTT+1)
}

## Set hyperparameters:

nuVal <- 2
sigsq.beta <- 1000
A.eps <- 1e5; A.u <- 1e5 
A.R <- 1e5;   sigsq.beta <- 1e5

## Set true values of parameters:

nuVal <- 2
beta0True <- 0.58
beta1True <- 1.89
SigmaL2True <- matrix(c(2.58,0.22,0.22,1.73),2,2)
sigsqEpsTrue <- 0.1

## Set user-defined functions:

fTrue <- function(x)
{
   return(1 - 2.6*pnorm(x,0.15,0.1) - (2.3*x - 0.07*x^2)
          + 0.5*(1-dnorm(x,0.8,0.07)))
}

inv.logit <- function(x) return(exp(x)/(exp(x)+1))

## Begin simulation studies, see Section 6:

for (isim in isimSTT:isimEND) 
{
   set.seed(isim)

   ## Set the number of groups m and within-group sample sizes n_i:
   
   m <- 100
   nVec <- NULL
   for (i in 1:m)
   {
      if (responseTypeVal == "Gaussian")
         nVec[i] <- sample(10:20,1,replace=TRUE)
      if (responseTypeVal == "Bernoulli")
         nVec[i] <- sample(90:110,1,replace=TRUE)
   }
   numObs <- sum(nVec) # Total number of observations

   ## Generate simulated data according to settings described in Section 6:

   x1 <- NULL; x2 <- NULL; y <- NULL
   idnum <- NULL; currStt <- 1
   reBlockInds <- vector("list", length=m)  
   uR <- mvrnorm(m,rep(0,2),SigmaL2True)
   for (i in 1:m)
   {  
      idnum <- c(idnum,rep(i,nVec[i]))
      
      x1Curr <- runif(nVec[i]);  x1 <- c(x1,x1Curr)
      x2Curr <- runif(nVec[i]); x2 <- c(x2,x2Curr)

      etaCurr <- ((beta0True + uR[i,1]) + (beta1True + uR[i,2])*x1Curr  
                   + fTrue(x2Curr))

      if (responseTypeVal=="Gaussian")
         yCurr <- rnorm(nVec[i],etaCurr,sqrt(sigsqEpsTrue))

      if (responseTypeVal=="Bernoulli")
         yCurr <- rbinom(nVec[i],rep(1,nVec[i]),1/(1+exp(-etaCurr)))
      
      y <- c(y,yCurr)
      
      currEnd <- currStt + length(yCurr) - 1
      reBlockInds[i] <- list(currStt:currEnd)
      currStt <- currEnd + 1
   }
      
   ## Set up general design matrices corresponding to the fixed effects
   ## vector and spline coefficients vector:

   XG <- x2
   
   numIntKnots <- 25 # Number of interior knots
   intKnots <-  quantile(unique(x2),seq(0,1,length=numIntKnots+2)[-c(1,numIntKnots+2)])
   range.x2 <- c(1.01*min(x2)-0.01*max(x2),1.01*max(x2)-0.01*min(x2))
   ZG <- ZOSull(x2,intKnots=intKnots,range.x=range.x2)
   
   ## Set up random group effects design matrices corresponding to the fixed
   ## effects vector and random group effects vector:

   XR <- cbind(rep(1,numObs),x1) 

   ZR <- matrix(0,numObs,2*m)
   for (i in 1:m)
   {
      indsCurr <- (1:numObs)[idnum==i]
      ZR[indsCurr,c((2*i-1),2*i)] <- cbind(rep(1,nVec[i]),x1[indsCurr])
   }

   X <- cbind(XR,XG); CG <- cbind(XR,XG,ZG)

   ## Set up dimension variables for design matrices:

   ncZG <- ncol(ZG); L <- length(ncZG)
   ncXR <- ncol(XR)
   ncX <- ncol(X); ncCG <- ncol(CG)
   
   ## Perform MCMC sampling via Stan:
   
   if (doMCMC)
   {
      MCMCfit <- mlevMCMC(responseTypeVal,compileCode,useSavedFit)
      timeMCMC[isim] <- MCMCfit$timeMCMC       
      
      ## Extract MCMC samples:

      MCMCsamples <- extract(MCMCfit$StanObj,permuted=FALSE)

      betaMCMC <- NULL
      for (j in 1:ncX)
      {
         charVar <- paste("beta[",as.character(j),"]",sep="") 
         betaMCMC <- rbind(betaMCMC,MCMCsamples[,1,
                           dimnames(MCMCsamples)$parameters==charVar])
      }

      uGMCMC <- NULL
      for (k in 1:ncol(ZG))
      {
         charVar <- paste("uG[",as.character(k),"]",sep="") 
         uGMCMC <- rbind(uGMCMC,MCMCsamples[,1,
                         dimnames(MCMCsamples)$parameters==charVar])
      }
      
      sigmaEpsMCMC <- MCMCsamples[,1,dimnames(MCMCsamples)$parameters=="sigmaEps"]
      Sigma11MCMC <- MCMCsamples[,1,dimnames(MCMCfit$StanObj)$parameters=="SigmaR[1,1]"]
      Sigma12MCMC <- MCMCsamples[,1,dimnames(MCMCfit$StanObj)$parameters=="SigmaR[1,2]"]
      Sigma22MCMC <- MCMCsamples[,1,dimnames(MCMCfit$StanObj)$parameters=="SigmaR[2,2]"]
      
      ## Extract MCMC samples for penalized spline function:

      x1g <- rep(mean(x1),ng)
      x2g <- seq(range.x2[1],range.x2[2],length=ng)
      Xg <- cbind(rep(1,ng),x2g) 
      ZGg <- ZOSull(x2g,range.x2,intKnots)
      fMCMC <- Xg%*%betaMCMC[c(1,3),] + ZGg%*%uGMCMC 
      fhatTrueg <- beta0True + fTrue(x2g)

      indQ1 <- length(x2g[x2g<=quantile(x2,0.2)])
      fQ1MCMC <- fMCMC[indQ1,]
      indQ2 <- length(x2g[x2g<=quantile(x2,0.4)])
      fQ2MCMC <- fMCMC[indQ2,]
      indQ3 <- length(x2g[x2g<=quantile(x2,0.6)])
      fQ3MCMC <- fMCMC[indQ3,]
      indQ4 <- length(x2g[x2g<=quantile(x2,0.8)])
      fQ4MCMC <- fMCMC[indQ4,]

      ## Produce compact graphical summary of MCMC outputs:

      if (responseTypeVal == "Gaussian")
      {
         parms <- list(cbind(fQ1MCMC,fQ2MCMC,fQ3MCMC,fQ4MCMC,
                       betaMCMC[2,],Sigma11MCMC,Sigma12MCMC,
                       Sigma22MCMC,sigmaEpsMCMC^2))

         parNamesVal <- list(expression(f (Q[1])),expression(f (Q[2])),
                             expression(f (Q[3])),expression(f (Q[4])),
                             expression(beta[x]),expression(Sigma[11]^{R}),
                             expression(Sigma[12]^{R}),expression(Sigma[22]^{R}),
                             expression(sigma[epsilon]^2))
         
         summMCMC(parms,parNames=parNamesVal,numerSummCex=1.1,KDEvertLine=FALSE,
                  addTruthToKDE=c(fhatTrueg[c(indQ1,indQ2,indQ3,indQ4)],beta1True,
                                  SigmaL2True[1,1],SigmaL2True[1,2],
                                  SigmaL2True[2,2],sigsqEpsTrue))
      }
      
      if (responseTypeVal == "Bernoulli")
      {
         parms <- list(cbind(fQ1MCMC,fQ2MCMC,fQ3MCMC,fQ4MCMC,
                             betaMCMC[2,],Sigma11MCMC,Sigma12MCMC,
                             Sigma22MCMC))

         parNamesVal <- list(expression(f (Q[1])),expression(f (Q[2])),
                             expression(f (Q[3])),expression(f (Q[4])),
                             expression(beta[x]),expression(Sigma[11]^{R}),
                             expression(Sigma[12]^{R}),expression(Sigma[22]^{R}))

         summMCMC(parms,parNames=parNamesVal,numerSummCex=1.1,KDEvertLine=FALSE,
                  addTruthToKDE=c(fhatTrueg[c(indQ1,indQ2,indQ3,indQ4)],beta1True,
                                  SigmaL2True[1,1],SigmaL2True[1,2],SigmaL2True[2,2]))
      }
      
      MCMC.plots[[isim]] <- recordPlot()
      if (createPDF) dev.off()
      if (userwait) wait()         
   }

   ## Perform mean field variational Bayes approximation:
   
   if (doMFVB)
   {
      timeInfo <- system.time(
                  MFVBfit <- mlevMFVB(y,XR,XG,ZR,ZG,reBlockInds,ncZG,
                                      responseType=responseTypeVal,
                                      doStreamlined=doStreamlinedVal,
                                      maxIter=maxIterVal))

      #saveRDS(MCMCfit,file=paste(getwd(),"/Data/MCMCfit.rds",sep=""))
      #MFVBfit <- readRDS(paste(getwd(),"/Data/MCMCfit.rds",sep=""))
      
      timeMFVB[isim] <- timeInfo[3]

      # Extract MFVB estimates:
      
      if (doStreamlinedVal)
      {
         mu.q.betauG <- MFVBfit$mu.q.betauG
         Sigma.q.betauG <- MFVBfit$Sigma.q.betauG
      }
      B.q.SigmauR <- MFVBfit$B.q.SigmaR 

      if (responseTypeVal=="Gaussian")
      {
         A.q.sigsq.eps <- MFVBfit$A.q.sigsq.eps
         B.q.sigsq.eps <- MFVBfit$B.q.sigsq.eps
      }
      logMLgrid <- MFVBfit$logMLgrid[MFVBfit$logMLgrid!=0]
   }

   ## Write MFVB and MCMC computation time results to file:

   outVec <- c(isim,timeMCMC[isim],timeMFVB[isim])
   write(outVec,filenameTime,n=length(outVec),append=TRUE)
   
   ## Figure 3 (plot 1): sucessive values of variational lower bound to monitor
   ## the convergence of the mean field variational Bayes algorithm:

   if (plotLowerBound)
   {
      par(mfrow=c(1,1),mar=c(5,5,1,1))
      if (!any(is.na(logMLgrid)))
         plot(1:length(logMLgrid),logMLgrid,type="b",bty="l",xlab="Iterations",
              ylab="Lower bound on marginal log-likelihood",cex.axis=cexVal,
              cex.lab=cexVal,lwd=lwdVal)
      logML.plots[[isim]] <- recordPlot()
      signDiffVec <- sign(diff(logMLgrid))
      if (any(signDiffVec!=1)) print("WARNING: Non-monotonic lower bound.")
      if (createPDF) dev.off()
      if (userwait) wait()
   }

   ## Figure 3 (plot 2): fitted function estimates and pointwise 95% credible sets
   ## for both MFVB and MCMC:

   if (plotFunctionFits)
   {
      par(mfrow=c(1,1),mar=c(5,5,1,1))
      obsCol <- "green4" ; colMCMC <- "blue"
      colMFVB <- "darkorange"

      x1g <- rep(mean(x1),ng)
      x2g <- seq(range.x2[1],range.x2[2],length=ng)
      ZGg <- ZOSull(x2g,intKnots=intKnots,range.x=range.x2)
      Xg <- cbind(rep(1,ng),x2g) 
      CGg <- cbind(Xg,ZGg)

      ## MCMC estimates:
      
      fhatMCMCg <- apply(fMCMC,1,mean)
      credLowerfMCMCg <- apply(fMCMC,1,quantile,0.025)
      credUpperfMCMCg <- apply(fMCMC,1,quantile,0.975)

      ## MFVB estimates:
      
      fhatMFVBg <- as.vector(CGg%*%mu.q.betauG[c(1,3:ncCG)])
      sdhatMFVBg <- as.vector(sqrt(diag(CGg%*%Sigma.q.betauG[c(1,3:ncCG),c(1,3:ncCG)]%*%t(CGg))))
      credLowerfMFVBg <- fhatMFVBg - qnorm(0.975)*sdhatMFVBg
      credUpperfMFVBg <- fhatMFVBg + qnorm(0.975)*sdhatMFVBg
      
      plot(0,0,type="n",bty="l",xlim=range(x2g),ylim=range(y),
           xlab="Continuous predictor",ylab="Function estimate",cex.lab=cexVal,
           cex.axis=cexVal,main="",cex.main=1.5)
      points(x2,y,col=obsCol,cex=0.7)
      
      if (responseTypeVal == "Gaussian") 
      {
         lines(x2g,fhatMCMCg,col=colMCMC,lwd=lwdVal,lty=1)
         lines(x2g,credLowerfMCMCg,col=colMCMC,lwd=lwdVal,lty=2)
         lines(x2g,credUpperfMCMCg,col=colMCMC,lwd=lwdVal,lty=2)
         
         lines(x2g,fhatMFVBg,col=colMFVB,lwd=lwdVal,lty=1)
         lines(x2g,credLowerfMFVBg,type="l",col=colMFVB,lty=2)
         lines(x2g,credUpperfMFVBg,type="l",col=colMFVB,lty=2)
         lines(x2g,fhatTrueg,col="red",lwd=lwdVal)
      }
      
      if (responseTypeVal == "Bernoulli")
      {
         lines(x2g,inv.logit(fhatMCMCg),col=colMCMC,lwd=lwdVal,lty=1)
         lines(x2g,inv.logit(credLowerfMCMCg),col=colMCMC,lwd=lwdVal,lty=2)
         lines(x2g,inv.logit(credUpperfMCMCg),col=colMCMC,lwd=lwdVal,lty=2)

         lines(x2g,inv.logit(fhatMFVBg),col=colMFVB,lwd=lwdVal,lty=1)
         lines(x2g,inv.logit(credLowerfMFVBg),col=colMFVB,lwd=lwdVal,lty=2)
         lines(x2g,inv.logit(credUpperfMFVBg),col=colMFVB,lwd=lwdVal,lty=2)
         lines(x2g,inv.logit(fhatTrueg),col="red",lwd=lwdVal)
      }

      legend("topright",legend=c("True","MFVB","MCMC"),bty="n",
             col=c("red",colMFVB,colMCMC),lwd=rep(lwdVal,2),cex=cexVal)
      
      fits.plots[[isim]] <- recordPlot()
      if (createPDF) dev.off()
      if (userwait) wait()
   }

   ## Figure 3: approximate posterior density functions obtained via MFVB and
   ## MCMC. The accuracy scores show the accuracy of MFVB approximation compared
   ## against a MCMC benchmark:

   if (doMFVBvsMCMC)
   {
      colVAval <- "darkorange" ; colMCval <- "blue"
      colTruthVal <- "green4" ; colAccVal <- "purple4"
      colBoxplots <- "lightblue"
   
      cex.mainVal <- 2.5 ; cex.accurVal <- 2.5; cex.labVal <- 2.5
      cex.axisVal <- 2.5 ; cexVal <- 2.5
      
      ## Obtain 95% credible sets of fitted function estimates at the quintiles:

      x1Quin <- quantile(x1,(1:4)/5)
      x2Quin <- quantile(x2,(1:4)/5)
      Xquin <- cbind(rep(1,4),x2Quin)
      Zquin <- ZOSull(x2Quin,intKnots=intKnots,range.x=range.x2) 

      fhatGibbsquin <- Xquin%*%betaMCMC[c(1,3),] + Zquin%*%uGMCMC
      QuinTrue <- fhatTrueg[c(indQ1,indQ2,indQ3,indQ4)]
      Cquin <- cbind(Xquin,Zquin) 

      fhatMFVBquin <- as.vector(Cquin%*%mu.q.betauG[c(1,3:30)])
      sdMFVBquin <- sqrt(diag(Cquin%*%Sigma.q.betauG[c(1,3:30),c(1,3:30)]%*%t(Cquin)))
      credLowMFVBquin <- fhatMFVBquin - qnorm(0.975)*sdMFVBquin 
      credUppMFVBquin <- fhatMFVBquin + qnorm(0.975)*sdMFVBquin 
 
      ## Compute accuracy scores for models parameters corresponding to model (31):
      
      par(mfrow=c(3,3),mai=c(0.35,0.35,0.5,0.05))
      accQuin <- rep(NA,4)
      accQuin[1] <- accVarApp(parVec=c(fhatMFVBquin[1],sdMFVBquin[1]^2),
                               parName=c(expression(paste("f(",Q[1],")"))),
                               MCsample=fhatGibbsquin[1,],type="Normal",parTrue=QuinTrue[1],
                               colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                               colAcc=colAccVal,cex.mainVal=cex.mainVal,
                               cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      legend("topleft",legend=c("True","MFVB","MCMC"),bty="n",
             col=c(colTruthVal,colVAval,colMCval),lwd=rep(lwdVal,3),cex=cexVal-1)

      accQuin[2] <- accVarApp(parVec=c(fhatMFVBquin[2],sdMFVBquin[2]^2),
                               parName=c(expression(paste("f(",Q[2],")"))),
                               MCsample=fhatGibbsquin[2,],type="Normal",parTrue=QuinTrue[2],
                               colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                               colAcc=colAccVal,cex.mainVal=cex.mainVal,
                               cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      accQuin[3] <- accVarApp(parVec=c(fhatMFVBquin[3],sdMFVBquin[3]^2),
                               parName=c(expression(paste("f(",Q[3],")"))),
                               MCsample=fhatGibbsquin[3,],type="Normal",parTrue=QuinTrue[3],
                               colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                               colAcc=colAccVal,cex.mainVal=cex.mainVal,
                               cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      accQuin[4] <- accVarApp(parVec=c(fhatMFVBquin[4],sdMFVBquin[4]^2),
                               parName=c(expression(paste("f(",Q[4],")"))),
                               MCsample=fhatGibbsquin[4,],type="Normal",parTrue=QuinTrue[4],
                               colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                               colAcc=colAccVal,cex.mainVal=cex.mainVal,
                               cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      accBetax <- accVarApp(c(mu.q.betauG[2],Sigma.q.betauG[2,2]),
                            expression(paste(q^{symbol("\052",cex=4)},"(",beta[x],")")),
                            betaMCMC[2,],"Normal",parTrue=beta1True,
                            colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                            colAcc=colAccVal,cex.mainVal=cex.mainVal,
                            cex.accur=cex.accurVal,cex.labVal=cex.labVal)
	    
      A.q.SigmauR11 <- ((nuVal + ncXR - 1 + m)- ncXR + 1)/2
      accSigmaR11 <- accVarApp(c(A.q.SigmauR11,B.q.SigmauR[1,1]/2),
                             expression(paste(q^{symbol("\052",cex=4)},"(",Sigma[11]^{R},")")),
                             Sigma11MCMC,"InvGamma",parTrue=SigmaL2True[1,1],
                             colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                             colAcc=colAccVal,cex.mainVal=cex.mainVal,
                             cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      accSigmaR22 <- accVarApp(c(A.q.SigmauR11,B.q.SigmauR[2,2]/2),
                             expression(paste(q^{symbol("\052",cex=4)},"(",Sigma[22]^{R},")")),
                             Sigma22MCMC,"InvGamma",parTrue=SigmaL2True[2,2],
                             colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                             colAcc=colAccVal,cex.mainVal=cex.mainVal,
                             cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      A.q.SigmaR12 <- (nuVal + ncXR - 1) + m
      accSigmaR12 <- accVarApp(list(A.q.SigmaR12,B.q.SigmauR,1,2),
                             expression(paste(q^{symbol("\052",cex=4)},"(",Sigma[12]^{R},")")),
                             Sigma12MCMC,"OffDiagWishart",parTrue=SigmaL2True[1,2],
                             colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                             colAcc=colAccVal,cex.mainVal=cex.mainVal,
                             cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      if (responseTypeVal=="Gaussian")
         accSigsqEps <- accVarApp(c(A.q.sigsq.eps,B.q.sigsq.eps),
                                    expression(paste(q^{symbol("\052",cex=4)},
                                               "(",sigma[epsilon]^2,")")),
                                    sigmaEpsMCMC^2,"InvGamma",parTrue=sigsqEpsTrue,
                                    colVA=colVAval,colMC=colMCval,colTruth=colTruthVal,
                                    colAcc=colAccVal,cex.mainVal=cex.mainVal,
                                    cex.accur=cex.accurVal,cex.labVal=cex.labVal)

      accu.plots[[isim]] <- recordPlot()
      if (userwait) wait()
      if (createPDF) dev.off()

      # Write accuracy results to file:

      if (responseTypeVal == "Gaussian")
      {
         outVec <- c(isim,accQuin[1],accQuin[2],accQuin[3],accQuin[4],
                     accBetax,accSigmaR11,accSigmaR22,accSigmaR12,accSigsqEps)
      }
      if (responseTypeVal == "Bernoulli")
      {
         outVec <- c(isim,accQuin[1],accQuin[2],accQuin[3],accQuin[4],
                     accBetax,accSigmaR11,accSigmaR22,accSigmaR12)
      }
      write(round(outVec),filenameAcc,ncolumns=length(outVec),append=TRUE)  

      ## Figure 4: Side-by-side boxplots of accuracy scores for MFVB approximation
      ## against MCMC over [100] runs:

      if (isim==isimEND)   
         source(paste(getwd(),"/Code/Functions/simuBoxplots.R",sep=""))
   }

   ## Calculate coverage percentage summaries of simulation result:
   
   if (doCoverSim)
   {
      # Set up the lower and upper bounds for penalized spline elements:
   
      inCredInt1 <- rep(0,4)
      for (idn in 1:4)
      {    
         if ((QuinTrue[idn]>=credLowMFVBquin[idn])&(QuinTrue[idn]<=credUppMFVBquin[idn]))
            inCredInt1[idn] <- 1
      }
   
      # Set up the lower and upper bounds for `betax':
   
      inCredInt2 <- 0
      credLowBeta1 <- mu.q.betauG[2] - qnorm(0.975)*sqrt(Sigma.q.betauG[2,2])
      credUppBeta1 <- mu.q.betauG[2] + qnorm(0.975)*sqrt(Sigma.q.betauG[2,2])
      if ((beta1True>=credLowBeta1)&(beta1True<=credUppBeta1)) inCredInt2 <- 1
         
      # Set up the lower and upper bounds for `MainDiagInverseWishart' elements:
   
      inCredInt3 <- rep(0,ncXR)
      for(idg in 1:ncXR)
      {
   	    credLowMainDiag <- 1/qgamma(0.975,A.q.SigmauR11,B.q.SigmauR[idg,idg]/2)
         credUppMainDiag <-1/qgamma(0.025,A.q.SigmauR11,B.q.SigmauR[idg,idg]/2)
         if ((SigmaL2True[idg,idg]>=credLowMainDiag)&(SigmaL2True[idg,idg]<=credUppMainDiag))
            inCredInt3[idg] <- 1
      } 
   
      # Set up the lower and upper bounds for `OffDiagInverseWishart' elements:
   
      inCredInt4 <- 0; inCredIntMat4 <- diag(ncXR)
      for (idw in 1:(ncXR-1))
      {
         for(jdw in (idw+1):ncXR)
         {
            nMC <- 10000
            xSamp <- rep(NA,nMC) 
   
            for(iMC in 1:nMC)
            {
               xSamp[iMC] <- solve(rwish(A.q.SigmaR12,solve(B.q.SigmauR)))[idw,jdw]
            }
   
            credLowOffDiag <- quantile(xSamp,0.025,names=FALSE)
            credUppOffDiag <- quantile(xSamp,0.975,names=FALSE)
   
            if ((SigmaL2True[idw,jdw]>=credLowOffDiag)&(SigmaL2True[idw,jdw]<=credUppOffDiag))
               inCredIntMat4[idw,jdw] <- 1  
         }
      } 
      inCredInt4 <- inCredIntMat4[1,2]
   
      if (responseTypeVal == "Gaussian") 
      {
         # Set up the lower and upper bounds for `sigsqeps':
   
         inCredInt5 <- 0
         credLowSigsqeps <- 1/qgamma(0.975,A.q.sigsq.eps,rate=B.q.sigsq.eps)
         credUppSigsqeps <- 1/qgamma(0.025,A.q.sigsq.eps,rate=B.q.sigsq.eps)
         if ((sigsqEpsTrue>=credLowSigsqeps)&(sigsqEpsTrue<=credUppSigsqeps))
            inCredInt5 <- 1
      }

      # Write coverage results to file:
      
      if (responseTypeVal == "Gaussian") 
         outVec <- c(isim,inCredInt1,inCredInt2,inCredInt3,inCredInt4,inCredInt5) 
      if (responseTypeVal == "Bernoulli") 
         outVec <- c(isim,inCredInt1,inCredInt2,inCredInt3,inCredInt4) 

      write(outVec,filenameCov,n=length(outVec),append=TRUE)
   
      if (isim==isimEND)   
         source(paste(getwd(),"/Code/Functions/simuCoverage.R",sep="")) 
   }
}

## Create pdf output files:

if (createPDF) 
{
   pdf(paste(output.path,"2LevModMCMC_",responseTypeVal,".pdf",sep=""),width=10)
   for (MCMC.plot in MCMC.plots)
   {
      replayPlot(MCMC.plot)
   }
   dev.off()

   pdf(paste(output.path,"2LevModlogML_",responseTypeVal,".pdf",sep=""),width=8)
   for (logML.plot in logML.plots)
   {
      replayPlot(logML.plot)
   }
   dev.off()

   pdf(paste(output.path,"2LevModfits_",responseTypeVal,".pdf",sep=""),width=8)
   for (fits.plot in fits.plots)
   {
      replayPlot(fits.plot)
   }
   dev.off()

   pdf(paste(output.path,"2LevModaccu_",responseTypeVal,".pdf",sep=""),width=8)
   for (accu.plot in accu.plots)
   {
      replayPlot(accu.plot)
   }
   dev.off()
}

# Read in txt file and output the mean and standard deviation of the
# total computation time:

if (doRecordTime)
{
   timedata <- read.table(filenameTime,header=T)
       
   cat("m (group size) = "); print(m)
   cat("ResponseType = "); print(responseTypeVal)
   cat("Streamlined? = "); print(doStreamlinedVal)
   cat("Mean MCMC = "); print(mean(timedata[,2]))
   cat("Std MCMC = "); print(sqrt(var(timedata[,2])))
   cat("Mean MFVB = "); print(mean(timedata[,3]))
   cat("Std MFVB = "); print(sqrt(var(timedata[,3])))
}	   
	   
######### End of mlevApproxInference.R ##########
