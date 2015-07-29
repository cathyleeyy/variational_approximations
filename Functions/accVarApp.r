## ----------------------------------------------------------------------------
## Name : accVarApp.r
## ----------------------------------------------------------------------------
## Authors: Matt P. Wand 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Assessment of the accuracy of a variational approximate
#               posterior density function
## ----------------------------------------------------------------------------
## Required Packages: KernSmooth
## ----------------------------------------------------------------------------
## Usage: accVarApp(parVec,parName,MCsample,type,plotDensities=TRUE,
##                  parTrue=NULL,nMC=10000,
##                  colVA="blue",colMC="darkorange",
##                  colTruth="seagreen",colAcc="DarkMagenta",
##                  xlimVal=NULL,gridSize=1001,lwdVal=2,ltyVA=1,ltyMC=1,
##                  cex.labVal=1.5,cex.mainVal=1.2,cex.accur=3,
##                  xlabVal="",ylabVal="",cex.axisVal=1.5)
##
## If type="Normal"            then parVec = (mean,variance)
## If type="InvGamma"          then parVec = (shape,rate)
## If type="InvGammaSqrt"      then parVec = (shape,rate)
## If type="Beta"              then parVec = (shape1,shape2)
## If type="Bernoulli"         then parVec = (mean)
## If type="Discrete"          then parVec = list(atoms,probabs)
## If type="Normal-Zero"       then parVec = (wtNormal,mean,variance)
## If type="NormalMix"         then parVec = list(wVec,meanVec,sdVec)
## If type="InvGammaMix"       then parVec = list(wVec,shapeVec,rateVec)
## If type="OffDiagWishart"    then parVec = list(shape,rateMatrix,rowNum,colNum)
## If type="scrF"              then parVec = (qq,rr,ss,tt)
## If type="scrG"              then parVec = (pp,qq,rr,ss)
## If type="scrJLN"            then parVec = (qq,rr,ss)
## If type="scrJplus"          then parVec = (pp,qq,rr)
## If type="scrJplusRecip"     then parVec = (pp,qq,rr)
## If type="scrJplusRecipMix"  then parVec = list(wVec,ppVec,qqVec,rrVec)
## If type="scrHrecip"         then parVec = (pp,qq,rr)
## ----------------------------------------------------------------------------

require("KernSmooth")

accVarApp <- function(parVec,parName,MCsample,type,plotDensities=TRUE,
                      parTrue=NULL,nMC=10000,
                      colVA="blue",colMC="darkorange",
                      colTruth="seagreen",colAcc="DarkMagenta",
                      xlimVal=NULL,gridSize=1001,lwdVal=2,ltyVA=1,ltyMC=1,
                      cex.labVal=1.5,cex.mainVal=1.2,cex.accur=3,
                      xlabVal="",ylabVal="",cex.axisVal=1.5)
{
   MClower <- quantile(MCsample,0.0005)
   MCupper <- quantile(MCsample,0.9995)

   if (type=="Normal")
   {
      VAlower <- parVec[1] - qnorm(0.9995)*sqrt(parVec[2])
      VAupper <- parVec[1] + qnorm(0.9995)*sqrt(parVec[2])
   }  

   if (type=="InvGamma")
   {
      VAlower <- 1/qgamma(0.9995,parVec[1],rate=parVec[2])
      VAupper <- 1/qgamma(0.0005,parVec[1],rate=parVec[2])
   }

   if (type=="InvGammaSqrt")
   {
      VAlower <- 1/sqrt(qgamma(0.9995,parVec[1],rate=parVec[2]))
      VAupper <- 1/sqrt(qgamma(0.0005,parVec[1],rate=parVec[2]))
   }

   if (type=="Beta")
   {
      VAlower <- qbeta(0.0005,parVec[1],parVec[2])
      VAupper <- qbeta(0.9995,parVec[1],parVec[2])
   }

   if (type=="Bernoulli")
   {
      VAlower <- -0.1
      VAupper <-  1.1
   }

   if (type=="Normal-Zero")
   {
      VAlower <- parVec[2] - qnorm(0.9995)*sqrt(parVec[3])
      VAupper <- parVec[2] + qnorm(0.9995)*sqrt(parVec[3])
   }  
   
   if (type=="Discrete")
   {
      VAlower <- 1.05*min(parVec[[1]]) - 0.05*max(parVec[[1]]) 
      VAupper <- 1.05*max(parVec[[1]]) - 0.05*min(parVec[[1]])      
   }
   
   if (type=="NormalMix")
   {
      VAlower <- qnorMix(0.0005,parVec[[1]],parVec[[2]],parVec[[3]])
      VAupper <- qnorMix(0.9995,parVec[[1]],parVec[[2]],parVec[[3]])
   }

   if (type=="InvGammaMix")
   {
      nMix <- length(parVec[[1]])
      xSamps <- matrix(NA,10000,nMix)
     
      for (k in 1:nMix)
           xSamps[,k] <- rgamma(10000,shape=parVec[[2]][k],rate=parVec[[3]][k])

      inds <- sample(1:nMix,10000,replace=TRUE,prob=parVec[[1]])

      xScrJplusMix <- rep(NA,10000)
      for (i in 1:10000) xScrJplusMix[i] <- xSamps[i,inds[i]]
                       
      VAlower <- 1/as.numeric(quantile(xScrJplusMix,0.9995))
      VAupper <- 1/as.numeric(quantile(xScrJplusMix,0.0005))              

   }

   if (type=="OffDiagWishart")
   {
      require("MCMCpack")

      Aval <- parVec[[1]]
      Bval <- parVec[[2]]
      rowNum <- parVec[[3]]
      colNum <- parVec[[4]]
       
      xSamp <- rep(NA,nMC)
      for (iMC in 1:nMC)
      {
          xSamp[iMC] <- solve(rwish(Aval,solve(Bval)))[rowNum,colNum]
      }

      VAlower <-  quantile(xSamp,0.0001)
      VAupper <-  quantile(xSamp,0.9999)
   }
     
   if (type=="scrF")
   {
      xScrF <- rScrF(100000,parVec[1],parVec[2],parVec[3],parVec[4])
      VAlower <- as.numeric(quantile(xScrF,0.0005))
      VAupper <- as.numeric(quantile(xScrF,0.9995))
   }
   
   if (type=="scrG")
   {
      VAlower <- qScrG(0.0005,"scrG",parVec)
      VAupper <- qScrG(0.9995,"scrG",parVec)
   }
   
   if (type=="scrJLN")
   {
      xScrJLN <- exp(rScrJ(100000,0.5*(parVec[1]+1),4*parVec[2],parVec[3])/2)
      VAlower <- as.numeric(quantile(xScrJLN,0.0005))
      VAupper <- as.numeric(quantile(xScrJLN,0.9995))
   }
   
   if (type=="scrJplus")
   {
      xScrJplus <- rScrJplus(100000,parVec[1],parVec[2],parVec[3])     
      VAlower <- as.numeric(quantile(xScrJplus,0.0005))
      VAupper <- as.numeric(quantile(xScrJplus,0.9995))
   }
   
   if (type=="scrJplusRecip")
   {
      xScrJplus <- rScrJplus(100000,parVec[1],parVec[2],parVec[3])     
      VAlower <- 1/as.numeric(quantile(xScrJplus,0.9995))
      VAupper <- 1/as.numeric(quantile(xScrJplus,0.0005))
   }

   if (type=="scrHrecip")
   {
      xScrH <- rScrH(100000,parVec[1],parVec[2],parVec[3])     
      VAlower <- 1/as.numeric(quantile(xScrH,0.9995))
      VAupper <- 1/as.numeric(quantile(xScrH,0.0005))
   }
   
   if (type=="scrJplusRecipMix")
   {
      nMix <- length(parVec[[1]])
      xSamps <- matrix(NA,10000,nMix)
     
      for (k in 1:nMix)
           xSamps[,k] <- rScrJplus(10000,parVec[[2]][k],parVec[[3]][k],parVec[[4]][k])

      inds <- sample(1:nMix,10000,replace=TRUE,prob=parVec[[1]])

      xScrJplusMix <- rep(NA,10000)
      for (i in 1:10000) xScrJplusMix[i] <- xSamps[i,inds[i]]
                       
      VAlower <- 1/as.numeric(quantile(xScrJplusMix,0.9995))
      VAupper <- 1/as.numeric(quantile(xScrJplusMix,0.0005))              

   }

   if (!is.null(xlimVal))
   {
      xgLow <- xlimVal[1] ; xgUpp <- xlimVal[2]
   }
   if (is.null(xlimVal))
   {
      xgLow <- min(MClower,VAlower) ; xgUpp <- max(MCupper,VAupper)
      xlimVal <- c(xgLow,xgUpp)
   }

   xg <- seq(xgLow,xgUpp,length=gridSize) 
      
   if (type=="Normal") VApostg <- dnorm(xg,parVec[1],sqrt(parVec[2]))
   if (type=="InvGamma")
      VApostg <- exp(parVec[1]*log(parVec[2])-lgamma(parVec[1])
                    -(parVec[1]+1)*log(xg)-parVec[2]/xg)
   if (type=="InvGammaSqrt")
      VApostg <- 2*exp(parVec[1]*log(parVec[2])-lgamma(parVec[1])
                    -(2*parVec[1]+1)*log(xg)-parVec[2]/(xg^2))
   if (type=="Beta") VApostg <- dbeta(xg,parVec[1],parVec[2])
   if (type=="NormalMix")
      VApostg <- dnorMix(xg,parVec[[1]],parVec[[2]],parVec[[3]])
   if (type=="InvGammaMix")
   {
      VApostg <- rep(0,length(xg))
      for (j in 1:length(parVec[[1]]))
      {
         newTerm <-  parVec[[1]][j]*exp(parVec[[2]][j]*log(parVec[[3]][j])
                                        - lgamma(parVec[[2]][j])
                                        - (parVec[[2]][j]+1)*log(xg)
                                        - parVec[[3]][j]/xg)
         VApostg <- VApostg + newTerm
      }   
   }
   if (type=="OffDiagWishart")
   {
      ## Obtain kernel density estimate based on Monte Carlo sample:

      dest <- bkde(xSamp,bandwidth=dpik(xSamp),range.x=range(xg),gridsize=gridSize)
      VApostg <- dest$y

   }  
   if (type=="Normal-Zero")
      VApostg <- parVec[1]*dnorm(xg,parVec[2],sqrt(parVec[3]))
   
   if (type=="scrF")
      VApostg <- exp(parVec[1]*log(1+xg^2)-parVec[2]*xg^2+
                     parVec[3]*xg*sqrt(1+xg^2)+parVec[4]*xg
                     -logScrF(0,parVec[1],parVec[2],parVec[3],parVec[4])$modulus) 
   if (type=="scrG")
      VApostg <- dScrG(xg,"scrG",parVec)
   if (type=="scrJLN")
      VApostg <- 2*exp(parVec[1]*log(xg)-parVec[2]*(log(xg)^2)-parVec[3]/(xg^2)
                     -logScrJ(0,0.5*(parVec[1]+1),0.25*parVec[2],0,parVec[3])$modulus)
   if (type=="scrJplus")
      VApostg <- exp(parVec[1]*log(xg)+parVec[2]*xg-parVec[3]*xg^2
                     -logScrJplusR(parVec[1],parVec[2],parVec[3]))
   if (type=="scrJplusRecip")
      VApostg <- exp(-(parVec[1]+2)*log(xg)+parVec[2]/xg-parVec[3]/(xg^2)
                     -logScrJplusR(parVec[1],parVec[2],parVec[3]))
   if (type=="scrHrecip")
      VApostg <- exp(-(parVec[1]+2)*log(xg)-parVec[2]/(xg^2)-log(parVec[3]+xg^2)
                     -logScrH(parVec[1],parVec[2],parVec[3]))
        
   if (type=="scrJplusRecipMix")
   {
      VApostg <- rep(0,length(xg))
      for (k in 1:length(parVec[[1]]))
         VApostg <- VApostg + parVec[[1]][k]*exp(
                    -(parVec[[2]][k]+2)*log(xg)+parVec[[3]][k]/xg-parVec[[4]][k]/(xg^2)
                    -logScrJplusR(parVec[[2]][k],parVec[[3]][k],parVec[[4]][k]))
   }


   if (any(type==c("Normal","InvGamma","InvGammaSqrt","Beta","NormalMix","InvGammaMix",
                   "OffDiagWishart","scrF","scrG","scrJLN","scrJplus","scrJplusRecip",
                   "scrHrecip","scrJplusRecipMix")))
   {
      ## Compute accuracy value:
     
      MCpostg <- bkde(MCsample,bandwidth=dpik(MCsample),
                      range.x=range(xg),gridsize=gridSize)$y
      accurVal <- round(100-50*trapint(xg,abs(MCpostg-VApostg)))
      
      if (plotDensities)
      {
         ylimVal <- range(c(VApostg,MCpostg))
         
         plot(xg,VApostg,type="n",bty="l",main=parName,cex.main=cex.mainVal,
              xlab=xlabVal,ylab=ylabVal,cex.lab=cex.labVal,
              ylim=ylimVal,cex.axis=cex.axisVal)
         lines(xg,VApostg,lwd=lwdVal,col=colVA,lty=ltyVA)
         lines(xg,MCpostg,lwd=lwdVal,col=colMC,lty=ltyMC)
   
         if (!is.null(parTrue)) lines(rep(parTrue,2),c(0,ylimVal[2]),col=colTruth)
         abline(h=0)
              
         ## Plot accuracy value:
   
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.1*min(ylimVal)+0.9*max(ylimVal),
              as.character(paste(accurVal,"%",sep="")),cex=cex.accur,col=colAcc)
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.3*min(ylimVal)+0.7*max(ylimVal),
              "accuracy",cex=0.5*cex.accur,col=colAcc)
      }      
   }

   if (type=="Bernoulli")
   {
      ## Plot accuracy value:

      muMC <- mean(MCsample)
      muVA <- parVec[1]
      accurVal <- round(100*(1 - abs(muMC - muVA)))

      if (plotDensities)
      {
     
         ylimVal <- c(0,1)
         xlims0 <- c(-0.05,0.05)
         xlims1 <- c(0.95,1.05)

         plot(0,0,type="n",xlim=c(-0.1,1.1),ylim=c(0,1),bty="l",main=parName,
              cex.main=cex.mainVal,cex.axis=cex.axisVal)

         x0VApoly <- list(x=c(xlims0,rev(xlims0)),y=c(rep(0,2),rep((1-muVA),2)))
         x0MCpoly <- list(x=c(xlims0,rev(xlims0)),y=c(rep(0,2),rep((1-muMC),2)))

         polygon(x0VApoly,col=colVA,density=5)
         polygon(x0MCpoly,col=colMC,density=5,angle=135)

         x1VApoly <- list(x=c(xlims1,rev(xlims1)),y=c(rep(0,2),rep(muVA,2)))
         x1MCpoly <- list(x=c(xlims1,rev(xlims1)),y=c(rep(0,2),rep(muMC,2)))

         polygon(x1VApoly,col=colMC,density=15)
         polygon(x1MCpoly,col=,density=15,angle=135)

         lines(c(xlims0[1],xlims1[2]),rep(0,2))

         text(0.1*min(xlimVal)+0.9*max(xlimVal),0.1*min(ylimVal)+0.9*max(ylimVal),
              as.character(paste(accurVal,"%",sep="")),cex=cex.accur,col=colAcc)
         text(0.1*min(xlimVal)+0.9*max(xlimVal),0.3*min(ylimVal)+0.7*max(ylimVal),
              "accuracy",cex=0.5*cex.accur,col=colAcc)
      }
         
    }

   if (type=="Normal-Zero")
   {
      ## Obtain numbers of unique and zero values in
      ## the MC sample:

      nMC <- length(MCsample)
      indsZero <- (1:nMC)[MCsample==0]
      nUqMC <- length(unique(MCsample))

      ## Obtain the weight values for both
      ## MC and VA approximations:
      
      MCwtZero <- length(indsZero)/nMC
      MCwtNorm <- 1 - MCwtZero
      VAwtNorm <- parVec[1]
      VAwtZero   <- 1 - VAwtNorm

      ## Compute the L1 distance for the discrete part:

      L1dct <- abs(MCwtZero-VAwtZero)
            
      if (plotDensities)
      {
         ## Determine if enough data for density estimation:

         if (nUqMC<10) densEst <- FALSE
         if (nUqMC>=10) densEst <- TRUE
         
         if (!densEst)
         {
            xlimVal <- c(0,1) ; ylimVal <- c(0,1)
            plot(0,0,type="n",bty="l",main=parName,cex.main=cex.mainVal,
                 xlab=xlabVal,ylab=ylabVal,cex.lab=cex.labVal,
                 ylim=ylimVal,cex.axis=cex.axisVal)
            text(0.5,0.5,"Not enough data for continuous comparison.",
                 col="green3")
            L1cts <- NA
         } 
         if (densEst)
         {
            MCsampleCts <- MCsample
            if (length(indsZero)>0)
               MCsampleCts <- MCsample[-indsZero]

            gsVal <- 1001
            h <- dpik(MCsampleCts,gridsize=gsVal)
            destObj <- bkde(MCsampleCts,bandwidth=h,
                            range.x=range(xg),gridsize=gsVal)
            
            MCpostg <- MCwtNorm*destObj$y
            xgFine <- destObj$x
   
            xlimVal <- range(xgFine)
            ylimVal <- range(c(VApostg,MCpostg))
            
            plot(xg,VApostg,type="n",bty="l",main=parName,cex.main=cex.mainVal,
                 xlab=xlabVal,ylab=ylabVal,cex.lab=cex.labVal,
                 ylim=ylimVal,cex.axis=cex.axisVal)
         
            lines(xg,VApostg,lwd=lwdVal,col=colVA,lty=ltyVA)
            lines(xgFine,MCpostg,lwd=lwdVal,col=colMC,lty=ltyMC)

            if (!is.null(parTrue)) lines(rep(parTrue,2),c(0,ylimVal[2]),col=colTruth)
            abline(h=0)

            L1cts <- trapint(xg,abs(MCpostg-VApostg))
            
         }

         ## Add `bar plot' showing the accuracy of the zero component:

         xpolLim <- c(0.85*min(xlimVal)+0.15*max(xlimVal),
                      0.70*min(xlimVal)+0.30*max(xlimVal))

         ypolLow <- 0.65*min(ylimVal)+0.35*max(ylimVal)
         ypolUpp <- 0.35*min(ylimVal)+0.65*max(ylimVal)
         
         framePolYvec <- c(ypolLow,ypolUpp)

         VApolHt <- ypolLow + VAwtZero*(ypolUpp - ypolLow)
         MCpolHt <- ypolLow + MCwtZero*(ypolUpp - ypolLow)

         framePoly <- list(x=c(xpolLim,rev(xpolLim)),
                           y=c(rep(ypolLow,2),rep(ypolUpp,2)))

         VApoly <- list(x=c(xpolLim,rev(xpolLim)),
                        y=c(rep(ypolLow,2),rep(VApolHt,2)))
         MCpoly <- list(x=c(xpolLim,rev(xpolLim)),
                        y=c(rep(ypolLow,2),rep(MCpolHt,2)))

         polygon(framePoly,col="white",border="forestgreen",lwd=2)
         polygon(VApoly,col=colVA,density=20,border=colVA)
         polygon(MCpoly,col=colMC,density=20,angle=135,border=colMC)
     
         ## Compute accuracy value: 
         
         if (!is.na(L1cts))
           accurVal <- round(100*(1 - 0.5*(L1dct + L1cts)))

         if (is.na(L1cts))
            accurVal <- round(100*(1 - 0.5*L1dct))
         
         ## Add accuracy value to plot:
         
         if (densEst) accCharStr <- as.character(paste(accurVal,"%",sep=""))
         if (!densEst) accCharStr <- as.character(paste(">",accurVal,"%",sep=""))
         
         ## Plot accuracy value:
   
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.1*min(ylimVal)+0.9*max(ylimVal),
              accCharStr,cex=cex.accur,col=colAcc)
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.3*min(ylimVal)+0.7*max(ylimVal),
                 "accuracy",cex=0.5*cex.accur,col=colAcc)
         
      }
   }   

   if (type=="Discrete")
   {
      VAatoms <- parVec[[1]] ; VAprobs <- parVec[[2]]
      VAsample <- sample(VAatoms,100000,replace=TRUE,prob=VAprobs) 

      atomHist <- intersect(signif(VAatoms,3),signif(unique(MCsample),3))

      if (length(atomHist)==0)
         atomHist <- VAatoms
      
      gapVal <- min(diff(atomHist))

      breaksLow <- min(atomHist)+0.5*gapVal
      breaksUpp <- max(atomHist)-0.5*gapVal
      
      breaksLowEnough <- FALSE
      while (!breaksLowEnough)
      {
         if (breaksLow < min(min(MCsample),min(VAsample)))
            breaksLowEnough <- TRUE
         breaksLow <- breaksLow - gapVal
      }

      breaksUppEnough <- FALSE
      while (!breaksUppEnough)
      {
         if (breaksUpp > max(max(MCsample),max(VAsample)))
            breaksUppEnough <- TRUE
         breaksUpp <- breaksUpp + gapVal
      }
      
      breaksHist <- seq(breaksLow,breaksUpp,by=gapVal)
      
      histObjMC <- hist(MCsample,breaks=breaksHist,plot=FALSE)
      histObjVA <- hist(VAsample,breaks=breaksHist,plot=FALSE)
      
      densMC <- histObjMC$counts/sum(histObjMC$counts)
      densVA <- histObjVA$counts/sum(histObjVA$counts)
      histObjMC$counts <- densMC
      histObjVA$counts <- densVA

      ylimVal <- range(c(densMC,densVA))
     
      if (plotDensities)
      {
         plot(histObjMC,density=20,col=colMC,border=colMC,
              xlim=c(xgLow,xgUpp),ylim=ylimVal,ylab="probability",
              xlab=xlabVal,main=parName,cex.main=cex.mainVal,
              cex.lab=cex.labVal,cex.axis=cex.axisVal)
                       
         plot(histObjVA,density=20,angle=135,col=colVA,
              border=colVA,ylab="probability",
              xlim=c(xgLow,xgUpp),ylim=ylimVal,add=plotDensities)
      }
         
      accurVal <- round(100 - 50*sum(abs(densMC - densVA)))
      
      if (plotDensities)
      {
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.1*min(ylimVal)+0.9*max(ylimVal),
              as.character(paste(accurVal,"%",sep="")),cex=cex.accur,col=colAcc)
         text(0.15*min(xlimVal)+0.85*max(xlimVal),0.3*min(ylimVal)+0.7*max(ylimVal),
              "accuracy",cex=0.5*cex.accur,col=colAcc)
         if (!is.null(parTrue)) lines(rep(parTrue,2),c(0,ylimVal[2]),col=colTruth)
      }
   }   

   return(accurVal)   
}

########### End of accVarApp.r ############

