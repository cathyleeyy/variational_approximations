## ----------------------------------------------------------------------------
## Name : summMCMC.r
## ----------------------------------------------------------------------------
## Authors: Matt P. Wand 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Summarises the Markov chain Monte Carlo (MCMC) output from a
##              Bayesian analysis
## ----------------------------------------------------------------------------
## Required Packages: KernSmooth
## ----------------------------------------------------------------------------
## Usage: summMCMC(xList,EPSfileName,PDFfileName,PNGfileName,
##                 plotInd=1,parNames,columnHeadings,colourVersion=TRUE,
##                 columnCols,credLevel=0.95,columnHeadCex=3,
##                 paletteNum=1,numerSummCex=1.5,BGRsttPos=10,
##                 BGRyRange=c(0.95,1.25),BGRtickPos=1.2,
##                 BGRlogTransf,BGRlogitTransf,KDExlim,KDEvertLine=TRUE,
##                 KDEvertLineCol="black",addTruthToKDE=NULL,
##                 addBaseLine=FALSE)
##
##         xList:         list of matrices, where each matrix corresponds to a
##                        different chain, and the columns of each matrix correspond
##                        to different parameters. The matrices each have dimension
##                        "numMCMC" by "numParms"; where "numMCMC" is the size of
##                        the MCMC sample and "numParms" is the number of parameters
##                        being summarised.
##
##         EPSfileName:   filename if the summary is to be saved as a (encapsulated)
##                        Postscript file. If this argument and PDFfileName are both
##                        not specified then the summary is printed to the screen.
##
##         PDFfileName:   filename if the summary is to be saved as a PDF file. If 
##                        this argument and EPSfileName are both not specified then  
##                        the summary is printed to the screen.
##
##         plotInd:       if "numChains" exceeds 1 then this indicates which chain is
##                        summarised in the non-BGR panels. The BGR panels are
##                        Brooks-Gelman-Rubin diagnostic plots are use all chains.  
##                        The default value is 1.
##
##         parNames:      list containing a vector of character strings for the 
##                        parameter names. The maximum length of the vector is 3.
##
##        columnHeadings: vector containing column headings. The default is:
##                        c("parameter","trace","lag 1","acf","BGR","density","summary")
## 
##         colourVersion: logical flag indicating if summary should be in colour.
##                        The default is TRUE.
##
##         columnCols:    vector containing colours for each column.
##                        The default is:
##                        ("darkmagenta", "green4","darkorange","dodgerblue",
##                         "darkgoldenrod1","red","navy")
##
##         credLevel:     number between 0 and 1 specifying the credible set level.
##                        The default is 0.95.
##
##         numerSummCex:  positive number specifying character expansion factor for the
##                        numerical summary (last column).
##
##         BGRsttPos:     starting position for the Brooks-Gelman-Rubin plots.
##                        The default value is 10.
##
##         BGRyRange:     vertical axis limits for the Brooks-Gelman-Rubin plots.
##                        The default value is  c(0.95,1.25).
##
##         BGRtickPos:    position of tick mark on vertical axis for the
##                        Brooks-Gelman-Rubin plots. The default value is 1.2.
##
##         BGRlogTransf:  vector containing indices of those parameters for which the
##                        Brooks-Gelman-Rubin plots should be done on a logarithmic scale.
##
##        BGRlogitTransf: vector containing indices of those parameters for which the
##                        Brooks-Gelman-Rubin plots should be done on a logit scale.
##
##         KDExlim:       list of vectors of length 2 specifying the horizontal axis
##                        limits for the kernel density estimates.
##          
##         KDEvertLine:   logical flag indicating if a vertical line at zero should be
##                        added to the kernel density estimates. The default value is TRUE.
##
##        KDEvertLineCol: colour of the vertical line at zero for kernel density estimates.
##                        The default value is "black".
##
##         addTruthToKDE: vector indicating `true' values of parameters.
##                        The default value is NULL. If addTruthToKDE is non-NULL
##                        then dashsed vertical lines corresponding to true values are added.
## ----------------------------------------------------------------------------

library(KernSmooth)

summMCMC <- function(xList,EPSfileName,PDFfileName,PNGfileName,
                     plotInd=1,parNames,columnHeadings,colourVersion=TRUE,
                     columnCols,credLevel=0.95,columnHeadCex=3,
                     paletteNum=1,numerSummCex=1.5,BGRsttPos=10,
                     BGRyRange=c(0.95,1.25),BGRtickPos=1.2,
                     BGRlogTransf,BGRlogitTransf,KDExlim,KDEvertLine=TRUE,
                     KDEvertLineCol="black",addTruthToKDE=NULL,
                     addBaseLine=FALSE)
{
   options(warn=-1)
   x <- lapply(xList,as.matrix)
   num.par <- ncol(x[[plotInd]])
   samp.size <- nrow(x[[plotInd]])
   num.chains <- length(x)

   ## Divert figure to EPS, PDF or PNG file if filename specified:

   if (!missing(EPSfileName))
       postscript(EPSfileName,horizontal=FALSE,width=11,height=(num.par+1))

   if (!missing(PDFfileName))
       pdf(PDFfileName,width=9,height=(num.par+1))

   if (!missing(PNGfileName))
      png(PNGfileName,width=630,height=70*(num.par+1))
            
   if (missing(PNGfileName)) op <- par()
   if (missing(columnHeadings))
      columnHeadings <- c("parameter","trace","lag 1","acf","BGR",
                           "density","summary")   

   if (colourVersion)
   {

      if (missing(columnCols))
      {
         if (paletteNum==1)
         {
            columnCols <- c("darkmagenta", "green4","darkorange","dodgerblue",
                            "darkgoldenrod1","red","navy")
            if (KDEvertLine&missing(KDEvertLineCol))
               KDEvertLineCol <- "DarkGreen"
         }
         if (paletteNum==2)
         {
            columnCols <- c("purple4","tomato","mediumblue",
                            "olivedrab4","DarkCyan","darkred","DarkGreen")
            if (KDEvertLine&missing(KDEvertLineCol))
               KDEvertLineCol <- "darkgoldenrod"
         }
      }
   }

   panels.per.par <- length(columnHeadings)
   if (num.chains==1) panels.per.par <- panels.per.par - 1

   if (!colourVersion)
   {
      columnCols <- rep(NA,num.par)
      for (i in 1:length(columnHeadings))
         columnCols[i] <- "black"
   }

   empty.panel <- function(addBaseLineFlag=addBaseLine)
   {
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",
           yaxt="n",xlab="",ylab="",bty="o")
      if (addBaseLineFlag) lines(c(-0.05,1.05),rep(-0.01,2),lwd=1)
      invisible()
   }
  
   par(mfrow=c((num.par+1),panels.per.par))

   headInds <- 1:7
   if (num.chains==1) headInds <- c(1:4,6,7)
   for (i in 1:panels.per.par)
   {
      par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
      empty.panel()
      text(0.5,0.5,columnHeadings[headInds[i]],cex=columnHeadCex,
           col=columnCols[headInds[i]])
   }

   for (j in 1:num.par)
   {
      ## Write the variable name:

      par(ann=F,mar=c(0,0,0,0),xaxt="n",yaxt="n")
    
      empty.panel(addBaseLineFlag=((j==num.par)&addBaseLine))
      if (length(parNames[[j]])==1)
         text(0.5,0.5,parNames[[j]][1],cex=3,col=columnCols[1])
      if (length(parNames[[j]])==2)
      {
         text(0.5,0.7,parNames[[j]][1],cex=3,col=columnCols[1])
         text(0.5,0.3,parNames[[j]][2],cex=3,col=columnCols[1])
      }
      if (length(parNames[[j]])==3)
      {
         text(0.5,0.8,parNames[[j]][1],cex=3,col=columnCols[1])
         text(0.5,0.5,parNames[[j]][2],cex=3,col=columnCols[1])
         text(0.5,0.2,parNames[[j]][3],cex=3,col=columnCols[1])
      }

      ## Do the trace plot:

      plot(x[[plotInd]][,j],xlab="",ylab="",type="l",col=columnCols[2])

      if ((!missing(PNGfileName))&(j==num.par)&addBaseLine)
      {
         xVec <- 1:length(x[[plotInd]][,j])
         yVec <- x[[plotInd]][,j]
         xLow <- 1.05*min(xVec) - 0.05*max(xVec)
         xUpp <- 1.05*max(xVec) - 0.05*min(xVec)
         yLow <- 1.01*min(yVec) - 0.01*max(yVec)
         lines(c(xLow,xUpp),rep(yLow,2),lwd=2,col="black")
      }    
      
      ## Do the lag 1 plot:

      plot(x[[plotInd]][1:(samp.size-1),j],
           x[[plotInd]][2:samp.size,j],xlab="",ylab="",type="n")

      points(x[[plotInd]][1:(samp.size-1),j],
             x[[plotInd]][2:samp.size,j],pch=1,cex=0.5,col=columnCols[3])

      
      if ((!missing(PNGfileName))&(j==num.par)&addBaseLine)
      {      
         xVec <- x[[plotInd]][1:(samp.size-1),j]
         yVec <- x[[plotInd]][2:samp.size,j]
         xLow <- 1.05*min(xVec) - 0.05*max(xVec)
         xUpp <- 1.05*max(xVec) - 0.05*min(xVec)
         yLow <- 1.01*min(yVec) - 0.01*max(yVec)
         lines(c(xLow,xUpp),rep(yLow,2),lwd=2,col="black")
      }    
      
      ## Do the ACF plot:

      ci.col.val <- "black"
      if (colourVersion) ci.col.val <- "blue"
      acf(x[[plotInd]][,j],lag.max=20,col=columnCols[4],
          lwd=2,ci.col=ci.col.val)

      if ((!missing(PNGfileName))&(j==num.par)&addBaseLine)
      {
         acfVec <- as.numeric(acf(x[[plotInd]][,j],lag.max=20,
                                  plot=FALSE)$acf)
         nACF <- length(x[[plotInd]][,j])
         CIlow <- -2/sqrt(nACF)
         yVec <- c(acfVec,CIlow)
         xVec <- 0:20
         xLow <- 1.05*min(xVec) - 0.05*max(xVec)
         xUpp <- 1.05*max(xVec) - 0.05*min(xVec)
         yLow <- 1.01*min(yVec) - 0.01*max(yVec)
         lines(c(xLow,xUpp),rep(yLow,2),lwd=2,col="black")
      }    
      
      if (num.chains>1) # Multiple chains.
      {

         ## Do the Gelman-Rubin plot:

         x.list <- list(x[[1]][,j])
            for (k in 2:num.chains)
               x.list <- c(x.list,list(x[[k]][,j]))

         if ((!missing(BGRlogTransf)))
         {
            if (any(j==BGRlogTransf))
            {
               x.list <- list(log(x[[1]][,j]))
                  for (k in 2:num.chains)
                     x.list <- c(x.list,list(log(x[[k]][,j])))
            }
         }
         if ((!missing(BGRlogitTransf)))
         {
            if (any(j==BGRlogitTransf))
            {
               x.list <- list(log(x[[1]][,j]))
                  for (k in 2:num.chains)
                     x.list <- c(x.list,list(logit(x[[k]][,j])))
            }
         }

         BGR.outp <- BGRinterval(x.list,BGRsttPos)

         ## Extract x-axis and y-axis components of `BGR.outp':
         
         BGR.x <- BGR.outp$x
         BGR.y <- BGR.outp$numer/BGR.outp$denom
 
         lb.wt <- 0.20*length(BGR.x)
         plot(0,0,type="n",xlim=c(-lb.wt,length(BGR.x)),ylim=BGRyRange,
              xlab="",ylab="", bty="l")
         lines(BGR.x,BGR.y,lwd=1,col=columnCols[5],err=-1)
         lines(rep(0,2),BGRyRange,err=-1)

         text(-0.75*lb.wt,BGRtickPos,as.character(BGRtickPos),srt=90,cex=0.8)
         text(-0.75*lb.wt,1,"1.00",srt=90,cex=0.8)
         abline(h=1,lwd=1)
         lines(c(-lb.wt/2,0),rep(BGRtickPos,2))
      }

      ## Do the density plot:

      h <- dpik(x[[plotInd]][,j])
      est <- bkde(x[[plotInd]][,j],bandwidth=h)

      if (!missing(KDExlim))
         xlim.val <-  KDExlim[[j]]   

      if (missing(KDExlim))
         xlim.val <-  range(est$x)

      lb.ht <- 0.2*(max(est$y)-min(est$y))
      ylim.val <- c(-lb.ht,1.1*max(est$y))
      plot(est,type="l",xlab="",ylab="",ylim= ylim.val,
           col=columnCols[6],xlim=xlim.val,lwd=2)
      lines(c(min(est$x),max(est$x)),rep(0,2))

      if (KDEvertLine)
         lines(rep(0,2),c(0,max(est$y)),lwd=3,err=-1,col=KDEvertLineCol)
      
      if (!is.null(addTruthToKDE))
         lines(rep(addTruthToKDE[j],2),c(0,max(est$y)),
               lwd=1,lty=2,,err=-1,col=KDEvertLineCol)
          
      x.labels <- pretty(x[[plotInd]][,j],n=3)
      for (i in 1:length(x.labels))
      {
         text(x.labels[i],-0.6*lb.ht,as.character(x.labels[i]),cex=0.8)
         lines(rep(x.labels[i],2),c(0,-0.1*lb.ht))
      }   

      if ((!missing(PNGfileName))&(j==num.par)&addBaseLine)
      {
         xVec <- xlim.val ; yVec <- ylim.val
         xLow <- 1.05*min(xVec) - 0.05*max(xVec)
         xUpp <- 1.05*max(xVec) - 0.05*min(xVec)
         yLow <- 1.01*min(yVec) - 0.01*max(yVec)
         lines(c(xLow,xUpp),rep(yLow,2),lwd=2,col="black")
      }    

      empty.panel(addBaseLine=(j==num.par))
      text(0.5,0.75,
      paste("posterior mean: ",
             as.character(signif(mean(x[[plotInd]][,j]),3)),sep=""),
             col=columnCols[7],cex=numerSummCex)
      text(0.5,0.5,paste(100*credLevel,"% credible interval: ",
           sep=""),col=columnCols[7],cex=numerSummCex)

      text(0.5,0.25,
      paste("(",
      as.character(signif(quantile(x[[plotInd]][,j],(1-credLevel)/2),3)),",",
      as.character(signif(quantile(x[[plotInd]][,j],(1+credLevel)/2),3)),")",
      sep=""),col=columnCols[7],cex=numerSummCex)
   }

   if (missing(PNGfileName)) par(op)
   invisible()
}

########## End of summMCMC.r ##########
