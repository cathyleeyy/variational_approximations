## ----------------------------------------------------------------------------
## Name : simuBoxplots.R
## ----------------------------------------------------------------------------
## Authors: Cathy Y. Y. Lee 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Required Packages: lattice
## ----------------------------------------------------------------------------
## Description: Obtain boxplots summaries of variational approximation
##              accuracy results
## ----------------------------------------------------------------------------

filename <- paste(output.path,"accuBoxplots_",responseTypeVal,".pdf",sep="")
if (createPDF) pdf(filename,width=18)

## Read in simulation results:

simRes <- read.table(filenameAcc,header=TRUE)[-1]
parNames <- names(simRes)

## Create data matrix for lattice graphics:

accVec <- as.vector(as.matrix(simRes))
parVec <- NULL
for (j in 1:ncol(simRes))
   parVec <- c(parVec,rep(as.character(parNames[j]),nrow(simRes)))
dataMat <- data.frame(accVec,parVec)
names(dataMat) <- c("Accuracy","parmFac")

## Obtain side-by-side boxplots:

tmp <- trellis.par.get("box.rectangle")
tmp$col <- "navy"
trellis.par.set("box.rectangle",tmp)

tmp <- trellis.par.get("box.umbrella")
tmp$col <- "navy"
trellis.par.set("box.umbrella",tmp)

tmp <- trellis.par.get("dot.symbol")
tmp$col <- "navy"
trellis.par.set("dot.symbol",tmp)

tmp <- trellis.par.get("plot.symbol")
tmp$col <- "navy"
trellis.par.set("plot.symbol",tmp)

tmp <- trellis.par.get("background")
tmp$col <- "white"
trellis.par.set("background",tmp)

tmp <- trellis.par.get("strip.background")
tmp$col <- c("honeydew","lavender")
trellis.par.set("strip.background",tmp)

tmp <- trellis.par.get("box.dot")
tmp$cex <- 0.5
trellis.par.set("box.dot",tmp)

numTicks <- 8
tmp <- trellis.par.get("par.xlab.text")
tmp$cex <- 1.8
trellis.par.set("par.xlab.text",tmp)

tmp <- trellis.par.get("par.ylab.text")
tmp$cex <- 2.8
trellis.par.set("par.ylab.text",tmp)

tmp <- trellis.par.get("add.text")
tmp$cex <- 2.65
trellis.par.set("add.text",tmp)

xlabVec <- c(expression(paste("f(",Q[1],")")),expression(paste("f(",Q[2],")")),
             expression(paste("f(",Q[3],")")),expression(paste("f(",Q[4],")")),
             expression(beta[x]),expression(Sigma[11]^{R}),expression(Sigma[22]^{R}),
             expression(Sigma[12]^{R}),expression(sigma[epsilon]^2))

dataMat$parmFac <- ordered(dataMat$parmFac,
                   levels=c("Quin1","Quin2","Quin3","Quin4","betax",
                            "SigmaR11","SigmaR22","SigmaR12","sigsqeps"))

plotObj <- bwplot(Accuracy~parmFac,data=dataMat,ylim=c(70,100),
                  scales=list(x=list(labels=xlabVec,cex=2.8),
                              y=list(labels=seq(70,100,5),cex=2.8)),cex.lab=3,
                  panel=function(...) 
                  {
                     hVals <- seq(70,100,by=5)
                     vVals <- 1:10
                     panel.abline(h=hVals,col="grey75")
                     panel.abline(v=vVals,col="grey75") 
                     panel.bwplot(fill="dodgerblue",...)
                  }
                  )
plot(plotObj)
if (createPDF) dev.off()

############ End of simuBoxplots.R ############

