## ----------------------------------------------------------------------------
## Name : ZOSull.r
## ----------------------------------------------------------------------------
## Authors: Matt P. Wand and John Ormerod 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Compute the Z matrix corresponding to O'Sullvan penalized
##              splines as defined by Wand & Ormerod (2008)
## ----------------------------------------------------------------------------
## Required Packages: splines
## ----------------------------------------------------------------------------
## Usage: ZOSull(x,range.x,intKnots,drv=0)
##
##        x:        vector containing the abscissae over which each column 
##                  of the Z matrix is computed. Typically, x is either a  
##                  univariate set of predictor data, or a plotting grid.
##
##        range.x:  vector of length 2 specifying the lower and upper anchor
##                  points for B-spline computation. It is required that
##                  range.x[1] is less than or equal to min(x), and range.x[2]
##                  is greater than or equal to max(x). The default value is
##                  c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x)).
##                     
##        intKnots: vector specifying the interior knot locations.
##                  The default value is as follows:
##                  numIntKnots <- min(length(unique(x)),35)
##                  intKnots <- quantile(unique(x),seq(0,1,length=
##                                       (numIntKnots+2))[-c(1,(numIntKnots+2))])
##
##        drv:      integer between 0 and 2 specifying the derivative of the 
##                  signal which is being estimated. The default value is 0.
## ----------------------------------------------------------------------------

library(splines)

ZOSull <- function(x,range.x,intKnots,drv=0)
{
   ## Check legality of range.x:
  
   if (!missing(range.x))
   {
      if (length(range.x)!=2) stop("range.x must be of length 2.")
      if (range.x[1]>range.x[2]) stop("range.x[1] exceeds range.x[1].")
      if (range.x[1]>min(x)) stop("range.x[1] must be <= than min(x).")
      if (range.x[2]<max(x)) stop("range.x[2] must be >= than max(x).")
   }
   
   if (drv>2) stop("splines not smooth enough for more than 2 derivatives")

   ## Set defaults for `range.x' and `intKnots':

   if (missing(range.x))
      range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))
   
   if (missing(intKnots))
   {
      numIntKnots <- min(length(unique(x)),35)
      intKnots <- quantile(unique(x),seq(0,1,length=
                  (numIntKnots+2))[-c(1,(numIntKnots+2))])
   }
   numIntKnots <- length(intKnots) 

   ## Obtain the penalty matrix:

   allKnots <- c(rep(range.x[1],4),intKnots,rep(range.x[2],4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- t(Bdd*wts)%*%Bdd     

   ## Use the spectral decomposition of Omega to obtain Z:

   svdOmega <- svd(Omega) 
   indsZ <- 1:(numIntKnots+2)
   UZ <- svdOmega$u[,indsZ] 
   LZ <- t(t(UZ)/sqrt(svdOmega$d[indsZ]))

   ## Perform stability check: 

   indsX <- (numIntKnots+3):(numIntKnots+4)
   UX <- svdOmega$u[,indsX]   
   L <- cbind(UX,LZ)
   stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
   if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
       print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")

   ## Obtain B and post-multiply by LZ matrix to get Z:

   B <- spline.des(allKnots,x,derivs=rep(drv,length(x)),
                     outer.ok=TRUE)$design  
   
   Z <- B%*%LZ

   ## Add the `range.x' and 'intKnots' as attributes
   ## of the return object:

   attr(Z,"range.x") <- range.x
   attr(Z,"intKnots") <- intKnots

   ## Return Z matrix with 2 attributes:

   return(Z)
}

########## End of ZOSull.r ##########

