## ----------------------------------------------------------------------------
## Name : mlevMFVB.r
## ----------------------------------------------------------------------------
## Authors: Cathy Y. Y. Lee and Matt P. Wand
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description:  Mean field variational Bayes fitting of two-level Bayesian
##               mixed effects models with Gaussian or Bernoulli response
## ----------------------------------------------------------------------------
## Required Packages: ---
## ----------------------------------------------------------------------------
## Usage: mlevMFVB(y,XR,XG,ZR=NULL,ZG,reBlockInds,ncZG,
##                 responseType="Gaussian",
##                 doStreamlined=TRUE,
##                 maxIter=35,useMatForDg=TRUE)
##
##        y:             Response vector
##
##        XR:            Random group effects design matrix corresponds to the 
##                       fixed effects vector betaR
##
##        XG:            General design matrix corresponds to the fixed effects
##                       vector betaG
##
##        ZR:            Random group effects design matrix corresponds to the 
##                       random group effects vector uR
##
##        ZG:            General design matrix corresponds to the spline 
##                       coefficients vector uG
##
##        reBlockInds:   Indices for random effects block
##
##        ncZG:          Number of columes in ZG
##
##        responseType:  Gaussian or Bernoulli response
##
##        doStreamlined: Naive or Streamlined approach
##
##        maxIter:       Maximum number of iterations for MFVB algorithms
## ----------------------------------------------------------------------------

mlevMFVB <- function(y,XR,XG,ZR=NULL,ZG,reBlockInds,ncZG,
                     responseType="Gaussian",
                     doStreamlined=TRUE,
                     maxIter=35,useMatForDg=TRUE)
{    
   ## Set adjusted response vector:

   if (responseType=="Gaussian") yAdj <- y
   if (responseType=="Bernoulli")
   {
      yAdj <- y - 0.5
      lambda <- function(x) 
      {
         nzi <- (1:length(x))[x!=0]
         ans <- rep(0.125,length(x))
         ans[nzi] <- tanh(x[nzi]/2)/(4*x[nzi])
         return(ans)
      }
      phi <- function(x) 
      {
         nzi <- (1:length(x))[x!=0]
         ans <- rep(0.125,length(x))
         ans[nzi] <- (x[nzi]/2) - log(1+exp(x[nzi])) + (0.25*x[nzi]*tanh(x[nzi]/2))
         return(ans)
      }
   }
   
   ## Obtain constant matrices and dimension variables:
  
   X <- cbind(XR,XG);  CG <- cbind(X,ZG)

   ncXR <- ncol(XR);   L <- length(ncZG)
   ncX <- ncol(X);     ncCG <- ncol(CG)
   
   m <- length(reBlockInds)
   nVec <- unlist(lapply(reBlockInds,length))
   numObs <- sum(nVec)
  
   if (!doStreamlined)
   {
      C <- cbind(CG,ZR);         CTC <- crossprod(C)
      CTy <- crossprod(C,yAdj);  ncC <- ncol(C)
   }

   if (doStreamlined)
   {
      CGTy <- crossprod(CG,yAdj) ; CGTCG <- crossprod(CG)

      ## Obtain the ZRTy list of matrices:

      ZRTy <- vector("list", length=m)
      for (i in 1:m)
      {
         ZRTy[[i]] <- crossprod(XR[reBlockInds[[i]],],yAdj[reBlockInds[[i]]])
      }
      
      ## Create lists of matrices required for MFVB:
      
      G <- vector("list",length=m) 
      H <- vector("list",length=m)
      zetaVec <- vector("list",length=m)
      kappaVec <- vector("list",length=m)
   }

   ## Do MFVB initialisations:
   
   mu.q.recip.sigsq.eps <- 1
   mu.q.recip.a.eps <- 1
   M.q.inv.SigmaR <- diag(ncXR)  
   mu.q.recip.aR <- rep(1,ncXR)
   mu.q.recip.sigsq.u <- rep(1,L)
   if (responseType=="Bernoulli")
      xiVec <- rep(1,length(y))

   ## Set up function for the lower bound on the marginal log-likelihood:

   logML <- function(nuVal,ncXR,ncX,ncZG,m,numObs,L,yAdj,C,
                     sigsq.beta,fitVecBetauG,fitVecuR,wtVec,phiVec,
                     xiVec,xiSqd,mu.q.betauG,Sigma.q.betauG,det.Sigma.q.betau,
                     B.q.SigmaR,B.q.sigsq.u,A.R,M.q.inv.SigmaR,
                     mu.q.recip.a.R,B.q.a.R,A.u,B.q.a.u,mu.q.recip.a.u,
                     mu.q.recip.sigsq.u,A.eps,B.q.sigsq.eps,B.q.a.eps,
                     mu.q.recip.a.eps,mu.q.recip.sigsq.eps)
   {
      ASigma <- nuVal + ncXR - 1
      logCqRA <- (0.5*ASigma*ncXR*log(2) + 0.25*ncXR*(ncXR-1)*log(pi) +
                  sum(lgamma(0.5*(ASigma+1-(1:ncXR)))))
      logCqRAm <- (0.5*(ASigma+m)*ncXR*log(2) + 0.25*ncXR*(ncXR-1)*log(pi) +
                  sum(lgamma(0.5*(ASigma+m+1-(1:ncXR)))))

      ans0 <- (0.5*ncXR*ASigma*log(2*nuVal)
               - (0.5*ncXR+L+1)*log(pi) - 0.5*ncX*log(sigsq.beta))

      if (responseType=="Bernoulli")
      {
         if (!doStreamlined)
         {
            ans1 <- (crossprod(yAdj,C%*%mu.q.betau)  
                    - sum((wtVec/2)*(xiVec^2)) + sum(phiVec))
         }
         if (doStreamlined)
         {
            ans1 <- (crossprod(yAdj,(fitVecBetauG + fitVecuR)) 
                    - sum((wtVec/2)*(xiVec^2)) + sum(phiVec))
         }
      }
       
      if (!doStreamlined)
         ans2 <- (-1/(2*sigsq.beta)*(sum(mu.q.betau[1:ncX]^2) + 
                  sum(diag(Sigma.q.betau[(1:ncX),(1:ncX)]))))

      if (doStreamlined)
         ans2 <- (-1/(2*sigsq.beta)*(sum(mu.q.betauG[1:ncX]^2) + 
                 sum(diag(Sigma.q.betauG[(1:ncX),(1:ncX)]))))

      ans3 <- (0.5*det.Sigma.q.betau + 0.5*(sum(ncZG)+m+ncX))
      ans4 <- (-logCqRA + logCqRAm -0.5*(ASigma+m)*determinant(B.q.SigmaR)$modulus)
      ans5 <- (sum(lgamma(0.5*(ncZG+1)))-0.5*sum((ncZG+1)*log(B.q.sigsq.u)))
      ans6 <- (-ncXR*log(A.R) + ncXR*lgamma(0.5*(nuVal+ncXR)))
      ans7 <- (sum(nuVal*diag(M.q.inv.SigmaR)*mu.q.recip.a.R))      
      ans8 <- (-0.5*(nuVal+ncXR)*sum(log(B.q.a.R)))
      ans9 <- (-L*log(A.u) - sum(log(B.q.a.u)))
      ans10 <- (sum(mu.q.recip.a.u*mu.q.recip.sigsq.u))
      
      if (responseType=="Gaussian")
      {
         ans11 <- (-0.5*(numObs+1)*log(B.q.sigsq.eps) + lgamma(0.5*(numObs+1))
                   - 0.5*numObs*log(2*pi))
         ans12 <- (-log(A.eps) - log(B.q.a.eps) + mu.q.recip.a.eps*mu.q.recip.sigsq.eps)
         ans <- ans0+ans2+ans3+ans4+ans5+ans6+ans7+ans8+ans9+ans10+ans11+ans12
         return(list(ans0,ans2,ans3,ans4,ans5,ans6,ans7,ans8,ans9,ans10,ans11,ans12,ans))
      }

      if (responseType=="Bernoulli")
      {
         ans <- ans0+ans1+ans2+ans3+ans4+ans5+ans6+ans7+ans8+ans9+ans10
         return(list(ans0,ans1,ans2,ans3,ans4,ans5,ans6,ans7,ans8,ans9,ans10,ans))
      }
   }
   
   logMLcurr <- -1e20
   converged <- FALSE ; itnum <- 0 
   tolerance <- 0.0000001
   logMLgrid <- NULL
   while (!converged) 
   {
      itnum <- itnum + 1

      ## Set current log(ML) value to previous one:
   
      if (itnum == 1) logMLprev <- logMLcurr
      if (itnum > 1)  logMLprev <- logMLcurr[[length(logMLcurr)]][1]
      
      ## Set weighed vectors:
      
      if (responseType=="Gaussian")
         wtVec <- rep(1,length(y))
      if (responseType=="Bernoulli")
      {
         wtVec <- 2*lambda(xiVec)
         phiVec <- phi(xiVec)
      }
      
      ridgeVec <- rep((1/sigsq.beta),ncX)
      for (ell in 1:L)
         ridgeVec <- c(ridgeVec,rep(mu.q.recip.sigsq.u[ell],ncZG[ell]))
      
      if (!doStreamlined)
      {
         ## Set up covariance matrix M.q.Sigma:
          
         if (ncXR == 1)
         {
             ridgeVec2 <- rep(M.q.inv.SigmaR,m)
             M.q.Sigma <- diag(c(ridgeVec,ridgeVec2))
         }
         
         if (ncXR > 1)
         {
            ridgeMat <- kronecker(diag(m),M.q.inv.SigmaR)
            M.q.Sigma <- adiag(diag(ridgeVec),ridgeMat)
         }
         
         ## Update q*(beta,u) parameters:
         
         Sigma.q.betau <- solve(mu.q.recip.sigsq.eps*crossprod(C,wtVec*C) + M.q.Sigma)
         mu.q.betau <- mu.q.recip.sigsq.eps*Sigma.q.betau%*%CTy

         ## Compute the determinant of Sigma.q.betau:
         
         det.Sigma.q.betau <- determinant(Sigma.q.betau)$modulus
      }

      if (doStreamlined)
      {      
         ## Update q*(beta,u) parameters:

         sVec <- rep(0,ncCG); Smat <- matrix(0,ncCG,ncCG)
         for (i in 1:m)
         {            
            G[[i]] <- (mu.q.recip.sigsq.eps*crossprod(CG[reBlockInds[[i]],],
                       wtVec[reBlockInds[[i]]]*XR[reBlockInds[[i]],]))
            H[[i]] <- (solve(mu.q.recip.sigsq.eps*crossprod(XR[reBlockInds[[i]],],
                       wtVec[reBlockInds[[i]]]*XR[reBlockInds[[i]],]) + M.q.inv.SigmaR))
            sVec <- sVec + as.vector(G[[i]]%*%H[[i]]%*%ZRTy[[i]])
            Smat <- Smat + G[[i]]%*%H[[i]]%*%t(G[[i]])
         }

         if (responseType=="Gaussian")
            Sigma.q.betauG <- (solve(mu.q.recip.sigsq.eps*CGTCG
                               + diag(ridgeVec) - Smat))
         
         if (responseType=="Bernoulli")
            Sigma.q.betauG <- (solve(mu.q.recip.sigsq.eps*crossprod(CG,wtVec*CG)
                               + diag(ridgeVec) - Smat))
         mu.q.betauG <- as.vector(mu.q.recip.sigsq.eps*Sigma.q.betauG%*%(CGTy - sVec))
      
         Sigma.q.uR <- vector("list",length=m)
         for (i in 1:m)
         {
            Sigma.q.uR[[i]] <- (H[[i]] + H[[i]]%*%crossprod(G[[i]],
                                (Sigma.q.betauG%*%G[[i]]%*%H[[i]])))
            zetaVec[[i]] <- H[[i]]%*%ZRTy[[i]]
         }
      
         mu.q.uR <- vector("list",length=m)
         for (i in 1:m)
         {
            kappaVec[[i]] <- as.vector(crossprod(G[[i]],Sigma.q.betauG%*%sVec))                
            mu.q.uR[[i]] <- as.vector(H[[i]]%*%(mu.q.recip.sigsq.eps*ZRTy[[i]] -
                                      crossprod(G[[i]],mu.q.betauG)))
         }
         
         ## Compute the determinant of Sigma.q.betau:

         detAMat <- NULL
         for (i in 1:m)
            detAMat <- c(detAMat,determinant(solve(H[[i]]))$modulus)
         
         detBMat <- matrix(0,ncCG,ncCG)
         for (i in 1:m)
            detBMat <- detBMat + G[[i]]%*%H[[i]]%*%t(G[[i]])

         detBMat <- mu.q.recip.sigsq.eps*crossprod(CG,(wtVec*CG)) + diag(ridgeVec) - detBMat         
         det.Sigma.q.betau <- -(sum(detAMat) + determinant(detBMat)$modulus)
                  
         ## Compute values of residual vectors:
      
         fitVecBetauG <- as.vector(CG%*%mu.q.betauG)
         fitVecuR <- NULL
         for (i in 1:m)
         {
            if (ncXR == 1) fitVecuR <- c(fitVecuR,rep(mu.q.uR[[i]],nVec[i]))
            if (ncXR > 1) fitVecuR <- rbind(fitVecuR, XR[reBlockInds[[i]],]%*%mu.q.uR[[i]])
         }
         
         if (responseType=="Gaussian") 
         {
            ## Update q*(sigsq.eps) parameter:
          
            residSS <- sum((y-fitVecBetauG-fitVecuR)^2)
            trTermTwo <- sum(diag(Sigma.q.betauG%*%CGTCG))
            for (i in 1:m)
               trTermTwo <- (trTermTwo + sum(diag(Sigma.q.uR[[i]]
                             %*%crossprod(XR[reBlockInds[[i]],]))))
            trTermThree <- 0
            for (i in 1:m)
               trTermThree <- (trTermThree + 
                               sum(diag(G[[i]]%*%H[[i]]%*%t(G[[i]])%*%Sigma.q.betauG)))
          }
       }

      if (responseType=="Gaussian") 
      {                        
         ## Update q*(1/sigsq.eps) parameters:

         A.q.sigsq.eps <- 0.5*(sum(nVec)+1)

         if (!doStreamlined)
            B.q.sigsq.eps <- (mu.q.recip.a.eps + 0.5*(sum((as.vector(y - C%*%mu.q.betau)^2)) + 
                              sum(diag(CTC%*%Sigma.q.betau))))
         if (doStreamlined)
            B.q.sigsq.eps <- (mu.q.recip.a.eps + 0.5*(residSS + trTermTwo 
                              - 2*(1/mu.q.recip.sigsq.eps)*trTermThree))

         mu.q.recip.sigsq.eps <- A.q.sigsq.eps/B.q.sigsq.eps

         ## Update q*(a.eps) parameter:

         B.q.a.eps <- mu.q.recip.sigsq.eps + (1/A.eps^2)
         mu.q.recip.a.eps <- 1/B.q.a.eps
      }

      ## Set mu.q.recip.sigsq.eps to unity and update xiVec:

      if (responseType=="Bernoulli")
      {
         A.q.sigsq.eps <- NA ; B.q.sigsq.eps <- NA
         mu.q.recip.sigsq.eps <- 1

         if (!doStreamlined)
         {
            EsqMat <- Sigma.q.betau + tcrossprod(mu.q.betau)
            xiVec <- sqrt(diag(C%*%EsqMat%*%t(C)))
         }
         
         if (doStreamlined)
         {
            if (useMatForDg)
               xiSqd <- diag(CG%*%(Sigma.q.betauG + tcrossprod(mu.q.betauG))%*%t(CG))
            if (!useMatForDg)
            {
               xiSqd <- rep(0,sum(nVec))
               for (i in 1:sum(nVec))
                  xiSqd[i] <- (crossprod(CG[i,],(Sigma.q.betauG
                                                + tcrossprod(mu.q.betauG)))%*%CG[i,])
            }

            for (i in 1:m)
            {
               EsqMatCurr <- (-Sigma.q.betauG)%*%G[[i]]%*%H[[i]]
                              + tcrossprod(mu.q.betauG,mu.q.uR[[i]])
               xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]] 
                      + 2*diag(CG[reBlockInds[[i]],]%*%EsqMatCurr%*%t(XR[reBlockInds[[i]],])))
               EsqMatCurr <- Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
               xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]] 
                        + diag(XR[reBlockInds[[i]],]%*%EsqMatCurr%*%t(XR[reBlockInds[[i]],])))
            }
            xiVec <- sqrt(xiSqd)
         }         
      }

      ## Update q*(a.R) parameters:

      B.q.a.R <- nuVal*diag(M.q.inv.SigmaR) + (1/A.R^2)
      mu.q.recip.a.R <- (0.5*(nuVal + ncXR))/B.q.a.R

      ## Update q*(SigmaR^{-1}) parameters:

      if (ncXR == 1) B.q.SigmaR <- 2*nuVal*mu.q.recip.a.R
      if (ncXR > 1) B.q.SigmaR <- 2*nuVal*diag(mu.q.recip.a.R)
      indsStt <- ncCG + 1
      for (i in 1:m)
      {
         if (!doStreamlined)
         { 
            indsEnd <- indsStt + ncXR - 1 ; inds <- indsStt:indsEnd
            B.q.SigmaR <- (B.q.SigmaR 
                           + Sigma.q.betau[inds,inds] 
                           + tcrossprod(mu.q.betau[inds]))
            indsStt <- indsStt + ncXR
         }

         if (doStreamlined)
            B.q.SigmaR <- B.q.SigmaR + Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
      }
      M.q.inv.SigmaR <- (nuVal + m + ncXR - 1)*solve(B.q.SigmaR)

      ## Update q*(a.u) parameters:

      B.q.a.u <- mu.q.recip.sigsq.u + (1/A.u^2)
      mu.q.recip.a.u <- 1/B.q.a.u

      ## Update q*(sigsq.u) parameters:
      
      indsStt <- ncX+1; A.q.sigsq.u <- NULL; B.q.sigsq.u <- NULL
      for (ell in 1:L)
      {
         if (!doStreamlined) 
         {
            mu.q.betauG <- mu.q.betau[1:ncCG]
            Sigma.q.betauG <- Sigma.q.betau[1:ncCG,1:ncCG]
         }
           
         indsEnd <- indsStt + ncZG[ell] - 1; inds <- indsStt:indsEnd
         A.q.sigsq.u[ell] <- ncZG[ell] + 1
         B.q.sigsq.u[ell] <- (2*mu.q.recip.a.u[ell] + sum(mu.q.betauG[inds]^2)
                              + sum(diag(Sigma.q.betauG[inds,inds])))
         indsStt <- indsEnd + 1                                              
      }  
      mu.q.recip.sigsq.u <- A.q.sigsq.u/B.q.sigsq.u
      
      ## Obtain current log(ML):
   
      logMLcurr <- logML(nuVal,ncXR,ncX,ncZG,m,numObs,L,yAdj,C,
                     sigsq.beta,fitVecBetauG,fitVecuR,wtVec,phiVec,
                     xiVec,xiSqd,mu.q.betauG,Sigma.q.betauG,det.Sigma.q.betau,
                     B.q.SigmaR,B.q.sigsq.u,A.R,M.q.inv.SigmaR,
                     mu.q.recip.a.R,B.q.a.R,A.u,B.q.a.u,mu.q.recip.a.u,
                     mu.q.recip.sigsq.u,A.eps,B.q.sigsq.eps,B.q.a.eps,
                     mu.q.recip.a.eps,mu.q.recip.sigsq.eps)
      logMLgrid[itnum] <- logMLcurr[[length(logMLcurr)]][1]

      ## Compute relative error:
   
      relErr <- abs((logMLcurr[[length(logMLcurr)]][1]/logMLprev)-1)

      ## Check `converged' conditions:
                    
      if (itnum >= maxIter) 
      {
         converged <- TRUE
         print("WARNING: maximum number of iterations exceeded.")
      }
   
      if (relErr < tolerance) converged <- TRUE 

      cat(itnum,sep="\n")
      print(cbind(itnum,(logMLcurr[[length(logMLcurr)]][1]-logMLprev)))
   }

   if (!doStreamlined)
   {
      if (responseType=="Gaussian")
         return(list(mu.q.betauG=mu.q.betau[1:ncCG],
                     Sigma.q.betauG=Sigma.q.betau[1:ncCG,1:ncCG],
                     A.q.sigsq.eps=A.q.sigsq.eps,B.q.sigsq.eps=B.q.sigsq.eps,
                     B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid))
      
      if (responseType=="Bernoulli")
         return(list(mu.q.betauG=mu.q.betau[1:ncCG],
                     Sigma.q.betauG=Sigma.q.betauG[1:ncCG,1:ncCG],
                     B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid))
   }
   
   if (doStreamlined)
   {
      if (responseType=="Gaussian")
         return(list(mu.q.betauG=mu.q.betauG,Sigma.q.betauG=Sigma.q.betauG,
                     A.q.sigsq.eps=A.q.sigsq.eps,B.q.sigsq.eps=B.q.sigsq.eps,
                     B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid,
                     mu.q.recip.sigsq.eps=mu.q.recip.sigsq.eps,
                     mu.q.recip.a.eps=mu.q.recip.a.eps,
                     M.q.inv.SigmaR=M.q.inv.SigmaR,
                     mu.q.recip.aR=mu.q.recip.aR,
                     mu.q.recip.sigsq.u=mu.q.recip.sigsq.u))

      if (responseType=="Bernoulli")
         return(list(mu.q.betauG=mu.q.betauG,Sigma.q.betauG=Sigma.q.betauG,
                     B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid,
                     M.q.inv.SigmaR=M.q.inv.SigmaR,
                     mu.q.recip.aR=mu.q.recip.aR,
                     mu.q.recip.sigsq.u=mu.q.recip.sigsq.u))
   }
}   

######### End of mlevMFVB.r ##########


