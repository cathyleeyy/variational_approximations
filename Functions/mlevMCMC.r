## ----------------------------------------------------------------------------
## Name : mlevMCMC.r
## ----------------------------------------------------------------------------
## Authors: Cathy Y. Y. Lee 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Obtain Markov chain Monte Carlo samples via Stan
## ----------------------------------------------------------------------------
## Required Packages: rstan
## ----------------------------------------------------------------------------
## Usage: mlevMCMC(responseTypeVal,compileCode,useSavedFit)
## 
##        responseTypeVal: Type of response distribution (Gaussian or Bernoulli)
##
##        complileCode:    TRUE or FALSE. If the same model is going to be used
##                         repeatedly, it is better to compile it just once. 
##
##        useSavedFit :    TRUE or FALSE. The Stan objects are saved in a file
##                         and reloaded for later use.   
## ----------------------------------------------------------------------------

mlevMCMC <- function(responseTypeVal,compileCode,useSavedFit)
{
   debugMode <- T
   
   if (responseTypeVal == "Gaussian")
   {    
      twoLevModel <-                  
      'data                               
      {
         int <lower=1> n;          int <lower=1> m;
         int <lower=1> numKnots;   int <lower=1> idnum[n];
         int <lower=1> ncX;        matrix[n,ncX] X;          
         vector[n] y;              matrix[n,numKnots] ZG;    
         real x[n];                real <lower=0> sigmaBeta; 
         real <lower=1> Aeps;      real <lower=1> Au;        
         real <lower=1> Ar;        vector[2] zeroVec;           
      }
      parameters                   
      {
         vector[ncX] beta;         matrix[2,m] uR;
         vector[numKnots] uG;      cov_matrix[2] SigmaR;  
         real <lower=0> sigmauG;   real <lower=0> sigmaEps;
         vector <lower=1e-10> [2] aR; 
      }
      transformed parameters
      {
         vector[n] fmean;          real mu[n];
         fmean <- X*beta + ZG*uG;
         for (i in 1:n)                                                  
         {             
            mu[i] <- (fmean[i] + uR[1,idnum[i]] + uR[2,idnum[i]]*x[i]);                         
         }
      }
      model 
      {
         matrix[2,2] rateForWish;
      
         y ~ normal(mu,sigmaEps);
      
         for (j in 1:m)
         {
            col(uR,j) ~ multi_normal(zeroVec,SigmaR);        
         }
      
         rateForWish[1,2] <- 0 ; rateForWish[2,1] <- 0 ; 
         for (r in 1:2)
         {
            aR[r] ~ inv_gamma(0.5,pow(Ar,-2));     
            rateForWish[r,r] <- 4/aR[r];
         }  
         SigmaR ~ inv_wishart(3,rateForWish);
      
         uG ~ normal(0,sigmauG);
         beta ~ normal(0,sigmaBeta);
         sigmaEps ~ cauchy(0,Aeps); 
         sigmauG ~ cauchy(0,Au);
      }'
      
      # Set up input data:
      
      allData <- list(n=numObs,m=m,X=X,y=y,ZG=ZG,ncX=ncX,
                      idnum=idnum,zeroVec=rep(0,2),x=x1,
                      numKnots=ncol(ZG),sigmaBeta=sqrt(sigsq.beta),
                      Aeps=A.eps,Au=A.u,Ar=A.R)
   }

   if (responseTypeVal == "Bernoulli")
   {
      twoLevModel <-                  
      'data                               
      {
         int <lower=1> n;           int <lower=1> m;
         int <lower=1> numKnots;    int <lower=1> idnum[n];
         int <lower=1> ncX;         matrix[n,ncX] X;            
         int<lower=0,upper=1> y[n]; matrix[n,numKnots] ZG;    
         real x[n];                 real <lower=0> sigmaBeta; 
         real <lower=1> Au;         real <lower=1> Ar;        
         vector[2] zeroVec;        
      }
      parameters                   
      {
         vector[ncX] beta;          matrix[2,m] uR;
         vector[numKnots] uG;       cov_matrix[2] SigmaR;  
         real <lower=0> sigmauG;    vector <lower=1e-10> [2] aR; 
      }
      transformed parameters
      {
         vector[n] fmean;           real mu[n];
         fmean <- X*beta + ZG*uG;
         for (i in 1:n)                                                  
         {             
            mu[i] <- fmean[i] + uR[1,idnum[i]] + uR[2,idnum[i]]*x[i];                         
         }
      }
      model 
      {
         matrix[2,2] rateForWish;
      
         y ~ bernoulli_logit(mu);
      
         for (j in 1:m)
         {
            col(uR,j) ~ multi_normal(zeroVec,SigmaR);        
         }
      
         rateForWish[1,2] <- 0 ; rateForWish[2,1] <- 0 ; 
         for (r in 1:2)
         {
            aR[r] ~ inv_gamma(0.5,pow(Ar,-2));     
            rateForWish[r,r] <- 4/aR[r];
         }  
         SigmaR ~ inv_wishart(3,rateForWish);
      
         uG ~ normal(0,sigmauG);
         beta ~ normal(0,sigmaBeta);
         sigmauG ~ cauchy(0,Au);
      }'
      
      # Set up input data:
      
      allData <- list(n=numObs,m=m,X=X,y=y,ZG=ZG,ncX=ncX,
                      idnum=idnum,zeroVec=rep(0,2),x=x1,
                      numKnots=ncol(ZG),sigmaBeta=sqrt(sigsq.beta),
                      Au=A.u,Ar=A.R)
   }

   filenameStan <- paste(getwd(),"/Data/StanObj_",responseTypeVal,isim,".txt",sep="")

   ## Compile code for model if required:

   if (compileCode)
      StanCompilObj <- stan(model_code=twoLevModel,data=allData,iter=2,chains=1)
   
   if (!useSavedFit)
   {
      cat("About to enter stan()... \n")
      timeInfo <- system.time(
      STANobj <- stan(model_code=twoLevModel,data=allData,warmup=nBurnin,
                      iter=(nBurnin+nIter),chains=1,thin=nThin,fit=StanCompilObj,
                      control=list(stepsize=0.001,adapt_delta=0.99))
      )
      timeMCMC <- timeInfo[3]       
      save(STANobj,file=filenameStan)
      cat("Finished with stan().\n")
   }

   ## Load the Stan object:

   if (useSavedFit)
   {
      load(filenameStan)
      timeMCMC <- NA
   }

   return(list(StanObj=STANobj,timeMCMC=timeMCMC))
}

########## End of mlevMCMC.r ##########
