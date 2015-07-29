Streamlined mean field variational Bayes algorithms
==============

A master R script, along with five user-written R functions, are included 
to illustrate fast and memory-efficient Bayesian fitting of longitudinal 
and multilevel models with Gaussian or Bernoulli response, via mean field 
variational Bayes approximations. The master script contains a setwd() 
command and this should point to the directory in which the R function 
and output files reside. It also contains a section that requires users
to specify Boolean flags for the type of simulation studies.

Summary of R files:
-------------------

**mlevAppoxInference.R**: Master R script

**ZOSull.r**:         R function for obtaining Z matrix of O'Sullivan spline basis functions

**mlevMCMC.r**:       R function for obtaining Markov chain Monte Carlo samples 

**mlevMFVB.r**:       R function for obtaining mean field variational Bayes estimates
  
**summMCMC.r**        R function for summarizing Markov chain Monte Carlo outputs

**accVarApp.r**:      R function for assessing accuracy of a variational approximate posterior density function

**simuBoxplots.R**:   R script for obtaining boxplot summaries of variational approximation accuracy results

**simuCoverage.R**:   R script for obtaining coverage percentage summaries of simulation results

The code has been written using R-3.2.0 (2015-04-16).

Platform: x86_64-apple-darwin13.4.0 (64-bit)

Attached base packages:

[1] splines stats graphics grDevices utils datasets methods base     

Other attached packages:

[1] MCMCpack_1.3-3   coda_0.17-1   KernSmooth_2.23-14   lattice_0.20-31   
[5] gdata_2.13.3     rstan_2.6.0   inline_0.3.14        Rcpp_0.11.5       
[9] Matrix_1.2-0     magic_1.5-6   abind_1.4-0          MASS_7.3-40       

Citation
-----------

The underlying algorithms are based on the following paper (under revision in Biometrical Journal)

http://matt-wand.utsacademics.info/LeeWand.pdf

    <div class="box">
        Foo Bar Baz
    </div>
    
@article{Lee2015,

  title={Streamlined mean field variational Bayes for longitudinal and multilevel data analysis},
  
  author={Lee, Cathy Y Y and Wand, Matt P},
  
  year={2015}
}


