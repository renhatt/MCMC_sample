IndepMH <-function(vx, dini, vhyper, cRep, cburn)
{
  
  # Sample size
  cn  = NROW(vx)
  
  # Set initial values for mu
  dmu = dini 
  
  # Set hyperparameter for mu
  dmu0  = vhyper[1]
  ds02  = vhyper[2]
  
  # Known parameters
  dnu     = 5.0
  dsigma2 = 4.0
  
  # Set the total number of MCMC iterations
  cN = cRep + cburn
  
  # Save MCMC samples for mu
  vmus   = matrix(0, cRep, 1)
  
  # Current value of mu
  dmu_c  = dmu
  
  # Compute mean and std dev of the proposal
  dmu1  = mean( vx )
  dsd1  = sqrt( dnu * dsigma2 / cn )
  
  # Compute the log posterior and proposal densites with dmu_c
  dlf_c = -0.5*(dmu_c - dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_c)^2/dsigma2) )
  dlh_c = -0.5*(dmu_c - dmu1)^2/dsd1^2
  
  # counter for acceptance in MH algorithm
  cacc = 0
  
  for(i in 1:cN)
  {
    # Propose mu
    dmu_p = rnorm( 1, mean = dmu1, sd = dsd1 )
    
    # Compute the log posterior
    dlf_p = -0.5*(dmu_p-dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_p)^2/dsigma2) )
    dlh_p = -0.5*(dmu_p - dmu1)^2/dsd1^2
    
    
    dratio = exp( dlf_p - dlf_c -dlh_p + dlh_c )
    if( runif(1) < dratio  ) 
    { 
      dmu_c = dmu_p;  
      if( i > cburn ) {cacc = cacc + 1;} 
    }
    
    dlf_c = -0.5*(dmu_c-dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_c)^2/dsigma2) )
    dlh_c = -0.5*(dmu_c - dmu1)^2/dsd1^2
    
    # Save MCMC samples
    if(i > cburn)
    {
      ck           = i - cburn
      vmus[ck]     = dmu_c
    }
  }
  cat("Acceptance Rate is", cacc/cRep*100, "%\n")   
  
  return(vmus)
}



# Generate data from t distribution
# X_i ~ i.i.d T(mu, dsigma2, dnu), i=1,...,cn 
# where E(X_i)= dmu, Var(X_i) = dnu/(dnu-2)*dsigma2
# and dnu and dsigma2 are assumed to be known.
#
dnu = 5.0; dsigma2 = 4.0; 
dmu = 5.0; cn = 100;
#
vx <- dmu + sqrt(dsigma2)*rt(cn, dnu)

# MCMC
# Set initial pameters mu
dini = 0.0
# priors mu ~ N(dmu0, ds02)
# Set hyperpameters (dmu0, ds02)
vhyper = c(0.0, 1000.0)

# Set number of iterations and burn-in periods
cRep = 10000; cburn = 0.1*cRep;

vmus <- IndepMH(vx, dini, vhyper, cRep, cburn)

# Summary Statistics
myname = expression(mu)
source('ReportMCMC.R')
ReportMCMC(vmus, vname = myname, soutfilename="EX10_R")



