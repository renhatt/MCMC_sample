library(numDeriv)

ARMH <-function(vx, vhyper, cRep, cburn)
{
  
  # Sample size
  cn  = NROW(vx)
  
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
  
  # Function to compute log posterior density
  fLogPost <-function(dpar)
  {
    dlpost = -0.5*(dpar - dmu0)^2/ds02 -
              0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dpar)^2/dsigma2) ) 
    return(dlpost)    
  }
  
  # Maximize log posterior
  result <- optim (0.0, fLogPost, method = "BFGS", 
                   control = list(fnscale = -1),  hessian =TRUE)

  # mu_hat, gradient, hessian (gradient = 0 since muhat is the maximizer)
  dmuhat = result$par
  # dhess  = result$hessian
  dlhat  = result$value
  
  # We can use numDeiv library to compute gradient and hessian
  dgrad  = grad(fLogPost, dmuhat)
  dhess  = hessian(fLogPost, dmuhat)   

  # Set initial values for mu
  dmu = dmuhat 
  
  # the proposal for AR step : N(dm, dv)
  dv = -1/dhess 
  dm = dmuhat + dv * dgrad  # dgrad = l'(mu_hat) = 0
    
  # Current value of mu
  dmu_c  = dmu
 
   
  # Compute the log posterior and proposal densites with dmu_c
  dlf_c = -0.5*(dmu_c - dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_c)^2/dsigma2) )
  dlh_c = dlhat  + dgrad * (dmu_c - dm) - 0.5*(dmu_c - dm)^2/dv
  
  # counters for the acceptance in AR step and MH step
  cARrep = 0
  cMHacc = 0
  
  for(i in 1:cN)
  {
    
    # AR step
    repeat
    {
      # Count the number of repetitions
      
      if( i > cburn ) cARrep = cARrep + 1

      # Propose mu
      dmu_p = rnorm( 1, mean = dm, sd = sqrt(dv) )
    
      # Compute the log posterior
      dlf_p = -0.5*(dmu_p- dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_p)^2/dsigma2) )
      dlh_p = dlhat   + dgrad * (dmu_p - dm) - 0.5*(dmu_p -dm)^2/dv
      dratio = exp( dlf_p - dlf_c -dlh_p + dlh_c )

      if( runif(1) < dratio  )  break
    }
    
    # MH step
    dratio =  exp( dlf_p - dlf_c - min(dlf_p, dlh_p)  + min(dlf_c, dlh_c) )
    
    if( runif(1) < dratio  ) 
    { 
      dmu_c = dmu_p;  
      if( i > cburn ) cMHacc = cMHacc + 1
    }
    
    dlf_c = -0.5*(dmu_c-dmu0)^2/ds02 - 0.5*(dnu + 1)*sum( log( 1 + 1/dnu*(vx - dmu_c)^2/dsigma2) )
    dlh_c = dlhat   + dgrad * (dmu_c - dm) -0.5*(dmu_c - dm)^2/dv
    
    # Save MCMC samples
    if(i > cburn)
    {
      ck           = i - cburn
      vmus[ck]     = dmu_c
    }
  }
  cat("Acceptance Rate in AR step is", cRep/cARrep*100, "%\n")   
  cat("Acceptance Rate in MH step is", cMHacc/cRep*100, "%\n")   
  
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
# priors mu ~ N(dmu0, ds02)
# Set hyperpameters (dmu0, ds02)
vhyper = c(0.0, 1000.0)

# Set number of iterations and burn-in periods
cRep = 10000; cburn = 0.1*cRep;

vmus <- ARMH(vx, vhyper, cRep, cburn)

# Summary Statistics
myname = expression(mu)
source('ReportMCMC.R')
ReportMCMC(vmus, vname = myname, soutfilename="EX11_R")



