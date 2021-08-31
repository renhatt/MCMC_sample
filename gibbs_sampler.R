RNormalGamma1 <-function(vx, vinit, vhyper, cRep, cburn)
{

  # Sample size
  cn  = NROW(vx)
  
  # Set initial values for mu & sigma2
  dmu = vinit[1]; dsigma2 = vinit[2]

  # Set hyperparameter for mu
  dmu0 = vhyper[1];
  dm0  = vhyper[2];
  
  # Set hyperparameter for sigma2
  dalpha0 = vhyper[3]; 
  dbeta0  = vhyper[4];
  
  # Set the total number of MCMC iterations
  cN = cRep + cburn
  
  # Save MCMC samples for mu and sigma
  vmus     = matrix(0, cRep, 1)
  vsigma2s = matrix(0, cRep, 1)
  vnewxs   = matrix(0, cRep, 1)
  
  # Compute xbar, m1, alpha1, beta1
  dxbar    = mean(vx)

  dm1      = dm0 + cn
  dmu1  = ( dm0*dmu0 + cn*dxbar ) / ( dm0 + cn ) 
  dalpha1  = dalpha0 + cn
  dbeta1   = dbeta0  + sum((vx-dxbar)^2) + cn*dm0*(dxbar-dmu0)^2 / (cn+dm0)
  
  for(i in 1:cN)
  {
    # Generate mu
    dmu   = rnorm( 1, mean = dmu1, sd = sqrt(dsigma2/dm1) )
    
    # Generate sigma2
    dsigma2 = 1/rgamma(1, shape = 0.5*( dalpha1+1 ), 
                           rate = 0.5*( dm1*(dmu-dmu1)^2 + dbeta1 ) )

    # Save MCMC samples
    if(i > cburn)
    {
      ck           = i - cburn
      vmus[ck]     = dmu
      vsigma2s[ck] = dsigma2
      vnewxs[ck]   = rnorm( 1, mean = dmu, sd = sqrt(dsigma2) )
    }

  }
     return(cbind(vmus, vsigma2s, vnewxs))
}

# x_i ~ i.i.d. N(mu,sigma2) (i=1,2,...,n)
# Set true values x_i ~ N(mu,sigma2)
set.seed(123)
dmu = 5.0; dsigma2 = 1.0; cn =100;
# Generate data
vx  = dmu + sqrt(dsigma2)*rnorm(cn);
#

# MCMC
# Set initial pameters (mu, sigma^2)
vini=c(0.0, 1.0);
# priors mu ~ N(mu0, sigma^2/m0), sigma2 ~ IG(dn0/2,dS0/2)
# Set hyperpameters (mu0, m0, alpha0, beta0)
vhyper = c(0.0, 1.0, 0.02, 0.02)
# Set number of iterations and burn-in periods
cRep = 10000; cburn = 0.1*cRep;
# 
mx <- RNormalGamma1(vx, vini, vhyper, cRep, cburn)
# Summary Statistics
myname = c(expression(mu), expression(sigma^2))
source('ReportMCMC.R')
ReportMCMC(mx[,1:2], vname = myname, soutfilename="EX7_R")
# Estimated density of new x
plot(density(mx[,3]), xlab="x", main="", col = "blue", lwd = 2, lty = 1)
curve(dnorm(x, mean=mean(mx[,1]), sd = sqrt(mean(mx[,2]))), add=TRUE, lwd = 1, lty =3)
legend("topright", legend = c("MCMC", "Plugged in"), col = c("blue", "black"), 
                   lwd = c( 2, 1), lty = c( 1, 3))


