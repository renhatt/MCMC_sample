ARsampling <-function()
{
  # Generate n samples by AR sampling
  cn    = 10000
  # Number of iterations
  citer = 0
  
  # Store samples
  vx = rep(0, cn)
  
  dc = sqrt(2*pi)*3^(3/2)*exp(-3/2)
  
  i = 1
  while(i <= cn)
  {
    dx     = rchisq(1, 1)
    dratio = dgamma( dx, 2, scale = 1 ) / (dc * dgamma( dx, 1/2, scale = 2 ) )
    if(runif(1) < dratio ) {  vx[i] = dx; i = i + 1;}
    citer = citer + 1

  }
  cat("Acceptance Rate is ", cn/citer*100, "%\n")   
  return(vx)
}

vx <- ARsampling()

# Estimated density of Gamma(2, 1)
plot(density(vx, from=0), xlab="x", main="", col = "blue", lwd = 2, lty = 1)
curve(dgamma(x, 2, scale = 1), add=TRUE, lwd = 1, lty =3)
legend("topright", legend = c("MCMC", "True"), col = c("blue", "black"), 
       lwd = c( 2, 1), lty = c( 1, 3))


