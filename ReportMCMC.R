#' Title
#' Report Summary statistics and plots for MCMC outputs
#'
#' @param mx
#' MCMC output
#' n x K matrix where n: number of samples, K: number of parameters
#'
#' @param dBm
#' Bandwidth to compute the variance of the sample mean
#' Also used for the lag of sample autocorrelation functions
#'
#' @param vname
#' labels for parameters
#'
#' @param soutfilename
#' filenames for output files
#'
#' @return
#' @export
#'
#' @examples
#' ReportMCMC(rnorm(1500))
#'
ReportMCMC <-function(mx, dBm=NULL, vname=NULL, soutfilename=NULL)
{

  cRep = NROW(mx)
  cdim = NCOL(mx)

  if(is.null(dBm)){dBm = 2*( floor( sqrt(cRep) ) +1 );}
  if(is.null(vname))
    {
      vname = "Param1"
      if(cdim > 1)
      {
        for(i in 2:cdim){vname = cbind(vname, paste("Param",i,sep=""));}
      }
    }
  if(is.null(soutfilename))
    {
      soutfilename = "Output"
    }

  library(freqdom) # spectral.density
  library(coda)    # autocorr.plot

  # Compute the estimate for the variance of sample mean of time series
  fTsvar <-function(vx,dbm)
    {
      dsp = as.double(spectral.density(vx,freq=c(0),q=dbm,weights="Parzen")$operators[1,,])
      # dsp = spec.ar(vx,plot=FALSE)$spec[1]
      # dsp = spectrum0.ar(vx)$spec      # using coda library
      # dsp = (mcse(vx)$se)^2*length(vx) # using mcmcse library
      return(dsp)
    }
  # Function to compute the p-value of Convergence Diagnostics by comparing
  # two means (for the first 10% and the last 50% of the series)
  # H0:two means are equal (convergence)
  fCD <-function(vx, dbm)
  {
    cT  = length(vx);  cn1 = floor(0.1*cT); cn2 = floor(0.5*cT);
    vx1 = vx[1:cn1]; vx2 = vx[cn2:cT];

    dx1bar = mean(vx1); dx2bar = mean(vx2);
    dvar1=fTsvar(vx1, floor( 2*sqrt(cn1) ) +1 )
    dvar2=fTsvar(vx2, floor( 2*sqrt(cn2) ) +1 )
    dz=(dx1bar - dx2bar)/sqrt(dvar1/cn1+dvar2/cn2)

    return(2*pnorm(abs(dz), lower.tail=FALSE));
  }

  ################################################################
  # Estimation results
  ################################################################
  result1        = matrix(0, cdim, 5)
  result2        = matrix(0, cdim, 4)
  vx             = matrix(0, cRep, 1)
  for(i in 1:cdim)
    {

     if( cdim==1 ){ vx = mx;}
     else         { vx = mx[,i];}

      # Compute Inefficiency Factor
      dIF           = fTsvar(vx, dBm) / var(vx)
      if(dIF <  1) {dIF = 1; }
      if(dIF < 10) {dIF.rounded = round(dIF, digits=1); }
      else         {dIF.rounded = round(dIF, digits=0); }
      # Summary Statistics
      result1[i,1]   = mean(vx)
      result1[i,2]   = sd(vx)
      result1[i,3:5] = quantile(vx, c(0.025, 0.5, 0.975))
      result2[i,1]   = round(cRep/dIF, digits=0) # ESS
      result2[i,2]   = dIF.rounded
      result2[i,3]   = round(fCD(vx,dBm),digits=3)
      result2[i,4]   = round(sum(vx>0)/cRep, digits=4)
          }
  colnames(result1) = c("Mean", "Std Dev", "95%L", "Median", "95%U")
  rownames(result1) = vname
  colnames(result2) = c("ESS", "IF", "CD", "Pr(+)")
  rownames(result2) = vname
  fname = paste(soutfilename,"_Result_it_",cRep,".txt",sep="")
  # Output to the file
  sink(fname);
  print(result1, digit=5); print(result2, digit=5);
  sink();
  # Output on screen
  print(result1, digit=5); print(result2, digit=5);

  if(cdim == 1){vdisp=c(1,1);}
  else if(cdim == 2){vdisp=c(2,1);}
  else if(cdim == 3){vdisp=c(3,1);}
  else if(cdim == 4){vdisp=c(2,2);}
  else              {vdisp=c(4,2);}

  defaultPar <- par(no.readonly = TRUE)
  ################################################################
  # Sample path
  ################################################################
  fname = paste(soutfilename,"_SamplePath_it_",cRep,".pdf",sep="")
  pdf(file=fname)
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
   {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    # plot thinned MCMC samples
    vx.thinned = vx[seq(1, length(vx), by = max(1,floor(length(vx)/1000)))]
    plot.ts(vx.thinned, xlab="Iteration", ylab=vname[i], main="")
  }
  dev.off()

  # Output on screen
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
  {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    # plot thinned MCMC samples
    vx.ind     = seq(1, length(vx), by = max(1,floor(length(vx)/1000)))
    vx.thinned = vx[vx.ind]
    plot(vx.ind, vx.thinned, xlab="Iteration", ylab=vname[i], main="", type="l")
  }

  ################################################################
  # Sample autocorrelation function
  ################################################################
  fname = paste(soutfilename,"_ACF_it_",cRep,".pdf",sep="")
  pdf(file=fname)
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
    {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    autocorr.plot(vx, lag.max=floor(length(vx)*0.05), auto.layout=FALSE, main=vname[i])
    }
  dev.off()

  # Output on screen
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
  {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    autocorr.plot(vx, lag.max=floor(length(vx)*0.05), auto.layout=FALSE, main=vname[i])
  }

  ################################################################
  # Estimated posterior density
  ################################################################
  fname = paste(soutfilename,"_Density_it_",cRep,".pdf",sep="")
  pdf(file=fname)
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
  {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    plot(density(vx), xlab=vname[i], main ="")
  }
  dev.off()
  # Output on screen
  par(mfrow=vdisp, plt = c(0.2,0.8,0.4,0.8))
  for(i in 1:cdim)
  {
    if( cdim==1 ){ vx = mx;}
    else         { vx = mx[,i];}
    plot(density(vx), xlab=vname[i], main ="")
  }
  par(defaultPar)
}
