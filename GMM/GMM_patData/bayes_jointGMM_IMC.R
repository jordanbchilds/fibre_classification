
library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
} else {
  imc_chan = args
}

cramp = colorRamp(c(rgb(0,0,1,0.20),rgb(1,0,0,0.2)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

myDarkGrey = rgb(169,169,159, max=255, alpha=100)
myGreen = rgb(25,90,0,max=255,alpha=50)
myYellow = rgb(225,200,50,max=255, alpha=50)
myBlue = cramp(1)
myRed = cramp(0)

classcols = function(classif){
  # A function using the cramp specified colours to generate rgb colour names
  # input: a number in [0,1]
  # output: rgb colour name
  rgbvals = cramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3],rgbvals[,4]))
}

priorpost = function(data, prior, posterior, classifs, ctrl=NULL,
                     title){
  # output: plots the prior and posterior regression lines and data
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  if( is.null(ctrl) ){
    x.lim = range(data[,1]) + c(-1,1)
    y.lim = range(data[,2]) + c(-1,1)
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    contour( kde2d(prior[,'compOne[1]'], prior[,'compOne[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    contour( kde2d(posterior[,'compOne[1]'], posterior[,'compOne[2]'], n=100), add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
  } else {
    x.lim = range( data[,1], ctrl[,1] ) + c(-1,1)
    y.lim = range( data[,2], ctrl[,2] ) + c(-1,1)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=myYellow, pch=20)
    contour( kde2d(prior[,'compOne[1]'], prior[,'compOne[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    
    points( data[,1], data[,2], col=classcols(classifs), pch=20)
    contour( kde2d(posterior[,'compOne[1]'], posterior[,'compOne[2]'], n=100), add=TRUE, nlevels=5 )
    title(main=title, line = -1, outer = TRUE)
  }
  
  par(op)
} 

priorpost_marginals = function(prior, posterior, data=NULL, title){
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  par(mfrow=c(2,2))
  ## mu_1
  # prior
  contour( kde2d(prior[,'mu[1,1]'], prior[,'mu[1,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[1,1]'], posterior[,'mu[1,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Posterior Density') )
  ## mu_2
  # prior
  contour( kde2d(prior[,'mu[2,1]'], prior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[2,1]'], posterior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Posterior Density') )
  title(main=title, line = -1, outer = TRUE)
  
  
  par(mfrow=c(2,3))
  ## tau_1
  # prior
  contour( kde2d( prior[,'tau[1,1,1]'], prior[,'tau[2,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[111]), ylab=expression(tau[221]), nlevels=5,
           main=expression(tau[1]~'Prior Density') )
  
  contour( kde2d( prior[,'tau[1,1,1]'], prior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[111]), ylab=expression(tau[121]), nlevels=5,
           main=expression(tau[1]~'Prior Density') )
  
  contour( kde2d( prior[,'tau[2,2,1]'], prior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[221]), ylab=expression(tau[121]), nlevels=5,
           main=expression(tau[1]~'Prior Density') )
  ## tau_1
  # posterior
  contour( kde2d( posterior[,'tau[1,1,1]'], posterior[,'tau[2,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[111]), ylab=expression(tau[221]), nlevels=5,
           main=expression(tau[1]~'Posterior Density') )
  
  contour( kde2d( posterior[,'tau[1,1,1]'], posterior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[111]), ylab=expression(tau[121]), nlevels=5,
           main=expression(tau[1]~'Posterior Density') )
  
  contour( kde2d( posterior[,'tau[2,2,1]'], posterior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[221]), ylab=expression(tau[121]), nlevels=5,
           main=expression(tau[1]~'Posterior Density') )
  title(main=title, line = -1, outer = TRUE)
  
  ## tau_2
  # prior
  contour( kde2d( prior[,'tau[1,1,2]'], prior[,'tau[2,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[112]), ylab=expression(tau[222]), nlevels=5,
           main=expression(tau[2]~'Prior Density') )
  
  contour( kde2d( prior[,'tau[1,1,2]'], prior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[112]), ylab=expression(tau[122]), nlevels=5,
           main=expression(tau[2]~'Prior Density') )
  
  contour( kde2d( prior[,'tau[2,2,2]'], prior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[222]), ylab=expression(tau[122]), nlevels=5,
           main=expression(tau[2]~'Prior Density') )
  ## tau_2
  # posterior
  contour( kde2d( posterior[,'tau[1,1,2]'], posterior[,'tau[2,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[112]), ylab=expression(tau[222]), nlevels=5,
           main=expression(tau[2]~'Posterior Density') )
  
  contour( kde2d( posterior[,'tau[1,1,2]'], posterior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[112]), ylab=expression(tau[122]), nlevels=5,
           main=expression(tau[2]~'Posterior Density') )
  
  contour( kde2d( posterior[,'tau[2,2,2]'], posterior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(tau[222]), ylab=expression(tau[122]), nlevels=5,
           main=expression(tau[2]~'Posterior Density') )
  title(main=title, line = -1, outer = TRUE)
  
  if( !is.null(data) ){
    par(mfrow=c(1,2))
    plot( density(posterior[,'probdiff']), cex.lab=2, cex.axis=1.5, xlim=c(0,1),
          xlab='probdiff', ylab='density', lwd=2, col='red', main='probdiff Density')
    lines( density(rbeta(5000, data$alpha, data$beta)), lwd=2, col='green')
    title(main=title, line = -1, outer = TRUE)
  }
  par(op)
}

MCMCplot = function( MCMCoutput, lag=20, title ){
  col.names = colnames(MCMCoutput[[1]])
  n.chains = length(MCMCoutput)
  
  par(mfrow=c(3,3), mar = c(5.5,5.5,3,3))
  if( n.chains==1 ){
    for(param in col.names){
      plot( 0:lag, autocorr(MCMCoutput[[1]][,param], lags=0:lag), 
            type='h', ylim=c(-1,1), xlab='Index', ylab='')
      plot( 1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type='l',
            ylab='', xlab='Iteration')
      plot(density(MCMCoutput[[1]][,param] ), main='', xlab='')
    }
  } else {
    for( param in col.names){
      plot( autocorr(MCMCoutput[[1]][,param], lags=0:lag), type='h', 
            xlab='Index', ylab='' )
      for(j in 2:n.chains) lines( autocorr(MCMCoutput[[j]][,param], lags=0:20), type='h', col=j)
      
      plot(1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type='l',
           ylab='', xlab='Iteration')
      for(j in 2:n.chains) lines(MCMCoutput[[j]][,param], type='l', col=j)
      
      plot(density(MCMCoutput[[1]][,param]) )
      for(j in 2:n.chains) lines(density(MCMCoutput[[j]][,param]), col=j )
    }
  }
  title(main=title, line = -1, outer = TRUE)
}

modelstring = "
model {
  for(i in 1:Nctrl ){ # fit to ctrl data
    Yctrl[i,] ~ dmnorm(mu[,1], tau[,,1] )
  }
  for(j in 1:Npat ){ # fit to patient data
    z[j] ~ dbern(probdiff)
    class[j] =   2 - z[j]
    Ypat[j,] ~ dmnorm(mu[,class[j]], tau[,,class[j]] )
  }
  
  # covariance matrix for component 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  # covariance matrix for component 2
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  probdiff ~ dbeta(alpha, beta)

  compOne ~ dmnorm(mu[,1], tau[,,1])
  compTwo~ dmnorm(mu[,2], tau[,,2])
}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
# dir.create(file.path("PNG"), showWarnings = FALSE)
dir.create(file.path("Time"), showWarnings = FALSE)

dir.create(file.path("Output/IMC_joint"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_joint"), showWarnings = FALSE)
# dir.create(file.path("PNG/IMC_joint"), showWarnings = FALSE)

dir.create(file.path("PDF/IMC_joint/MCMC"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_joint/classifs"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_joint/marginals"), showWarnings = FALSE)


# burn-in, chain length, thinning lag
MCMCUpdates = 2000
MCMCUpdates_Report = 5000
MCMCUpdates_Thin = 1
n.chains = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

mitochan = "VDAC1"

# removing unwanted info 
imcDat = imc_data[imc_data$channel %in% c(imc_chan, mitochan), ]

froot = gsub('.RAW.txt', '', fulldat)

# getting the ranges of the axis

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P", sbj, value = TRUE)

for( chan in imc_chan ){
  for( pat in pts){
    outroot = paste( froot, pat, chan, sep='__')
    posterior_file = file.path("Output/IMC_joint", paste0(outroot, "__POSTERIOR.txt") )
  
    if( !file.exists(posterior_file)){
    
      ## CONTROL DATA
      control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
      Xctrl = log(control$value[control$channel==mitochan])
      Yctrl = log(control$value[control$channel==chan])
      XY_ctrl = cbind( Xctrl, Yctrl )
      Nctrl = nrow(XY_ctrl)
      
      ## PATIENT DATA
      patient = imcDat[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"), ] 
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      XY_pat = cbind(Xpat, Ypat)
      Npat = nrow(XY_pat)
      
      ## PRIORS
      mu1_mean = c(1,1.5)
      mu2_mean = 1.5*mu1_mean
      mu1_prec = solve( matrix(c(0.2,0.1,0.1,0.2), ncol=2, nrow=2, byrow=TRUE) )
      mu2_prec = solve( 5*diag(2) )
      
      U_1 = matrix( c(10,7,7,10), ncol=2, nrow=2, byrow=TRUE)
      n_1 = 50
      U_2 = 3*diag(2)
      n_2 = 20
      
      alpha = 1
      beta = 1
      
      data = list(Yctrl=XY_ctrl, Nctrl=Nctrl, Ypat=XY_pat, Npat=Npat,
                  mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                  mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                  U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)

      data_priorpred = data
      data_priorpred$Yctrl = NULL
      data_priorpred$Ypat = NULL
      data_priorpred$Nctrl = 0
      data_priorpred$Npat = 0 
      
      model = jags.model(textConnection(modelstring), data=data, n.chains=n.chains) 
      
      model_priorpred = jags.model(textConnection(modelstring), data=data_priorpred) 
      
      update(model, n.iter=MCMCUpdates)
      
      converge = coda.samples(model=model, variable.names=c("mu","tau","z","probdiff", "compOne", "compTwo"),
                                  n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output = coda.samples(model=model, variable.names=c("mu","tau","z","probdiff", "compOne", "compTwo"),
                                n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output_priorpred = coda.samples(model=model_priorpred,
                                          variable.names=c("mu", "tau","z","probdiff", "compOne", "compTwo"),
                                          n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      MCMCoutput = output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                             "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                             "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                             "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]",
                             "compTwo[2]")]
      
      posterior = as.data.frame(output[[1]])
      prior = as.data.frame(output_priorpred[[1]])
      
      classifs = colMeans( posterior[, grepl('z', colnames(posterior))] )
      colnames(posterior) = colnames(output[[1]])
      colnames(prior) = colnames(output_priorpred[[1]])
      
      if( pat=='CTRL'){
        pdf(file.path("PDF/IMC_joint/classifs", paste0(paste(outroot, pat, sep='__'), ".pdf")), width=14,height=8.5)
        priorpost( data=data$Yctrl, prior=prior, posterior=posterior,
                   classifs=classifs, title=paste(froot, pat, chan, sep='__'))
        dev.off()
        pdf(file.path("PDF/IMC_joint/marginals", paste0(paste(outroot, pat, sep='__'), ".pdf")), width=14, height=8.5)
        priorpost_marginals(prior=prior, posterior=posterior, 
                            title=paste(froot, pat , chan, sep='__'))
        dev.off()
      } else { 
        pdf(file.path("PDF/IMC_joint/classifs", paste0(paste(outroot, pat, sep="__"), ".pdf")), width=14,height=8.5)
        priorpost( data=data$Ypat, prior=prior, posterior=posterior, ctrl=data$Yctrl,
                   classifs=classifs, title=paste(froot, pat, chan, sep='__'))
        dev.off()
        pdf(file.path("PDF/IMC_joint/marginals", paste0(paste(outroot, pat ,sep='__'), ".pdf")), width=14, height=8.5)
        priorpost_marginals(prior=prior, posterior=posterior, data=data,
                            title=paste(froot, pat, chan, sep='__'))
        dev.off()
      }
      write.table(as.numeric(classifs),file.path("Output/IMC_joint",paste0(outroot,"__CLASS.txt")),
                  row.names=FALSE,quote=FALSE,col.names=FALSE)
      
      pdf(file.path("PDF/IMC_joint/MCMC", paste0(paste(outroot, pat, sep='__'), '.pdf')), width=14, height=8.5)
      MCMCplot( MCMCoutput, title=paste(froot, pat, chan, sep='__'))
      dev.off()
      
      write.table(posterior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]",
                               "compTwo[2]")],
                  posterior_file, row.names=FALSE, quote=FALSE)
    }
  }
}








