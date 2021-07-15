
library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

myDarkGrey = rgb(169,169,159, max=255, alpha=100)
myGreen = rgb(25,90,0,max=255,alpha=50)
myYellow = rgb(225,200,50,max=255, alpha=50)

cramp = colorRamp(c(rgb(1,0,0,0.2), rgb(0,0,1,0.20)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

classcols = function(classif){
  # A function using the cramp specified colours to generate rgb colour names
  # input: a number in [0,1]
  # output: rgb colour name
  rgbvals = cramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3],rgbvals[,4]))
}

comparedensities=function(priorvec, posteriorvec, xlab="", main="", xlim=-99){
  # output: figure of prior and posterior densities (one figure)
  d1 = density(priorvec)
  d2 = density(posteriorvec)
  xvals = c(d1$x,d2$x)
  yvals = c(d1$y,d2$y)
  if(length(xlim)==1){
    xlim = range(xvals)
  }
  plot(density(posteriorvec),lwd=3,xlim=xlim,ylim=c(0,max(yvals)),col="green",main=main,xlab=xlab,cex.lab=2,cex.axis=1.5)
  points(density(priorvec),lwd=3,col="black",type="l")
}

priorpost = function(ctrl_data, ctrl_prior, ctrl_posterior, 
                     pat_data=NULL, pat_prior=NULL, pat_posterior=NULL, 
                     class_posterior=NULL, classifs=NULL, output_mcmc=NULL, title){
  # output: plots the prior and posterior regression lines and data
  
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  if( is.null(pat_data) ){ # plot ctrl prior and posterior
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Control Prior')
    contour( kde2d(ctrl_prior[,'Y_syn[1]'], ctrl_prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Control Posterior')
    contour( kde2d(ctrl_posterior[,'Y_syn[1]'], ctrl_posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
    
    prior = ctrl_prior
    posterior = ctrl_posterior
  } else { # plot ctrl prior-posterior and patient prior-posterior
    class_posterior = colMeans( pat_posterior[,grepl("z",colnames(pat_posterior))])
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myGreen, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Contorl Prior')
    contour( kde2d(ctrl_prior[,'Y_syn[1]'], ctrl_prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myGreen, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Control Posterior')
    contour( kde2d(ctrl_posterior[,'Y_syn[1]'], ctrl_posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Patient Prior',
          ylim=(range(c(ctrl_data$Y[,2], pat_data$Y[,2]))+c(-1,1)), xlim=(range(c(ctrl_data$Y[,1], pat_data$Y[,1]))+c(-1,1)) )
    points(  pat_data$Y[,1], pat_data$Y[,2], col=myGreen,  pch=20 )
    contour( kde2d(pat_prior[,'Y_syn[1]'], pat_prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Patient Posterior',
          ylim=(range(c(ctrl_data$Y[,2], pat_data$Y[,2]))+c(-1,1)), xlim=(range(c(ctrl_data$Y[,1], pat_data$Y[,1]))+c(-1,1)) )
    points( pat_data$Y[,1], pat_data$Y[,2], col=classcols(class_posterior), pch=20 )
    contour( kde2d(pat_posterior[,'Y_syn[1]'], pat_posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
    
    title(main=title, line = -1, outer = TRUE)
    
    prior = pat_prior
    posterior = pat_posterior
  }
  
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
  
  if( !is.null(pat_data) ){
    par(mfrow=c(1,2))
    plot( density(posterior[,'probdiff']), cex.lab=2, cex.axis=1.5, xlim=c(0,1),
          xlab='probdiff', ylab='density', lwd=2, col='red', main='probdiff Density')
    lines( density(rbeta(5000,pat_data$alpha_p, pat_data$beta_p)), lwd=2, col='green')
    title(main=title, line = -1, outer = TRUE)
  }
  if( !is.null(output_mcmc) ){
    par(mfrow=c(2,3))
    plot(output_mcmc[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                        "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                        "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                        "probdiff", "Y_syn[1]", "Y_syn[2]")])
    par(mfrow=c(1,1))
  }
  
  par(op)
} 

priorpost_pred = function(Yctrl, Ypat, posterior, prior, classifs, title){
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  
  x.lim = range(c(Yctrl[,1], Ypat[,1])) + c(-1,1)
  y.lim = range(c(Yctrl[,2], Ypat[,2])) + c(-1,1)
  
  # Y_syn_prior = matrix(0, nrow=2*nrow(prior), ncol=2)
  # Y_syn_prior[,1] = c(prior[,'Y_syn_g1[1]'], prior[,'Y_syn_g2[1]'])
  # Y_syn_prior[,2] = c(prior[,'Y_syn_g1[2]'], prior[,'Y_syn_g2[2]'])
  # 
  plot(Yctrl[,1], Yctrl[,2], pch=20, col=myDarkGrey, cex.lab=2, cex.axis=1.5,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Prior Predictive',
       xlim=x.lim, ylim=y.lim)
  points(Ypat[,1], Ypat[,2], pch=20, col=myGreen)
  contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
  
  # Y_syn_post = matrix(0, nrow=2*nrow(prior), ncol=2)
  # Y_syn_post[,1] = c(posterior[,'Y_syn_g1[1]'], posterior[,'Y_syn_g2[1]'])
  # Y_syn_post[,2] = c(posterior[,'Y_syn_g1[2]'], posterior[,'Y_syn_g2[2]'])
  
  plot(Yctrl[,1], Yctrl[,2], pch=20, col=myDarkGrey, cex.lab=2, cex.axis=1.5,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Posterior Predictive',
       xlim=x.lim, ylim=y.lim)
  points( Ypat[,1], Ypat[,2], col=classcols(classifs), pch=20)
  contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
  
  title(main=title, line = -1, outer = TRUE)
  
  # plot( prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], pch=20, col=myDarkGrey)
  # contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
  # 
  # plot( posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], pch=20, col=myDarkGrey)
  # contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
  # 
  # title(main=title, line = -1, outer = TRUE)
  
  par(op)
}

modelstring = "
model {
  for(i in 1:Nctrl ){ # fit to ctrl data
    Yctrl[i,] ~ dmnorm(mu[,1], tau[,,1] )
  }
  
  for(j in 1:Npat ){ # fit to patient data
    z[j] ~ dbern(probdiff)
    class[j] =  2 - z[j]
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

  # posterior distribution
  z_syn ~ dbern(probdiff)
  class_syn = z_syn + 1
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
  
  pred_1 ~ dmnorm(mu[,1], tau[,,1])
  pred_2 ~ dmnorm(mu[,2], tau[,,2])
}
"
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_joint"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_joint"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_joint"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
MCMCUpdates = 2000
MCMCUpdates_Report = 5000
MCMCUpdates_Thin = 1
n.chains = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

colnames(imc_data)

unique(imc_data$channel)

# removing unwanted info 
imc_chan = c('SDHA','OSCP', 'VDAC1', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imc_data_1.0 = imc_data[imc_data$channel %in% imc_chan, ]

imcDat=imc_data_1.0

mitochan = "VDAC1"

froot = gsub('.RAW.txt', '', fulldat)

# getting the ranges of the axis
imc_lims = list()
for(ch in imc_chan){ 
  imc_lims[[ch]] = quantile(log(imcDat$value[imcDat$channel==ch]),
                            c(0.001,0.999), na.rm=TRUE)
}

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P", sbj, value = TRUE)
pts = c('P01',  'P07')

for( chan in imc_chan[-which(imc_chan == 'VDAC1')]){
  for( pat in pts){
    outroot = paste( froot, pat, chan, sep='__')
    posterior_file = file.path("Output/Output_joint", paste0(outroot, "__POSTERIOR.txt") )
  
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
      mu1_prec = solve( matrix(c(1,0.7,0.7,1), ncol=2, nrow=2, byrow=TRUE) )
      mu2_prec = solve( 5*diag(2) )
      
      U_1 = matrix( c(1,0.5,0.5,1), ncol=2, nrow=2, byrow=TRUE)
      n_1 = 10
      U_2 = 3*diag(2)
      n_2 = 3
      
      alpha = 1
      beta = 1
      
      data = list(Yctrl=XY_ctrl, Nctrl=Nctrl, Ypat=XY_pat, Npat=Npat,
                  mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                  mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                  U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
      
      data$Yctrl = NULL
      data$Nctrl = 0
      
      data_priorpred = data
      data_priorpred$Yctrl = NULL
      data_priorpred$Ypat = NULL
      data_priorpred$Nctrl = 0
      data_priorpred$Npat = 0 
      
      model = jags.model(textConnection(modelstring), data=data, n.chains=n.chains) 
      
      model_priorpred = jags.model(textConnection(modelstring), data=data_priorpred) 
      
      update(model, n.iter=MCMCUpdates)
      
      converge = coda.samples(model=model, variable.names=c("mu","tau","Y_syn","z","probdiff", "pred_1", "pred_2"),
                                  n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output = coda.samples(model=model, variable.names=c("mu", "tau","Y_syn","z","probdiff", "pred_1", "pred_2"),
                                n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output_priorpred = coda.samples(model=model_priorpred,
                                          variable.names=c("mu", "tau","Y_syn","z","probdiff", "pred_1", "pred_2"),
                                          n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      plot(output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                     "probdiff", "Y_syn[1]", "Y_syn[2]")] )
      
      posterior = as.data.frame(output[[1]])
      prior = as.data.frame(output_priorpred[[1]])
      
      classifs = colMeans( posterior[, grepl('z', colnames(posterior))] )
      colnames(posterior) = colnames(output[[1]])
      colnames(prior) = colnames(output_priorpred[[1]])
      
      #predpsumm_pat=summary(output_pat_priorpred)
      pdf(file.path("PDF/PDF_joint",paste0(outroot,".pdf")),width=14,height=8.5)
      
      priorpost_pred(Yctrl=XY_ctrl, Ypat=XY_pat, posterior=posterior, prior=prior, 
                     classifs=classifs, title=paste(froot, pat)  )
      dev.off()
      
      write.table(as.numeric(classifs),file.path("Output/Output_joint",paste0(outroot,"__CLASS.txt")),
                  row.names=FALSE,quote=FALSE,col.names=FALSE)
      
      write.table(posterior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "probdiff", "Y_syn[1]", "Y_syn[2]")],
                  posterior_file, row.names=FALSE, quote=FALSE)
  }else{ # if file exists load previous data
      class_pat_file = file.path("Output/Output_joint", paste0(outroot, "__CLASS.txt"))
    }
  }
}










