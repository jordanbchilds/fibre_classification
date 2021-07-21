library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

myDarkGrey = rgb(169,169,159, max=255, alpha=50)
myGreen = rgb(25,90,0,max=255,alpha=50)
myYellow = rgb(225,200,50,max=255, alpha=50)

cramp = colorRamp(c(rgb(1,0,0,0.25),rgb(0,0,1,0.25)),alpha=TRUE)
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
  
  # par(mfrow=c(2,2))
  # ## mu_1
  # # prior
  # contour( kde2d(prior[,'mu[1,1]'], prior[,'mu[1,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
  #          main=expression(mu[1]~'Prior Density') )
  # # posterior 
  # contour( kde2d(posterior[,'mu[1,1]'], posterior[,'mu[1,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
  #          main=expression(mu[1]~'Posterior Density') )
  # ## mu_2
  # # prior
  # contour( kde2d(prior[,'mu[2,1]'], prior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
  #          main=expression(mu[2]~'Prior Density') ) 
  # # posterior
  # contour( kde2d(posterior[,'mu[2,1]'], posterior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
  #          main=expression(mu[2]~'Posterior Density') )
  # title(main=title, line = -1, outer = TRUE)
  # 
  # 
  # par(mfrow=c(2,3))
  # ## tau_1
  # # prior
  # contour( kde2d( prior[,'tau[1,1,1]'], prior[,'tau[2,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[111]), ylab=expression(tau[221]), nlevels=5,
  #          main=expression(tau[1]~'Prior Density') )
  # 
  # contour( kde2d( prior[,'tau[1,1,1]'], prior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[111]), ylab=expression(tau[121]), nlevels=5,
  #          main=expression(tau[1]~'Prior Density') )
  # 
  # contour( kde2d( prior[,'tau[2,2,1]'], prior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[221]), ylab=expression(tau[121]), nlevels=5,
  #          main=expression(tau[1]~'Prior Density') )
  # ## tau_1
  # # posterior
  # contour( kde2d( posterior[,'tau[1,1,1]'], posterior[,'tau[2,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[111]), ylab=expression(tau[221]), nlevels=5,
  #          main=expression(tau[1]~'Posterior Density') )
  # 
  # contour( kde2d( posterior[,'tau[1,1,1]'], posterior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[111]), ylab=expression(tau[121]), nlevels=5,
  #          main=expression(tau[1]~'Posterior Density') )
  # 
  # contour( kde2d( posterior[,'tau[2,2,1]'], posterior[,'tau[1,2,1]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[221]), ylab=expression(tau[121]), nlevels=5,
  #          main=expression(tau[1]~'Posterior Density') )
  # title(main=title, line = -1, outer = TRUE)
  # 
  # ## tau_2
  # # prior
  # contour( kde2d( prior[,'tau[1,1,2]'], prior[,'tau[2,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[112]), ylab=expression(tau[222]), nlevels=5,
  #          main=expression(tau[2]~'Prior Density') )
  # 
  # contour( kde2d( prior[,'tau[1,1,2]'], prior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[112]), ylab=expression(tau[122]), nlevels=5,
  #          main=expression(tau[2]~'Prior Density') )
  # 
  # contour( kde2d( prior[,'tau[2,2,2]'], prior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[222]), ylab=expression(tau[122]), nlevels=5,
  #          main=expression(tau[2]~'Prior Density') )
  # ## tau_2
  # # posterior
  # contour( kde2d( posterior[,'tau[1,1,2]'], posterior[,'tau[2,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[112]), ylab=expression(tau[222]), nlevels=5,
  #          main=expression(tau[2]~'Posterior Density') )
  # 
  # contour( kde2d( posterior[,'tau[1,1,2]'], posterior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[112]), ylab=expression(tau[122]), nlevels=5,
  #          main=expression(tau[2]~'Posterior Density') )
  # 
  # contour( kde2d( posterior[,'tau[2,2,2]'], posterior[,'tau[1,2,2]'], n=100), cex.lab=2, cex.axis=1.5,
  #          xlab=expression(tau[222]), ylab=expression(tau[122]), nlevels=5,
  #          main=expression(tau[2]~'Posterior Density') )
  # title(main=title, line = -1, outer = TRUE)
  # 
  if( !is.null(pat_data) ){
    par(mfrow=c(1,2))
    plot( density(posterior[,'probdiff']), cex.lab=2, cex.axis=1.5, xlim=c(0,1),
          xlab='probdiff', ylab='density', lwd=2, col='red', main='probdiff Density')
    lines( density(rbeta(5000,pat_data$alpha_p, pat_data$beta_p)), lwd=2, col='green')
    title(main=title, line = -1, outer = TRUE)
  }
  # if( !is.null(output_mcmc) ){
  #   par(mfrow=c(2,3))
  #   plot(output_mcmc[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
  #                       "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
  #                       "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
  #                       "probdiff", "Y_syn[1]", "Y_syn[2]")])
  #   par(mfrow=c(1,1))
  # }

  par(op)
} 

modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dbern(probdiff)
    class[i] =  2 - z[i]
    Y[i,] ~ dmnorm(mu[,class[i]], tau[,,class[i]] )
  }
  
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha_p, beta_p)
  probdiff = ifelse( pi==1, p, 0) # probability of being 'like-control'
  
  # posterior distribution
  z_syn ~ dbern(probdiff)
  class_syn = 2 - z_syn 
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
}
"
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_IMC"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_IMC"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_IMC"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
MCMCUpdates = 2000
MCMCUpdates_Report = 5000
MCMCUpdates_Thin = 3

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

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


# seperating the control patients 
for( chan in imc_chan[-which(imc_chan == 'VDAC1')]){
  outroot_ctrl = paste( froot, 'CTRL', chan, sep='__')
  posterior_ctrl_file = file.path("Output/Output_IMC",paste0(outroot_ctrl,"__POSTERIOR.txt"))
  
  # data frame for mean intensity 
  control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
  
  # get data for controls 
  Xctrl = log(control$value[control$channel==mitochan])
  Yctrl = log(control$value[control$channel==chan])
  XY_ctrl = cbind( Xctrl, Yctrl )
  
  if(!file.exists(posterior_ctrl_file)){
    
    # define prior parameters
    mu1_mean = c(4,4)
    mu2_mean = c(5,5)
    mu1_prec = 0.25*diag(2)
    mu2_prec = 0.25*diag(2)
    
    U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_1 = 2
    U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_2 = 2
    alpha_p = 1
    beta_p = 1
    pi = 0
    
    # Bayesian inference using JAGS
    # Assume a prior centred around the true values (assume good estimate from control data)
    
    Nctrl = nrow(XY_ctrl)
    # parameter list for RJAGS
    data_ctrl = list(Y=XY_ctrl, N=Nctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                     mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                     U_1=U_1, U_2=U_2, alpha_p=alpha_p, beta_p=beta_p, pi=pi)
    
    data_ctrl_priorpred = data_ctrl # same parameters used for prior prediction RJAGS code
    data_ctrl_priorpred$Y = NULL # removes for prior prediction RJAGS 
    data_ctrl_priorpred$N = 0 # N: number of observed control points, removes for prior prediction
    
    # run the JAGS model for prior prediction and control parameter inference
    model_ctrl=jags.model(textConnection(modelstring), data=data_ctrl) # no initial vals given -> draw from prior
    
    model_ctrl_priorpred=jags.model(textConnection(modelstring), data=data_ctrl_priorpred) # no initial vals given -> draw from prior
    
    update(model_ctrl, n.iter=MCMCUpdates)
    
    output_ctrl=coda.samples(model=model_ctrl,variable.names=c("mu","tau","Y_syn","z","probdiff"),
                             n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    output_ctrl_priorpred=coda.samples(model=model_ctrl_priorpred,variable.names=c("mu","tau","Y_syn", "z", "probdiff"),
                                       n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])
    
    colnames(posterior_ctrl) = colnames(output_ctrl[[1]])
    colnames(prior_ctrl) = colnames(output_ctrl_priorpred[[1]])
    
    par(mfrow=c(2,3))
    plot(output_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                        "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                        "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                        "probdiff", "Y_syn[1]", "Y_syn[2]")])
    par(mfrow=c(1,1))
    
    classifs_ctrl = colMeans(posterior_ctrl[, grepl('z', colnames(posterior_ctrl))])
    
    # prior and posterior prediction for control data
    predpsumm_ctrl=summary(output_ctrl_priorpred)
    ctrlroot = paste(froot,"CONTROL",chan,sep="__") 
    
    pdf(file.path("PDF/PDF_IMC",paste0(ctrlroot,".pdf")),width=14,height=7)
    priorpost(ctrl_prior=prior_ctrl, ctrl_posterior=posterior_ctrl, ctrl_data=data_ctrl, 
              title=paste(froot,"CTRL"))
    dev.off()
    
    
    write.table(posterior_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                  "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                  "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                  "probdiff", "Y_syn[1]", "Y_syn[2]")],
                posterior_ctrl_file,row.names=FALSE,quote=FALSE)
    
  } else {
    posterior_ctrl = read.delim(posterior_ctrl_file, sep=" ",stringsAsFactors=FALSE)
    
    colnames(posterior_ctrl) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                 "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                 "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                 "probdiff", "Y_syn[1]", "Y_syn[2]")
  }
  ###
  ### prior specification for patient data
  ###
  
 
  # define the expected value of the patient prior (prec_pred) be the mean of the control
  # posterior
  prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')]),
                      nrow=2, ncol=2, byrow=TRUE )
  
  prec_pred_inv = solve( prec_pred )
  n_1 = 2 # degrees of freedom
  
  # define prior parameter
  U_1 = prec_pred_inv*n_1
  n_2 = 2
  U_2 = solve( matrix( c(2,0,0,2), nrow=2, ncol=2, byrow=TRUE) )*n_2
  
  mu1_mean = colMeans( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])
  mu1_prec = solve( var( posterior_ctrl[,c('mu[1,1]','mu[2,1]')]) )
  
  mu2_mean = mu1_mean
  mu2_prec = mu1_prec/100
  
  alpha_p = 1
  beta_p = 1
  pi = 1
  
  for(pat in pts){ # loop through patients
    outroot = paste(froot,pat,chan,sep="__")
    
    patient = imcDat[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"), ] 
    
    posterior_file = file.path("Output/Output_IMC",paste0(outroot,"__POSTERIOR.txt"))
    
    if(!file.exists(posterior_file)){ # regression for mitochondrial disease patients
      # Block off file from analysis
      file.create(posterior_file)
      
      op = par(mfrow=c(2,3) ) 
      
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      XY_pat = cbind(Xpat, Ypat)
      
      # Bayesian inference using JAGS                                                                                                                                                  
      # Assume a prior centered around the estimate from control data
      Npat = nrow(XY_pat)
      
      data_pat = list(Y=XY_pat, N=Npat, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                      mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                      U_1=U_1, U_2=U_2, alpha_p=alpha_p, beta_p=beta_p, pi=pi)
      
      data_pat_priorpred = data_pat
      data_pat_priorpred$Y = NULL
      data_pat_priorpred$N = 0
      
      model_pat=jags.model(textConnection(modelstring), data=data_pat, n.chains=1) 
      
      model_pat_priorpred=jags.model(textConnection(modelstring), data=data_pat_priorpred) 
      update(model_pat,n.iter=MCMCUpdates)
      
      converge_pat = coda.samples(model=model_pat,variable.names=c("mu","tau","Y_syn","z","probdiff"),
                                  n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      output_pat = coda.samples(model=model_pat,variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      output_pat_priorpred = coda.samples(model=model_pat_priorpred,
                                          variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                          n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      plot(output_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                         "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                         "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                         "probdiff", "Y_syn[1]", "Y_syn[2]")] )
      
      posterior_pat = as.data.frame(output_pat[[1]])
      prior_pat = as.data.frame(output_pat_priorpred[[1]])
      
      classifs_pat = colMeans( posterior_pat[, grepl('z', colnames(posterior_pat))] )
      colnames(posterior_pat) = colnames(output_pat[[1]])
      colnames(prior_pat) = colnames(output_pat_priorpred[[1]])

      #predpsumm_pat=summary(output_pat_priorpred)
      pdf(file.path("PDF/PDF_IMC",paste0(outroot,".pdf")),width=14,height=8.5)
      priorpost(ctrl_data=data_ctrl, ctrl_prior=prior_ctrl, ctrl_posterior=posterior_ctrl, 
                pat_prior=prior_pat, pat_posterior=posterior_pat, 
                pat_data=data_pat, classifs=classifs_pat, output_mcmc=output_pat, title=paste(froot,pat)  )
      dev.off()
      
      write.table(as.numeric(classifs_pat),file.path("Output/Output_IMC",paste0(outroot,"__CLASS.txt")),row.names=FALSE,quote=FALSE,col.names=FALSE)
      write.table(posterior_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                   "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                   "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                   "probdiff", "Y_syn[1]", "Y_syn[2]")],posterior_file,row.names=FALSE,quote=FALSE)
    }else{ # if file exists load previous data
      
      class_pat_file = file.path("Output/Output_IMC", paste0(outroot, "__CLASS.txt"))
    }
  }
}


