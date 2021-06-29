library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

DarkGrey = rgb(169,169,159, max=255, alpha=50)

cramp = colorRamp(c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)),alpha=TRUE)
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

priorpost = function(ctrl_data=NULL, prior, posterior, data, class_posterior=NULL, classifs, title){
  # output: plots the prior and posterior regression lines and data
  
  op = par(mfrow=c(1,2)) #,mar = c(5.5,5.5,3,3))
  #for(param in c("mu","tau","probdiff", "")) {
  #  xlimp = -99
  #  if(param=="probdiff") xlimp = c(0.0,1.0)
  #  comparedensities(prior[[param]],posterior[[param]], xlab=param, xlim=xlimp)
  #}
  if(is.null(class_posterior)){
    class_posterior = colMeans( posterior[,grepl("z",colnames(posterior))])
  }
  if( is.null(ctrl_data) ){
    plot( data$Y[,1], data$Y[,2], col='darkgrey', pch=19, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"))
    contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot( data$Y[,1], data$Y[,2], col=classcols(class_posterior), pch=19, cex.lab=2,cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"),ylab=paste0("log(",chan,")"))
    contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
    title(main=title, line = -1, outer = TRUE)
    
  } else {
    plot( ctrl_data[,1], ctrl_data[,2], col='darkgrey', pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
          ylim=(range(c(ctrl_data[,2], data$Y[,2]))+c(-1,1)), xlim=(range(c(ctrl_data[,2], data$Y[,2]))+c(-1,1)) )
    points( data$Y[,1], data$Y[,2], col='green', pch=19 )
    contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
    
    plot( ctrl_data[,1], ctrl_data[,2], col=DarkGrey, pch=19, cex.lab=2,cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"),ylab=paste0("log(",chan,")"))
    points( data$Y[,1], data$Y[,2], col=classcols(class_posterior), pch=20 )
    contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
    title(main=title, line = -1, outer = TRUE)
  }
  
  ## classification and full posterior distributions
  plot( data$Y[,1], data$Y[,2], col=classcols(class_posterior), pch=20, cex.lab=2,cex.axis=1.5,
        xlab=paste0("log(",mitochan,")"),ylab=paste0("log(",chan,")"))
  contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5)
  title(main=title, line = -1, outer = TRUE)
  
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
  
  if( !is.null(ctrl_data) ){
    par(mfrow=c(1,2))
    plot( density(posterior[,'probdiff']), cex.lab=2, cex.axis=1.5,
          xlab='probdiff', ylab='density', lwd=2, col='red', main='probdiff Density')
    lines( density(rbeta(5000,data_pat$alpha_p, data_pat$beta_p)), lwd=2, col='green')
    title(main=title, line = -1, outer = TRUE)
  }
  par(op)
} 

modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dbern(probdiff)
    class[i] =  z[i] + 1
    Y[i,] ~ dmnorm(mu[,class[i]], tau[,,class[i]] )
  }
  
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha_p, beta_p)
  probdiff = ifelse( pi==1, p, 1)
  
  # posterior distribution
  z_syn ~ dbern(probdiff)
  class_syn = z_syn + 1
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
}
"

dir.create(file.path("Output_IMC"), showWarnings = FALSE)
dir.create(file.path("PDF_IMC"), showWarnings = FALSE)
dir.create(file.path("PNG_IMC"), showWarnings = FALSE)

# 12,000 burn-in
MCMCUpdates = 2000
# 5,000 posterior draws after burn-in
MCMCUpdates_Report = 5000
# thin with lag-10- > 5,000 draws from posterior
MCMCUpdates_Thin = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping/IMC.RAW.txt"), stringsAsFactors=FALSE)

colnames(imc_data)

unique(imc_data$channel)

# removing unwanted info 
imc_chan = c('SDHA','OSCP', 'VDAC1', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
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
pts = NULL # what's the point in this?
pts = grep("P", sbj, value = TRUE) 


# seperating the control patients 
for( chan in imc_chan[-which(imc_chan == 'VDAC1')]){
  
  outroot_ctrl = paste( froot, 'CTRL', chan, sep='__')
  posterior_ctrl_file = file.path("Output_IMC",paste0(outroot_ctrl,"__POSTERIOR.txt"))
  
  # data frame for mean intensity 
  control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
  
  # get data for controls 
  Xctrl = log(control$value[control$channel==mitochan])
  Yctrl = log(control$value[control$channel==chan])
  XY_ctrl = cbind( Xctrl, Yctrl )
  
  if(!file.exists(posterior_ctrl_file)){
    
    # define prior parameters
    mu1_mean = c(0,0)
    mu2_mean = c(0,0)
    mu1_prec = 0.25*diag(2)
    mu2_prec = 0.25*diag(2)
    U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_1 = 2
    U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_2 = 2
    alpha_p = 2
    beta_p = 2
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
    
    summ_ctrl = summary(output_ctrl)
    classifs_ctrl = summ_ctrl$statistics[grepl("z",rownames(summ_ctrl$statistics)),"Mean"]
    
    # print MCMC output
    #print(summ_ctrl)
    #plot(output_ctrl)
    #autocorr.plot(output_ctrl)
    #pairs(as.matrix(output_ctrl))
    #crosscorr.plot(output_ctrl)
    
    # prior and posterior prediction for control data
    predpsumm_ctrl=summary(output_ctrl_priorpred)
    ctrlroot = paste(froot,"CONTROL",chan,sep="__") 
    
    pdf(file.path("PDF_IMC",paste0(ctrlroot,".pdf")),width=14,height=7)
    
    priorpost(prior=prior_ctrl, posterior=posterior_ctrl, data=data_ctrl, 
              classifs=classifs_ctrl, title=paste(froot,"CTRL", chan, sep='__'))
    #title(paste(froot,"CTRL"), line = -1, outer = TRUE)
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
  
  n_1 = 6 # degrees of freedom
  # define the expected value of the patient prior (prec_pred) be the mean of the control
  # posterior
  prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[1,2,1]','tau[2,2,1]')]),
                      nrow=2, ncol=2, byrow=TRUE)
  
  # increase the covariance between 'x' and 'y', keep variances the same
  Sigma = solve(prec_pred)
  delta = matrix(c(1,-0.9,-0.9,1), ncol=2, nrow=2, byrow=TRUE)
  # re-define the expectation of the prior
  prec_pred = solve( Sigma + delta )
  # define prior parameter
  U_1 = prec_pred/n_1
  n_2 = 3
  U_2 = solve( matrix( c(2,0,0,2), nrow=2, ncol=2, byrow=TRUE) )/n_2
  
  mu1_mean = colMeans( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])
  mu1_prec = solve( matrix( c(1,0,0,1), ncol=2, nrow=2, byrow=TRUE) )
  
  mu2_mean = mu1_mean/2
  mu2_prec = solve( matrix( c(5,0,0,5), ncol=2, nrow=2, byrow=TRUE) )
  
  alpha_p = 2
  beta_p = 2
  pi = 1
  
  for(pat in pts){ # loop through patients
    
    outroot = paste(froot,pat,chan,sep="__")
    
    patient = imcDat[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"), ] 
    
    posterior_file = file.path("Output_IMC",paste0(outroot,"__POSTERIOR.txt"))
    
    if(!file.exists(posterior_file)){ # regression for mitochondrial disease patients
      # Block off file from analysis
      file.create(posterior_file)
      
      op = par(mfrow=c(2,3) ) #, mar = c(5.5,5.5,3,3))
      
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
      
      class_posterior_pat = posterior_pat[, grepl('z', colnames(posterior_pat))]
      colnames(posterior_pat) = colnames(output_pat[[1]])
      colnames(prior_pat) = colnames(output_pat_priorpred[[1]])
      
      summ_pat = summary(output_pat)
      #classifs_pat = summ_pat$statistics[grepl("z",rownames(summ_pat$statistics)),"Mean"]
      classifs_pat = colMeans(class_posterior_pat)
      
      #print(summ_pat)
      #plot(converge_pat)
      #autocorr.plot(converge_pat)
      #pairs(as.matrix(converge_pat))
      #crosscorr.plot(converge_pat)
      
      #predpsumm_pat=summary(output_pat_priorpred)
      pdf(file.path("PDF",paste0(outroot,".pdf")),width=14,height=8.5)
      priorpost(ctrl_data=XY_ctrl, prior=prior_pat, posterior=posterior_pat, 
                data=data_pat, classifs=classifs_pat, title=paste(froot,pat) )
      # title(paste(froot,pat), line = -1, outer = TRUE)
      dev.off()
      write.table(as.numeric(classifs_pat),file.path("Output",paste0(outroot,"__CLASS.txt")),row.names=FALSE,quote=FALSE,col.names=FALSE)
      write.table(posterior_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                   "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                   "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                   "probdiff", "Y_syn[1]", "Y_syn[2]")],posterior_file,row.names=FALSE,quote=FALSE)
    }else{ # if file exists load previous data
      
      class_pat_file = file.path("Output_IMC", paste0(outroot, "__CLASS.txt"))
    }
  }
}

