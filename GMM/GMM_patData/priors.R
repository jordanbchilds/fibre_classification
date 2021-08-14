
library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

myDarkGrey = rgb(169,169,159, max=255, alpha=50)
myGreen = rgb(25,90,0,max=255,alpha=50)
myYellow = rgb(225,200,50,max=255, alpha=50)
myBlue = rgb(0,0,255, max=255 ,alpha=80)
myRed = rgb(255,0,0, max=255, alpha=80)

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
                     pat_data=NULL, pat_prior=NULL, title){
  # output: plots the prior and posterior regression lines and data
  
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  if( is.null(pat_data) ){ # plot ctrl prior and posterior
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Control Prior')
    contour( kde2d(ctrl_prior[,'compOne[1]'], ctrl_prior[,'compOne[2]'], n=100), col=myBlue, add=TRUE, nlevels=5 )
    contour( kde2d(ctrl_prior[,'compTwo[1]'], ctrl_prior[,'compTwo[2]'], n=100), col=myRed, add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
    
  } else { # plot ctrl prior-posterior and patient prior-posterior

    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Contorl Prior')
    contour( kde2d(ctrl_prior[,'compOne[1]'], ctrl_prior[,'compOne[2]'], n=100), col=myBlue, add=TRUE, nlevels=5 )
    contour( kde2d(ctrl_prior[,'compTwo[1]'], ctrl_prior[,'compTwo[2]'], n=100), col=myRed, add=TRUE, nlevels=5 )
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Control Posterior')
    contour( kde2d(ctrl_posterior[,'compOne[1]'], ctrl_posterior[,'compOne[2]'], n=100), col=myBlue, add=TRUE, nlevels=5 )
    contour( kde2d(ctrl_posterior[,'compTwo[1]'], ctrl_posterior[,'compTwo[2]'], n=100), col=myRed, add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
    
    plot( ctrl_data$Y[,1], ctrl_data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
          xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Patient Prior',
          ylim=(range(c(ctrl_data$Y[,2], pat_data$Y[,2]))+c(-1,1)), xlim=(range(c(ctrl_data$Y[,1], pat_data$Y[,1]))+c(-1,1)) )
    points(  pat_data$Y[,1], pat_data$Y[,2], col=myYellow,  pch=20 )
    contour( kde2d(pat_prior[,'compOne[1]'], pat_prior[,'compOne[2]'], n=100), col=myBlue, add=TRUE, nlevels=5)
    contour( kde2d(pat_prior[,'compTwo[1]'], pat_prior[,'compTwo[2]'], n=100), col=myRed, add=TRUE, nlevels=5)

    title(main=title, line = -1, outer = TRUE)
  }
  par(op)
} 

modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dbern(probdiff)
    class[i] =  z[i] + 1
    Y[i,1:2] ~ dmnorm(mu[,class[i]], tau[,,class[i]] )
  }
  
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha, beta)
  probdiff = ifelse( pi==1, p, 0) 
  
  # posterior distribution
  compOne ~ dmnorm(mu[,1], tau[,,1])
  compTwo ~ dmnorm(mu[,2], tau[,,2])
}
"

dir.create(file.path("Priors"), showWarnings = FALSE)
dir.create(file.path("Priors/IMC"), showWarnings = FALSE)

MCMCUpdates = 2000
MCMCUpdates_Report = 5000
MCMCUpdates_Thin = 1
n.chains = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

mitochan = "VDAC1"
imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')

# removing unwanted info 
imcDat = imc_data[imc_data$channel %in% c(imc_chan, mitochan), ]

froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P", sbj, value = TRUE)

# seperating the control patients 
# for( chan in imc_chan){
for( chan in c('NDUFB8')){
  outroot_ctrl = paste( froot, 'CTRL', chan, sep='__')
  posterior_ctrl_file = file.path("Output/IMC",paste0(outroot_ctrl,"__POSTERIOR.txt"))
  
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
    
    n_1 = 2
    n_2 = 2
    U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))*n_2
    U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))*n_2

    alpha = 1
    beta = 1
    pi = 0
    
    # Bayesian inference using JAGS
    # Assume a prior centred around the true values (assume good estimate from control data)
    
    Nctrl = nrow(XY_ctrl)
    # parameter list for RJAGS
    data_ctrl = list(Y=XY_ctrl, N=Nctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                     mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                     U_1=U_1, U_2=U_2, alpha=alpha, beta=beta, pi=pi)
    
    data_ctrl_priorpred = data_ctrl # same parameters used for prior prediction RJAGS code
    data_ctrl_priorpred$Y = NULL # removes for prior prediction RJAGS 
    data_ctrl_priorpred$N = 0 # N: number of observed control points, removes for prior prediction
    
    # run the JAGS model for prior prediction and control parameter inference
    model_ctrl=jags.model(textConnection(modelstring), data=data_ctrl) # no initial vals given -> draw from prior
    
    model_ctrl_priorpred=jags.model(textConnection(modelstring), data=data_ctrl_priorpred) # no initial vals given -> draw from prior
    
    update(model_ctrl, n.iter=MCMCUpdates)
    
    output_ctrl=coda.samples(model=model_ctrl,variable.names=c("mu","tau","z", "probdiff", "compOne", "compTwo"),
                             n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    output_ctrl_priorpred=coda.samples(model=model_ctrl_priorpred,variable.names=c("mu","tau","z", "probdiff", "compOne", "compTwo"),
                                       n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])
    
    
    
    colnames(posterior_ctrl) = colnames(output_ctrl[[1]])
    colnames(prior_ctrl) = colnames(output_ctrl_priorpred[[1]])
    
    
    MCMCoutput = output_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                "probdiff")]
    
    
    classifs_ctrl = colMeans(posterior_ctrl[, grepl('z', colnames(posterior_ctrl))])
    
    # prior and posterior prediction for control data
    ctrl_title = paste(froot, "CTRL")
    
    pdf(file.path("Priors/IMC", paste("CTRL", chan, sep="__" )), width=14,height=7)
    priorpost(ctrl_prior=prior_ctrl, ctrl_posterior=posterior_ctrl, ctrl_data=data_ctrl, 
              title=ctrl_title)
    dev.off()

    
  } 

  # define the expected value of the patient prior (prec_pred) be the mean of the control
  # posterior
  prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')]),
                      nrow=2, ncol=2, byrow=TRUE )
  
  prec_pred_inv = solve( prec_pred )
  n_1 = 2 # degrees of freedom
  
  # define prior parameter
  U_1 = prec_pred_inv*n_1
  n_2 = 2
  U_2 = solve( matrix( c(5,0,0,5), nrow=2, ncol=2, byrow=TRUE) )*n_2
  
  mu1_mean = colMeans( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])
  mu1_prec = solve( var( posterior_ctrl[,c('mu[1,1]','mu[2,1]')]) )
  
  mu2_mean = mu1_mean
  mu2_prec = solve( matrix(c(4,0,0,4), nrow=2, ncol=2, byrow=TRUE))
  
  alpha = 1
  beta = 1
  pi = 1
  
  for(pat in pts[1]){ # P01 
    outroot = paste(froot,pat,chan,sep="__")
    patient = imcDat[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"), ] 
    if(TRUE){ # regression for mitochondrial disease patients
      # Block off file from analysis
      
      op = par(mfrow=c(2,3) ) 
      
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      XY_pat = cbind(Xpat, Ypat)
      
      # Bayesian inference using JAGS                                                                                                                                                  
      # Assume a prior centered around the estimate from control data
      Npat = nrow(XY_pat)
      
      data_pat_priorpred = list(Y=NULL, N=0, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                      mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                      U_1=U_1, U_2=U_2, alpha=alpha, beta=beta, pi=pi)
      

      model_pat_priorpred=jags.model(textConnection(modelstring), data=data_pat_priorpred) 

      output_pat_priorpred = coda.samples(model=model_pat_priorpred,
                                          variable.names=c("mu","tau","z","probdiff","compOne","compTwo"),
                                          n.iter=10000,thin=MCMCUpdates_Thin)

      prior_pat = as.data.frame(output_pat_priorpred[[1]])

      colnames(prior_pat) = colnames(output_pat_priorpred[[1]])
      
      pat_title = paste(froot, pat)
      
      pdf(file.path("Priors/IMC",paste(outroot, ".pdf")), width=14,height=8.5)
      priorpost(ctrl_data=data_ctrl, ctrl_prior=prior_ctrl, ctrl_posterior=posterior_ctrl, 
                pat_prior=prior_pat, pat_data=data_pat_priorpred, title=pat_title  )
      dev.off()
  
      
      }
  }
}

## mu_1 ctrl posterior vs pat prior
par(mfrow=c(1,1))
contour( kde2d(posterior_ctrl[,'mu[1,1]'], posterior_ctrl[,'mu[2,1]'], n=20), col='green', nlevels=5, main=expression(mu[1]~'Ctrl Post Vs Pat Prior'))
contour( kde2d(prior_pat[,'mu[1,1]'], prior_pat[,'mu[2,1]'], n=20), add=TRUE, col='blue', nlevels=5 )
legend('topleft', legend=c('ctrl posterior', 'pat prior'), lty=1, col=c('green', 'blue'))

contour( kde2d(posterior_ctrl[,'mu[1,2]'], posterior_ctrl[,'mu[2,2]'], n=20), col='green', nlevels=5, main=expression(mu[1]~'Ctrl Post Vs Pat Prior'))
contour( kde2d(prior_pat[,'mu[1,2]'], prior_pat[,'mu[2,2]'], n=20), add=TRUE, col='blue', nlevels=5 )
legend('topleft', legend=c('ctrl posterior', 'pat prior'), lty=1, col=c('green', 'blue'))


wishes = rWishart(n=5000, Sigma=prec_pred/n_1, df=n_1)
par(mfrow=c(2,3))
plot( density( prior_pat[,'tau[1,1,1]']), col='green', main=expression(T[1]), xlab=expression(tau[11]))
lines( density( posterior_ctrl[,'tau[1,1,1]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)
plot( density( prior_pat[,'tau[1,2,1]']), col='green', main=expression(T[1]), xlab=expression(tau[12]))
lines( density( posterior_ctrl[,'tau[1,2,1]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)
plot( density( prior_pat[,'tau[2,2,1]']), col='green', main=expression(T[1]), xlab=expression(tau[22]))
lines( density( posterior_ctrl[,'tau[2,2,1]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)

plot( density( prior_pat[,'tau[1,1,2]']), col='green', main=expression(T[2]), xlab=expression(tau[11]))
lines( density( posterior_ctrl[,'tau[1,1,2]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)
plot( density( prior_pat[,'tau[1,2,2]']), col='green', main=expression(T[2]), xlab=expression(tau[12]))
lines( density( posterior_ctrl[,'tau[1,2,2]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)
plot( density( prior_pat[,'tau[2,2,2]']), col='green', main=expression(T[2]), xlab=expression(tau[22]))
lines( density( posterior_ctrl[,'tau[2,2,2]']), col='blue')
legend('topright', legend=c('ctrl post', 'pat prior'), col=c('blue', 'green'), lty=1)
par(mfrow=c(1,1))
