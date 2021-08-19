library(rjags)
library(ggplot2)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

myDarkGrey = rgb(169,169,159, max=255, alpha=50)
myGreen = rgb(0,255,0,max=255,alpha=50)
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

priorpred = function(output, data_ctrl, data_pat){
  par(mfrow=c(1,2))
  contourOne = percentiles(output[,"compOne[1]"], output[,"compOne[2]"])
  contourTwo = percentiles(output[,"compTwo[1]"], output[,"compTwo[2]"])
  plot(data_ctrl[,1], data_ctrl[,2], pch=20, col=myDarkGrey, 
       xlab=paste0("log(",mitochan,")"), ylab=paste("log(",chan,")"), 
       main='Component One', xlim=c(0,10), ylim=c(0,10))
  points(data_pat[,1], data_pat[,2], pch=20, col=myYellow)
  contour( contourOne$dens, levels=contourOne$levels, labels=contourOne$probs,
           lwd=2, col="blue", add=TRUE)
  
  plot(data_ctrl[,1], data_ctrl[,2], pch=20, col=myDarkGrey, 
       xlab=paste0("log(",mitochan,")"), ylab=paste("log(",chan,")"), 
       main='Component Two', xlim=c(0,10), ylim=c(0,10))
  points(data_pat[,1], data_pat[,2], pch=20, col=myYellow)
  contour( contourTwo$dens, levels=contourTwo$levels, labels=contourTwo$probs,
           lwd=2, col='red', add=TRUE)
}

prior_marg = function(output){
  
  par(mfrow=c(2,2))
  plot( density(output[,"mu[1,1]"]), lwd=2, 
        xlab='', ylab='', main=expression(mu[11]))
  plot( density(output[,"mu[2,1]"]), lwd=2, 
        xlab='', ylab='', main=expression(mu[21]))
  plot( density(output[,'mu[1,2]']), lwd=2, 
        xlab='', ylab='', main=expression(mu[12]))
  plot(density(output[,'mu[2,2]']), lwd=2, 
       xlab='', ylab='', expression(mu[22]))
  
  par(mfrow=c(2,2))
  plot( density(output[,'tau[1,1,1]']), lwd=2, 
        xlab='', ylab='', main=expression(tau[111]))
  plot( density(output[,'tau[1,2,1]']), lwd=2,
        xlab='', ylab='', main=expression(tau[121]))
  plot( density(output[,'tau[2,2,1]']), lwd=2, 
        xlab='', ylab='', main=expression(tau[221]))
  
  par(mfrow=c(2,2))
  plot( density(output[,'tau[1,1,2]']), lwd=2, 
        xlab='', ylab='', main=expression(tau[112]))
  plot( density(output[,'tau[1,2,2]']), lwd=2,
        xlab='', ylab='', main=expression(tau[122]))
  plot( density(output[,'tau[2,2,2]']), lwd=2, 
        xlab='', ylab='', main=expression(tau[222]))
  
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

percentiles = function(xdat, ydat, probs=c(0.95, 0.5, 0.1)){
  dens = kde2d(xdat, ydat, n=200); ## estimate the z counts
  dx = diff(dens$x[1:2])
  dy = diff(dens$y[1:2])
  sz = sort(dens$z)
  c1 = cumsum(sz) * dx * dy
  levs = sapply(probs, function(x) {
    approx(c1, sz, xout = 1 - x)$y
  })
  return( list(dens=dens, levels=levs, probs=probs))
}

priorpost_pred = function(data, ctrl_posterior, pat_prior){
  # output: plots the prior and posterior regression lines and data
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  plot(data$Yctrl[,1], data$Yctrl[,2], pch=20, col=myDarkGrey, 
       main=" Control Predictive Posterior",
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"))
  contours = percentiles(ctrl_posterior[,"compOne[1]"], ctrl_posterior[,"compOne[2]"])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
  
  xRange = range(data$Yctrl[,1])
  yRange = range(data$Yctrl[,2])
  
  plot(data$Ypat[,1], data$Ypat[,2], pch=20, col=myYellow, 
       main="Patient Predictive Prior",
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xRange, ylim=yRange)
  contours = percentiles(pat_prior[,"compOne[1]"], pat_prior[,"compOne[2]"])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
  
  par(op)
} 

postprior_params = function(ctrl_posterior, pat_prior){
  
  par(mfrow=c(1,2))
  plot( density(ctrl_posterior[,"mu[1,1]"]), ylab='', xlab='', 
        main=expression(mu[11]), lwd=2)
  lines( density( pat_prior[,"mu[1,1]"]), lwd=2, col='blue')
  legend('topleft', col=c('black','blue'), lty=1, legend=c("ctrl post", "pat prior"))
  
  plot( density(ctrl_posterior[,"mu[2,1]"]), ylab='', xlab='', 
        main=expression(mu[21]), lwd=2)
  lines( density( pat_prior[,"mu[2,1]"]), lwd=2, col='blue')
  legend('topleft', col=c('black','blue'), lty=1, legend=c("ctrl post", "pat prior"))
  
  par(mfrow=c(2,2))
  plot( density(ctrl_posterior[,"tau[1,1,1]"]), ylab='', xlab='',
        main=expression(tau[11]), lwd=2)
  lines(density(pat_prior[,"tau[1,1,1]"]), lwd=2, col='blue')
  legend('topleft', legend=c("ctrl post", "pat prior"), lty=1, col=c('black', 'blue'))
  
  plot( density(ctrl_posterior[,"tau[1,2,1]"]), ylab='', xlab='',
        main=expression(tau[12]), lwd=2)
  lines(density(pat_prior[,"tau[1,2,1]"]), lwd=2, col='blue')
  legend('topleft', legend=c("ctrl post", "pat prior"), lty=1, col=c('black', 'blue'))
  
  plot( density(ctrl_posterior[,"tau[2,2,1]"]), ylab='', xlab='',
        main=expression(tau[22]), lwd=2)
  lines(density(pat_prior[,"tau[2,2,1]"]), lwd=2, col='blue')
  legend('topleft', legend=c("ctrl post", "pat prior"), lty=1, col=c('black', 'blue'))
  par(mfrow=c(1,1))
}

modelstring = "
model {
  for(i in 1:N){
    Y[i,] ~ dmnorm( mu[,1], tau[,,1])
  }
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha, beta)
  probdiff = ifelse( pi==1, p, 1)
  
  # posterior distribution
  compOne ~ dmnorm(mu[,1], tau[,,1])
  compTwo ~ dmnorm(mu[,2], tau[,,2])
}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_patPrior"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior"), showWarnings = FALSE)

dir.create(file.path("PDF/PDF_patPrior/Joint"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior/IMC"), showWarnings = FALSE)


dir.create(file.path("PDF/PDF_patPrior/IMC/Predictive"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior/IMC/Marginals"), showWarnings = FALSE)

dir.create(file.path("PDF/PDF_patPrior/Joint/Predictive"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior/Joint/Marginals"), showWarnings = FALSE)


dir.create(file.path("PNG/PNG_patPrior"), showWarnings = FALSE)


# 12,000 burn-in
MCMCUpdates = 2000
# 5,000 posterior draws after burn-in
MCMCUpdates_Report = 10000
# thin with lag-10- > 5,000 draws from posterior
MCMCUpdates_Thin = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping/IMC.RAW.txt"), stringsAsFactors=FALSE)

mitochan = "VDAC1"

# removing unwanted info 
imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imcDat = imc_data[imc_data$channel %in% c(imc_chan, mitochan), ]

froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = c('P01')

for( chan in imc_chan ){
    outroot = paste( froot, chan, sep='__')
    posterior_file = file.path("Output/Output_patPrior", paste0(outroot, "__POSTERIOR.txt") )
    
    ctrlDat = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
    
    # get data for controls 
    Yctrl = log( cbind(ctrlDat$value[ctrlDat$channel==mitochan], ctrlDat$value[ctrlDat$channel==chan])) 
    Nctrl = nrow(Yctrl)
    
    patDat = imcDat[(imcDat$patient_id==pts)&(imcDat$type=='mean intensity'), ]
    Ypat = log( cbind(patDat$value[patDat$channel==mitochan], patDat$value[patDat$channel==chan]))
    
    data = list(Yctrl=Yctrl, Ypat=Ypat)
    
    if( !file.exists(posterior_file )){
      # define prior parameters
      mu1_mean = c(4,4)
      mu2_mean = c(5,5)
      mu1_prec = 0.25*diag(2)
      mu2_prec = 0.25*diag(2)
      
      U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
      n_1 = 2
      U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
      n_2 = 2
      alpha = 1
      beta = 1
      pi = 0
      
      # Bayesian inference using JAGS
      # Assume a prior centred around the true values (assume good estimate from control data)
      
      # parameter list for RJAGS
      data_ctrl = list(Y=Yctrl, N=Nctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                       mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                       U_1=U_1, U_2=U_2, alpha=alpha, beta=beta, pi=pi)
      
      data_ctrl_priorpred = data_ctrl # same parameters used for prior prediction RJAGS code
      data_ctrl_priorpred$Y = NULL # removes for prior prediction RJAGS 
      data_ctrl_priorpred$N = 0 # N: number of observed control points, removes for prior prediction
      
      # run the JAGS model for prior prediction and control parameter inference
      model_ctrl=jags.model(textConnection(modelstring), data=data_ctrl) # no initial vals given -> draw from prior
      model_ctrl_priorpred=jags.model(textConnection(modelstring), data=data_ctrl_priorpred) # no initial vals given -> draw from prior
      update(model_ctrl, n.iter=MCMCUpdates)
      output_ctrl=coda.samples(model=model_ctrl,variable.names=c("mu","tau","z","probdiff", "compOne", "compTwo"),
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
                                  "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]",
                                  "compTwo[2]")]
      
      classifs_ctrl = colMeans(posterior_ctrl[, grepl('z', colnames(posterior_ctrl))])
      
      
      
      write.table(posterior_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                    "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                    "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                    "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]", 
                                    "compTwo[2]")],
                  posterior_file,row.names=FALSE,quote=FALSE)
      
    }else{ # if file exists load previous data
      filePath = file.path( "Output/Output_patPrior", paste("IMC", chan, 'POSTERIOR.txt', sep="__"))
      posterior_ctrl = read.delim(filePath, sep=" ", header=TRUE )
      colnames(posterior_ctrl) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                   "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                   "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                   "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]", 
                                   "compTwo[2]")
      
    }
    
    ###
    ### patient prior
    ###
    
    prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')]),
                        nrow=2, ncol=2, byrow=TRUE )
    
    n_1 = 560 # degrees of freedom
    U_1 = solve( prec_pred )/n_1
    n_2 = 10
    U_2 = solve( matrix( c(2,0,0,2), nrow=2, ncol=2, byrow=TRUE) )/n_2
    
    mu1_mean = colMeans( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])
    mu1_prec = solve( var( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])*10 )
    
    mu2_mean = mu1_mean
    mu2_prec = mu1_prec/100
    
    alpha = 1
    beta = 1
    pi = 1
    
    data_pat_priorpred = list(Y=NULL, N=0, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                    mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                    U_1=U_1, U_2=U_2, alpha=alpha, beta=beta, pi=pi)
    
    model_pat_priorpred=jags.model(textConnection(modelstring), data=data_pat_priorpred) 
    output_pat_priorpred = coda.samples(model=model_pat_priorpred,
                                        variable.names=c("mu", "tau", "z", "probdiff", "compOne", "compTwo"),
                                        n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    prior_pat = as.data.frame(output_pat_priorpred[[1]])
    colnames(prior_pat) = colnames(output_pat_priorpred[[1]])
    
    pdf(file.path("PDF/PDF_patPrior/IMC/Predictive",paste0(outroot,".pdf")),width=14,height=8.5)
    priorpost_pred(data=data, ctrl_posterior=posterior_ctrl, pat_prior=prior_pat)
    dev.off()
    
    pdf(file.path("PDF/PDF_patPrior/IMC/Marginals", paste0(outroot, ".pdf")), width=14, height=8.5)
    postprior_params(ctrl_posterior=posterior_ctrl, pat_prior=prior_pat)
    dev.off()
}

ndufb8_output = read.delim("Output/Output_patPrior/IMC__NDUFB8__POSTERIOR.txt", header=TRUE, sep=" ")
colnames(ndufb8_output) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                     "probdiff", "compOne[1]", "compOne[2]", "compTwo[1]", 
                     "compTwo[2]")

par(mfrow=c(1,1))
contours_one = percentiles(ndufb8_output[,"compOne[1]"], ndufb8_output[,"compOne[2]"])
contours_two = percentiles(ndufb8_output[,"compTwo[1]"], ndufb8_output[,"compTwo[2]"])
contour( contours_one$dens, levels=contours_one$levels, labels=contours_one$probs, col='blue', lwd=2)

#######################################
######## JOINT MODEL 
#######################################

## prior predictive 
for(chan in imc_chan){
  
  ### priors
  mu1_mean = c(1,1.5)
  mu2_mean = 2*mu1_mean
  mu1_prec = solve( matrix(c(0.2,0.1,0.1,0.2), ncol=2, nrow=2, byrow=TRUE))
  mu2_prec = solve( 5*diag(2) )
  
  U_1 = matrix( c(10,7,7,10), ncol=2, nrow=2, byrow=TRUE)
  n_1 = 10
  U_2 = U_1/5
  n_2 = 5
  
  alpha = 1
  beta = 1
  pi = 1 
  
  data_lst = list(mu1_mean=mu1_mean, mu1_prec=mu1_prec, mu2_mean=mu2_mean, 
                  mu2_prec=mu2_prec, U_1=U_1, n_1=n_1, 
                  U_2=U_2, n_2=n_2, alpha=alpha, beta=beta, pi=pi)
  data_lst$Y = NULL
  data_lst$N = 0
  
  for(pat in pts){
    ### control data
    ctrlDat = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
    
    # get data for controls 
    Yctrl = log( cbind(ctrlDat$value[ctrlDat$channel==mitochan], ctrlDat$value[ctrlDat$channel==chan])) 
    Nctrl = nrow(Yctrl)
    
    patDat = imcDat[(imcDat$patient_id==pat)&(imcDat$type=='mean intensity'),]
    Ypat = log( cbind( patDat$value[patDat$channel==mitochan], patDat$value[patDat$channel==chan]))

    model_priorpred = jags.model(textConnection(modelstring), data=data_lst) # no initial vals given -> draw from prior
    
    output_priorpred=coda.samples(model=model_priorpred,variable.names=c("mu","tau","z", "probdiff", "compOne", "compTwo"),
                                       n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    priorpred_chan = output_priorpred[[1]]
    
    pdf(file.path("PDF/PDF_patPrior/Joint/Predictive", paste0(paste(froot, chan, sep="__"),".pdf")), width=14, height=8.5 )
    priorpred(output=priorpred_chan, data_ctrl=Yctrl, data_pat=Ypat)
    dev.off()
    
    pdf(file.path("PDF/PDF_patPrior/Joint/Marginals", paste0(paste(froot, chan, sep="__"),".pdf" )), width=14, height=8.5)
    prior_marg( priorpred_chan )
    dev.off()
  }
}




























