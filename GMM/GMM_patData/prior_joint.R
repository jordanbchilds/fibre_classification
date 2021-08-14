library(rjags)
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

priorpost_den = function(ctrl_data, prior, title){
  # output: plots the prior and posterior regression lines and data
  x.lim = range(ctrl_data$Yctrl[,1]) + c(-1,1)
  y.lim = range(ctrl_data$Yctrl[,2]) + c(-1,1)
  
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  plot( prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], pch=20, col=myDarkGrey, cex.lab=2, cex.axis=1.5,
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Draws From Prior Predictive')
  contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )

    plot( ctrl_data$Yctrl[,1], ctrl_data$Yctrl[,2], pch=20, col=myGreen, cex.lab=2, cex.axis=1.5,
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Prior Predictive and Control Data',
        xlim=x.lim, ylim=y.lim)
  contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
  par(op)
} 

modelstring = "
model {
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  probdiff ~ dbeta(alpha, beta)
  
  # posterior distribution
  compOne ~ dmnorm( mu[,1], tau[,,1] )
  compTwo ~ dmnorm( mu[,2], tau[,,2] )
}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_patPrior"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_patPrior"), showWarnings = FALSE)

dir.create(file.path("Output/Output_jointPrior"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_jointPrior"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_jointPrior"), showWarnings = FALSE)

# 12,000 burn-in
MCMCUpdates = 2000
# 5,000 posterior draws after burn-in
MCMCUpdates_Report = 10000
# thin with lag-10- > 5,000 draws from posterior
MCMCUpdates_Thin = 1
n.chains = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping/IMC.RAW.txt"), stringsAsFactors=FALSE)


# removing unwanted info 
imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
mitochan = "VDAC1"

imcDat = imc_data[imc_data$channel %in% c(imc_chan, mitochan), ]

froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = c('P01')


# seperating the control patients 
# for( chan in imc_chan[-which(imc_chan == 'VDAC1')]){
for( chan in imc_chan){
    outroot = paste( froot, chan, sep='__')
    posterior_file = file.path("Output/Output_jointPrior", paste0(outroot, "__POSTERIOR.txt") )
    if( !file.exists(posterior_file)){
      
      control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
      Xctrl = log(control$value[control$channel==mitochan])
      Yctrl = log(control$value[control$channel==chan])
      XY_ctrl = cbind( Xctrl, Yctrl )
      
      ## PRIORS
      mu1_mean = colMeans(XY_ctrl)
      mu2_mean = mu1_mean
      mu1_prec = solve( matrix(c(0.2,0.1,0.1,0.2), ncol=2, nrow=2, byrow=TRUE) )
      mu2_prec = solve( 4*diag(2) )
      
      U_1 = matrix( c(10,7,7,10), ncol=2, nrow=2, byrow=TRUE)
      n_1 = 50
      U_2 = 3*diag(2)
      n_2 = 20
      
      alpha = 1
      beta = 1
      
      data = list(Yctrl=XY_ctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                  mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                  U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
      
      model = jags.model(textConnection(modelstring), data=data, n.chains=n.chains) 
      
      update(model, n.iter=MCMCUpdates)
      
      converge = coda.samples(model=model, variable.names=c("mu","tau","z","compOne","compTwo"),
                              n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output = coda.samples(model=model, variable.names=c("mu","tau","z","compOne","compTwo"),
                            n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
    
      plot(output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                     "compOne[1]","compOne[1]","compTwo[1]","compTwo[2]")] )
      
      prior = as.data.frame(output[[1]])

      colnames(prior) = colnames(output[[1]])
      
      #predpsumm_pat=summary(output_pat_priorpred)
      pdf(file.path("PDF/PDF_jointPrior",paste0(outroot,".pdf")),width=14,height=8.5)
      
      priorpost_den(ctrl_data=data, prior=prior, title=paste(froot, pat)  )
      dev.off()
      
      write.table(prior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "compOne[1]","compOne[2]","compTwo[1]","compTwo[2]")],
                  posterior_file, row.names=FALSE, quote=FALSE)
    }else{ # if file exists load previous data
      class_pat_file = file.path("Output/Output_jointPrior", paste0(outroot, "__CLASS.txt"))
    }
}



