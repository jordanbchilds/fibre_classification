library(rjags)
library(R2jags)
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

dens_plot = function( ctrl_data, prior, title ){
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,
       xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
       main="Component One", xlim=c(0,4), ylim=c(0,4))
  contour_one = percentiles(prior[,"compOne[1]"], prior[,"compOne[2]"])
  contour(contour_one$dens, levels=contour_one$levels, labels=contour_one$probs,
          col='blue', lwd=2, add=TRUE)
  
  title(main=title, line = -1, outer = TRUE)
  
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

plot_marg = function(comp){
  xlim11 = range(prior.0[,paste0("tau[1,1,",comp,"]")], prior.1[,paste0("tau[1,1,",comp,"]")] )
  xlim12 = range(prior.0[,paste0("tau[1,2,",comp,"]")], prior.1[,paste0("tau[1,2,",comp,"]")] )
  xlim22 = range(prior.0[,paste0("tau[2,2,",comp,"]")], prior.1[,paste0("tau[2,2,",comp,"]")] )
  
  par(mfrow=c(1,3))
  plot( density( prior.0[,paste0("tau[1,1,",comp,"]")] ), xlim=xlim11,
        xlab='', ylab='', main=expression(tau[11]) )
  lines( density( prior.1[,paste0("tau[1,1,",comp,"]")] ), col="blue")
  plot( density( prior.0[,paste0("tau[1,2,",comp,"]")] ), xlim=xlim12,
        xlab='', ylab='', main=expression(tau[12]) )
  lines( density( prior.1[,paste0("tau[1,2,",comp,"]")] ), col="blue")
  plot( density( prior.0[,paste0("tau[2,2,",comp,"]")] ), xlim=xlim22,
        xlab='', ylab='', main=expression(tau[22]) )
  lines( density( prior.1[,paste0("tau[2,2,",comp,"]")] ), col="blue")
  title(main=expression(T[1]), line=-1, outer=TRUE)
  
  par(mfrow=c(1,2))
  plot( density( prior.0[,paste0("mu[1,",comp,"]")] ), xlab='', ylab='', main=expression(mu[1]) )
  lines( density( prior.1[,paste0("mu[1,",comp,"]")] ), col="blue")
  plot( density( prior.0[,paste0("mu[2,",comp,"]")] ), xlab='', ylab='', main=expression(mu[2]) )
  lines( density( prior.1[,paste0("mu[2,",comp,"]")] ), col="blue")
  title(main=expression(mu[1]), line=-1, outer=TRUE)
  
  par(mfrow=c(1,2))
  plot(XY_ctrl[,1], XY_ctrl[,2], pch=20, col=myDarkGrey,
       xlab="log(VDAC1)", ylab="log([protein])", xlim=c(-1,4), ylim=c(-1,4))
  perc.0 = percentiles(prior.0[,paste0("mu[1,",comp,"]")], prior.0[,paste0("mu[2,",comp,"]")], probs=c(0.95))
  contour(perc.0$dens, levels=perc.0$levels, labels=perc.0$probs,
          add=TRUE)
  plot(XY_ctrl[,1], XY_ctrl[,2], pch=20, col=myDarkGrey,
       xlab="log(VDAC1)", ylab="log([protein])", xlim=c(-1,4), ylim=c(-1,4))
  perc.1 = percentiles(prior.1[,paste0("mu[1,",comp,"]")], prior.1[,paste0("mu[2,",comp,"]")], probs=c(0.95))
  contour(perc.1$dens, levels=perc.1$levels, labels=perc.1$probs, col="blue",
          add=TRUE)
  title(main=expression(mu[1]), line=-1,  outer=TRUE)
  
  comp = ifelse(comp==1, "compOne", "compTwo")
  par(mfrow=c(1,2))
  plot(XY_ctrl[,1], XY_ctrl[,2], pch=20, col=myDarkGrey,
       xlab="log(VDAC1)", ylab="log([protein])", xlim=c(-5,10), ylim=c(-5,10))
  dens.0 = percentiles( prior.0[,paste0(comp,"[1]")], prior.0[,paste0(comp,"[2]")], probs=c(0.95))
  contour( dens.0$dens, levels=dens.0$levels, labels=dens.0$probs, add=TRUE)
  
  plot(XY_ctrl[,1], XY_ctrl[,2], pch=20, col=myDarkGrey,
       xlab="log(VDAC1)", ylab="log([protein])", xlim=c(-5,10), ylim=c(-5,10))
  dens.1 = percentiles( prior.1[,paste0(comp,"[1]")], prior.1[,paste0(comp,"[2]")], probs=c(0.95))
  contour( dens.1$dens, levels=dens.1$levels, labels=dens.1$probs, 
           col="blue", add=TRUE)
  title(main="Predictive Posterior", line=-1, outer=TRUE)
  
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
  compOne ~ dmnorm(mu[,1], tau[,,1])
  compTwo ~ dmnorm(mu[,2], tau[,,2])
}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)

dir.create(file.path("Output/Output_patPrior"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_patPrior"), showWarnings = FALSE)

dir.create(file.path("Output/Output_jointPrior"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_jointPrior"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_jointPrior"), showWarnings = FALSE)

# 12,000 burn-in
MCMCBurnin = 0
# 5,000 posterior draws after burn-in
MCMCUpdates = 10000
# thin with lag-10- > 5,000 draws from posterior
MCMCThin = 1
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
imc_chan = c("MTCO1")

for( chan in imc_chan){
    outroot = paste( froot, chan, sep='__')
    posterior_file = file.path("Output/Output_jointPrior", paste0(outroot, "__POSTERIOR.txt") )
    
    control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
    Xctrl = log(control$value[control$channel==mitochan])
    Yctrl = log(control$value[control$channel==chan])
    XY_ctrl = cbind( Xctrl, Yctrl )
      
      ## PRIORS
      mu1_mean = 2*colMeans(XY_ctrl)
      mu2_mean = c(1,2)
      mu1_prec = solve( matrix( c(1,1,1,1.1), ncol=2, nrow=2, byrow=TRUE) )
      mu2_prec = 0.1*diag(2)
      
      n_1 = 10
      U_1 = matrix( c(0.5,0.5,0.5,1), ncol=2, nrow=2, byrow=TRUE)/n_1
      n_2 = 4
      U_2 = 10*diag(2)/n_2

      
      alpha = 1
      beta = 1
      
      data_prior.0 = list(mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                       mu2_mean=mu1_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                       U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
      
      data_prior.1 = list(mu1_mean=mu1_mean, mu1_prec=solve(matrix(c(0.7,0.7,0.7,0.8), ncol=2,nrow=2)), 
                          mu2_mean=mu1_mean, mu2_prec=100*diag(2), 
                          n_1=20, n_2=2, U_1=U_1*10/20, U_2=10*diag(2)/2,
                          alpha=alpha, beta=beta)
      
      model_priorpred.0 = jags(data=data_prior.0, parameters.to.save=c("mu","tau","compOne","compTwo"),
                                  model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdates, 
                                  n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
      
      model_priorpred.1 = jags(data=data_prior.1, parameters.to.save=c("mu","tau","compOne","compTwo"),
                             model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdates, 
                             n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
      
      output.0 = as.mcmc(model_priorpred.0)
      output.1 = as.mcmc(model_priorpred.1)
      
      prior.0 = as.data.frame(output.0[[1]])
      prior.1 = as.data.frame(output.1[[1]])

      colnames(prior.0) = colnames(output.0[[1]])
      colnames(prior.1) = colnames(output.1[[1]])
      
      #predpsumm_pat=summary(output_pat_priorpred)
      pdf(file.path("PDF/PDF_jointPrior",paste0(chan,"__PRIOR.pdf")), width=14,height=8.5)
      dens_plot(ctrl_data=XY_ctrl, prior=prior.1, title=paste0(chan, "__PRIOR")   )
      dev.off()
      
      write.table(prior.1[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "compOne[1]", "compOne[1]","compTwo[1]","compTwo[2]")],
                  posterior_file, row.names=FALSE, quote=FALSE)
}

plot_marg(2)

plot(XY_ctrl[,1], XY_ctrl[,2], pch=20, col=myDarkGrey,
     xlab="log(VDAC1)", ylab="log(MTCO1)")
tt = percentiles( prior.1[,"compTwo[1]"], prior.1[,"compTwo[2]"], probs=c(0.1))
contour(tt$dens, levels=tt$levels, labels=tt$probs, nlevels=5)

dens = kde2d(prior.1[,"compTwo[1]"], prior.1[,"compTwo[2]"], n=100); ## estimate the z counts
dx = diff(dens$x[1:2])
dy = diff(dens$y[1:2])
sz = sort(dens$z)
c1 = cumsum(sz) * dx * dy
levs = sapply(c(0.1), function(x) {
  approx(c1, sz, xout = 1 - x)$y
})


#### testing BUGS Wishart specification
wishart.string="
model{
  tau ~ dwish(V, df)
}"

V = matrix(c(10,0,0,10), ncol=2, nrow=2, byrow=TRUE)
df = 10
data_wish = list(V=V, df=df)

model_wish= jags(data=data_wish, parameters.to.save=c("tau"),
                       model.file=textConnection(wishart.string), n.chains=1, n.iter=10000, 
                       n.thin=1, n.burnin=0, DIC=FALSE, progress.bar="text")

output_wish = as.mcmc(model_wish)
BUGSwishart = as.data.frame(output_wish[[1]])
colnames(BUGSwishart) = colnames(output_wish[[1]])

Rwishart = rWishart(n=10000, df=df, Sigma=solve(V) )
par(mfrow=c(1,3))
plot( density(BUGSwishart[,"tau[1,1]"]), lwd=2, xlab='', ylab='', 
      main=expression(tau[11]))
lines( density(Rwishart[1,1,]), col="blue", lwd=2)
plot( density(BUGSwishart[,"tau[1,2]"]), lwd=2, xlab='', ylab='', 
      main=expression(tau[12]))
lines( density(Rwishart[1,2,]), col="blue", lwd=2)
plot( density(BUGSwishart[,"tau[2,2]"]), lwd=2, xlab='', ylab='', 
      main=expression(tau[22]))
lines( density(Rwishart[2,2,]), col="blue", lwd=2)







