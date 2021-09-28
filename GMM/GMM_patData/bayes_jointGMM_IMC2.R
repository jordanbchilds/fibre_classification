library(loo)
library(rjags)
library(R2jags)
library(MASS)
library(parallel)
source("../BootStrapping/parseData.R", local = TRUE)

cramp = colorRamp(c(rgb(1,0,0,0.2),rgb(0,0,1,0.20)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

myDarkGrey = rgb(169,169,159, max=255, alpha=50)
myGreen = rgb(25,90,0,max=255)
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
    contours = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    title(main=title, line = -1, outer = TRUE)
  } else {
    x.lim = range( data[,1], ctrl[,1] ) + c(-1,1)
    y.lim = range( data[,2], ctrl[,2] ) + c(-1,1)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=myYellow, pch=20)
    contours = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    
    points( data[,1], data[,2], col=classcols(classifs), pch=20)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    title(main=title, line = -1, outer = TRUE)
  }
  
  par(op)
} 

priorpost_marginals = function(prior, posterior, title){
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  par(mfrow=c(2,2))
  ## mu_1
  # prior
  contour( kde2d(prior[,'mu[1,1]'], prior[,'mu[2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[1,1]'], posterior[,'mu[2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Posterior Density') )
  ## mu_2
  # prior
  contour( kde2d(prior[,'mu[1,2]'], prior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[1,2]'], posterior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Posterior Density') )
  title(main=title, line = -1, outer = TRUE)
  
  
  par(mfrow=c(1,3))
  ## tau_1
  # prior
  plot( density(prior[,"tau[1,1,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[11]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,1]"]), col=myGreen, lwd=2)
  
  plot( density(prior[,"tau[1,2,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[12]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,1]"]), lwd=2, col=myGreen)
  
  plot( density(prior[,"tau[2,2,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[22]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[2,2,1]"]), lwd=2, col=myGreen)
  
  ## tau_1
  # posterior
  plot( density(prior[,"tau[1,1,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[11]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,2]"]), col=myGreen, lwd=2)
  plot( density(prior[,"tau[1,2,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[12]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[1,2,2]"]), col=myGreen, lwd=2)
  plot( density(prior[,"tau[2,2,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[22]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[2,2,2]"]), col=myGreen, lwd=2)
  title(main=title, line = -1, outer = TRUE)
  
  par(mfrow=c(1,2))
  plot( density(posterior[,"probdiff_pat"]), cex.lab=2, cex.axis=1.5, xlim=c(0,1),          
        xlab="probdiff_pat", ylab="", lwd=2, col="red", main="probdiff_pat Density")
  lines( density(prior[,"probdiff_pat"]), lwd=2, col="green")
  title(main=title, line = -1, outer = TRUE)
  
  par(op)
}

component_densities = function(ctrl_data, pat_data, posterior, prior, 
                               classifs_pat, title, chan ){
  with( inf_data, {
    par(mfrow=c(2,2))
    plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
         main="Component One", xlim=c(-1,6), ylim=c(-1,6))
    points( pat_data[,1], pat_data[,2], pch=20, col=myYellow)
    prior_one = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour(prior_one$dens, levels=prior_one$levels, labels=prior_one$probs,
            col='blue', lwd=2, add=TRUE)
    
    plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
         main="Component Two", xlim=c(-1,6), ylim=c(-1,6))
    points( pat_data[,1], pat_data[,2], pch=20, col=myYellow)
    prior_two = percentiles(prior[,"predTwo[1]"], prior[,"predTwo[2]"])
    contour(prior_two$dens, levels=prior_two$levels, labels=prior_two$probs, 
            col='red', lwd=2, add=TRUE)
    
    plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
         main="Component One", xlim=c(-1,6), ylim=c(-1,6))
    points( pat_data[,1], pat_data[,2], pch=20, col=classcols(classifs_pat))
    post_one = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour(post_one$dens, levels=post_one$levels, labels=post_one$probs,
            col="blue", lwd=2, add=TRUE)
    
    plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
         main="Component Two", xlim=c(-1,6), ylim=c(-1,6))
    points( pat_data[,1], pat_data[,2], pch=20, col=classcols(classifs_pat))
    post_two = percentiles(posterior[,"predTwo[1]"], posterior[,"predTwo[2]"])
    contour(post_two$dens, levels=post_two$levels, labels=post_two$probs, 
            col="red", lwd=2, add=TRUE)
    
    title(main=title, line = -1, outer = TRUE)
  })
}

MCMCplot = function( MCMCoutput, lag=20, title ){
  col.names = colnames(MCMCoutput[[1]])
  n.chains = length(MCMCoutput)
  
  par(mfrow=c(3,3), mar = c(5.5,5.5,4,4))
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
      plot( autocorr(MCMCoutput[[1]][,param], lags=0:lag), 
            type='h', ylim=c(-1,1), xlab='Index', ylab='' )
      for(j in 2:n.chains) lines( autocorr(MCMCoutput[[j]][,param], lags=0:20), type='h', col=j)
      
      plot(1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type='l',
           ylab='', xlab='Iteration')
      for(j in 2:n.chains) lines(MCMCoutput[[j]][,param], type='l', col=j)
      
      plot(density(MCMCoutput[[1]][,param]), main="", xlab="" )
      for(j in 2:n.chains) lines(density(MCMCoutput[[j]][,param]), col=j )
      
      title(main=title, line = -1, outer = TRUE)
    }
  }
}

inf_data = list()

inf_data$modelstring = "
model {
  # fit to ctrl data
  for( i in 1:Nctrl ){
    Yctrl[i,] ~ dmnorm( mu[,1], tau[,,1] )
    loglik[i] = logdensity.mnorm(Yctrl[i,], mu[,1], tau[,,1])
  }
  # fit to patient data
  for( j in 1:Npat ){ 
    z[j] ~ dbern(probdiff_pat)
    class[j] = 2 - z[j]
    Ypat[j,] ~ dmnorm(mu[,class[j]], tau[,,class[j]] )
    loglik[Nctrl+j] = logdensity.mnorm(Ypat[j,], mu[,class[j]], tau[,,class[j]])
  }
  # component one prior
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  # component two prior
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  probdiff_pat ~ dbeta(alpha_pat, beta_pat)
  
  # predictive distribution
  predOne ~ dmnorm(mu[,1], tau[,,1])
  predTwo ~ dmnorm(mu[,2], tau[,,2])
}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("Time"), showWarnings = FALSE)
dir.create(file.path("Information_Criteria"), showWarnings=FALSE)


dir.create(file.path("Output/IMC_joint2"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_joint2"), showWarnings = FALSE)
dir.create(file.path("Information_Criteria/IMC_joint2"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
inf_data$MCMCBurnin = 2000
inf_data$MCMCUpdate = 3000 + inf_data$MCMCBurnin
inf_data$MCMCThin = 1
inf_data$n.chains = 2

data_file = "IMC_data.txt"

if( !file.exists(data_file) ){
  url = "https://raw.githubusercontent.com/CnrLwlss/Warren_2019/master/shiny/dat.txt"
  data = read.csv(url,sep="\t", stringsAsFactors=FALSE)
  write.table(data, file=outfile, row.names=FALSE, quote=FALSE, sep="\t")
}else{
  data = read.delim(data_file, sep="\t",stringsAsFactors=FALSE)
}

inf_data$imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
inf_data$mitochan = "VDAC1"

# removing unwanted info 
inf_data$imcDat = data[data$channel %in% c(inf_data$imc_chan, inf_data$mitochan), ]

inf_data$froot = "IMC"

sbj = sort(unique(inf_data$imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
inf_data$pts = grep("P", sbj, value = TRUE)

inference = function(chan_pat){
  with(as.list(c(inf_data, chan_pat)), {
    outroot = paste( froot, chan, pat, sep='__')
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
    mutation_type = imcDat$mutation_type[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"),]
    
    mu1_mean = c(mean(Xctrl), mean(Yctrl))
    mu1_prec = solve( matrix(c(0.05, 0.068, 0.068, 0.1), ncol=2, nrow=2, byrow=TRUE) ) # correlation of ~96.167%
    
    mu2_mean = mu1_mean 
    mu2_prec = 1*diag(2) 
    
    n_1 = 2000
    U_1 = matrix(c(0.2,0.339,0.339,0.6), nrow=2,ncol=2)*n_1 # correlation of ~95% 
    n_2 = 200
    U_2 = matrix(c(5,1,1,5),nrow=2,ncol=2)*n_2

    alpha_pat = 1
    beta_pat = 1
    
    data = list(Yctrl=XY_ctrl, Nctrl=Nctrl, Ypat=XY_pat, Npat=Npat,
                mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                mu2_mean=mu2_mean, mu2_prec=mu2_prec, 
                n_1=n_1, n_2=n_2, U_1=U_1, U_2=U_2,
                alpha_pat=alpha_pat, beta_pat=beta_pat)
    
    data_priorpred = data
    data_priorpred$Yctrl = NULL
    data_priorpred$Ypat = NULL
    data_priorpred$Nctrl = 0
    data_priorpred$Npat = 0 
    
    model_jags = jags(data=data, parameters.to.save=c("mu","tau","z","probdiff_pat","predOne","predTwo","loglik"),
                      model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                      n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")
    
    model_priorpred_jags = jags(data=data_priorpred, parameters.to.save=c("mu","tau","probdiff_pat","predOne","predTwo"),
                                model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                                n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
    
    DIC = model_jags$BUGSoutput$DIC
    WAIC = waic(model_jags$BUGSoutput$sims.list$loglik)
    
    output = as.mcmc(model_jags)
    output_priorpred = as.mcmc(model_priorpred_jags)
    
    MCMCoutput = output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                           "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                           "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                           "probdiff_pat", "predOne[1]", "predOne[2]", "predTwo[1]",
                           "predTwo[2]")]
    
    posterior = as.data.frame(output[[1]])
    prior = as.data.frame(output_priorpred[[1]])
    
    zpat = colnames(posterior[,grep("z", colnames(posterior))])
    zpat.split = strsplit(zpat, split="")
    zpat.vec = double(length(zpat.split))
    for(i in seq_along(zpat.split)){
      rr = zpat.split[[i]][ !zpat.split[[i]] %in% c("z","[","]") ]
      zpat.vec[i] = as.numeric(paste(rr, collapse=""))
    }
    names(zpat.vec) = zpat
    zpat.vec = sort(zpat.vec)
    
    class_pat = posterior[, names(zpat.vec)]  
    classifs_pat = colMeans(class_pat)

    colnames(posterior) = colnames(output[[1]])
    colnames(prior) = colnames(output_priorpred[[1]])
    
    return( list(
      "plot_comp" = function() { 
        component_densities(ctrl_data=XY_ctrl, pat_data=XY_pat, posterior=posterior, 
                            prior=prior, classifs_pat=classifs_pat,
                            chan=chan, title=paste(froot, chan, pat, sep="__") ) },
      "plot_mcmc" = function() {
        MCMCplot( MCMCoutput, title=paste(froot, chan, pat, sep="__") ) },
      "plot_marg" = function() {
        priorpost_marginals(prior, posterior, title=paste(froot, chan, pat, sep="__")) },
      "output" = MCMCoutput[[1]],
      "classifs_pat" = classifs_pat,
      "DIC" = DIC,
      "WAIC" = WAIC,
      "channel" = chan,
      "patient" = pat) 
    )
  })
}

chanpat_list = list()
for(chan in inf_data$imc_chan){
  for(pat in inf_data$pts){
    chan_pat = paste(chan, pat, sep="_")
    chanpat_list[[chan_pat]] = list(chan=chan, pat=pat)
  }
}

cl  = makeCluster(21) 
clusterExport(cl, c("inference", "chanpat_list", "inf_data"))
clusterEvalQ(cl, {
  library("R2jags")
  library("loo")
})

time = system.time({
  inference_out = parLapply(cl, chanpat_list, inference)
})

stopCluster(cl)

pdf(paste0("PDF/IMC_joint2/predictive.pdf"), width=10, height=8.5)
for(chan_pat in inference_out){
  chan_pat[["plot_comp"]]()
}
dev.off()

pdf(paste0("PDF/IMC_joint2/marginal.pdf"), width=10, height=8.5)
for(chan_pat in inference_out){
  chan_pat[["plot_marg"]]()
}
dev.off()

pdf(paste0("PDF/IMC_joint2/mcmc.pdf"), width=10, height=8.5)
for(chan_pat in inference_out){
  chan_pat[["plot_mcmc"]]()
}
dev.off()

######
# mcmc posterior draws 
######
for(chan_pat in inference_out){
  write.table(chan_pat[["output"]],
              file.path("Output/IMC_joint2", paste(chan_pat$channel, chan_pat$patient, "POSTERIOR.txt", sep="__") ),
              row.names=FALSE, quote=FALSE )
  
}

######
# classifications
######
for(chan_pat in inference_out){
  write.table(chan_pat[["classifs_pat"]],
              file.path("Output/IMC_joint2", paste(chan_pat$channel, chan_pat$patient, "CLASS.txt", sep="__") ),
              row.names=FALSE, quote=FALSE, col.names=FALSE )
  
}

######
# time taken for inference
######
time_df = data.frame(time=time[3])
write.table(time_df, file=file.path("Time/IMC_joint2.txt") )

######
# DIC
######
DIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(inf_data$imc_chan), 
                dimnames=list(inf_data$pts, inf_data$imc_chan))
for(chan_pat in inference_out){
  DIC_df[chan_pat$patient, chan_pat$channel]= chan_pat[["DIC"]]
}
DICpath = "Information_Criteria/IMC_joint2/DIC.txt"
write.table(DIC_df, file=DICpath, row.names=T, quote=FALSE, col.names=T)

######
# WAIC estimate and SE
######
WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
WAIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(inf_data$imc_chan), 
                 dimnames=list(inf_data$pts, inf_data$imc_chan))
WAICse_df = WAIC_df
for(chan_pat in inference_out){
  WAIC_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic","Estimate"]
  WAICse_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic", "SE"]
}
WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
WAICse_path = "Information_Criteria/IMC_joint2/WAICse.txt"
write.table(WAIC_df, file=WAICpath, row.names=TRUE, quote=FALSE, col.names=TRUE)
write.table(WAICse_df, file=WAICse_path, row.names=TRUE, quote=FALSE, col.names=TRUE)
