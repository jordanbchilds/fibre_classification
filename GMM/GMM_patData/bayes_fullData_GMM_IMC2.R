library(loo)
library(rjags)
library(R2jags)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)

cramp = colorRamp(c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)),alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

myDarkGrey = rgb(169,169,159, max=255, alpha=50)
myGreen = rgb(0,255,0,max=255,alpha=50)
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
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    
    title(main=title, line = -1, outer = TRUE)
  } else {
    x.lim = range( data[,1], ctrl[,1] ) + c(-1,1)
    y.lim = range( data[,2], ctrl[,2] ) + c(-1,1)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=myYellow, pch=20)
    contours = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=classcols(classifs), pch=20)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    title(main=title, line = -1, outer = TRUE)
    
    # densities
    # plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
    #      xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
    #      xlim=x.lim, ylim=y.lim)
    # points( data[,1], data[,2], col=classcols(classifs), pch=20)
    # contours = percetniles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    # contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    # contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    # contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2)
    # 
    # title(main=title, line = -1, outer = TRUE)
  }
  
  par(op)
} 

comp_lines = function(mcmc, comp, col="black"){
  dens = percentiles(mcmc[,paste0(comp,"[1]")], mcmc[,paste0(comp,"[1]")])
  contour(dens$dens, levels=dens$levels, labels=dens$probs, 
          lwd=3, col=col)
}

priorpost_marginals = function(prior, posterior, data, title){
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
  
  par(mfrow=c(2,5))
  for(i in 1:length(pts)){
    pdiff = paste0( "pi[",i+1,"]" )
    plot( density(posterior[,pdiff]), cex.lab=2, cex.axis=1.5, xlim=c(0,1),
          xlab=pdiff, ylab='density', lwd=2, col='red', main=paste(pts[i], pdiff, 'Density'))
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

inf_data = list()
inf_data$modelstring = "
model {
  # fit model to control data
  for(i in 1:N[1]){
    z[i] = 1
    Y[i,] ~ dmnorm(mu[,z[i]],tau[,,z[i]])
    loglik[i] = logdensity.mnorm(Y[i,], mu[,z[i]], tau[,,z[i]])
  }
  # fir model to pateint data
  for(j in 2:length(N)){
    for(k in 1:N[j]){
      Y[pat_index[j]+k-1,] ~ dmnorm(mu[,z[pat_index[j]+k-1]], tau[,,z[pat_index[j]+k-1]])
      z[pat_index[j]+k-1] ~ dcat(q[j,])
      loglik[pat_index[j]+k-1] = logdensity.mnorm(Y[pat_index[j]+k-1,], mu[,z[pat_index[j]+k-1]], tau[,,z[pat_index[j]+k-1]]) 
    }
  }  
  # component one prior
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  # component two prior
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  # classifcation
  for(l in 1:length(N)){ 
    pi[l] ~ dbeta(alpha, beta) 
    q[l,1] = pi[l]; q[l,2] = 1 - pi[l]
  }
  # predictive distribution
  predOne ~ dmnorm(mu[,1], tau[,,1])
  predTwo ~ dmnorm(mu[,2], tau[,,2])
}
"
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("Time"), showWarnings = FALSE)

dir.create(file.path("Output/IMC_allData2"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_allData2"), showWarnings = FALSE)

# dir.create(file.path("PDF/IMC_allData2/MCMC"), showWarnings = FALSE)
# dir.create(file.path("PDF/IMC_allData2/classifs"), showWarnings = FALSE)
# dir.create(file.path("PDF/IMC_allData2/marginals"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_allData2/components"), showWarnings = FALSE)
# dir.create(file.path("PDF/IMC_allData2/components/pat_singular"), showWarnings = FALSE)
# dir.create(file.path("PDF/IMC_allData2/components/pat_joined"), showWarnings = FALSE)

dir.create(file.path("Information_Criteria"), showWarnings=FALSE)
dir.create(file.path("Information_Criteria/IMC_allData2"), showWarnings = FALSE)
dir.create(file.path("Information_Criteria/IMC_allData2/WAIC"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
inf_data$MCMCBurnin = 2000
inf_data$MCMCUpdate = 3000 + inf_data$MCMCBurnin
inf_data$MCMCThin = 1
inf_data$n.chains = 2

## tests for RJAGS
fulldat = 'IMC.RAW.txt'
imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imc_chan = c("SDHA", "OSCP")
inf_data$mitochan = "VDAC1"

# removing unwanted info 
inf_data$imcDat = imc_data[imc_data$channel %in% c(imc_chan, inf_data$mitochan), ]

inf_data$froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(inf_data$imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
inf_data$pts = grep("P", sbj, value = TRUE)[1:2]

# sorts the dataset into the correct for inference
chan_data = function(chan, mitochan, pts, imcDat){
  control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
  Xctrl = log(control$value[control$channel==mitochan])
  Yctrl = log(control$value[control$channel==chan])
  Nctrl = length(Yctrl)
  data_ctrl_lst = list(Y=cbind(Xctrl, Yctrl))
  
  Xchan = Xctrl
  Ychan = Yctrl
  patient_id = rep('ctrl', Nctrl)
  
  data_ctrl = data.frame(Xctrl, Yctrl,  rep('ctrl', Nctrl))
  colnames(data_ctrl) = c(mitochan, chan, 'patient')
  
  N = double(10) # store the number of observations per patient 
  N[1] = Nctrl
  
  for( j in 1:length(pts) ){
    # all the patient data for chan
    patient = imcDat[(imcDat$patient_id==pts[j])&(imcDat$type=="mean intensity"), ] 
    # patient data
    Xpat = log(patient$value[patient$channel==mitochan])
    Ypat = log(patient$value[patient$channel==chan]) 
    Npat = length(Xpat)
    # add patient data to data matrix
    Xchan = c(Xchan, Xpat)
    Ychan = c(Ychan, Ypat)
    patient_id = c(patient_id, rep(pts[j], Npat) )
    
    N[j+1] = Npat
  }
  
  data_chan = data.frame(Xchan, Ychan, patient_id)
  colnames(data_chan) = c(mitochan, chan, "patient")
  
  Ychan = data_chan[,c(mitochan, chan)]
  
  ctrl_pts = c("ctrl", pts)
  # row index for each change in patient
  pat_index = double(length(ctrl_pts))
  for(i in 1:length(ctrl_pts)) pat_index[i] = min(which(data_chan[,"patient"]==ctrl_pts[i]))
  
  return(list(
    "data" = Ychan, 
    "pat_index" = pat_index,
    "N" = N,
    "ctrlMean" = c(mean(Xctrl), mean(Yctrl))
  ))
}

inference = function(chan){
  with(as.list(c(inf_data, chan)), {
    outroot = paste( froot, chan, sep='__')
    
    chan_data = chan_data(chan, mitochan, pts, imcDat)
    
    N = chan_data$N
    data = chan_data$data
    pat_index = chan_data$pat_index
    
    ### prior specification
    mu1_mean = 1.5*chan_data$ctrlMean
    mu1_prec = solve( matrix(c(0.1,0.125,0.125,0.2), ncol=2, nrow=2, byrow=TRUE) )
    mu2_mean = mu1_mean
    mu2_prec = 0.5*diag(2) 
    
    n_1 = 100
    U_1 = matrix(c(0.3,0.48,0.48,0.9), nrow=2,ncol=2)*n_1
    n_2 = 50
    U_2 = 2.5*diag(2)*n_2
    
    alpha = 1
    beta = 1

    data = list(Y=data, N=N, pat_index=pat_index,
                mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
    
    data_priorpred = data
    data_priorpred$Y = NULL
    data_priorpred$N = 0
    
    model_jags = jags(data=data, parameters.to.save=c("mu","tau","z","pi","predOne","predTwo","loglik"),
                      model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                      n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")
    
    model_priorpred_jags = jags(data=data_priorpred, parameters.to.save=c("mu","tau", "predOne","predTwo"),
                                model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                                n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
    
    DIC = model_jags$BUGSoutput$DIC
    WAIC = waic(model_jags$BUGSoutput$sims.list$loglik)
    
    output = as.mcmc(model_jags)
    output_priorpred = as.mcmc(model_priorpred_jags)
    
    MCMCoutput = output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                           "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                           "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                           "pi[2]","pi[3]","pi[4]","pi[5]","pi[6]","pi[7]","pi[8]",
                           "pi[9]","pi[10]", "predOne[1]", "predOne[2]", "predTwo[1]",
                           "predTwo[2]")]
    
    posterior = as.data.frame(output[[1]])
    prior = as.data.frame(output_priorpred[[1]])
    
    tt = colnames(posterior[,grep("z", colnames(posterior))])
    tt.split = strsplit(tt, split="")
    tt.vec = double(length(tt.split))
    for(i in seq_along(tt.split)){
      rr = tt.split[[i]][ !tt.split[[i]] %in% c("z","[","]") ]
      tt.vec[i] = as.numeric(paste(rr, collapse=""))
    }
    names(tt.vec) = tt 
    tt.vec = sort(tt.vec)
    
    class_posterior = posterior[, names(tt.vec)]  
    classifs = colMeans(class_posterior)
    
    colnames(posterior) = colnames(output[[1]])
    colnames(prior) = colnames(output_priorpred[[1]])
    
    return( list(
      "plot_mcmc" = function(){
        MCMCplot( MCMCoutput, title=paste(froot, chan, sep="__") ) },
      "plot_marg" = function(){
        priorpost_marginals(prior, posterior, title=paste(froot, chan, sep="__")) },
      "prior_predOne" = function(){
        comp_lines(prior, comp="predOne")
      },
      "prior_predTwo" = function(){
        comp_lines(prior, comp="predTwo")
      },
      "post_predOne" = function(){
        comp_lines(posterior, comp="predOne")
      },
      "post_predTwo" = function(){
        comp_lines(posterior, comp="predTwo")
      },
      "output" = MCMCoutput[[1]],
      "classifs" = classifs,
      "DIC" = DIC,
      "WAIC" = WAIC,
      "channel" = chan
      ) 
    )
  })
}

for(i in 1:1){
chan_list = as.list(imc_chan)
names(chan_list) = imc_chan

cl  = makeCluster(4) 
clusterExport(cl, c("inference", "chan_list", "inf_data", "chan_data"))
clusterEvalQ(cl, {
  library("R2jags")
  library("loo")
})

time = system.time({
  inference_out = parLapply(cl, chan_list, inference)
})

stopCluster(cl)

tt_out = lapply(chan_list, inference)

data = lapply(chan_list, chan_data, imc_data$mitochan, inf_data$pts, inf_data$imcDat )

pdf(paste0("PDF/IMC_joint2/pred_allData.pdf"), width=10, height=8.5)
par(mfrow=c(2,2))
for(chan in chan_list){
  chan_data = data[[chan]]
  ctrl_data = chan_data[1:data[[chan]][["N"]], ]
  pat_data = chan_data[(data[[chan]][["N"]]+1):nrow(chan_data), ]
  classifs_pat = inference_out[[chan]][["classifs"]][(data[[chan]][["N"]]+1):nrow(chan_data)]
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, 
       main="", xlab=paste0("log(",mitochan,")"), ylab=psate0("log(",chan,")"),
       xlim=c(-1,6), ylim=c(-1,6))
  points(pat_data[,1], pat_data[,2], pch=20, col=myYellow)
  inference_out[[chan]][["prior_predOne"]]()
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, 
       main="", xlab=paste0("log(",mitochan,")"), ylab=psate0("log(",chan,")"),
       xlim=c(-1,6), ylim=c(-1,6))
  points(pat_data[,1], pat_data[,2], pch=20, col=myYellow)
  inference_out[[chan]][["prior_predTwo"]]()
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, 
       main="", xlab=paste0("log(",mitochan,")"), ylab=psate0("log(",chan,")"),
       xlim=c(-1,6), ylim=c(-1,6))
  points(pat_data[,1], pat_data[,2], pch=20, col=classcols(classifs_pat))
  inference_out[[chan]][["post_predOne"]]()
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, 
       main="", xlab=paste0("log(",mitochan,")"), ylab=psate0("log(",chan,")"),
       xlim=c(-1,6), ylim=c(-1,6))
  points(pat_data[,1], pat_data[,2], pch=20, col=classcols(classifs_pat))
  inference_out[[chan]][["post_predTwo"]]()
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
  write.table(chan_pat[["classifs"]],
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
DIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(imc_chan), 
                dimnames=list(inf_data$pts, imc_chan))
for(chan_pat in inference_out){
  DIC_df[chan_pat$patient, chan_pat$channel]= chan_pat[["DIC"]]
}
DICpath = "Information_Criteria/IMC_joint2/DIC.txt"
write.table(DIC_df, file=DICpath, row.names=T, quote=FALSE, col.names=T)

######
# WAIC estimate and SE
######
WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
WAIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(imc_chan), 
                 dimnames=list(inf_data$pts, imc_chan))
WAICse_df = WAIC_df
for(chan_pat in inference_out){
  WAIC_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic","Estimate"]
  WAICse_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic", "SE"]
}
WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
WAICse_path = "Information_Criteria/IMC_joint2/WAICse.txt"
write.table(WAIC_df, file=WAICpath, row.names=TRUE, quote=FALSE, col.names=TRUE)
write.table(WAICse_df, file=WAICse_path, row.names=TRUE, quote=FALSE, col.names=TRUE)

}

time = system.time({
  for( chan in imc_chan ){
    outroot = paste(froot, chan, sep='__')
    posterior_file = file.path("Output/IMC_allData2", paste0(outroot, "__POSTERIOR.txt") )
    # dataset with only the current protein and VDAC1
    data_chan = imcDat[(imcDat$channel==chan)|(imcDat$channel==mitochan),]
    # control data for chan
    control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
    Xctrl = log(control$value[control$channel==mitochan])
    Yctrl = log(control$value[control$channel==chan])
    Nctrl = length(Yctrl)
    data_ctrl_lst = list(Y=cbind(Xctrl, Yctrl))
    
    Xchan = Xctrl
    Ychan = Yctrl
    patient_id = rep('ctrl', Nctrl)
    
    data_ctrl = data.frame(Xctrl, Yctrl,  rep('ctrl', Nctrl))
    colnames(data_ctrl) = c(mitochan, chan, 'patient')
    
    N = double(10) # store the number of observations per patient 
    N[1] = Nctrl
    
    for( j in 1:length(pts) ){
      # all the patient data for chan
      patient = imcDat[(imcDat$patient_id==pts[j])&(imcDat$type=="mean intensity"), ] 
      # patient data
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      Npat = length(Xpat)
      # add patient data to data matrix
      Xchan = c(Xchan, Xpat)
      Ychan = c(Ychan, Ypat)
      patient_id = c(patient_id, rep(paste(pts[j]), Npat) )
      
      N[j+1] = Npat
    }
    
    data_chan = data.frame(Xchan, Ychan, patient_id)
    colnames(data_chan) = c(mitochan, chan, 'patient')
    
    if( !file.exists(posterior_file) ){
      # make data frame from data matrix
      Ychan = data_chan[,c(mitochan, chan)]
      
      ctrl_pts = c("ctrl", pts)
      # row index for each change in patient
      pat_index = double(length(ctrl_pts))
      for(i in 1:length(ctrl_pts)) pat_index[i] = min(which(data_chan[,'patient']==ctrl_pts[i]))
      
      ## PRIORS
      mu1_mean = c(mean(Xctrl), mean(Yctrl))
      mu2_mean = c(2,2)
      mu1_prec = solve( matrix(c(0.3,0.3,0.3,0.5), ncol=2, nrow=2, byrow=TRUE) )
      mu2_prec = solve( diag(2) )
      
      n_1 = 50
      U_1 = matrix( c(0.4,0.4,0.4,0.5), ncol=2, nrow=2, byrow=TRUE)/n_1
      n_2 = 5
      U_2 = 10*diag(2)/n_1
      
      alpha = 1
      beta = 1
      
      data = list(Y=Ychan, N=N, pat_index=pat_index,
                  mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                  mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                  U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
      
      data_priorpred = data
      data_priorpred$Y = NULL
      data_priorpred$N = 0
      
      model_jags = jags(data=data, parameters.to.save=c("mu","tau","pi","z","predOne","predTwo","loglik"),
                        model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                        n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")
      
      model_priorpred_jags = jags(data=data_priorpred, parameters.to.save=c("mu","tau","predOne","predTwo"),
                                  model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                                  n.thin=MCMCThin, n.burnin=MCMCBurnin, progress.bar="text", DIC=FALSE)
      
      DIC_df[,chan] = model_jags$BUGSoutput$DIC
      WAIC_lst[[chan]] = waic(model_jags$BUGSoutput$sims.list$loglik)$estimates
      
      output = as.mcmc(model_jags)
      output_priorpred = as.mcmc(model_priorpred_jags)
      
      MCMCoutput = output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                             "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                             "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                             "pi[2]","pi[3]","pi[4]","pi[5]","pi[6]","pi[7]","pi[8]",
                             "pi[9]","pi[10]","predOne[1]","predOne[2]", 
                             "predTwo[1]","predTwo[2]")]
      
      posterior = as.data.frame(output[[1]])
      prior = as.data.frame(output_priorpred[[1]])
      
      tt = colnames(posterior[,grep("z", colnames(posterior))])
      tt.split = strsplit(tt, split="")
      tt.vec = double(length(tt.split))
      for(i in 1:length(tt.split)){
        rr = tt.split[[i]][ !tt.split[[i]] %in% c("z","[","]") ]
        tt.vec[i] = as.numeric(paste(rr, collapse=""))
      }
      names(tt.vec) = tt 
      tt.vec = sort(tt.vec)
      
      class_posterior = posterior[, names(tt.vec)]
      
      classifs_all = colMeans( class_posterior )
      classifs_probs = classifs_all - 1
      
      colnames(posterior) = colnames(output[[1]])
      colnames(prior) = colnames(output_priorpred[[1]])
      
      write.table(posterior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "pi[2]","pi[3]","pi[4]","pi[5]","pi[6]","pi[7]","pi[8]",
                               "pi[9]","pi[10]","predOne[1]","predOne[2]", 
                               "predTwo[1]", "predTwo[2]")],
                  file=posterior_file, row.names=FALSE, quote=FALSE)
      
    } else {
      
    }
    
    pat_ind = c(0,N)
    pat_ind = cumsum(pat_ind)
    # plots for each patient
    for(i in 1:length(N)) {
      pat = ctrl_pts[i]
      outroot_pat = paste0(outroot, "__", pat)
      data_pat = data_chan[(pat_ind[i]+1):pat_ind[i+1], ]
      classifs = classifs_probs[(pat_ind[i]+1):pat_ind[i+1]] 
      
      class_filePath = file.path("PDF/IMC_allData2/classifs", paste0(outroot_pat, "__CLASSIF.pdf"))
      if( pat=='ctrl'){
        pdf(class_filePath, width=14,height=8.5)
        priorpost( data=data_pat, prior=prior, posterior=posterior,
                   classifs=classifs, title=paste(froot, pat, chan, sep='__'))
        dev.off()
      } else { 
        data_pat_lst = list(Y=cbind(data_pat[,c(mitochan, chan)]))
        pdf(class_filePath, width=14,height=8.5)
        priorpost( data=data_pat, prior=prior, posterior=posterior, ctrl=data_ctrl,
                   classifs=classifs, title=paste(froot, pat, chan, sep='__'))
        dev.off()
        
        write.table(as.numeric(classifs),file.path("Output/IMC_allData2",paste0(outroot_pat, "__CLASS.txt")),
                    row.names=FALSE,quote=FALSE,col.names=FALSE)
        
        comp_filePath = file.path("PDF/IMC_allData2/components/pat_singular", paste0(outroot_pat, "__COMPS.pdf"))
        pdf(comp_filePath, width=14, height=8.5)
        component_densities(ctrl_data=data_ctrl_lst, pat_data=data_pat_lst,
                            pat_posterior=posterior, classifs=classifs_probs,
                            title=paste(froot, chan, sep="__"))
        dev.off()
      }
    }
    
    mcmc_filePath = file.path("PDF/IMC_allData2/MCMC", paste0(outroot, "__MCMC.pdf"))
    pdf(mcmc_filePath, width=14,height=8.5)
    MCMCplot(MCMCoutput, title=paste(froot, chan, sep='__'))
    dev.off()
    
    marg_filePath = file.path("PDF/IMC_allData2/marginals", paste0(outroot, "__MARG.pdf"))
    pdf(marg_filePath, width=14, height=8.5)
    priorpost_marginals(prior=prior, posterior=posterior, data=data,
                        title=paste(froot, pat, chan, sep='__'))
    dev.off()
    
    comp_alldata_path = file.path("PDF/IMC_allData2/components/pat_joined", paste0(outroot, "__COMP.pdf"))
    pdf(comp_alldata_path, width=14, height=8.5)
    comp_dens_allData(data=data_chan, Nctrl=Nctrl, posterior=posterior,
                      classifs=classifs_probs, title=paste(froot, chan, sep='__'))
    dev.off()
  }
})

time_df = data.frame(time=time[3])
write.table(time_df, file=paste("Time/fullData2", imc_chan, sep="__"))

DICpath = file.path("Information_Criteria/IMC_allData2/DIC.txt")
write.table(DIC_df, file=DICpath)

WAICpath = "Information_Criteria/IMC_joint2/WAIC"
chan_pat = names(WAIC_lst)
for(i in 1:length(WAIC_lst)){ 
  write.table(WAIC_lst[[i]][[1]], 
              file=file.path(WAICpath, paste0(chan_pat[i], ".txt") ))
}


