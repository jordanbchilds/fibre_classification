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

priorpost_marginals = function(prior, posterior, title){
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  par(mfrow=c(2,2))
  ## mu_1
  # prior
  contour( kde2d(prior[,"mu[1,1]"], prior[,"mu[2,1]"], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~"Prior Density") )
  # posterior
  contour( kde2d(posterior[,"mu[1,1]"], posterior[,"mu[2,1]"], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~"Posterior Density") )
  ## mu_2
  # prior
  contour( kde2d(prior[,"mu[1,2]"], prior[,"mu[2,2]"], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~"Prior Density") )
  # posterior
  contour( kde2d(posterior[,"mu[1,2]"], posterior[,"mu[2,2]"], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~"Posterior Density") )
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
  plot( density(posterior[,"probdiff"]), cex.lab=2, cex.axis=1.5, xlim=c(0,1),          
        xlab="probdiff", ylab="", lwd=2, col="red", main="probdiff Density")
  lines( density(prior[,"probdiff"]), lwd=2, col="green")
  title(main=title, line = -1, outer = TRUE)
  
  par(op)
}

MCMCplot = function( MCMCoutput, lag=20, title ){
  col.names = colnames(MCMCoutput[[1]])
  n.chains = length(MCMCoutput)
  
  par(mfrow=c(3,3), mar = c(5.5,5.5,4,4))
  if( n.chains==1 ){
    for(param in col.names){
      plot( 0:lag, autocorr(MCMCoutput[[1]][,param], lags=0:lag), 
            type="h", ylim=c(-1,1), xlab="Index", ylab="")
      plot( 1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type="l",
            ylab="", xlab="Iteration")
      plot(density(MCMCoutput[[1]][,param] ), main="", xlab="")
    }
  } else {
    for( param in col.names){
      plot( autocorr(MCMCoutput[[1]][,param], lags=0:lag), 
            type="h", ylim=c(-1,1), xlab="Index", ylab="" )
      for(j in 2:n.chains) lines( autocorr(MCMCoutput[[j]][,param], lags=0:20), type="h", col=j)
      
      plot(1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type="l",
           ylab="", xlab="Iteration")
      for(j in 2:n.chains) lines(MCMCoutput[[j]][,param], type="l", col=j)
      
      plot(density(MCMCoutput[[1]][,param]), main="", xlab="" )
      for(j in 2:n.chains) lines(density(MCMCoutput[[j]][,param]), col=j )
      
      title(main=title, line = -1, outer = TRUE)
    }
  }
}

ctrl_plot = function(ctrl_data, prior, posterior, chan, title){
  par(mfrow=c(1,2))
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, 
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Prior Predictive", cex.lab=1.2, cex.main=1.4)
  prior_dens = percentiles(prior[,"exp_predOne[1]"], prior[,"exp_predOne[2]"])
  contour(prior_dens$dens, levels=prior_dens$levels, labels=prior_dens$probs, 
          lwd=2, col="blue", add=TRUE)
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey,        
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Posterior Predictive", cex.lab=1.2, cex.main=1.4)
  post_dens = percentiles(posterior[,"exp_predOne[1]"], posterior[,"exp_predOne[2]"])
  contour(post_dens$dens, levels=post_dens$levels, labels=post_dens$probs, 
          lwd=2, col="blue", add=TRUE)
  
  title(main=title, outer=TRUE, line=-1)
}

pat_plot = function(ctrl_data, pat_data, prior, posterior, classifs, chan, pat, title){
  x.lim = range(ctrl_data[,1], pat_data[,1])+c(0,1)
  y.lim = range(ctrl_data[,2], pat_data[,2])+c(0,1)
  
  par(mfrow=c(2,2))
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, xlim=x.lim, ylim=y.lim,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"),
       main="Prior Component One", cex.lab=1.2, cex.main=1.4)
  points(pat_data[,1], pat_data[,2], col=myYellow, pch=20)
  priorOne_dens = percentiles(prior[,"exp_predOne[1]"], prior[,"exp_predOne[2]"])
  contour(priorOne_dens$dens, levels=priorOne_dens$levels, labels=priorOne_dens$probs,
          lwd=2, col="blue", add=TRUE)
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, xlim=x.lim, ylim=y.lim,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"),
       main="Prior Component Two", cex.lab=1.2, cex.main=1.4)
  points(pat_data[,1], pat_data[,2], col=myYellow, pch=20)
  priorTwo_dens = percentiles(prior[,"exp_predTwo[1]"], prior[,"exp_predTwo[2]"])
  contour(priorTwo_dens$dens, levels=priorTwo_dens$levels, labels=priorTwo_dens$probs,
          lwd=2, col="red", add=TRUE)
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, xlim=x.lim, ylim=y.lim,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"),
       main="Posterior Component One", cex.lab=1.2, cex.main=1.4)
  points(pat_data[,1], pat_data[,2], col=classcols(classifs), pch=20)
  postOne_dens = percentiles(posterior[,"exp_predOne[1]"], posterior[,"exp_predOne[2]"])
  contour(postOne_dens$dens, levels=postOne_dens$levels, labels=postOne_dens$probs,
          lwd=2, col="blue", add=TRUE)
  
  plot(ctrl_data[,1], ctrl_data[,2], pch=20, col=myDarkGrey, xlim=x.lim, ylim=y.lim,
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"),
       main="Posterior Component Two", cex.lab=1.2, cex.main=1.4)
  points(pat_data[,1], pat_data[,2], col=classcols(classifs), pch=20)
  postTwo_dens = percentiles(posterior[,"exp_predTwo[1]"], posterior[,"exp_predTwo[2]"])
  contour(postTwo_dens$dens, levels=postTwo_dens$levels, labels=postTwo_dens$probs,
          lwd=2, col="red", add=TRUE)
  
  title(main=title, outer=TRUE, line=-1)
}

inf_data = list()

inf_data$modelstring = "
model {
  for(i in 1:N){
    class[i] ~ dbern(probdiff)
    z[i] = 2 - class[i]
    Y[i,] ~ dmnorm( mu[,z[i]], tau[,,z[i]])
    loglik[i] = logdensity.mnorm(Y[i,], mu[,z[i]], tau[,,z[i]])
  }
  # priors
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  pi ~ dbeta(alpha, beta)
  probdiff = ifelse( p==1, pi, 1) 

  # predictive distribution
  predOne ~ dmnorm(mu[,1], tau[,,1])
  predTwo ~ dmnorm(mu[,2], tau[,,2])
  
  exp_predOne[1] = exp(predOne[1])
  exp_predOne[2] = exp(predOne[2])
  exp_predTwo[1] = exp(predTwo[1])
  exp_predTwo[2] = exp(predTwo[2])

}"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("Time"), showWarnings = FALSE)
dir.create(file.path("Information_Criteria"), showWarnings=FALSE)

dir.create(file.path("Output/IMC_logNorm"), showWarnings = FALSE)
dir.create(file.path("PDF/IMC_logNorm"), showWarnings = FALSE)
dir.create(file.path("Information_Criteria/IMC_logNorm"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
inf_data$MCMCBurnin = 2000
inf_data$MCMCUpdate = 3000 + inf_data$MCMCBurnin
inf_data$MCMCThin = 1
inf_data$n.chains = 2

fulldat = "IMC.RAW.txt"

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

inf_data$imc_chan = c("SDHA","OSCP", "GRIM19", "MTCO1", "NDUFB8", "COX4+4L2", "UqCRC2")
inf_data$mitochan = "VDAC1"

# removing unwanted info 
inf_data$imcDat = imc_data[imc_data$channel %in% c(inf_data$imc_chan, inf_data$mitochan), ]
inf_data$froot = gsub(".RAW.txt", "", fulldat)

# getting the ranges of the axis

sbj = sort(unique(inf_data$imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
inf_data$pts = grep("P", sbj, value = TRUE)

inference = function( chan ){
  with(as.list(c(inf_data, chan)), {
    output_list = list()
    outroot = paste( froot, chan, sep="__")
    ## CONTROL DATA
    control = imcDat[(imcDat$patient_type=="control")&(imcDat$type=="mean intensity"), ]
    Xctrl = log(control$value[control$channel==mitochan])
    Yctrl = log(control$value[control$channel==chan])
    XY_ctrl = cbind( Xctrl, Yctrl )
    logXY_ctrl = log(XY_ctrl)
    Nctrl = nrow(XY_ctrl)
    
    mu1_mean = colMeans(logXY_ctrl)
    mu1_prec = solve( matrix(c(0.5,0,0,0.5), ncol=2, nrow=2, byrow=TRUE) ) 
    
    mu2_mean = mu1_mean 
    mu2_prec = 0.5*diag(2) 
    
    n_1 = 5
    U_1 = matrix(c(1,0,0,1), nrow=2,ncol=2)*n_1 
    n_2 = 50
    U_2 = matrix(c(5,1,1,5),nrow=2,ncol=2)*n_2
    
    alpha = 1
    beta = 1
    p = 1
    
    data_ctrl = list(Y=logXY_ctrl, N=Nctrl,
                     mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                     mu2_mean=mu2_mean, mu2_prec=mu2_prec, 
                     n_1=n_1, U_1=U_1,
                     n_2=n_2, U_2=U_2, 
                     alpha=alpha, beta=beta, p=p)
    
    data_ctrl_priorpred = data_ctrl
    data_ctrl_priorpred$Y = NULL
    data_ctrl_priorpred$N = 0
    
    model_ctrl = jags(data=data_ctrl, parameters.to.save=c("mu","tau","exp_predOne","exp_predTwo","loglik"),
                      model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                      n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")
    
    model_ctrl_priorpred = jags(data=data_ctrl_priorpred, parameters.to.save=c("mu","tau","probdiff","exp_predOne","exp_predTwo"),
                                model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                                n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
    
    output_ctrl = as.mcmc(model_ctrl)
    output_ctrl_priorpred = as.mcmc(model_ctrl_priorpred)
    
    MCMC_ctrl = output_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                     "exp_predOne[1]", "exp_predOne[2]", 
                                     "exp_predTwo[1]", "exp_predTwo[2]")]
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])
    
    colnames(posterior_ctrl) = colnames(output_ctrl[[1]])
    colnames(prior_ctrl) = colnames(output_ctrl_priorpred[[1]])
    
    posterior_ctrl_filepath = file.path("Output/IMC_logNorm", paste0("IMC__", chan,".txt"))
    
    write.table(posterior_ctrl, file=posterior_ctrl_filepath, sep=" ")
    
    output_list[["ctrl_plotter"]] = function(){
      ctrl_plot(ctrl_data=XY_ctrl, prior=prior_ctrl, posterior=posterior_ctrl,
                chan=chan, title=paste("IMC", chan, sep=" "))
    }
    output_list[["ctrl_mcmc"]] = function(){
      MCMCplot(MCMC_ctrl, title=paste("IMC", chan, sep=" "))
    }
    output_list[["ctrl_marg"]] = function(){
      priorpost_marginals(prior=prior_ctrl, posterior=posterior_ctrl, 
                          title=paste("IMC", chan, sep=" "))
    }
    
    ####
    ## inference for patient data
    ####
    
    prec_pred = matrix( colMeans(posterior_ctrl[,c("tau[1,1,1]", "tau[1,2,1]", "tau[2,1,1]","tau[2,2,1]")]),
                        nrow=2, ncol=2, byrow=TRUE )
    
    var_pred = solve(prec_pred)
    # n1 = var(posterior_ctrl[,"tau[1,1,1]"]) / (2*(var_pred[1,1]^2))
    # n2 = var(posterior_ctrl[,"tau[1,2,1]"]) / (var_pred[1,2]^2 + var_pred[1,1]*var_pred[2,2])
    # n3 = var(posterior_ctrl[,"tau[2,2,1]"]) / (2*(var_pred[2,2]^2))
    # 
    # output_list[["df_est"]] = c(n1,n2,n3)
    
    mu1_pred = colMeans(posterior_ctrl[,c("mu[1,1]", "mu[2,1]")])
    mu1_var = 10*var(posterior_ctrl[,c("mu[1,1]","mu[2,1]")])
    mu1_prec = solve(mu1_var)
    
    mu2_prec = 2*diag(2) 
    mu2_mean = mu1_mean 
    
    
    n_1 = 100
    U_1 = var_pred*n_1    
    n_2 = 50
    U_2 = matrix(c(1,0,0,1),nrow=2,ncol=2)*n_2
    
    alpha = 1
    beta = 1
    p=0
    
    for( pat in pts ){
      ## PATIENT DATA
      patient = imcDat[(imcDat$patient_id==pat)&(imcDat$type=="mean intensity"), ] 
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      XY_pat = cbind(Xpat, Ypat)
      logXY_pat = log(XY_pat)
      Npat = nrow(XY_pat)
      
      data_pat = list(Y=logXY_pat, N=Npat,
                      mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                      mu2_mean=mu2_mean, mu2_prec=mu2_prec, 
                      n_1=n_1, U_1=U_1,
                      n_2=n_2, U_2=U_2, 
                      alpha=alpha, beta=beta, p=p)
      
      data_pat_priorpred = data_pat
      data_pat_priorpred$Y = NULL
      data_pat_priorpred$N = 0
      
      model_pat = jags(data=data_pat, parameters.to.save=c("mu","tau", "z","probdiff","exp_predOne","exp_predTwo","loglik"),
                       model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                       n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")
      
      model_pat_priorpred = jags(data=data_pat_priorpred, parameters.to.save=c("mu","tau","probdiff","exp_predOne","exp_predTwo"),
                                 model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
                                 n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=FALSE, progress.bar="text")
      
      output_pat = as.mcmc(model_pat)
      output_pat_priorpred = as.mcmc(model_pat_priorpred)
      
      posterior_pat = as.data.frame(output_pat[[1]])
      prior_pat = as.data.frame(output_pat_priorpred[[1]])
      
      MCMC_pat = output_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                               "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                               "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                               "exp_predOne[1]", "exp_predOne[2]", 
                               "exp_predTwo[1]", "exp_predTwo[2]")]
      
      tt = colnames(posterior_pat[,grep("z", colnames(posterior_pat))])
      tt.split = strsplit(tt, split="")
      tt.vec = double(length(tt.split))
      for(i in seq_along(tt.split)){
        rr = tt.split[[i]][ !tt.split[[i]] %in% c("z","[","]") ]
        tt.vec[i] = as.numeric(paste(rr, collapse=""))
      }
      names(tt.vec) = tt 
      tt.vec = sort(tt.vec)
      
      class_posterior = posterior_pat[, names(tt.vec)]  
      classifs = colMeans(class_posterior)
      
      posterior_pat_filepath = file.path("Output/IMC_logNorm", paste0("IMC__", chan, "__", pat, ".txt"))
      write.table(posterior_pat, file=posterior_pat_filepath, sep=" ")
      
      classifs_filepath = file.path("Output/IMC_logNorm", paste0("IMC__", chan, "__", pat, ".txt"))
      write.table(posterior_pat, file=classifs_filepath, sep=" ")
      
      output_list[[paste0(pat, "_plot")]] = function(){
        pat_plot(pat_data=XY_pat, ctrl_data=XY_ctrl, 
                 prior=prior_pat, posterior=posterior_pat, classifs=classifs, 
                 chan=chan, pat=pat, title=paste("IMC", chan, pat, sep="  "))
      }
      output_list[[paste0(pat,"_logPlot")]] = function(){
        pat_plot(pat_data=logXY_pat, ctrl_data=logXY_ctrl, prior=log(prior_pat),
                 posterior=log(posterior_pat), classifs=classifs, 
                 chan=chan, pat=pat, title=paste("IMC", chan, pat, sep="  "))
      }
      output_list[[paste0(pat, "_mcmc")]] = function(){
        MCMCplot(MCMC_pat, title=paste("IMC", chan, pat, sep="  "))
      }
      output_list[[paste0(pat, "_marg")]] = function(){
        priorpost_marginals(prior=prior_pat, posterior=posterior_pat, 
                            title=paste("IMC", chan, pat, sep="  "))
      }
    }
    return(output_list)
  })
}

chan_list = as.list(inf_data$imc_chan)
names(chan_list) = inf_data$imc_chan

cl  = makeCluster(7)
clusterExport(cl, c("inference", "chan_list", "inf_data"))
clusterEvalQ(cl, {
  library("R2jags")
  library("loo")
})

time = system.time({
  inference_out = parLapply(cl, chan_list, inference)
})

stopCluster(cl)

# CTRL prior and posterior
pdf("PDF/IMC_logNorm/ctrl_PRED.pdf", width=10, height=8.5)
for(chan in names(inference_out)){
  inference_out[[chan]][["ctrl_plotter"]]()
}
dev.off()

# PAT log Normal model
pdf("PDF/IMC_logNorm/pat_PRED.pdf", width=14, height=8.5)
for(chan in names(inference_out)){
  for(pat in inf_data$pts){
    inference_out[[chan]][[paste0(pat, "_plot")]]()
  }
}
dev.off()

# PAT Normal on logged (twice) data
pdf("PDF/IMC_logNorm/logpat_PRED.pdf", width=14, height=8.5)
for(chan in names(inference_out)){
  for(pat in inf_data$pts){
    inference_out[[chan]][[paste0(pat, "_logPlot")]]()
  }
}
dev.off()

# MCMC output
pdf("PDF/IMC_logNorm/MCMC.pdf", width=14, height=8.5)
for(chan in names(inference_out)){
  inference_out[[chan]][["ctrl_mcmc"]]()
  for(pat in inf_data$pts){
    inference_out[[chan]][[paste0(pat, "_mcmc")]]()
  }
}
dev.off()

# marginal distributions
pdf("PDF/IMC_logNorm/marginals.pdf", width=14, height=8.5)
for(chan in names(inference_out)){
  inference_out[[chan]][["ctrl_marg"]]()
  for(pat in inf_data$pts){
    inference_out[[chan]][[paste0(pat, "_marg")]]()
  }
}
dev.off()





























# ######
# # time taken for inference
# ######
# time_df = data.frame(time=time[3])
# write.table(time_df, file=file.path("Time/IMC_joint2.txt") )
# 
# ######
# # DIC
# ######
# DIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(imc_chan), 
#                 dimnames=list(inf_data$pts, imc_chan))
# for(chan_pat in inference_out){
#   DIC_df[chan_pat$patient, chan_pat$channel]= chan_pat[["DIC"]]
# }
# DICpath = "Information_Criteria/IMC_joint2/DIC.txt"
# write.table(DIC_df, file=DICpath, row.names=T, quote=FALSE, col.names=T)
# 
# ######
# # WAIC estimate and SE
# ######
# WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
# WAIC_df = matrix(NA, nrow=length(inf_data$pts), ncol=length(imc_chan), 
#                  dimnames=list(inf_data$pts, imc_chan))
# WAICse_df = WAIC_df
# for(chan_pat in inference_out){
#   WAIC_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic","Estimate"]
#   WAICse_df[chan_pat$patient, chan_pat$channel] = chan_pat[["WAIC"]][[1]]["waic", "SE"]
# }
# WAICpath = "Information_Criteria/IMC_joint2/WAIC.txt"
# WAICse_path = "Information_Criteria/IMC_joint2/WAICse.txt"
# write.table(WAIC_df, file=WAICpath, row.names=TRUE, quote=FALSE, col.names=TRUE)
# write.table(WAICse_df, file=WAICse_path, row.names=TRUE, quote=FALSE, col.names=TRUE)
# 
# 
# 
