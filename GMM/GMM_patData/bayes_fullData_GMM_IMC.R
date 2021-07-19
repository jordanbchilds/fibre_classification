######################
#### A MODEL FOR ALL THE DATA, AT ONCE, IN ONE GO, BECASUE WHY NOT?
######################
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

priorpost = function(data, prior, posterior, classifs, output_mcmc=NULL, title, marginals=FALSE){
  # output: plots the prior and posterior regression lines and data
  x.lim = range(c(data$Yctrl[,1], data$Ypat[,1])) + c(-1,1)
  y.lim = range(c(data$Yctrl[,2], data$Ypat[,2])) + c(-1,1)
  
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  
  plot( data$Yctrl[,1], data$Yctrl[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Prior Predictive',
        xlim=x.lim, ylim=y.lim)
  points( data$Ypat[,1], data$Ypat[,2], col=myGreen, pch=20)
  contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
  
  plot( data$Y[,1], data$Y[,2], col=myDarkGrey, pch=20, cex.lab=2, cex.axis=1.5,
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), main='Posterior Predictive',
        xlim=x.lim, ylim=y.lim)
  points( data$Ypat[,1], data$Ypat[,2], col=classcols(classifs), pch=20)
  contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
  
  title(main=title, line = -1, outer = TRUE)
  
  if( marginals ){
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
    
    par(mfrow=c(1,2))
    plot( density(posterior[,'probdiff']), cex.lab=2, cex.axis=1.5, xlim=c(0,1),
          xlab='probdiff', ylab='density', lwd=2, col='red', main='probdiff Density')
    lines( density(rbeta(5000, data$alpha, data$beta)), lwd=2, col='green')
    title(main=title, line = -1, outer = TRUE)
  }
  if( !is.null(output_mcmc) ){
    par(mfrow=c(2,3))
    plot(output_mcmc[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                        "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                        "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                        "probdiff", "Y_syn[1]", "Y_syn[2]")])
  }
  par(op)
} 

modelstring = "
model {
  for(i in 1:length(N)){
    probdiff[i] = ifelse(i==1, 0, p)
    
    for(j in 1:N[i]){
    
      z[ pat_index[i]+j-1 ] ~ dbern( probdiff[i] )
      comp[pat_index[i]+j-1] = 2 - z[pat_index[i]+j-1] 
      Y[pat_index[i]+j-1, 1:2] ~ dmnorm(mu[,comp[pat_index[i]+j-1]], tau[,,comp[pat_index[i]+j-1]])
    }
  }

  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha, beta)

  # posterior distribution
  z_syn ~ dbern(p)
  class_syn = 2 - z_syn 
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])

}
"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_allData"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_allData"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_allData"), showWarnings = FALSE)

dir.create(file.path("PDF/PDF_allData/MCMC_output"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_allData/classification"), showWarnings = FALSE)


## tests for RJAGS
fulldat = 'IMC.RAW.txt'
imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

# removing unwanted info 
imc_chan = c('SDHA','OSCP', 'VDAC1', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imc_data_1.0 = imc_data[imc_data$channel %in% imc_chan, ]

imcDat=imc_data_1.0

mitochan = "VDAC1"

froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P", sbj, value = TRUE)

MCMCUpdates = 2000
MCMCUpdates_Report = 3000
MCMCUpdates_Thin = 1
n.chains = 1

# for( chan in imc_chan[-which(imc_chan == mitochan)]){
for( chan in c('NDUFB8')){
  outroot = paste( froot, chan, sep='__')
  posterior_file = file.path("Output/Output_allData", paste0(outroot, "__POSTERIOR.txt") )
  
  # dataset with only the current protein and VDAC1
  data_chan = imcDat[(imcDat$channel==chan)|(imcDat$channel==mitochan),]
  
  # control data for chan
  control = imcDat[(imcDat$patient_type=='control')&(imcDat$type=='mean intensity'), ]
  Xctrl = log(control$value[control$channel==mitochan])
  Yctrl = log(control$value[control$channel==chan])
  Nctrl = length(Yctrl)
  data_ctrl = cbind( Xctrl, Yctrl, rep('control', Nctrl))
  colnames(data_ctrl) = c(paste(mitochan), paste(chan), 'patient')
  
  data_chan = data_ctrl# data frame for chan (add patient data to bottom
  
  N = double(10) # store the number of observations per patient 
  N[1] = Nctrl
  
  for( j in 1:length(pts) ){
    # all the patient data for chan
    patient = imcDat[(imcDat$patient_id==pts[j])&(imcDat$type=="mean intensity"), ] 

    Xpat = log(patient$value[patient$channel==mitochan])
    Ypat = log(patient$value[patient$channel==chan]) 
    Npat = length(Xpat)
    data_pat = cbind(Xpat, Ypat, rep(paste(pts[j]), Npat))
    
    # data_pat = data.frame(Xpat, Ypat, rep(pts[j], Npat))
    colnames(data_pat) = c(paste(mitochan), paste(chan), 'patient')

    # combining with previous data
    data_chan = rbind(data_chan, data_pat)
    
    N[j+1] = Npat
  }
  # saving output for mcmc 
  Ychan = data_chan[,c(paste(mitochan), paste(chan))]
  pat_id = data_chan[,'patient']
  is_control = as.numeric( pat_id=="control" )
  
  # row index for each change in patient
  con_pts = c('control', pts)
  pat_index = double(length(con_pts))
  for(i in 1:length(con_pts)) pat_index[i] = min(which(data_chan[,'patient']==con_pts[i]))
  
  ## PRIORS
  mu1_mean = c(1,1.5)
  mu2_mean = 1.5*mu1_mean
  mu1_prec = solve( matrix(c(0.2,0.1,0.1,0.2), ncol=2, nrow=2, byrow=TRUE) )
  mu2_prec = solve( 5*diag(2) )
  
  U_1 = matrix( c(10,7,7,10), ncol=2, nrow=2, byrow=TRUE)
  n_1 = 50
  U_2 = 3*diag(2)
  n_2 = 20
  
  alpha = 1
  beta = 1
  
  
  data = list(Y=Ychan, N=N, pat_index=pat_index,
              mu1_mean=mu1_mean, mu1_prec=mu1_prec,
              mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
              U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
  
  data_priorpred = data
  data_priorpred$Y = NULL
  data_priorpred$N = 0

  model = jags.model(textConnection(modelstring), data=data, n.chains=n.chains) 
  
  model_priorpred = jags.model(textConnection(modelstring), data=data_priorpred) 
  
  update(model, n.iter=MCMCUpdates)
  
  converge = coda.samples(model=model, variable.names=c("mu","tau","Y_syn","z","probdiff"),
                          n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
  
  output = coda.samples(model=model, variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                        n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
  
  output_priorpred = coda.samples(model=model_priorpred,
                                  variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                  n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
  
  MCMCoutput = output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                         "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                         "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                         "probdiff[1]", "probdiff[2]","probdiff[3]","probdiff[4]",
                         "probdiff[5]","probdiff[6]","probdiff[7]","probdiff[8]",
                         "probdiff[9]","probdiff[10]","Y_syn[1]", "Y_syn[2]")]
  
  posterior = as.data.frame(output[[1]])
  prior = as.data.frame(output_priorpred[[1]])
  
  classifs = colMeans( posterior[, grepl('z', colnames(posterior))] )
  colnames(posterior) = colnames(output[[1]])
  colnames(prior) = colnames(output_priorpred[[1]])
  
  pat_ind = c(0,N)
  pat_ind = cumsum(pat_ind)
  
  for(i in 1:length(N)){
    
    pat = c('CTRL', pts)[i]
    class_pat = classifs[(pat_ind[i]+1):pat_ind[i+1]]
    
    write.table(as.numeric(classifs),file.path("Output/Output_allData",paste0(outroot, pts[i], "CLASS.txt", sep='__')),
                row.names=FALSE,quote=FALSE,col.names=FALSE)    
    
    data_pat = data_chan[(pat_ind[i]+1):pat_ind[i+1], ]
    
    pdf(file.path("PDF/PDF_fullData", paste0(outroot,".pdf")), width=14,height=8.5)
    priorpost( data=data_pat, prior=prior, posterior=posterior, 
               classifs=class_pat, title=paste(froot, pat[i], chan))
    dev.off()
  }
  
  write.table(posterior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                           "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                           "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                           "probdiff[1]","probdiff[2]","probdiff[3]","probdiff[4]",
                           "probdiff[5]","probdiff[6]","probdiff[7]","probdiff[8]",
                           "probdiff[9]","probdiff[10]","Y_syn[1]", "Y_syn[2]")],
              posterior_file, row.names=FALSE, quote=FALSE)

}



