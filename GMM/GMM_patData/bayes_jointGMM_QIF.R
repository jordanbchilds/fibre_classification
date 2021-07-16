#################################################################
#########              THE CONOR PROCESS 
#################################################################
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
  for(i in 1:Nctrl){
    Yctrl[i,] ~ dmnorm(mu[,1], tau[,,1] )
  }
  for(i in 1:Npat){
    z[i] ~ dbern(probdiff)
    class[i] = 2 - z[i]
    Ypat[i,] ~ dmnorm(mu[,class[i]], tau[,,class[i]])
  }

  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  probdiff ~ dbeta(alpha, beta)

  # posterior distribution
  z_syn ~ dbern(probdiff)
  class_syn = z_syn + 1
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
  
  pred_1 ~ dmnorm(mu[,1], tau[,,1])
  pred_2 ~ dmnorm(mu[,2], tau[,,2])
}
"
# create files to store raw output and figures
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output/Output_jointQIF"),  showWarnings = FALSE)

dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_jointQIF"), showWarnings = FALSE)

dir.create(file.path("PNG"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_jointQIF"), showWarnings = FALSE)

# all datasets
fulldats = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

# 10,000 burn-in
MCMCUpdates = 2000
# 10,000 posterior draws after burn-in
MCMCUpdates_Report = 10000
# thin with lag-10 - > 1,000 draws from posterior
MCMCUpdates_Thin = 1
# number of inferences
n.chains = 1

# choose which datasets to estimate parameters for
# for(fulldat in fulldats){
# for(fulldat in c("IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt")){
for(fulldat in c("IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt")){
  # removes '.RAW.txt' from fulldat, stores as froot
  froot = gsub(".RAW.txt","",fulldat)
  
  mchans = c("Ch1","Ch2","Ch3","Ch4","Area")
  if(grepl(".TCF.", fulldat)){# attention confocal - swap channels
    chans=c("LAMA1","VDAC1","MTCO1","NDUFB8","Area") # swapped MTCO1 and VDAC1 for confocal
  }else{
    chans=c("LAMA1","MTCO1","VDAC1","NDUFB8","Area") # standard for CD7
  }
  
  # names each element of chans with elements of mchans
  names(chans) = mchans
  
  dat = read.delim(file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)
  # maps {Ch1,Ch2,Ch3,Ch4,Area} to {LAMA1,MTCO1,VDAC1,NDUFB8,Area}
  dat$Channel = chans[dat$Channel] 
  
  write.table(dat,gsub(".RAW", ".RAW_ren", fulldat),row.names=FALSE,quote=FALSE,sep="\t")
  
  cord = c("NDUFB8","MTCO1","VDAC1") # unsure 
  chlabs = c("CI","CIV","OMM") # unsure
  names(chlabs) = cord # unsure
  mitochan = "VDAC1"
  correctnpc = TRUE
  updatechans = FALSE
  
  # dat = getData(gsub(".RAW", ".RAW_ren", fulldat), cord,
  #               mitochan = mitochan,updatechans = updatechans, correctnpc = correctnpc)
  # dat$fn = gsub("_.0", "", dat$filename)
  # dat$pch = paste(dat$fn,dat$ch,sep="_")
  
  correctnpc = FALSE
  
  dat = getData(gsub(".RAW", ".RAW_ren", fulldat), c(cord,"Area"),
                mitochan = mitochan,updatechans = updatechans, correctnpc = correctnpc)
  dat$fn = gsub("_.0", "", dat$filename)
  dat$pch = paste(dat$fn,dat$ch,sep="_")
  
  # if not correcting, NPC need to be removed
  if(!correctnpc) {
    dat = dat[grepl("OXPHOS",dat$filename),]
    dat$filename = gsub(" OXPHOS","",dat$filename)
    dat$fn = gsub(" OXPHOS","",dat$fn)
    dat$pch = gsub("OXPHOS_","",dat$pch)
  }
  
  
  # Get plot axis ranges
  lims = list()
  for(ch in cord){
    lims[[ch]] = quantile(log(dat$value[dat$Channel==ch]),c(0.001,0.999),na.rm=TRUE)
  }
  # Merge different regions of the same section
  dat$fn = gsub("_R1","",dat$fn)
  dat$fn = gsub("_R2","",dat$fn)
  
  # grabbing and seperating ctrls and patients
  sbj = sort(unique(dat$fn))
  crl = grep("C._H", sbj, value = TRUE)
  pts = grep("P.", sbj, value = TRUE) 
  
  for(chan in c("NDUFB8","MTCO1")){
    for(pat in pts){
      # froot: data name 
      outroot = paste(froot,pat,chan,sep="__")
      # saves posterior draws in "Output" file
      posterior_file = file.path("Output/Output_jointQIF",paste0(outroot,"__POSTERIOR.txt"))
      if(!file.exists(posterior_file)){ # regression for mitochondrial disease patients
        
      file.create(posterior_file)
      
      # dataframe for mean intensity of Control data 
      control = dat[(dat$fn%in%crl)&(dat$type=="Mean intensity"),]
      patient = dat[(dat$fn==pat)&(dat$type=="Mean intensity"),] 
      
      # control data
      Xctrl = log(control$value[control$channel==mitochan])
      Yctrl = log(control$value[control$channel==chan])
      XY_ctrl = cbind( Xctrl, Yctrl ) 
      Nctrl = nrow(XY_ctrl)
      
      # patient data
      
      Xpat = log(patient$value[patient$channel==mitochan])
      Ypat = log(patient$value[patient$channel==chan]) 
      XY_pat = cbind(Xpat, Ypat)
      Npat = nrow(XY_pat)
      
      # define prior parameters
      
      mu1_mean = c(1,1.5)
      mu2_mean = mu1_mean*1.5
      mu1_prec = solve( matrix(c(0.2,0.1,0.1,0.2), ncol=2, nrow=2, byrow=TRUE))
      mu2_prec = solve( 5*diag(2) )
      
      U_1 = matrix(c(10,7,7,10), ncol=2, nrow=2, byrow=TRUE)
      n_1 = 50
      U_2 = matrix(c(3,0,0,3), ncol=2, nrow=2, byrow=TRUE)
      n_2 = 20
      alpha = 1
      beta = 1
      
      data = list(Yctrl=XY_ctrl, Nctrl=Nctrl, Ypat=XY_pat, Npat=Npat,
                  mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                  mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                  U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)
      
      data_priorpred = data # same parameters used for prior prediction RJAGS code
      data_priorpred$Yctrl = NULL
      data_priorpred$Ypat = NULL
      # removes for prior prediction RJAGS 
      data_priorpred$Nctrl = 0
      data_priorpred$Npat = 0 
      
      op = par(mfrow=c(2,3), mar = c(5.5,5.5,3,3))
      
      
      model = jags.model(textConnection(modelstring), data=data, n.chains=n.chains) 
      
      model_priorpred = jags.model(textConnection(modelstring), data=data_priorpred) 
      
      update(model, n.iter=MCMCUpdates)
      
      converge = coda.samples(model=model, variable.names=c("mu","tau","Y_syn","z","probdiff"),
                                  n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)
      
      output = coda.samples(model=model, variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      output_priorpred = coda.samples(model=model_priorpred,
                                          variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                          n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      plot(output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                         "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                         "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                         "probdiff", "Y_syn[1]", "Y_syn[2]")] )
      
      posterior = as.data.frame(output[[1]])
      prior = as.data.frame(output_priorpred[[1]])
      
      class_posterior = posterior[, grepl('z', colnames(posterior))]
      classifs = colMeans(class_posterior)
      
      
      colnames(posterior) = colnames(output[[1]])
      colnames(prior) = colnames(output_priorpred[[1]])
      
      pdf(file.path("PDF/PDF_jointQIF", paste0(outroot,".pdf")),width=14,height=8.5)
      
      priorpost(data=data, prior=prior, posterior=posterior,
                classifs=classifs, title=paste(froot, pat) )
      
      dev.off()
      write.table(as.numeric(classifs),file.path("Output/Output_jointQIF",paste0(outroot,"__CLASS.txt")),row.names=FALSE,quote=FALSE,col.names=FALSE)
      write.table(posterior[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                   "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                   "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                   "probdiff", "Y_syn[1]", "Y_syn[2]")],posterior_file,row.names=FALSE,quote=FALSE)
      
      } else {
        
        class_file = file.path("Output/Output_jointQIF", paste0(outroot, "__CLASS.txt"))
        posterior = read.delim(posterior_file, sep=" ",stringsAsFactors=FALSE)
        colnames(posterior) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                    "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                    "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                    "probdiff", "Y_syn[1]", "Y_syn[2]")
      }
    } # pats
  } # chans
} # fulldats
