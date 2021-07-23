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


priorpost = function(data, prior, posterior, classifs, ctrl=NULL,
                     output_mcmc=NULL, title, marginals=FALSE){
  # output: plots the prior and posterior regression lines and data
  
  
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  if( is.null(ctrl) ){
    x.lim = range(data[,1]) + c(-1,1)
    y.lim = range(data[,2]) + c(-1,1)
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    title(main=title, line = -1, outer = TRUE)
  } else {
    x.lim = range( data[,1], ctrl[,1] ) + c(-1,1)
    y.lim = range( data[,2], ctrl[,2] ) + c(-1,1)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=myYellow, pch=20)
    contour( kde2d(prior[,'Y_syn[1]'], prior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=classcols(classifs), pch=20)
    contour( kde2d(posterior[,'Y_syn[1]'], posterior[,'Y_syn[2]'], n=100), add=TRUE, nlevels=5 )
    title(main=title, line = -1, outer = TRUE)
  }
  
  par(op)
} 

MCMCplot = function( MCMCoutput, lag=20){
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
      plot( autocorr(MCMcoutput[[1]][,param], lags=0:lag), type='h', 
            xlab='Index', ylab='' )
      for(j in 2:n.chains) lines( autocorr(MCMCoutput[[j]][,param], lags=0:20), type='h', col=j)
      
      plot(1:nrow(MCMCoutput[[1]]), MCMCoutput[[1]][,param], main=param, type='l',
           ylab='', xlab='Iteration')
      for(j in 2:n.chains) lines(MCMCoutput[[j]][,param], type='l', col=j)
      
      plot(density(MCMCoutput[[1]][,param]) )
      for(j in 2:n.chains) lines(density(MCMCoutput[[j]][,param]), col=j )
    }
  }
  
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

dir.create(file.path("Output/QIF_allData"), showWarnings = FALSE)
dir.create(file.path("PDF/QIF_allData"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_allData"), showWarnings = FALSE)

dir.create(file.path("PDF/QIF_allData/MCMC_output"), showWarnings = FALSE)
dir.create(file.path("PDF/QIF_allData/classification"), showWarnings = FALSE)


# all datasets
fulldats = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

MCMCUpdates = 2000
MCMCUpdates_Report = 5000
MCMCUpdates_Thin = 1
n.chains = 3


# choose which datasets to estimate parameters for
# for(fulldat in fulldats){
# for(fulldat in c("IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt")){
for(fulldat in c("IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt")) {
  
  #fulldat = c("IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt")
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
    Xchan = vector()
    Ychan = vector()
    N = double(length(c(crl, pts)))
    
    for(i in 1:length(crl)){
      Xctrl = dat[(dat$fn%in%crl[i])&(dat$type=="Mean intensity")&(dat$channel==mitochan),'value']
      Yctrl = dat[(dat$fn%in%crl[i])&(dat$type=="Mean intensity")&(dat$channel==chan),'value']
      N[i] = length(Xctrl)
      
      Xchan = c(Xchan, Xctrl)
      Ychan = c(Ychan, Yctrl)
    }

    for( j in 1:length(pts) ){

        # dataframe for mean intensity of Control data 
      Xpat = dat[(dat$fn==pts[j])&(dat$type=="Mean intensity")&(dat$channel==mitochan),'value'] 
      Ypat = dat[(dat$fn==pts[j])&(dat$type=="Mean intensity")&(dat$channel==chan),'value'] 
        
      Xchan = c(Xchan, Xpat )
      Ychan = c(Ychan, Ypat )
      N[j+length(crl)] = length(Xpat)
    }
    
    pat_id = rep(c(crl, pts), N)
    data_chan = data.frame(log(Xchan), log(Ychan), pat_id)
    colnames(data_chan) = c(mitochan, chan, 'fn')
    
    is_control = 
    
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
    
    
  }
} # fulldats



