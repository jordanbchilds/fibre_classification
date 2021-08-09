library(rjags)
library(MASS)
library(twosamples)

source("../BootStrapping/parseData.R", local = TRUE)

myDarkGrey = rgb(169,169,159, max=255, alpha=20)
myGreen = rgb(25,90,0,max=255,alpha=50)
myYellow = rgb(225,200,50,max=255, alpha=50)
myBlue = rgb(0,0,255,max-255, alpha=50)

cramp = colorRamp(c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)),alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

plotWishart = function(output){
  if(length(dim(output))==3){
    par(mfrow=c(2,2))
    plot( output[1,1,], output[1,2,], pch=20, col=myDarkGrey, 
          xlab=expression(tau[11]), ylab=expression(tau[12]))
    plot( output[2,2,], output[1,2,], pch=20, col=myDarkGrey, 
          xlab=expression(tau[22]), ylab=expression(tau[12]))
    plot( output[1,1,], output[2,2,], pch=20, col=myDarkGrey, 
          xlab=expression(tau[11]), ylab=expression(tau[22]))
  } else {
    par(mfrow=c(2,2))
    plot( output[,"tau[1,1,1]"], output[,"tau[1,2,1]"], pch=20, col=myDarkGrey, 
          xlab=expression(tau[11]), ylab=expression(tau[12]))
    plot( output[,"tau[2,2,1]"], output[,"tau[1,2,1]"], pch=20, col=myDarkGrey, 
          xlab=expression(tau[22]), ylab=expression(tau[12]))
    plot( output[,"tau[1,1,1]"], output[,"tau[2,2,1]"], pch=20, col=myDarkGrey, 
          xlab=expression(tau[11]), ylab=expression(tau[22]))
  }
}

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

modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dbern(probdiff)
    class[i] =  2 - z[i]
    Y[i,1:2] ~ dmnorm(mu[1:2,class[i]], tau[1:2,1:2,class[i]] )
  }
  
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  p ~ dbeta(alpha_p, beta_p)
  probdiff = ifelse( pi==1, p, 1)
}
"
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

dir.create(file.path("Output/Output_tauPost"), showWarnings = FALSE)
dir.create(file.path("PDF/PDF_tauPost"), showWarnings = FALSE)
dir.create(file.path("PNG/PNG_tauPost"), showWarnings = FALSE)

# burn-in, chain length, thinning lag
MCMCUpdates = 2000
MCMCUpdates_Report = 10000
MCMCUpdates_Thin = 1

fulldat = 'IMC.RAW.txt'

imc_data = read.delim( file.path("../BootStrapping", fulldat), stringsAsFactors=FALSE)

mitochan = "VDAC1"
# imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imc_chan = "UqCRC2"

# removing unwanted info 
# imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
imcDat= imc_data[imc_data$channel %in% c(imc_chan, mitochan), ]

froot = gsub('.RAW.txt', '', fulldat)

sbj = sort(unique(imcDat$patient_id))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P", sbj, value = TRUE)
pts = c("P01")


# seperating the control patients 
for( chan in imc_chan ) {
  outroot_ctrl = paste( froot, 'CTRL', chan, sep='__')
  posterior_ctrl_file = file.path("Output/Output_IMC",paste0(outroot_ctrl,"__POSTERIOR.txt"))
  
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
    
    U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_1 = 2
    U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
    n_2 = 2
    alpha_p = 1
    beta_p = 1
    pi = 0
    
    # Bayesian inference using JAGS
    # Assume a prior centred around the true values (assume good estimate from control data)
    
    Nctrl = nrow(XY_ctrl)
    # parameter list for RJAGS
    data_ctrl = list(Y=XY_ctrl, N=Nctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                     mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                     U_1=U_1, U_2=U_2, alpha_p=alpha_p, beta_p=beta_p, pi=pi)
    
    data_ctrl_priorpred = data_ctrl # same parameters used for prior prediction RJAGS code
    data_ctrl_priorpred$Y = NULL # removes for prior prediction RJAGS 
    data_ctrl_priorpred$N = 0 # N: number of observed control points, removes for prior prediction
    
    # run the JAGS model for prior prediction and control parameter inference
    model_ctrl=jags.model(textConnection(modelstring), data=data_ctrl) # no initial vals given -> draw from prior
    
    model_ctrl_priorpred=jags.model(textConnection(modelstring), data=data_ctrl_priorpred) # no initial vals given -> draw from prior
    
    update(model_ctrl, n.iter=MCMCUpdates)
    
    output_ctrl=coda.samples(model=model_ctrl,variable.names=c("tau"),
                             n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    output_ctrl_priorpred=coda.samples(model=model_ctrl_priorpred,variable.names=c("tau"),
                                       n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])
    
    colnames(posterior_ctrl) = colnames(output_ctrl[[1]])
    colnames(prior_ctrl) = colnames(output_ctrl_priorpred[[1]])
    
  } 
  # define the expected value of the patient prior (prec_pred) be the mean of the control
  # posterior
  prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')]),
                      nrow=2, ncol=2, byrow=TRUE )
  
  prec_pred_inv = solve( prec_pred )
  
  var = var( posterior_ctrl[c('tau[1,1,1]', 'tau[1,2,1]','tau[2,2,1]')])
  d11 = var[1,1]/( 2*prec_pred_inv[1,1]^2 )
  d12 = var[2,2]/( prec_pred_inv[1,2]^2 + prec_pred_inv[1,1]*prec_pred_inv[2,2] )
  d22 = var[3,3]/( 2*prec_pred_inv[2,2]^2 )
  
  df_vec = seq(500, 650, by=20)
  KS_test11 = double(length(df_vec))
  KS_test22 = KS_test11
  KS_test12 = KS_test11
  
  # AD_test11 = double(length(df_vec))
  # AD_test22 = AD_test11
  # AD_test12 = AD_test11
  # 
  # CvM_test11 = double(length(df_vec))
  # CvM_test22 = CvM_test11
  # CvM_test12 = CvM_test11
  
  for( i in 1:length(df_vec)){
    wishart = rWishart(n=MCMCUpdates_Report, df=df_vec[i], Sigma=prec_pred/df_vec[i] )
    
    KS_test11[i] = ks.test(posterior_ctrl[,'tau[1,1,1]'], wishart[1,1,])$p.value
    KS_test22[i] = ks.test(posterior_ctrl[,'tau[2,2,1]'], wishart[2,2,])$p.value
    KS_test12[i] = ks.test(posterior_ctrl[,'tau[1,2,1]'], wishart[1,2,])$p.value
    
    # AD_test11[i] = ad_test(posterior_ctrl[,'tau[1,1,1]'], wishart[1,1,])[2]
    # AD_test12[i] = ad_test(posterior_ctrl[,'tau[1,2,1]'], wishart[1,2,])[2]
    # AD_test22[i] = ad_test(posterior_ctrl[,'tau[2,2,1]'], wishart[2,2,])[2]
    # 
    # CvM_test11[i] = cvm_test(posterior_ctrl[,'tau[1,1,1]'], wishart[1,1,])[2]
    # CvM_test12[i] = cvm_test(posterior_ctrl[,'tau[1,2,1]'], wishart[1,2,])[2]
    # CvM_test22[i] = cvm_test(posterior_ctrl[,'tau[2,2,1]'], wishart[2,2,])[2]
  }
  
  par(mfrow=c(3,1))
  plot(df_vec, KS_test11, type='b', col='black',
       xlab='Degrees of Freedom', ylab='p value', main=paste(chan))
  # lines(df_vec, AD_test11, type='b', col='blue')
  # lines(df_vec, CvM_test11, type='b', col='red')
  
  plot(df_vec, KS_test12, type='b', col='black',
       xlab='Degrees of Freedom', ylab='p value')
  # lines(df_vec, AD_test12, type='b', col='blue')
  # lines(df_vec, CvM_test12, type='b', col='red')
  
  plot(df_vec, KS_test22, type='b', 
       xlab='Degrees of Freedom', ylab='p value')
  # lines(df_vec, AD_test22, type='b', col='blue')
  # lines(df_vec, CvM_test22, type='b', col='red')
}

d=100
rWish = rWishart(n=10000, df=d, Sigma=prec_pred/d )
par(mfrow=c(3,1))
plot(density(rWish[1,1,]), main=expression(tau[11] ~ ' Density'))
lines(density(posterior_ctrl[,'tau[1,1,1]']), col='red')
plot(density(rWish[1,2,]), main=expression(tau[12]~' Density'))
lines(density(posterior_ctrl[,'tau[1,2,1]']), col='red') 
plot(density(rWish[2,2,]), main=expression(tau[22]~' Density'))
lines(density(posterior_ctrl[,'tau[2,2,1]']), col='red')

plotWishart(output=rWish)
plotWishart(output=posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')])





## NDUFB8 has an optimal d in (500,700), best tested was 560
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

## SDHA has an optimal d in (500,700), best tested was 560
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

## OSCP has an optimal d in (500,700), best tested was 560
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

## GRIM19 has an optimal d in (500,700), best tested was 560
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

##  MTCO1has an optimal d in (500,700), best tested was 580
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

##  COX4+4L2 has an optimal d in (500,700), best tested was 580
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar

##  UqCRC2 has an optimal d in (500,700), best tested was 580
### this however is not clean cut as one d can varying effects on different elements of T
### many values of d in this range have densities that look very similar





