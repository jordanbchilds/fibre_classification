# Jags-MultivariateNormal.R
# John Kruschke, November 2015 - June 2017.
# For further info, see:
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

#--------------------------------------------------------------------------
# Load the data:
#--------------------------------------------------------------------------
alldat = read.delim("IFdata.txt",sep="\t",stringsAsFactors=FALSE) # must have file in curr. work. dir.
# Discard repeated occurences of controls
alldat = alldat[match(unique(alldat$cellid),alldat$cellid),]
# Discard zeros
alldat = alldat[(alldat$raw_CIV!=0)&(alldat$raw_porin!=0)&(alldat$raw_CI!=0),]
alldat = alldat[!is.na(alldat$sno),]

#args = commandArgs(trailingOnly=TRUE)
#if (length(args)!=0) {
#  stop("Exactly one patient ID (e.g. P11) required.n", call.=FALSE)
#} else {
#  pat = args[1]
#}
pat = "P11"

dat = alldat[alldat$sno==pat,]
datc = alldat[(alldat$controls=="control")&(alldat$Batch%in%unique(dat$Batch)),] 

# y must have named columns, with no missing values!
y = dat[,c("raw_porin","raw_CI","raw_CIV")]
colnames(y) = gsub("raw_","",colnames(y))
#y = y[sample(1:length(y[,1]),100,replace=FALSE),]

yc = datc[,c("raw_porin","raw_CI","raw_CIV")]
colnames(yc) = gsub("raw_","",colnames(yc))
if(length(yc[,1])>length(y[,1])){
 yc = yc[sample(1:length(yc[,1]),length(y[,1]),replace=FALSE),]
}

meany = apply(y,2,mean)
tauy = apply(y,2,function(x) 0.5)

meanyc = apply(yc,2,mean)
tauyc = apply(yc,2,function(x) 0.5)

#----------------------------------------------------------------------------
# The rest can remain unchanged, except for the specification of difference of
# correlations at the very end.
#----------------------------------------------------------------------------

# Load some functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.
# Install the ellipse package if not already:
want = c("ellipse")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

#Rscal = ncol(y)
#Rmat = diag(apply(yc,2,var))

Rscal = ncol(y)+2
Rmat = diag(5,3)
delta = row(Rmat) - col(Rmat)
Rmat[delta < 0 | delta > 0] = 0.5

# Assemble data for sending to JAGS:
dataList = list(
  y = y ,
  yc = yc,
  Ntotal =  nrow(y) ,
  Ntotalc =  nrow(yc) ,
  Nvar = ncol(y) ,
  # For prior on means
  mu_est = meany,
  tau_est = tauy,
  # For wishart (dwish) prior on inverse covariance matrix:
  Rscal = Rscal ,  # for dwish prior
  Rmat = Rmat,
  # Beta prior for theta, probability of being in mixture
  alpha_pi = 2,
  beta_pi = 2
)

library(mvtnorm)
#library(matlib)
pdf("PriorHeterogeneity.pdf",width=12,height=4)
for(jj in 1:20){

InvCovMat = rWishart(1,dataList$Rscal,dataList$Rmat)[,,1]
#CovMat = solve(InvCovMat) %*% InvCovMat 
CovMat = solve(InvCovMat)
mu = matrix(0,1,dataList$Nvar)
sigma = matrix(0,1,dataList$Nvar)
rho = matrix(0,dataList$Nvar,dataList$Nvar)
for(varIdx in 1:dataList$Nvar){
  mu[varIdx] = dnorm(1,dataList$mu_est[varIdx],1/sqrt(dataList$tau_est[varIdx]))
  sigma[varIdx] = sqrt(CovMat[varIdx,varIdx])
  for(varIdx2 in 1:dataList$Nvar){
     rho[varIdx,varIdx2] = CovMat[varIdx,varIdx2]/(sigma[varIdx]*sigma[varIdx2])
  }
}
print(rho) 
print(sigma)
print(mu)

# Preparation -- define useful functions:
library(ellipse)
expandRange = function( x , exMult=0.2 ) {
  lowVal = min(x)
  highVal = max(x)
  wid = max(x)-min(x)
  return( c( lowVal - exMult*wid , highVal + exMult*wid ) )
}

op=par(mfrow=c(1,3))
for ( varIdx1 in 1:(dataList$Nvar-1) ) {
  for ( varIdx2 in (varIdx1+1):dataList$Nvar ) {
    print(c(varIdx1,varIdx2))
    ellipseLevel = 0.90
    #plot( y[,c(varIdx1,varIdx2)] , type="n" ,
    #      xlim=c(0,max(y[,varIdx1])) , ylim=c(0,max(y[,varIdx2])) ,
    #      xlab=colnames(y)[varIdx1] , ylab=colnames(y)[varIdx2] ,
    #      main=bquote("Data with posterior "*.(ellipseLevel)*" level contour") )
    plot(ellipse(rho[varIdx2,varIdx1],scale = c(sigma[varIdx1],sigma[varIdx2]),centre=c(mu[varIdx1],mu[varIdx2]),level=ellipseLevel),
col="red",type="l",xlab=colnames(y)[varIdx1] , ylab=colnames(y)[varIdx2],xlim=c(-2,2),ylim=c(-2,2))
    points(yc[,c(varIdx1,varIdx2)]-meanyc[c(varIdx1,varIdx2)],pch=16,col=rgb(0,0,0,0.2))
    }
}
par(op)
}

dev.off()