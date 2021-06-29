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
myData = read.delim("IFdata.txt",sep="\t",stringsAsFactors=FALSE) # must have file in curr. work. dir.
# y must have named columns, with no missing values!
y = myData[myData$sno=="C01",c("raw_porin","raw_CI","raw_CIV")]
colnames(y) = gsub("raw_","",colnames(y))
meany = apply(y,2,mean)
tauy = apply(y,2,function(x) 0.5)

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

# Assemble data for sending to JAGS:
dataList = list(
  y = y ,
  Ntotal =  nrow(y) ,
  Nvar = ncol(y) ,
  # For prior on means
  mu_est = meany,
  tau_est = tauy,
  # For wishart (dwish) prior on inverse covariance matrix:
  Rscal = ncol(y) ,  # for dwish prior
  Rmat = diag(apply(y,2,var))
)

# Define the model:
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i,1:Nvar] ~ dmnorm( mu[1:Nvar] , InvCovMat[1:Nvar,1:Nvar] )
  }
  for ( varIdx in 1:Nvar ) { mu[varIdx] ~ dnorm( mu_est[varIdx] , tau_est[varIdx] ) }
  InvCovMat ~ dwish( Rmat[1:Nvar,1:Nvar] , Rscal )
  # Convert invCovMat to sd and correlation:
  CovMat <- inverse( InvCovMat )
  for ( varIdx in 1:Nvar ) { sigma[varIdx] <- sqrt(CovMat[varIdx,varIdx]) }
  for ( varIdx1 in 1:Nvar ) { for ( varIdx2 in 1:Nvar ) {
    rho[varIdx1,varIdx2] <- ( CovMat[varIdx1,varIdx2]
                               / (sigma[varIdx1]*sigma[varIdx2]) )
  } }
}
" # close quote for modelString
writeLines( modelString , con="Jags-MultivariateNormal-model.txt" )

# Run the chains:
nChain = 3
nAdapt = 500
nBurnIn = 500
nThin = 10
nStepToSave = 20000
require(rjags)
jagsModel = jags.model( file="Jags-MultivariateNormal-model.txt" ,
                        data=dataList , n.chains=nChain , n.adapt=nAdapt )
update( jagsModel , n.iter=nBurnIn )
codaSamples = coda.samples( jagsModel ,
                            variable.names=c("mu","sigma","rho") ,
                            n.iter=nStepToSave/nChain*nThin , thin=nThin )

# Convergence diagnostics:
parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=codaSamples , parName=parName )
}

# Examine the posterior distribution:
mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
Nvar = ncol(y)
# Create subsequence of steps through chain for plotting:
stepVec = floor(seq(1,chainLength,length=20))

# Make plots of posterior distribution:

# Preparation -- define useful functions:
library(ellipse)
expandRange = function( x , exMult=0.2 ) {
  lowVal = min(x)
  highVal = max(x)
  wid = max(x)-min(x)
  return( c( lowVal - exMult*wid , highVal + exMult*wid ) )
}

for ( varIdx in 1:Nvar ) {
  openGraph(width=7,height=3.5)
  par( mar=c(3.5,3,2,1) , mgp=c(2.0,0.7,0) )
  layout(matrix(1:2,nrow=1))
  # Marginal posterior on means:
  plotPost( mcmcMat[ , paste0("mu[",varIdx,"]") ] ,
            xlab=paste0("mu[",varIdx,"]") ,
            main=paste( "Mean of" , colnames(y)[varIdx] ) )
  # Marginal posterior on standard deviations:
  plotPost( mcmcMat[ , paste0("sigma[",varIdx,"]") ] ,
            xlab=paste0("sigma[",varIdx,"]") ,
            main=paste( "SD of" , colnames(y)[varIdx] ) )
}

for ( varIdx1 in 1:(Nvar-1) ) {
  for ( varIdx2 in (varIdx1+1):Nvar ) {
    openGraph(width=7,height=3.5)
    par( mar=c(3.5,3,2,1) , mgp=c(2.0,0.7,0) )
    layout(matrix(1:2,nrow=1))
    # Marginal posterior on correlation coefficient
    plotPost( mcmcMat[ , paste0("rho[",varIdx1,",",varIdx2,"]") ] ,
              xlab=paste0("rho[",varIdx1,",",varIdx2,"]") ,
              main=paste( "Corr. of" , colnames(y)[varIdx1] ,
                          "and" , colnames(y)[varIdx2] ) )
    # Data with posterior ellipse
    ellipseLevel = 0.90
    plot( y[,c(varIdx1,varIdx2)] , type="n" ,
          xlim=expandRange(y[,varIdx1]) , ylim=expandRange(y[,varIdx2]) ,
          xlab=colnames(y)[varIdx1] , ylab=colnames(y)[varIdx2] ,
          main=bquote("Data with posterior "*.(ellipseLevel)*" level contour") )
    # Posterior ellipses:
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[",varIdx1,"]") ,
                                         paste0("sigma[",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[",varIdx1,"]") ,
                                          paste0("mu[",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col="skyblue" , lwd=1 )
    }
    # replot data:
    points( y[,c(varIdx1,varIdx2)],pch=16,col=rgb(0,0,0,0.05) )
  }
}

# Show data descriptives on console:
cor( y )
apply(y,2,mean)
apply(y,2,sd)

##### Try to sketch prior
library(mvtnorm)
library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

Rscal = ncol(y)  # for Wishart prior
Rmat = diag(apply(y,2,var))

getsamp = function(){
 Q = rWishart(1,Rscal,Rmat)[,,1]
 sigma = solve(Q) %*% Q 
 return(rmvnorm(1,dnorm(meany,sqrt(1/tauy)),sigma)[1,])
}

prior_sim = t(replicate(1000,getsamp()))

k = kde2d(prior_sim[,"porin"],prior_sim[,"CI"],n=400)
image(k,col=r)
points(y[,"porin"],y[,"CI"])



