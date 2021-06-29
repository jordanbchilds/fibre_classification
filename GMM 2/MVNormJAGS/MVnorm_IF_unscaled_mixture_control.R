# Jags-MultivariateNormal.R
# John Kruschke, November 2015 - June 2017.
# For further info, see:
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# Optional generic preliminaries:
#graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Careful! This clears all of R's memory!

#--------------------------------------------------------------------------
# Load the data:
#--------------------------------------------------------------------------
alldat = read.delim("IFdata.txt",sep="\t",stringsAsFactors=FALSE) # must have file in curr. work. dir.
# Discard repeated occurences of controls
alldat = alldat[match(unique(alldat$cellid),alldat$cellid),]
# Discard zeros
alldat = alldat[(alldat$raw_CIV!=0)&(alldat$raw_porin!=0)&(alldat$raw_CI!=0),]
alldat = alldat[!is.na(alldat$sno),]

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Exactly one patient ID (e.g. P11) required.n", call.=FALSE)
} else {
  pat = args[1]
}
print(args)
#pat = "P19"

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


# Install the ellipse package if not already:
want = c("ellipse")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

#Rscal = ncol(y)
#Rmat = diag(apply(yc,2,var))

Rscal = ncol(y)+2
Rmat = diag(5,3)

#delta = row(Rmat) - col(Rmat)
#Rmat[delta < 0 | delta > 0] = 0.5

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

# Define the model:
modelString = "
model {
  pi ~ dbeta(alpha_pi, beta_pi)

  for(k in 1:2){
     InvCovMat[k,1:Nvar,1:Nvar] ~ dwish(Rmat[1:Nvar,1:Nvar] , Rscal)
     CovMat[k,1:Nvar,1:Nvar] <- inverse(InvCovMat[k,1:Nvar,1:Nvar])
     for ( varIdx in 1:Nvar ) { 
       mu[k,varIdx] ~ dnorm(mu_est[varIdx] , tau_est[varIdx])
       sigma[k,varIdx] <- sqrt(CovMat[k,varIdx,varIdx])
       for ( varIdx2 in 1:Nvar ) {
         rho[k,varIdx,varIdx2] <- (CovMat[k,varIdx,varIdx2]/(sigma[k,varIdx]*sigma[k,varIdx2]))
       }
     }
  }

  for ( i in 1:Ntotal ) {
    model_index[i] ~ dbern(pi) #ddirch & dcat for more than 2 clusters?
    y[i,1:Nvar] ~ dmnorm( mu[model_index[i]+1,1:Nvar] , InvCovMat[model_index[i]+1,1:Nvar,1:Nvar] )
  }

  for ( i in 1:Ntotalc ) {
    yc[i,1:Nvar] ~ dmnorm( mu[2,1:Nvar] , InvCovMat[2,1:Nvar,1:Nvar] )
  }
}
" # close quote for modelString
writeLines( modelString , con="Jags-MultivariateNormal-model.txt" )

# Run the chains:
nChain = 1
nAdapt = 5000
nBurnIn = 17500000
nThin = 100
nStepToSave = 5000
require(rjags)
jagsModel = jags.model( file="Jags-MultivariateNormal-model.txt" ,
                        data=dataList , n.chains=nChain , n.adapt=nAdapt )
update( jagsModel , n.iter=nBurnIn )

codaSamples = coda.samples( jagsModel ,
                            variable.names=c("model_index","mu","sigma","rho") ,
                            n.iter=nStepToSave/nChain*nThin , thin=nThin )

# Save output
save(codaSamples, file = paste(pat,".MCMC",sep="_"))
load(paste(pat,".MCMC",sep="_"))

source("DBDA2E-utilities.R") # Must be in R's current working directory.

# Convergence diagnostics:
parameterNames = varnames(codaSamples) # get all parameter names 
parameterNames = parameterNames[!grepl("model_index",parameterNames)] # too many indices to look at diagnostics for those
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
# Load some functions used below:

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
  plotPost( mcmcMat[ , paste0("mu[1,",varIdx,"]") ] ,
            xlab=paste0("mu[1,",varIdx,"]") ,
            main=paste( "Mean of" , colnames(y)[varIdx] ) )
  # Marginal posterior on standard deviations:
  plotPost( mcmcMat[ , paste0("sigma[1,",varIdx,"]") ] ,
            xlab=paste0("sigma[1,",varIdx,"]") ,
            main=paste( "SD of" , colnames(y)[varIdx] ) )
}

for ( varIdx1 in 1:(Nvar-1) ) {
  for ( varIdx2 in (varIdx1+1):Nvar ) {
    openGraph(width=7,height=3.5)
    par( mar=c(3.5,3,2,1) , mgp=c(2.0,0.7,0) )
    layout(matrix(1:2,nrow=1))
    # Marginal posterior on correlation coefficient
    plotPost( mcmcMat[ , paste0("rho[1,",varIdx1,",",varIdx2,"]") ] ,
              xlab=paste0("rho[1,",varIdx1,",",varIdx2,"]") ,
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
                                paste0("rho[1,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[1,",varIdx1,"]") ,
                                         paste0("sigma[1,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[1,",varIdx1,"]") ,
                                          paste0("mu[1,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col="skyblue" , lwd=1 )
    }
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[2,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[2,",varIdx1,"]") ,
                                         paste0("sigma[2,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[2,",varIdx1,"]") ,
                                          paste0("mu[2,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col="pink" , lwd=1 )
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

#prior_sim = t(replicate(1000,getsamp()))

#k = kde2d(prior_sim[,"porin"],prior_sim[,"CI"],n=400)
#image(k,col=r)
#points(y[,"porin"],y[,"CI"])

head(mcmcMat)
zerobiggest = (mcmcMat[,"mu[1,2]"]<=mcmcMat[,"mu[2,2]"])
mcmcMat[!zerobiggest,1:length(y[,1])] = as.numeric(xor(mcmcMat[!zerobiggest,1:length(y[,1])],1))
probbiggest = apply(mcmcMat[,1:length(y[,1])],2,sum)/dim(mcmcMat)[1]
y$probupper = probbiggest

cramp = colorRamp(c("red","yellow","blue"),space="Lab")
getcol = function(x,alph=1){
 chans = cramp(x)/255.0
 return(rgb(chans[1],chans[2],chans[3],alph))
}
y$col = sapply(y$probupper,getcol,0.3)

varids=1:3
names(varids)=c("porin","CI","CIV")

pdf(paste(pat,"report.pdf",sep="_"),width=8,height=8)

op = par(mfrow=c(2,2),mai=c(0.75,0.75,0.1,0.1))
for(ch in c("CI","CIV")){
    axrng = range(c(y$porin,yc$porin,y[[ch]],yc[[ch]]))
    plot(NULL,xlim=axrng,ylim=axrng,xlab="Porin",ylab=ch)
    varIdx1=1
    varIdx2=varids[ch]
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[1,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[1,",varIdx1,"]") ,
                                         paste0("sigma[1,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[1,",varIdx1,"]") ,
                                          paste0("mu[1,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col=rgb(1,0,0,0.3) , lwd=1 )
    }
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[2,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[2,",varIdx1,"]") ,
                                         paste0("sigma[2,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[2,",varIdx1,"]") ,
                                          paste0("mu[2,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col=rgb(0,0,1,0.3) , lwd=1 )
    }
    points(yc$porin,yc[[ch]],pch=16,col=rgb(0,0,0,0.1),cex=0.75)
    points(y$porin,y[[ch]],pch=16,col=y$col,cex=0.75)
}
    axrng = range(c(y$CI,yc$CI,y$CIV,yc$CIV))
    plot(NULL,xlim=axrng,ylim=axrng,xlab="CI",ylab="CIV")
    varIdx1=2
    varIdx2=3
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[1,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[1,",varIdx1,"]") ,
                                         paste0("sigma[1,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[1,",varIdx1,"]") ,
                                          paste0("mu[1,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col=rgb(1,0,0,0.3) , lwd=1 )
    }
    for ( stepIdx in stepVec ) {
      points( ellipse( mcmcMat[ stepIdx ,
                                paste0("rho[2,",varIdx1,",",varIdx2,"]") ] ,
                       scale=mcmcMat[ stepIdx ,
                                      c( paste0("sigma[2,",varIdx1,"]") ,
                                         paste0("sigma[2,",varIdx2,"]") ) ] ,
                       centre=mcmcMat[ stepIdx ,
                                       c( paste0("mu[2,",varIdx1,"]") ,
                                          paste0("mu[2,",varIdx2,"]") ) ] ,
                       level=ellipseLevel ) ,
              type="l" , col=rgb(0,0,1,0.3) , lwd=1 )
    }
    points(yc$CI,yc$CIV,pch=16,col=rgb(0,0,0,0.1),cex=0.75)
    points(y$CI,y$CIV,pch=16,col=y$col,cex=0.75)

    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, pat,cex = 7.5, col = "black", family="sans", font=2, adj=0.5)

par(op)

dev.off()


