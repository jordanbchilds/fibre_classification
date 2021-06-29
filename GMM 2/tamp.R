library(rjags)
library(beanplot)
library(MASS)

DarkGrey = rgb(169,169,159, max=255, alpha=50)

plot_Wishart = function( array, inv=FALSE ){
  N = dim( array )
  if(inv){
    for(i in 1:N[3]){
      array[,,i] = solve( array[,,i] )
    }
  }
  sig11 = array[1,1,]
  sig22 = array[2,2,]
  sig12 = array[1,2,]
  
  var_lim = c(0, max(c(sig11, sig22) ) )
  cov_lim = c( min(sig12), max(sig12) ) 
  
  
  den_11_22 = kde2d(sig11, sig22, n=100)
  den_11_12 = kde2d(sig11, sig12, n=100)
  den_22_12 = kde2d(sig22, sig12, n=100)
  
  par(mfrow=c(2,2))
  plot(sig11, sig22, pch=20, col=DarkGrey, main=expression(sigma[11]~'vs.'~sigma[22]),
       xlab=expression(sigma[11]), ylab=expression(sigma[22]), ylim=var_lim, xlim=var_lim)
  contour(den_11_22, add=TRUE, nlevels=20)
  
  plot(sig11, sig12, pch=20, col=DarkGrey, main=expression(sigma[11]~'vs.'~sigma[12]),
       xlab=expression(sigma[11]), ylab=expression(sigma[12]), ylim=cov_lim, xlim=var_lim)
  contour(den_11_12, add=TRUE, nlevels=20)
  
  plot(sig22, sig12, pch=20, col=DarkGrey, main=expression(sigma[22]~'vs.'~sigma[12]),
       xlab=expression(sigma[22]), ylab=expression(sigma[12]), ylim=cov_lim, xlim=var_lim)
  contour(den_22_12, add=TRUE, nlevels=20)
  
  par(mfrow=c(1,1))
}

# modelstring_temp = "
# model {
#   for(i in 1:N){
#     z[i] ~ dbern(probdiff)
#     class[i] =  z[i] + 1
#     Y[i,] ~ dmnorm(mu[,class[i]], tau[,,class[i]] )
#   }
#   
#   # construsting covariance matrix for group 1
#   tau[1:2,1:2,1] ~ dwish(U_1, n_1)
#   mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
#   
#   tau[1:2,1:2,2] ~ dwish(U_2, n_2)
#   mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
#   
#   # classification
#   p ~ dbeta(alpha_p, beta_p)
#   probdiff = ifelse( pi==1, p, 1)
#   
#   # posterior distribution
#   z_syn ~ dbern(probdiff)
#   class_syn = z_syn + 1
#   Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
# }
# "

modelstring_temp = "
model {
  p ~ dbeta(alpha_p, beta_p)
  
  z ~ dbern(p)
  class =  z + 1
  Y ~ dmnorm( mu[,class], tau[,,class] )
  
  mu1 ~ dmnorm(mu1_mean, mu1_prec)
  mu2 ~ dmnorm(mu2_mean, mu2_prec)
  mu[1:2,1] = mu1
  mu[1:2,2] = mu2
  
  tau1 ~ dwish( U_1, df_1 )
  tau2~ dwish( U_2, df_2 )
  tau[1:2,1:2,1] = tau1
  tau[1:2,1:2,2] = tau2
}
"
n_1 = 8 # degrees of freedom
# define the expected value of the patient prior (prec_pred) be the mean of the control
# posterior
prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[2,1,1]','tau[2,2,1]')]),
                    nrow=2, ncol=2, byrow=TRUE)
# increase the covariance between 'x' and 'y', keep variances the same
Sigma = solve(prec_pred)
delta_vec = c(2,-1.5,-1.5,2)
delta = matrix(delta_vec, nrow=2, ncol=2, byrow=TRUE)
# delta is positive definite
# if M and N are PD then M+N is also PD
# (the same holds for M,N PSD)

# 
# # re-define the expectation of the prior
prec_pred = solve( Sigma + delta )
# define prior parameter
U_1 = prec_pred/n_1
n_2 = 6
U_2 = solve( matrix(c(2,0,0,2),nrow=2,ncol=2,byrow=TRUE) )/n_2

plot_Wishart( rWishart(10000, df=n_1, Sigma=U_1), inv=TRUE )
plot_Wishart( rWishart(10000, df=n_2, Sigma=U_2), inv=TRUE )

mu1_mean = colMeans( posterior_ctrl[,c('mu[1,2]','mu[2,2]')])
mu1_prec = solve( var( posterior_ctrl[,c('mu[1,2]','mu[2,2]')] ) )

mu2_mean = mu1_mean/2
mu2_prec = solve( matrix( c(4,0,0,4), nrow=2, ncol=2, byrow=TRUE))

alpha_p = 2
beta_p = 2
pi = 1

data_temp = list(mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                mu2_mean=mu2_mean, mu2_prec=mu2_prec, df_1=n_1, df_2=n_1,
                U_1=U_1, U_2=U_1, alpha_p=alpha_p, beta_p=beta_p)

model_temp = jags.model( textConnection(modelstring_temp), data=data_temp) 

update(model_temp, n.iter=MCMCUpdates)

output_temp = coda.samples( model=model_temp, variable.names=c('Y'),
                         n.iter=100000, thin=MCMCUpdates_Thin)

plot(output_temp[,c('Y[1]','Y[2]')])

output_temp_df = as.data.frame(output_temp[[1]])

par(mfrow=c(1,1))
contour( kde2d(output_temp_df[,'Y[1]'], output_temp_df[,'Y[2]'], n=100),
         nlevels=20, xlab=paste('n_1=',n_1,', n_2=',n_2), main=paste('delta=',delta_vec) )



















