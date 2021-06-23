modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dcat(q[])
    Y[i,] ~ dmnorm(mu[,z[i]], tau[,(2*z[i]-1):(2*z[i])] )
  }
  
  # construsting covariance matrix for group 1
  for(i in 1:3){ taus_1[i] ~ dgamma(alpha_tau, beta_tau) }
  taus_1s = sort(taus_1)
  tau[1,2] = taus_1s[1]
  tau[2,1] = tau[1,2]
  p1 ~ dunif(0,1)
  tau[1,1] = ifelse( p1<=0.5, taus_1s[2], taus_1s[3] )
  tau[2,2] = ifelse( p1>0.5, taus_1s[3], taus_1s[2])
  
  # constructing covarinace matrix fro group 2
  for(i in 1:3){ taus_2[i] ~ dgamma(alpha_tau, beta_tau) }
  taus_2s = sort(taus_2)
  tau[1,4] = taus_2s[1]
  tau[2,3] = tau[1,2]
  p2 ~ dunif(0,1)
  tau[1,3] = ifelse( p2<=0.5, taus_2s[2], taus_2s[3] )
  tau[2,4] = ifelse( p2>0.5, taus_2s[3], taus_2s[2])
  
  # constructing means 
  mu_1 ~ dmnorm( mu_mean, mu_prec)
  mu_2 ~ dmnorm( mu_mean, mu_prec)
  mu[1:2,1] = mu_1
  mu[1:2,2] = mu_2
  
  # classification
  pi ~ dbeta(alpha_pi, beta_pi)
  q[1] = pi
  q[2] = 1-pi

}
"
### SYNTHETIC DATA ####
N = 500
mu_x = 5
sd_x = 2.5
intercept = 0.1
slope_normal = 1.0
slope_deficient = 0.0
sd_err = 1.6

x_normal = rnorm(N/2,mu_x,sd_x)
#x_deficient = runif(N/2, mu_x-2*sd_X, mu_X+2*sd_x)
x_deficient = seq(mu_x-2*sd_x, mu_x+2*sd_x, length.out=N/2)
y_normal = slope_normal*x_normal + intercept + rnorm(N/2,0,sd_err)
y_deficient = slope_deficient*x_deficient + intercept + rnorm(N/2,0,sd_err)
#x = c(x_normal,y_normal)
x = c(x_normal, x_deficient)
y = c(y_normal, y_deficient)
Y = cbind(x,y)

#### PRIOR DENSITIES ####
par(mfrow=c(1,1))
pi_supp = seq(0,1,length.out=200)
pi_den = dbeta(pi_supp, 2,2)
plot(pi_supp, pi_den, type='l', lwd=2, xlab=expression(pi), ylab='Density', main=expression('Prior Density, '~ pi))

tau_supp = seq(0,10, length.out=200)
tau_den = dgamma(tau_supp, 3,1)
plot(tau_supp, tau_den, type='l', lwd=2, xlab=expression(tau), ylab='Density', main=expression('Prior Density, '~ tau))


#### INFERENCE ####
data_input = list(Y=Y, N=N, mu_mean=c(0,0), mu_prec=0.25*diag(2), alpha_tau=3, beta_tau=1, alpha_pi=2, beta_pi=2)

MCMCUpdates = 5000
# 10,000 posterior draws after burn-in
MCMCUpdates_Report = 5000
# thin with lag-10 - > 1,000 draws from posterior
MCMCUpdates_Thin = 0

model_0=jags.model(textConnection(modelstring),data=data_input, n.chains=3)

update(model_0, n.iter=MCMCUpdates)

output_0=coda.samples(model=model_0,variable.names=c('tau_1', 'tau_2','mu_1','mu_1', 'z', 'pi'),
                      n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)





