cramp = colorRamp(c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)),alpha=TRUE)
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

DarkGrey = rgb(169,169,159, max=255, alpha=50)

plot_Wishart = function( array ){
  N = dim( array )
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

### SYNTHETIC DATA ####
N = 600
mu_x1 = 5
mu_x2 = 10
sd_x = 0.9
intercept = 0.1
slope_normal = 1.0
slope_deficient = 0.0
sd_err_norm = 0.8
sd_err_def = 2.5

x_normal_g1 = rnorm(N/4, mu_x1, sd_x)
x_normal_g2 = rnorm(N/4, mu_x2, sd_x)

x_deficient_g1 = seq(mu_x1-2*sd_x, mu_x1+2*sd_x, length.out=N/4)
x_deficient_g2 = seq(mu_x2-2*sd_x, mu_x2+2*sd_x, length.out=N/4)

y_normal_g1 = slope_normal*x_normal_g1 + intercept + rnorm(N/4,0,sd_err_norm)
y_normal_g2 = slope_normal*x_normal_g2 + intercept + rnorm(N/4,0,sd_err_norm)

y_deficient_g1 = slope_deficient*x_deficient_g1 + intercept + rnorm(N/4,0,sd_err_def)
y_deficient_g2 = slope_deficient*x_deficient_g2 + intercept + rnorm(N/4,0,sd_err_def)

x = c(x_normal_g1, x_normal_g2, x_deficient_g1, x_deficient_g2)
y = c(y_normal_g1, y_normal_g2, y_deficient_g1, y_deficient_g2)
Y = cbind(x,y)

par(mfrow=c(1,1))
plot(x,y, col=c(rep('black', N/2), rep('red',N/2)), pch=20,
     main='Simulated Data')



###--------------------------------------------------
###### MODEL MODEL MODEL MODEL #######
###--------------------------------------------------

modelstring = "
model {
  for(i in 1:N){
    z[i] ~ dbern(pi)
    class[i] =  z[i] + 1
    Y[i,] ~ dmnorm(mu[,class[i]], tau[,,class[i]] )
  }
  
  # construsting covariance matrix for group 1
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  pi~ dbeta(alpha_pi, beta_pi)
  
  z_syn ~ dbern(pi)
  class_syn = z_syn + 1
  Y_syn ~ dmnorm(mu[,class_syn], tau[,,class_syn])
}
"

#### PRIOR DENSITIES ####
par(mfrow=c(1,1))
mu_rand = mvrnorm(n=10000, mu=c(0,0), Sigma=solve(0.25*diag(2)))
plot(mu_rand[,1], mu_rand[,2], col=DarkGrey, pch=20, xlab=expression(mu[11]),
     ylab=expression(mu[12]), main=expression('Prior'~mu[1]))
contour(kde2d(mu_rand[,1],mu_rand[,2], n=100), add=TRUE, nlevels=15)

mu_rand = mvrnorm(n=10000, mu=c(0,0), Sigma=solve(0.25*diag(2)))
plot(mu_rand[,1], mu_rand[,2], col=DarkGrey, pch=20, xlab=expression(mu[21]),
     ylab=expression(mu[22]), main=expression('Prior'~mu[2]))
contour(kde2d(mu_rand[,1],mu_rand[,2], n=100), add=TRUE, nlevels=15)

prec1_rand = rWishart(10000, df=2, solve(matrix(c(1.1,1,1,1.1),ncol=2,nrow=2,byrow=TRUE) ) )
plot_wishart(prec1_rand)

prec2_rand = rWishart(10000, df=2, solve(matrix(c(1,0,0,2),ncol=2,nrow=2,byrow=TRUE) ) )
plot_wishart(prec2_rand)


#### INFERENCE ####
data_input = list(Y=Y, N=N, mu1_mean=c(0,0), mu1_prec=0.25*diag(2), 
                  mu2_mean=c(0,0), mu2_prec=0.25*diag(2), n_1=2, n_2=2,
                  U_1=solve(matrix(c(2,0,0,2)/2,ncol=2,nrow=2,byrow=TRUE)),
                  U_2=solve(matrix(c(2,0,0,2)/2,ncol=2,nrow=2,byrow=TRUE)),
                  alpha_pi=2, beta_pi=2)

MCMCUpdates = 5000
# 10,000 posterior draws after burn-in
MCMCUpdates_Report = 5000
# thin with lag-10 - > 1,000 draws from posterior
MCMCUpdates_Thin = 1

model = jags.model(textConnection(modelstring), data=data_input, n.chains=3)

update(model, n.iter=MCMCUpdates)

output = coda.samples(model=model,variable.names=c('z', 'pi', 'mu', 'tau', 'z_syn', 'Y_syn'),
                      n.iter=MCMCUpdates_Report, thin=MCMCUpdates_Thin)

plot(output[,c('mu[1,1]', 'mu[1,2]', 'mu[2,1]','mu[2,2]')])
plot(output[,c('tau[1,1,1]','tau[1,2,1]','tau[2,2,1]')])
plot(output[,c('tau[1,1,2]','tau[1,2,2]','tau[2,2,2]')])

output_comb = rbind(output[[1]], output[[2]], output[[3]])

class_posterior = output_comb[,16:515]

class_mean = colMeans(class_posterior)
class_prob = class_mean 

y_syn = output_comb[,c('Y_syn[1]','Y_syn[2]')]

par(mfrow=c(1,1))
plot(x,y, col=classcols(class_prob), pch=19)
contour( kde2d(y_syn[,1], y_syn[,2], n=100), add=TRUE, nlevels=15)


#################################################################
#########              WISHART DISTRIBUTION 
#################################################################
N_d = 10000

## changing the degrees of freedom

n2 = rWishart(n=N_d, df=2, diag(2))
plot_Wishart(n2)

n5 = rWishart(n=N_d, df=5, diag(2))
plot_Wishart(n5)

n20 = rWishart(n=N_d, df=20, diag(2))
plot_Wishart(n20)

sig_0.5 = rWishart(n=N_d, df=2, 0.5*diag(2))
plot_Wishart(sig_0.5)

sig_1.5 = rWishart(n=N_d, df=2, 1.5*diag(2))
plot_Wishart(sig_1.5)

n5_sig0.5 = rWishart(n=N_d, df=5, 05*diag(2))
plot_Wishart(n5_sig0.5)

#################################################################
#########         COVARIANCE vs. PRECISION MATRICES
#################################################################

prec = matrix(c(2.532165283,	-0.061085756,	-0.061085756,	0.885847557),
              ncol=2, nrow=2, byrow=TRUE)
delta = matrix(c(0,0.1,0.1,0), ncol=2, nrow=2, byrow=TRUE)

cov = solve(prec)
solve(prec+delta)

solve(cov+delta)



