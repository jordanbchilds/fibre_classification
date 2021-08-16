library(MASS)
library(R2jags)

myDarkGrey = rgb(169,169,159, max=255, alpha=50)

modelstring = "
model{
  for(i in 1:N){
    Y[i,] ~ dmnorm(mu[,comp[i]], tau[,,comp[i]])
    comp[i] = 2 - z[i]
    z[i] ~ dbern(pi)
  }
  pi ~ dbeta(alpha,beta)
  # comp one
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  tau[1:2,1:2,1] ~ dwish(U_1, d_1)
  # comp two
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  tau[1:2,1:2,2] ~ dwish(U_2, d_2)
  # predective
  predOne ~ dmnorm(mu[,1], tau[,,1])
  predTwo ~ dmnorm(mu[,2], tau[,,2])
}
"

n1 = 1000
mu1 = c(3,3)
Sigma1 = matrix(c(5,5,5,5.5), nrow=2, ncol=2)
T1 = solve(Sigma1)

n2 = 1000
mu2 = c(5,0)
Sigma2 = matrix(c(4,0.3,0.3,0.3), nrow=2, ncol=2)
T2 = solve(Sigma2)

syn.1  = mvrnorm(n=n1, mu=mu1, Sigma=Sigma1)
xbar1 = colMeans(syn.1)
syn.2 = mvrnorm(n=n2, mu=mu2, Sigma=Sigma2)
xbar2 = colMeans(syn.2)

syn = rbind(syn.1, syn.2)

plot(syn, xlab='', ylab='', col=rep(c("blue","red"), c(n1,n2)))

MCMCBurnin = 2000
MCMCUpdate = 5000 + MCMCBurnin
MCMCThin = 1
n.chains = 1

mu1_mean = mu1
mu2_mean = mu2
mu1_prec = solve( 1*diag(2) )
mu2_prec = solve( 1*diag(2) )
d_1 = 50
U_1 = Sigma1*d_1
d_2 = 50
U_2 = Sigma2*d_2
alpha = 1
beta = 1

data = list(Y=syn, N=(n1+n2),
            mu1_mean=mu1_mean, mu1_prec=mu1_prec,
            mu2_mean=mu2_mean, mu2_prec=mu2_prec, d_1=d_1, d_2=d_2,
            U_1=U_1, U_2=U_2, alpha=alpha, beta=beta)

model = jags(data=data, parameters.to.save=c("mu","tau","pi","z","predOne","predTwo"),
             model.file=textConnection(modelstring), n.chains=n.chains, n.iter=MCMCUpdate, 
             n.thin=MCMCThin, n.burnin=MCMCBurnin, DIC=TRUE, progress.bar="text")

output = as.mcmc(model)

# traceplot(output[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
#                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
#                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
#                     "pi")])

posterior = output[[1]][,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                      "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                      "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                      "pi")]

postpred = output[[1]][,c("predOne[1]","predOne[2]","predTwo[1]","predTwo[2]")]

classifs_prob = colMeans( output[[1]][, grepl('z', colnames(posterior))] )
classifs = round(classifs_prob)

par(mfrow=c(1,2))
plot(syn, col=myDarkGrey, pch=20,
     xlab=expression(X[1]), ylab=expression(X[2]))
contour( kde2d(postpred[,"predOne[1]"], postpred[,"predOne[2]"], n=100),
         col="blue", lwd=2, add=TRUE)

plot(syn, col=myDarkGrey, pch=20,
     xlab=expression(X[1]), ylab=expression(X[2]))
contour( kde2d(postpred[,"predTwo[1]"], postpred[,"predTwo[2]"], n=100),
         col="red", lwd=2, add=TRUE)

data = cbind(syn, classifs)

# mu theoretical posterior parameters
S.1 = solve(mu1_prec+n1*T1)
M.1 = S.1%*%(mu1_prec%*%mu1_mean + n1*T1%*%xbar1)

mu1.rand = mvrnorm(n=10000, mu=M.1, Sigma=S.1 )
contour(kde2d(mu1.rand[,1], mu1.rand[,2], n=100), nlevels=5, col='blue', 
        main=expression(mu[1]), xlab=expression(mu[11]), ylab=expression(mu[21]),
        cex.main=2)
contour(kde2d(posterior[,"mu[1,1]"], posterior[,"mu[2,1]"], n=100), nlevels=5,
        col='red', add=TRUE)

S.2 = solve(mu2_prec+n2*T2)
M.2 = S.2%*%(mu2_prec%*%mu2_mean + n2*T2%*%xbar2)

mu2.rand = mvrnorm(n=10000, mu=M.2, Sigma=S.2 )
contour(kde2d(mu2.rand[,1], mu2.rand[,2], n=100), nlevels=5, col='blue',
        main=expression(mu[2]), xlab=expression(mu[12]), ylab=expression(mu[22]),
        cex.main=2)
contour(kde2d(posterior[,"mu[1,2]"], posterior[,"mu[2,2]"], n=100), nlevels=5,
        col='red', add=TRUE)

# tau theoretical posterior parameters
sum1 = 0 
for(i in 1:n1){ sum1 = sum1 + (syn.1[i,]-mu1)%*%t(syn.1[i,]-mu1) }
V.1 = sum1 + U_1
d.1 = n1 + d_1 + 2

T1.rand = rWishart(n=5000, Sigma=solve(V.1), df=d.1)
par(mfrow=c(1,3))
plot(density(T1.rand[1,1,]), xlab=expression(tau[11]), ylab='', 
     main=expression(T[1]), col="blue")
lines(density(posterior[,"tau[1,1,1]"]), col="red")
plot(density(T1.rand[1,2,]), xlab=expression(tau[12]), ylab='', 
     main=expression(T[1]), col="blue")
lines(density(posterior[,"tau[1,2,1]"]), col="red")
plot(density(T1.rand[2,2,]), xlab=expression(tau[22]), ylab='', 
     main=expression(T[1]), col="blue")
lines(density(posterior[,"tau[2,2,1]"]), col="red")

sum2 = 0 
for(i in 1:n2){ sum2 = sum2 + (syn.2[i,]-mu2)%*%t(syn.2[i,]-mu2) }
V.2 = sum2 + U_2
d.2 = n2 + d_2 + 2

T2.rand = rWishart(n=5000, Sigma=solve(V.2), df=d.2)
par(mfrow=c(1,3))
plot(density(T2.rand[1,1,]), xlab=expression(tau[11]), ylab='', 
     main=expression(T[2]), col="blue")
lines(density(posterior[,"tau[1,1,2]"]), col="red")
plot(density(T2.rand[1,2,]), xlab=expression(tau[12]), ylab='', 
     main=expression(T[2]), col="blue")
lines(density(posterior[,"tau[1,2,2]"]), col="red")
plot(density(T2.rand[2,2,]), xlab=expression(tau[22]), ylab='', 
     main=expression(T[2]), col="blue")
lines(density(posterior[,"tau[2,2,2]"]), col="red")


# pi theoretical posterior pararmeters
alpha.1 = n1+alpha
beta.1 = n2+beta

par(mfrow=c(1,1))
pi.rand = rbeta(n=5000, alpha.1, beta.1)
plot(density(pi.rand), xlab='p', ylab='', col="blue", xlim=c(0,1))
lines(density(posterior[,"pi"]), col="red")











