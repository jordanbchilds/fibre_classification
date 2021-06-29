# Make some synthetic data, like patient with v-shaped 2Dmito scatterplot
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

# Fit linear model to "normal" cells
mod = lm(y~x, data=data.frame(x=x_normal,y=y_normal))

# Make some predictions for plotting, using regularly spaced x values
xsyn = seq(min(x),max(x),length.out=100)
pred_syn = predict(mod, newdata = data.frame(x=xsyn), se.fit=TRUE, interval = "prediction",na.action=na.omit, level=0.95)$fit
mid_syn = pred_syn[,1]
up_syn = pred_syn[,3]
low_syn = pred_syn[,2]

# Make some predictions based on the data
pred_dat = predict(mod, newdata=data.frame(x=x), se.fit=TRUE,  interval = "prediction",na.action=na.omit, level=0.95)$fit
mid_dat = pred_dat[,1]
up_dat = pred_dat[,3]
low_dat = pred_dat[,2]
z = y - mid_dat
#z = (y - mid_dat)/((up_dat-low_dat)/(2*qnorm(0.975)))

# Show scatterplot and z-score distribution
op = par(mfrow=c(2,2))
 plot(x,y,pch=16,col=rgb(1,0,0,0.2),cex.lab=1.55,cex.axis=2)
 points(xsyn,mid_syn,lwd=3,type="l")
 points(xsyn,up_syn,lwd=3,type="l",lty=2)
 points(xsyn,low_syn,lwd=3,type="l",lty=2)

 plot(x,y-mid_dat,pch=16,col=rgb(1,0,0,0.2),cex.lab=1.55,cex.axis=2)

 plot(density(z),main="z-scores",lwd=3,cex.lab=1.55,cex.main=2,cex.axis=2)
par(op)
