setwd("../MOBILAB/Statistik-Hilfe GIUB/Markus_VulnerabilitätFahrhabe")
load("Loss/Loss_Data.RData")
# contains 384 by 5 dataframes "Contents" and "Structure"

source("Global_LossAnalysisLF.R")

plot( Structure$Loss, Contents$Loss, xlab = "Structure", ylab = "Content", main = "Absolute Loss" )
points( Structure$DoL[c(145,249,328)], Contents$DoL[c(145,249,328)], pch = 4, col = c('red','red','blue') )
abline( 0, 1, lty = 4, col = 'grey' )

##########
#==== Fit Box-Cox transformation
##########

p <- 2

# Execute only one of the three possibilities.
#==== All data points ====#
y.Loss <- Contents$Loss
x.Loss <- Structure$Loss

modelname <- "Full"
#==== End ====#


#==== W/o outliers ====#
idx.outl <- 145

y.Loss <- Contents$Loss[-idx.outl]
x.Loss <- Structure$Loss[-idx.outl]

modelname <- "noOutl"

PTBS.mle( y = Contents$Loss, x = Structure$Loss, interval = c(-3,1) )  # as comparison
#==== End ====#


#==== W/o outliers nor lev. pts ====#
idx.remove <- c(145,249) # based on complete data set

y.Loss <- Contents$Loss[-idx.remove]
x.Loss <- Structure$Loss[-idx.remove]

modelname <- "notOutLev"

PTBS.mle( y = Contents$Loss, x = Structure$Loss, interval = c(-3,1) ) # comparison
PTBS.mle( y = Contents$Loss[-idx.outl], x = Structure$Loss[-idx.outl], interval = c(-3,1) ) # comparison
#==== End ====#

n.obs <- length(y.Loss)

sep.lam.global <- TRUE
if(sep.lam.global) { modelname <- paste( modelname, "_seplam", sep = "" ) }

PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = sep.lam.global )

CI.lambda.Loss <- PTBS.CI.lambda( y = y.Loss, x = x.Loss, interval = c(-3,1), sep.lam = sep.lam.global )
CI.lambda.Loss
#==# Full: CI = (0.0556,0.1814)
#==# W/o outliers: CI = ()
#==# W/o extr.Cdist: CI = ()

mle.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = sep.lam.global )
par.Loss <- unlist(mle.Loss[c(2,3,1)])
s2.hat.Loss <- mle.Loss$sigma2.hat * n.obs/(n.obs-p)
# Q: not clear whether n-p or n-p-1 ?! 
#--> I would say n-p because the rank of I - H.lambda is still n-p even though lambda is estimated
covML.Loss <- solve( -optimHess( par.Loss, fn = PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ) ) )


# Estimation results
cbind( par.Loss, sqrt(diag(covML.Loss)) )
mle.Loss$ell.opt
2*( -mle.Loss$ell.opt + p+2+sep.lam.global )

# p-value LRT seplam vs. single lambda (in case of sep.lam.global = TRUE)
1 - pchisq( 2*( mle.Loss$ell.opt - PTBS.mle( y = y.Loss, x = x.Loss, sep.lam = FALSE )$ell.opt ), df = 1 )


# Transformed quantities
y.lambda.Loss <- BC.transform( mle.Loss$lambda.hat[1], y.Loss )
x.lambda.Loss <- BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], x.Loss )
X.lambda.Loss <- cbind( 1, x.lambda.Loss )

## Covariance matrix of beta based on regression (#3 below):
covbeta.regr.Loss <- s2.hat.Loss * solve(crossprod(X.lambda.Loss))  # much smaller than ML ones
# var(sigma2.hat) = 2 sigma2.hat^2/n, var(s2.hat) = 2 s2.hat^2/n

# theta.init.Loss <- if(sep.lam.global) { c( 7.5, 0.2, 2.9, 0.15, 0.2 ) } else { c( 2, 1, 4, 0.1 ) + c(384-n.obs,0,0,0) } #EVENTUELL
# theta.init.Loss <- if(sep.lam.global) { c( 8, 0.05, 3.1, 0.2, 0.38 ) } else { c( 2, 1, 4, 0.1 ) + c(384-n.obs,0,0,0) } # w/o outlier
# opt.Loss <- optim( theta.init.Loss, fn = PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ), hessian = TRUE ) # ev. 3 for sigma2.init
# #opt.Loss <- optim( c( theta.init.Loss, if(sep.lam.global){ 0.1 } else { NULL } ), fn = PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ), hessian = TRUE )
# opt.Loss$val - mle.Loss$ell.opt
# covML.Loss <- solve( -opt.Loss$hessian )

## !! sigma2.hat is now correlated because of the estimated lambda:
#cov2cor( covML.Loss )

## Bootstrap covariance matrix (not needed)
# mle.boot.Loss <- t(apply( boot.idx, 1, function(star){ unlist( PTBS.mle( y = y.Loss[star], x = x.Loss[star], interval = range(lam.seq.Loss) )[c(2,3,1)] ) } )) # boot.idx is R by n.obs
# covBoot.Loss <- var( mle.boot.Loss )


###### Covariance relations (can be skipped)
opt.lamfix.Loss <- optim( theta.init.Loss[1:3], fn = regr.llkhd, y = y.lambda.Loss, X = X.lambda.Loss, control = list( fnscale = -1 ), hessian = TRUE )
# Estimating sigma instead of sigma2 (via sigma = TRUE in regr.llkhd) is not so helpful bcs true distn of sigma.hat is only half-normal and not normal (and therefore certainly not equal the asymptotic one).
sqrt(diag(covML.Loss[1:3,1:3]/solve( -opt.lamfix.Loss$hessian ))) # increase factor of se

# J(theta.hat)^{-1} vs. regression (sigma2.hat * (X^T X)^{-1})
solve( -opt.lamfix.Loss$hessian )[1:2,1:2]/( opt.lamfix.Loss$par[3] * solve( crossprod( X.lambda.Loss ) ) )
# Under lamfix: covML.beta = sigma2.hat * (X^T X)^{-1}  (thus w/ biased sigma2.hat)

# estimated lambda vs. fixed lambda (both J(theta.hat)^{-1}):
covML.Loss[1:3,1:3]/solve( -opt.lamfix.Loss$hessian ) 
# --> with estimated lambda much larger variance of sigma2
## having a different lambda for x and y largely increases the variance of beta.hat !!

# cov(beta.hat): s2.hat * (X^T X)^{-1} (regression) vs. J(theta.hat)^{-1} for lambda fixed
covbeta.regr.Loss/solve( -opt.lamfix.Loss$hessian )[1:2,1:2]  # roughly n.obs/(n.obs-p)
# Thus J(theta.hat)^{-1}[1:2,1:2] = sigma2*(X^t X)^{-1} equals the exact covariance of beta.hat in the Gaussian case.
# Thus covbeta.regr = s2.hat * (X^t X)^{-1} equals n.obs/(n.obs-p) * J(theta.hat)^{-1}[1:2,1:2]
# --> Use covbeta.MLadj (?)

# cov(beta.hat) for lambda fixed: J(theta.hat)^{-1} vs. regression version for theta.hat(opt):
solve( -opt.lamfix.Loss$hessian )[1:2,1:2]/( n.obs/(n.obs-p) * opt.lamfix.Loss$par[3] * solve(crossprod( cbind( 1, BC.transform( opt.Loss$par[4+sep.lam.global], x.Loss ) ) )) )

opt.lamfix.Loss$par[3] * solve( crossprod( cbind( 1, BC.transform( opt.Loss$par[4+sep.lam.global], x.Loss ) ) ) )/solve( -opt.lamfix.Loss$hessian )[1:2,1:2]
# --> sigma2.hat(lamfix) {X(lamfix)^t X(lamfix)}^{-1} = J(lamfix)^{-1}
######

#==== End Fit ====#


##===
# Plot confidence region for lambda
##===

##===== same lambda =====##
proflambda.Loss <- sapply( lam.seq.Loss, PTBS.profllkhd.lambda, y = y.Loss, x = x.Loss, intercept = TRUE )

## Plot profile log-likelihood for lambda
plot( lam.seq.Loss, proflambda.Loss - CI.lambda.Loss[[4]], xlab = expression(lambda), ylab = "Profile log likelihood", type = 'l', xlim = c(-1,1), ylim = c(-200,0) )
abline( h = - qchisq(0.95,1)/2, lty = 2 )
abline( v = CI.lambda.Loss[1:3], lty = 3 )
##===== End =====##


##===== separate lambdas =====##
# Profile likelihood array
lam1.seq.Loss <- seq( -0.07, 0.25, 0.001 )
lam2.seq.Loss <- seq( 0.1, 0.5, 0.001 )
lam.arr.Loss <- expand.grid( lam1 = lam1.seq.Loss, lam2 = lam2.seq.Loss )

proflambda.Loss <- matrix( apply( lam.arr.Loss, 1, PTBS.profllkhd.lambda, y = y.Loss, x = x.Loss ), nr = length(lam1.seq.Loss) )

## Plot confidence region for (lambda.y, lambda.x)
image( x = lam1.seq.Loss, y = lam2.seq.Loss, pmax( proflambda.Loss - (mle.Loss$ell.opt - qchisq(0.95, df=2)/2), 0 ), xlab = expression(lambda[y]), ylab = expression(lambda[x]), breaks = c(-0.0001,0.0001,seq( 0.05, 3, 0.05 )), col = c('white',topo.colors(60)), main = "Absolute loss, 95% CR lambdas" )
abline( v = mle.Loss$lambda.hat[1], col = 'red' )
abline( h = mle.Loss$lambda.hat[2], col = 'red' )
abline( v = unlist(CI.lambda.Loss)[c(3,5)], col = 'orange', lty = 2 )
abline( h = unlist(CI.lambda.Loss)[c(4,6)], col = 'orange', lty = 2 )
# abline( v = mle.Loss$lambda.hat[1] + c(1,-1)*qnorm(0.975)*sqrt(diag(covML.Loss)[4]), col = 'red', lty = 3 )
# abline( h = mle.Loss$lambda.hat[2] + c(1,-1)*qnorm(0.975)*sqrt(diag(covML.Loss)[5]), col = 'red', lty = 3 )
abline( 0, 1, col = 'darkgrey', lty = 2 )
box()
savePlot( paste( "FiguresLF/Loss_", modelname, "_ConfReg", sep = '' ), type = "pdf" )
# identity line NOT through CR --> lambdas are different
##===== End =====##



##===================##
##=== Diagnostics ===##
##===================##

## Involved quantities
H.lambda.Loss <- X.lambda.Loss %*% solve( crossprod( X.lambda.Loss ) ) %*% t(X.lambda.Loss)
#sum( diag( H.lambda.Loss ) ) # tr(H) = p = 2
# Fitted values
yhat.lambda.Loss <- c( H.lambda.Loss %*% y.lambda.Loss )
# Residuals:
e.Loss <- y.lambda.Loss - yhat.lambda.Loss
# Standardised residuals:
r.Loss <- e.Loss/sqrt(s2.hat.Loss*(1-diag(H.lambda.Loss)))
# Deletion residuals (Davison 2008, p.395) = "studentized residuals" in Weisberg 2005 (p.196):
del.resid.Loss <- sqrt( (n.obs-p-1)/(n.obs-p-r.Loss^2) )*r.Loss

## Outliers?
idx.tail <- rev( tail( order( abs(del.resid.Loss) ) ) )
2*( 1 - pt( abs(del.resid.Loss[idx.tail[1]]), df = n.obs - p - 1 ) )
0.05/n.obs  # Bonferroni correction
2*( 1 - pt( abs(del.resid.Loss[idx.tail[2:3]]), df = n.obs - p - 1 ) )
#==# Full: No outlier?! (also with two lambdas)
#==# W/o outlier: No outliers anymore.
idx.outl <- idx.tail[1]

# Since
# which( abs(r.Loss) > 2 ) # large residual
# which( diag(H.lambda.Loss) > 2*p/n.obs ) # high leverage
# deserve attention,
# which( Cook.Loss > 8/(n.obs - 2*p) ) # ~both
# is worth a closer look.

lrg.resid <- which( abs(r.Loss) > 2 )
length( lrg.resid )
#==# Full: 18 (... with two lambdas)
#==# W/o outlier: 19, sep.lam 18

# High leverage (Davison 2008, p.394 top)
idx.highlever <- which( diag(H.lambda.Loss) > 2*p/n.obs )
length( idx.highlever ) # quite a lot
#==# Full: 26 (24 with two lambdas)
#==# W/o outlier: 24
#==# W/o outl nor lev pts: 34
#==# W/o outl nor lev pts and lambda = 0.25: 33

Cook.Loss <- (r.Loss^2)*diag(H.lambda.Loss)/(2*(1 - diag(H.lambda.Loss)))
extr.Cdist <- which( Cook.Loss > 8/(n.obs-2*p) )
extr.Cdist
#==# Full: 51  69 145 224 359; with sep.lam c(51, 69, 145), 200 is close to limit bcs very high leverage
#==# W/o outlier: 248
#==# W/o outl. nor lev. pts: none
#==# W/o outl. nor lev. pts and lambda = 0.25: 326

any(lrg.resid %in% idx.highlever)
sum(lrg.resid %in% idx.highlever)
extr.Cdist %in% lrg.resid
extr.Cdist %in% idx.highlever
#==# Full: w/ sep.lam all three extr.Cdist have large residuals, 69 also high leverage
#==# W/o outlier: all extr.Cdist have large resdiuals, the last two also high leverage


# Tukey's one degree of freedom for non-additivity (Davison 2008, p.391)
X.delta.Loss <- cbind( X.lambda.Loss, yhat.lambda.Loss^2 )
lm.delta.Loss <- my.lm( y.lambda.Loss, X.delta.Loss )
SS.delta.Loss <- sum(lm.delta.Loss$e^2)
1 - pf( c( crossprod(e.Loss)- SS.delta.Loss )/lm.delta.Loss$s2, df1 = 1, df2 = lm.delta.Loss$n - lm.delta.Loss$p ) # equivalent to t-test for delta
#==# Full: 0.00003885 (arrgh!!) --> Maybe because of the outlier??
#==# W/o outlier: 2.6e-05 (arrgh!!) --> Apparently NOT because of the outlier...!!
#==# W/o outl. nor lev. pts: 0.01488 (arrgh!!) --> Apparently NOT because of outlier nor leverage point...!!
#==# W/o outl. nor lev. pts and lambda = 0.25: 0.0658 (petit ouf)


#======
#== Plots
#======

## Residuals
x11( width = 10.7, height = 10.5 )
par(mfrow = c(2,2))
dev.set(2)
# Normality of standardised residuals
qqnorm( r.Loss )
abline( 0, 1, col = 'red' )

# Std residuals vs fitted values: Fit doesn't seem too good!!
plot( yhat.lambda.Loss, r.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
abline( h = 0, col = 'red' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( yhat.lambda.Loss[extr.Cdist], r.Loss[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.Loss[idx.highlever], r.Loss[idx.highlever], col = 'blue' )
# quite obvious pattern

# Deletion or studentized residuals (plot not interesting)
#plot( del.resid.Loss, ylab = "Deletion residuals", main = "Absolute loss, transformed scale" )
#points( lrg.resid, del.resid.Loss[lrg.resid], col = 'sandybrown', pch = 19 )
# High leverage points not special in this plot

# Absolute std residuals vs. fitted values
plot( yhat.lambda.Loss, abs(r.Loss), xlab = "Fitted values", ylab = "", main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
points( yhat.lambda.Loss[extr.Cdist], abs(r.Loss)[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.Loss[idx.highlever], abs(r.Loss)[idx.highlever], col = 'blue' )
abline( h = 2, col = 'sandybrown', lty = 2 )
#
plot( yhat.lambda.Loss, sqrt(abs(r.Loss)), xlab = "Fitted values", ylab = "", main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
points( yhat.lambda.Loss[extr.Cdist], sqrt(abs(r.Loss[extr.Cdist])), col = 'red', pch = 19 )
points( yhat.lambda.Loss[idx.highlever], sqrt(abs(r.Loss[idx.highlever])), col = 'blue' ) 
abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
#text( -1.5, 2, labels = expression( sqrt(group("|","Standardised residuals","|")) ) )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )

# Running variance along fitted values for block sizes 20, 26, 32
x11( width = 9.8, height = 6.7 )
dev.set(5)
plot( sort(yhat.lambda.Loss)[9+seq(1,365,4)], sapply( seq(1,365,4), function(i){ var( r.Loss[order(yhat.lambda.Loss)][i+0:19], na.rm = TRUE ) } ), type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) ) # block size 20
lines( sort(yhat.lambda.Loss)[12+seq(1,357,4)], c( sapply( seq(1,353,4), function(i){ var( r.Loss[order(yhat.lambda.Loss)][i+0:25], na.rm = TRUE ) } ), var( r.Loss[order(yhat.lambda.Loss)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
lines( sort(yhat.lambda.Loss)[15+seq(1,353,4)], sapply( seq(1,353,4), function(i){ var( r.Loss[order(yhat.lambda.Loss)][i+0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
abline( h = 1, col = 'darkgrey', lty = 2 )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

# Compare residual plots for the two models (ONLY if separate lambda)
x11( width = 6.6, height = 7 )
dev.set(6)
#
plot( yhat.lambda.Loss, r.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, separate lambda" )
abline( h = 0, col = 'red' )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_compare", sep = '' ), type = "pdf" )


## Leverage etc.

# Cook's distance (plot not interesting)
#plot( Cook.Loss, xlab = "Index", ylab = "Cook's distance", main = "Absolute loss, transformed scale" )
#abline( h = 8/(n.obs - 2*p), col = 'red' )
# High leverage points not special in this plot

x11( width = 13.2, height = 7 )
par(mfrow = c(1,2))
dev.set(3)
# Distinguish outliers from leverage points (Davison 2008, p.395)
plot( diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)), Cook.Loss, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", main = "Absolute loss, transformed scale" )
points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[extr.Cdist], Cook.Loss[extr.Cdist], col = 'red', pch = 19 )
points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[idx.highlever], Cook.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
abline( h = 8/(n.obs - 2*p), col = 'red' )
points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[lrg.resid], Cook.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
abline( 0, 4/p, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
abline( v = 2*p/(n.obs - 2*p), col = 'blue', lty = 3 )

## Residuals vs. leverage (like plot.lm but with better level curves)
plot( diag(H.lambda.Loss), r.Loss, xlab = "Leverage", ylab = "Standardised residuals", main = "Absolute loss, transformed scale", xlim = c(0,max(diag(H.lambda.Loss))) ) # xlim = c(0,0.03)
#plot( diag(H.lambda.Loss), r.Loss, xlab = "Leverage", ylab = "Standardised residuals", main = "Absolute loss, transformed scale", xlim = c(0,0.03), ylim = c(-3,3) )
abline( h = 0, col = 'grey' )
abline( v = 2*p/n.obs, col = 'blue' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( diag(H.lambda.Loss)[extr.Cdist], r.Loss[extr.Cdist], col = 'red', pch = 19 )
points( diag(H.lambda.Loss)[idx.highlever], r.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
points( diag(H.lambda.Loss)[lrg.resid], r.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p)) * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'orange', lty = 3 )
# r^2 h / (p*(1-h)) = C
# r^2 = C * p*(1-h)/h
savePlot( paste( "FiguresLF/Loss_", modelname, "_leverage", sep = "" ), type = "pdf" )

head(rev(order(Cook.Loss)))


# Regression equation (response vs. covariate)
x11()
par(mfrow = c(1,1))
dev.set(4)
plot( x.lambda.Loss, y.lambda.Loss, xlab = "Structure", ylab = "Contents", main = "Absolute loss, transformed scale" )
abline( mle.Loss$beta.hat, col = 'red' )
points( x.lambda.Loss[extr.Cdist], y.lambda.Loss[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.Loss[idx.highlever], y.lambda.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
points( x.lambda.Loss[lrg.resid], y.lambda.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
lines( lowess( x = x.lambda.Loss, y = y.lambda.Loss ), col = 'green' )
#==# Full: --> 
#==# W/o outlier: --> Still the leverage point left.
#==# W/o outl nor lev pts: Nothing left.
savePlot( paste( "FiguresLF/Loss_", modelname, "_oview", sep = "" ), type = "pdf" )

# Output of the non-additivity model:
lm.delta.Loss$coefs
lm.delta.Loss$se.coefs
plot( c(X.delta.Loss %*% lm.delta.Loss$coefs), lm.delta.Loss$e/sqrt( lm.delta.Loss$s2 * (1-lm.delta.Loss$levs) ), xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, with additivity term" ) 
abline( h = 0, col = 'red' )
# Pattern of residuals basically the same as before. But this is not meaningful because the non-linear dependence doesn't need to have form of yhat^2.

## Try with x.lambda^2 as additional covariate
X.sq.lambda.Loss <- cbind( X.lambda.Loss, x.lambda.Loss^2 )
lm.sq.Loss <- my.lm( y.lambda.Loss, X.sq.lambda.Loss )
SS.sq.lambda.Loss <- sum(lm.sq.Loss$e^2)
lm.sq.Loss$coefs
lm.sq.Loss$se.coefs
plot( c(X.sq.lambda.Loss %*% lm.sq.Loss$coefs), lm.sq.Loss$e/sqrt( lm.sq.Loss$s2 * (1-lm.sq.Loss$levs) ), xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, with squared covariate" )
abline( h = 0, col = 'red' )
# Not really helpful, as before.

# Again the overall plot (response vs. covariate)
plot( x.lambda.Loss, y.lambda.Loss, xlab = "Structure", ylab = "Contents", main = "Absolute loss, transformed scale" )
abline( 0, 1, lty = 4, col = 'grey' )
abline( mle.Loss$beta.hat, col = 'red' )
lines( lowess( x = x.lambda.Loss, y = y.lambda.Loss ), col = 'blue', lwd = 1.2 )
#points( x.lambda.Loss[idx.outl], y.lambda.Loss[idx.outl], col = 'red', pch = 19 )
points( x.lambda.Loss[extr.Cdist], y.lambda.Loss[extr.Cdist], col = 'red', pch = 19 )
rug( x.lambda.Loss, quiet = TRUE, ticksize = 0.02 )
abline( v = c(15.5,16,16.5), col = 'blue', lty = 2)
#--> Try a model with two piecewise linear functions (1 knot)

#=====
# Fit piecewise linear function  !!!--only for full data--!!! (not good idea)
#=====
modelname <- paste( modelname, "_split", sep = '' )
p.split <- p+1

# Optim is not very stable, in particular in the subsequent steps:
opt.knot.Loss <- optim( c( 5, 0.5, 0.5, 4, 0.15, 15 ), fn = PTBS.llkhd.knot, y = y.Loss, x = x.Loss, control = list( fnscale = -1, maxit = 10000 ), hessian = TRUE ) # full data
opt.knot.Loss <- 
optim( c( 8, 0.3, 0.5, 3, 0.1, 17 ), fn = PTBS.llkhd.knot, y = y.Loss, x = x.Loss, control = list( fnscale = -1, maxit = 10000 ), hessian = F ) # w/o outlier c( 8, 0.3, 0.5, 3, 0.1, 17 )
opt.knot.Loss$par
sqrt( diag(solve( -opt.knot.Loss$hessian )) )

# Use profiling instead
### Estimate lambda again
modelname <- paste( modelname, "_lamgen", sep = '' )
prof.knot.Loss <- sapply( seq(12,21,0.01), function(kn){ optimize( f = PTBS.profllkhd.knot.lambda, interval = c(-3,1), y = y.Loss, x = x.Loss, knot.x = kn, intercept = TRUE, maximum = TRUE )$objective } )
plot( seq(12,21,0.01), prof.knot.Loss )

knot.Loss <- optimize( f = PTBS.profllkhd.knot.kn, interval = c(13.5,14.5), y = y.Loss, x = x.Loss, intercept = TRUE, intval.lam = c(-3,1), maximum = TRUE )$maximum
abline( v = knot.Loss, col = 'red' )

mle.split.Loss <- PTBS.knot.mle( y = y.Loss, x = x.Loss, knot = knot.Loss, interval = c(-3,1) )

PTBS.knot.CI.lambda( y = y.Loss, x = x.Loss, knot = knot.Loss, interval = c(-3,1) )
# CI = (0.06709,0.1067) --> no overlap
plot( lam.seq.Loss, sapply( lam.seq.Loss, PTBS.profllkhd.knot.lambda, y = y.Loss, x = x.Loss, knot = knot.Loss, intercept = TRUE ), type = 'l', xlim = c(-0.5,0.5) )
abline( v = mle.split.Loss$lambda.hat, col = 'red' )
abline( v = unlist( PTBS.knot.CI.lambda( y = y.Loss, x = x.Loss, knot = knot.Loss, interval = c(-3,1) )[2:3] ), col = 'red', lty = 2 )
# quite flat!!

CI.lambda.Loss
# overlap here!
## Nevertheless keep same lambda for both cases??

covML.split.Loss <- solve( -optimHess( unlist(mle.split.Loss[c(2,3,1)]), PTBS.llkhd.knot, theta = c(rep(NA,5), knot.Loss), y = y.Loss, x = x.Loss, control = list( maxit = 5000, fnscale = -1 ) ) ) # knot considered as fixed

### Keep lambda fixed (not sure this is so useful) !!! reset lrg.resid, idx.highlever and extr.Cdist to non-split model !!!
lrg.resid <- which( abs(r.Loss) > 2 )
idx.highlever <- which( diag(H.lambda.Loss) > 2*p/n.obs )
extr.Cdist <- which( Cook.Loss > 8/(n.obs-2*p) )
modelname <- paste( unlist( strsplit( modelname, "_lam" ) )[1], "_lamfix", sep = '' )

prof.knot.Loss <- sapply( seq(12,21,0.01), function(kn){ PTBS.profllkhd.knot.lambda( lambda = mle.Loss$lambda.hat, y = y.Loss, x = x.Loss, knot.x = kn, intercept = TRUE ) } )
plot( seq(12,21,0.01), prof.knot.Loss )

knot.Loss <- optimize( f = PTBS.profllkhd.knot.kn, interval = c(17,19), y = y.Loss, x = x.Loss, lambda = mle.Loss$lambda.hat, intercept = TRUE, intval.lam = c(-3,1), maximum = TRUE )$maximum
mle.split.Loss <- PTBS.knot.mle( y = y.Loss, x = x.Loss, lambda = mle.Loss$lambda.hat, knot = knot.Loss )
abline( v = knot.Loss, col = 'red' )

covML.split.Loss <- solve( -optimHess( unlist(mle.split.Loss[c(2,3)]), PTBS.llkhd.knot, theta = c(rep(NA,4), mle.split.Loss$lambda.hat, knot.Loss), y = y.Loss, x = x.Loss, control = list( maxit = 5000, fnscale = -1 ) ) ) # knot considered as fixed, lambda fixed as well
###

x.lambda.split.Loss <- BC.transform( mle.split.Loss$lambda.hat, x.Loss )

y.lambda.split.Loss <- BC.transform( mle.split.Loss$lambda.hat, y.Loss )
X.lambda.split.Loss <- cbind( 1, pmin( x.lambda.split.Loss, knot.Loss ), pmax( x.lambda.split.Loss - knot.Loss, 0 ) )

H.lambda.split.Loss <- X.lambda.split.Loss %*% solve( crossprod(X.lambda.split.Loss) ) %*% t(X.lambda.split.Loss)

cbind( unlist( mle.split.Loss[c(2,3,1)] ), sqrt( diag(covML.split.Loss) ) )
mle.split.Loss$ell.opt

# Model comparison
2*( -mle.split.Loss$ell.opt + p.split+2+grepl("lamgen",modelname) )
2*( -mle.Loss$ell.opt + p+2 ) # as comparison
# AIC suggests that split-model is better
1 - pchisq( 2*(mle.split.Loss$ell.opt - mle.Loss$ell.opt), df = 1 )  # LRT (also valid if lambda estimated)

# Fitted values:
yhat.lambda.split.Loss <- c( H.lambda.split.Loss %*% y.lambda.split.Loss )
# Residuals:
e.split.Loss <- y.lambda.split.Loss - yhat.lambda.split.Loss
# Standardised residuals:
r.split.Loss <- e.split.Loss/sqrt( (n.obs/(n.obs-p.split))*mle.split.Loss$sigma2.hat * (1-diag(H.lambda.split.Loss)) )
# Deletion residuals (Davison 2008, p.395) = "studentized residuals" in Weisberg 2005 (p.196):
del.resid.split.Loss <- sqrt( (n.obs-p.split-1)/(n.obs-p.split-r.split.Loss^2) )*r.split.Loss

# Fit model with two different sigma2 but same knot:
pars.curr <- unlist(mle.split.Loss[c(2,3,3,if(grepl("lamgen",modelname)){ 1 } else { NULL })])
##ADAPT values!!!!
intvl <- list( c(-5,0), c(0,1.5), c(0,2.5), c(0,2), c(0,2), if(grepl("lamgen",modelname)){ c(-3,1) } else { NULL } ) # last element is present even if NULL, but this doesn't matter
R <- length(pars.curr)
ell.vec <- rep(NA,R)
repeat{
	for( r in 1:R ) {
		pars <- pars.curr
		pars[r] <- NA
		opt.r <- optimize( f = PTBS.llkhd.knot, interval = intvl[[r]], theta = c( pars, if(grepl("lamfix",modelname)){ mle.split.Loss$lambda.hat } else { NULL }, knot.Loss ), y = y.Loss, x = x.Loss, var.split = TRUE, maximum = TRUE )
		pars.curr[r] <- opt.r$maximum
		ell.vec[r] <- opt.r$objective
	}
	if( max( abs( combn( R, 2, function(i){ diff(ell.vec[i]) } ) ) ) < 1e-07 ) {
		break
	}
}
ell.vec
pars.curr
1 - pchisq( 2*(ell.vec[R] - mle.split.Loss$ell.opt), df = 1 )
# Two different sigma2 doesn't look like an improvement.
# W/o outlier for lambda fixed, the p-value is 

# Other diagnostics for split model

## Outliers?
idx.tail <- rev( tail( order( abs(del.resid.split.Loss) ) ) )
2*( 1 - pt( abs(del.resid.split.Loss[idx.tail[1]]), df = n.obs - p.split - 1 ) )
0.05/n.obs  # Bonferroni correction
2*( 1 - pt( abs(del.resid.split.Loss[idx.tail[2:3]]), df = n.obs - p.split - 1 ) )
idx.tail[1]
#==# Full: Even the most extreme one (145) is (just!) no outlier.
#==# W/o outlier: No outliers anymore.
#==# W/o outl nor lev pts: No outliers anymore.

lrg.resid <- which( abs(r.split.Loss) > 2 )
length( lrg.resid )
#==# Full: 19 if lambda general and ... fixed
#==# W/o outlier: ... if lambda general, ... if lambda fixed

# High leverage (Davison 2008, p.394 top)
idx.highlever <- which( diag(H.lambda.split.Loss) > 2*p.split/n.obs )
length( idx.highlever ) # quite a lot
#==# Full: 37
#==# W/o outlier:
#==# W/o outl nor lev pts: 33
#==# W/o outl nor lev pts and lambda = 0.25: 33

Cook.split.Loss <- (r.split.Loss^2)*diag(H.lambda.split.Loss)/(p.split*(1 - diag(H.lambda.split.Loss)))
extr.Cdist <- which( Cook.split.Loss > 8/(n.obs-2*p.split) )
extr.Cdist
#==# Full: 69, 359; 145 is just below the limit
#==# W/o outlier: 248
#==# W/o outl. nor lev. pts: none
#==# W/o outl. nor lev. pts and lambda = 0.25: 326

any(lrg.resid %in% idx.highlever)
extr.Cdist %in% lrg.resid
extr.Cdist %in% idx.highlever
#==# Full: 
#==# W/o outlier:

# Tukey's one degree of freedom for non-additivity (Davison 2008, p.391)
X.delta.split.Loss <- cbind( X.lambda.split.Loss, yhat.lambda.split.Loss^2 )
H.delta.split.Loss <- X.delta.split.Loss %*% solve( crossprod(X.delta.split.Loss) ) %*% t(X.delta.split.Loss)
SS.delta.split.Loss <- c( y.lambda.split.Loss %*% ( diag(n.obs) - H.delta.split.Loss ) %*% y.lambda.split.Loss )
1 - pf( c( crossprod(e.split.Loss) - SS.delta.split.Loss )/( SS.delta.split.Loss/(n.obs - p.split - 1) ), df1 = 1, df2 = n.obs-p.split-1 )  # equivalent to t-test for delta
#==# Full: 0.7695 (ouf!)


#=====
#== Plots
#=====
plot( x.lambda.split.Loss, y.lambda.split.Loss, xlab = "Structure", ylab = "Contents", main = paste( "Absolute loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( 0, 1, lty = 4, col = 'grey' )
lines( lowess( x = x.lambda.split.Loss, y = y.lambda.split.Loss ), col = 'green' )
abline( c( solve( crossprod( cbind( 1, x.lambda.split.Loss ) ) ) %*% t( cbind( 1, x.lambda.split.Loss ) ) %*% y.lambda.split.Loss ), col = 'red' ) # linear fit for current lambda
#points( x.lambda.split.Loss[idx.outl], y.lambda.split.Loss[idx.outl], col = 'red', pch = 19 )
points( x.lambda.split.Loss[extr.Cdist], y.lambda.split.Loss[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.split.Loss[idx.highlever], y.lambda.split.Loss[idx.highlever], col = 'blue' )
points( x.lambda.split.Loss[lrg.resid], y.lambda.split.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
rug( x.lambda.split.Loss, quiet = TRUE, ticksize = 0.02 )
lines( xlam.seq.Loss, c( cbind( 1, pmin( xlam.seq.Loss, knot.Loss ), pmax( xlam.seq.Loss - knot.Loss, 0 ) ) %*% mle.split.Loss$beta.hat ), col = 'green4' )
abline( v = knot.Loss, col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/Loss_", modelname, "_oview1", sep = '' ), type = "pdf" )

# Compare residual plots for the two models
x11( width = 13.2, height = 7 )
par( mfrow = c(1,2) )
dev.set(3)
#
plot( yhat.lambda.Loss, r.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, initial model" )
abline( h = 0, col = 'red' )
abline( v = c( mle.Loss$beta.hat %*% c( 1, BC.transform( mle.Loss$lambda.hat, BC.backtransform( mle.split.Loss$lambda.hat, knot.Loss ) ) ) ), col = 'green4', lty = 3 )
#
plot( yhat.lambda.split.Loss, r.split.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = paste( "Absolute loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) ) 
abline( h = 0, col = 'red' )
abline( v = c( mle.split.Loss$beta.hat[1:2] %*% c(1,knot.Loss) ), col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_compare", sep = '' ), type = "pdf" )


## Residuals
par(mfrow = c(2,2))
dev.set(2)
# Normality of standardised residuals
qqnorm( r.split.Loss )
abline( 0, 1, col = 'red' )

# Std residuals vs fitted values:
plot( yhat.lambda.split.Loss, r.split.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, transformed scale" )
abline( h = 0, col = 'red' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( yhat.lambda.split.Loss[extr.Cdist], r.split.Loss[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.split.Loss[idx.highlever], r.split.Loss[idx.highlever], col = 'blue' )
# quite obvious pattern
#==# Full: 
#==# W/o outl nor lev pts:
# High leverage not very interesting since just largest and smallest values of x...

# Deletion or studentized residuals (plot not interesting)
#plot( del.resid.split.Loss, ylab = "Deletion residuals", main = "Absolute loss, transformed scale" )
# High leverage points not special in this plot

# Absolute std residuals vs. fitted values (see heteroscedasticity)
plot( yhat.lambda.split.Loss, abs(r.split.Loss), xlab = "Fitted values", ylab = "", main = "Absolute loss, transformed scale" )
title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
points( yhat.lambda.split.Loss[extr.Cdist], abs(r.split.Loss)[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.split.Loss[idx.highlever], abs(r.split.Loss)[idx.highlever], col = 'blue' )
abline( h = 2, col = 'sandybrown', lty = 2 )
#
plot( yhat.lambda.split.Loss, sqrt(abs(r.split.Loss)), xlab = "Fitted values", ylab = "", main = "Absolute loss, transformed scale" )
title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
points( yhat.lambda.split.Loss[extr.Cdist], sqrt(abs(r.split.Loss[extr.Cdist])), col = 'red', pch = 19 )
points( yhat.lambda.split.Loss[idx.highlever], sqrt(abs(r.split.Loss[idx.highlever])), col = 'blue' )
abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
#text( -1.5, 2, labels = expression( sqrt(group("|","Standardised residuals","|")) ) )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )

#==# Full: fit and homoskedasticity reasonable (lambda general)
#==# W/o outlier: fit & homosk. not too bad, smaller variance for high-leverage points with large x?

# Running variance along fitted values for block sizes 20, 26, 32
x11( width = 9.8, height = 6.7 )
dev.set(5)
plot( sort(yhat.lambda.split.Loss)[9+seq(1,365,4)], sapply( seq(1,365,4), function(i){ var( r.split.Loss[order(yhat.lambda.split.Loss)][i+0:19], na.rm = TRUE ) } ), type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = "Absolute loss, transformed scale" ) # block size 20
lines( sort(yhat.lambda.split.Loss)[12+seq(1,357,4)], c( sapply( seq(1,353,4), function(i){ var( r.split.Loss[order(yhat.lambda.split.Loss)][i+0:25], na.rm = TRUE ) } ), var( r.split.Loss[order(yhat.lambda.split.Loss)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
lines( sort(yhat.lambda.split.Loss)[15+seq(1,353,4)], sapply( seq(1,353,4), function(i){ var( r.split.Loss[order(yhat.lambda.split.Loss)][i+0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
abline( h = 1, col = 'darkgrey', lty = 2 )
savePlot( paste( "FiguresLF/Loss_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

## Leverage etc.
par(mrow = c(1,2))
dev.set(3)
# Cook's distance (plot not interesting)
#plot( Cook.split.Loss, xlab = "Index", ylab = "Cook's distance" )
#abline( h = 8/(n.obs - 2*p.split), col = 'red' )
# High leverage points not special in this plot

# Distinguish outliers from leverage points (Davison 2008, p.395)
plot( diag(H.lambda.split.Loss)/(1-diag(H.lambda.split.Loss)), Cook.split.Loss, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", main = "Absolute loss, transformed scale" )
points( (diag(H.lambda.split.Loss)/(1-diag(H.lambda.split.Loss)))[extr.Cdist], Cook.split.Loss[extr.Cdist], col = 'red', pch = 19 )
points( (diag(H.lambda.split.Loss)/(1-diag(H.lambda.split.Loss)))[idx.highlever], Cook.split.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
abline( h = 8/(n.obs - 2*p.split), col = 'red' )
points( (diag(H.lambda.split.Loss)/(1-diag(H.lambda.split.Loss)))[lrg.resid], Cook.split.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
abline( 0, 4/p.split, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
abline( v = 2*p.split/(n.obs - 2*p.split), col = 'blue', lty = 3 )

# Residuals vs. leverage (like plot.lm but with better level curves)
plot( diag(H.lambda.split.Loss), r.split.Loss, xlab = "Leverage", ylab = "Standardised residuals", main = "Absolute loss, transformed scale", xlim = c(0,max(diag(H.lambda.split.Loss))) ) # xlim = c(0,0.025)
#plot( diag(H.lambda.split.Loss), r.split.Loss, xlab = "Leverage", ylab = "Standardised residuals", main = "Absolute loss, transformed scale", xlim = c(0,0.025), ylim = c(-3,3) )
abline( h = 0, col = 'grey' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
abline( v = 2*p.split/n.obs, col = 'blue' )
points( diag(H.lambda.split.Loss)[extr.Cdist], r.split.Loss[extr.Cdist], col = 'red', pch = 19 )
points( diag(H.lambda.split.Loss)[idx.highlever], r.split.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
points( diag(H.lambda.split.Loss)[lrg.resid], r.split.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p.split)) * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'orange', lty = 3 )
# r^2 h / (p.split*(1-h)) = C
# r^2 = C * p.split*(1-h)/h
savePlot( paste( "FiguresLF/Loss_", modelname, "_leverage", sep = "" ), type = "pdf" )

head(rev(order(Cook.split.Loss)))

dev.set(4)
plot( x.lambda.split.Loss, y.lambda.split.Loss, xlab = "Structure", ylab = "Contents", main = paste( "Absolute loss, split model",  tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( h = 0, lty = 2, col = 'grey' )
abline( v = 0, lty = 2, col = 'grey' )
abline( 0, 1, lty = 4, col = 'grey' )
lines( lowess( x = x.lambda.split.Loss, y = y.lambda.split.Loss ), col = 'green' )
abline( c( solve( crossprod( cbind( 1, x.lambda.split.Loss ) ) ) %*% t( cbind( 1, x.lambda.split.Loss ) ) %*% y.lambda.split.Loss ), col = 'red' ) # linear fit for current lambda
points( x.lambda.split.Loss[extr.Cdist], y.lambda.split.Loss[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.split.Loss[idx.highlever], y.lambda.split.Loss[idx.highlever], col = 'blue' )
points( x.lambda.split.Loss[lrg.resid], y.lambda.split.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
rug( x.lambda.split.Loss, quiet = TRUE, ticksize = 0.02 )
lines( xlam.seq.Loss, c( cbind( 1, pmin( xlam.seq.Loss, knot.Loss ), pmax( xlam.seq.Loss - knot.Loss, 0 ) ) %*% mle.split.Loss$beta.hat ), col = 'green4' )
abline( v = knot.Loss, col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/Loss_", modelname, "_oview2", sep = '' ), type = "pdf" )

pars.curr <- c( 8, 0.2, 0.6, 4, 0.2, 15 )
R <- length(pars.init)
intvl <- list( c(4,10), c(0,1), c(0,2), c(1,7), c(-3,1), c(10,20) ) 
ell.vec <- rep(NA,R)

repeat{
	for( r in 1:R ) {
		pars <- pars.curr
		pars[r] <- NA
		opt.r <- optimize( f = PTBS.llkhd.knot, interval = intvl[[r]], theta = pars, y = y.Loss, x = x.Loss, maximum = TRUE )
		pars.curr[r] <- opt.r$maximum
		ell.vec[r] <- opt.r$objective
	}
	if( max( abs( combn( R, 2, function(i){ diff(ell.vec[i]) } ) ) ) < 1e-07 ) {
		break
	}
}
ell.vec
pars.curr
# This is also sensitive to the starting values...!

#====== End piecewise linear function =======#



######
# CI and PI on transformed scale
######

# There are three possible scenarios under which standard errors of prediction x^T*beta.hat can be computed:
# 1. Delta method on full covariance matrix: Account for estimation of lambda. This is not proper because uncertainty of lambda in y.lambda is not accounted for.
# 3. Regression-based while treating lambda = lambda.hat as fixed. Thus the uncertainty of estimating lambda is contained nowhere, not even in the std errors of beta.
# 2. Use std errors for beta from full covariance matrix (thus estimation of lambda is accounted for in se.beta), but then ignore that lambda was estimated when appyling the delta method for var(x0^T beta.hat).
# 4. As 2 but bootstrap-based.
# 5. Detour via original scale by computing confidence bands for y.hat (p-quantile) (accounting for uncertainty in lambda), which are then transformed again.

# --> I have two (very similar) sets of estimated parameters (mle and opt). Based on my calculations, I could also code the complete observed information matrix for mle, but currently I haven't done it. Therefore I have to used the opt-estimates for full uncertainty.
# --> By applying optimHess to mle.Loss I get the Hessian and thus the joint covariance at mle.DoL. Thus I do not need the opt version anymore. Even though 2 seems better than 3, I keep 2, 3, and 5.

fitted.mle.seq.Loss <- c( cbind( 1, xlam.seq.Loss ) %*% mle.Loss$beta.hat )

# (2) Assuming lambda = lambda.hat known, based on mle
var.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ c(1,x) %*% covML.Loss[1:p,1:p] %*% c(1,x) } )
var.predict.seq.Loss <- var.fitted.seq.Loss + s2.hat.Loss

# (3) Cov matrix of beta.hat from the regression of y.lambda on X.lambda with lambda = lambda.hat fixed
varRegr.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ c(1,x) %*% covbeta.regr.Loss %*% c(1,x) } )

# (5) Detour via original scale to account for estimated lambda in y as well (but not sure it is properly done)
varOrig.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(xlam){ grad.yhat( theta = par.Loss, x.lam = c(1,xlam), sep.lam = sep.lam.global ) %*% covML.Loss %*% grad.yhat( theta = par.Loss, x.lam = c(1,xlam), sep.lam = sep.lam.global ) } )

apply( t( sapply( xlam.seq.Loss, function(xlam){ grad.yhat( theta = par.Loss, x.lam = c(1,xlam), sep.lam = sep.lam.global ) } ) ), 2, range )

# (1) Including uncertainty on lambda.hat, based on opt
varFull.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ grad.fitted( par.Loss, c(1,x) ) %*% covML.Loss %*% grad.fitted( par.Loss, c(1,x) ) } ) # sigma2 has no influence on this
varFull.predict.seq.Loss <- varFull.fitted.seq.Loss + s2.hat.Loss

# (4) As 1 but with cov(beta) based on non-parametric bootstrap 
# varBoot.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ grad.fitted( par.Loss, c(1,x) ) %*% covBoot.Loss %*% grad.fitted( par.Loss, c(1,x) ) } )
# varBoot.predict.seq.Loss <- varBoot.fitted.seq.Loss + s2.hat.Loss

x11()
plot( xlam.seq.Loss, varOrig.fitted.seq.Loss, type = 'l', xlab = "x(lambda)", ylab = "Variance of X(lambda)*beta.hat", main = "Absolute loss, transformed scale" )
#lines( xlam.seq.Loss, varBoot.fitted.seq.Loss, col = 'dodgerblue' )
lines( xlam.seq.Loss, var.fitted.seq.Loss, col = 'blue' )
lines( xlam.seq.Loss, varRegr.fitted.seq.Loss, col = 'cyan1' )
dev.off()

x11()
# Transformed scale: y.lambda vs. X.lambda*beta
plot( x.lambda.Loss, y.lambda.Loss, xlab = "Structure", ylab = "Contents", main = paste( "Absolute loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
# if appropriate:
points( BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], Structure$Loss[idx.outl] ), BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[idx.outl] ), pch = 4 ) # ev. pch = 8
points( BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], Structure$Loss[idx.remove] ), BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[idx.remove] ), pch = 4 )
#
abline( mle.Loss$beta.hat, col = 'red' )
#abline( opt.Loss$par[1:p], col = 'lightpink2' )

# 95% Confidence interval (mle, lambda fixed --> 3)
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ), col = 'orange', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ), col = 'orange', lty = 2 )
## much smaller than ML ones for lower X.lambda

# 95% Prediction interval (mle, lambda fixed --> 3)
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ), col = 'cyan1', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ), col = 'cyan1', lty = 2 )

# 95% Confidence interval (mle, lambda known --> 2)  ---> use qnorm??
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), col = 'red', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), col = 'red', lty = 2 )

# 95% Prediction interval (mle, lambda known --> 2)
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ), col = 'blue', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ), col = 'blue', lty = 2 )
# --> Practically no difference of prediction intervals because sigma2 dominates the variance.

# 95% Confidence interval (opt, via original scale --> 5)
lines( xlam.seq.Loss, BC.transform( mle.Loss$lambda.hat[1], BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss) ), col = 'mediumorchid1', lty = 2 ) # sienna1
lines( xlam.seq.Loss, BC.transform( mle.Loss$lambda.hat[1], BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) + qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss) ), col = 'mediumorchid1', lty = 2 )
# The detour via the original scale is closest to 3, whereas 2 is too wide for smaller x-values (!!??)
# Prediction interval accounting for estimation uncertainty is too complex for this case
lines( lowess( x = x.lambda.Loss, y = y.lambda.Loss ), col = 'green' )

# 95% Confidence interval (opt, lambda unknown --> 1)
lines( xlam.seq.Loss, fitted.opt.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ), col = 'lightpink2', lty = 2 )
lines( xlam.seq.Loss, fitted.opt.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ), col = 'lightpink2', lty = 2 )

# 95% Prediction interval (opt, lambda unknown --> 1)
lines( xlam.seq.Loss, fitted.opt.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ), col = 'deepskyblue', lty = 2 )
lines( xlam.seq.Loss, fitted.opt.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ), col = 'deepskyblue', lty = 2 )

# 95% Confidence interval (mle, lambda unknown --> 4)
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.Loss ), col = 'chocolate3', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.Loss ), col = 'chocolate3', lty = 2 )  # use qnorm here???

# 95% Prediction interval (mle, lambda unknown --> 4)
lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.Loss ), col = 'dodgerblue', lty = 2 )
lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.Loss ), col = 'dodgerblue', lty = 2 )

## Either 2 or 3 seem reasonable since very close. 1 and 4 are not plausible because uncertainty of lambda is accounted for in x, but not in y. And the detour via the original scale is very close to 2 or 3.

( BC.transform( opt.Loss$par[4], BC.backtransform( opt.Loss$par[4], fitted.opt.seq.Loss ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss) ) - ( fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ) )[-(1:12)]
########


########
# Backtransformed to original scale
########

x.seq.Loss <- BC.backtransform( mle.Loss$lambda.hat[1+sep.lam.global], xlam.seq.Loss )

## Conditional mean on original scale, estimator by Taylor 1986
condmean.seq.Loss <- condmean.orig( theta = par.Loss, x = cbind(1,xlam.seq.Loss), sep.lam = sep.lam.global )
## use sigma2 or s2 here???
var.condmean.Loss <- sapply( xlam.seq.Loss, function(x){ c( grad.condmean.orig( theta = par.Loss, x = c(1,x), sep.lam = sep.lam.global ) %*% covML.Loss %*% grad.condmean.orig( theta = par.Loss, x = c(1,x), sep.lam = sep.lam.global ) ) } )


plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Contents", main = paste( "Absolute loss", ifelse( sep.lam.global, " seplam", "" ), ", original scale", sep = '' ) )
# if appropriate:
points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4 ) # ev. pch = 8
points( Structure$Loss[idx.remove], Contents$Loss[idx.remove], pch = 4 )
#
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ), col = 'red' ) # fitted median 
#lines( x.seq.Loss, BC.backtransform( opt.Loss$par[p+2], fitted.opt.seq.Loss ), col = 'lightpink2' )

# Backtransformed 95% Confidence interval (mle, lambda fixed --> 3)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ) ), col = 'orange', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ) ), col = 'orange', lty = 2 )

# Backtransformed 95% Prediction interval (mle, lambda fixed --> 3)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ) ), col = 'cyan1', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ) ), col = 'cyan1', lty = 2 )

# Backtransformed 95% Confidence interval (mle, lambda known --> 2)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), col = 'red', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), col = 'red', lty = 2 )

# Backtransformed 95% Prediction interval (mle, lambda known --> 2)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), col = 'blue', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), col = 'blue', lty = 2 )

# Backtransformed 95% Confidence interval (opt, var of back-transformed fitted value --> 5) 
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss), col = 'sienna1', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) + qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss), col = 'sienna1', lty = 2 )

# Backtransformed 95% Confidence interval (opt, lambda unknown --> 1)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ) ), col = 'lightpink2', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ) ), col = 'lightpink2', lty = 2 )  # use qnorm here???

# Backtransformed 95% Prediction interval (opt, lambda unknown --> 1)
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ) ), col = 'deepskyblue', lty = 2 )
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ) ), col = 'deepskyblue', lty = 2 )

# Conditional mean on original scale
lines( x.seq.Loss, condmean.seq.Loss, col = 'green3', lwd = 1 )
lines( x.seq.Loss, condmean.seq.Loss - qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )
lines( x.seq.Loss, condmean.seq.Loss + qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )

# Smearing estimate of the mean on original scale (Duan, 1983)
lines( x.seq.Loss, sapply( xlam.seq.Loss, function(x){ mean( BC.backtransform( mle.Loss$lambda.hat, sum( c(1,x)*mle.Loss$beta.hat ) + e.Loss ) ) } ), col = 'mediumpurple1' )
# Basically identical with condmean.seq (relative to mean of Contents$Loss)
max( abs( condmean.seq.Loss - sapply( xlam.seq.Loss, function(x){ mean( BC.backtransform( mle.Loss$lambda.hat[1], sum( c(1,x)*mle.Loss$beta.hat ) + e.Loss ) ) } ) ) )/mean(range(Contents$Loss))



## "Confidence interval" for the mean on original scale: Use 95% CI bounds as x instead of the fitted values (what Markus had done...)
lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), col = 'magenta', lty = 2 )
lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), col = 'magenta', lty = 2 )

# "Prediction interval" for the mean on original scale (obtained as before)
lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), col = 'skyblue', lty = 2 )
lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), col = 'skyblue', lty = 2 )


#####
#°°° TBS model (named "other") °°°#
#####
# loglik.other( c(3,0.5,5,0.25), y = y.Loss, x = x.Loss )
# ( opt.other.Loss <- optim( c(1200,0.4,5,0.1), loglik.other, y = y.Loss, x = x.Loss, control = list( fnscale = -1, maxit = 10000 ), hessian = TRUE ) )
# sqrt( diag( solve( -opt.other.Loss$hessian ) ) )

plot( lam.seq.Loss, sapply( lam.seq.Loss, TBS.profllkhd.lambda, y = y.Loss, x = x.Loss, intercept = TRUE, beta.init = c(1000,1)) )
mle.other.Loss <- TBS.mle( y = y.Loss, x = x.Loss, beta.init = c(1000,0.5), interval = c(-3,1) )
covML.other.Loss <- solve( -optimHess( unlist( mle.other.Loss[c(2,3,1)] ), loglik.other, y = y.Loss, x = x.Loss, control = list( fnscale = -1 ) ) )
optim( unlist( mle.other.Loss[c(2,3,1)] ), loglik.other, y = y.Loss, x = x.Loss, control = list( fnscale = -1 ) )
sqrt( diag( covML.other.Loss ) )
TBS.CI.lambda( y.Loss, x.Loss, interval = c(-3,1), beta.init = c(1000,0.5) )

yhat.Loss <- c( cbind( 1, x.Loss ) %*% mle.other.Loss$beta.hat )
y.lambda.Loss <- BC.transform( mle.other.Loss$lambda.hat, y.Loss )
yhat.lambda.Loss <- BC.transform( mle.other.Loss$lambda.hat, yhat.Loss )
e.Loss <- y.lambda.Loss - yhat.lambda.Loss
x.lambda.Loss <- BC.transform( mle.other.Loss$lambda.hat, x.Loss )

plot( qnorm( 1:n.obs/(n.obs+1) ), sort( e.Loss/sqrt(mle.other.Loss$sigma2.hat) ) )
abline( 0, 1, col = 'red' )

plot( yhat.lambda.Loss, e.Loss, xlab = "Fitted values", ylab = "Residuals", main = "Absolute loss" ) # /median( abs(e.Loss) )
abline( h = 0, col = 'red' )
lines( lowess( x = yhat.lambda.Loss, y = e.Loss ), col = 'green' )

fitted.seq.Loss <- c( cbind( 1, x.seq.Loss ) %*% mle.other.Loss$beta.hat )

plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Content", main = "Absolute loss, original scale" )
points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4 )
abline( mle.other.Loss$beta.hat, col = 'red' )
lines( lowess( x = x.Loss, y = y.Loss ), col = 'green' )
# 95% CI of x^T beta.hat (which is the fitted median):
lines( x.seq.Loss, fitted.seq.Loss - qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ c( c(1,x) %*% covML.other.Loss[1:p,1:p] %*% c(1,x) ) } ) ), col = 'red', lty = 2 )
lines( x.seq.Loss, fitted.seq.Loss + qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ c( c(1,x) %*% covML.other.Loss[1:p,1:p] %*% c(1,x) ) } ) ), col = 'red', lty = 2 )
# Estimated quantiles (2.5% and 97.5%) --> these are an approximate 95% PI:
lines( x.seq.Loss, ( fitted.seq.Loss^mle.other.Loss$lambda.hat + qnorm(0.025) * mle.other.Loss$lambda.hat * sqrt(mle.other.Loss$sigma2.hat) )^(1/mle.other.Loss$lambda.hat), col = 'blue', lty = 2 )
lines( x.seq.Loss, ( fitted.seq.Loss^mle.other.Loss$lambda.hat + qnorm(0.975) * mle.other.Loss$lambda.hat * sqrt(mle.other.Loss$sigma2.hat) )^(1/mle.other.Loss$lambda.hat), col = 'blue', lty = 2 )

# Conditional mean:
condmean.other.seq.Loss <- TBS.condmean.orig( theta = unlist(mle.other.Loss[c(2,3,1)]), x.beta = fitted.seq.Loss )
var.condmean.other.Loss <- sapply( x.seq.Loss, function(x){ xbeta <- c( c(1,x) %*% mle.other.Loss$beta.hat ); grad.mean <- c( ( 1 + mle.other.Loss$sigma2.hat * (1 - mle.other.Loss$lambda.hat) * (1 - 2*mle.other.Loss$lambda.hat)/( 2*xbeta^(2*mle.other.Loss$lambda.hat) ) ) * c(1,x),  (1 - mle.other.Loss$lambda.hat)/( 2*xbeta^(2*mle.other.Loss$lambda.hat - 1) ), -mle.other.Loss$sigma2.hat/( 2*xbeta^(2*mle.other.Loss$lambda.hat - 1) ) * ( 1 + 2*(1 - mle.other.Loss$lambda.hat) * log(xbeta) ) ); c( grad.mean %*% covML.other.Loss %*% grad.mean ) } )

lines( x.seq.Loss, condmean.other.seq.Loss, col = 'green3' )
# 95% CI:
lines( x.seq.Loss, condmean.other.seq.Loss - qnorm(0.975) * sqrt( var.condmean.other.Loss ), col = 'green3', lty = 2 )
lines( x.seq.Loss, condmean.other.seq.Loss + qnorm(0.975) * sqrt( var.condmean.other.Loss ), col = 'green3', lty = 2 )

# Quantity that should be small for Taylor appx to be valid
range( mle.other.Loss$lambda.hat*sqrt(mle.other.Loss$sigma2.hat)/c(cbind(1,x.Loss)%*% mle.other.Loss$beta.hat)^mle.other.Loss$lambda.hat )

# Plot of TBS model on transformed scale:
plot( BC.transform( mle.other.Loss$lambda.hat, x.Loss ), y.lambda.Loss, xlab = "Structure", ylab = "Content", main = "Absolute loss, transformed scale" )
lines( BC.transform( mle.other.Loss$lambda.hat, sort(x.Loss) ), BC.transform( mle.other.Loss$lambda.hat, yhat.Loss[order(x.Loss)] ), col = 'red' )
lines( lowess( x = BC.transform( mle.other.Loss$lambda.hat, x.Loss ), y = y.lambda.Loss ), col = 'green' )
abline( 0, 1, col = 'orange' )
lines( BC.transform( mle.other.Loss$lambda.hat, x.seq.Loss ), BC.transform( mle.other.Loss$lambda.hat, c( cbind( 1, x.seq.Loss ) %*% mle.other.Loss$beta.hat ) ) - qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ xbeta <- c( c(1,x) %*% mle.other.Loss$beta.hat ); grad.fit <- c( xbeta^(mle.other.Loss$lambda.hat - 1) * c(1,x), 0, ( log(xbeta) * xbeta^mle.other.Loss$lambda.hat )/mle.other.Loss$lambda.hat - BC.transform( mle.other.Loss$lambda.hat, xbeta )/mle.other.Loss$lambda.hat ); c( grad.fit %*% covML.other.Loss %*% grad.fit ) } ) ), col = 'skyblue' )
lines( BC.transform( mle.other.Loss$lambda.hat, x.seq.Loss ), BC.transform( mle.other.Loss$lambda.hat, c( cbind( 1, x.seq.Loss ) %*% mle.other.Loss$beta.hat ) ) + qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ xbeta <- c( c(1,x) %*% mle.other.Loss$beta.hat ); grad.fit <- c( xbeta^(mle.other.Loss$lambda.hat - 1) * c(1,x), 0, ( log(xbeta) * xbeta^mle.other.Loss$lambda.hat )/mle.other.Loss$lambda.hat - BC.transform( mle.other.Loss$lambda.hat, xbeta )/mle.other.Loss$lambda.hat ); c( grad.fit %*% covML.other.Loss %*% grad.fit ) } ) ), col = 'skyblue' )

#lines( x.seq.Loss, BC.backtransform( mle.other.Loss$lambda.hat, BC.transform( mle.other.Loss$lambda.hat, c( cbind( 1, x.seq.Loss ) %*% mle.other.Loss$beta.hat ) ) - qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ xbeta <- c( c(1,x) %*% mle.other.Loss$beta.hat ); grad.fit <- c( xbeta^(mle.other.Loss$lambda.hat - 1) * c(1,x), 0, ( log(xbeta) * xbeta^mle.other.Loss$lambda.hat )/mle.other.Loss$lambda.hat - BC.transform( mle.other.Loss$lambda.hat, xbeta )/mle.other.Loss$lambda.hat ); grad.fit %*% covML.other.Loss %*% grad.fit } ) ) ), col = 'skyblue' )
#lines( x.seq.Loss, BC.backtransform( mle.other.Loss$lambda.hat, BC.transform( mle.other.Loss$lambda.hat, c( cbind( 1, x.seq.Loss ) %*% mle.other.Loss$beta.hat ) ) + qnorm(0.975) * sqrt( sapply( x.seq.Loss, function(x){ xbeta <- c( c(1,x) %*% mle.other.Loss$beta.hat ); grad.fit <- c( xbeta^(mle.other.Loss$lambda.hat - 1) * c(1,x), 0, ( log(xbeta) * xbeta^mle.other.Loss$lambda.hat )/mle.other.Loss$lambda.hat - BC.transform( mle.other.Loss$lambda.hat, xbeta )/mle.other.Loss$lambda.hat ); grad.fit %*% covML.other.Loss %*% grad.fit } ) ) ), col = 'skyblue' )
#°°° End TBS model °°°#
#######


########
# Alternative estimators for lambda
########

## lambda transforming to a symmetric distribution (p.135 Carroll and Ruppert 1984a)
plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), T.skew, x = x.Loss, y = y.Loss, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Absolute loss" )
abline( h = 0, col = 'red' )
uniroot( f = T.skew, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE )
# robust: lambda_sk = -0.088 or around -1.85, w/o outlier -0.0868

plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), T.skew, x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Absolute loss" )
abline( h = 0, col = 'red' )
uniroot( f = T.skew, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE )
# less robust: lambda_sk = 0.1433 or between -1.75 and -1.95, w/o outlier 0.215


## lambda transforming to a homoskedastic distribution (p.135 Carroll and Ruppert 1984a)
plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), H.hesk, x = x.Loss, y = y.Loss, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Absolute loss" )
abline( h = 0, col = 'red' )
uniroot( f = H.hesk, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE )
# robust: lambda_hesk = 0.101 or between -1.75 and -2, w/o outlier 0.0995

plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), H.hesk, x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE ), type = 'l', xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Absolute loss" )
abline( h = 0, col = 'red' )
uniroot( f = H.hesk, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE )
# less robust: lambda_hesk = 0.0632 or between -1.75 and -2, w/o outlier 0.0558

##############


# Overview plots
par( mfrow = c(1,2) )
plot( Structure$Loss, Contents$Loss, xlab = "Structure", ylab = "Contents", main = "Absolute loss, original scale" )
points( Structure$Loss[idx.highlever], Contents$Loss[idx.highlever], col = 'blue', pch = 19 )
points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], col = 'red', pch = 19 )

plot( Structure$DoL, Contents$DoL, xlab = "Structure", ylab = "Contents", main = "Relative loss, original scale" )
points( Structure$DoL[idx.highlever], Contents$DoL[idx.highlever], col = 'blue', pch = 19 )
points( Structure$DoL[idx.outl], Contents$DoL[idx.outl], col = 'red', pch = 19 )


###################
# Other things in the paper
###################

### Table 1
# First line:
tapply( Contents$Loss, Contents$Region, sum )/( tapply( Contents$Loss, Contents$Region, sum ) + tapply( Structure$Loss, Structure$Region, sum ) )
# Second line:
tapply( Contents$Loss/(Contents$Loss + Structure$Loss), Contents$Region, mean )
tapply( Contents$Loss/(Contents$Loss + Structure$Loss), Contents$Region, median )
# Third line:
tapply( Contents$DoL[-145]/Structure$DoL[-145], Contents$Region[-145], mean ) ## Pb for TI
tapply( Contents$DoL/Structure$DoL, Contents$Region, median )


## Breusch-Paggan test (trials, not complete!!)
e.Loss^2/mle.Loss$sigma2.hat 
1- pf( ( sum( c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% e.Loss^2 )^2 ) - sum( (( diag(n.obs) - H.lambda.Loss ) %*% e.Loss^2)^2 ) )/( sum( (( diag(n.obs) - H.lambda.Loss ) %*% e.Loss^2)^2 )/(n.obs-1) ), df1 = 1, df2 = n.obs - 1 )



sum( c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% (e.Loss^2/mle.Loss$sigma2.hat)^2 )^2 )
sum( (( diag(n.obs) - H.lambda.Loss ) %*% (e.Loss^2/mle.Loss$sigma2.hat)^2)^2 )

/( sum( (( diag(n.obs) - H.lambda.Loss ) %*% (e.Loss^2/mle.Loss$sigma2.hat)^2)^2 )/(n.obs-1) ), df = 1 )

c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% e.Loss^2 ) - 
c( matrix( 1/n.obs, n.obs, n.obs ) %*% e.Loss^2 ) - mean( e.Loss^2 )

solve( crossprod( X.lambda.Loss ) ) %*% t(X.lambda.Loss) %*% e.Loss

1 - pf( ( sum( (y.lambda.Loss - mean(y.lambda.Loss))^2 ) - sum( e.Loss^2 ) )/( sum( e.Loss^2 )/(n.obs - p) ), df1 = 1, df2 = n.obs - p )

### Correlations between y.lambda and x.lambda
cor.test( y.lambda.Loss, x.lambda.Loss )
cor.test( y.lambda.Loss, x.lambda.Loss, method = "spearman" )
cor.test( y.lambda.Loss, x.lambda.Loss, method = "kendall" )

### Adjusted R^2:
1 - ( sum(e.Loss^2)/(n.obs-p) )/( sum( ( y.lambda.Loss - mean(y.lambda.Loss) )^2 )/(n.obs-1) )

