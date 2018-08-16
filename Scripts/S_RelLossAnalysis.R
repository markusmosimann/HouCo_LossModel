setwd("../MOBILAB/Statistik-Hilfe GIUB/Markus_VulnerabilitätFahrhabe")
load("Loss/Loss_Data.RData")
# contains 384 by 5 dataframes "Contents" and "Structure"

source("Global_LossAnalysisLF.R")

plot( Structure$DoL, Contents$DoL, xlab = "Structure", ylab = "Content", main = "Relative Loss" )
points( Structure$DoL[c(145,249,328)], Contents$DoL[c(145,249,328)], pch = 4, col = c('red','red','blue') )
abline( h = 1, lty = 2, col = 'grey' )
abline( v = 1, lty = 2, col = 'grey' )
abline( 0, 1, lty = 4, col = 'grey' )

##########
#==== Fit Box-Cox transformation
##########

p <- 2

# Execute only one of the three possibilities.
#==== All data points ====#
y.DoL <- Contents$DoL
x.DoL <- Structure$DoL

modelname <- "Full"
#==== End ====#


#==== W/o outliers ====#  <--- Good thing to do.
idx.outl <- 145

y.DoL <- Contents$DoL[-idx.outl]
x.DoL <- Structure$DoL[-idx.outl]

modelname <- "noOutl"

PTBS.mle( y = Contents$DoL, x = Structure$DoL )  # as comparison
#==== End ====#


#==== W/o outliers nor lev. pts ====#
idx.remove <- c(145,249) # based on complete data set

y.DoL <- Contents$DoL[-idx.remove]
x.DoL <- Structure$DoL[-idx.remove]

modelname <- "notOutLev"

PTBS.mle( y = Contents$DoL, x = Structure$DoL ) # comparison
PTBS.mle( y = Contents$DoL[-idx.outl], x = Structure$DoL[-idx.outl] ) # comparison
#==== End ====#

n.obs <- length(y.DoL)

sep.lam.global <- FALSE ## Falls gewünscht auf TRUE stellen
if(sep.lam.global) { modelname <- paste( modelname, "_seplam", sep = "" ) }

PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global )

CI.lambda.DoL <- PTBS.CI.lambda( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global )
CI.lambda.DoL
#==# Full: CI = (0.1352,0.2579)
#==# W/o outliers: CI = (0.144,0.265), lambda also slightly shifted to the right
#==# W/o extr.Cdist: CI = (0.151,0.272) even more shifted to the right (as lambda itself)

mle.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global )
# hmmm, seems much better w/ sep.lam
par.DoL <- unlist(mle.DoL[c(2,3,1)])
s2.hat.DoL <- mle.DoL$sigma2.hat * n.obs/(n.obs-p)
# Q: not clear whether n-p or n-p-1 ?!
#--> I would say n-p because the rank of I - H.lambda is still n-p even though lambda is estimated
covML.DoL <- solve( -optimHess( par.DoL, PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ) ) )


# Estimation results
cbind( par.DoL, sqrt(diag(covML.DoL)) )
mle.DoL$ell.opt
2*( -mle.DoL$ell.opt + p+2+sep.lam.global )

# p-value LRT seplam vs. single lambda (in case of sep.lam.global = TRUE)
1 - pchisq( 2*( mle.DoL$ell.opt - PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = FALSE )$ell.opt ), df = 1 )


# Transformed quantities
y.lambda.DoL <- BC.transform( mle.DoL$lambda.hat[1], y.DoL )
x.lambda.DoL <- BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], x.DoL )
X.lambda.DoL <- cbind( 1, x.lambda.DoL )


## Covariance matrix of beta based on regression (#3 below):
covbeta.regr.DoL <- s2.hat.DoL * solve(crossprod(X.lambda.DoL))

# theta.init.DoL <- if(sep.lam.global) { c(0, 1, 0.4, 0.3, 0.5) } else {  c( 0, 1, 0.6, 0.5 ) }
# opt.DoL <- optim( theta.init.DoL, fn = PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ), hessian = TRUE )
# opt.DoL$val - mle.DoL$ell.opt
# covML.DoL <- solve( -opt.DoL$hessian )

## !! sigma2.hat is now correlated because of the estimated lambda:
#cov2cor( covML.DoL )

## Bootstrap covariance matrix (not needed)
# mle.boot.DoL <- t(apply( boot.idx, 1, function(star){ unlist( PTBS.mle( y = y.DoL[star], x = x.DoL[star] )[c(2,3,1)] ) } )) # boot.idx is R by n.obs
# covBoot.DoL <- var( mle.boot.DoL )


###### Covariance relations (can be skipped)
opt.lamfix.DoL <- optim( theta.init.DoL[1:3], fn = regr.llkhd, y = y.lambda.DoL, X = X.lambda.DoL, control = list( fnscale = -1 ), hessian = TRUE )
# Estimating sigma instead of sigma2 (via sigma = TRUE in regr.llkhd) is not so helpful bcs true distn of sigma.hat is only half-normal and not normal (and therefore certainly not equal the asymptotic one).
sqrt(diag(covML.DoL[1:3,1:3]/solve( -opt.lamfix.DoL$hessian ))) # increase factors of se

# J(theta.hat)^{-1} vs. regression (sigma2.hat * (X^T X)^{-1})
solve( -opt.lamfix.DoL$hessian )[1:2,1:2]/( opt.lamfix.DoL$par[3] * solve( crossprod( X.lambda.DoL ) ) )
# Under lamfix: covML.beta = sigma2.hat * (X^T X)^{-1}  (thus w/ biased sigma2.hat)

# estimated lambda vs. fixed lambda (both J(theta.hat)^{-1}):
covML.DoL[1:3,1:3]/solve( -opt.lamfix.DoL$hessian ) 
# --> with estimated lambda much larger variance of sigma2
## having a different lambda for x and y largely increases the variance of beta.hat !!

# cov(beta.hat): s2.hat * (X^T X)^{-1} (regression) vs. J(theta.hat)^{-1} for lambda fixed
covbeta.regr.DoL/solve( -opt.lamfix.DoL$hessian )[1:2,1:2]  # roughly n.obs/(n.obs-p)
# Thus J(theta.hat)^{-1}[1:2,1:2] = sigma2*(X^t X)^{-1} equals the exact covariance of beta.hat in the Gaussian case.
# Thus covbeta.regr = s2.hat * (X^t X)^{-1} equals n.obs/(n.obs-p) * J(theta.hat)^{-1}[1:2,1:2]
# --> Use covbeta.MLadj (?)

# cov(beta.hat) for lambda fixed: J(theta.hat)^{-1} vs. regression version for theta.hat(opt):
solve( -opt.lamfix.DoL$hessian )[1:2,1:2]/( n.obs/(n.obs-p) * opt.lamfix.DoL$par[3] * solve(crossprod( cbind( 1, BC.transform( opt.DoL$par[4+sep.lam.global], x.DoL ) ) )) )

opt.lamfix.DoL$par[3] * solve( crossprod( cbind( 1, BC.transform( opt.DoL$par[4+sep.lam.global], x.DoL ) ) ) )/solve( -opt.lamfix.DoL$hessian )[1:2,1:2]
# --> sigma2.hat(lamfix) {X(lamfix)^t X(lamfix)}^{-1} = J(lamfix)^{-1}
###### End skip

#==== End Fit ====#

#### To be skipped
cbind( estim = c(unlist(mle.DoL)[c(2:4,1)], CI.lambda.DoL$loglik.opt), se = c( sqrt(diag(covbeta.regr.DoL)),rep(NA,3) ) )
cbind( estim = c(opt.DoL$par, opt.DoL$val), se = c( sqrt(diag(covML.DoL)),NA ) )

# One of the three lines at a time:
yhat.lambda.DoL.Full <- yhat.lambda.DoL
yhat.lambda.DoL.noOutl <- yhat.lambda.DoL
yhat.lambda.DoL.noOutLev <- yhat.lambda.DoL

yhat.mat.DoL <- cbind( Full = yhat.lambda.DoL.Full[-idx.remove], noOutl = yhat.lambda.DoL.noOutl[-248], noOutLev = yhat.lambda.DoL.noOutLev )
head(yhat.mat.DoL)

# Absolute difference between fitted values of the three models (on the n.obs-2 common observations):
yhat.diff.mat.DoL <- abs( yhat.mat.DoL %*% cbind( c(1,-1,0), c(1,0,-1), c(0,1,-1) ) )
apply( yhat.diff.mat.DoL, 2, range )
# largest absolute difference is ~0.1 (between Full and noOutLev), for other two combinations the largest absolute difference is ~0.05
#### End skip


##===
# Plot confidence region for lambda
##===

##===== same lambda =====##
proflambda.DoL <- sapply( lam.seq.DoL, PTBS.profllkhd.lambda, y = y.DoL, x = x.DoL, intercept = TRUE )

## Plot profile log-likelihood for lambda
plot( lam.seq.DoL, proflambda.DoL - CI.lambda.DoL[[4]], xlab = expression(lambda), ylab = "Profile log likelihood", type = 'l', xlim = c(-1,1), ylim = c(-200,0) )
abline( h = - qchisq(0.95,1)/2, lty = 2 )
abline( v = CI.lambda.DoL[1:3], lty = 3 )
##===== End =====##


##===== separate lambdas =====##
# Profile likelihood array
lam1.seq.DoL <- seq( 0.001, 0.35, 0.001 )
lam2.seq.DoL <- seq( 0.08, 0.5, 0.001 )
lam.arr.DoL <- expand.grid( lam1 = lam1.seq.DoL, lam2 = lam2.seq.DoL )

proflambda.DoL <- matrix( apply( lam.arr.DoL, 1, PTBS.profllkhd.lambda, y = y.DoL, x = x.DoL ), nr = length(lam1.seq.DoL) )

## Plot confidence region for (lambda.y, lambda.x)
image( x = lam1.seq.DoL, y = lam2.seq.DoL, pmax( proflambda.DoL - (mle.DoL$ell.opt - qchisq(0.95, df=2)/2), 0 ), xlab = expression(lambda[y]), ylab = expression(lambda[x]), breaks = c(-0.0001,0.0001,seq( 0.05, 3, 0.05 )), col = c('white',topo.colors(60)), main = "Relative loss, 95% CR lambdas" )
abline( v = mle.DoL$lambda.hat[1], col = 'red' )
abline( h = mle.DoL$lambda.hat[2], col = 'red' )
abline( v = unlist(CI.lambda.DoL)[c(3,5)], col = 'orange', lty = 2 )
abline( h = unlist(CI.lambda.DoL)[c(4,6)], col = 'orange', lty = 2 )
abline( 0, 1, col = 'darkgrey', lty = 2 )
box()
savePlot( paste( "FiguresLF/DoL_", modelname, "_ConfReg", sep = '' ), type = "pdf" )
# identity line NOT through CR --> lambdas are different
##===== End =====##



##===================##
##=== Diagnostics ===##
##===================##

#x.seq.DoL <- BC.backtransform( 0.2, seq( -5.75, 0.3, 0.025 ) )

## Involved quantities
H.lambda.DoL <- X.lambda.DoL %*% solve( crossprod( X.lambda.DoL ) ) %*% t(X.lambda.DoL)
#sum( diag( H.lambda.DoL ) ) # tr(H) = p = 2
# Fitted values:
yhat.lambda.DoL <- c( H.lambda.DoL %*% y.lambda.DoL )
# Residuals:
e.DoL <- y.lambda.DoL - yhat.lambda.DoL
# Standardised residuals:
r.DoL <- e.DoL/sqrt(s2.hat.DoL*(1-diag(H.lambda.DoL)))
# Deletion residuals (Davison 2008, p.395) = "studentized residuals" in Weisberg 2005 (p.196):
del.resid.DoL <- sqrt( (n.obs-p-1)/(n.obs-p-r.DoL^2) )*r.DoL

## Outliers?
idx.tail <- rev( tail( order( abs(del.resid.DoL) ) ) )
2*( 1 - pt( abs(del.resid.DoL[idx.tail[1]]), df = n.obs - p - 1 ) )
0.05/n.obs  # Bonferroni correction
2*( 1 - pt( abs(del.resid.DoL[idx.tail[2:3]]), df = n.obs - p - 1 ) )
#==# Full: The most extreme one (145) is an outlier, the following not anymore.
#==# W/o outlier: No outliers anymore.
idx.outl <- idx.tail[1]

# Since
# which( abs(r.DoL) > 2 ) # large residual
# which( diag(H.lambda.DoL) > 2*p/n.obs ) # high leverage
# deserve attention,
# which( Cook.DoL > 8/(n.obs - 2*p) ) # ~both
# is worth a closer look.

lrg.resid <- which( abs(r.DoL) > 2 )
length( lrg.resid )
#==# Full: 15
#==# W/o outlier: 16, sep.lam 18

# High leverage (Davison 2008, p.394 top)
idx.highlever <- which( diag(H.lambda.DoL) > 2*p/n.obs )
length( idx.highlever ) # quite a lot
#==# Full: 34
#==# W/o outlier: 33, sep.lam 34
#==# W/o outl nor lev pts:

Cook.DoL <- (r.DoL^2)*diag(H.lambda.DoL)/(p*(1 - diag(H.lambda.DoL)))
extr.Cdist <- which( Cook.DoL > 8/(n.obs-2*p) )
extr.Cdist
#==# Full: 145, 249
#==# W/o outlier: 248; 327 which has high leverage but no large residual is close to limit; outside if sep.lam
#==# W/o outl. nor lev. pts: 326

any(lrg.resid %in% idx.highlever)
extr.Cdist %in% lrg.resid
extr.Cdist %in% idx.highlever
#==# Full: lrg.resid and highlever disjoint, extr.Cdist both have large residuals but not high leverages
#==# W/o outlier: lrg.resid and highlever disjoint, extr.Cdist both has a large residual but no high leverage, one each with sep.lam


## Tukey's one degree of freedom for non-additivity (Davison 2008, p.391)
X.delta.DoL <- cbind( X.lambda.DoL, yhat.lambda.DoL^2 )
lm.delta.DoL <- my.lm( y.lambda.DoL, X.delta.DoL )
SS.delta.DoL <- sum(lm.delta.DoL$e^2)
1 - pf( c( crossprod(e.DoL) - SS.delta.DoL )/lm.delta.DoL$s2, df1 = 1, df2 = lm.delta.DoL$n - lm.delta.DoL$p )  # equivalent to t-test for delta
#==# Full: 0.01173 (arrgh!!) --> Maybe because of the outlier??
#==# W/o outlier: 0.01836 (arrgh!!) --> Apparently NOT because of the outlier...!!
#==# W/o outl. nor lev. pts: 0.0219 (arrgh!!) --> Apparently NOT because of outlier nor leverage point...!!
#==# W/o outl. nor lev. pts and lambda = 0.25: 0.0658 (petit ouf)


#======
#== Plots
#======

## Residuals
x11( width = 10.7, height = 10.5 )
par(mfrow = c(2,2))
dev.set(2)
# Normality of standardised residuals
qqnorm( r.DoL )
abline( 0, 1, col = 'red' )

# Std residuals vs fitted values: Fit doesn't seem too good (see below)!!
plot( yhat.lambda.DoL, r.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
abline( h = 0, col = 'red' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( yhat.lambda.DoL[extr.Cdist], r.DoL[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.DoL[idx.highlever], r.DoL[idx.highlever], col = 'blue' )
# quite obvious pattern
#==# W/o outl nor lev pts: Extreme Cook's distance has also high leverage
# High leverage not very interesting since just largest and smallest values of x...

# Deletion or studentized residuals (plot not interesting)
#plot( del.resid.DoL, ylab = "Deletion residuals", main = "Relative loss, transformed scale" )
#points( lrg.resid, del.resid.DoL[lrg.resid], col = 'sandybrown', pch = 19 )
# High leverage points not special in this plot

# Absolute std residuals vs. fitted values (see heteroscedasticity)
plot( yhat.lambda.DoL, abs(r.DoL), xlab = "Fitted values", ylab = "", main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
points( yhat.lambda.DoL[extr.Cdist], abs(r.DoL)[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.DoL[idx.highlever], abs(r.DoL)[idx.highlever], col = 'blue' )
abline( h = 2, col = 'sandybrown', lty = 2 )
#
plot( yhat.lambda.DoL, sqrt(abs(r.DoL)), xlab = "Fitted values", ylab = "", main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
points( yhat.lambda.DoL[extr.Cdist], sqrt(abs(r.DoL[extr.Cdist])), col = 'red', pch = 19 )
points( yhat.lambda.DoL[idx.highlever], sqrt(abs(r.DoL[idx.highlever])), col = 'blue' )
abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
#text( -1.5, 2, labels = expression( sqrt(group("|","Standardised residuals","|")) ) )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )

#==# Full: Fit maybe not that good, but apart from outlier quite homoskedastic (maybe smaller variance for large x-values?).
#==# W/o outlier: Still some problem of fit. Smaller variance for high leverage points with large x??

# Running variance along fitted values for block sizes 20, 26, 32
x11( width = 9.8, height = 6.7 )
dev.set(5)
plot( sort(yhat.lambda.DoL)[9+seq(1,365,4)], sapply( seq(1,365,4), function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:19], na.rm = TRUE ) } ), type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) ) # block size 20
lines( sort(yhat.lambda.DoL)[12+seq(1,357,4)], c( sapply( seq(1,353,4), function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:25], na.rm = TRUE ) } ), var( r.DoL[order(yhat.lambda.DoL)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
lines( sort(yhat.lambda.DoL)[15+seq(1,353,4)], sapply( seq(1,353,4), function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
abline( h = 1, col = 'darkgrey', lty = 2 )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

# Compare residual plots for the two models (ONLY if separate lambda)
x11( width = 6.6, height = 7 )
dev.set(6)
#
plot( yhat.lambda.DoL, r.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = "Relative loss, separate lambda" )
abline( h = 0, col = 'red' )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_compare", sep = '' ), type = "pdf" )


## Leverage etc.

# Cook's distance (plot not interesting)
#plot( Cook.DoL, xlab = "Index", ylab = "Cook's distance", main = "Relative loss, transformed scale" )
#abline( h = 8/(n.obs - 2*p), col = 'red' )
# High leverage points not special in this plot

x11( width = 13.2, height = 7 )
par(mfrow = c(1,2))
dev.set(3)
# Distinguish outliers from leverage points (Davison 2008, p.395)
plot( diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)), Cook.DoL, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", main = "Relative loss, transformed scale" )
points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[extr.Cdist], Cook.DoL[extr.Cdist], col = 'red', pch = 19 )
points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[idx.highlever], Cook.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
abline( h = 8/(n.obs - 2*p), col = 'red' )
points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[lrg.resid], Cook.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
abline( 0, 4/p, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
abline( v = 2*p/(n.obs - 2*p), col = 'blue' )

# Residuals vs. leverage (like plot.lm but with better level curves)
plot( diag(H.lambda.DoL), r.DoL, xlab = "Leverage", ylab = "Standardised residuals", main = "Relative loss, transformed scale", xlim = c(0,max(diag(H.lambda.DoL))) ) # xlim = c(0,0.025)
#plot( diag(H.lambda.DoL), r.DoL, xlab = "Leverage", ylab = "Standardised residuals", main = "Relative loss, transformed scale", xlim = c(0,0.025), ylim = c(-3,3) )
abline( h = 0, col = 'grey' )
abline( v = 2*p/n.obs, col = 'blue' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( diag(H.lambda.DoL)[extr.Cdist], r.DoL[extr.Cdist], col = 'red', pch = 19 )
points( diag(H.lambda.DoL)[idx.highlever], r.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
points( diag(H.lambda.DoL)[lrg.resid], r.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p)) * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'orange', lty = 3 )
# r^2 h / (p*(1-h)) = C
# r^2 = C * p*(1-h)/h
savePlot( paste( "FiguresLF/DoL_", modelname, "_leverage", sep = "" ), type = "pdf" )

head(rev(order(Cook.DoL)))


## Regression equation (response vs. covariate)
x11()
par(mfrow = c(1,1))
dev.set(4)
plot( x.lambda.DoL, y.lambda.DoL, xlab = "Structure", ylab = "Contents", main = "Relative loss, transformed scale" )
abline( mle.DoL$beta.hat, col = 'red' )
points( x.lambda.DoL[extr.Cdist], y.lambda.DoL[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.DoL[idx.highlever], y.lambda.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
points( x.lambda.DoL[lrg.resid], y.lambda.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
lines( lowess( x = x.lambda.DoL, y = y.lambda.DoL ), col = 'green' )
#==# Full: --> The first red point is also an outlier, but I don't know why the second is a leverage point.
#==# W/o outlier: --> Still the (not really explainable) leverage point left. Point with smallest x has no high leverage anymore.
#==# W/o outl nor lev pts: Nothing left.
savePlot( paste( "FiguresLF/DoL_", modelname, "_oview", sep = "" ), type = "pdf" )


# Output of non-additivity model:
lm.delta.DoL$coefs
lm.delta.DoL$se.coefs
plot( c(X.delta.DoL %*% lm.delta.DoL$coefs), lm.delta.DoL$e/sqrt( lm.delta.DoL$s2 * (1-lm.delta.DoL$levs) ), xlab = "Fitted values", ylab = "Standardised residuals", main = "Relative loss, with additivity term" ) 
abline( h = 0, col = 'red' )
# Pattern of residuals basically the same as before. But this is not meaningful because the non-linear dependence doesn't need to have form of yhat^2.

## Try with x.lambda^2 as additional covariate
X.sq.lambda.DoL <- cbind( X.lambda.DoL, x.lambda.DoL^2 )
lm.sq.DoL <- my.lm( y.lambda.DoL, X.sq.lambda.DoL )
SS.sq.lambda.DoL <- sum(lm.sq.DoL$e^2)
lm.sq.DoL$coefs
lm.sq.DoL$se.coefs
plot( c(X.sq.lambda.DoL %*% lm.sq.DoL$coefs), lm.sq.DoL$e/sqrt( lm.sq.DoL$s2 * (1-lm.sq.DoL$levs) ), xlab = "Fitted values", ylab = "Standardised residuals", main = "Relative loss, with squared covariate" )
abline( h = 0, col = 'red' )
# Not really helpful, as before.

# Again the overall plot (response vs. covariate)
plot( x.lambda.DoL, y.lambda.DoL, xlab = "Structure", ylab = "Contents", main = "Relative loss, transformed scale" )
abline( h = 0, lty = 2, col = 'grey' )
abline( v = 0, lty = 2, col = 'grey' )
abline( 0, 1, lty = 4, col = 'grey' )
abline( mle.DoL$beta.hat, col = 'red' )
lines( lowess( x = x.lambda.DoL, y = y.lambda.DoL ), col = 'blue', lwd = 1.2 )
#points( x.lambda.DoL[idx.outl], y.lambda.DoL[idx.outl], col = 'red', pch = 19 )
points( x.lambda.DoL[extr.Cdist], y.lambda.DoL[extr.Cdist], col = 'red', pch = 19 )
rug( x.lambda.DoL, quiet = TRUE, ticksize = 0.02 )
abline( v = c(-2.7,-2.75,-2.8), col = 'blue', lty = 2)
#--> Try a model with two piecewise linear functions (1 knot)

#=====
# Fit piecewise linear function (not good idea)
#=====
modelname <- paste( modelname, "_split", sep = '' )
p.split <- p+1

# Optim is not very stable, in particular in the subsequent steps:
opt.knot.DoL <- optim( c( -1, 0.5, 1, 0.6, 0.4, -2.5 ), fn = PTBS.llkhd.knot, y = y.DoL, x = x.DoL, control = list( fnscale = -1, maxit = 10000 ), hessian = TRUE )
opt.knot.DoL$par
sqrt( diag(solve( -opt.knot.DoL$hessian )) )
#==# Full: suggests taking knot at about -2.8, se = 0.04
#==# W/o outlier: suggests taking knot at about -2.8, se = 0.055

# Use profiling instead
### Estimate lambda again
modelname <- paste( modelname, "_lamgen", sep = '' )
prof.knot.DoL <- sapply( seq(-4,-1.5,0.01), function(kn){ optimize( f = PTBS.profllkhd.knot.lambda, interval = c(-3,3), y = y.DoL, x = x.DoL, knot.x = kn, intercept = TRUE, maximum = TRUE )$objective } )
plot( seq(-4,-1.5,0.01), prof.knot.DoL )

knot.DoL <- optimize( f = PTBS.profllkhd.knot.kn, interval = c(-4,-1.5), y = y.DoL, x = x.DoL, intercept = TRUE, maximum = TRUE )$maximum
abline( v = knot.DoL, col = 'red' )

mle.split.DoL <- PTBS.knot.mle( y = y.DoL, x = x.DoL, knot = knot.DoL )

PTBS.knot.CI.lambda( y = y.DoL, x = x.DoL, knot = knot.DoL )
#==# Full: CI = (0.1166,0.2118) --> overlap
#==# W/o outlier: CI = (0.1253,0.2211) --> overlap
## Thus could keep same lambda for both cases.

covML.split.DoL <- solve( -optimHess( unlist(mle.split.DoL[c(2,3,1)]), PTBS.llkhd.knot, theta = c(rep(NA,5), knot.DoL), y = y.DoL, x = x.DoL, control = list( maxit = 5000, fnscale = -1 ) ) ) # knot considered as fixed

### Keep lambda fixed (not sure this is so useful) !!! reset lrg.resid, idx.highlever and extr.Cdist to non-split model !!!
lrg.resid <- which( abs(r.DoL) > 2 )
idx.highlever <- which( diag(H.lambda.DoL) > 2*p/n.obs )
extr.Cdist <- which( Cook.DoL > 8/(n.obs-2*p) )
modelname <- paste( unlist( strsplit( modelname, "_lam" ) )[1], "_lamfix", sep = '' )

prof.knot.DoL <- sapply( seq(-4,-1.5,0.01), function(kn){ PTBS.profllkhd.knot.lambda( lambda = mle.DoL$lambda.hat, y = y.DoL, x = x.DoL, knot.x = kn, intercept = TRUE ) } )
plot( seq(-4,-1.5,0.01), prof.knot.DoL )

knot.DoL <- optimize( f = PTBS.profllkhd.knot.kn, interval = c(-4,-1.5), y = y.DoL, x = x.DoL, lambda = mle.DoL$lambda.hat, intercept = TRUE, maximum = TRUE )$maximum
mle.split.DoL <- PTBS.knot.mle( y = y.DoL, x = x.DoL, lambda = mle.DoL$lambda.hat, knot = knot.DoL )
abline( v = knot.DoL, col = 'red' )

covML.split.DoL <- solve( -optimHess( unlist(mle.split.DoL[c(2,3)]), PTBS.llkhd.knot, theta = c(rep(NA,4), mle.split.DoL$lambda.hat, knot.DoL), y = y.DoL, x = x.DoL, control = list( maxit = 5000, fnscale = -1 ) ) ) # knot considered as fixed, lambda fixed as well
###

x.lambda.split.DoL <- BC.transform( mle.split.DoL$lambda.hat, x.DoL )

y.lambda.split.DoL <- BC.transform( mle.split.DoL$lambda.hat, y.DoL )
X.lambda.split.DoL <- cbind( 1, pmin( x.lambda.split.DoL, knot.DoL ), pmax( x.lambda.split.DoL - knot.DoL, 0 ) )

H.lambda.split.DoL <- X.lambda.split.DoL %*% solve( crossprod(X.lambda.split.DoL) ) %*% t(X.lambda.split.DoL)

cbind( unlist(mle.split.DoL[c(2,3,1)]), sqrt( diag(covML.split.DoL) ) ) # sqrt( n.obs*mle.split.DoL$sigma2.hat/(n.obs - p.split) * diag( solve( crossprod(X.lambda.split.DoL) ) ) )
mle.split.DoL$ell.opt

# Model comparison
2*( -mle.split.DoL$ell.opt + p.split+2+grepl("lamgen",modelname) )
2*( -mle.DoL$ell.opt + p+2 ) # as comparison
# AIC suggests that split-model is better
1 - pchisq( 2*(mle.split.DoL$ell.opt - mle.DoL$ell.opt), df = 1 )  # LRT (also valid if lambda estimated)
#==# Full: 0.00076 (lambda fixed)
#==# W/o outlier: 0.0013 (lambda fixed)

# Fitted values:
yhat.lambda.split.DoL <- c( H.lambda.split.DoL %*% y.lambda.split.DoL )
# Residuals:
e.split.DoL <- y.lambda.split.DoL - yhat.lambda.split.DoL
# Standardised residuals:
r.split.DoL <- e.split.DoL/sqrt( (n.obs/(n.obs-p.split))*mle.split.DoL$sigma2.hat * (1-diag(H.lambda.split.DoL)) )
# Deletion or studentized residuals:
del.resid.split.DoL <- sqrt( (n.obs - p.split - 1)/(n.obs - p.split - r.split.DoL^2) )*r.split.DoL

# Fit model with two different sigma2 but same knot:
# pars.curr <- c( -1.2669356, 0.4459105, 0.9109834, 0.2871712, 0.3222869, 0.1670669) # full, lamgen
pars.curr <- unlist(mle.split.DoL[c(2,3,3,if(grepl("lamgen",modelname)){ 1 } else { NULL })])
intvl <- list( c(-5,0), c(0,1.5), c(0,2.5), c(0,2), c(0,2), if(grepl("lamgen",modelname)){ c(-3,3) } else { NULL } ) # last element is present even if NULL, but this doesn't matter
R <- length(pars.curr)
ell.vec <- rep(NA,R)
repeat{
	for( r in 1:R ) {
		pars <- pars.curr
		pars[r] <- NA
		opt.r <- optimize( f = PTBS.llkhd.knot, interval = intvl[[r]], theta = c( pars, if(grepl("lamfix",modelname)){ mle.split.DoL$lambda.hat } else { NULL }, knot.DoL ), y = y.DoL, x = x.DoL, var.split = TRUE, maximum = TRUE )
		pars.curr[r] <- opt.r$maximum
		ell.vec[r] <- opt.r$objective
	}
	if( max( abs( combn( R, 2, function(i){ diff(ell.vec[i]) } ) ) ) < 1e-07 ) {
		break
	}
}
ell.vec
pars.curr
1 - pchisq( 2*(ell.vec[R] - mle.split.DoL$ell.opt), df = 1 )
# Two different sigma2 doesn't look like an improvement.
# W/o outlier for lambda fixed, the p-value is 0.040 !!!


# Other diagnostics for split model

## Outliers?
idx.tail <- rev( tail( order( abs(del.resid.split.DoL) ) ) )
2*( 1 - pt( abs(del.resid.split.DoL[idx.tail[1]]), df = n.obs - p.split - 1 ) )
0.05/n.obs  # Bonferroni correction
2*( 1 - pt( abs(del.resid.split.DoL[idx.tail[2:3]]), df = n.obs - p.split - 1 ) )
idx.tail[1]
#==# Full: The most extreme one (145) is an outlier, the following not anymore (lambda general and fixed)
#==# W/o outlier: No outliers anymore (lambda general).
#==# W/o outl nor lev pts: No outliers anymore.

lrg.resid <- which( abs(r.split.DoL) > 2 )
length( lrg.resid )
#==# Full: 17 if lambda general and fixed
#==# W/o outlier: 18 if lambda general, 20 if lambda fixed

# High leverage (Davison 2008, p.394 top)
idx.highlever <- which( diag(H.lambda.split.DoL) > 2*p.split/n.obs )
length( idx.highlever ) # quite a lot
#==# Full: 47 if lambda general, 49 if lambda fixed
#==# W/o outlier: 47 if lambda general, 49 if lambda fixed
#==# W/o outl nor lev pts: 33
#==# W/o outl nor lev pts and lambda = 0.25: 33

Cook.split.DoL <- (r.split.DoL^2)*diag(H.lambda.split.DoL)/(p.split*(1 - diag(H.lambda.split.DoL)))
extr.Cdist <- which( Cook.split.DoL > 8/(n.obs-2*p.split) )
extr.Cdist
#==# Full: 145 only, but 328 (which has high leverage) is very close to limit for lambda general; 328 also is outside if lambda is fixed
#==# W/o outlier: 327 (both lambda general and fixed)
#==# W/o outl. nor lev. pts: none
#==# W/o outl. nor lev. pts and lambda = 0.25: 326

any(lrg.resid %in% idx.highlever)
extr.Cdist %in% lrg.resid
extr.Cdist %in% idx.highlever
#==# Full: lrg.resid and highlever disjoint (both lambda general and fixed); extr.Cdist has large residual but not high leverage if lambda general, if lambda fixed extr.Cdist is one each
#==# W/o outlier: lrg.resid and highlever disjoint, extr.Cdist has high leverage but no large residual (lambda general and fixed)

# Tukey's one degree of freedom for non-additivity (Davison 2008, p.391)
X.delta.split.DoL <- cbind( X.lambda.split.DoL, yhat.lambda.split.DoL^2 )
H.delta.split.DoL <- X.delta.split.DoL %*% solve( crossprod(X.delta.split.DoL) ) %*% t(X.delta.split.DoL)
SS.delta.split.DoL <- c( y.lambda.split.DoL %*% ( diag(n.obs) - H.delta.split.DoL ) %*% y.lambda.split.DoL )
1 - pf( c( crossprod(e.split.DoL) - SS.delta.split.DoL )/( SS.delta.split.DoL/(n.obs - p.split - 1) ), df1 = 1, df2 = n.obs-p.split-1 )  # equivalent to t-test for delta
#==# Full: 0.4282 for lambda general, 0.3939 for lambda fixed (ouf!)
#==# W/o outlier: 0.4026 for lambda general, 0.3736 for lambda fixed (ouf!)


#=====
#== Plots
#=====
plot( x.lambda.split.DoL, y.lambda.split.DoL, xlab = "Structure", ylab = "Contents", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( h = 0, lty = 2, col = 'grey' )
abline( v = 0, lty = 2, col = 'grey' )
abline( 0, 1, lty = 4, col = 'grey' )
lines( lowess( x = x.lambda.split.DoL, y = y.lambda.split.DoL ), col = 'green' )
abline( c( solve( crossprod( cbind( 1, x.lambda.split.DoL ) ) ) %*% t( cbind( 1, x.lambda.split.DoL ) ) %*% y.lambda.split.DoL ), col = 'red' ) # linear fit for current lambda
#points( x.lambda.split.DoL[idx.outl], y.lambda.split.DoL[idx.outl], col = 'red', pch = 19 )
points( x.lambda.split.DoL[extr.Cdist], y.lambda.split.DoL[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.split.DoL[idx.highlever], y.lambda.split.DoL[idx.highlever], col = 'blue' )
points( x.lambda.split.DoL[lrg.resid], y.lambda.split.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
rug( x.lambda.split.DoL, quiet = TRUE, ticksize = 0.02 )
lines( xlam.seq.DoL, c( cbind( 1, pmin( xlam.seq.DoL, knot.DoL ), pmax( xlam.seq.DoL - knot.DoL, 0 ) ) %*% mle.split.DoL$beta.hat ), col = 'green4' )
abline( v = knot.DoL, col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/DoL_", modelname, "_oview1", sep = '' ), type = "pdf" )

# Compare residual plots for the two models
x11( width = 13.2, height = 7 )
par( mfrow = c(1,2) )
dev.set(3)
#
plot( yhat.lambda.DoL, r.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = "Relative loss, initial model" )
abline( h = 0, col = 'red' )
abline( v = c( mle.DoL$beta.hat %*% c( 1, BC.transform( mle.DoL$lambda.hat, BC.backtransform( mle.split.DoL$lambda.hat, knot.DoL ) ) ) ), col = 'green4', lty = 3 )
#
plot( yhat.lambda.split.DoL, r.split.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( h = 0, col = 'red' )
abline( v = c( mle.split.DoL$beta.hat[1:2] %*% c(1,knot.DoL) ), col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_compare", sep = '' ), type = "pdf" )
#==# Full: Split model looks slightly better (both lambda general and fixed)
#==# W/o outlier: Split model looks slightly better (lambda general and fixed)


## Residuals
par(mfrow = c(2,2))
dev.set(2)
# Normality of standardised residuals
qqnorm( r.split.DoL )
abline( 0, 1, col = 'red' )

# Std residuals vs fitted values:
plot( yhat.lambda.split.DoL, r.split.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( h = 0, col = 'red' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
points( yhat.lambda.split.DoL[extr.Cdist], r.split.DoL[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.split.DoL[idx.highlever], r.split.DoL[idx.highlever], col = 'blue' )
# quite obvious pattern
#==# Full: much more high leverages for small yhat (and thus x) if lambda general
#==# W/o outl nor lev pts: Extreme Cook's distance has also high leverage
# High leverage not very interesting since just largest and smallest values of x...

# Deletion or studentized residuals (plot not interesting)
#plot( del.resid.split.DoL, ylab = "Deletion residuals", main = "Relative loss, transformed scale" )
# High leverage points not special in this plot

# Absolute std residuals vs. fitted values (see heteroscedasticity)
plot( yhat.lambda.split.DoL, abs(r.split.DoL), xlab = "Fitted values", ylab = "", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
points( yhat.lambda.split.DoL[extr.Cdist], abs(r.split.DoL)[extr.Cdist], col = 'red', pch = 19 )
points( yhat.lambda.split.DoL[idx.highlever], abs(r.split.DoL)[idx.highlever], col = 'blue' )
abline( h = 2, col = 'sandybrown', lty = 2 )
#
plot( yhat.lambda.split.DoL, sqrt(abs(r.split.DoL)), xlab = "Fitted values", ylab = "", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
points( yhat.lambda.split.DoL[extr.Cdist], sqrt(abs(r.split.DoL[extr.Cdist])), col = 'red', pch = 19 )
points( yhat.lambda.split.DoL[idx.highlever], sqrt(abs(r.split.DoL[idx.highlever])), col = 'blue' )
abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
#text( -1.5, 2, labels = expression( sqrt(group("|","Standardised residuals","|")) ) )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )

#==# Full: fit and homoskedasticity reasonable (lambda general)
#==# W/o outlier: fit & homosk. not too bad, smaller variance for high-leverage points with large x?

# Running variance along fitted values for block sizes 20, 26, 32
x11( width = 9.8, height = 6.7 )
dev.set(5)
plot( sort(yhat.lambda.split.DoL)[9+seq(1,365,4)], sapply( seq(1,365,4), function(i){ var( r.split.DoL[order(yhat.lambda.split.DoL)][i+0:19], na.rm = TRUE ) } ), type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = "Relative loss, transformed scale" ) # block size 20
lines( sort(yhat.lambda.split.DoL)[12+seq(1,357,4)], c( sapply( seq(1,353,4), function(i){ var( r.split.DoL[order(yhat.lambda.split.DoL)][i+0:25], na.rm = TRUE ) } ), var( r.split.DoL[order(yhat.lambda.split.DoL)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
lines( sort(yhat.lambda.split.DoL)[15+seq(1,353,4)], sapply( seq(1,353,4), function(i){ var( r.split.DoL[order(yhat.lambda.split.DoL)][i+0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
abline( h = 1, col = 'darkgrey', lty = 2 )
savePlot( paste( "FiguresLF/DoL_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

## Leverage etc.
par(mrow = c(1,2))
dev.set(3)
# Cook's distance (plot not interesting)
#plot( Cook.split.DoL, xlab = "Index", ylab = "Cook's distance" )
#abline( h = 8/(n.obs - 2*p.split), col = 'red' )
# High leverage points not special in this plot

# Distinguish outliers from leverage points (Davison 2008, p.395)
plot( diag(H.lambda.split.DoL)/(1-diag(H.lambda.split.DoL)), Cook.split.DoL, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
points( (diag(H.lambda.split.DoL)/(1-diag(H.lambda.split.DoL)))[extr.Cdist], Cook.split.DoL[extr.Cdist], col = 'red', pch = 19 )
points( (diag(H.lambda.split.DoL)/(1-diag(H.lambda.split.DoL)))[idx.highlever], Cook.split.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
abline( h = 8/(n.obs - 2*p.split), col = 'red' )
points( (diag(H.lambda.split.DoL)/(1-diag(H.lambda.split.DoL)))[lrg.resid], Cook.split.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
abline( 0, 4/p.split, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
abline( v = 2*p.split/(n.obs - 2*p.split), col = 'blue' )

# Residuals vs. leverage (like plot.lm but with better level curves)
plot( diag(H.lambda.split.DoL), r.split.DoL, xlab = "Leverage", ylab = "Standardised residuals", main = paste( "Relative loss, split model", tail( unlist( strsplit( modelname, "_" ) ), 1 ) ), xlim = c(0,max(diag(H.lambda.split.DoL))) ) # xlim = c(0,0.025)
#plot( diag(H.lambda.split.DoL), r.split.DoL, xlab = "Leverage", ylab = "Standardised residuals", main = "Relative loss, transformed scale", xlim = c(0,0.025), ylim = c(-3,3) )
abline( h = 0, col = 'grey' )
abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
abline( v = 2*p.split/n.obs, col = 'blue' )
points( diag(H.lambda.split.DoL)[extr.Cdist], r.split.DoL[extr.Cdist], col = 'red', pch = 19 )
points( diag(H.lambda.split.DoL)[idx.highlever], r.split.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
points( diag(H.lambda.split.DoL)[lrg.resid], r.split.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p.split)) * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p.split*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'orange', lty = 3 )
# r^2 h / (p.split*(1-h)) = C
# r^2 = C * p.split*(1-h)/h
savePlot( paste( "FiguresLF/DoL_", modelname, "_leverage", sep = "" ), type = "pdf" )

head(rev(order(Cook.split.DoL)))

dev.set(4)
plot( x.lambda.split.DoL, y.lambda.split.DoL, xlab = "Structure", ylab = "Contents", main = paste( "Relative loss, split model",  tail( unlist( strsplit( modelname, "_" ) ), 1 ) ) )
abline( h = 0, lty = 2, col = 'grey' )
abline( v = 0, lty = 2, col = 'grey' )
abline( 0, 1, lty = 4, col = 'grey' )
lines( lowess( x = x.lambda.split.DoL, y = y.lambda.split.DoL ), col = 'green' )
abline( c( solve( crossprod( cbind( 1, x.lambda.split.DoL ) ) ) %*% t( cbind( 1, x.lambda.split.DoL ) ) %*% y.lambda.split.DoL ), col = 'red' ) # linear fit for current lambda
points( x.lambda.split.DoL[extr.Cdist], y.lambda.split.DoL[extr.Cdist], col = 'red', pch = 19 )
points( x.lambda.split.DoL[idx.highlever], y.lambda.split.DoL[idx.highlever], col = 'blue' )
points( x.lambda.split.DoL[lrg.resid], y.lambda.split.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
rug( x.lambda.split.DoL, quiet = TRUE, ticksize = 0.02 )
lines( xlam.seq.DoL, c( cbind( 1, pmin( xlam.seq.DoL, knot.DoL ), pmax( xlam.seq.DoL - knot.DoL, 0 ) ) %*% mle.split.DoL$beta.hat ), col = 'green4' )
abline( v = knot.DoL, col = 'green4', lty = 3 )
savePlot( paste( "FiguresLF/DoL_", modelname, "_oview2", sep = '' ), type = "pdf" )

#====== End piecewise linear function =======#


######
# CI and PI on transformed scale
######

# There are three possible scenarios under which standard errors of prediction x^T*beta.hat can be computed:
# 1. Delta method on full covariance matrix: Account for estimation of lambda. This is not proper because uncertainty of lambda in y.lambda is not accounted for.
# 2. Use std errors for beta from full covariance matrix (thus estimation of lambda is accounted for in se.beta), but then ignore that lambda was estimated when appyling the delta method for var(x0^T beta.hat).
# 3. Regression-based while treating lambda = lambda.hat as fixed. Thus the uncertainty of estimating lambda is contained nowhere, not even in the std errors of beta.
# 4. As 2 but bootstrap-based.
# 5. Detour via original scale by computing confidence bands for y.hat (p-quantile) (accounting for uncertainty in lambda), which are then transformed again.

# --> I have two (very similar) sets of estimated parameters (mle and opt). Based on my calculations, I could also code the complete observed information matrix for mle, but currently I haven't done it. Therefore I have to used the opt-estimates for full uncertainty.
# --> By applying optimHess to mle.DoL I get the Hessian and thus the joint covariance at mle.DoL. Thus I do not need the opt version anymore. Since 2 seems better than 3, I keep 2 and 5.

fitted.mle.seq.DoL <- c( cbind( 1, xlam.seq.DoL ) %*% mle.DoL$beta.hat )
#fitted.opt.seq.DoL <- c( cbind( 1, xlam.seq.DoL ) %*% opt.DoL$par[1:p] )
#max( abs( fitted.mle.seq.DoL - fitted.opt.seq.DoL ) )

# (2) Assuming lambda = lambda.hat known, based on mle
var.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ c(1,x) %*% covML.DoL[1:p,1:p] %*% c(1,x) } )
var.predict.seq.DoL <- var.fitted.seq.DoL + s2.hat.DoL

# (3) Cov matrix of beta.hat from the regression of y.lambda on X.lambda with lambda = lambda.hat fixed
varRegr.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ c(1,x) %*% covbeta.regr.DoL %*% c(1,x) } )

# (5) Detour via original scale to account for estimated lambda in y as well (but not sure it is properly done)
varOrig.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(xlam){ grad.yhat( theta = par.DoL, x.lam = c(1,xlam), sep.lam = sep.lam.global ) %*% covML.DoL %*% grad.yhat( theta = par.DoL, x.lam = c(1,xlam), sep.lam = sep.lam.global ) } )

# (1) Including uncertainty on lambda.hat in x.lambda, based on opt
varFull.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ grad.fitted( par.DoL, c(1,x) ) %*% covML.DoL %*% grad.fitted( par.DoL, c(1,x) ) } ) # sigma2 has no influence on this
varFull.predict.seq.DoL <- varFull.fitted.seq.DoL + s2.hat.DoL

# (4) As 1 but with cov(beta) based on non-parametric bootstrap 
varBoot.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ grad.fitted( par.DoL, c(1,x) ) %*% covBoot.DoL %*% grad.fitted( par.DoL, c(1,x) ) } )
varBoot.predict.seq.DoL <- varBoot.fitted.seq.DoL + s2.hat.DoL

x11()
plot( xlam.seq.DoL, varOrig.fitted.seq.DoL, type = 'l', xlab = "x(lambda)", ylab = "Variance of X(lambda)*beta.hat", main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
#lines( xlam.seq.DoL, varBoot.fitted.seq.DoL, col = 'dodgerblue' )
lines( xlam.seq.DoL, var.fitted.seq.DoL, col = 'blue' )
lines( xlam.seq.DoL, varRegr.fitted.seq.DoL, col = 'cyan1' )
dev.off()

x11()
# Transformed scale: y.lambda vs. X.lambda*beta
plot( x.lambda.DoL, y.lambda.DoL, xlab = "Structure", ylab = "Contents", main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
polygon( x = BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], c( 1.5, 1, 1, 0, 0, 1.5 ) ), y = BC.transform( mle.DoL$lambda.hat[1], c( 0, 0, 1, 1, 1.5, 1.5 ) ), density = 5, angle = -45, col = "lightgrey")
# if appropriate:
points( BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], Structure$DoL[idx.outl] ), BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[idx.outl] ), pch = 4 ) # ev. pch = 8
points( BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], Structure$DoL[idx.remove] ), BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[idx.remove] ), pch = 4 )
#
abline( mle.DoL$beta.hat, col = 'red' )
#abline( opt.DoL$par[1:p], col = 'lightpink2' )

# 95% Confidence interval (mle, lambda fixed --> 3)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ), col = 'orange', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ), col = 'orange', lty = 2 )

# 95% Prediction interval (mle, lambda fixed --> 3)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ), col = 'cyan1', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ), col = 'cyan1', lty = 2 )

# 95% Confidence interval (mle, lambda known --> 2)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ), col = 'red', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ), col = 'red', lty = 2 )

# 95% Prediction interval (mle, lambda known --> 2)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ), col = 'blue', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ), col = 'blue', lty = 2 )

# 95% Confidence interval (opt, via original scale --> 5)
lines( xlam.seq.DoL, BC.transform( mle.DoL$lambda.hat[1], BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL) ), col = 'sienna1', lty = 2 )
lines( xlam.seq.DoL, BC.transform( mle.DoL$lambda.hat[1], BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) + qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL) ), col = 'sienna1', lty = 2 )
# Prediction interval accounting for estimation uncertainty is too complex for this case
lines( lowess( x = x.lambda.DoL, y = y.lambda.DoL ), col = 'green' )
box()

# 95% Confidence interval (opt, lambda unknown --> 1)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ), col = 'lightpink2', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ), col = 'lightpink2', lty = 2 )  # use qnorm here???

# 95% Prediction interval (opt, lambda unknown --> 1)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ), col = 'deepskyblue', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ), col = 'deepskyblue', lty = 2 )

# 95% Confidence interval (mle, lambda unknown --> 4)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.DoL ), col = 'chocolate3', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.DoL ), col = 'chocolate3', lty = 2 )  # use qnorm here???

# 95% Prediction interval (mle, lambda unknown --> 4)
lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.DoL ), col = 'dodgerblue', lty = 2 )
lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.DoL ), col = 'dodgerblue', lty = 2 )
# --> Much less difference for prediction intervals because sigma2 dominates the variance.


## Either 2 or 3 seem reasonable since very close. 1 and 4 are not plausible because uncertainty of lambda is accounted for in x, but not in y. And the detour via the original scale is very close to 2 or 3.

( BC.transform( opt.DoL$par[4], BC.backtransform( opt.DoL$par[4], fitted.opt.seq.DoL ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL) ) - ( fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ) )
covML.DoL[1:2,1:2]/covbeta.regr.DoL
########


########
# Backtransformed to original scale (!! only for full data !!)
########

x.seq.DoL <- BC.backtransform( mle.DoL$lambda.hat[1+sep.lam.global], xlam.seq.DoL )

## Conditional mean on original scale, estimator by Taylor 1986
condmean.seq.DoL <- condmean.orig( theta = par.DoL, x = cbind(1,xlam.seq.DoL), sep.lam = sep.lam.global )
## use sigma2 or s2 here???
var.condmean.DoL <- sapply( xlam.seq.DoL, function(x){ c( grad.condmean.orig( theta = par.DoL, x = c(1,x), sep.lam = sep.lam.global ) %*% covML.DoL %*% grad.condmean.orig( theta = par.DoL, x = c(1,x), sep.lam = sep.lam.global ) ) } )


plot( x.DoL, y.DoL, xlab = "Structure", ylab = "Contents", main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", original scale", sep = '' ) )
polygon(x = c( 1.5, 1, 1, -1, -1, 1.5 ), y =  c( -1, -1, 1, 1, 1.5, 1.5 ), density = 5, angle = -45, col = "lightgrey")
# if appropriate:
points( Structure$DoL[idx.outl], Contents$DoL[idx.outl], pch = 4 ) # ev. pch = 8
points( Structure$DoL[idx.remove], Contents$DoL[idx.remove], pch = 4 )
#
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ), col = 'red' ) # fitted median 
#lines( x.seq.DoL, BC.backtransform( opt.DoL$par[p+2], fitted.opt.seq.DoL ), col = 'lightpink2' )

# Backtransformed 95% Confidence interval (mle, lambda fixed --> 3)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ) ), col = 'orange', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ) ), col = 'orange', lty = 2 )

# Backtransformed 95% Prediction interval (mle, lambda fixed --> 3)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ) ), col = 'cyan1', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ) ), col = 'cyan1', lty = 2 )

# Backtransformed 95% Confidence interval (mle, lambda known --> 2)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), col = 'red', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), col = 'red', lty = 2 )

# Backtransformed 95% Prediction interval (mle, lambda known --> 2)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), col = 'blue', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), col = 'blue', lty = 2 )

# Backtransformed 95% Confidence interval (opt, var of back-transformed fitted value --> 5) 
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) - qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL), col = 'sienna1', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) + qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL), col = 'sienna1', lty = 2 )

# Backtransformed 95% Confidence interval (opt, lambda unknown --> 1)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ) ), col = 'lightpink2', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ) ), col = 'lightpink2', lty = 2 )  # use qnorm here???

# Backtransformed 95% Prediction interval (opt, lambda unknown --> 1)
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ) ), col = 'deepskyblue', lty = 2 )
lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ) ), col = 'deepskyblue', lty = 2 )

# Conditional mean on original scale
lines( x.seq.DoL, condmean.seq.DoL, col = 'green3', lwd = 1 )
lines( x.seq.DoL, condmean.seq.DoL - qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
lines( x.seq.DoL, condmean.seq.DoL + qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
box()

# Smearing estimate of the mean on original scale (Duan, 1983)
lines( x.seq.DoL, sapply( xlam.seq.DoL, function(x){ mean( BC.backtransform( mle.DoL$lambda.hat[1], sum( c(1,x)*mle.DoL$beta.hat ) + e.DoL ) ) } ), col = 'mediumpurple1' )
# Basically identical with condmean.seq
max( abs( condmean.seq.DoL - sapply( xlam.seq.DoL, function(x){ mean( BC.backtransform( mle.DoL$lambda.hat[1], sum( c(1,x)*mle.DoL$beta.hat ) + e.DoL ) ) } ) ), na.rm = TRUE )



## "Confidence interval" for the mean on original scale: Use 95% CI bounds as x instead of the fitted values (what Markus had done...)
lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), col = 'magenta', lty = 2 )
lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), col = 'magenta', lty = 2 )

# "Prediction interval" for the mean on original scale (obtained as before)
lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), col = 'skyblue', lty = 2 )
lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), col = 'skyblue', lty = 2 )


#####
#°°° TBS model (named "other") °°°#
#####
# loglik.other( c(0.1,1,0.3,0.25), y = y.DoL, x = x.DoL )
# opt.other.DoL <- optim( c(0.1,1,0.4,0.25), loglik.other, y = y.DoL, x = x.DoL, control = list( fnscale = -1 ), hessian = TRUE )
# sqrt( diag( solve( -opt.other.DoL$hessian ) ) )
# TBS.CI.lambda( y.DoL, x.DoL, beta.init = c(0.1,1) )

# y.lambda.DoL <- BC.transform( opt.other.DoL$par[p+2], y.DoL )
# yhat.lambda.DoL <- BC.transform( opt.other.DoL$par[p+2], c( cbind( 1, x.DoL ) %*% opt.other.DoL$par[1:p] ) )

plot( lam.seq.DoL, sapply( lam.seq.DoL, TBS.profllkhd.lambda, y = y.DoL, x = x.DoL, intercept = TRUE, beta.init = c(0.1,1) ) )
mle.other.DoL <- TBS.mle( y = y.DoL, x = x.DoL, beta.init = c(0.1,1) )
covML.other.DoL <- solve( -optimHess( unlist( mle.other.DoL[c(2,3,1)] ), loglik.other, y = y.DoL, x = x.DoL, control = list( fnscale = -1 ) ) )
TBS.CI.lambda( y.DoL, x.DoL, beta.init = c(0.1,1) )

y.lambda.DoL <- BC.transform( mle.other.DoL$lambda.hat, y.DoL )
yhat.lambda.DoL <- BC.transform( mle.other.DoL$lambda.hat, c( cbind( 1, x.DoL ) %*% mle.other.DoL$beta.hat ) )
e.DoL <- y.lambda.DoL - yhat.lambda.DoL

plot( qnorm( 1:n.obs/(n.obs+1) ), sort( e.DoL/sqrt(mle.other.DoL$sigma2.hat) ) )
abline( 0, 1, col = 'red' )

plot( yhat.lambda.DoL, e.DoL, xlab = "Fitted values", ylab = "Residuals", main = "Relative loss" ) # /median( abs(e.DoL) )
abline( h = 0, col = 'red' )
lines( lowess( x = yhat.lambda.DoL, y = e.DoL ), col = 'green' )

fitted.seq.DoL <- c( cbind( 1, x.seq.DoL ) %*% mle.other.DoL$beta.hat )

plot( x.DoL, y.DoL, xlab = "Structure", ylab = "Content", main = "Relative loss, original scale" )
points( Structure$DoL[idx.outl], Contents$DoL[idx.outl], pch = 4 )
abline( mle.other.DoL$beta.hat, col = 'red' )
lines( lowess( x = x.DoL, y = y.DoL ), col = 'green' )
# 95% CI of x^T beta.hat (which is the fitted median):
lines( x.seq.DoL, fitted.seq.DoL - qnorm(0.975) * sqrt( sapply( x.seq.DoL, function(x){ c( c(1,x) %*% covML.other.DoL[1:p,1:p] %*% c(1,x) ) } ) ), col = 'red', lty = 2 )
lines( x.seq.DoL, fitted.seq.DoL + qnorm(0.975) * sqrt( sapply( x.seq.DoL, function(x){ c( c(1,x) %*% covML.other.DoL[1:p,1:p] %*% c(1,x) ) } ) ), col = 'red', lty = 2 )
# Estimated quantiles (2.5% and 97.5%) --> these are an approximate 95% PI:
lines( x.seq.DoL, fitted.seq.DoL * ( 1 + qnorm(0.025) * mle.other.DoL$lambda.hat * sqrt(mle.other.DoL$sigma2.hat)/fitted.seq.DoL^mle.other.DoL$lambda.hat )^(1/mle.other.DoL$lambda.hat), col = 'blue', lty = 2 )
lines( x.seq.DoL, fitted.seq.DoL * ( 1 + qnorm(0.975) * mle.other.DoL$lambda.hat * sqrt(mle.other.DoL$sigma2.hat)/fitted.seq.DoL^mle.other.DoL$lambda.hat )^(1/mle.other.DoL$lambda.hat), col = 'blue', lty = 2 )

# Conditional mean:
condmean.other.seq.DoL <- TBS.condmean.orig( theta = unlist(mle.other.DoL[c(2,3,1)]), x.beta = fitted.seq.DoL )
var.condmean.other.DoL <- sapply( x.seq.DoL, function(x){ xbeta <- c( c(1,x) %*% mle.other.DoL$beta.hat ); grad.mean <- c( ( 1 + mle.other.DoL$sigma2.hat * (1 - mle.other.DoL$lambda.hat) * (1 - 2*mle.other.DoL$lambda.hat)/( 2*xbeta^(2*mle.other.DoL$lambda.hat) ) ) * c(1,x),  (1 - mle.other.DoL$lambda.hat)/( 2*xbeta^(2*mle.other.DoL$lambda.hat - 1) ), -mle.other.DoL$sigma2.hat/( 2*xbeta^(2*mle.other.DoL$lambda.hat - 1) ) * ( 1 + 2*(1 - mle.other.DoL$lambda.hat) * log(xbeta) ) ); c( grad.mean %*% covML.other.DoL %*% grad.mean ) } )

lines( x.seq.DoL, condmean.other.seq.DoL, col = 'green3' )
# 95% CI:
lines( x.seq.DoL, condmean.other.seq.DoL - qnorm(0.975) * sqrt( var.condmean.other.DoL ), col = 'green3', lty = 2 )
lines( x.seq.DoL, condmean.other.seq.DoL + qnorm(0.975) * sqrt( var.condmean.other.DoL ), col = 'green3', lty = 2 )

# Quantity that should be small for Taylor appx to be valid
range( mle.other.DoL$lambda.hat*sqrt(mle.other.DoL$sigma2.hat)/c(cbind(1,x.DoL)%*% mle.other.DoL$beta.hat)^mle.other.DoL$lambda.hat )  # not so optimal...
#°°° End TBS model °°°#
#######



########
# Alternative estimators for lambda
########

#### Sensitivity to removed data points
## lambda transforming to a symmetric distribution (p.135 Carroll and Ruppert 1984a)
plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), T.skew, x = x.DoL, y = y.DoL, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Relative loss" )
abline( h = 0, col = 'red' )
uniroot( f = T.skew, interval = c(-2,2), x = x.DoL, y = y.DoL, intercept = TRUE )
# robust: lambda_sk = 0.344, w/o outlier 0.339

plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), T.skew, x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Relative loss" )
abline( h = 0, col = 'red' )
uniroot( f = T.skew, interval = c(-2,2), x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE )
# less robust: lambda_sk = 0.179

## lambda transforming to a homoskedastic distribution (p.135 Carroll and Ruppert 1984a)
plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), H.hesk, x = x.DoL, y = y.DoL, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Relative loss" )
abline( h = 0, col = 'red' )
uniroot( f = H.hesk, interval = c(-1,1), x = x.DoL, y = y.DoL, intercept = TRUE )
# robust: lambda_hesk = 0.136, w/o outlier 0.130

### DOESN'T WORK because fitted values negative (thus can't take their log)
plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), H.hesk, x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE ), type = 'l', xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Relative loss" )
abline( h = 0, col = 'red' )
uniroot( f = H.hesk, interval = c(-1,1), x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE )
# less robust: lambda_hesk = 0.294 (???? bcs of log(fitted^2))

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


### Breusch-Paggan test (trials, not complete!!)
e.DoL^2/mle.DoL$sigma2.hat 
1- pf( ( sum( c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% e.DoL^2 )^2 ) - sum( (( diag(n.obs) - H.lambda.DoL ) %*% e.DoL^2)^2 ) )/( sum( (( diag(n.obs) - H.lambda.DoL ) %*% e.DoL^2)^2 )/(n.obs-1) ), df1 = 1, df2 = n.obs - 1 )


sum( c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% (e.DoL^2/mle.DoL$sigma2.hat)^2 )^2 )
sum( (( diag(n.obs) - H.lambda.DoL ) %*% (e.DoL^2/mle.DoL$sigma2.hat)^2)^2 )

/( sum( (( diag(n.obs) - H.lambda.DoL ) %*% (e.DoL^2/mle.DoL$sigma2.hat)^2)^2 )/(n.obs-1) ), df = 1 )

c( ( diag(n.obs) - matrix( 1/n.obs, n.obs, n.obs ) ) %*% e.DoL^2 ) - 
c( matrix( 1/n.obs, n.obs, n.obs ) %*% e.DoL^2 ) - mean( e.DoL^2 )

solve( crossprod( X.lambda.DoL ) ) %*% t(X.lambda.DoL) %*% e.DoL

1 - pf( ( sum( (y.lambda.DoL - mean(y.lambda.DoL))^2 ) - sum( e.DoL^2 ) )/( sum( e.DoL^2 )/(n.obs - p) ), df1 = 1, df2 = n.obs - p )

### Correlations between y.lambda and x.lambda
cor.test( y.lambda.DoL, x.lambda.DoL )
cor.test( y.lambda.DoL, x.lambda.DoL, method = "spearman" )
cor.test( y.lambda.DoL, x.lambda.DoL, method = "kendall" )

### Adjusted R^2:
1 - ( sum(e.DoL^2)/(n.obs-p) )/( sum( ( y.lambda.DoL - mean(y.lambda.DoL) )^2 )/(n.obs-1) )


