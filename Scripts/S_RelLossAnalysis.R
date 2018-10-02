############################################## Relative loss model #####################################################
# Overview:
  plot( Structure$DoL, Contents$DoL, xlab = "Structure", ylab = "Content", main = "Relative Loss" )
  points( Structure$DoL[c(145,249,328)], Contents$DoL[c(145,249,328)], pch = 4, col = c('red','red','blue') )
  abline( h = 1, lty = 2, col = 'grey' )
  abline( v = 1, lty = 2, col = 'grey' )
  abline( 0, 1, lty = 4, col = 'grey' )

########################################################################################################################
#==== Fit Box-Cox transformation =======================================================================================
########################################################################################################################

p <- 2 # (???)

# 1. Choose model type: ================================================================================================
  # Execute only one of the three possibilities:
    #==== 1: All data points ====#
      y.DoL <- Contents$DoL
      x.DoL <- Structure$DoL
      
      modelname <- "Full"
      
    #==== 3: W/o outliers nor lev. pts ====#
      idx.remove <- c(145,249) # based on complete data set
      
      y.DoL <- Contents$DoL[-idx.remove]
      x.DoL <- Structure$DoL[-idx.remove]
    
      modelname <- "notOutLev"
    
    #==== 2: W/o outliers ====#  <--- suggested & used in paper.
      idx.outl <- 145
      
      y.DoL <- Contents$DoL[-idx.outl]
      x.DoL <- Structure$DoL[-idx.outl]
      
      modelname <- "noOutl"
      
  #==== End ====#
  
  n.obs <- length(y.DoL)
  
  sep.lam.global <- F # Set TRUE if lambda paramaters should be estimated separately -> if T, PTBS seplam model
  if(sep.lam.global) { modelname <- paste( modelname, "_seplam", sep = "" ) } # Change modelname according to decision
  
  # Comparison of models:
    PTBS.mle( y = Contents$DoL, x = Structure$DoL, sep.lam = sep.lam.global) # Full
    PTBS.mle( y = Contents$DoL[-idx.outl], x = Structure$DoL[-idx.outl], sep.lam = sep.lam.global ) # rmv Outliers
    # remove Outliers & leverage points:
    PTBS.mle( y = Contents$DoL[-idx.remove], x = Structure$DoL[-idx.remove], sep.lam = sep.lam.global )
    # no big differences, just remove outliers!
  
  PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global)

#=======================================================================================================================
# 2. Model Fit: ========================================================================================================
  # Confidence intervals:
    CI.lambda.DoL <- PTBS.CI.lambda( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global )
    CI.lambda.DoL
    #==# Full: CI = (0.1352,0.2579)
    #==# W/o outliers: CI = (0.144,0.265), lambda also slightly shifted to the right
    #==# W/o extr.Cdist: CI = (0.151,0.272) even more shifted to the right (as lambda itself)

  # fitted parameters:
    mle.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = sep.lam.global )
    # fit seems better with sep.lam = T
    par.DoL <- unlist(mle.DoL[c(2,3,1)])
    s2.hat.DoL <- mle.DoL$sigma2.hat * n.obs/(n.obs-p)
  # Covariance:
    covML.DoL <- solve( -optimHess( par.DoL, PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = sep.lam.global, 
                                    control = list( maxit = 5000, fnscale = -1 ) ) )
    
  # Estimation results
    cbind( par.DoL, sqrt(diag(covML.DoL)) )
    mle.DoL$ell.opt
    2*( -mle.DoL$ell.opt + p + 2 + sep.lam.global )

  # p-value LRT seplam vs. single lambda (in case of sep.lam.global = TRUE)
    # 1 - pchisq( 2*( mle.DoL$ell.opt - PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = F )$ell.opt ), df = 1 )


  # Transformed quantities
    y.lambda.DoL <- BC.transform( mle.DoL$lambda.hat[1], y.DoL )
    x.lambda.DoL <- BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], x.DoL )
    X.lambda.DoL <- cbind( 1, x.lambda.DoL )
    
  # Covariance matrix of beta based on regression (#3 below):
    covbeta.regr.DoL <- s2.hat.DoL * solve(crossprod(X.lambda.DoL))
    
    theta.init.DoL <- if(sep.lam.global) { c(0, 1, 0.4, 0.3, 0.5) } else {  c( 0, 1, 0.6, 0.5 ) }
    opt.DoL <- optim( theta.init.DoL, fn = PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = sep.lam.global, 
                      control = list( maxit = 5000, fnscale = -1 ), hessian = TRUE )
    opt.DoL$val - mle.DoL$ell.opt
    covML.DoL <- solve( -opt.DoL$hessian )
  
  ## !! sigma2.hat is now correlated because of the estimated lambda:
    cov2cor( covML.DoL )

#=======================================================================================================================
# 3. Plots: ============================================================================================================

##=======================================##
##== Plot confidence region for lambda ==##
##=======================================##
if(!sep.lam.global){ 
    ##===== same lambda =====##
    x11()
    proflambda.DoL <- sapply( lam.seq.DoL, PTBS.profllkhd.lambda, y = y.DoL, x = x.DoL, intercept = TRUE )
    
    ## Plot profile log-likelihood for lambda
    plot( lam.seq.DoL, proflambda.DoL - CI.lambda.DoL[[4]], xlab = expression(lambda), ylab = "Profile log likelihood", 
          type = 'l', xlim = c(-1,1), ylim = c(-200,0) )
    abline( h = - qchisq(0.95,1)/2, lty = 2 )
    abline( v = CI.lambda.DoL[1:3], lty = 3 )
    savePlot( paste( "Figures/DoL_", modelname, "_proflambda", sep = '' ), type = "pdf" )
    ##===== End =====##
    
} else { 
    ##===== separate lambdas =====##
    # Profile likelihood array
    x11()
    lam1.seq.DoL <- seq( 0.001, 0.35, 0.001 )
    lam2.seq.DoL <- seq( 0.08, 0.5, 0.001 )
    lam.arr.DoL <- expand.grid( lam1 = lam1.seq.DoL, lam2 = lam2.seq.DoL )
    
    proflambda.DoL <- matrix( apply( lam.arr.DoL, 1, PTBS.profllkhd.lambda, y = y.DoL, x = x.DoL ), 
                              nr = length(lam1.seq.DoL) )
    
    # Plot confidence region for (lambda.y, lambda.x)
    image( x = lam1.seq.DoL, y = lam2.seq.DoL, pmax( proflambda.DoL - (mle.DoL$ell.opt - qchisq(0.95, df=2)/2), 0 ),
           xlab = expression(lambda[y]), ylab = expression(lambda[x]), breaks = c(-0.0001,0.0001,seq( 0.05, 3, 0.05 )), 
           col = c('white',topo.colors(60)), main = "Relative loss, 95% CR lambdas" )
    abline( v = mle.DoL$lambda.hat[1], col = 'red' )
    abline( h = mle.DoL$lambda.hat[2], col = 'red' )
    abline( v = unlist(CI.lambda.DoL)[c(3,5)], col = 'orange', lty = 2 )
    abline( h = unlist(CI.lambda.DoL)[c(4,6)], col = 'orange', lty = 2 )
    abline( 0, 1, col = 'darkgrey', lty = 2 )
    box()
    savePlot( paste( "Figures/DoL_", modelname, "_ConfReg", sep = '' ), type = "pdf" )
    # identity line NOT through CR --> lambdas are different
    ##===== End =====##
}

##===================##
##=== Diagnostics ===##
##===================##

x.seq.DoL <- BC.backtransform( 0.2, seq( -5.75, 0.3, 0.025 ) )

# Involved quantities
  H.lambda.DoL <- X.lambda.DoL %*% solve( crossprod( X.lambda.DoL ) ) %*% t(X.lambda.DoL)
  # sum( diag( H.lambda.DoL ) ) # tr(H) = p = 2
# Fitted values:
  yhat.lambda.DoL <- c( H.lambda.DoL %*% y.lambda.DoL )
# Residuals:
  e.DoL <- y.lambda.DoL - yhat.lambda.DoL
# Standardised residuals:
  r.DoL <- e.DoL/sqrt(s2.hat.DoL*(1-diag(H.lambda.DoL)))
# Deletion residuals (Davison 2008, p.395) = "studentized residuals" in Weisberg 2005 (p.196):
  del.resid.DoL <- sqrt( (n.obs-p-1)/(n.obs-p-r.DoL^2) )*r.DoL

# Outliers
  if(modelname=="Full"){
    idx.tail <- rev( tail( order( abs(del.resid.DoL) ) ) )
    0.05/n.obs  # Bonferroni correction. Smaller p-values from t-test indicate outliers:
    2*( 1 - pt( abs(del.resid.DoL[idx.tail[1:3]]), df = n.obs - p - 1 ) )
    #==# Full: The most extreme one (145) is an outlier, the following not anymore.
    #==# W/o outlier: No outliers anymore.
    idx.outl <- idx.tail[1]
    min(abs(del.resid.DoL))
  }
  
# large residuals
  lrg.resid <- which( abs(r.DoL) > 2 )
  length( lrg.resid )
  #==# Full: 15
  #==# W/o outlier: 16, sep.lam 18

# High leverage (Davison 2008, p.394 top)
  idx.highlever <- which( diag(H.lambda.DoL) > 2*p/n.obs )
  length( idx.highlever ) # quite a lot
  #==# Full: 34
  #==# W/o outlier: 33, sep.lam 34

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
  #==# W/o outlier: lrg.resid and highlever disjoint, extr.Cdist both has a large residual but no high leverage, 
  #==#  one each with sep.lam


  ##===================##
  ##=== Diag. Plots ===##
  ##===================##
  

# Residuals
  x11( width = 10.7, height = 10.5 )
  par(mfrow = c(2,2))
# Normality of standardised residuals
  qqnorm( r.DoL )
  abline( 0, 1, col = 'red' )

# Std residuals vs fitted values
  plot( yhat.lambda.DoL, r.DoL, xlab = "Fitted values", ylab = "Standardised residuals", 
        main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  abline( h = 0, col = 'red' )
  abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
  points( yhat.lambda.DoL[extr.Cdist], r.DoL[extr.Cdist], col = 'red', pch = 19 )
  points( yhat.lambda.DoL[idx.highlever], r.DoL[idx.highlever], col = 'blue' )
  #==# quite obvious pattern
  #==# W/o outl nor lev pts: Extreme Cook's distance has also high leverage
  #==# High leverage not very interesting since just largest and smallest values of x...

# Absolute std residuals vs. fitted values (see heteroscedasticity)
  plot( yhat.lambda.DoL, abs(r.DoL), xlab = "Fitted values", ylab = "", 
        main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
  points( yhat.lambda.DoL[extr.Cdist], abs(r.DoL)[extr.Cdist], col = 'red', pch = 19 )
  points( yhat.lambda.DoL[idx.highlever], abs(r.DoL)[idx.highlever], col = 'blue' )
  abline( h = 2, col = 'sandybrown', lty = 2 )
  
# Sqrt of absolute std residuals vs. fittet values
  plot( yhat.lambda.DoL, sqrt(abs(r.DoL)), xlab = "Fitted values", ylab = "", 
        main = paste( "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
  points( yhat.lambda.DoL[extr.Cdist], sqrt(abs(r.DoL[extr.Cdist])), col = 'red', pch = 19 )
  points( yhat.lambda.DoL[idx.highlever], sqrt(abs(r.DoL[idx.highlever])), col = 'blue' )
  abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
  savePlot( paste( "Figures/DoL_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )
  #==# Full: Fit maybe not that good, but apart from outlier quite homoskedastic (smaller variance for large x-values).
  #==# Some problem of fit, but apart from outlier quite homoskedastic. Smaller variance for high leverage points
  #==#  with large x!

# Running variance along fitted values for block sizes 20, 26, 32
  x11( width = 9.8, height = 6.7 )
  plot( sort(yhat.lambda.DoL)[9+seq(1,365,4)], sapply( seq(1,365,4), 
    function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:19], na.rm = TRUE ) } ),
    type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = paste(
      "Relative loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' )) # block size 20
  lines( sort(yhat.lambda.DoL)[12+seq(1,357,4)], c( sapply( 
    seq(1,353,4), function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:25], na.rm = TRUE ) } ), 
    var( r.DoL[order(yhat.lambda.DoL)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
  lines( sort(yhat.lambda.DoL)[15+seq(1,353,4)], sapply(seq(1,353,4), 
              function(i){ var( r.DoL[order(yhat.lambda.DoL)][i+0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
  abline( h = 1, col = 'darkgrey', lty = 2 )
  savePlot( paste( "Figures/DoL_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

# Compare residual plots for the two models (ONLY if separate lambda)
  # x11( width = 6.6, height = 7 )
  # plot( yhat.lambda.DoL, r.DoL, xlab = "Fitted values", ylab = "Standardised residuals", main = "Relative loss, 
  #       separate lambda" )
  # abline( h = 0, col = 'red' )
  # savePlot( paste( "Figures/DoL_", modelname, "_resid_compare", sep = '' ), type = "pdf" )

# Leverage vs. Cook's distance vs. high residuals -> master pieces of plots
  x11( width = 13.2, height = 7 )
  par(mfrow = c(1,2))
  # Distinguish outliers from leverage points (Davison 2008, p.395)
  plot( diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)), Cook.DoL, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", 
        main = "Relative loss, transformed scale" )
  points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[extr.Cdist], Cook.DoL[extr.Cdist], col = 'red', pch = 19 )
  points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[idx.highlever], Cook.DoL[idx.highlever], 
          col = 'blue' ) # quite obvious pattern
  abline( h = 8/(n.obs - 2*p), col = 'red' )
  points( (diag(H.lambda.DoL)/(1-diag(H.lambda.DoL)))[lrg.resid], Cook.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
  abline( 0, 4/p, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
  abline( v = 2*p/(n.obs - 2*p), col = 'blue' )
  
# Residuals vs. leverage (like plot.lm but with better level curves)
  plot( diag(H.lambda.DoL), r.DoL, xlab = "Leverage", ylab = "Standardised residuals", 
        main = "Relative loss, transformed scale", xlim = c(0,max(diag(H.lambda.DoL))) ) # xlim = c(0,0.025)
  abline( h = 0, col = 'grey' )
  abline( v = 2*p/n.obs, col = 'blue' )
  abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
  points( diag(H.lambda.DoL)[extr.Cdist], r.DoL[extr.Cdist], col = 'red', pch = 19 )
  points( diag(H.lambda.DoL)[idx.highlever], r.DoL[idx.highlever], col = 'blue' ) # quite obvious pattern
  points( diag(H.lambda.DoL)[lrg.resid], r.DoL[lrg.resid], pch = 3, col = 'sandybrown' )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p)) * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
  savePlot( paste( "Figures/DoL_", modelname, "_leverage", sep = "" ), type = "pdf" )

   # We also tried a model with two piecewise linear functions (1 knot) -> not good idea
  
##======================================##
##=== CI and PI on transformed scale ===##
##======================================##

# There are three possible scenarios under which standard errors of prediction x^T*beta.hat can be computed:
# 1. Delta method on full covariance matrix: Account for estimation of lambda. This is not proper because uncertainty of
#    lambda in y.lambda is not accounted for.
# 2. Use std errors for beta from full covariance matrix (thus estimation of lambda is accounted for in se.beta), but 
#    then ignore that lambda was estimated when appyling the delta method for var(x0^T beta.hat).
# 3. Regression-based while treating lambda = lambda.hat as fixed. Thus the uncertainty of estimating lambda is 
#    contained nowhere, not even in the std errors of beta.
# 4. As 2 but bootstrap-based (not good)
# 5. Detour via original scale by computing confidence bands for y.hat (p-quantile) (accounting for uncertainty in 
#    lambda), which are then transformed again.

# --> I have two (very similar) sets of estimated parameters (mle and opt). Based on my calculations, I could also code 
#     the complete observed information matrix for mle, but currently I haven't done it. Therefore I have to used the 
#     opt-estimates for full uncertainty.
# --> By applying optimHess to mle.DoL I get the Hessian and thus the joint covariance at mle.DoL. Thus I do not need 
#     the opt version anymore. Since 2 seems better than 3, I keep 2 and 5.

fitted.mle.seq.DoL <- c( cbind( 1, xlam.seq.DoL ) %*% mle.DoL$beta.hat )

# (2) Assuming lambda = lambda.hat known, based on mle
  var.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ c(1,x) %*% covML.DoL[1:p,1:p] %*% c(1,x) } )
  var.predict.seq.DoL <- var.fitted.seq.DoL + s2.hat.DoL

# (3) Cov matrix of beta.hat from the regression of y.lambda on X.lambda with lambda = lambda.hat fixed
  varRegr.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ c(1,x) %*% covbeta.regr.DoL %*% c(1,x) } )

# (5) Detour via original scale to account for estimated lambda in y as well (but not sure it is properly done)
  varOrig.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(xlam){ 
    grad.yhat( theta = par.DoL, x.lam = c(1,xlam), sep.lam = sep.lam.global ) %*% 
      covML.DoL %*% grad.yhat( theta = par.DoL, x.lam = c(1,xlam), sep.lam = sep.lam.global ) } )

# (1) Including uncertainty on lambda.hat in x.lambda, based on opt
  # varFull.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ grad.fitted( par.DoL, c(1,x) ) %*% covML.DoL %*% 
  #     grad.fitted( par.DoL, c(1,x) ) } ) # sigma2 has no influence on this
  # varFull.predict.seq.DoL <- varFull.fitted.seq.DoL + s2.hat.DoL

# (4) As 1 but with cov(beta) based on non-parametric bootstrap (not good)
  # varBoot.fitted.seq.DoL <- sapply( xlam.seq.DoL, function(x){ 
  #   grad.fitted( par.DoL, c(1,x) ) %*% covBoot.DoL %*% grad.fitted( par.DoL, c(1,x) ) } )
  # varBoot.predict.seq.DoL <- varBoot.fitted.seq.DoL + s2.hat.DoL

  # x11()
  # plot( xlam.seq.DoL, varOrig.fitted.seq.DoL, type = 'l', xlab = "x(lambda)", ylab = "Variance of X(lambda)*beta.hat", 
  #       main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
  # # lines( xlam.seq.DoL, varBoot.fitted.seq.DoL, col = 'dodgerblue' )
  # lines( xlam.seq.DoL, var.fitted.seq.DoL, col = 'blue' )
  # lines( xlam.seq.DoL, varRegr.fitted.seq.DoL, col = 'cyan1' )
  # dev.off()

  x11()
# Transformed scale: y.lambda vs. X.lambda*beta
  plot( x.lambda.DoL, y.lambda.DoL, xlab = "Structure", ylab = "Contents", 
        main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
  polygon( x = BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], c( 1.5, 1, 1, 0, 0, 1.5 ) ), 
           y = BC.transform( mle.DoL$lambda.hat[1], c( 0, 0, 1, 1, 1.5, 1.5 ) ), density = 5, angle = -45, 
           col = "lightgrey")
  # if appropriate:
  points( BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], Structure$DoL[idx.outl] ),
          BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[idx.outl] ), pch = 4 ) # ev. pch = 8
  # points( BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], Structure$DoL[idx.remove] ),
  #         BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[idx.remove] ), pch = 4 )
  #
  abline( mle.DoL$beta.hat, col = 'red' )

# Method 3:
  # 95% Confidence interval (mle, lambda fixed --> 3)
  lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ), col = 'orange', 
         lty = 2 )
  lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.DoL ), col = 'orange', 
         lty = 2 )
  # 95% Prediction interval (mle, lambda fixed --> 3)
  lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ), 
         col = 'cyan1', lty = 2 )
  lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ), 
         col = 'cyan1', lty = 2 )

# Method 2:
  # 95% Confidence interval (mle, lambda known --> 2)
  lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ), col = 'red', lty = 2 )
  lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ), col = 'red', lty = 2 )
  # 95% Prediction interval (mle, lambda known --> 2)
  lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ), col = 'blue', lty = 2 )
  lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ), col = 'blue', lty = 2 )

# Method 5:
  # 95% Confidence interval (opt, via original scale --> 5)
  lines( xlam.seq.DoL, BC.transform( mle.DoL$lambda.hat[1], BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) - 
                                       qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL) ), col = 'sienna1', lty = 2 )
  lines( xlam.seq.DoL, BC.transform( mle.DoL$lambda.hat[1], BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) + 
                                       qnorm(0.975) * sqrt(varOrig.fitted.seq.DoL) ), col = 'sienna1', lty = 2 )
  # Prediction interval accounting for estimation uncertainty is too complex for this case -> neglect
  lines( lowess( x = x.lambda.DoL, y = y.lambda.DoL ), col = 'green' )
  box()

# Method 1:
  # # 95% Confidence interval (opt, lambda unknown --> 1)
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ), 
  #        col = 'lightpink2', lty = 2 )
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.DoL ), 
  #        col = 'lightpink2', lty = 2 )  # use qnorm here???
  # # 95% Prediction interval (opt, lambda unknown --> 1)
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ), col = 'deepskyblue', lty = 2 )
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.DoL ), col = 'deepskyblue', lty = 2 )

# Method 4 (neglected):
  # # 95% Confidence interval (mle, lambda unknown --> 4)
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.DoL ), 
  #        col = 'chocolate3', lty = 2 )
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.DoL ), 
  #        col = 'chocolate3', lty = 2 )  # use qnorm here???
  # # 95% Prediction interval (mle, lambda unknown --> 4)
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.DoL ), 
  #        col = 'dodgerblue', lty = 2 )
  # lines( xlam.seq.DoL, fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.DoL ), 
  #        col = 'dodgerblue', lty = 2 )
  # # --> Much less difference for prediction intervals because sigma2 dominates the variance.


# Either 2 or 3 seem reasonable since very close. 1 and 4 are not plausible because uncertainty of lambda is accounted 
# for x, but not for y. And the detour via the original scale is very close to 2 or 3.


##====================================================================##
##=== Backtransformed to original scale (!! only for full data !!) ===##
##====================================================================##
  x.seq.DoL <- BC.backtransform( mle.DoL$lambda.hat[1 + sep.lam.global], xlam.seq.DoL )
  
# Conditional mean on original scale, estimator by Taylor 1986
  condmean.seq.DoL <- condmean.orig( theta = par.DoL, x = cbind(1,xlam.seq.DoL), sep.lam = sep.lam.global )
  # use sigma2 or s2 here???
  var.condmean.DoL <- sapply( xlam.seq.DoL, function(x){
    c( grad.condmean.orig( theta = par.DoL, x = c(1,x), sep.lam = sep.lam.global ) %*% covML.DoL %*% 
         grad.condmean.orig( theta = par.DoL, x = c(1,x), sep.lam = sep.lam.global ) ) } )


  plot( x.DoL, y.DoL, xlab = "Structure", ylab = "Contents", 
        main = paste( "Relative loss", ifelse( sep.lam.global, " seplam", "" ), ", original scale", sep = '' ) )
  polygon(x = c( 1.5, 1, 1, -1, -1, 1.5 ), y =  c( -1, -1, 1, 1, 1.5, 1.5 ), density = 5, angle = -45, 
          col = "lightgrey")
  
# if appropriate:
  points( Structure$DoL[idx.outl], Contents$DoL[idx.outl], pch = 4 ) # ev. pch = 8
  # points( Structure$DoL[idx.remove], Contents$DoL[idx.remove], pch = 4 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ), col = 'red' ) # fitted median 
  # lines( x.seq.DoL, BC.backtransform( opt.DoL$par[p+2], fitted.opt.seq.DoL ), col = 'lightpink2' )

# Method 3:
  # Backtransformed 95% Confidence interval (mle, lambda fixed --> 3)
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                        sqrt( varRegr.fitted.seq.DoL ) ), col = 'orange', lty = 2 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                        sqrt( varRegr.fitted.seq.DoL ) ), col = 'orange', lty = 2 )
  # Backtransformed 95% Prediction interval (mle, lambda fixed --> 3)
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                        sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ) ), col = 'cyan1', lty = 2 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                        sqrt( s2.hat.DoL + varRegr.fitted.seq.DoL ) ), col = 'cyan1', lty = 2 )

# Method 2:   
  # Backtransformed 95% Confidence interval (mle, lambda known --> 2)
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.fitted.seq.DoL ) ), col = 'red', lty = 2 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.fitted.seq.DoL ) ), col = 'red', lty = 2 )
# Backtransformed 95% Prediction interval (mle, lambda known --> 2)
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) ), col = 'blue', lty = 2 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) ), col = 'blue', lty = 2 )

# Method 5:  
  # Backtransformed 95% Confidence interval (opt, var of back-transformed fitted value --> 5) 
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) - qnorm(0.975) * 
           sqrt(varOrig.fitted.seq.DoL), col = 'sienna1', lty = 2 )
  lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ) + qnorm(0.975) * 
           sqrt(varOrig.fitted.seq.DoL), col = 'sienna1', lty = 2 )

# # Method 1:
#   # Backtransformed 95% Confidence interval (opt, lambda unknown --> 1)
#   lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p ) *
#                                         sqrt( varFull.fitted.seq.DoL ) ), col = 'lightpink2', lty = 2 )
#   lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p ) *
#                                         sqrt( varFull.fitted.seq.DoL ) ), col = 'lightpink2', lty = 2 )  # use qnorm here???
# # Backtransformed 95% Prediction interval (opt, lambda unknown --> 1)
#   lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p ) * 
#                                         sqrt( varFull.predict.seq.DoL ) ), col = 'deepskyblue', lty = 2 )
#   lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p ) * 
#                                         sqrt( varFull.predict.seq.DoL ) ), col = 'deepskyblue', lty = 2 )

  
##==========================================##
##=== Conditional mean on original scale ===##
##==========================================##
  
  lines( x.seq.DoL, condmean.seq.DoL, col = 'green3', lwd = 1 )
  lines( x.seq.DoL, condmean.seq.DoL - qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
  lines( x.seq.DoL, condmean.seq.DoL + qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
  box()

  #<------------------------------------------------- Experimental ------------------------------------------------------> 
# # Smearing estimate of the mean on original scale (Duan, 1983)
#   lines( x.seq.DoL, sapply( xlam.seq.DoL, function(x){ 
#     mean( BC.backtransform( mle.DoL$lambda.hat[1], sum( c(1,x)*mle.DoL$beta.hat ) + e.DoL ) ) } ), 
#     col = 'mediumpurple1' )
#   # Basically identical with condmean.seq
#   max( abs( condmean.seq.DoL - sapply( xlam.seq.DoL, function(x){ 
#     mean( BC.backtransform( mle.DoL$lambda.hat[1], sum( c(1,x)*mle.DoL$beta.hat ) + e.DoL ) ) } ) ), na.rm = TRUE )
# 
#   
# # Apply to confidence interval 
#   # "Confidence interval" for the mean on orig. scale: Use 95% CI bounds as x instead of the fitted values (1st try)
#   lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL - 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), 
#          col = 'magenta', lty = 2 )
#   lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL + 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL ) ), 
#          col = 'magenta', lty = 2 )
#   # "Prediction interval" for the mean on original scale (obtained as before)
#   lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL - 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), 
#          col = 'skyblue', lty = 2 )
#   lines( x.seq.DoL, condmean.orig( theta = unlist(mle.DoL)[c(2:4,1)], x.beta = fitted.mle.seq.DoL + 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL ) ), 
#          col = 'skyblue', lty = 2 )
#<------------------------------------------------- Experimental ------------------------------------------------------> 

##==========================================##
##=== Alternative estimators for lambda ====##
##==========================================##
  
# Sensitivity to removed data points
  # lambda transforming to a symmetric distribution (p.135 Carroll and Ruppert 1984a)
  plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), T.skew, x = x.DoL, y = y.DoL, intercept = TRUE ), type = 'l', 
        xlab = "lambda", ylab = "Skewness measure", main = "Relative loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = T.skew, interval = c(-2,2), x = x.DoL, y = y.DoL, intercept = TRUE )
  # robust: lambda_sk = 0.344, w/o outlier 0.339

  plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), T.skew, x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE ), 
        type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Relative loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = T.skew, interval = c(-2,2), x = x.DoL, y = y.DoL, intercept = TRUE, robust = FALSE )
  # less robust: lambda_sk = 0.179

# lambda transforming to a homoskedastic distribution (p.135 Carroll and Ruppert 1984a)
  plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), H.hesk, x = x.DoL, y = y.DoL, intercept = TRUE ), type = 'l', 
        xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Relative loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = H.hesk, interval = c(-1,1), x = x.DoL, y = y.DoL, intercept = TRUE )
# robust: lambda_hesk = 0.136, w/o outlier 0.130


# Adjusted R^2:
  adj_R2.DoL <- 1 - ( sum(e.DoL^2)/(n.obs-p) )/( sum( ( y.lambda.DoL - mean(y.lambda.DoL) )^2 )/(n.obs-1) )


