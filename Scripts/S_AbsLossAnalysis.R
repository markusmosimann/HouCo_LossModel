############################################## Absolute loss model #####################################################
# Overview:
  plot( Structure$Loss, Contents$Loss, xlab = "Structure", ylab = "Content", main = "Absolute Loss" )
  points( Structure$Loss[c(145,249,328)], Contents$Loss[c(145,249,328)], pch = 4, col = c('red','red','blue') )
  abline( 0, 1, lty = 4, col = 'grey' )

########################################################################################################################
#==== Fit Box-Cox transformation ======================================================================================#
########################################################################################################################

p <- 2 # (???)

# 1. Choose model type: ================================================================================================
  # Execute only one of the three possibilities:
    #==== 1: All data points ====#
      y.Loss <- Contents$Loss
      x.Loss <- Structure$Loss
      
      modelname <- "Full"
      
    #==== 3: W/o outliers nor lev. pts ====#
      idx.remove <- c(145,249) # based on complete data set
      
      y.Loss <- Contents$Loss[-idx.remove]
      x.Loss <- Structure$Loss[-idx.remove]
    
      modelname <- "notOutLev"
    
    #==== 2: W/o outliers ====#  <--- suggested & used in paper.
      idx.outl <- 145
      
      y.Loss <- Contents$Loss[-idx.outl]
      x.Loss <- Structure$Loss[-idx.outl]
      
      modelname <- "noOutl"
    
  #==== End ====#
  
  n.obs <- length(y.Loss)
  
  sep.lam.global <- F # Set TRUE if lambda paramaters should be estimated separately -> if T, PTBS seplam model
  if(sep.lam.global) { modelname <- paste( modelname, "_seplam", sep = "" ) } # Change modelname according to decision
  
  # Comparison of models:
    PTBS.mle( y = Contents$Loss, x = Structure$Loss, sep.lam = sep.lam.global) # Full
    PTBS.mle( y = Contents$Loss[-idx.outl], x = Structure$Loss[-idx.outl], sep.lam = sep.lam.global ) # rmv Outliers
    # remove Outliers & leverage points:
    PTBS.mle( y = Contents$Loss[-idx.remove], x = Structure$Loss[-idx.remove], sep.lam = sep.lam.global )
    # no big differences, just remove outliers!
  
  PTBS.mle( y = y.Loss, x = x.Loss, sep.lam = sep.lam.global)

#=======================================================================================================================
# 2. Model Fit: ========================================================================================================
  # Confidence intervals:
    CI.lambda.Loss <- PTBS.CI.lambda( y = y.Loss, x = x.Loss, sep.lam = sep.lam.global )
    CI.lambda.Loss
    #==# Full: CI = (0.0556,0.1814)
    #==# W/o outliers: CI = (0.0675, 0.1930)
    #==# W/o extr.Cdist: CI = ()

  # fitted parameters:
    mle.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, sep.lam = sep.lam.global )
    par.Loss <- unlist(mle.Loss[c(2,3,1)])
    s2.hat.Loss <- mle.Loss$sigma2.hat * n.obs/(n.obs-p)
  # Covariance:
    covML.Loss <- solve( -optimHess( par.Loss, PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = sep.lam.global, 
                                    control = list( maxit = 5000, fnscale = -1 ) ) )
    
  # Estimation results
    cbind( par.Loss, sqrt(diag(covML.Loss)) )
    mle.Loss$ell.opt
    2*( -mle.Loss$ell.opt + p + 2 + sep.lam.global )

  # p-value LRT seplam vs. single lambda (in case of sep.lam.global = TRUE)
    # 1 - pchisq( 2*( mle.Loss$ell.opt - PTBS.mle( y = y.Loss, x = x.Loss, sep.lam = F )$ell.opt ), df = 1 )


  # Transformed quantities
    y.lambda.Loss <- BC.transform( mle.Loss$lambda.hat[1], y.Loss )
    x.lambda.Loss <- BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], x.Loss )
    X.lambda.Loss <- cbind( 1, x.lambda.Loss )
    
  # Covariance matrix of beta based on regression (#3 below):
    covbeta.regr.Loss <- s2.hat.Loss * solve(crossprod(X.lambda.Loss))
    
    theta.init.Loss <- if(sep.lam.global) { c( 8, 0.05, 3.1, 0.2, 0.38 ) } else {c( 2, 1, 4, 0.1 ) + c(384-n.obs,0,0,0)}
    opt.Loss <- optim( c( theta.init.Loss, if(sep.lam.global){ 0.1 } else { NULL } ), fn = PTBS.llkhd, y = y.Loss, 
                       x = x.Loss, sep.lam = sep.lam.global, control = list( maxit = 5000, fnscale = -1 ), 
                       hessian = TRUE )
    opt.Loss$val - mle.Loss$ell.opt
    covML.Loss <- solve( -opt.Loss$hessian )
  
  ## !! sigma2.hat is now correlated because of the estimated lambda:
    cov2cor( covML.Loss )

#=======================================================================================================================
# 3. Plots: ============================================================================================================

##=======================================##
##== Plot confidence region for lambda ==##
##=======================================##
if(!sep.lam.global){ 
    ##===== same lambda =====##
    x11()
    proflambda.Loss <- sapply( lam.seq.Loss, PTBS.profllkhd.lambda, y = y.Loss, x = x.Loss, intercept = TRUE )
    
    ## Plot profile log-likelihood for lambda
    plot( lam.seq.Loss, proflambda.Loss - CI.lambda.Loss[[4]], xlab = expression(lambda), ylab = "Profile log likelihood", 
          type = 'l', xlim = c(-1,1), ylim = c(-200,0) )
    abline( h = - qchisq(0.95,1)/2, lty = 2 )
    abline( v = CI.lambda.Loss[1:3], lty = 3)
    savePlot( paste( "Figures/Loss_", modelname, "_proflambda", sep = '' ), type = "pdf" )
    ##===== End =====##
    
} else { 
    ##===== separate lambdas =====##
    # Profile likelihood array
    x11()
    lam1.seq.Loss <- seq( -0.07, 0.25, 0.001 )
    lam2.seq.Loss <- seq( 0.1, 0.5, 0.001 )
    lam.arr.Loss <- expand.grid( lam1 = lam1.seq.Loss, lam2 = lam2.seq.Loss )
    
    proflambda.Loss <- matrix( apply( lam.arr.Loss, 1, PTBS.profllkhd.lambda, y = y.Loss, x = x.Loss ), 
                              nr = length(lam1.seq.Loss) )
    
    # Plot confidence region for (lambda.y, lambda.x)
    image( x = lam1.seq.Loss, y = lam2.seq.Loss, pmax( proflambda.Loss - (mle.Loss$ell.opt - qchisq(0.95, df=2)/2), 0 ),
           xlab = expression(lambda[y]), ylab = expression(lambda[x]), breaks = c(-0.0001,0.0001,seq( 0.05, 3, 0.05 )), 
           col = c('white',topo.colors(60)), main = "Absolute loss, 95% CR lambdas" )
    abline( v = mle.Loss$lambda.hat[1], col = 'red' )
    abline( h = mle.Loss$lambda.hat[2], col = 'red' )
    abline( v = unlist(CI.lambda.Loss)[c(3,5)], col = 'orange', lty = 2 )
    abline( h = unlist(CI.lambda.Loss)[c(4,6)], col = 'orange', lty = 2 )
    abline( 0, 1, col = 'darkgrey', lty = 2 )
    box()
    savePlot( paste( "Figures/Loss_", modelname, "_ConfReg", sep = '' ), type = "pdf" )
    # identity line NOT through CR --> lambdas are different
    ##===== End =====##
}

##===================##
##=== Diagnostics ===##
##===================##

# Involved quantities
  H.lambda.Loss <- X.lambda.Loss %*% solve( crossprod( X.lambda.Loss ) ) %*% t(X.lambda.Loss)
  # sum( diag( H.lambda.Loss ) ) # tr(H) = p = 2
# Fitted values:
  yhat.lambda.Loss <- c( H.lambda.Loss %*% y.lambda.Loss )
# Residuals:
  e.Loss <- y.lambda.Loss - yhat.lambda.Loss
# Standardised residuals:
  r.Loss <- e.Loss/sqrt(s2.hat.Loss*(1-diag(H.lambda.Loss)))
# Deletion residuals (Davison 2008, p.395) = "studentized residuals" in Weisberg 2005 (p.196):
  del.resid.Loss <- sqrt( (n.obs-p-1)/(n.obs-p-r.Loss^2) )*r.Loss

# Outliers
  if(modelname=="Full"){
    idx.tail <- rev( tail( order( abs(del.resid.Loss) ) ) )
    0.05/n.obs  # Bonferroni correction. Smaller p-values from t-test indicate outliers:
    2*( 1 - pt( abs(del.resid.Loss[idx.tail[1:3]]), df = n.obs - p - 1 ) )
    #==# Full: No outlier?! (also with two lambdas)
    #==# W/o outlier: No outliers anymore.
    idx.outl <- idx.tail[1]
    min(abs(del.resid.Loss))
  }
  
# large residuals
  lrg.resid <- which( abs(r.Loss) > 2 )
  length( lrg.resid )
  #==# W/o outlier: 19, sep.lam 18

# High leverage (Davison 2008, p.394 top)
  idx.highlever <- which( diag(H.lambda.Loss) > 2*p/n.obs )
  length( idx.highlever ) # quite a lot
  #==# Full: 26 (24 with two lambdas)
  #==# W/o outlier: 24
  #==# W/o outl nor lev pts: 34
  #==# W/o outl nor lev pts and lambda = 0.25: 33 (24 with two lambdas)

  Cook.Loss <- (r.Loss^2)*diag(H.lambda.Loss)/(p*(1 - diag(H.lambda.Loss)))
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


  ##===================##
  ##=== Diag. Plots ===##
  ##===================##
  

# Residuals
  x11( width = 10.7, height = 10.5 )
  par(mfrow = c(2,2))
# Normality of standardised residuals
  qqnorm( r.Loss )
  abline( 0, 1, col = 'red' )

# Std residuals vs fitted values
  plot( yhat.lambda.Loss, r.Loss, xlab = "Fitted values", ylab = "Standardised residuals", 
        main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  abline( h = 0, col = 'red' )
  abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
  points( yhat.lambda.Loss[extr.Cdist], r.Loss[extr.Cdist], col = 'red', pch = 19 )
  points( yhat.lambda.Loss[idx.highlever], r.Loss[idx.highlever], col = 'blue' )
  #==# quite obvious pattern
  #==# High leverage not very interesting since just largest and smallest values of x...

# Absolute std residuals vs. fitted values (see heteroscedasticity)
  plot( yhat.lambda.Loss, abs(r.Loss), xlab = "Fitted values", ylab = "", 
        main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  title( ylab = expression( group("|","Standardised residuals","|") ), mgp = c(2.5,1,0) )
  points( yhat.lambda.Loss[extr.Cdist], abs(r.Loss)[extr.Cdist], col = 'red', pch = 19 )
  points( yhat.lambda.Loss[idx.highlever], abs(r.Loss)[idx.highlever], col = 'blue' )
  abline( h = 2, col = 'sandybrown', lty = 2 )
  
# Sqrt of absolute std residuals vs. fittet values
  plot( yhat.lambda.Loss, sqrt(abs(r.Loss)), xlab = "Fitted values", ylab = "", 
        main = paste( "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' ) )
  title( ylab = expression( sqrt(group("|","Standardised residuals","|") ) ), mgp = c(2.5,1,0) )
  points( yhat.lambda.Loss[extr.Cdist], sqrt(abs(r.Loss[extr.Cdist])), col = 'red', pch = 19 )
  points( yhat.lambda.Loss[idx.highlever], sqrt(abs(r.Loss[idx.highlever])), col = 'blue' )
  abline( h = sqrt(2), col = 'sandybrown', lty = 2 )
  savePlot( paste( "Figures/Loss_", modelname, "_resid_diagn", sep = "" ), type = "pdf" )
  #==# Full: Fit maybe not that good, but apart from outlier quite homoskedastic (smaller variance for large x-values).
  #==# Some problem of fit, but apart from outlier quite homoskedastic. Smaller variance for high leverage points
  #==#  with large x!

# Running variance along fitted values for block sizes 20, 26, 32
  x11( width = 9.8, height = 6.7 )
  plot( sort(yhat.lambda.Loss)[9 + seq(1, 365, 4)], sapply( seq(1, 365, 4), 
    function(i){ var( r.Loss[order(yhat.lambda.Loss)][i + 0:19], na.rm = TRUE ) } ),
    type = 'l', xlab = "Fitted values", ylab = "Running variance of standardised resdiuals", main = paste(
      "Absolute loss", ifelse(sep.lam.global, " seplam", ""), ", transformed scale", sep = '' )) # block size 20
  lines( sort(yhat.lambda.Loss)[ 12 + seq(1, 357, 4)], c( sapply( 
    seq(1,353,4), function(i){ var( r.Loss[order(yhat.lambda.Loss)][i + 0:25], na.rm = TRUE ) } ), 
    var( r.Loss[order(yhat.lambda.Loss)][357:n.obs], na.rm = TRUE ) ), col = 'blue' ) # block size 26
  lines( sort(yhat.lambda.Loss)[15 + seq(1,353,4)], sapply(seq(1,353,4), 
              function(i){ var( r.Loss[order(yhat.lambda.Loss)][i + 0:31], na.rm = TRUE ) } ), col = 'red' ) # block size 32
  abline( h = 1, col = 'darkgrey', lty = 2 )
  savePlot( paste( "Figures/Loss_", modelname, "_resid_runvar", sep = "" ), type = "pdf" )

# Compare residual plots for the two models (ONLY if separate lambda)
  # x11( width = 6.6, height = 7 )
  # plot( yhat.lambda.Loss, r.Loss, xlab = "Fitted values", ylab = "Standardised residuals", main = "Absolute loss, 
  #       separate lambda" )
  # abline( h = 0, col = 'red' )
  # savePlot( paste( "Figures/Loss_", modelname, "_resid_compare", sep = '' ), type = "pdf" )

# Leverage vs. Cook's distance vs. high residuals -> master pieces of plots
  x11( width = 13.2, height = 7 )
  par(mfrow = c(1, 2))
  # Distinguish outliers from leverage points (Davison 2008, p.395)
  plot( diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)), Cook.Loss, xlab = "h_{ii}/(1-h_{ii})", ylab = "Cook's distance", 
        main = "Absolute loss, transformed scale" )
  points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[extr.Cdist], Cook.Loss[extr.Cdist], col = 'red', pch = 19 )
  points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[idx.highlever], Cook.Loss[idx.highlever], 
          col = 'blue' ) # quite obvious pattern
  abline( h = 8/(n.obs - 2*p), col = 'red' )
  points( (diag(H.lambda.Loss)/(1-diag(H.lambda.Loss)))[lrg.resid], Cook.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
  abline( 0, 4/p, col = 'sandybrown', lty = 2 ) # 4 h/(p*(1-h)) = C
  abline( v = 2*p/(n.obs - 2*p), col = 'blue' )
  
# Residuals vs. leverage (like plot.lm but with better level curves)
  plot( diag(H.lambda.Loss), r.Loss, xlab = "Leverage", ylab = "Standardised residuals", 
        main = "Absolute loss, transformed scale", xlim = c(0,max(diag(H.lambda.Loss))) ) # xlim = c(0,0.025)
  abline( h = 0, col = 'grey' )
  abline( v = 2*p/n.obs, col = 'blue' )
  abline( h = c(-2,2), col = 'sandybrown', lty = 2 )
  points( diag(H.lambda.Loss)[extr.Cdist], r.Loss[extr.Cdist], col = 'red', pch = 19 )
  points( diag(H.lambda.Loss)[idx.highlever], r.Loss[idx.highlever], col = 'blue' ) # quite obvious pattern
  points( diag(H.lambda.Loss)[lrg.resid], r.Loss[lrg.resid], pch = 3, col = 'sandybrown' )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( (8/(n.obs-2*p)) * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red' )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.015 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 2 )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.01 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
  lines( c(h.seq,NA) %*% cbind(1,1), c(sqrt( 0.005 * p*(1-h.seq)/h.seq ),NA) %*% cbind(1,-1), col = 'red', lty = 3 )
  savePlot( paste( "Figures/Loss_", modelname, "_leverage", sep = "" ), type = "pdf" )

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
# --> By applying optimHess to mle.Loss I get the Hessian and thus the joint covariance at mle.Loss. Thus I do not need 
#     the opt version anymore. Since 2 seems better than 3, I keep 2 and 5.

fitted.mle.seq.Loss <- c( cbind( 1, xlam.seq.Loss ) %*% mle.Loss$beta.hat )

# (2) Assuming lambda = lambda.hat known, based on mle
  var.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ c(1,x) %*% covML.Loss[1:p,1:p] %*% c(1,x) } )
  var.predict.seq.Loss <- var.fitted.seq.Loss + s2.hat.Loss

# (3) Cov matrix of beta.hat from the regression of y.lambda on X.lambda with lambda = lambda.hat fixed
  varRegr.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ c(1,x) %*% covbeta.regr.Loss %*% c(1,x) } )

# (5) Detour via original scale to account for estimated lambda in y as well (but not sure it is properly done)
  varOrig.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(xlam){ 
    grad.yhat( theta = par.Loss, x.lam = c(1,xlam), sep.lam = sep.lam.global ) %*% 
      covML.Loss %*% grad.yhat( theta = par.Loss, x.lam = c(1,xlam), sep.lam = sep.lam.global ) } )

# (1) Including uncertainty on lambda.hat in x.lambda, based on opt
  # varFull.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ grad.fitted( par.Loss, c(1,x) ) %*% covML.Loss %*% 
  #     grad.fitted( par.Loss, c(1,x) ) } ) # sigma2 has no influence on this
  # varFull.predict.seq.Loss <- varFull.fitted.seq.Loss + s2.hat.Loss

# (4) As 1 but with cov(beta) based on non-parametric bootstrap (not good)
  # varBoot.fitted.seq.Loss <- sapply( xlam.seq.Loss, function(x){ 
  #   grad.fitted( par.Loss, c(1,x) ) %*% covBoot.Loss %*% grad.fitted( par.Loss, c(1,x) ) } )
  # varBoot.predict.seq.Loss <- varBoot.fitted.seq.Loss + s2.hat.Loss

  # x11()
  # plot( xlam.seq.Loss, varOrig.fitted.seq.Loss, type = 'l', xlab = "x(lambda)", ylab = "Variance of X(lambda)*beta.hat", 
  #       main = paste( "Absolute loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
  # # lines( xlam.seq.Loss, varBoot.fitted.seq.Loss, col = 'dodgerblue' )
  # lines( xlam.seq.Loss, var.fitted.seq.Loss, col = 'blue' )
  # lines( xlam.seq.Loss, varRegr.fitted.seq.Loss, col = 'cyan1' )
  # dev.off()

  x11()
# Transformed scale: y.lambda vs. X.lambda*beta
  plot( x.lambda.Loss, y.lambda.Loss, xlab = "Structure", ylab = "Contents", 
        main = paste( "Absolute loss", ifelse( sep.lam.global, " seplam", "" ), ", transformed scale", sep = '' ) )
  # if appropriate:
  # points( BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], Structure$Loss[idx.outl] ), 
  #         BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[idx.outl] ), pch = 4 ) # ev. pch = 8
  # points( BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], Structure$Loss[idx.remove] ), 
  #         BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[idx.remove] ), pch = 4 )
  #
  abline( mle.Loss$beta.hat, col = 'red' )

# Method 3:
  # 95% Confidence interval (mle, lambda fixed --> 3)
  lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ), col = 'orange', 
         lty = 2 )
  lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varRegr.fitted.seq.Loss ), col = 'orange', 
         lty = 2 )
  # 95% Prediction interval (mle, lambda fixed --> 3)
  lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ), 
         col = 'cyan1', lty = 2 )
  lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ), 
         col = 'cyan1', lty = 2 )

# Method 2:
  # 95% Confidence interval (mle, lambda known --> 2)
  lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), col = 'red', lty = 2 )
  lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), col = 'red', lty = 2 )
  # 95% Prediction interval (mle, lambda known --> 2)
  lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ), col = 'blue', lty = 2 )
  lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ), col = 'blue', lty = 2 )

# Method 5:
  # 95% Confidence interval (opt, via original scale --> 5)
  lines( xlam.seq.Loss, BC.transform( mle.Loss$lambda.hat[1], BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) - 
                                       qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss) ), col = 'sienna1', lty = 2 )
  lines( xlam.seq.Loss, BC.transform( mle.Loss$lambda.hat[1], BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) + 
                                       qnorm(0.975) * sqrt(varOrig.fitted.seq.Loss) ), col = 'sienna1', lty = 2 )
  # Prediction interval accounting for estimation uncertainty is too complex for this case -> neglect
  lines( lowess( x = x.lambda.Loss, y = y.lambda.Loss ), col = 'green' )
  box()

# Method 1:
  # # 95% Confidence interval (opt, lambda unknown --> 1)
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ), 
  #        col = 'lightpink2', lty = 2 )
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.fitted.seq.Loss ), 
  #        col = 'lightpink2', lty = 2 )  # use qnorm here???
  # # 95% Prediction interval (opt, lambda unknown --> 1)
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ), col = 'deepskyblue', lty = 2 )
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varFull.predict.seq.Loss ), col = 'deepskyblue', lty = 2 )

# Method 4 (neglected):
  # # 95% Confidence interval (mle, lambda unknown --> 4)
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.Loss ), 
  #        col = 'chocolate3', lty = 2 )
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varBoot.fitted.seq.Loss ), 
  #        col = 'chocolate3', lty = 2 )  # use qnorm here???
  # # 95% Prediction interval (mle, lambda unknown --> 4)
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.Loss ), 
  #        col = 'dodgerblue', lty = 2 )
  # lines( xlam.seq.Loss, fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( varBoot.predict.seq.Loss ), 
  #        col = 'dodgerblue', lty = 2 )
  # # --> Much less difference for prediction intervals because sigma2 dominates the variance.


# Either 2 or 3 seem reasonable since very close. 1 and 4 are not plausible because uncertainty of lambda is accounted 
# for x, but not for y. And the detour via the original scale is very close to 2 or 3.


##====================================================================##
##=== Backtransformed to original scale (!! only for full data !!) ===##
##====================================================================##
  x.seq.Loss <- BC.backtransform( mle.Loss$lambda.hat[1 + sep.lam.global], xlam.seq.Loss )
  
# Conditional mean on original scale, estimator by Taylor 1986
  condmean.seq.Loss <- condmean.orig( theta = par.Loss, x = cbind(1,xlam.seq.Loss), sep.lam = sep.lam.global )
  # use sigma2 or s2 here???
  var.condmean.Loss <- sapply( xlam.seq.Loss, function(x){
    c( grad.condmean.orig( theta = par.Loss, x = c(1,x), sep.lam = sep.lam.global ) %*% covML.Loss %*% 
         grad.condmean.orig( theta = par.Loss, x = c(1,x), sep.lam = sep.lam.global ) ) } )


  plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Contents", 
        main = paste( "Absolute loss", ifelse( sep.lam.global, " seplam", "" ), ", original scale", sep = '' ) )
  polygon(x = c( 1.5, 1, 1, -1, -1, 1.5 ), y =  c( -1, -1, 1, 1, 1.5, 1.5 ), density = 5, angle = -45, 
          col = "lightgrey")
  
# if appropriate:
  points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4 ) # ev. pch = 8
  # points( Structure$Loss[idx.remove], Contents$Loss[idx.remove], pch = 4 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ), col = 'red' ) # fitted median 
  # lines( x.seq.Loss, BC.backtransform( opt.Loss$par[p+2], fitted.opt.seq.Loss ), col = 'lightpink2' )

# Method 3:
  # Backtransformed 95% Confidence interval (mle, lambda fixed --> 3)
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                        sqrt( varRegr.fitted.seq.Loss ) ), col = 'orange', lty = 2 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                        sqrt( varRegr.fitted.seq.Loss ) ), col = 'orange', lty = 2 )
  # Backtransformed 95% Prediction interval (mle, lambda fixed --> 3)
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                        sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ) ), col = 'cyan1', lty = 2 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                        sqrt( s2.hat.Loss + varRegr.fitted.seq.Loss ) ), col = 'cyan1', lty = 2 )

# Method 2:   
  # Backtransformed 95% Confidence interval (mle, lambda known --> 2)
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.fitted.seq.Loss ) ), col = 'red', lty = 2 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.fitted.seq.Loss ) ), col = 'red', lty = 2 )
# Backtransformed 95% Prediction interval (mle, lambda known --> 2)
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) ), col = 'blue', lty = 2 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) ), col = 'blue', lty = 2 )

# Method 5:  
  # Backtransformed 95% Confidence interval (opt, var of back-transformed fitted value --> 5) 
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) - qnorm(0.975) * 
           sqrt(varOrig.fitted.seq.Loss), col = 'sienna1', lty = 2 )
  lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ) + qnorm(0.975) * 
           sqrt(varOrig.fitted.seq.Loss), col = 'sienna1', lty = 2 )

# # Method 1:
#   # Backtransformed 95% Confidence interval (opt, lambda unknown --> 1)
#   lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p ) *
#                                         sqrt( varFull.fitted.seq.Loss ) ), col = 'lightpink2', lty = 2 )
#   lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p ) *
#                                         sqrt( varFull.fitted.seq.Loss ) ), col = 'lightpink2', lty = 2 )  # use qnorm here???
# # Backtransformed 95% Prediction interval (opt, lambda unknown --> 1)
#   lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p ) * 
#                                         sqrt( varFull.predict.seq.Loss ) ), col = 'deepskyblue', lty = 2 )
#   lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p ) * 
#                                         sqrt( varFull.predict.seq.Loss ) ), col = 'deepskyblue', lty = 2 )

  
##==========================================##
##=== Conditional mean on original scale ===##
##==========================================##
  
  lines( x.seq.Loss, condmean.seq.Loss, col = 'green3', lwd = 1 )
  lines( x.seq.Loss, condmean.seq.Loss - qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )
  lines( x.seq.Loss, condmean.seq.Loss + qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )
  box()

  #<------------------------------------------------- Experimental ------------------------------------------------------> 
# # Smearing estimate of the mean on original scale (Duan, 1983)
#   lines( x.seq.Loss, sapply( xlam.seq.Loss, function(x){ 
#     mean( BC.backtransform( mle.Loss$lambda.hat[1], sum( c(1,x)*mle.Loss$beta.hat ) + e.Loss ) ) } ), 
#     col = 'mediumpurple1' )
#   # Basically identical with condmean.seq
#   max( abs( condmean.seq.Loss - sapply( xlam.seq.Loss, function(x){ 
#     mean( BC.backtransform( mle.Loss$lambda.hat[1], sum( c(1,x)*mle.Loss$beta.hat ) + e.Loss ) ) } ) ), na.rm = TRUE )
# 
#   
# # Apply to confidence interval 
#   # "Confidence interval" for the mean on orig. scale: Use 95% CI bounds as x instead of the fitted values (1st try)
#   lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss - 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), 
#          col = 'magenta', lty = 2 )
#   lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss + 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ) ), 
#          col = 'magenta', lty = 2 )
#   # "Prediction interval" for the mean on original scale (obtained as before)
#   lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss - 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), 
#          col = 'skyblue', lty = 2 )
#   lines( x.seq.Loss, condmean.orig( theta = unlist(mle.Loss)[c(2:4,1)], x.beta = fitted.mle.seq.Loss + 
#                                      qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss ) ), 
#          col = 'skyblue', lty = 2 )
#<------------------------------------------------- Experimental ------------------------------------------------------> 

##==========================================##
##=== Alternative estimators for lambda ====##
##==========================================##
  
# Sensitivity to removed data points
  # lambda transforming to a symmetric distribution (p.135 Carroll and Ruppert 1984a)
  plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), T.skew, x = x.Loss, y = y.Loss, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Absolute loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = T.skew, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE )
  # robust: lambda_sk = 0.344, w/o outlier 0.339

  plot( seq(-2,1,0.02), sapply( seq(-2,1,0.02), T.skew, x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE ), type = 'l', xlab = "lambda", ylab = "Skewness measure", main = "Absolute loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = T.skew, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE, robust = FALSE )
  # less robust: lambda_sk = 0.179

# lambda transforming to a homoskedastic distribution (p.135 Carroll and Ruppert 1984a)
  plot( seq(-2,2,0.02), sapply( seq(-2,2,0.02), H.hesk, x = x.Loss, y = y.Loss, intercept = TRUE ), type = 'l', xlab = "lambda", ylab = "Heteroscedasticity measure", main = "Absolute loss" )
  abline( h = 0, col = 'red' )
  uniroot( f = H.hesk, interval = c(-1,1), x = x.Loss, y = y.Loss, intercept = TRUE )
# robust: lambda_hesk = 0.136, w/o outlier 0.130


# Adjusted R^2:
  adj_R2.Loss <- 1 - ( sum(e.Loss^2)/(n.obs-p) )/( sum( ( y.lambda.Loss - mean(y.lambda.Loss) )^2 )/(n.obs-1) )
  
  
  
  
  
  

