########################################################################################################################
########################################### Analysis of variance wrt cantons ###########################################
########################################################################################################################

# Observations to omit
  which( Contents$Region %in% c("GE","AI") )
  idx.omit <- c( idx.outl, which( Contents$Region %in% c("GE","AI") ) )

  region.aov <- Contents$Region[-idx.omit]
  
  Region.mat.aov <- Region.mat[,-which(regionnames %in% c("GE","AI"))][-idx.omit,]
  n.obs.aov <- nrow(Region.mat.aov)
  
  n.obs.aov.cant <- as.vector( table( region.aov ) )
  
  pairdiff.mat <- rbind( SZ.OW = c(-1,1,0,0,0), TI.OW = c(-1,0,1,0,0), UR.OW = c(-1,0,0,1,0), VS.OW = c(-1,0,0,0,1), 
                         TI.SZ = c(0,-1,1,0,0), UR.SZ = c(0,-1,0,1,0), VS.SZ = c(0,-1,0,0,1), 
                         UR.TI = c(0,0,-1,1,0), VS.TI = c(0,0,-1,0,1), 
                         VS.UR = c(0,0,0,-1,1) )


  
##=======================================##
##==== Relative loss ====================##
##=======================================##

# Fit the three models
  mle.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = FALSE )
  mle.seplam.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = TRUE )
  mle.other.DoL <- TBS.mle( y = y.DoL, x = x.DoL, beta.init = c(0.1,1) )

# LRT sep.lam vs. single lambda
  1 - pchisq( -2*( mle.DoL$ell - mle.seplam.DoL$ell ), df = 1 )
  
# AIC:
  -2*mle.DoL$ell + 2*length(unlist(mle.DoL[1:3])) # -1877.458
  -2*mle.seplam.DoL$ell + 2*length(unlist(mle.seplam.DoL[1:3])) # -1882.11
  -2*mle.other.DoL$ell + 2*length(unlist(mle.other.DoL[1:3])) # -1883.873
  
# BIC:
  -2*mle.DoL$ell + log(n.obs)*length(unlist(mle.DoL[1:3])) # -1861.666
  -2*mle.seplam.DoL$ell + log(n.obs)*length(unlist(mle.seplam.DoL[1:3])) # -1862.37
  -2*mle.other.DoL$ell + log(n.obs)*length(unlist(mle.other.DoL[1:3])) # -1868.081


# Redefine quantities (including xy.lambda) with appropriate observations omitted
# Original scale:
  y.aov.DoL <- Contents$DoL[-idx.omit]
  x.aov.DoL <- Structure$DoL[-idx.omit]
  
  X.canton.aov.DoL <- cbind( Region.mat.aov, Region.mat.aov * x.aov.DoL )
  colnames(X.canton.aov.DoL)[6:10] <- paste( 'x', colnames(X.canton.aov.DoL)[1:5], sep = ':' )

#(A) Take lambda as fixed:
  # (A.1) PTBS model (need to redefine xy.lambda first --> execute whole block)
    # Transformed scale (lambda fixed to its fit over the whole data set)
      y.lambda.aov.DoL <- BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[-idx.omit] )
      x.lambda.aov.DoL <- BC.transform( mle.DoL$lambda.hat[1+0], Structure$DoL[-idx.omit] )
  
    #X.lambda.canton.DoL <- cbind( Region.mat, Region.mat * x.lambda.DoL )
      X.lambda.canton.aov.DoL <- cbind( Region.mat.aov, Region.mat.aov * x.lambda.aov.DoL )
      colnames(X.lambda.canton.aov.DoL)[6:10] <- paste( 'x.lam', colnames(X.lambda.canton.aov.DoL)[1:5], sep = ':' )
    
    # Intercept and slope by canton
      canton.both.DoL <- my.lm( y.lambda.aov.DoL, X.lambda.canton.aov.DoL )
      names( canton.both.DoL$coefs ) <- colnames(X.lambda.canton.aov.DoL)
    # Only intercept by canton
      canton.intcpt.DoL <- my.lm( y.lambda.aov.DoL, cbind( Region.mat.aov, x.lam = x.lambda.aov.DoL ) )
      names( canton.intcpt.DoL$coefs ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      canton.slope.DoL <- my.lm( y.lambda.aov.DoL, cbind( intcpt = 1, X.lambda.canton.aov.DoL[,6:10] ) )
      names( canton.slope.DoL$coefs ) <- c( 'one', colnames(X.lambda.canton.aov.DoL[,6:10]) )
    # Nothing by canton
      canton.no.DoL <- my.lm( y.lambda.aov.DoL, cbind( intcpt = 1, x.lam = x.lambda.aov.DoL ) )
      names( canton.no.DoL$coefs ) <- c('one','x')
    

  # (A.2) PTBS model with seplam (need to redefine xy.lambda first --> execute whole block)
    # Transformed scale (lambda fixed to its fit over the whole data set)
      y.lambda.aov.DoL <- BC.transform( mle.seplam.DoL$lambda.hat[1], Contents$DoL[-idx.omit] )
      x.lambda.aov.DoL <- BC.transform( mle.seplam.DoL$lambda.hat[1+1], Structure$DoL[-idx.omit] )
    # X.lambda.canton.DoL <- cbind( Region.mat, Region.mat * x.lambda.DoL )
      X.lambda.canton.aov.DoL <- cbind( Region.mat.aov, Region.mat.aov * x.lambda.aov.DoL )
      colnames(X.lambda.canton.aov.DoL)[6:10] <- paste( 'x.lam', colnames(X.lambda.canton.aov.DoL)[1:5], sep = ':' )
    # Intercept and slope by canton
      canton.both.seplam.DoL <- my.lm( y.lambda.aov.DoL, X.lambda.canton.aov.DoL )
      names( canton.both.seplam.DoL$coefs ) <- colnames(X.lambda.canton.aov.DoL)
    # Only intercept by canton
      canton.intcpt.seplam.DoL <- my.lm( y.lambda.aov.DoL, cbind( Region.mat.aov, x.lam = x.lambda.aov.DoL ) )
      names( canton.intcpt.seplam.DoL$coefs ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      canton.slope.seplam.DoL <- my.lm( y.lambda.aov.DoL, cbind( intcpt = 1, X.lambda.canton.aov.DoL[,6:10] ) )
      names( canton.slope.seplam.DoL$coefs ) <- c( 'one', colnames(X.lambda.canton.aov.DoL[,6:10]) )
    # Nothing by canton
      canton.no.seplam.DoL <- my.lm( y.lambda.aov.DoL, cbind( intcpt = 1, x.lam = x.lambda.aov.DoL ) )
      names( canton.no.seplam.DoL$coefs ) <- c('one','x')
    

  # (A.3) TBS model (also lambda fixed; model is on original scale)
    # Intercept and slope by canton
      canton.both.other.DoL <- TBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL, lambda = mle.other.DoL$lambda.hat, 
                                        beta.init = rep( c(0.1,1), each = 5 ), intercept = FALSE )
     names( canton.both.other.DoL$beta.hat ) <- colnames(X.canton.aov.DoL)
    # Only intercept by canton
      canton.intcpt.other.DoL <- TBS.mle( y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), 
                                          lambda = mle.other.DoL$lambda.hat, beta.init = c(rep(0.1,5),1), intercept = FALSE )
      names( canton.intcpt.other.DoL$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      canton.slope.other.DoL <- TBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], lambda = mle.other.DoL$lambda.hat, 
                                         beta.init = c(0.1,rep(1,5)), intercept = TRUE )
      names( canton.slope.other.DoL$beta.hat ) <- c( 'one', colnames(X.canton.aov.DoL[,6:10]) )
    # Nothing by canton
      canton.no.other.DoL <- TBS.mle( y = y.aov.DoL, x = x.aov.DoL, lambda = mle.other.DoL$lambda.hat, 
                                      beta.init = c(0.1,1), intercept = TRUE )
      names( canton.no.other.DoL$beta.hat ) <- c('one','x')

  #=== Comparisons (lambda fixed)
  # For regression models use F-tests corresponding to anova( lm.simple, lm.complex )
    RSS.no.DoL <- sum(canton.no.DoL$e^2)
    RSS.intcpt.DoL <- sum(canton.intcpt.DoL$e^2)
    RSS.slope.DoL <- sum(canton.slope.DoL$e^2)
    RSS.both.DoL <- sum(canton.both.DoL$e^2)
    
    RSS.no.seplam.DoL <- sum(canton.no.seplam.DoL$e^2)
    RSS.intcpt.seplam.DoL <- sum(canton.intcpt.seplam.DoL$e^2)
    RSS.slope.seplam.DoL <- sum(canton.slope.seplam.DoL$e^2)
    RSS.both.seplam.DoL <- sum(canton.both.seplam.DoL$e^2)

  # (a) both vs. nothing
    1 - pf( ( ( RSS.no.DoL - RSS.both.DoL )/(canton.both.DoL$p - canton.no.DoL$p) )/canton.both.DoL$s2, 
            df1 = canton.both.DoL$p - canton.no.DoL$p, df2 = canton.both.DoL$n - canton.both.DoL$p )
    1 - pf( ( ( RSS.no.seplam.DoL - RSS.both.seplam.DoL )/(canton.both.seplam.DoL$p - canton.no.seplam.DoL$p) )/
              canton.both.seplam.DoL$s2, df1 = canton.both.seplam.DoL$p - canton.no.seplam.DoL$p, df2 = canton.both.seplam.DoL$n - canton.both.seplam.DoL$p )
    1 - pchisq( -2*( canton.no.other.DoL$ell - canton.both.other.DoL$ell ), df = ncol(X.canton.aov.DoL) - 2 )
    # just not signif. (0.057); seplam: not signif.; TBS: not signif.
  # (b) both vs. intercept
    1 - pf( ( ( RSS.intcpt.DoL - RSS.both.DoL )/(canton.both.DoL$p - canton.intcpt.DoL$p) )/canton.both.DoL$s2,
            df1 = canton.both.DoL$p - canton.intcpt.DoL$p, df2 = canton.both.DoL$n - canton.both.DoL$p )
    1 - pf( ( ( RSS.intcpt.seplam.DoL - RSS.both.seplam.DoL )/(canton.both.seplam.DoL$p - canton.intcpt.seplam.DoL$p) )/
              canton.both.seplam.DoL$s2, df1 = canton.both.seplam.DoL$p - canton.intcpt.seplam.DoL$p, 
            df2 = canton.both.seplam.DoL$n - canton.both.seplam.DoL$p )
    1 - pchisq( -2*( canton.intcpt.other.DoL$ell - canton.both.other.DoL$ell ), df = ncol(X.canton.aov.DoL) - 
                  (ncol(Region.mat.aov)+1) )
    # not signif. (neither for seplam nor TBS)
  # (c) both vs. slope
    1 - pf( ( ( RSS.slope.DoL - RSS.both.DoL )/(canton.both.DoL$p - canton.slope.DoL$p) )/canton.both.DoL$s2, 
            df1 = canton.both.DoL$p - canton.slope.DoL$p, df2 = canton.both.DoL$n - canton.both.DoL$p )
    1 - pf( ( ( RSS.slope.seplam.DoL - RSS.both.seplam.DoL )/(canton.both.seplam.DoL$p - canton.slope.seplam.DoL$p) )/
              canton.both.seplam.DoL$s2, df1 = canton.both.seplam.DoL$p - canton.slope.seplam.DoL$p, 
            df2 = canton.both.seplam.DoL$n - canton.both.seplam.DoL$p )
    1 - pchisq( -2*( canton.slope.other.DoL$ell - canton.both.other.DoL$ell ), df = ncol(X.canton.aov.DoL) - 
                  (ncol(Region.mat.aov)+1) )
    # not signif. (neither for seplam nor TBS)
  # (d) intercept vs. nothing 
    1 - pf( ( ( RSS.no.DoL - RSS.intcpt.DoL )/(canton.intcpt.DoL$p - canton.no.DoL$p) )/canton.intcpt.DoL$s2, 
            df1 = canton.intcpt.DoL$p - canton.no.DoL$p, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
    1 - pf( ( ( RSS.no.seplam.DoL - RSS.intcpt.seplam.DoL )/(canton.intcpt.seplam.DoL$p - canton.no.seplam.DoL$p) )/
              canton.intcpt.seplam.DoL$s2, df1 = canton.intcpt.seplam.DoL$p - canton.no.seplam.DoL$p, 
            df2 = canton.intcpt.seplam.DoL$n - canton.intcpt.seplam.DoL$p )
    1 - pchisq( -2*( canton.no.other.DoL$ell - canton.intcpt.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
    # signif. (also for seplam); not signif. for TBS
  # (e) slope vs. nothing
    1 - pf( ( ( RSS.no.DoL - RSS.slope.DoL )/(canton.slope.DoL$p - canton.no.DoL$p) )/canton.slope.DoL$s2, 
            df1 = canton.slope.DoL$p - canton.no.DoL$p, df2 = canton.slope.DoL$n - canton.slope.DoL$p )
    1 - pf( ( ( RSS.no.seplam.DoL - RSS.slope.seplam.DoL )/(canton.slope.seplam.DoL$p - canton.no.seplam.DoL$p) )/
              canton.slope.seplam.DoL$s2, df1 = canton.slope.seplam.DoL$p - canton.no.seplam.DoL$p, 
            df2 = canton.slope.seplam.DoL$n - canton.slope.seplam.DoL$p )
    1 - pchisq( -2*( canton.no.other.DoL$ell - canton.slope.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
    # just not signif. (0.061); just not signif. for seplam (0.059); not signif. for TBS (0.089)
  


# (B) Fit new PTBS/TBS model:
  # (B.1) PTBS model
    # Intercept and slope by canton
      mle.canton.both.DoL <- PTBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL, intercept = 1:5, dummy.zeros = TRUE, 
                                       sep.lam = FALSE )
      names( mle.canton.both.DoL$beta.hat ) <- colnames(X.canton.aov.DoL)
    # Only intercept by canton
      mle.canton.intcpt.DoL <- PTBS.mle( y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), intercept = 1:5, 
                                         sep.lam = FALSE )
      names( mle.canton.intcpt.DoL$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.DoL <- PTBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], intercept = TRUE, dummy.zeros = TRUE, 
                                        sep.lam = FALSE ) # ! dummy.zeros only works if all(x.aov.DoL != 0) !
      names( mle.canton.slope.DoL$beta.hat ) <- c( 'one', colnames(X.canton.aov.DoL[, 6:10]) )
    # Nothing by canton
      mle.canton.no.DoL <- PTBS.mle( y = y.aov.DoL, x = x.aov.DoL, intercept = TRUE, sep.lam = FALSE )
      names(mle.canton.no.DoL$beta.hat) <- c('one','x')

    # Covariance matrix of coefficients:
      covML.canton.both.DoL <- solve( -optimHess( unlist(mle.canton.both.DoL[c(2,3,1)]), PTBS.llkhd, y = y.aov.DoL, 
                                                  x = X.canton.aov.DoL, intercept = 1:5, dummy.zeros = TRUE, 
                                                  sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.intcpt.DoL <- solve( -optimHess( unlist(mle.canton.intcpt.DoL[c(2,3,1)]), PTBS.llkhd, y = y.aov.DoL, 
                                                    x = cbind( Region.mat.aov, x.aov.DoL ), intercept = 1:5, 
                                                    sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.slope.DoL <- solve( -optimHess( unlist(mle.canton.slope.DoL[c(2,3,1)]), PTBS.llkhd, y = y.aov.DoL, 
                                                   x = X.canton.aov.DoL[,6:10], intercept = TRUE, dummy.zeros = TRUE, 
                                                   sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.no.DoL <- solve( -optimHess( unlist(mle.canton.no.DoL[c(2,3,1)]), PTBS.llkhd, y = y.aov.DoL, 
                                                x = x.aov.DoL, intercept = TRUE, sep.lam = FALSE, 
                                                control = list( maxit = 5000, fnscale = -1 ) ) )

    # Standard errors of coefficients:
      se.canton.both.DoL <- sqrt( diag(covML.canton.both.DoL)[1:length(mle.canton.both.DoL$beta.hat)] )
      se.canton.intcpt.DoL <- sqrt( diag(covML.canton.intcpt.DoL)[1:length(mle.canton.intcpt.DoL$beta.hat)] )
      se.canton.slope.DoL <- sqrt( diag(covML.canton.slope.DoL)[1:length(mle.canton.slope.DoL$beta.hat)] )
      se.canton.no.DoL <- sqrt( diag(covML.canton.no.DoL)[1:p] )


  # (B.2) PTBS model with seplam
    # Intercept and slope by canton
      mle.canton.both.seplam.DoL <- PTBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL, intercept = 1:5, dummy.zeros = TRUE, 
                                              sep.lam = TRUE )
      names( mle.canton.both.seplam.DoL$beta.hat ) <- colnames(X.canton.aov.DoL)
    # Only intercept by canton
      mle.canton.intcpt.seplam.DoL <- PTBS.mle( y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), intercept = 1:5, 
                                                sep.lam = TRUE )
      names( mle.canton.intcpt.seplam.DoL$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.seplam.DoL <- PTBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], intercept = TRUE, 
                                               dummy.zeros = TRUE, sep.lam = TRUE ) # ! dummy.zeros only works if all(x.aov.DoL != 0) !
      names( mle.canton.slope.seplam.DoL$beta.hat ) <- c( 'one', colnames(X.canton.aov.DoL[,6:10]) )
    # Nothing by canton
      mle.canton.no.seplam.DoL <- PTBS.mle( y = y.aov.DoL, x = x.aov.DoL, intercept = TRUE, sep.lam = TRUE )
      names(mle.canton.no.seplam.DoL$beta.hat) <- c('one','x')

    # Covariance matrix of coefficients:
      covML.canton.both.seplam.DoL <- solve( -optimHess( unlist(mle.canton.both.seplam.DoL[c(2,3,1)]), PTBS.llkhd, 
                                                         y = y.aov.DoL, x = X.canton.aov.DoL, intercept = 1:5, 
                                                         dummy.zeros = TRUE, sep.lam = TRUE, 
                                                         control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.intcpt.seplam.DoL <- solve( -optimHess( unlist(mle.canton.intcpt.seplam.DoL[c(2,3,1)]), PTBS.llkhd, 
                                                           y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), 
                                                           intercept = 1:5, sep.lam = TRUE, 
                                                           control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.slope.seplam.DoL <- solve( -optimHess( unlist(mle.canton.slope.seplam.DoL[c(2,3,1)]), PTBS.llkhd, 
                                                          y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], intercept = TRUE, 
                                                          dummy.zeros = TRUE, sep.lam = TRUE, 
                                                          control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.no.seplam.DoL <- solve( -optimHess( unlist(mle.canton.no.seplam.DoL[c(2,3,1)]), PTBS.llkhd, 
                                                       y = y.aov.DoL, x = x.aov.DoL, intercept = TRUE, sep.lam = TRUE, 
                                                       control = list( maxit = 5000, fnscale = -1 ) ) )
    
    # Standard errors of coefficients:
      se.canton.both.seplam.DoL <- sqrt( diag(covML.canton.both.seplam.DoL)[1:length(mle.canton.both.seplam.DoL$beta.hat)] )
      se.canton.intcpt.seplam.DoL <- sqrt( diag(covML.canton.intcpt.seplam.DoL)[1:length(mle.canton.intcpt.seplam.DoL$beta.hat)] )
      se.canton.slope.seplam.DoL <- sqrt( diag(covML.canton.slope.seplam.DoL)[1:length(mle.canton.slope.seplam.DoL$beta.hat)] )
      se.canton.no.seplam.DoL <- sqrt( diag(covML.canton.no.seplam.DoL)[1:p] )

      
  # (B.3) TBS model
    # Both by canton
      mle.canton.both.other.DoL <- TBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL, beta.init =rep( c(0.1,1), each = 5 ), 
                                            intercept = FALSE )
      names( mle.canton.both.other.DoL$beta.hat ) <- colnames(X.canton.aov.DoL)
    # Only intercept by canton
      mle.canton.intcpt.other.DoL <- TBS.mle( y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), 
                                              beta.init = c( rep(0.1,5), 1 ), intercept = FALSE )
      names( mle.canton.intcpt.other.DoL$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.other.DoL <- TBS.mle( y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], beta.init = c( 0.1, rep(1,5) ), 
                                             intercept = TRUE )
      names( mle.canton.slope.other.DoL$beta.hat ) <- c( 'one', colnames(X.canton.aov.DoL[,6:10]) )
    # Nothing by canton
      mle.canton.no.other.DoL <- TBS.mle( y = y.aov.DoL, x = x.aov.DoL, beta.init = c(0.1,1), intercept = TRUE )
      names(mle.canton.no.other.DoL$beta.hat) <- c('one','x')

    # Covariance matrix of coefficients:
      covML.canton.both.other.DoL <- solve( -optimHess( unlist(mle.canton.both.other.DoL[c(2,3,1)]), loglik.other, 
                                                        y = y.aov.DoL, x = X.canton.aov.DoL, intercept = FALSE, 
                                                        control = list( fnscale = -1 ) ) )
      covML.canton.intcpt.other.DoL <- solve( -optimHess( unlist(mle.canton.intcpt.other.DoL[c(2,3,1)]), loglik.other, 
                                                          y = y.aov.DoL, x = cbind( Region.mat.aov, x.aov.DoL ), 
                                                          intercept = FALSE, control = list( fnscale = -1 ) ) )
      covML.canton.slope.other.DoL <- solve( -optimHess( unlist(mle.canton.slope.other.DoL[c(2,3,1)]), loglik.other, 
                                                         y = y.aov.DoL, x = X.canton.aov.DoL[,6:10], intercept = TRUE, 
                                                         control = list( fnscale = -1 ) ) )
      covML.canton.no.other.DoL <- solve( -optimHess( unlist(mle.canton.no.other.DoL[c(2,3,1)]), loglik.other, 
                                                      y = y.aov.DoL, x = x.aov.DoL, intercept = TRUE, 
                                                      control = list( fnscale = -1 ) ) )
    
    # Standard errors of coefficients:
    se.canton.both.other.DoL <- sqrt( diag(covML.canton.both.other.DoL)[1:length(mle.canton.both.other.DoL$beta.hat)] )
    se.canton.intcpt.other.DoL <- sqrt( diag(covML.canton.intcpt.other.DoL)[1:length(mle.canton.intcpt.other.DoL$beta.hat)] )
    se.canton.slope.other.DoL <- sqrt( diag(covML.canton.slope.other.DoL)[1:length(mle.canton.slope.other.DoL$beta.hat)] )
    se.canton.no.other.DoL <- sqrt( diag(covML.canton.no.other.DoL)[1:p] )
    

    #=== Comparisons (newly fitted models)
    #(a) both vs. nothing
      1 - pchisq( -2*( mle.canton.no.DoL$ell - mle.canton.both.DoL$ell ), df = ncol(X.canton.aov.DoL) - 2 )
      1 - pchisq( -2*( mle.canton.no.seplam.DoL$ell - mle.canton.both.seplam.DoL$ell ), df = ncol(X.canton.aov.DoL) - 2 )
      1 - pchisq( -2*( mle.canton.no.other.DoL$ell - mle.canton.both.other.DoL$ell ), df = ncol(X.canton.aov.DoL) - 2 )
      # just signif. (0.048); seplam: not signif.; TBS: not signif.
    #(b) both vs. intercept
      1 - pchisq( -2*( mle.canton.intcpt.DoL$ell - mle.canton.both.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.intcpt.seplam.DoL$ell - mle.canton.both.seplam.DoL$ell ),df = ncol(Region.mat.aov) - 1)
      1 - pchisq( -2*( mle.canton.intcpt.other.DoL$ell - mle.canton.both.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      # not signif. (neither for seplam nor TBS)
    #(c) both vs. slope
      1 - pchisq( -2*( mle.canton.slope.DoL$ell - mle.canton.both.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.slope.seplam.DoL$ell - mle.canton.both.seplam.DoL$ell ),df = ncol(Region.mat.aov) - 1)
      1 - pchisq( -2*( mle.canton.slope.other.DoL$ell - mle.canton.both.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      # not signif. (neither for seplam nor TBS)
    #(d) intercept vs. nothing 
      1 - pchisq( -2*( mle.canton.no.DoL$ell - mle.canton.intcpt.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.no.seplam.DoL$ell - mle.canton.intcpt.seplam.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.no.other.DoL$ell - mle.canton.intcpt.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      # signif. (also for seplam); not signif. for TBS
    #(e) slope vs. nothing 
      1 - pchisq( -2*( mle.canton.no.DoL$ell - mle.canton.slope.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.no.seplam.DoL$ell - mle.canton.slope.seplam.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      1 - pchisq( -2*( mle.canton.no.other.DoL$ell - mle.canton.slope.other.DoL$ell ), df = ncol(Region.mat.aov) - 1 )
      # just not signif. (0.057) (also for seplam: 0.055); not signif. for TBS


  # Overview tables of the individual model coefficients (only for new estimations B; not adapted for the three models)
    # both
      cbind( estim = mle.canton.both.DoL$beta.hat, se = se.canton.both.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.both.DoL$beta.hat/se.canton.both.DoL) ) ) )
      cbind( estim = mle.canton.both.seplam.DoL$beta.hat, se = se.canton.both.seplam.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.both.seplam.DoL$beta.hat/se.canton.both.seplam.DoL) ) ) )
      cbind( estim = mle.canton.both.other.DoL$beta.hat, se = se.canton.both.other.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.both.other.DoL$beta.hat/se.canton.both.other.DoL) ) ) )
    # intercept
      cbind( estim = mle.canton.intcpt.DoL$beta.hat, se = se.canton.intcpt.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.DoL$beta.hat/se.canton.intcpt.DoL) ) ) )
      cbind( estim = mle.canton.intcpt.seplam.DoL$beta.hat, se = se.canton.intcpt.seplam.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.seplam.DoL$beta.hat/se.canton.intcpt.seplam.DoL) ) ) )
      cbind( estim = mle.canton.intcpt.other.DoL$beta.hat, se = se.canton.intcpt.other.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.other.DoL$beta.hat/se.canton.intcpt.other.DoL) ) ) )
    # slope
      cbind( estim = mle.canton.slope.DoL$beta.hat, se = se.canton.slope.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.slope.DoL$beta.hat/se.canton.slope.DoL) ) ) )
      cbind( estim = mle.canton.slope.seplam.DoL$beta.hat, se = se.canton.slope.seplam.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.slope.seplam.DoL$beta.hat/se.canton.slope.seplam.DoL) ) ) )
      cbind( estim = mle.canton.slope.other.DoL$beta.hat, se = se.canton.slope.other.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.slope.other.DoL$beta.hat/se.canton.slope.other.DoL) ) ) )
    # nothing
      cbind( estim = mle.canton.no.DoL$beta.hat, se = se.canton.no.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.no.DoL$beta.hat/se.canton.no.DoL) ) ) )
      cbind( estim = mle.canton.no.seplam.DoL$beta.hat, se = se.canton.no.seplam.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.no.seplam.DoL$beta.hat/se.canton.no.seplam.DoL) ) ) )
      cbind( estim = mle.canton.no.other.DoL$beta.hat, se = se.canton.no.other.DoL, 
             pval = 2*( 1 - pnorm( abs(mle.canton.no.other.DoL$beta.hat/se.canton.no.other.DoL) ) ) )


  # Differences between coefficients
    # Only intercept:
      diffs.intcpt.DoL <- pairdiff.mat %*% mle.canton.intcpt.DoL$beta.hat[1:5]
      se.diffs.intcpt.DoL <- sqrt( diag( pairdiff.mat %*% covML.canton.intcpt.DoL[1:5,1:5] %*% t(pairdiff.mat) ) )
      diffs.intcpt.seplam.DoL <- pairdiff.mat %*% mle.canton.intcpt.seplam.DoL$beta.hat[1:5]
      se.diffs.intcpt.seplam.DoL <- sqrt( diag( pairdiff.mat %*% covML.canton.intcpt.seplam.DoL[1:5,1:5] %*% 
                                                  t(pairdiff.mat) ) )
      diffs.intcpt.other.DoL <- pairdiff.mat %*% mle.canton.intcpt.other.DoL$beta.hat[1:5]
      se.diffs.intcpt.other.DoL <- sqrt( diag( pairdiff.mat %*% covML.canton.intcpt.other.DoL[1:5,1:5] %*% 
                                                 t(pairdiff.mat) ) )

    # Overview tables !!! p-values with caution !!!
      cbind( estim = c(diffs.intcpt.DoL), se = se.diffs.intcpt.DoL, 
             pval = 2*( 1 - pnorm( abs(c(diffs.intcpt.DoL))/se.diffs.intcpt.DoL ) ) )
      # VS-OW and VS-TI largest; UR-OW is also significant but only a bit more than half as large as the two others; 
      #   SZ-OW is of similar value as UR-OW but just not significant

      cbind( estim = c(diffs.intcpt.seplam.DoL), se = se.diffs.intcpt.seplam.DoL, 
             pval = 2*( 1 - pnorm( abs(c(diffs.intcpt.seplam.DoL))/se.diffs.intcpt.seplam.DoL ) ) )
      # VS-OW and VS-TI largest; UR-OW is also significant but only a bit more than half as large as the two others; 
      #   other differences of similar magnitude are not significant (inlcluding SZ-OW)
  
      cbind( estim = c(diffs.intcpt.other.DoL), se = se.diffs.intcpt.other.DoL, 
             pval = 2*( 1 - pnorm( abs(c(diffs.intcpt.other.DoL))/se.diffs.intcpt.other.DoL ) ) )
      # TBS: UR-TI and VS-TI are significant and also the largest two differences (together with VS-OW which is not 
      #   significant), but all differences are very small (abs <= 0.013). Thus there are two significant differences in 
      #   the medians of y, but overall canton-specific intercepts do not seem appropriate by the LRT.
      # Moreover, not so sure that considering differences of betas makes sense for this model!!

  # Only slope (by curiosity):
    diffs.slope.DoL <- pairdiff.mat %*% mle.canton.slope.DoL$beta.hat[2:6]
    rownames(diffs.slope.DoL) <- sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ), 
                                                 function(cn){ paste( 'x', cn, sep = ':' ) } ), paste, collapse = "." )
    se.diffs.slope.DoL <- sqrt( diag( pairdiff.mat %*% covML.canton.slope.DoL[2:6,2:6] %*% t(pairdiff.mat) ) )
    names(se.diffs.slope.DoL) <- sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ), 
                                                 function(cn){ paste( 'x', cn, sep = ':' ) } ), paste, collapse = "." )
    diffs.slope.seplam.DoL <- pairdiff.mat %*% mle.canton.slope.seplam.DoL$beta.hat[2:6]
    rownames(diffs.slope.seplam.DoL) <- sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ),
                                                  function(cn){ paste( 'x', cn, sep = ':' ) } ), paste, collapse = "." )
    se.diffs.slope.seplam.DoL <- sqrt(diag(pairdiff.mat %*% covML.canton.slope.seplam.DoL[2:6,2:6] %*% t(pairdiff.mat)))
    names(se.diffs.slope.seplam.DoL) <- sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ),
                                                  function(cn){ paste( 'x', cn, sep = ':' ) } ), paste, collapse = "." )

  # Overview tables !!! p-values with caution !!!
    cbind( estim = c(diffs.slope.DoL), se = se.diffs.slope.DoL, 
           pval = 2*( 1 - pnorm( abs(c(diffs.slope.DoL))/se.diffs.slope.DoL ) ) )
    # VS-OW and VS-TI largest; UR-TI is also significant but only a bit more than half as large as the two others; 
    #   UR-OW is of similar value as UR-TI but not significant
  
    cbind( estim = c(diffs.slope.seplam.DoL), se = se.diffs.slope.seplam.DoL, 
           pval = 2*( 1 - pnorm( abs(c(diffs.slope.seplam.DoL))/se.diffs.slope.seplam.DoL ) ) )
    # VS-OW and VS-TI largest; UR-TI and UR.OW only a bit more than half as large as the two others but both not signif.
  
    # Both (only relevant for PTBS model with single lambda):
      diffs.both.DoL <- ( diag(2) %x% pairdiff.mat ) %*% mle.canton.both.DoL$beta.hat
      rownames(diffs.both.DoL) <- c( rownames(pairdiff.mat), 
                                     sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ), 
                                                     function(cn){paste( 'x', cn, sep = ':' )}), paste, collapse = "."))
      se.diffs.both.DoL <- sqrt( diag( (diag(2) %x% pairdiff.mat) %*% covML.canton.both.DoL[1:10,1:10] %*% 
                                         t(diag(2) %x% pairdiff.mat) ) )
      names(se.diffs.both.DoL) <- c( rownames(pairdiff.mat), 
                                     sapply( lapply( strsplit( rownames(pairdiff.mat), ".", fixed = TRUE ), 
                                                     function(cn){paste( 'x', cn, sep = ':' )}), paste, collapse = "."))
      # HOWEVER: Cannot really compare lines with different slopes AND intercepts!!


#======================================================================================================================#
#===== Which differences are really significant?? (use a posteriori tests that are adjusted for multiple testing) =====#
#======================================================================================================================#

## Possible methods
# ( Fisher LSD not suitable for a posteriori tests )
# (i) Fixed range test (LSR) = Tukey HSD = Range-STP: se(single mean)*Q_{k,sum(n_i)-k} (increase of type II error); 
      # not exact for ranges of >=3 means with unequal sample size!! JMH: use 1/n.H = mean(1/n_i, i involved)
# ( Scheff: se(diff)*sqrt( (k-1) F_{k-1,N-k} ) ) too conservative
# (ii) Shortest siginificant range test (SSR) = Student-Newman-Keuls: ordered coefs, se(single mean)*Q_{k:2,sum(n_i)-k}
# (iii) SS-STP (Sum of Squares simult. test procedure), SS > (k-1) * MS.error * F_{k-1,N-k} (always exact, also for 
      # unequal sample sizes)
# Duncan: ordered, test at level alpha_p = 1 - (1-alpha)^{p-1} for p-2 separating means, 
      # critical value r_{p,nu,alpha} = Q_{p,nu,(1-alpha)^{p-1}} if p=2, r_{p,nu,alpha} = max( Q_{p,nu,(1-alpha)^{p-1}},
      # r_{p-1,nu,alpha} ) --> QUITE FAULTY ON WIKIPEDIA!!!

# (i) Tukey HSD = Range-STP = LSR: se(single mean)*Q_{k,sum(n_i)-k}; not exact for ranges of >=3 means with unequal 
      # sample size!! JMH: use 1/n.H = mean(1/n_i, i involved)
# (ii) Student-Newman-Keuls = SSR: ordered coefs, se(single mean)*Q_{k:2,sum(n_i)-k}; also not exact with unequal 
      # sample sizes!!

#(i) Tukey HSD:
  # Single critical value based on artificial sample size 1/mean(1/n_i):
    sqrt( canton.intcpt.DoL$s2 * mean(1/n.obs.aov.cant) ) * qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - 
                                                                      canton.intcpt.DoL$p )
    abs(diffs.intcpt.DoL) > sqrt( canton.intcpt.DoL$s2 * mean(1/n.obs.aov.cant) ) * 
          qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
    # VS.OW and VS.TI are significant, others not
  
  # "Worst case scenario" for artificial sample size 1/mean(1/range(n_i))
    sqrt( canton.intcpt.DoL$s2 * mean(1/range(n.obs.aov.cant)) ) * qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - 
                                                                             canton.intcpt.DoL$p )
    abs(diffs.intcpt.DoL) > sqrt( canton.intcpt.DoL$s2 * mean(1/range(n.obs.aov.cant)) ) * 
      qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
    # Nothing significant!
  
  # Critical values for pairs w/ different sample sizes (Tukey-Kramer method according to Wikipedia):
    qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * sqrt(canton.intcpt.DoL$s2) * 
      sqrt( combn( 5, 2, FUN = function(pair){ mean(1/n.obs.aov.cant[pair]) } ) )
    abs(c(diffs.intcpt.DoL)) > qtukey( 0.95, nmeans = 5, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * 
      sqrt(canton.intcpt.DoL$s2) * sqrt( combn( 5, 2, FUN = function(pair){ mean(1/n.obs.aov.cant[pair]) } ) )
    # Nothing significant!
  
  # Scheff confidence intervals (which are more conservative than Tukey):
    cbind( c(diffs.intcpt.DoL) - sqrt( 4 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) ) * 
             sqrt( canton.intcpt.DoL$s2 * combn( 5, 2, FUN = function(pair){ sum(1/n.obs.aov.cant[pair]) } ) ), 
           c(diffs.intcpt.DoL) + sqrt( 4 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) ) * 
             sqrt( canton.intcpt.DoL$s2 * combn( 5, 2, FUN = function(pair){ sum(1/n.obs.aov.cant[pair]) } ) ) )
    # All CIs contain zero, thus no significant differences!

# (ii) Shortest significant range test:
  # In order of decreasing abs(diffs.intcpt.DoL)
  # Critical values (not ordered):
    combn( 5, 2, FUN = function(pair){ qtukey( 0.95, nmeans = diff(pair)+1, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * 
        sqrt( mean(1/n.obs.aov.cant[pair]) ) } ) * sqrt(canton.intcpt.DoL$s2)
  
  # Test for ordered coefficients
    combn( 5, 2, function(pair){ diff( sort(canton.intcpt.DoL$coefs[1:5])[pair] ) } ) > 
      combn( 5, 2, FUN = function(pair){ qtukey( 0.95, nmeans = diff(pair)+1, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) *
          sqrt( mean(1/n.obs.aov.cant[order(canton.intcpt.DoL$coefs[1:5])][pair]) ) } ) * sqrt(canton.intcpt.DoL$s2)
    # Largest difference (4th entry) is not significant!! --> don't go further.

# (iii) SS-STP: implementation probably not meaningful!
  # Critical value:
    4 * canton.intcpt.DoL$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
  # SS for all pairs:
    combn( 5, 2, FUN = function(pair){ sum( (tapply( y.lambda.aov.DoL - canton.intcpt.DoL$coefs[6] * x.lambda.aov.DoL, 
                                                     region.aov, sum )^2/n.obs.aov.cant)[pair] ) - 
        sum( tapply( y.lambda.aov.DoL, region.aov, sum )[pair] )^2/sum(n.obs.aov.cant[pair]) } ) # this is probably wrong!!!
    # TI.SZ, UR.SZ, VS.SZ not significant, all others well (?!
    
    anova( lm( y.lambda.aov.DoL ~ x.lambda.aov.DoL ) )
    anova( lm( y.lambda.aov.DoL ~ x.lambda.aov.DoL + region.aov ) )
    (n.obs.aov-1)*var(y.lambda.aov.DoL) # SS.tot
    
# I need SS.tot = RSS(canton.no.DoL) = 91.527
# SS.tot(intcpt) = SS.mod(intcpt) + SS.error(intcpt) = SS.x + SS.region + SS.error(intcpt)
# SS.tot(regr) = SS.mod(regr) + SS.error(regr) = SS.x + SS.error(regr)
# SS.error(regr) = SS.region + SS.error(intcpt)
# SS.mod(intcpt) - SS.mod(regr)
  sum( n.obs.aov.cant * tapply( y.lambda.aov.DoL, region.aov, mean )^2 ) - n.obs.aov * mean(y.lambda.aov.DoL)^2 + 
    ( sum(x.lambda.aov.DoL*y.lambda.aov.DoL) - sum( n.obs.aov.cant * tapply( x.lambda.aov.DoL, region.aov, mean ) * 
                                                      tapply( y.lambda.aov.DoL, region.aov, mean ) ) )^2/
    ( sum(x.lambda.aov.DoL^2) - sum(n.obs.aov.cant * tapply( x.lambda.aov.DoL, region.aov, mean )^2) ) - 
    ( sum(x.lambda.aov.DoL*y.lambda.aov.DoL) - n.obs.aov * mean(x.lambda.aov.DoL) * mean(y.lambda.aov.DoL) )^2/
    ( sum(x.lambda.aov.DoL^2) - n.obs.aov * mean(x.lambda.aov.DoL)^2 )
  # This sum cannot be decomposed into k components due to the presence of the covariate x!!!

# Thus try SS-STP test with residuals from regression with single intercept.
# Critical value:
  4 * canton.intcpt.DoL$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
# SS for all pairs:
  combn( 5, 2, FUN = function(pair){ sum( (tapply( canton.no.DoL$e, region.aov, sum )^2/n.obs.aov.cant)[pair] ) - 
      sum( tapply( canton.no.DoL$e, region.aov, sum )[pair] )^2/sum(n.obs.aov.cant[pair]) } )
# Nothing significant

SS.STP.m <- function( group.list, y.i.bar, n.i ) {
	sum( sapply( group.list, function(grp){ sum( (n.i*y.i.bar)[grp] )^2/sum(n.i[grp]) } ) ) - 
    sum(( n.i * y.i.bar )[unlist(group.list)])^2/sum(n.i[unlist(group.list)])
}

mu.i.hat.srt <- sort(tapply( canton.no.DoL$e, region.aov, mean ))
n.i.srt <- n.obs.aov.cant[order(tapply( canton.no.DoL$e, region.aov, mean ))]

# Critical value:
  4 * my.lm( canton.no.DoL$e, Region.mat.aov )$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p )

# cantons are different:
  SS.STP.m( 1:5, mu.i.hat.srt, n.i.srt ) # not all equal: yes!

# Three groups: {OW,TI}, {SZ,UR}, {VS}
  SS.STP.m( list( 1, 2:3, 4:5 ), mu.i.hat.srt, n.i.srt ) # {VS},{UR,SZ},{TI,OW}: yes!


H0.mat <- rbind( c(1,0,-1,0,0), c(0,1,0,-1,0) )

# sum(canton.no.DoL$e^2) - sum( ( canton.no.DoL$e - Region.mat.aov %*% ( tapply( canton.no.DoL$e, region.aov, mean ) -
#   c( solve(crossprod(Region.mat.aov)) %*% t(H0.mat) %*% solve( H0.mat %*% solve(crossprod(Region.mat.aov)) %*% t(H0.mat) ) %*%
       # H0.mat %*% tapply( canton.no.DoL$e, region.aov, mean ) ) ) )^2 )
# SS.mod.corr(H0) = SS.tot.corr - SS.err(H0) = SS.err(all equal) - SS.err(H_0: OW=TI,SZ=UR)

# t( H0.mat %*% tapply( canton.no.DoL$e, region.aov, mean ) ) %*% solve( H0.mat %*% solve(crossprod(Region.mat.aov)) %*%
  # t(H0.mat) ) %*% H0.mat %*% tapply( canton.no.DoL$e, region.aov, mean )

# {VS} is different from {SZ,UR}
SS.STP.m( list( c(2,4), 5 ), mu.i.hat.srt, n.i.srt ) # No

# {VS} is different from {OW,TI}
SS.STP.m( list( c(1,3), 5 ), mu.i.hat.srt, n.i.srt ) # No

# {VS} is different from the other 4
SS.STP.m( list( 1:4, 5 ), mu.i.hat.srt, n.i.srt ) # No

# {OW,TI} different from {SZ,UR,VS}
SS.STP.m( list( c(1,3), c(2,4,5) ), mu.i.hat.srt, n.i.srt ) # Just not.

# {OW,TI} different from {SZ,UR}
SS.STP.m( list( c(1,3), c(2,4) ), mu.i.hat.srt, n.i.srt ) # No

# Pairwise differences:
combn( 5, 2, FUN = function(pair){ SS.STP.m( pair, mu.i.hat.srt, n.i.srt ) } ) # None significant

# Again DoL for residuals:
4 * my.lm( canton.no.DoL$e, Region.mat.aov )$s2 *(canton.intcpt.DoL$n - (canton.intcpt.DoL$p - 1))/
  (canton.intcpt.DoL$n - canton.intcpt.DoL$p) * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p )
    # = 2.291624; Correction in s2 to account for the fact the the residuals have df = n-1.

SS.STP.m( 1:5, mu.i.hat.srt, n.i.srt ) # not all equal: yes!
SS.STP.m( list( 1, 5 ), mu.i.hat.srt, n.i.srt ) # VS-OW: no
SS.STP.m( list( 1, 2:5 ), mu.i.hat.srt, n.i.srt ) # {VS},{UR,SZ,TI,OW}: no
SS.STP.m( list( 1:2, 3:5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR},{SZ,TI,OW}: no
SS.STP.m( list( 1:3, 4:5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR,SZ},{TI,OW}: Just not (p-val = 0.060).
SS.STP.m( list( 1:4, 5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR,SZ,TI},{OW}: no
SS.STP.m( list( 2:3, 4:5), mu.i.hat.srt, n.i.srt ) # {UR,SZ},{TI,OW}: no
SS.STP.m( list( 1, 2:3, 4:5 ), mu.i.hat.srt, n.i.srt ) # {VS},{UR,SZ},{TI,OW}: yes!

SS.STP.m( list( 2, 5 ), mu.i.hat.srt, n.i.srt ) # UR-OW: no
SS.STP.m( list( 3, 5 ), mu.i.hat.srt, n.i.srt ) # SZ-OW: no
SS.STP.m( list( 4, 5 ), mu.i.hat.srt, n.i.srt ) # TI-OW: no
SS.STP.m( list( 1, 4 ), mu.i.hat.srt, n.i.srt ) # VS-TI: no

1 - pf( SS.STP.m( list( 1:3, 4:5 ), mu.i.hat.srt, n.i.srt ), df1 = 4, df2 = canton.intcpt.DoL$n - 1 - 2 )

library(partitions)
test.SSm.partitions.r <- sapply( apply( setparts(5), 2, function(setidx){ lapply( 1:max(setidx), function(i){
  which(setidx == i) } ) } ), function(grl){ SS.STP.m( grl, mu.i.hat.srt, n.i.srt ) } )
# Critical values with adapted error sum in denominator:
crit.SSm.partitions.r <- sapply( apply( setparts(5), 2, function(setidx){ lapply( 1:max(setidx), function(i){
  which(setidx == i) } ) } ), function(grl){ 4 * ( ( sum( canton.no.DoL$e^2 ) - sum( sapply( grl, function(grp){
    sum( c(n.i.srt*mu.i.hat.srt)[grp] )^2/sum(n.i.srt[grp]) } ) ) )/(n.obs.aov - 1 - length(grl)) ) * 
      qf( 0.95, df1 = 4, df2 = n.obs.aov - 1 - length(grl) ) } ) 
# this seems correct (at least for the computation, for the statistics I don't know...!)

( sum( canton.no.DoL$e^2 ) - sum( sapply( list(1,2,3,4,5), function(grp){ sum( c(n.i.srt*mu.i.hat.srt)[grp] )^2/
    sum(n.i.srt[grp]) } ) ) )/(n.obs.aov - 1 - length(list(1,2,3,4,5)))

setparts(5)[,which( test.SSm.partitions.r > 4 * my.lm( canton.no.DoL$e, Region.mat.aov )$s2 * 
                      qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) )]

# Rejected (i.e., reasonable separation) partitions are: 
  # {VS,UR},{SZ},{TI,OW}
  # {VS,SZ},{UR},{TI,OW}
  # {VS},{UR,SZ},{TI,OW} <-- most reasonable grouping
  # {VS,UR},{SZ},{TI},{OW}
  # {VS,SZ},{UR},{TI},{OW}
  # {VS},{UR,SZ},{TI},{OW}
  # {VS},{UR},{SZ},{TI,OW}
  # {VS},{UR},{SZ},{TI},{OW} (All five different)
  setparts(5)[,which( test.SSm.partitions.r > crit.SSm.partitions.r )]
  # same result as before, thus non-adjustment of MSE does not matter


## SS-STP for model with intercept by canton
mu.i.hat.srt <- sort(canton.intcpt.DoL$coefs[1:5])
n.i.srt <- n.obs.aov.cant[order(canton.intcpt.DoL$coefs[1:5])]

# Critical value
4 * canton.intcpt.DoL$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) # = 2.2911

SS.STP.m( 1:5, mu.i.hat.srt, n.i.srt ) # not all equal: yes!
SS.STP.m( list( 1, 5 ), mu.i.hat.srt, n.i.srt ) # VS-OW: no
SS.STP.m( list( 1, 2:5 ), mu.i.hat.srt, n.i.srt ) # {VS},{UR,SZ,TI,OW}: no
SS.STP.m( list( 1:2, 3:5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR},{SZ,TI,OW}: no
SS.STP.m( list( 1:3, 4:5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR,SZ},{TI,OW}: yes (just)
SS.STP.m( list( 1:4, 5 ), mu.i.hat.srt, n.i.srt ) # {VS,UR,SZ,TI},{OW}: no
SS.STP.m( list( 2:3, 4:5), mu.i.hat.srt, n.i.srt ) # {UR,SZ},{TI,OW}: no
SS.STP.m( list( 1, 2:3, 4:5 ), mu.i.hat.srt, n.i.srt ) # {VS},{UR,SZ},{TI,OW}: yes!

SS.STP.m( list( 2, 5 ), mu.i.hat.srt, n.i.srt ) # UR-OW: no
SS.STP.m( list( 3, 5 ), mu.i.hat.srt, n.i.srt ) # SZ-OW: no
SS.STP.m( list( 4, 5 ), mu.i.hat.srt, n.i.srt ) # TI-OW: no
SS.STP.m( list( 1, 4 ), mu.i.hat.srt, n.i.srt ) # VS-TI: no

test.SSm.partitions <- sapply( apply( setparts(5), 2, function(setidx){ lapply( 1:max(setidx), function(i){
  which(setidx == i) } ) } ), function(grl){ SS.STP.m( grl, mu.i.hat.srt, n.i.srt ) } )
# Critical values with adapted error sum in denominator:
crit.SSm.partitions <- sapply( apply( setparts(5), 2, function(setidx){ lapply( 1:max(setidx), function(i){
  which(setidx == i) } ) } ), function(grl){ 4 * ( ( sum( ( y.lambda.aov.DoL - canton.intcpt.DoL$coefs[6] * 
                                                              x.lambda.aov.DoL )^2 ) - sum( sapply( grl, function(grp){
                                                                sum( c(n.i.srt*mu.i.hat.srt)[grp] )^2/sum(n.i.srt[grp]) } ) ) )/
                                                     (n.obs.aov - length(grl) - 1) ) * qf( 0.95, df1 = 4, df2 = n.obs.aov - 
                                                                                             length(grl) - 1 ) } ) 
# this seems correct (for computation, statistically I have no idea...)

sum( sapply( list(1,2,3,4,5 ), function(grp){ sum( c(n.i.srt*mu.i.hat.srt)[grp] )^2/sum(n.i.srt[grp]) } ) )
 
# SS.error(H_0: mu12 = mu3) - SS.error(mu12,mu3) = SS.error(H_0: mu1 = mu2 = mu3) - SS.error(H_0: mu1 = mu2)

setparts(5)[,which( test.SSm.partitions > 4 * canton.intcpt.DoL$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - 
                                                                           canton.intcpt.DoL$p ) )] # 10 partitions
# Rejected (i.e., reasonable separation) partitions are: 
  # {VS,UR,SZ},{TI,OW}
  # {VS,UR,SZ},{TI},{OW}
  # {VS,UR},{SZ},{TI,OW}
  # {VS,SZ},{UR},{TI,OW}
  # {VS},{UR,SZ},{TI,OW} <-- most reasonable grouping
  # {VS,UR},{SZ},{TI},{OW}
  # {VS,SZ},{UR},{TI},{OW}
  # {VS},{UR,SZ},{TI},{OW}
  # {VS},{UR},{SZ},{TI,OW}
  # {VS},{UR},{SZ},{TI},{OW} (All five different)
  setparts(5)[,which( test.SSm.partitions > crit.SSm.partitions )] # same result as before --> non-adjustment of critical value does not matter
  
  test.SSm.partitions[which( test.SSm.partitions > 4 * canton.intcpt.DoL$s2 * qf( 0.95, df1 = 4, df2 = canton.intcpt.DoL$n - 
                                                                                    canton.intcpt.DoL$p ) )]
  my.lm( canton.no.DoL$e, Region.mat.aov )$se.coefs

# SS.mod.corr(intcpt)
sum( n.obs.aov.cant * tapply( y.lambda.aov.DoL, region.aov, mean )^2 ) + 
  ( sum(x.lambda.aov.DoL*y.lambda.aov.DoL) - sum( n.obs.aov.cant * tapply( x.lambda.aov.DoL, region.aov, mean ) * 
                                                    tapply( y.lambda.aov.DoL, region.aov, mean ) ) )^2/
  ( sum(x.lambda.aov.DoL^2) - sum(n.obs.aov.cant * tapply( x.lambda.aov.DoL, region.aov, mean )^2) ) - 
  n.obs.aov * mean(y.lambda.aov.DoL)^2
# SS.mod.corr(regr)
( sum(x.lambda.aov.DoL*y.lambda.aov.DoL) - n.obs.aov * mean(x.lambda.aov.DoL) * mean(y.lambda.aov.DoL) )^2/
  ( sum(x.lambda.aov.DoL^2) - n.obs.aov * mean(x.lambda.aov.DoL)^2 )

summary( lm( canton.no.DoL$e ~ region.aov ) )
anova( lm( canton.no.DoL$e ~ region.aov ) )

cbind( pairdiff.mat %*% my.lm( canton.no.DoL$e, Region.mat.aov )$coefs, diffs.intcpt.DoL ) 
# differences for anova of residuals are very similar but weaker (except for UR-SZ and UR-TI)

# Newman-Keuls for residuals
  # Sorted means:
    sort(tapply( canton.no.DoL$e, region.aov, mean ))
  # Pooled variance:
    sum( (n.obs.aov.cant-1) * tapply( canton.no.DoL$e, region.aov, var ) )/sum(n.obs.aov.cant-1)
  
  # Critical value pattern for "pair" (!! diff(pair) must be in sorted case !!)
    # sqrt( sum( (n.obs.aov.cant-1) * tapply( canton.no.DoL$e, region.aov, var ) )/sum(n.obs.aov.cant-1) ) * 
    #   qtukey( 0.95, diff(pair)+1, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * sqrt( mean(1/n.obs.aov.cant[pair]) ) 
  
  # Critical value for range of 5 means:
    sqrt( sum( (n.obs.aov.cant-1) * tapply( canton.no.DoL$e, region.aov, var ) )/sum(n.obs.aov.cant-1) ) * 
      qtukey( 0.95, 5, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * sqrt( mean(1/n.obs.aov.cant[c(1,5)]) )
  # OW-VS (largest difference) 
    sort( abs( pairdiff.mat %*% tapply( canton.no.DoL$e, region.aov, mean ) ), decreasing = TRUE )[1] 
    # This is not significant!!! Thus no more tests...
  
  # Nevertheless test of all differences (due to unequal sample sizes):
    combn( 5, 2, function(pair){ diff( sort(tapply( canton.no.DoL$e, region.aov, mean ))[pair] ) } ) > 
      combn( 5, 2, function(pair){ sqrt( sum( (n.obs.aov.cant-1) * tapply( canton.no.DoL$e, region.aov, var ) )/
                                           sum(n.obs.aov.cant-1) ) * 
          qtukey( 0.95, diff(pair)+1, df = canton.intcpt.DoL$n - canton.intcpt.DoL$p ) * 
          sqrt( mean(1/n.obs.aov.cant[order(tapply( canton.no.DoL$e, region.aov, mean ))][pair]) ) } )
  
  # Test of homogeneity of variances:
    1 - pchisq( ( sum(n.obs.aov.cant-1) * log( sum( (n.obs.aov.cant-1) * tapply( canton.no.DoL$e, region.aov, var ) )/
                                                 sum(n.obs.aov.cant-1) ) - 
                    sum( (n.obs.aov.cant-1) * log(tapply( canton.no.DoL$e, region.aov, var )) ) )/
                  ( 1 + ( sum( 1/(n.obs.aov.cant-1) ) - 1/sum(n.obs.aov.cant-1) )/(3*4) ), df = 4 )
  # Variances not equal...
  
  # Check normality of samples (because previous test is sensitive to lack of normality, in which case it might also reject the homogeneity of variances):
    plot( qnorm( (1:n.obs.aov.cant[1])/(n.obs.aov.cant[1] + 1) ), 
          sort( ( canton.no.DoL$e[as.logical(Region.mat.aov[,1])] - 
                    mean(canton.no.DoL$e[as.logical(Region.mat.aov[,1])]) )/
                  sd(canton.no.DoL$e[as.logical(Region.mat.aov[,1])]) ), 
          ylab = "Observed standardized", xlab = "Theoretical N(0,1) quantiles" )
    abline( 0, 1, col = 'red' )
    
    plot( qnorm( (1:n.obs.aov.cant[2])/(n.obs.aov.cant[2] + 1) ), 
          sort( ( canton.no.DoL$e[as.logical(Region.mat.aov[,2])] - 
                    mean(canton.no.DoL$e[as.logical(Region.mat.aov[,2])]) )/
                  sd(canton.no.DoL$e[as.logical(Region.mat.aov[,2])]) ),
          ylab = "Observed standardized", xlab = "Theoretical N(0,1) quantiles" )
    abline( 0, 1, col = 'red' )
    
    plot( qnorm( (1:n.obs.aov.cant[3])/(n.obs.aov.cant[3] + 1) ), 
          sort( ( canton.no.DoL$e[as.logical(Region.mat.aov[,3])] -
                    mean(canton.no.DoL$e[as.logical(Region.mat.aov[,3])]) )/
                  sd(canton.no.DoL$e[as.logical(Region.mat.aov[,3])]) ), 
          ylab = "Observed standardized", xlab = "Theoretical N(0,1) quantiles" )
    abline( 0, 1, col = 'red' )
    
    plot( qnorm( (1:n.obs.aov.cant[4])/(n.obs.aov.cant[4] + 1) ), 
          sort( ( canton.no.DoL$e[as.logical(Region.mat.aov[,4])] - 
                    mean(canton.no.DoL$e[as.logical(Region.mat.aov[,4])]) )/
                  sd(canton.no.DoL$e[as.logical(Region.mat.aov[,4])]) ), 
          ylab = "Observed standardized", xlab = "Theoretical N(0,1) quantiles" )
    abline( 0, 1, col = 'red' )
    
    plot( qnorm( (1:n.obs.aov.cant[5])/(n.obs.aov.cant[5] + 1) ), 
          sort( ( canton.no.DoL$e[as.logical(Region.mat.aov[,5])] - 
                    mean(canton.no.DoL$e[as.logical(Region.mat.aov[,5])]) )/
                  sd(canton.no.DoL$e[as.logical(Region.mat.aov[,5])]) ), 
          ylab = "Observed standardized", xlab = "Theoretical N(0,1) quantiles" )
    abline( 0, 1, col = 'red' )
    # Last two (UR and VS) slightly questionable but not terribly bad.
  
  # F_max test for equal variances
  library(SuppDists)
  1 - pmaxFratio( max(tapply( canton.no.DoL$e, region.aov, var ))/min(tapply( canton.no.DoL$e, region.aov, var )), 
                  df = 18, k = 5 ) # S&R suggest to use lesser df
  # Not significant (with larger df = 76 it would however, and with mean df 47 as well)!


# Box 13.2 Sokal and Rohlf 1961 (approximate test of equality of means when variances are heterogeneous)
  w <- n.obs.aov.cant/tapply( canton.no.DoL$e, region.aov, var )
  CT.w <- sum( w * tapply( canton.no.DoL$e, region.aov, mean ) )^2/sum(w)
  1 - pf( ( (sum( w * tapply( canton.no.DoL$e, region.aov, mean )^2 ) - CT.w)/4 )/( 1 + (2*3)/(5^2-1) * sum( (1 - w/sum(w))^2/(n.obs.aov.cant-1) ) ), df1 = 4, df2 = (5^2-1)/(3*sum( (1 - w/sum(w))^2/(n.obs.aov.cant-1) )) ) # p-val = 0.0146 --> Reject equality of means
  # Checked with data from book that correct

# Kruskal-Wallis test
  rank(canton.no.DoL$e) 
  min( diff(sort(canton.no.DoL$e)) ) # no ties, o/w need correction factor which I haven't implemented
  tapply( rank(canton.no.DoL$e), region.aov, sum )
  1 - pchisq( 12/( n.obs.aov * (n.obs.aov + 1) ) * sum( tapply( rank(canton.no.DoL$e), region.aov, sum )^2/
                                                          n.obs.aov.cant ) - 3*(n.obs.aov + 1), df = 4 ) # not all equal
  # Checked with book example

# Dunn test (?)
# Critical value:
  combn( 5, 2, function(pair){ diff( sort(tapply( rank(canton.no.DoL$e), region.aov, mean ))[pair] ) } ) > combn( 5, 2, function(pair){ 2.807 * sqrt( n.obs.aov * (n.obs.aov + 1)/12 * mean(1/n.obs.aov.cant[order(tapply( rank(canton.no.DoL$e), region.aov, mean ))][pair]) ) } )
  sort(tapply( rank(canton.no.DoL$e), region.aov, mean ))
  # VS-TI, VS-OW, and SZ-OW seem to be significant -- but no adjustment for multiple comparisons, thus actually they 
  #   might not be significant!!
  library(dunn.test)
  dunn.test( canton.no.DoL$e, g = region.aov, altp = TRUE, method = "none" )
  dunn.test( canton.no.DoL$e, g = region.aov, altp = TRUE, method = "by"  ) # with any adjustment method, nothing is significant anymore...


# For combined model not even worth trying because the artificial overall estimate of sigma2 is much larger than previously!
# Use an artificial MS.error based on my estimates from the combined model (suppose se(beta.hat[i]) is se(y.bar[i]):
  sum( (n.obs.aov.cant-1) * n.obs.aov.cant * canton.intcpt.DoL$se.coefs[1:5]^2 )/sum(n.obs.aov.cant-1)
  
  tapply( canton.no.DoL$e, region.aov, sd )
  sqrt(n.obs.aov.cant) * canton.intcpt.DoL$se.coefs[1:5]
  
  canton.intcpt.DoL$s2
  sum( (n.obs.aov.cant-1) * n.obs.aov.cant * canton.intcpt.DoL$se.coefs[1:5]^2 )/sum(n.obs.aov.cant-1)
  sum( (n.obs.aov.cant-1) * n.obs.aov.cant * my.lm( canton.no.DoL$e, Region.mat.aov )$se.coefs^2 )/sum(n.obs.aov.cant-1)
  
  cbind( canton.intcpt.DoL$se.coefs[1:5], sqrt(canton.intcpt.DoL$s2/n.obs.aov.cant) )


##=======================================##
##==== Absolute loss ====================##
##=======================================##

# Fit the three models
  mle.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = FALSE )
  mle.seplam.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = TRUE )
  mle.other.Loss <- TBS.mle( y = y.Loss, x = x.Loss, beta.init = c(1000,0.5), interval = c(-3,1) )

# LRT sep.lam vs. single lambda
  1 - pchisq( -2*( mle.Loss$ell - mle.seplam.Loss$ell ), df = 1 )
# AIC:
  -2*mle.Loss$ell + 2*length(unlist(mle.Loss[1:3])) # 6841.861
  -2*mle.seplam.Loss$ell + 2*length(unlist(mle.seplam.Loss[1:3])) # 6825.102
  -2*mle.other.Loss$ell + 2*length(unlist(mle.other.Loss[1:3])) # 6821.25
# BIC:
  -2*mle.Loss$ell + log(n.obs)*length(unlist(mle.Loss[1:3])) # 6857.653
  -2*mle.seplam.Loss$ell + log(n.obs)*length(unlist(mle.seplam.Loss[1:3])) # 6844.842
  -2*mle.other.Loss$ell + log(n.obs)*length(unlist(mle.other.Loss[1:3])) # 6837.042


# Redefine quantities (including xy.lambda) with appropriate observations omitted
# Original scale:
  y.aov.Loss <- Contents$Loss[-idx.omit]
  x.aov.Loss <- Structure$Loss[-idx.omit]
  
  X.canton.aov.Loss <- cbind( Region.mat.aov, Region.mat.aov * x.aov.Loss )
  colnames(X.canton.aov.Loss)[6:10] <- paste( 'x', colnames(X.canton.aov.Loss)[1:5], sep = ':' )

#(A) Take lambda as fixed:
  # (A.1) PTBS model (need to redefine xy.lambda first --> execute whole block)
    # Transformed scale (lambda fixed to its fit over the whole data set)
      y.lambda.aov.Loss <- BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[-idx.omit] )
      x.lambda.aov.Loss <- BC.transform( mle.Loss$lambda.hat[1+0], Structure$Loss[-idx.omit] )

      # X.lambda.canton.Loss <- cbind( Region.mat, Region.mat * x.lambda.Loss )
      X.lambda.canton.aov.Loss <- cbind( Region.mat.aov, Region.mat.aov * x.lambda.aov.Loss )
      colnames(X.lambda.canton.aov.Loss)[6:10] <- paste( 'x.lam', colnames(X.lambda.canton.aov.Loss)[1:5], sep = ':' )

    # Intercept and slope by canton
      canton.both.Loss <- my.lm( y.lambda.aov.Loss, X.lambda.canton.aov.Loss )
      names( canton.both.Loss$coefs ) <- colnames(X.lambda.canton.aov.Loss)
    # Only intercept by canton
      canton.intcpt.Loss <- my.lm( y.lambda.aov.Loss, cbind( Region.mat.aov, x.lam = x.lambda.aov.Loss ) )
      names( canton.intcpt.Loss$coefs ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      canton.slope.Loss <- my.lm( y.lambda.aov.Loss, cbind( intcpt = 1, X.lambda.canton.aov.Loss[,6:10] ) )
      names( canton.slope.Loss$coefs ) <- c( 'one', colnames(X.lambda.canton.aov.Loss[,6:10]) )
    # Nothing by canton
      canton.no.Loss <- my.lm( y.lambda.aov.Loss, cbind( intcpt = 1, x.lam = x.lambda.aov.Loss ) )
      names( canton.no.Loss$coefs ) <- c('one','x')

# (A.2) PTBS model with seplam (need to redefine xy.lambda first --> execute whole block)
    # Transformed scale (lambda fixed to its fit over the whole data set)
      y.lambda.aov.Loss <- BC.transform( mle.seplam.Loss$lambda.hat[1], Contents$Loss[-idx.omit] )
      x.lambda.aov.Loss <- BC.transform( mle.seplam.Loss$lambda.hat[1+1], Structure$Loss[-idx.omit] )

    #X.lambda.canton.Loss <- cbind( Region.mat, Region.mat * x.lambda.Loss )
      X.lambda.canton.aov.Loss <- cbind( Region.mat.aov, Region.mat.aov * x.lambda.aov.Loss )
      colnames(X.lambda.canton.aov.Loss)[6:10] <- paste( 'x.lam', colnames(X.lambda.canton.aov.Loss)[1:5], sep = ':' )
    
    # Intercept and slope by canton
      canton.both.seplam.Loss <- my.lm( y.lambda.aov.Loss, X.lambda.canton.aov.Loss )
      names( canton.both.seplam.Loss$coefs ) <- colnames(X.lambda.canton.aov.Loss)
    # Only intercept by canton
      canton.intcpt.seplam.Loss <- my.lm( y.lambda.aov.Loss, cbind( Region.mat.aov, x.lam = x.lambda.aov.Loss ) )
      names( canton.intcpt.seplam.Loss$coefs ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      canton.slope.seplam.Loss <- my.lm( y.lambda.aov.Loss, cbind( intcpt = 1, X.lambda.canton.aov.Loss[,6:10] ) )
      names( canton.slope.seplam.Loss$coefs ) <- c( 'one', colnames(X.lambda.canton.aov.Loss[,6:10]) )
    # Nothing by canton
      canton.no.seplam.Loss <- my.lm( y.lambda.aov.Loss, cbind( intcpt = 1, x.lam = x.lambda.aov.Loss ) )
      names( canton.no.seplam.Loss$coefs ) <- c('one','x')

# (A.3) TBS model (also lambda fixed; model is on original scale)
  # Nothing by canton
    canton.no.other.Loss <- TBS.mle( y = y.aov.Loss, x = x.aov.Loss, lambda = mle.other.Loss$lambda.hat, 
                                     beta.init = c(1000,0.5), intercept = TRUE )
    names( canton.no.other.Loss$beta.hat ) <- c('one','x')
  # Only intercept by canton
    canton.intcpt.other.Loss <- TBS.mle( y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), 
                                         lambda = mle.other.Loss$lambda.hat, beta.init = rep( canton.no.other.Loss$beta.hat, 
                                                                                              c(5,1) ), intercept = FALSE )
    names( canton.intcpt.other.Loss$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
  # Only slope by canton
    canton.slope.other.Loss <- TBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss[,6:10], lambda = mle.other.Loss$lambda.hat, 
                                        beta.init = rep( canton.no.other.Loss$beta.hat, c(1,5) ), intercept = TRUE )
    names( canton.slope.other.Loss$beta.hat ) <- c( 'one', colnames(X.canton.aov.Loss[,6:10]) )
  # Intercept and slope by canton
    canton.both.other.Loss <- TBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss, lambda = mle.other.Loss$lambda.hat, 
                                       beta.init = c(canton.intcpt.other.Loss$beta.hat[1:5], 
                                                     canton.slope.other.Loss$beta.hat[2:6]), intercept = FALSE )
    names( canton.both.other.Loss$beta.hat ) <- colnames(X.canton.aov.Loss)
 

#=== Comparisons (lambda fixed)
## For regression models use F-tests corresponding to anova( lm.simple, lm.complex )
  RSS.no.Loss <- sum(canton.no.Loss$e^2)
  RSS.intcpt.Loss <- sum(canton.intcpt.Loss$e^2)
  RSS.slope.Loss <- sum(canton.slope.Loss$e^2)
  RSS.both.Loss <- sum(canton.both.Loss$e^2)
  
  RSS.no.seplam.Loss <- sum(canton.no.seplam.Loss$e^2)
  RSS.intcpt.seplam.Loss <- sum(canton.intcpt.seplam.Loss$e^2)
  RSS.slope.seplam.Loss <- sum(canton.slope.seplam.Loss$e^2)
  RSS.both.seplam.Loss <- sum(canton.both.seplam.Loss$e^2)

#(a) both vs. nothing
  1 - pf( ( ( RSS.no.Loss - RSS.both.Loss )/(canton.both.Loss$p - canton.no.Loss$p) )/canton.both.Loss$s2, 
          df1 = canton.both.Loss$p - canton.no.Loss$p, df2 = canton.both.Loss$n - canton.both.Loss$p )
  1 - pf( ( ( RSS.no.seplam.Loss - RSS.both.seplam.Loss )/(canton.both.seplam.Loss$p - canton.no.seplam.Loss$p) )/
            canton.both.seplam.Loss$s2, df1 = canton.both.seplam.Loss$p - canton.no.seplam.Loss$p, 
          df2 = canton.both.seplam.Loss$n - canton.both.seplam.Loss$p )
  1 - pchisq( -2*( canton.no.other.Loss$ell - canton.both.other.Loss$ell ), df = ncol(X.canton.aov.Loss) - 2 )
  # not signif. (neither for seplam nor TBS)
#(b) both vs. intercept
  1 - pf( ( ( RSS.intcpt.Loss - RSS.both.Loss )/(canton.both.Loss$p - canton.intcpt.Loss$p) )/canton.both.Loss$s2, 
          df1 = canton.both.Loss$p - canton.intcpt.Loss$p, df2 = canton.both.Loss$n - canton.both.Loss$p )
  1 - pf( ( ( RSS.intcpt.seplam.Loss - RSS.both.seplam.Loss )/(canton.both.seplam.Loss$p - canton.intcpt.seplam.Loss$p) )/
            canton.both.seplam.Loss$s2, df1 = canton.both.seplam.Loss$p - canton.intcpt.seplam.Loss$p, 
          df2 = canton.both.seplam.Loss$n - canton.both.seplam.Loss$p )
  1 - pchisq( -2*( canton.intcpt.other.Loss$ell - canton.both.other.Loss$ell ), df = ncol(X.canton.aov.Loss) - 
                (ncol(Region.mat.aov)+1) )
  # not signif. (neither for seplam nor TBS)
#(c) both vs. slope
  1 - pf( ( ( RSS.slope.Loss - RSS.both.Loss )/(canton.both.Loss$p - canton.slope.Loss$p) )/canton.both.Loss$s2, 
          df1 = canton.both.Loss$p - canton.slope.Loss$p, df2 = canton.both.Loss$n - canton.both.Loss$p )
  1 - pf( ( ( RSS.slope.seplam.Loss - RSS.both.seplam.Loss )/(canton.both.seplam.Loss$p - canton.slope.seplam.Loss$p) )/
            canton.both.seplam.Loss$s2, df1 = canton.both.seplam.Loss$p - canton.slope.seplam.Loss$p, 
          df2 = canton.both.seplam.Loss$n - canton.both.seplam.Loss$p )
  1 - pchisq( -2*( canton.slope.other.Loss$ell - canton.both.other.Loss$ell ), df = ncol(X.canton.aov.Loss) - 
                (ncol(Region.mat.aov)+1) )
  # not signif. (neither for seplam nor TBS)
#(d) intercept vs. nothing 
  1 - pf( ( ( RSS.no.Loss - RSS.intcpt.Loss )/(canton.intcpt.Loss$p - canton.no.Loss$p) )/canton.intcpt.Loss$s2, 
          df1 = canton.intcpt.Loss$p - canton.no.Loss$p, df2 = canton.intcpt.Loss$n - canton.intcpt.Loss$p )
  1 - pf( ( ( RSS.no.seplam.Loss - RSS.intcpt.seplam.Loss )/(canton.intcpt.seplam.Loss$p - canton.no.seplam.Loss$p) )/
            canton.intcpt.seplam.Loss$s2, df1 = canton.intcpt.seplam.Loss$p - canton.no.seplam.Loss$p, 
          df2 = canton.intcpt.seplam.Loss$n - canton.intcpt.seplam.Loss$p )
  1 - pchisq( -2*( canton.no.other.Loss$ell - canton.intcpt.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
  # not signif. (neither for seplam nor TBS)
#(e) slope vs. nothing
  1 - pf( ( ( RSS.no.Loss - RSS.slope.Loss )/(canton.slope.Loss$p - canton.no.Loss$p) )/canton.slope.Loss$s2, 
          df1 = canton.slope.Loss$p - canton.no.Loss$p, df2 = canton.slope.Loss$n - canton.slope.Loss$p )
  1 - pf( ( ( RSS.no.seplam.Loss - RSS.slope.seplam.Loss )/(canton.slope.seplam.Loss$p - canton.no.seplam.Loss$p) )/
            canton.slope.seplam.Loss$s2, df1 = canton.slope.seplam.Loss$p - canton.no.seplam.Loss$p, 
          df2 = canton.slope.seplam.Loss$n - canton.slope.seplam.Loss$p )
  1 - pchisq( -2*( canton.no.other.Loss$ell - canton.slope.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
  # not signif. (neither for seplam nor TBS)



#(B) Fit new PTBS/TBS model:
  # (B.1) PTBS model
    # Intercept and slope by canton
      mle.canton.both.Loss <- PTBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss, intercept = 1:5, dummy.zeros = TRUE, 
                                        sep.lam = FALSE, interval = c(-3,1) )
      names( mle.canton.both.Loss$beta.hat ) <- colnames(X.canton.aov.Loss)
    # Only intercept by canton
      mle.canton.intcpt.Loss <- PTBS.mle( y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), intercept = 1:5, 
                                          sep.lam = FALSE, interval = c(-3,1) )
      names( mle.canton.intcpt.Loss$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.Loss <- PTBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss[,6:10], intercept = TRUE,
                                         dummy.zeros = TRUE, sep.lam = FALSE, interval = c(-3,1) ) 
      # ! dummy.zeros only works if all(x.aov.Loss != 0) !
      names( mle.canton.slope.Loss$beta.hat ) <- c( 'one', colnames(X.canton.aov.Loss[,6:10]) )
    # Nothing by canton
      mle.canton.no.Loss <- PTBS.mle( y = y.aov.Loss, x = x.aov.Loss, intercept = TRUE, sep.lam = FALSE, interval = c(-3,1) )
      names(mle.canton.no.Loss$beta.hat) <- c('one','x')

    # Covariance matrix of coefficients:
      covML.canton.both.Loss <- solve( -optimHess( unlist(mle.canton.both.Loss[c(2,3,1)]), PTBS.llkhd, y = y.aov.Loss, 
                                                   x = X.canton.aov.Loss, intercept = 1:5, dummy.zeros = TRUE, 
                                                   sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.intcpt.Loss <- solve( -optimHess( unlist(mle.canton.intcpt.Loss[c(2,3,1)]), PTBS.llkhd, y = y.aov.Loss,
                                                     x = cbind( Region.mat.aov, x.aov.Loss ), intercept = 1:5, 
                                                     sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.slope.Loss <- solve( -optimHess( unlist(mle.canton.slope.Loss[c(2,3,1)]), PTBS.llkhd, y = y.aov.Loss, 
                                                    x = X.canton.aov.Loss[,6:10], intercept = TRUE, dummy.zeros = TRUE, 
                                                    sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.no.Loss <- solve( -optimHess( unlist(mle.canton.no.Loss[c(2,3,1)]), PTBS.llkhd, y = y.aov.Loss, 
                                                 x = x.aov.Loss, intercept = TRUE, sep.lam = FALSE,
                                                 control = list( maxit = 5000, fnscale = -1 ) ) )
    
    # Standard errors of coefficients:
      se.canton.both.Loss <- sqrt( diag(covML.canton.both.Loss)[1:length(mle.canton.both.Loss$beta.hat)] )
      se.canton.intcpt.Loss <- sqrt( diag(covML.canton.intcpt.Loss)[1:length(mle.canton.intcpt.Loss$beta.hat)] )
      se.canton.slope.Loss <- sqrt( diag(covML.canton.slope.Loss)[1:length(mle.canton.slope.Loss$beta.hat)] )
      se.canton.no.Loss <- sqrt( diag(covML.canton.no.Loss)[1:p] )


  #(B.2) PTBS model with seplam
    # Intercept and slope by canton
      mle.canton.both.seplam.Loss <- PTBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss, intercept = 1:5, 
                                               dummy.zeros = TRUE, sep.lam = TRUE )
      names( mle.canton.both.seplam.Loss$beta.hat ) <- colnames(X.canton.aov.Loss)
    # Only intercept by canton
      mle.canton.intcpt.seplam.Loss <- PTBS.mle( y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), 
                                                 intercept = 1:5, sep.lam = TRUE )
      names( mle.canton.intcpt.seplam.Loss$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.seplam.Loss <- PTBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss[,6:10], intercept = TRUE, 
                                                dummy.zeros = TRUE, sep.lam = TRUE ) # ! dummy.zeros only works if all(x.aov.Loss != 0) !
      names( mle.canton.slope.seplam.Loss$beta.hat ) <- c( 'one', colnames(X.canton.aov.Loss[,6:10]) )
    # Nothing by canton
      mle.canton.no.seplam.Loss <- PTBS.mle( y = y.aov.Loss, x = x.aov.Loss, intercept = TRUE, sep.lam = TRUE )
      names(mle.canton.no.seplam.Loss$beta.hat) <- c('one','x')

    # Covariance matrix of coefficients:
      covML.canton.both.seplam.Loss <- solve( -optimHess( unlist(mle.canton.both.seplam.Loss[c(2,3,1)]), PTBS.llkhd,
                                                          y = y.aov.Loss, x = X.canton.aov.Loss, intercept = 1:5, 
                                                          dummy.zeros = TRUE, sep.lam = TRUE, 
                                                          control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.intcpt.seplam.Loss <- solve( -optimHess( unlist(mle.canton.intcpt.seplam.Loss[c(2,3,1)]), PTBS.llkhd, 
                                                            y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), 
                                                            intercept = 1:5, sep.lam = TRUE, 
                                                            control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.slope.seplam.Loss <- solve( -optimHess( unlist(mle.canton.slope.seplam.Loss[c(2,3,1)]), PTBS.llkhd, 
                                                           y = y.aov.Loss, x = X.canton.aov.Loss[,6:10],
                                                           intercept = TRUE, dummy.zeros = TRUE, sep.lam = TRUE, 
                                                           control = list( maxit = 5000, fnscale = -1 ) ) )
      covML.canton.no.seplam.Loss <- solve( -optimHess( unlist(mle.canton.no.seplam.Loss[c(2,3,1)]), PTBS.llkhd, 
                                                        y = y.aov.Loss, x = x.aov.Loss, intercept = TRUE, sep.lam = TRUE, 
                                                        control = list( maxit = 5000, fnscale = -1 ) ) )
    
    # Standard errors of coefficients:
      se.canton.both.seplam.Loss <- sqrt( diag(covML.canton.both.seplam.Loss)[1:length(mle.canton.both.seplam.Loss$beta.hat)] )
      se.canton.intcpt.seplam.Loss <- sqrt( diag(covML.canton.intcpt.seplam.Loss)[1:length(mle.canton.intcpt.seplam.Loss$beta.hat)] )
      se.canton.slope.seplam.Loss <- sqrt( diag(covML.canton.slope.seplam.Loss)[1:length(mle.canton.slope.seplam.Loss$beta.hat)] )
      se.canton.no.seplam.Loss <- sqrt( diag(covML.canton.no.seplam.Loss)[1:p] )


  # (B.3) TBS model
    # Nothing by canton
      mle.canton.no.other.Loss <- TBS.mle( y = y.aov.Loss, x = x.aov.Loss, beta.init = c(1000,0.5), intercept = TRUE )
      names(mle.canton.no.other.Loss$beta.hat) <- c('one','x')
    # Only intercept by canton
      mle.canton.intcpt.other.Loss <- TBS.mle( y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), 
                                               beta.init = rep( mle.canton.no.other.Loss$beta.hat, c(5,1) ), 
                                               intercept = FALSE )
      names( mle.canton.intcpt.other.Loss$beta.hat ) <- c( colnames(Region.mat.aov), 'x' )
    # Only slope by canton
      mle.canton.slope.other.Loss <- TBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss[,6:10], 
                                              beta.init = rep( mle.canton.no.other.Loss$beta.hat, c(1,5) ), 
                                              intercept = TRUE ) 
      names( mle.canton.slope.other.Loss$beta.hat ) <- c( 'one', colnames(X.canton.aov.Loss[,6:10]) )
    # Both by canton
      mle.canton.both.other.Loss <- TBS.mle( y = y.aov.Loss, x = X.canton.aov.Loss, 
                                             beta.init = rep( mle.canton.no.other.Loss$beta.hat, each = 5 ), 
                                             intercept = FALSE )
      names( mle.canton.both.other.Loss$beta.hat ) <- colnames(X.canton.aov.Loss)
    
    # Covariance matrix of coefficients:
      covML.canton.both.other.Loss <- solve( -optimHess( unlist(mle.canton.both.other.Loss[c(2,3,1)]), loglik.other, 
                                                         y = y.aov.Loss, x = X.canton.aov.Loss, intercept = FALSE, 
                                                         control = list( fnscale = -1 ) ) )
      covML.canton.intcpt.other.Loss <- solve( -optimHess( unlist(mle.canton.intcpt.other.Loss[c(2,3,1)]), loglik.other, 
                                                           y = y.aov.Loss, x = cbind( Region.mat.aov, x.aov.Loss ), 
                                                           intercept = FALSE, control = list( fnscale = -1 ) ) )
      covML.canton.slope.other.Loss <- solve( -optimHess( unlist(mle.canton.slope.other.Loss[c(2,3,1)]), loglik.other,
                                                          y = y.aov.Loss, x = X.canton.aov.Loss[,6:10], intercept = TRUE, 
                                                          control = list( fnscale = -1 ) ) )
      covML.canton.no.other.Loss <- solve( -optimHess( unlist(mle.canton.no.other.Loss[c(2,3,1)]), loglik.other, 
                                                       y = y.aov.Loss, x = x.aov.Loss, intercept = TRUE, 
                                                       control = list( fnscale = -1 ) ) )
    
    # Standard errors of coefficients:
      se.canton.both.other.Loss <- sqrt( diag(covML.canton.both.other.Loss)[1:length(mle.canton.both.other.Loss$beta.hat)] )
      se.canton.intcpt.other.Loss <- sqrt( diag(covML.canton.intcpt.other.Loss)[1:length(mle.canton.intcpt.other.Loss$beta.hat)] )
      se.canton.slope.other.Loss <- sqrt( diag(covML.canton.slope.other.Loss)[1:length(mle.canton.slope.other.Loss$beta.hat)] )
      se.canton.no.other.Loss <- sqrt( diag(covML.canton.no.other.Loss)[1:p] )


#=== Comparisons (newly fitted models)
  # (a) both vs. nothing
    1 - pchisq( -2*( mle.canton.no.Loss$ell - mle.canton.both.Loss$ell ), df = ncol(X.canton.aov.Loss) - 2 )
    1 - pchisq( -2*( mle.canton.no.seplam.Loss$ell - mle.canton.both.seplam.Loss$ell ), df = ncol(X.canton.aov.Loss) - 2 )
    1 - pchisq( -2*( mle.canton.no.other.Loss$ell - mle.canton.both.other.Loss$ell ), df = ncol(X.canton.aov.Loss) - 2 )
    # not signif. (neither for seplam nor TBS)
  # (b) both vs. intercept
    1 - pchisq( -2*( mle.canton.intcpt.Loss$ell - mle.canton.both.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.intcpt.seplam.Loss$ell - mle.canton.both.seplam.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.intcpt.other.Loss$ell - mle.canton.both.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    # not signif. (neither for seplam nor TBS)
  # (c) both vs. slope
    1 - pchisq( -2*( mle.canton.slope.Loss$ell - mle.canton.both.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.slope.seplam.Loss$ell - mle.canton.both.seplam.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.slope.other.Loss$ell - mle.canton.both.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    # not signif. (neither for seplam nor TBS)
  # (d) intercept vs. nothing 
    1 - pchisq( -2*( mle.canton.no.Loss$ell - mle.canton.intcpt.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.no.seplam.Loss$ell - mle.canton.intcpt.seplam.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.no.other.Loss$ell - mle.canton.intcpt.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    # not signif. (neither for seplam nor TBS)
  # (e) slope vs. nothing 
    1 - pchisq( -2*( mle.canton.no.Loss$ell - mle.canton.slope.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.no.seplam.Loss$ell - mle.canton.slope.seplam.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    1 - pchisq( -2*( mle.canton.no.other.Loss$ell - mle.canton.slope.other.Loss$ell ), df = ncol(Region.mat.aov) - 1 )
    # not signif. (neither for seplam nor TBS)


# Overview tables of the individual model coefficients (only for new estimations B; not adapted for the three models)
  # both
    cbind( estim = mle.canton.both.Loss$beta.hat, se = se.canton.both.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.both.Loss$beta.hat/se.canton.both.Loss) ) ) )
    cbind( estim = mle.canton.both.seplam.Loss$beta.hat, se = se.canton.both.seplam.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.both.seplam.Loss$beta.hat/se.canton.both.seplam.Loss) ) ) )
    cbind( estim = mle.canton.both.other.Loss$beta.hat, se = se.canton.both.other.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.both.other.Loss$beta.hat/se.canton.both.other.Loss) ) ) )
  # intercept
    cbind( estim = mle.canton.intcpt.Loss$beta.hat, se = se.canton.intcpt.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.Loss$beta.hat/se.canton.intcpt.Loss) ) ) )
    cbind( estim = mle.canton.intcpt.seplam.Loss$beta.hat, se = se.canton.intcpt.seplam.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.seplam.Loss$beta.hat/se.canton.intcpt.seplam.Loss) ) ) )
    cbind( estim = mle.canton.intcpt.other.Loss$beta.hat, se = se.canton.intcpt.other.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.intcpt.other.Loss$beta.hat/se.canton.intcpt.other.Loss) ) ) )
  # slope
    cbind( estim = mle.canton.slope.Loss$beta.hat, se = se.canton.slope.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.slope.Loss$beta.hat/se.canton.slope.Loss) ) ) )
    cbind( estim = mle.canton.slope.seplam.Loss$beta.hat, se = se.canton.slope.seplam.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.slope.seplam.Loss$beta.hat/se.canton.slope.seplam.Loss) ) ) )
    cbind( estim = mle.canton.slope.other.Loss$beta.hat, se = se.canton.slope.other.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.slope.other.Loss$beta.hat/se.canton.slope.other.Loss) ) ) )
  # nothing
    cbind( estim = mle.canton.no.Loss$beta.hat, se = se.canton.no.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.no.Loss$beta.hat/se.canton.no.Loss) ) ) )
    cbind( estim = mle.canton.no.seplam.Loss$beta.hat, se = se.canton.no.seplam.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.no.seplam.Loss$beta.hat/se.canton.no.seplam.Loss) ) ) )
    cbind( estim = mle.canton.no.other.Loss$beta.hat, se = se.canton.no.other.Loss, 
           pval = 2*( 1 - pnorm( abs(mle.canton.no.other.Loss$beta.hat/se.canton.no.other.Loss) ) ) )



