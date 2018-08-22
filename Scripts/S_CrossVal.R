############################################## Cross-Validation ########################################################
# =====================================================================================================================#
# =========== IF S_RelLossAnalysis.R & S_AbsLossAnalysis.R ar skipped, activate the the content of this box ===========#

#   idx.outl <- 145
#   
#   y.Loss <- Contents$Loss[-idx.outl]
#   x.Loss <- Structure$Loss[-idx.outl]
# 
#   y.DoL <- Contents$DoL[-idx.outl]
#   x.DoL <- Structure$DoL[-idx.outl]
# 
#   n.obs <- length(y.Loss)
#   p <- 2
#   
#   # Parameter estimation, PTBS same lambda:
#     mle.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = FALSE )
#     mle.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = FALSE )
# =====================================================================================================================#

# PTBS sep.lam:
  mle.seplam.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = TRUE )
  s2.hat.seplam.Loss <- mle.seplam.Loss$sigma2.hat * n.obs/(n.obs-p)

# TBS:
  mle.other.Loss <- TBS.mle( y = y.Loss, x = x.Loss, beta.init = c(1000,0.5), interval = c(-3,1) )
  s2.hat.other.Loss <- mle.other.Loss$sigma2.hat * n.obs/(n.obs-p)

# Insurance sum of contents:
  y.InSum <- Contents$InSum[-idx.outl]

# Parameter estimation:
  # PTBS sep.lam:
  mle.seplam.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = TRUE )
  s2.hat.seplam.DoL <- mle.seplam.DoL$sigma2.hat * n.obs/(n.obs-p)
  
  # TBS:
    mle.other.DoL <- TBS.mle( y = y.DoL, x = x.DoL, beta.init = c(0.1,1) )
    s2.hat.other.DoL <- mle.other.DoL$sigma2.hat * n.obs/(n.obs-p)

# === Functions for CV =================================================================================================
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
err.measures <- function( pred, target ) { # all relative quantities are in percent
	diffs <- pred - target
	c( rmspe = sqrt( mean( diffs^2 ) ), bias = mean( diffs ), rel.tot.bias = 100 * mean(diffs)/mean(target), 
	   rel.bias = 100 * ( mean( pred/target ) - 1 ), mae = mean( abs(diffs) ), rel.tot.ae = 100 * 
	     mean( abs(diffs) )/mean(target), rel.ae = 100 * mean( abs( 1 - pred/target ) ) )  
}

err.measures.med <- function( pred, target ) { # all relative quantities are in percent
	diffs <- pred - target
	c( med.bias = median( diffs ), med.rel.bias = 100 * median( diffs/target ), med.ae = median( abs(diffs) ), 
	   med.rel.ae = 100 * median( abs( diffs/target ) ) ) 
}

pred.models.comb <- function( x, mle, mle.seplam, mle.other, type ) {
	# TBS model:
  	fitted.TBS <- c( cbind( 1, x ) %*% mle.other$beta.hat )
  	# Median prediction
  	  TBS.med <- fitted.TBS
  	# Mean prediction
  	  TBS.mean <- TBS.condmean.orig( theta = unlist(mle.other[c(2,3,1)]), x.beta = fitted.TBS )
  	  names(TBS.mean) <- NULL
	# PTBS model
  	# Same lambda:
    	# Median prediction:
      	PTBS.med <- BC.backtransform( mle$lambda.hat[1], c( cbind( 1, BC.transform( mle$lambda.hat[1+0], x ) ) %*% 
      	                                                      mle$beta.hat ) )
    	# Mean prediction:
      	PTBS.mean <- condmean.orig( theta = unlist(mle[c(2,3,1)]), x = cbind( 1, BC.transform( mle$lambda.hat[1+0], x)),
      	                            sep.lam = FALSE )
    	  names(PTBS.mean) <- NULL
	# Separate lambda
  	# Median prediction:
  	  PTBS.seplam.med <- BC.backtransform( mle.seplam$lambda.hat[1], c( 
  	    cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ) %*% mle.seplam$beta.hat ) )
  	# Mean prediction:
    	PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.seplam[c(2,3,1)]), 
    	                                   x = cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ), sep.lam = TRUE )
    	names(PTBS.seplam.mean) <- NULL	
	  # Predictions (matrix n.obs by |{error measures}| )
  	  ret <- cbind( TBS.med, PTBS.med, PTBS.seplam.med, TBS.mean, PTBS.mean, PTBS.seplam.mean )
  	  colnames(ret) <- paste( type, colnames(ret), sep = '.' )
  	  ret
}

pbias.models.comb.transf <- function( x, mle, mle.seplam, mle.other, y, type ) {
	# TBS model:
  	pbias.TBS.tr <- BC.transform( mle.other$lambda.hat, c( cbind( 1, x ) %*% mle.other$beta.hat ) ) - 
  	  BC.transform( mle.other$lambda.hat, y )
	# PTBS model:
  	# Same Lambda
	  pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle$lambda.hat[1+0], x ) ) %*% mle$beta.hat ) - 
	    BC.transform( mle$lambda.hat[1], y )
  	# Separate lambda
	    pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ) %*% mle.seplam$beta.hat ) - 
	      BC.transform( mle.seplam$lambda.hat[1], y )
  	# Predictions (matrix n.obs by |{error measures}| )
    	ret <- cbind( pbias.TBS.tr, pbias.PTBS.tr, pbias.PTBS.seplam.tr )
    	colnames(ret) <- paste( type, c( "TBS.tr", "PTBS.tr", "PTBS.seplam.tr" ), sep = '.' )
    	ret
}
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ==================================================================================================================== #

# Mean approximation by Taylor series (ok, i.e. these quantities are rather small):
  summary( mle.Loss$lambda.hat * sqrt(s2.hat.Loss)/( 1 + mle.Loss$lambda.hat * 
                                                       c( cbind( 1, BC.transform( mle.Loss$lambda.hat, x.Loss ) ) %*% 
                                                            mle.Loss$beta.hat ) ) )
  
  summary( mle.seplam.Loss$lambda.hat[1] * sqrt(s2.hat.seplam.Loss)/
             ( 1 + mle.seplam.Loss$lambda.hat[1] * c( cbind( 1, BC.transform(mle.seplam.Loss$lambda.hat[2], x.Loss)) %*%
                                                        mle.seplam.Loss$beta.hat ) ) )
  
  summary( mle.other.Loss$lambda.hat * sqrt(s2.hat.other.Loss)/c( cbind( 1, x.Loss ) %*% 
                                                                    mle.other.Loss$beta.hat )^mle.other.Loss$lambda.hat)
  
  summary( mle.DoL$lambda.hat * sqrt(s2.hat.DoL)/
             ( 1 + mle.DoL$lambda.hat * c( cbind( 1, BC.transform( mle.DoL$lambda.hat, x.DoL ) ) %*% mle.DoL$beta.hat)))
  
  summary( mle.seplam.DoL$lambda.hat[1] * sqrt(s2.hat.seplam.DoL)/
             ( 1 + mle.seplam.DoL$lambda.hat[1] *c( cbind( 1, BC.transform( mle.seplam.DoL$lambda.hat[2], x.DoL ) ) %*% 
                                                      mle.seplam.DoL$beta.hat ) ) )
  
  summary( mle.other.DoL$lambda.hat * sqrt(s2.hat.other.DoL)/c( cbind( 1, x.DoL ) %*% 
                                                                  mle.other.DoL$beta.hat )^mle.other.DoL$lambda.hat )

  
##=======================================##
##==== Absolute loss ====================##
##=======================================##

# Leave-one-out cross validation
  pred.fullCV.Loss <- t( sapply( 1:n.obs, function(k){ # Predictions (matrix n.obs by |{model:pred}|)
  	# Fit model to y[-k], giving estimated parameters theta.fit.
  	# TBS model:
    	mle.k <- TBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
    	# Predict y[k] by yhat( theta.fit, x[k] ).
    	  fitted.k <- c( cbind( 1, x.Loss[k] ) %*% mle.k$beta.hat )
    	# Median prediction
    	  TBS.med <- fitted.k
    	# Mean prediction
    	  TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k )
  	    names(TBS.mean) <- NULL
  	# PTBS model:
      # Same lambda:
  	    mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
      	# Median prediction:
      	  PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ) 
      	                                                        %*% mle.k$beta.hat ) )
      	# Mean prediction:
        	PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), 
        	                            x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ), sep.lam = FALSE )
        	names(PTBS.mean) <- NULL
  	  # Separate lambda
  	    mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
      	# Median prediction:
      	  PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( 
      	    cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[k] ) ) %*% mle.k$beta.hat ) )
      	# Mean prediction:
      	  PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( 
      	    mle.k$lambda.hat[1+1], x.Loss[k] ) ), sep.lam = TRUE )
      	  names(PTBS.seplam.mean) <- NULL
      	  c( Loss.TBS.med = TBS.med, Loss.PTBS.med = PTBS.med, Loss.PTBS.seplam.med = PTBS.seplam.med, 
      	     Loss.TBS.mean = TBS.mean, Loss.PTBS.mean = PTBS.mean, Loss.PTBS.seplam.mean = PTBS.seplam.mean )
  } ) )
  
# Apparent error (prediction of used data under model):
  pred.mod.Loss <- pred.models.comb( x = x.Loss, mle = mle.Loss, mle.seplam = mle.seplam.Loss, 
                                   mle.other = mle.other.Loss, type = "Loss" ) # n.obs by |{model:pred}|
  
# Aggregate prediction errors (|{model:pred}| by |{err.measures}|)
  Delta.CV.Loss <- t( apply( pred.fullCV.Loss, 2, err.measures, target = y.Loss ) )
  Delta.CV.med.Loss <- t( apply( pred.fullCV.Loss, 2, err.measures.med, target = y.Loss ) )

  
  
  plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Contents", main = "Absolute loss, original scale" )
  points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4 ) # ev. pch = 8
  # PTBS same lambda:
    # median:
      lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ), col = 'red' ) # fitted median
      points( x.Loss, pred.fullCV.Loss[,2], col = 'red', pch = 3)
    # conditional mean:
      lines( x.seq.Loss, condmean.seq.Loss, col = 'green3', lwd = 1)
      points( x.Loss, pred.fullCV.Loss[,5], col = 'green3', pch = 3 )
  # TBS:
    # median:
      lines( x.seq.Loss, c( cbind(1,x.seq.Loss) %*% mle.other.Loss$beta.hat ), col = 'blue' )
      points( x.Loss, pred.fullCV.Loss[,1], col = 'blue', pch = 3 )
    # mean:
      points( x.Loss, pred.fullCV.Loss[,4], col = 'purple', pch = 3 )
  # PTBS sep.lam:
      # median:
        lines( x.Loss, pred.fullCV.Loss[,3], col = 'salmon', pch = 4 )
      # mean:
        points( x.Loss, pred.fullCV.Loss[,6], col = 'skyblue', pch = 4 )
  
## =========================================================##
## ==== Prediction on transformed scale ====================##
## =========================================================##
    
    # Fit model to y[-k], giving estimated parameters theta.fit.
    pbias.transf.fullCV.Loss <- t( sapply( 1:n.obs, function(k){
    	# TBS model
      	mle.k <- TBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
      	pbias.TBS.tr <- BC.transform( mle.k$lambda.hat, c( cbind( 1, x.Loss[k] ) %*% mle.k$beta.hat ) ) -
      	  BC.transform( mle.k$lambda.hat, y.Loss[k] )
    	# PTBS model
    	  mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
    	  pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ) %*% mle.k$beta.hat ) -
    	    BC.transform( mle.k$lambda.hat[1], y.Loss[k] )
    	# Separate lambda
      	mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
      	pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[k] ) ) %*% mle.k$beta.hat ) -
      	  BC.transform( mle.k$lambda.hat[1], y.Loss[k] )
      	cbind( CV.TBS.tr.Loss = pbias.TBS.tr, CV.PTBS.tr.Loss = pbias.PTBS.tr, CV.PTBS.seplam.Loss = pbias.PTBS.seplam.tr )
    } ) )

    s2.hat.Loss <- mle.Loss$sigma2.hat * n.obs/(n.obs-p-1)
    covML.Loss <- solve( -optimHess( unlist(mle.Loss[c(2,3,1)]), PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = FALSE,
                                     control = list( maxit = 5000, fnscale = -1 ) ) )
    var.pred.PTBS.transf.Loss <- sapply( BC.transform( mle.Loss$lambda.hat[1], x.Loss ), function(x){ c(1,x) %*%
                                          covML.Loss[1:p,1:p] %*% c(1,x) } ) + s2.hat.Loss

    s2.hat.seplam.Loss <- mle.seplam.Loss$sigma2.hat * n.obs/(n.obs-p-2)
    covML.seplam.Loss <- solve( -optimHess( unlist(mle.seplam.Loss[c(2,3,1)]), PTBS.llkhd, y = y.Loss, x = x.Loss, 
                                            sep.lam = TRUE, control = list( maxit = 5000, fnscale = -1 ) ) )
    var.pred.PTBS.seplam.transf.Loss <- sapply( BC.transform( mle.seplam.Loss$lambda.hat[2], x.Loss ), 
                                                function(x){ c(1,x) %*% covML.seplam.Loss[1:p,1:p] %*% c(1,x) } 
                                                ) + s2.hat.seplam.Loss

    s2.hat.other.Loss <- mle.other.Loss$sigma2.hat * n.obs/(n.obs-p-1)
    covML.other.Loss <- solve( -optimHess( unlist( mle.other.Loss[c(2,3,1)] ), loglik.other, y = y.Loss, x = x.Loss, 
                                           control = list( fnscale = -1 ) ) )
    var.pred.TBS.transf.Loss <- sapply( x.Loss, function(x){ c( c( 1, x ) %*% mle.other.Loss$beta.hat )^
        (2*(mle.other.Loss$lambda.hat - 1)) * c( c( 1, x ) %*% covML.other.Loss[1:p,1:p] %*% c(1, x) ) } 
        ) + s2.hat.other.Loss
    # With regression version of var(beta.hat) (gives similar result):
    sapply( x.Loss, function(x){ c( c( 1, x ) %*% mle.other.Loss$beta.hat )^(2*(mle.other.Loss$lambda.hat - 1)) * 
        c( c( 1, x ) %*% ( s2.hat.other.Loss * solve( crossprod( diag( c( cbind( 1, x.Loss ) %*% mle.other.Loss$beta.hat )^(mle.other.Loss$lambda.hat - 1) ) %*% cbind( 1, x.Loss ) ) ) ) %*% c(1, x) ) } ) + s2.hat.other.Loss

    var.pred.transf.Loss <- cbind( TBS = var.pred.TBS.transf.Loss, PTBS = var.pred.PTBS.transf.Loss, PTBS.seplam = var.pred.PTBS.seplam.transf.Loss )

    par(mfrow=c(3,1))
    for(i in 1:3) {
    	boxplot( t( pbias.transf.fullCV.Loss[ranking,i] ), outlty = 1, outpch = NA, xlab = "Rank of target loss", 
    	         ylab = "Prediction bias [transf.CHF]", main = dimnames(pbias.transf.fullCV.Loss)[i] )
    	lines( -qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.Loss[ranking,i]), col = 'blue', lty = 2 )
    	lines( qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.Loss[ranking,i]), col = 'blue', lty = 2 )
    	abline( h = 0, col = 'blue' )
    	rug( which( ( apply( pbias.transf.fullCV.Loss, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * 
    	                sqrt(var.pred.transf.Loss[,i]) | apply( pbias.transf.fullCV.Loss, 1:2, quantile, probs = 0.975 ) 
    	              < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss[,i]) )[ranking,i] ), col = 'red' )
    	box()
    }


##=======================================##
##==== Relative loss ====================##
##=======================================##
    

# Leave-one-out cross validation
  # Predictions (matrix n.obs by |{model:pred}|)
    pred.fullCV.DoL <- t( sapply( 1:n.obs, function(k){
      # Fit model to y[-k], giving estimated parameters theta.fit.
    	# TBS model:
    	    mle.k <- TBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, beta.init = c(0.1,1) ) 
    	    # Predict y[k] by yhat( theta.fit, x[k] ).
    	      fitted.k <- c( cbind( 1, x.DoL[k] ) %*% mle.k$beta.hat )
    	    # Median prediction
    	      TBS.med <- fitted.k * y.InSum[k]
    	    # Mean prediction
          	TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) * y.InSum[k]
          	names(TBS.mean) <- NULL
    	# PTBS model:
          # Same lambda:
    	      mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = FALSE )
      	    # Median prediction:
      	      PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], 
      	                                    c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ) %*% 
      	                                         mle.k$beta.hat ) ) * y.InSum[k]
      	    # Mean prediction:
            	PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), 
            	                            x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ), 
            	                            sep.lam = FALSE ) * y.InSum[k]
            	names(PTBS.mean) <- NULL
    	    # Separate lambda:
          	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = TRUE )
          	# Median prediction:
          	  PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( 
          	                      cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ) %*% 
          	                        mle.k$beta.hat ) ) * y.InSum[k]
          	# Mean prediction:
          	  PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), 
          	                                     x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ), 
          	                                     sep.lam = TRUE ) * y.InSum[k]
          	  names(PTBS.seplam.mean) <- NULL
    	c( DoL.TBS.med = TBS.med, DoL.PTBS.med = PTBS.med, DoL.PTBS.seplam.med = PTBS.seplam.med, DoL.TBS.mean = TBS.mean, 
    	   DoL.PTBS.mean = PTBS.mean, DoL.PTBS.seplam.mean = PTBS.seplam.mean )
  } ) ) # n.obs by 6
  
  # Aggregate prediction errors (|{model:pred}| by |{err.measures}|)
    Delta.CV.DoL <- t( apply( pred.fullCV.DoL, 2, err.measures, target = y.DoL * y.InSum ) )
    Delta.CV.med.DoL <- t( apply( pred.fullCV.DoL, 2, err.measures.med, target = y.DoL * y.InSum ) )

  # Apparent error (prediction of used data under model):
    pred.mod.DoL <- pred.models.comb( x = x.DoL, mle = mle.DoL, mle.seplam = mle.seplam.DoL, mle.other = mle.other.DoL, 
                                      type = "DoL" ) * y.InSum  # n.obs by |{model:pred}|
    

  # Individual biases for mean prediction (resamples do or don't contain obs j)
    # x11()
    ranking <- order(y.Loss)
    par( mfcol = c(3, 2))
    for(i in 1:6) {
    	boxplot( 1:n.obs, (pred.fullCV.DoL[,i] - y.DoL * y.InSum)[ranking], outlty = 1, outpch = NA, xlab = "Rank of target value", 
    	         ylab = "Prediction bias [CHF]", main = dimnames(pred.fullCV.DoL)[[2]][i] )
    	points( 1:n.obs, (pred.fullCV.DoL[,i] - y.DoL * y.InSum)[ranking], col = 'red', pch = 3, cex = 0.8 )
    	abline( h = seq(-140000,320000,10000), col = 'grey', lty = 3 )
    	abline( h = 0, col = 'blue' )
    	rug( which( apply( pred.fullCV.DoL - y.DoL * y.InSum, 1:2, function(e){ min(e) <= 0 & max(e) >= 0 } )[ranking,i] ), col = 'green3' )
    	box()
    }
    # savePlot( paste( "Figures/CV_", modelname, "_individualBias", sep = '' ), type = "pdf" )

##=======================================##
##==== Prediction on transformed scale===##
##=======================================##
    
# Prediction bias on transformed scale (same for mean and median), n.obs by |models| by R
# Fit model to y[-k], giving estimated parameters theta.fit.
  pbias.transf.fullCV.DoL <- t( sapply( 1:n.obs, function(k){
  	# TBS model
  	  mle.k <- TBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, beta.init = c(0.1,1) )
  	  pbias.TBS.tr <- BC.transform( mle.k$lambda.hat, c( cbind( 1, x.DoL[k] ) %*% mle.k$beta.hat ) ) - 
  	    BC.transform( mle.k$lambda.hat, y.DoL[k] )
  	# PTBS model
  	  # Same Lambda
  	    mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = FALSE )
  	    pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ) %*% mle.k$beta.hat ) -
  	      BC.transform( mle.k$lambda.hat[1], y.DoL[k] )
  	  # Separate lambda
      	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = TRUE )
      	pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ) %*% mle.k$beta.hat ) - 
      	  BC.transform( mle.k$lambda.hat[1], y.DoL[k] )
  	cbind( CV.TBS.tr.DoL = pbias.TBS.tr, CV.PTBS.tr.DoL = pbias.PTBS.tr, CV.PTBS.seplam.DoL = pbias.PTBS.seplam.tr )
  } ) )
  
  
  s2.hat.DoL <- mle.DoL$sigma2.hat * n.obs/(n.obs-p-1)
  covML.DoL <- solve( -optimHess( unlist(mle.DoL[c(2,3,1)]), PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = FALSE, 
                                  control = list( maxit = 5000, fnscale = -1 ) ) )
  var.pred.PTBS.transf.DoL <- sapply( BC.transform( mle.DoL$lambda.hat[1], x.DoL ), function(x){ c(1,x) %*% 
      covML.DoL[1:p,1:p] %*% c(1,x) } ) + s2.hat.DoL
  
  s2.hat.seplam.DoL <- mle.seplam.DoL$sigma2.hat * n.obs/(n.obs-p-2)
  covML.seplam.DoL <- solve( -optimHess( unlist(mle.seplam.DoL[c(2,3,1)]), PTBS.llkhd, y = y.DoL, x = x.DoL, 
                                         sep.lam = TRUE, control = list( maxit = 5000, fnscale = -1 ) ) )
  var.pred.PTBS.seplam.transf.DoL <- sapply( BC.transform( mle.seplam.DoL$lambda.hat[2], x.DoL ), 
                                             function(x){ c(1,x) %*% covML.seplam.DoL[1:p,1:p] %*% c(1,x) } 
                                             ) + s2.hat.seplam.DoL
  s2.hat.other.DoL <- mle.other.DoL$sigma2.hat * n.obs/(n.obs-p-1)
  covML.other.DoL <- solve( -optimHess( unlist( mle.other.DoL[c(2,3,1)] ), loglik.other, y = y.DoL, x = x.DoL, 
                                        control = list( fnscale = -1 ) ) )
  var.pred.TBS.transf.DoL <- sapply( x.DoL, function(x){ c( c( 1, x ) %*% 
                                                              mle.other.DoL$beta.hat )^(2*(mle.other.DoL$lambda.hat - 1)) * 
      c( c( 1, x ) %*% covML.other.DoL[1:p,1:p] %*% c(1, x) ) } ) + s2.hat.other.DoL

# With regression version of var(beta.hat) (gives similar result):
#   sapply( x.DoL, function(x){ c( c( 1, x ) %*% mle.other.DoL$beta.hat )^(2*(mle.other.DoL$lambda.hat - 1)) * 
#       c( c( 1, x ) %*% ( s2.hat.other.DoL * solve( crossprod( diag( c( cbind( 1, x.DoL ) %*% mle.other.DoL$beta.hat )^
#                                                                       (mle.other.DoL$lambda.hat - 1) ) %*% 
#                                                                 cbind( 1, x.DoL ) ) ) ) %*% c(1, x) ) } 
#         ) + s2.hat.other.DoL
  
  var.pred.transf.DoL <- cbind( TBS = var.pred.TBS.transf.DoL, PTBS = var.pred.PTBS.transf.DoL, PTBS.seplam = var.pred.PTBS.seplam.transf.DoL )
  
  par(mfrow=c(3,1))
  for(i in 1:3) {
    boxplot( t( pbias.transf.fullCV.DoL[ranking,i] ), outlty = 1, outpch = NA, xlab = "Rank of target loss", 
             ylab = "Prediction bias [transf.CHF]", main = dimnames(pbias.transf.fullCV.DoL)[i] )
    lines( -qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.DoL[ranking,i]), col = 'blue', lty = 2 )
    lines( qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.DoL[ranking,i]), col = 'blue', lty = 2 )
    abline( h = 0, col = 'blue' )
    rug( which( ( apply( pbias.transf.fullCV.DoL, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * 
                    sqrt(var.pred.transf.DoL[,i]) | apply( pbias.transf.fullCV.DoL, 1:2, quantile, probs = 0.975 ) 
                  < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL[,i]) )[ranking,i] ), col = 'red' )
    box()
  }
  

##=======================================##
##============ Summary plots ============##
##=======================================##

library(abind)
# BIAS:
  se.bias <- cbind( mod = apply( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) 
                                 - y.Loss, 2, sd ), 
                    CV = apply( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], 
                                       pred.fullCV.DoL[,4:6] ) - y.Loss, 2, sd ) )
  apply( se.bias, 1, diff )

# Relative BIAS:
  se.rel.bias <- cbind( mod = apply( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], 
                                                  pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss, 2, sd ), 
                        CV = apply( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], 
                                                 pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss, 2, sd ) )
  apply( se.rel.bias, 1, diff )

# Mean absolute error:
  se.mae <- cbind( mod = apply( abs( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], 
                                            pred.mod.DoL[,4:6] ) - y.Loss ), 2, sd ), 
                   CV = apply( abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], 
                                           pred.fullCV.DoL[,4:6] ) - y.Loss ), 2, sd ) )
  apply( se.mae, 1, diff )

  
# Relative absolute error:
  se.rel.ae <- cbind( mod = apply( abs( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], 
                                                     pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss ), 2, sd ), 
                      CV = apply( abs( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], 
                                                    pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss ),
                                  2, sd ) )
  apply( se.rel.ae, 1, diff ) # all positive up to one case
  
# high variability across the sample

# bring together:
  se.preds <- abind( se.rmspe = matrix( NA, 12, 2 ), se.bias = se.bias, se.rel.tot.bias = matrix( NA, 12, 2 ), 
                     se.rel.bias = se.rel.bias, se.mae = se.mae, se.rel.tot.ae = matrix( NA, 12, 2 ), 
                     se.rel.ae = se.rel.ae, along = 3 )
  
  err.name <- c( "Root mean square prediction error [CHF]", "Bias [CHF]", "Relative total bias [%]", "Relative bias [%]", 
                 "Mean absolute error [CHF]", "Relative total absolute error [%]", "Relative absolute error [%]" )


# Summary plots for prediction on original scale:
  x11(width = 12, height = 9)
  par( mfrow = c(2,2), mar = c( 9.5, 4, 3.5, 4 ), cex.axis=0.9)
  for( i in c(2,4,5,7) ) {
  	plot( 1:12, c( Delta.CV.Loss[1:3,i], Delta.CV.DoL[1:3,i], 
  	                                                                  Delta.CV.Loss[4:6,i], Delta.CV.DoL[4:6,i]
  	                                                                  #, 
  	                                                                  # Delta.B.Loss[1:3,i], Delta.B.DoL[1:3,i], 
  	                                                                  # Delta.B.Loss[4:6,i], Delta.B.DoL[4:6,i], 
  	                                                                  # Delta.632.Loss[1:3,i], Delta.632.DoL[1:3,i], 
  	                                                                  # Delta.632.Loss[4:6,i], Delta.632.DoL[4:6,i], 
  	                                                                  # Delta.app.Loss[1:3,i], Delta.app.DoL[1:3,i], 
  	                                                                  # Delta.app.Loss[4:6,i], Delta.app.DoL[4:6,i] 
  	                                                                  ), 
  	      xlab = "", ylab = err.name[i], main = "Aggregate prediction error and variability", 
  	      pch = rep( c(4,1,0,6), each = 12 ), xaxt = 'n', mgp = c(2.5,1,0) )
    axis( 1, at = 1:12, labels = c( rownames(Delta.CV.Loss)[1:3], rownames(Delta.CV.DoL)[1:3], rownames(Delta.CV.Loss)[4:6],
                                    rownames(Delta.CV.DoL)[4:6] ), las = 2 )
  	abline( h = 0, col = 'red' )
  	abline( v = seq(0.5,12.5,1), lty = 3 )
  	if(i!=1) {
  		par(new = TRUE,  cex.axis=0.9, mar = c( 9.5, 4, 3.5, 4 ))
  		plot( 1:12, se.preds[, 'CV' , i], axes = FALSE, xlab = "", pch = 18,
  		      ylab = "", type = 'b', lty = 6, col = 'blue' )
  		axis( 4, at = pretty(range(se.preds[, , i])), col = 'blue', col.axis = 'blue' )
  		mtext( paste( "Standard deviation in the sample [", ifelse( i==2 | i==5, "CHF", "%"), "]", sep = '' ), side = 4, 
  		       line = 2.5, col= 'blue', cex = par("cex") )
  	}
  }
  savePlot( paste( "Figures/CV_", modelname, "_ModMetrics", sep = '' ), type = "pdf" )
  

# General patterns in aggregate errors:
  # Bias smaller for mean than for median; se smaller for DoL than for Loss
  # Absolute error smaller for DoL than for Loss and for median than for mean; se smaller for DoL than for Loss
  # Relative bias and relative absolute error smaller for median than for mean and smaller for Loss than for DoL; se same pattern
  # RMSPE smaller for DoL than for Loss; no se available

# Ranks of the models in terms of different prediction errors
  # Bias: 11, 12, 10, 7, 9, 8, 1, 5, 3, 6, 2, 4
  # Relative bias: 1, 3, 1, 4, 6, 5, 7, 9, 7, 10, 12, 11
  # MAE: 7, 9, 9 , 2, 1, 3, 8, 11, 12, 5, 3, 6
  # Relative AE: 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11
  # RMSE: 8, 11, 10, 3, 2, 3, 7, 9, 12, 6, 1, 5
  err.ranks <- rbind( bias = c( 11, 12, 10, 7, 9, 8, 1, 5, 3, 6, 2, 4 ), 
                      rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), 
                      mae = c( 7, 9.5, 9.5 , 2, 1, 3.5, 8, 11, 12, 5, 3.5, 6 ), 
                      rel.ae = c( 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11 ), 
                      rmse = c( 8, 11, 10, 3.5, 2, 3.5, 7, 9, 12, 6, 1, 5 ) )
  
  apply( err.ranks, 2, sum ) # 4
  apply( err.ranks == apply( err.ranks, 1, min ), 2, sum ) # 1
  apply( err.ranks == apply( err.ranks, 1, max ), 2, sum ) # not 9, 11, 2
  apply( err.ranks, 2, max ) # 4
  apply( err.ranks, 2, min ) # 1, 3, 5, 7, 11

# Ranks of the models in terms of the std deviation for different prediction errors
  # Bias: 7.5, 11, 10, 4, 2, 3, 7.5, 9, 12, 6, 1, 5
  # Relative bias: 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11
  # MAE: 9, 11.5, 10, 3.5, 2, 3.5, 7, 8, 11.5, 6, 1, 5
  # Relative AE: 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11
  se.ranks <- rbind( se.bias = c( 7.5, 11, 10, 4, 2, 3, 7.5, 9, 12, 6, 1, 5 ), 
                     se.rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), 
                     se.mae = c( 9, 11.5, 10, 3.5, 2, 3.5, 7, 8, 11.5, 6, 1, 5 ), 
                     se.rel.ae = c( 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ) )
  apply( se.ranks, 2, sum ) # 4(, 5)
  apply( se.ranks == apply( se.ranks, 1, min ), 2, sum ) # 1, 3, 11
  apply( se.ranks == apply( se.ranks, 1, max ), 2, sum ) # not 9, 11, 2
  apply( se.ranks, 2, max ) # 4
  apply( se.ranks, 2, min ) # 1, 3, 11, not 7:9

# Combined ranks:
  floor( err.ranks[1:4,] ) + floor( se.ranks )
  
  comb.ranks <- rbind( bias = c( 10, 12, 11, 5, 5, 5, 2, 8, 9, 7, 1, 3 ), rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), mae = c( 8, 11, 9.5, 3, 1, 4, 7, 9.5, 12, 5.5, 2, 5.5 ), rel.ae = c( 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11 ) )
  apply( comb.ranks, 2, sum ) # 4 followed by 5, 6, 1
  
  apply( comb.ranks, 2, sum ) # 4(, 5)
  apply( comb.ranks == apply( comb.ranks, 1, min ), 2, sum ) # 1, 3, 5, 11
  apply( comb.ranks == apply( comb.ranks, 1, max ), 2, sum ) # not 11, 9, 2
  apply( comb.ranks, 2, max ) # 4, 6
  apply( comb.ranks, 2, min ) # 1, 3, 5, 11, not 7:8


# TBS Loss mean  (best in terms of bias for Loss, has smallest mean bias but larger se than other models) [#7 overall]
  pdf(file = paste( "Figures/CV_", modelname, "_bestperformers.pdf", sep = '' ))
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.Loss[,4]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.Loss[,4]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.Loss[,4]-y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,4]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "TBS model for absolute loss, mean prediction - best in absolute bias", outer = TRUE, line = 0.8 )

# PTBS DoL mean (overall best in terms of bias when accounting also for se, best in terms of RMSPE, 
#   co-best in terms of absolute error) [#11 overall]
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.DoL[,5]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.DoL[,5]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.DoL[,5]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,5]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "PTBS model for relative loss, mean prediction - best in RMSPE", outer = TRUE, line = 0.8 )

  
  
# PTBS DoL median (overall best in terms of absolute error, practically together with PTBS DoL mean) [#5 overall]
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.DoL[,2]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.DoL[,2]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.DoL[,2]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,2]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "PTBS model for relative loss, median prediction - best in abs. error", outer = TRUE, line = 0.8 )

# Best in terms of relative bias and relative aboslute error: TBS median Loss and PTBS.seplam median Loss
# TBS Loss median [#1 overall]
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.Loss[,1]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.Loss[,1]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.Loss[,1]-y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,1]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "TBS model for absolute loss, median prediction - best in rel. bias and rel. abs. error", 
         outer = TRUE, line = 0.8 )
  
# PTBS.seplam Loss median [#3 overall]
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.Loss[,3]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.Loss[,3]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.Loss[,3]-y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,3]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "PTBS model with separate lambda for absolute loss - 2nd best in rel. bias and rel. abs. error",
         outer = TRUE, line = 0.8 )

# TBS median DoL [#4 overall] (best in terms of overall ranks)
  par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
  plot( 1:n.obs, (pred.fullCV.DoL[,1]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", 
        mgp = c(2.5,0.8,0), pch = 20, cex = 0.8)
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, 100*((pred.fullCV.DoL[,1]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  plot( 1:n.obs, abs(pred.fullCV.DoL[,1]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", 
        ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,1]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", 
        ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
  abline( h = 0, col = 'blue' )
  abline( h = 100, col = 'blue', lty = 2 )
  title( main = "TBS model for relative loss, median prediction - best in terms of overall ranks", outer = TRUE, line = 0.8 )
dev.off()

# Boxplots
  x11(width = 11, height = 9)
  par( mfrow = c(2,2), mar = c(3,4,0,1), oma = c(5,0,1,0), cex.lab=1.1, cex.axis=0.9)
  boxplot( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - 
             y.Loss, ylab = "Bias [CHF]", xaxt="n")
  abline( h = 0, col = 'blue' )
  lines( apply( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) -
                  y.Loss, 2, mean ), col = 'red', lwd = 1.2, lty = 5 )
  boxplot( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) -
                   y.Loss )/y.Loss, ylab = "Relative bias [%]", xaxt="n" )
  abline( h = 0, col = 'blue' )
  par( mar=c(5,4, 0, 1))
  boxplot( abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) -
                  y.Loss ), ylab = "Mean absolute error [CHF]", las = 3)
  abline( h = 0, col = 'blue' )
  boxplot( 100*abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6])
                    - y.Loss )/y.Loss, ylab = "Relative absolute error [%]", las = 3)
  abline( h = 100, lty = 2 )
  lines( 1:12, c( Delta.CV.Loss[1:3,i], Delta.CV.DoL[1:3,i], Delta.CV.Loss[4:6,i], Delta.CV.DoL[4:6,i] ), col = 'red', 
         lty = 5, lwd = 1.2 )
  savePlot( paste( "Figures/CV_", modelname, "_boxplots", sep = '' ), type = "pdf" )

