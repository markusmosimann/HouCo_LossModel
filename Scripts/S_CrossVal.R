# K = 3: m <- c( 128, 128, 127 )
# K = 5: m <- c( rep(77,3), rep(76,2) )
# K = 8: m <- c( rep(48,7), 47 )
# K = 10: m <- c( rep(39,3), rep(38,7) )
# K = 15: m <- c( rep(26,8), rep(25,7) )
# K = 20: m <- c( rep(20,3), rep(19,17) )
# K = 50: m <- c( rep(8,33), rep(7,17) )
# K = 100: m <- c( rep(4,83), rep(3,17) )

setwd("../MOBILAB/Statistik-Hilfe GIUB/Markus_VulnerabilitätFahrhabe")
# load("Working_20180724.RData")
load("Loss/Loss_Data.RData")
# contains 384 by 5 dataframes "Contents" and "Structure"

source("Global_LossAnalysisLF.R")

idx.outl <- 145

y.Loss <- Contents$Loss[-idx.outl]
x.Loss <- Structure$Loss[-idx.outl]

n.obs <- length(y.Loss)
p <- 2

mle.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = FALSE )
mle.seplam.Loss <- PTBS.mle( y = y.Loss, x = x.Loss, interval = range(lam.seq.Loss), sep.lam = TRUE )
mle.other.Loss <- TBS.mle( y = y.Loss, x = x.Loss, beta.init = c(1000,0.5), interval = c(-3,1) )

y.DoL <- Contents$DoL[-idx.outl]
x.DoL <- Structure$DoL[-idx.outl]

y.InSum <- Contents$InSum[-idx.outl]

mle.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = FALSE )
mle.seplam.DoL <- PTBS.mle( y = y.DoL, x = x.DoL, sep.lam = TRUE )
mle.other.DoL <- TBS.mle( y = y.DoL, x = x.DoL, beta.init = c(0.1,1) )


err.measures <- function( pred, target ) {
	diffs <- pred - target
	c( rmspe = sqrt( mean( diffs^2 ) ), bias = mean( diffs ), rel.tot.bias = 100 * mean(diffs)/mean(target), rel.bias = 100 * ( mean( pred/target ) - 1 ), mae = mean( abs(diffs) ), rel.tot.ae = 100 * mean( abs(diffs) )/mean(target), rel.ae = 100 * mean( abs( 1 - pred/target ) ) )  # all relative quantities are in percent
	# mspe = mean( diffs^2 ) is numerically too large
}

err.measures.med <- function( pred, target ) {
	diffs <- pred - target
	c( med.bias = median( diffs ), med.rel.bias = 100 * median( diffs/target ), med.ae = median( abs(diffs) ), med.rel.ae = 100 * median( abs( diffs/target ) ) )  # all relative quantities are in percent
}

pred.models.comb <- function( x, mle, mle.seplam, mle.other, type ) {
	### TBS model
	fitted.TBS <- c( cbind( 1, x ) %*% mle.other$beta.hat )
	# Median prediction
	TBS.med <- fitted.TBS
	# Mean prediction
	TBS.mean <- TBS.condmean.orig( theta = unlist(mle.other[c(2,3,1)]), x.beta = fitted.TBS )
	names(TBS.mean) <- NULL
	### PTBS model
	# Median prediction:
	PTBS.med <- BC.backtransform( mle$lambda.hat[1], c( cbind( 1, BC.transform( mle$lambda.hat[1+0], x ) ) %*% mle$beta.hat ) )
	# Mean prediction:
	PTBS.mean <- condmean.orig( theta = unlist(mle[c(2,3,1)]), x = cbind( 1, BC.transform( mle$lambda.hat[1+0], x ) ), sep.lam = FALSE )
	names(PTBS.mean) <- NULL
	## Separate lambda
	# Median prediction:
	PTBS.seplam.med <- BC.backtransform( mle.seplam$lambda.hat[1], c( cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ) %*% mle.seplam$beta.hat ) )
	# Mean prediction:
	PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.seplam[c(2,3,1)]), x = cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ), sep.lam = TRUE )
	names(PTBS.seplam.mean) <- NULL	
	# Predictions (matrix n.obs by |{error measures}| )
	ret <- cbind( TBS.med, PTBS.med, PTBS.seplam.med, TBS.mean, PTBS.mean, PTBS.seplam.mean )
	#cbind( Loss.TBS.med = TBS.med, Loss.PTBS.med = PTBS.med, Loss.PTBS.seplam.med = PTBS.seplam.med, Loss.TBS.mean = TBS.mean, Loss.PTBS.mean = PTBS.mean, Loss.PTBS.seplam.mean = PTBS.seplam.mean )
	colnames(ret) <- paste( type, colnames(ret), sep = '.' )
	ret
}

pbias.models.comb.transf <- function( x, mle, mle.seplam, mle.other, y, type ) {
	### TBS model
	pbias.TBS.tr <- BC.transform( mle.other$lambda.hat, c( cbind( 1, x ) %*% mle.other$beta.hat ) ) - BC.transform( mle.other$lambda.hat, y )
	### PTBS model
	pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle$lambda.hat[1+0], x ) ) %*% mle$beta.hat ) - BC.transform( mle$lambda.hat[1], y )
	## Separate lambda
	pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.seplam$lambda.hat[1+1], x ) ) %*% mle.seplam$beta.hat ) - BC.transform( mle.seplam$lambda.hat[1], y )
	# Predictions (matrix n.obs by |{error measures}| )
	ret <- cbind( pbias.TBS.tr, pbias.PTBS.tr, pbias.PTBS.seplam.tr )
	colnames(ret) <- paste( type, c( "TBS.tr", "PTBS.tr", "PTBS.seplam.tr" ), sep = '.' )
	ret
}


## Mean approximation by Taylor series more or less ok (i.e. these quantities are rather small):
summary( mle.Loss$lambda.hat * sqrt(s2.hat.Loss)/( 1 + mle.Loss$lambda.hat * c( cbind( 1, BC.transform( mle.Loss$lambda.hat, x.Loss ) ) %*% mle.Loss$beta.hat ) ) )
summary( mle.seplam.Loss$lambda.hat[1] * sqrt(s2.hat.seplam.Loss)/( 1 + mle.seplam.Loss$lambda.hat[1] * c( cbind( 1, BC.transform( mle.seplam.Loss$lambda.hat[2], x.Loss ) ) %*% mle.seplam.Loss$beta.hat ) ) )
summary( mle.other.Loss$lambda.hat * sqrt(s2.hat.other.Loss)/c( cbind( 1, x.Loss ) %*% mle.other.Loss$beta.hat )^mle.other.Loss$lambda.hat )

summary( mle.DoL$lambda.hat * sqrt(s2.hat.DoL)/( 1 + mle.DoL$lambda.hat * c( cbind( 1, BC.transform( mle.DoL$lambda.hat, x.DoL ) ) %*% mle.DoL$beta.hat ) ) )
summary( mle.seplam.DoL$lambda.hat[1] * sqrt(s2.hat.seplam.DoL)/( 1 + mle.seplam.DoL$lambda.hat[1] * c( cbind( 1, BC.transform( mle.seplam.DoL$lambda.hat[2], x.DoL ) ) %*% mle.seplam.DoL$beta.hat ) ) )
summary( mle.other.DoL$lambda.hat * sqrt(s2.hat.other.DoL)/c( cbind( 1, x.DoL ) %*% mle.other.DoL$beta.hat )^mle.other.DoL$lambda.hat )


##########
# Absolute loss
##########

### K-fold cross-validation (!! long computation time --> I have saved the results)
set.seed(9343029)
seeds <- sample( 1:6000000, 5000 )
K.vec <- c(3,5,8,10,15,20,50,100)
names(K.vec) <- K.vec
nb.K <- length(K.vec)
m.list <- list( '3' = c( 128, 128, 127 ), '5' = c( rep(77,3), rep(76,2) ), '8' = c( rep(48,7), 47 ), '10' = c( rep(39,3), rep(38,7) ), '15' = c( rep(26,8), rep(25,7) ), '20' = c( rep(20,3), rep(19,17) ), '50' = c( rep(8,33), rep(7,17) ), '100' = c( rep(4,83), rep(3,17) ) )
n.runs <- 100


pred.array.Loss <- array( dim = c( nb.K, n.obs, 6, n.runs ), dimnames = list( names(K.vec), NULL, c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" ), NULL ) )

for( K in K.vec ) {
	print( paste( "K =", K ) )
	K.nam <- as.character(K)
	it <- which(K.vec %in% K.nam)
	# n.runs * K groups, of sizes m.list[[]]
	for( r in 1:n.runs ) {
		set.seed(seeds[(it-1)*n.runs + r]) # avoid same seed for all K
		groups <- split( sample( n.obs, n.obs, replace = FALSE ), f = rep( 1:K, m.list[[K.nam]] ) )
		for( k in 1:K ) {
		#=======#
			# Fit model to y[-groups[[k]]], giving estimated parameters theta.fit.
			### TBS model
			mle.k <- TBS.mle( y = y.Loss[-groups[[k]]], x = x.Loss[-groups[[k]]], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) ) 
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			fitted.k <- c( cbind( 1, x.Loss[groups[[k]]] ) %*% mle.k$beta.hat )
			TBS.med <- fitted.k # predicted median
			TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) # predicted mean
			names(TBS.mean) <- NULL
			### PTBS model
			mle.k <- PTBS.mle( y = y.Loss[-groups[[k]]], x = x.Loss[-groups[[k]]], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[groups[[k]]] ) ) %*% mle.k$beta.hat ) ) # median prediction
			PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[groups[[k]]] ) ), sep.lam = FALSE ) # mean prediction
			names(PTBS.mean) <- NULL
			## Separate lambda
			mle.k <- PTBS.mle( y = y.Loss[-groups[[k]]], x = x.Loss[-groups[[k]]], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[groups[[k]]] ) ) %*% mle.k$beta.hat ) ) # median prediction
			PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[groups[[k]]] ) ), sep.lam = TRUE ) # mean prediction
			names(PTBS.seplam.mean) <- NULL
		#=======#
			pred.array.Loss[it,groups[[k]],,r] <- cbind( Loss.TBS.med = TBS.med, Loss.PTBS.med = PTBS.med, Loss.PTBS.seplam.med = PTBS.seplam.med, Loss.TBS.mean = TBS.mean, Loss.PTBS.mean = PTBS.mean, Loss.PTBS.seplam.mean = PTBS.seplam.mean )
		}
	}
}
pred.err.Kfold.Loss <- aperm( apply( pred.array.Loss, c(1,3:4), err.measures, target = y.Loss ), c(2,1,3:4) )

dimnames( pred.err.Kfold.Loss )[[1]] <- names(K.vec)
dimnames( pred.err.Kfold.Loss )[[3]] <- c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" )

dimnames( pred.array.Loss )[[1]] <- names(K.vec)
dimnames( pred.array.Loss )[[3]] <- c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" )

save( list = c( "pred.array.Loss", "pred.err.Kfold.Loss" ), file = "pred_KfoldCV_Loss.RData" )

par( mfrow = c(2,4) )
for( i in c(2:4,1,5:7) ) {
	boxplot( t( pred.err.Kfold.Loss[,i,1,] ), ylab = err.name[i], xlab = "K-fold", main = "Absolute loss" )
	abline( h = Delta.CV.Loss[1,i], col = 'blue' )
}
# variance not increasing with K

pred.err.byfold.Loss <- array( dim = c( nb.K, 7, 6, n.runs ), dimnames = list( names(K.vec), NULL, c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" ), NULL )  )
# pred.array.Loss: c( nb.K, n.obs, 6, n.runs )
for( K in K.vec ) {
	K.nam <- as.character(K)
	it <- which(K.vec %in% K.nam)
	# n.runs * K groups, of sizes m.list[[]]
	for( r in 1:n.runs ) {
		#print( paste( "n.runs =", r ) ) 
		set.seed(seeds[(it-1)*n.runs + r]) # avoid same seed for all K
		groups <- split( sample( n.obs, n.obs, replace = FALSE ), f = rep( 1:K, m.list[[K.nam]] ) )
		pred.err.byfold.Loss[it,,,r] <- apply( vapply( X = groups, FUN = function(idx){ apply( pred.array.Loss[it,idx,,r], 2, err.measures, target = y.Loss[idx] ) }, FUN.VALUE = array( 0.1, dim = c(7,6) ) ), 1:2, mean )
	}
}

dev.set(2)
par( mfrow = c(2,4) )
for( i in c(2:4,1,5:7) ) {
	boxplot( t( pred.err.byfold.Loss[,i,1,] ), ylab = err.name[i], xlab = "K-fold", main = "Absolute loss" )
	abline( h = Delta.CV.Loss[1,i], col = 'blue' )
}
# now the variance increases with K, but this is not consistent with single value for K = n.obs!!
###

### Leave-one-out cross validation
# Predictions (matrix n.obs by |{model:pred}|)
pred.fullCV.Loss <- t( sapply( 1:n.obs, function(k){
	# Fit model to y[-k], giving estimated parameters theta.fit.
	### TBS model
	mle.k <- TBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
	# Predict y[k] by yhat( theta.fit, x[k] ).
	fitted.k <- c( cbind( 1, x.Loss[k] ) %*% mle.k$beta.hat )
	# Median prediction
	TBS.med <- fitted.k
	# Mean prediction
	TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k )
	names(TBS.mean) <- NULL
	### PTBS model
	mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
	# Median prediction:
	PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ) %*% mle.k$beta.hat ) )
	# Mean prediction:
	PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ), sep.lam = FALSE )
	names(PTBS.mean) <- NULL
	## Separate lambda
	mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
	# Median prediction:
	PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[k] ) ) %*% mle.k$beta.hat ) )
	# Mean prediction:
	PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[k] ) ), sep.lam = TRUE )
	names(PTBS.seplam.mean) <- NULL
	c( Loss.TBS.med = TBS.med, Loss.PTBS.med = PTBS.med, Loss.PTBS.seplam.med = PTBS.seplam.med, Loss.TBS.mean = TBS.mean, Loss.PTBS.mean = PTBS.mean, Loss.PTBS.seplam.mean = PTBS.seplam.mean )
} ) )

# Aggregate prediction errors (|{model:pred}| by |{err.measures}|)
Delta.CV.Loss <- t( apply( pred.fullCV.Loss, 2, err.measures, target = y.Loss ) )
Delta.CV.med.Loss <- t( apply( pred.fullCV.Loss, 2, err.measures.med, target = y.Loss ) )
###

### Apparent error (prediction of used data under model):
pred.mod.Loss <- pred.models.comb( x = x.Loss, mle = mle.Loss, mle.seplam = mle.seplam.Loss, mle.other = mle.other.Loss, type = "Loss" ) # n.obs by |{model:pred}|

# Aggregate prediction errors
Delta.app.Loss <- t( apply( pred.mod.Loss, 2, err.measures, target = y.Loss ) )
Delta.app.med.Loss <- t( apply( pred.mod.Loss, 2, err.measures.med, target = y.Loss ) )
###


### Bootstrap
R <- 500
set.seed(572042)
resamples <- sapply( 1:R, function(r){ sample( n.obs, n.obs, replace = TRUE ) } ) # n.obs by R
j.in <- apply( resamples, 2, function(star){ 1:n.obs %in% star } ) # n.obs by R

# Prediction of all data points under model fits from resamples (n.obs by |{model:pred}| by R)
pred.resamp.Loss <- array( apply( resamples, 2, function(star){ 
	### TBS model
	mle.other.star <- TBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
	### PTBS model
	mle.star <- PTBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
	## Separate lambda
	mle.seplam.star <- PTBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )

	pred.models.comb( x = x.Loss, mle = mle.star, mle.seplam = mle.seplam.star, mle.other = mle.other.star, type = "Loss" )
} ), dim = c(n.obs,6,R), dimnames = list( NULL, paste( "Loss", c( "TBS.med", "PTBS.med", "PTBS.seplam.med", "TBS.mean", "PTBS.mean", "PTBS.seplam.mean" ), sep = '.' ), NULL ) )
# Check that array is correctly constructed:
# array( apply( matrix( 1:8, nr = 4, nc = 2 ), 2, function(x){ cbind( x, x^2, x^3 ) } ), dim = c(4,3,2), dimnames = list( NULL, c('A','B','C'), NULL ) )

# Prediction errors wrt observed data (|{model:pred}| by |{error measures}| by R)
pred.err.resamp.Loss <- aperm( apply( pred.resamp.Loss, 2:3, err.measures, target = y.Loss ), c(2,1,3) )
pred.err.med.resamp.Loss <- aperm( apply( pred.resamp.Loss, 2:3, err.measures.med, target = y.Loss ), c(2,1,3) )

# Prediction errors wrt resampled data (|{model:pred}| by |{error measures}| by R)
pred.mod.err.resamp.Loss <- vapply( 1:R, function(r){ t( apply( pred.resamp.Loss[resamples[,r],,r], 2, err.measures, target = y.Loss[resamples[,r]] ) ) }, FUN.VALUE = array( 0.1, dim = c(6,7) ) )
pred.mod.err.med.resamp.Loss <- vapply( 1:R, function(r){ t( apply( pred.resamp.Loss[resamples[,r],,r], 2, err.measures.med, target = y.Loss[resamples[,r]] ) ) }, FUN.VALUE = array( 0.1, dim = c(6,4) ) )

# Bootstrap aggregate prediction error (|{model:pred}| by |{error measures}|)
Delta.B.Loss <- apply( pred.err.resamp.Loss - pred.mod.err.resamp.Loss, 1:2, mean ) + Delta.app.Loss
Delta.B.med.Loss <- apply( pred.err.med.resamp.Loss - pred.mod.err.med.resamp.Loss, 1:2, mean ) + Delta.app.med.Loss

## Two different interpretations of leave-one-out bootstrap error:
# (1) Average of prediction errors for resamples without j:
pred.err.BCV.Loss <- vapply( 1:R, function(r){ t( apply( pred.resamp.Loss[!j.in[,r],,r], 2, err.measures, target = y.Loss[!j.in[,r]] ) ) }, FUN.VALUE = matrix( 0.1, nr = 6, nc = 7 ) ) # |{model:pred}| by |{error measures}| by R
pred.err.med.BCV.Loss <- vapply( 1:R, function(r){ t( apply( pred.resamp.Loss[!j.in[,r],,r], 2, err.measures.med, target = y.Loss[!j.in[,r]] ) ) }, FUN.VALUE = matrix( 0.1, nr = 6, nc = 4 ) ) # |{model:pred}| by |{error measures}| by R
Delta.BCV1.Loss <- apply( pred.err.BCV.Loss, 1:2, mean )
Delta.BCV1.med.Loss <- apply( pred.err.med.BCV.Loss, 1:2, mean )

# (2) Error of average predictions over resamples without j
# Average prediction for each j over resamples without j:
pred.smooth.out.Loss <- t( sapply( 1:n.obs, function(j){ apply( pred.resamp.Loss[j,,!j.in[j,]], 1, mean ) } ) ) # n.obs by |{model:pred}|
Delta.BCV2.Loss <- t( apply( pred.smooth.out.Loss, 2, err.measures, target = y.Loss ) )
Delta.BCV2.med.Loss <- t( apply( pred.smooth.out.Loss, 2, err.measures.med, target = y.Loss ) )

Delta.BCV1.Loss/Delta.BCV2.Loss # largest relative discrepancies for bias and rel.tot.bias for mean prediction
Delta.BCV1.med.Loss/Delta.BCV2.med.Loss # largest relative discrepancies for bias and relative bias under median prediction

Delta.632.Loss <- (1 - exp(-1)) * Delta.BCV2.Loss + exp(-1) * Delta.app.Loss
Delta.632.med.Loss <- (1 - exp(-1)) * Delta.BCV2.med.Loss + exp(-1) * Delta.app.med.Loss


## Approach not useful/necessary(?) because R either infinite or small
gamma.Loss <- t( apply( pred.mod.Loss, 2, function(pred){ err.measures( pred[rep( 1:n.obs, n.obs )], target = y.Loss[rep( 1:n.obs, each = n.obs )] ) } ) )

R.Loss <- (Delta.BCV1.Loss - Delta.app.Loss)/(gamma.Loss - Delta.app.Loss)
w.hat.Loss <- (1 - exp(-1))/(1 - exp(-1)*R.Loss)

Delta.632plus.Loss <- w.hat.Loss * Delta.BCV2.Loss + (1 - w.hat.Loss) * Delta.app.Loss
###

ranking <- order(y.Loss)

dev.set(3)
# Individual biases for mean prediction (resamples do or don't contain obs j)
par( mfcol = c(3,2) )
for(i in 1:6) {
	boxplot( t( (pred.resamp.Loss - y.Loss)[ranking,i,] ), outlty = 1, outpch = NA, xlab = "Rank of target value", ylab = "Prediction bias [CHF]", main = dimnames(pred.resamp.Loss)[[2]][i] )
	points( 1:n.obs, (pred.fullCV.Loss[,i] - y.Loss)[ranking], col = 'red', pch = 3, cex = 0.8 )
	abline( h = seq(-140000,320000,10000), col = 'grey', lty = 3 )
	abline( h = 0, col = 'blue' )
	rug( which( apply( pred.resamp.Loss - y.Loss, 1:2, function(e){ min(e) <= 0 & max(e) >= 0 } )[ranking,i] ), col = 'green3' )
	box()
}
# Only resamples without obs j: 
#sapply( 1:n.obs, function(j){ (pred.resamp.Loss[j,i,!j.in[j,]] - y.Loss[j]) } )[ranking]

apply( apply( pred.resamp.Loss - y.Loss, 1:2, function(e){ min(e) > 0 | max(e) < 0 } )[,4:6], 2, sum ) # oops, over 300 predictions are biased!!!


plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Contents", main = "Absolute loss, original scale" )
points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4 ) # ev. pch = 8
lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ), col = 'red' ) # fitted median 
points( x.Loss, pred.fullCV.Loss[,2], col = 'red', pch = 3 )
lines( x.seq.Loss, condmean.seq.Loss, col = 'green3', lwd = 1 )
points( x.Loss, pred.fullCV.Loss[,5], col = 'green3', pch = 3 )
lines( x.seq.Loss, c( cbind(1,x.seq.Loss) %*% mle.other.Loss$beta.hat ), col = 'blue' )
points( x.Loss, pred.fullCV.Loss[,1], col = 'blue', pch = 3 )

points( x.Loss, pred.fullCV.Loss[,4], col = 'purple', pch = 3 )
points( x.Loss, pred.fullCV.Loss[,3], col = 'salmon', pch = 4 )
points( x.Loss, pred.fullCV.Loss[,6], col = 'skyblue', pch = 4 )


#====
# Prediction on transformed scale
#====

# Prediction bias on transformed scale (same for mean an median), n.obs by |models| by R
pbias.transf.resamp.Loss <- array( apply( resamples, 2, function(star){ 
	### TBS model
	mle.other.star <- TBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
	### PTBS model
	mle.star <- PTBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
	## Separate lambda
	mle.seplam.star <- PTBS.mle( y = y.Loss[star], x = x.Loss[star], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
	
	pbias.models.comb.transf( x = x.Loss, mle = mle.star, mle.seplam = mle.seplam.star, mle.other = mle.other.star, y = y.Loss, type = "Loss" )
} ), dim = c(n.obs,3,R), dimnames = list( NULL, paste( "Loss", c( "TBS.tr", "PTBS.tr", "PTBS.seplam.tr" ), sep = '.' ), NULL ) )

pbias.transf.fullCV.Loss <- t( sapply( 1:n.obs, function(k){
	# Fit model to y[-k], giving estimated parameters theta.fit.
	### TBS model
	mle.k <- TBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, beta.init = c(1000,0.5), interval = c(-3,1) )
	pbias.TBS.tr <- BC.transform( mle.k$lambda.hat, c( cbind( 1, x.Loss[k] ) %*% mle.k$beta.hat ) ) - BC.transform( mle.k$lambda.hat, y.Loss[k] )
	### PTBS model
	mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = FALSE )
	pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.Loss[k] ) ) %*% mle.k$beta.hat ) - BC.transform( mle.k$lambda.hat[1], y.Loss[k] )
	## Separate lambda
	mle.k <- PTBS.mle( y = y.Loss[-k], x = x.Loss[-k], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
	pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.Loss[k] ) ) %*% mle.k$beta.hat ) - BC.transform( mle.k$lambda.hat[1], y.Loss[k] )
	cbind( CV.TBS.tr.Loss = pbias.TBS.tr, CV.PTBS.tr.Loss = pbias.PTBS.tr, CV.PTBS.seplam.Loss = pbias.PTBS.seplam.tr )
} ) )

s2.hat.Loss <- mle.Loss$sigma2.hat * n.obs/(n.obs-p-1)
covML.Loss <- solve( -optimHess( unlist(mle.Loss[c(2,3,1)]), PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
var.pred.PTBS.transf.Loss <- sapply( BC.transform( mle.Loss$lambda.hat[1], x.Loss ), function(x){ c(1,x) %*% covML.Loss[1:p,1:p] %*% c(1,x) } ) + s2.hat.Loss

s2.hat.seplam.Loss <- mle.seplam.Loss$sigma2.hat * n.obs/(n.obs-p-2)
covML.seplam.Loss <- solve( -optimHess( unlist(mle.seplam.Loss[c(2,3,1)]), PTBS.llkhd, y = y.Loss, x = x.Loss, sep.lam = TRUE, control = list( maxit = 5000, fnscale = -1 ) ) )
var.pred.PTBS.seplam.transf.Loss <- sapply( BC.transform( mle.seplam.Loss$lambda.hat[2], x.Loss ), function(x){ c(1,x) %*% covML.seplam.Loss[1:p,1:p] %*% c(1,x) } ) + s2.hat.seplam.Loss

s2.hat.other.Loss <- mle.other.Loss$sigma2.hat * n.obs/(n.obs-p-1)
covML.other.Loss <- solve( -optimHess( unlist( mle.other.Loss[c(2,3,1)] ), loglik.other, y = y.Loss, x = x.Loss, control = list( fnscale = -1 ) ) )
var.pred.TBS.transf.Loss <- sapply( x.Loss, function(x){ c( c( 1, x ) %*% mle.other.Loss$beta.hat )^(2*(mle.other.Loss$lambda.hat - 1)) * c( c( 1, x ) %*% covML.other.Loss[1:p,1:p] %*% c(1, x) ) } ) + s2.hat.other.Loss
# With regression version of var(beta.hat) (gives similar result):
#sapply( x.Loss, function(x){ c( c( 1, x ) %*% mle.other.Loss$beta.hat )^(2*(mle.other.Loss$lambda.hat - 1)) * c( c( 1, x ) %*% ( s2.hat.other.Loss * solve( crossprod( diag( c( cbind( 1, x.Loss ) %*% mle.other.Loss$beta.hat )^(mle.other.Loss$lambda.hat - 1) ) %*% cbind( 1, x.Loss ) ) ) ) %*% c(1, x) ) } ) + s2.hat.other.Loss

var.pred.transf.Loss <- cbind( TBS = var.pred.TBS.transf.Loss, PTBS = var.pred.PTBS.transf.Loss, PTBS.seplam = var.pred.PTBS.seplam.transf.Loss )

par(mfrow=c(3,1))
for(i in 1:3) {
	boxplot( t( pbias.transf.resamp.Loss[ranking,i,] ), outlty = 1, outpch = NA, xlab = "Rank of target loss", ylab = "Prediction bias [transf.CHF]", main = dimnames(pbias.transf.resamp.Loss)[[2]][i] )
	points( 1:n.obs, pbias.transf.fullCV.Loss[ranking,i], col = 'red', pch = 3, cex = 0.8 )
	lines( -qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.Loss[ranking,i]), col = 'blue', lty = 2 )
	lines( qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.Loss[ranking,i]), col = 'blue', lty = 2 )
#	abline( h = seq(-15,20,5), col = 'grey', lty = 3 )
	abline( h = 0, col = 'blue' )
	rug( which( ( apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss[,i]) | apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss[,i]) )[ranking,i] ), col = 'red' )
	box()
}

apply( apply( pbias.transf.resamp.Loss, 1:2, function(e){ min(e) > 0 | max(e) < 0 } ), 2, sum ) # only very few less!!!
sum( apply( apply( pbias.transf.resamp.Loss, 1:2, function(e){ min(e) <= 0 & max(e) >= 0 } ), 1, all ) )  # predictions for 50 observation are unbiased under all three models

apply( apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss[,i]) | apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss), 2, sum ) # no obs. outside prediction interval
apply( apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss[,i]) | apply( pbias.transf.resamp.Loss, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.Loss), 2, which )
# none
#======


##########
# Relative loss
##########

### K-fold cross-validation (!! long computation time --> I have saved the results)
set.seed(9343029)
seeds <- sample( 1:6000000, 5000 )
K.vec <- c(3,5,8,10,15,20,50,100)
names(K.vec) <- K.vec
nb.K <- length(K.vec)
m.list <- list( '3' = c( 128, 128, 127 ), '5' = c( rep(77,3), rep(76,2) ), '8' = c( rep(48,7), 47 ), '10' = c( rep(39,3), rep(38,7) ), '15' = c( rep(26,8), rep(25,7) ), '20' = c( rep(20,3), rep(19,17) ), '50' = c( rep(8,33), rep(7,17) ), '100' = c( rep(4,83), rep(3,17) ) )
n.runs <- 100


pred.array.DoL <- array( dim = c( nb.K, n.obs, 6, n.runs ), dimnames = list( names(K.vec), NULL, c( "DoL.TBS.med", "DoL.PTBS.med", "DoL.PTBS.seplam.med", "DoL.TBS.mean", "DoL.PTBS.mean", "DoL.PTBS.seplam.mean" ), NULL ) )

for( K in K.vec ) {
	print( paste( "K =", K ) )
	K.nam <- as.character(K)
	it <- which(K.vec %in% K.nam)
	# n.runs * K groups, of sizes m.list[[]]
	for( r in 1:n.runs ) {
		set.seed(seeds[(it-1)*n.runs + r]) # avoid same seed for all K
		groups <- split( sample( n.obs, n.obs, replace = FALSE ), f = rep( 1:K, m.list[[K.nam]] ) )
		for( k in 1:K ) {
		#=======#
			# Fit model to y[-groups[[k]]], giving estimated parameters theta.fit.
			### TBS model
			mle.k <- TBS.mle( y = y.DoL[-groups[[k]]], x = x.DoL[-groups[[k]]], intercept = TRUE, beta.init = c(0.1,1) ) 
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			fitted.k <- c( cbind( 1, x.DoL[groups[[k]]] ) %*% mle.k$beta.hat )
			TBS.med <- fitted.k * y.InSum[groups[[k]]] # predicted median
			TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) * y.InSum[groups[[k]]] # predicted mean
			names(TBS.mean) <- NULL
			### PTBS model
			mle.k <- PTBS.mle( y = y.DoL[-groups[[k]]], x = x.DoL[-groups[[k]]], intercept = TRUE, sep.lam = FALSE )
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[groups[[k]]] ) ) %*% mle.k$beta.hat ) ) * y.InSum[groups[[k]]] # median prediction
			PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[groups[[k]]] ) ), sep.lam = FALSE ) * y.InSum[groups[[k]]] # mean prediction
			names(PTBS.mean) <- NULL
			## Separate lambda
			mle.k <- PTBS.mle( y = y.DoL[-groups[[k]]], x = x.DoL[-groups[[k]]], intercept = TRUE, sep.lam = TRUE )
			# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
			PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[groups[[k]]] ) ) %*% mle.k$beta.hat ) ) * y.InSum[groups[[k]]] # median prediction
			PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[groups[[k]]] ) ), sep.lam = TRUE ) * y.InSum[groups[[k]]] # mean prediction
			names(PTBS.seplam.mean) <- NULL
		#=======#
			pred.array.DoL[it,groups[[k]],,r] <- cbind( DoL.TBS.med = TBS.med, DoL.PTBS.med = PTBS.med, DoL.PTBS.seplam.med = PTBS.seplam.med, DoL.TBS.mean = TBS.mean, DoL.PTBS.mean = PTBS.mean, DoL.PTBS.seplam.mean = PTBS.seplam.mean )
		}
	}
}
pred.err.Kfold.DoL <- aperm( apply( pred.array.DoL, c(1,3:4), err.measures, target = y.DoL * y.InSum ), c(2,1,3:4) )

save( list = c( "pred.array.DoL", "pred.err.Kfold.DoL" ), file = "pred_KfoldCV_DoL.RData" )

par( mfrow = c(2,4) )
for( i in c(2:4,1,5:7) ) {
	boxplot( t( pred.err.Kfold.DoL[,i,1,] ), ylab = err.name[i], xlab = "K-fold", main = "Relative loss" )
	abline( h = Delta.CV.DoL[1,i], col = 'blue' )
}
# variance not increasing with K

pred.err.byfold.DoL <- array( dim = c( nb.K, 7, 6, n.runs ), dimnames = list( names(K.vec), NULL, c( "DoL.TBS.med", "DoL.PTBS.med", "DoL.PTBS.seplam.med", "DoL.TBS.mean", "DoL.PTBS.mean", "DoL.PTBS.seplam.mean" ), NULL )  )
# pred.array.DoL: c( nb.K, n.obs, 6, n.runs )
for( K in K.vec ) {
	K.nam <- as.character(K)
	it <- which(K.vec %in% K.nam)
	# n.runs * K groups, of sizes m.list[[]]
	for( r in 1:n.runs ) {
		#print( paste( "n.runs =", r ) ) 
		set.seed(seeds[(it-1)*n.runs + r]) # avoid same seed for all K
		groups <- split( sample( n.obs, n.obs, replace = FALSE ), f = rep( 1:K, m.list[[K.nam]] ) )
		pred.err.byfold.DoL[it,,,r] <- apply( vapply( X = groups, FUN = function(idx){ apply( pred.array.DoL[it,idx,,r], 2, err.measures, target = y.DoL[idx] * y.InSum[idx] ) }, FUN.VALUE = array( 0.1, dim = c(7,6) ) ), 1:2, mean )
	}
}

dev.set(2)
par( mfrow = c(2,4) )
for( i in c(2:4,1,5:7) ) {
	boxplot( t( pred.err.byfold.DoL[,i,1,] ), ylab = err.name[i], xlab = "K-fold", main = "Relative loss" )
	abline( h = Delta.CV.DoL[1,i], col = 'blue' )
}
# now the variance increases with K, but this is not consistent with single value for K = n.obs!!
###

### Leave-one-out cross validation
# Predictions (matrix n.obs by |{model:pred}|)
pred.fullCV.DoL <- t( sapply( 1:n.obs, function(k){
	# Fit model to y[-k], giving estimated parameters theta.fit.
	### TBS model
	mle.k <- TBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, beta.init = c(0.1,1) ) 
	# Predict y[k] by yhat( theta.fit, x[k] ).
	fitted.k <- c( cbind( 1, x.DoL[k] ) %*% mle.k$beta.hat )
	# Median prediction
	TBS.med <- fitted.k * y.InSum[k]
	# Mean prediction
	TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) * y.InSum[k]
	names(TBS.mean) <- NULL
	### PTBS model
	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = FALSE )
	# Median prediction:
	PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ) %*% mle.k$beta.hat ) ) * y.InSum[k]
	# Mean prediction:
	PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ), sep.lam = FALSE ) * y.InSum[k]
	names(PTBS.mean) <- NULL
	## Separate lambda
	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = TRUE )
	# Median prediction:
	PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ) %*% mle.k$beta.hat ) ) * y.InSum[k]
	# Mean prediction:
	PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ), sep.lam = TRUE ) * y.InSum[k]
	names(PTBS.seplam.mean) <- NULL
	c( DoL.TBS.med = TBS.med, DoL.PTBS.med = PTBS.med, DoL.PTBS.seplam.med = PTBS.seplam.med, DoL.TBS.mean = TBS.mean, DoL.PTBS.mean = PTBS.mean, DoL.PTBS.seplam.mean = PTBS.seplam.mean )
} ) ) # n.obs by 6

# Aggregate prediction errors (|{model:pred}| by |{err.measures}|)
Delta.CV.DoL <- t( apply( pred.fullCV.DoL, 2, err.measures, target = y.DoL * y.InSum ) )
Delta.CV.med.DoL <- t( apply( pred.fullCV.DoL, 2, err.measures.med, target = y.DoL * y.InSum ) )
###

### Apparent error (prediction of used data under model):
pred.mod.DoL <- pred.models.comb( x = x.DoL, mle = mle.DoL, mle.seplam = mle.seplam.DoL, mle.other = mle.other.DoL, type = "DoL" ) * y.InSum  # n.obs by |{model:pred}|

# Aggregate prediction errors
Delta.app.DoL <- t( apply( pred.mod.DoL, 2, err.measures, target = y.DoL * y.InSum ) )
Delta.app.med.DoL <- t( apply( pred.mod.DoL, 2, err.measures.med, target = y.DoL * y.InSum ) )
###


### Bootstrap
R <- 500
set.seed(572042)
resamples <- sapply( 1:R, function(r){ sample( n.obs, n.obs, replace = TRUE ) } ) # n.obs by R
j.in <- apply( resamples, 2, function(star){ 1:n.obs %in% star } ) # n.obs by R

# Prediction of all data points under model fits from resamples (n.obs by |{model:pred}| by R)
pred.resamp.DoL <- array( apply( resamples, 2, function(star){ 
	### TBS model
	mle.other.star <- TBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, beta.init = c(0.1,1) ) 
	### PTBS model
	mle.star <- PTBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, sep.lam = FALSE )
	## Separate lambda
	mle.seplam.star <- PTBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, sep.lam = TRUE )

	pred.models.comb( x = x.DoL, mle = mle.star, mle.seplam = mle.seplam.star, mle.other = mle.other.star, type = "DoL" ) * y.InSum
} ), dim = c(n.obs,6,R), dimnames = list( NULL, paste( "DoL", c( "TBS.med", "PTBS.med", "PTBS.seplam.med", "TBS.mean", "PTBS.mean", "PTBS.seplam.mean" ), sep = '.' ), NULL ) )

# Prediction errors wrt observed data (|{model:pred}| by |{error measures}| by R)
pred.err.resamp.DoL <- aperm( apply( pred.resamp.DoL, 2:3, err.measures, target = y.DoL * y.InSum ), c(2,1,3) )
pred.err.med.resamp.DoL <- aperm( apply( pred.resamp.DoL, 2:3, err.measures.med, target = y.DoL * y.InSum ), c(2,1,3) )

# Prediction errors wrt resampled data (|{model:pred}| by |{error measures}| by R)
pred.mod.err.resamp.DoL <- vapply( 1:R, function(r){ t( apply( pred.resamp.DoL[resamples[,r],,r], 2, err.measures, target = (y.DoL * y.InSum)[resamples[,r]] ) ) }, FUN.VALUE = array( 0.1, dim = c(6,7) ) )
pred.mod.err.med.resamp.DoL <- vapply( 1:R, function(r){ t( apply( pred.resamp.DoL[resamples[,r],,r], 2, err.measures.med, target = (y.DoL * y.InSum)[resamples[,r]] ) ) }, FUN.VALUE = array( 0.1, dim = c(6,4) ) )

# Bootstrap aggregate prediction error (|{model:pred}| by |{error measures}|)
Delta.B.DoL <- apply( pred.err.resamp.DoL - pred.mod.err.resamp.DoL, 1:2, mean ) + Delta.app.DoL
Delta.B.med.DoL <- apply( pred.err.med.resamp.DoL - pred.mod.err.med.resamp.DoL, 1:2, mean ) + Delta.app.med.DoL


## Two different interpretations of leave-one-out bootstrap error:
# (1) Average of prediction errors for resamples without j:
pred.err.BCV.DoL <- vapply( 1:R, function(r){ t( apply( pred.resamp.DoL[!j.in[,r],,r], 2, err.measures, target = (y.DoL * y.InSum)[!j.in[,r]] ) ) }, FUN.VALUE = matrix( 0.1, nr = 6, nc = 7 ) ) # |{model:pred}| by |{error measures}| by R
pred.err.med.BCV.DoL <- vapply( 1:R, function(r){ t( apply( pred.resamp.DoL[!j.in[,r],,r], 2, err.measures.med, target = (y.DoL * y.InSum)[!j.in[,r]] ) ) }, FUN.VALUE = matrix( 0.1, nr = 6, nc = 4 ) ) # |{model:pred}| by |{error measures med}| by R
Delta.BCV1.DoL <- apply( pred.err.BCV.DoL, 1:2, mean )
Delta.BCV1.med.DoL <- apply( pred.err.med.BCV.DoL, 1:2, mean )

# (2) Error of average predictions over resamples without j
# Average prediction for each j over resamples without j:
pred.smooth.out.DoL <- t( sapply( 1:n.obs, function(j){ apply( pred.resamp.DoL[j,,!j.in[j,]], 1, mean ) } ) ) # n.obs by |{model:pred}|
Delta.BCV2.DoL <- t( apply( pred.smooth.out.DoL, 2, err.measures, target = y.DoL * y.InSum ) )
Delta.BCV2.med.DoL <- t( apply( pred.smooth.out.DoL, 2, err.measures.med, target = y.DoL * y.InSum ) )

Delta.BCV1.DoL/Delta.BCV2.DoL # largest relative discrepancies for bias and rel.tot.bias for mean prediction

Delta.BCV1.med.DoL/Delta.BCV2.med.DoL # larger differences for med.bias and med.rel.bias

Delta.632.DoL <- (1 - exp(-1)) * Delta.BCV2.DoL + exp(-1) * Delta.app.DoL # largest discrepancies for bias and rel.tot.bias under PTBS mean estimation
Delta.632.med.DoL <- (1 - exp(-1)) * Delta.BCV2.med.DoL + exp(-1) * Delta.app.med.DoL
###

## Approach not useful/necessary(?) because R either infinite or small
gamma.DoL <- t( apply( pred.mod.DoL, 2, function(pred){ err.measures( pred[rep( 1:n.obs, n.obs )], target = (y.DoL * y.InSum)[rep( 1:n.obs, each = n.obs )] ) } ) )

R.DoL <- (Delta.BCV1.DoL - Delta.app.DoL)/(gamma.DoL - Delta.app.DoL)
w.hat.DoL <- (1 - exp(-1))/(1 - exp(-1)*R.DoL)

Delta.632plus.DoL <- w.hat.DoL * Delta.BCV2.DoL + (1 - w.hat.DoL) * Delta.app.DoL
###


# 
# Plots of individual bias on original scale
dev.set(3)
# Individual biases for mean prediction (resamples do or don't contain obs j)
par( mfcol = c(3,2) )
for(i in 1:6) {
	boxplot( t( (pred.resamp.DoL - y.DoL * y.InSum)[ranking,i,] ), outlty = 1, outpch = NA, xlab = "Rank of target value", ylab = "Prediction bias [CHF]", main = dimnames(pred.resamp.DoL)[[2]][i] )
	points( 1:n.obs, (pred.fullCV.DoL[,i] - y.DoL * y.InSum)[ranking], col = 'red', pch = 3, cex = 0.8 )
	abline( h = seq(-140000,320000,10000), col = 'grey', lty = 3 )
	abline( h = 0, col = 'blue' )
	rug( which( apply( pred.resamp.DoL - y.DoL * y.InSum, 1:2, function(e){ min(e) <= 0 & max(e) >= 0 } )[ranking,i] ), col = 'green3' )
	box()
}
# Only resamples without obs j: 
#sapply( 1:n.obs, function(j){ pred.resamp.DoL[j,i,!j.in[j,]] - y.DoL[j] * y.InSum[j] } )[ranking]

apply( apply( pred.resamp.DoL - y.DoL * y.InSum, 1:2, function(e){ min(e) > 0 | max(e) < 0 } )[,4:6], 2, sum ) # oops, over 300 predictions are biased!!!


#====
# Prediction on transformed scale
#====

# Prediction bias on transformed scale (same for mean and median), n.obs by |models| by R
pbias.transf.resamp.DoL <- array( apply( resamples, 2, function(star){ 
	### TBS model
	mle.other.star <- TBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, beta.init = c(0.1,1) )
	### PTBS model
	mle.star <- PTBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, sep.lam = FALSE )
	## Separate lambda
	mle.seplam.star <- PTBS.mle( y = y.DoL[star], x = x.DoL[star], intercept = TRUE, sep.lam = TRUE )
	
	pbias.models.comb.transf( x = x.DoL, mle = mle.star, mle.seplam = mle.seplam.star, mle.other = mle.other.star, y = y.DoL, type = "DoL" )
} ), dim = c(n.obs,3,R), dimnames = list( NULL, paste( "DoL", c( "TBS.tr", "PTBS.tr", "PTBS.seplam.tr" ), sep = '.' ), NULL ) ) # on transformed scale

pbias.transf.fullCV.DoL <- t( sapply( 1:n.obs, function(k){
	# Fit model to y[-k], giving estimated parameters theta.fit.
	### TBS model
	mle.k <- TBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, beta.init = c(0.1,1) )
	pbias.TBS.tr <- BC.transform( mle.k$lambda.hat, c( cbind( 1, x.DoL[k] ) %*% mle.k$beta.hat ) ) - BC.transform( mle.k$lambda.hat, y.DoL[k] )
	### PTBS model
	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = FALSE )
	pbias.PTBS.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.DoL[k] ) ) %*% mle.k$beta.hat ) - BC.transform( mle.k$lambda.hat[1], y.DoL[k] )
	## Separate lambda
	mle.k <- PTBS.mle( y = y.DoL[-k], x = x.DoL[-k], intercept = TRUE, sep.lam = TRUE )
	pbias.PTBS.seplam.tr <- c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.DoL[k] ) ) %*% mle.k$beta.hat ) - BC.transform( mle.k$lambda.hat[1], y.DoL[k] )
	cbind( CV.TBS.tr.DoL = pbias.TBS.tr, CV.PTBS.tr.DoL = pbias.PTBS.tr, CV.PTBS.seplam.DoL = pbias.PTBS.seplam.tr )
} ) )


s2.hat.DoL <- mle.DoL$sigma2.hat * n.obs/(n.obs-p-1)
covML.DoL <- solve( -optimHess( unlist(mle.DoL[c(2,3,1)]), PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = FALSE, control = list( maxit = 5000, fnscale = -1 ) ) )
var.pred.PTBS.transf.DoL <- sapply( BC.transform( mle.DoL$lambda.hat[1], x.DoL ), function(x){ c(1,x) %*% covML.DoL[1:p,1:p] %*% c(1,x) } ) + s2.hat.DoL

s2.hat.seplam.DoL <- mle.seplam.DoL$sigma2.hat * n.obs/(n.obs-p-2)
covML.seplam.DoL <- solve( -optimHess( unlist(mle.seplam.DoL[c(2,3,1)]), PTBS.llkhd, y = y.DoL, x = x.DoL, sep.lam = TRUE, control = list( maxit = 5000, fnscale = -1 ) ) )
var.pred.PTBS.seplam.transf.DoL <- sapply( BC.transform( mle.seplam.DoL$lambda.hat[2], x.DoL ), function(x){ c(1,x) %*% covML.seplam.DoL[1:p,1:p] %*% c(1,x) } ) + s2.hat.seplam.DoL

s2.hat.other.DoL <- mle.other.DoL$sigma2.hat * n.obs/(n.obs-p-1)
covML.other.DoL <- solve( -optimHess( unlist( mle.other.DoL[c(2,3,1)] ), loglik.other, y = y.DoL, x = x.DoL, control = list( fnscale = -1 ) ) )
var.pred.TBS.transf.DoL <- sapply( x.DoL, function(x){ c( c( 1, x ) %*% mle.other.DoL$beta.hat )^(2*(mle.other.DoL$lambda.hat - 1)) * c( c( 1, x ) %*% covML.other.DoL[1:p,1:p] %*% c(1, x) ) } ) + s2.hat.other.DoL
# With regression version of var(beta.hat) (gives similar result):
#sapply( x.DoL, function(x){ c( c( 1, x ) %*% mle.other.DoL$beta.hat )^(2*(mle.other.DoL$lambda.hat - 1)) * c( c( 1, x ) %*% ( s2.hat.other.DoL * solve( crossprod( diag( c( cbind( 1, x.DoL ) %*% mle.other.DoL$beta.hat )^(mle.other.DoL$lambda.hat - 1) ) %*% cbind( 1, x.DoL ) ) ) ) %*% c(1, x) ) } ) + s2.hat.other.DoL

var.pred.transf.DoL <- cbind( TBS = var.pred.TBS.transf.DoL, PTBS = var.pred.PTBS.transf.DoL, PTBS.seplam = var.pred.PTBS.seplam.transf.DoL )

par(mfrow=c(3,1))
for(i in 1:3) {
	boxplot( t( pbias.transf.resamp.DoL[ranking,i,] ), outlty = 1, outpch = NA, xlab = "Rank of target loss", ylab = "Prediction bias [-]", main = dimnames(pbias.transf.resamp.DoL)[[2]][i] )
	points( 1:n.obs, pbias.transf.fullCV.DoL[ranking,i], col = 'red', pch = 3, cex = 0.8 )
	lines( -qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.DoL[ranking,i]), col = 'blue', lty = 2 )
	lines( qt( 0.975, df = n.obs - p - 1 ) * sqrt(var.pred.transf.DoL[ranking,i]), col = 'blue', lty = 2 )
#	abline( h = seq(-5,6,1), col = 'grey', lty = 3 )
	abline( h = 0, col = 'blue' )
	rug( which( ( apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL[,i]) | apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL[,i]) )[ranking,i] ), col = 'red' )
	box()
}

apply( apply( pbias.transf.resamp.DoL, 1:2, function(e){ min(e) > 0 | max(e) < 0 } ), 2, sum ) # only very few less!!!
sum( apply( apply( pbias.transf.resamp.DoL, 1:2, function(e){ min(e) <= 0 & max(e) >= 0 } ), 1, all ) )  # predictions for 31 observation are unbiased under all three models

apply( apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL[,i]) | apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL), 2, sum )
apply( apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.025 ) > qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL[,i]) | apply( pbias.transf.resamp.DoL, 1:2, quantile, probs = 0.975 ) < -qt( 0.975, df = n.obs - p ) * sqrt(var.pred.transf.DoL), 2, which )
# c( 43, 69, 125, 218, 233, 248 )
#======


######
# Summary plots
######

library(abind)
se.bias <- cbind( mod = apply( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss, 2, sd ), resamp = apply( apply( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss, 2:3, sd ), 1, mean ), CV = apply( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss, 2, sd ) )
#apply( se.bias, 1, diff ) # all positive up to one case

se.rel.bias <- cbind( mod = apply( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss, 2, sd ), resamp = apply( apply( 100*( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss )/y.Loss, 2:3, sd ), 1, mean ), CV = apply( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss, 2, sd ) )
#apply( se.rel.bias, 1, diff ) # all positive up to one case (not the same as before)

se.mae <- cbind( mod = apply( abs( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss ), 2, sd ), resamp = apply( apply( abs( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss ), 2:3, sd ), 1, mean ), CV = apply( abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss ), 2, sd ) )
#apply( se.mae, 1, diff ) # all positive except DoL.TBS

se.rel.ae <- cbind( mod = apply( abs( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss ), 2, sd ), resamp = apply( apply( abs( 100*( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss )/y.Loss ), 2:3, sd ), 1, mean ), CV = apply( abs( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss ), 2, sd ) )
#apply( se.rel.ae, 1, diff ) # all positive up to one case
#--> huge variability across the sample

se.preds <- abind( se.rmspe = matrix( NA, 12, 3 ), se.bias = se.bias, se.rel.tot.bias = matrix( NA, 12, 3 ), se.rel.bias = se.rel.bias, se.mae = se.mae, se.rel.tot.ae = matrix( NA, 12, 3 ), se.rel.ae = se.rel.ae, along = 3 )


err.name <- c( "Root mean square prediction error [CHF]", "Bias [CHF]", "Relative total bias [%]", "Relative bias [%]", "Mean absolute error [CHF]", "Relative total absolute error [%]", "Relative absolute error [%]" )
# Summary plots for prediction on original scale:
par( mfrow = c(2,3), mar = c( 9, 4, 3.5, 4 ) )
for( i in c(2,4,1,5,7) ) {
	plot( rep( 1:12, 4 ) + rep( c(-0.3,-0.1,0.1,0.3), each = 12 ), c( Delta.CV.Loss[1:3,i], Delta.CV.DoL[1:3,i], Delta.CV.Loss[4:6,i], Delta.CV.DoL[4:6,i], Delta.B.Loss[1:3,i], Delta.B.DoL[1:3,i], Delta.B.Loss[4:6,i], Delta.B.DoL[4:6,i], Delta.632.Loss[1:3,i], Delta.632.DoL[1:3,i], Delta.632.Loss[4:6,i], Delta.632.DoL[4:6,i], Delta.app.Loss[1:3,i], Delta.app.DoL[1:3,i], Delta.app.Loss[4:6,i], Delta.app.DoL[4:6,i] ), xlab = "", ylab = err.name[i], main = "Aggregate prediction error and variability", pch = rep( c(4,1,0,6), each = 12 ), xaxt = 'n', mgp = c(2.5,1,0) )
	axis( 1, at = 1:12, labels = c( rownames(Delta.CV.Loss)[1:3], rownames(Delta.CV.DoL)[1:3], rownames(Delta.CV.Loss)[4:6], rownames(Delta.CV.DoL)[4:6] ), las = 2 )
	#axis( 1, at = 1:12, labels = c( sapply( strsplit( rownames(Delta.CV.Loss)[1:6], ".me" ), '[', 1 ), sapply( strsplit( rownames(Delta.CV.DoL)[1:6], ".me" ), '[', 1 ) ), las = 2 )
	abline( h = 0, col = 'red' )
	abline( v = seq(0.5,12.5,1), lty = 3 )
	if(i!=1) {
		par(new = TRUE)
		plot( rep( 1:12, each = 4 ) + rep( c(-0.3,0,0.3,NA), 12 ), c( t(cbind(se.preds[,,i],NA)) ), axes = FALSE, xlab = "", ylab = "", type = 'l', lty = 6, col = 'blue' )
		axis( 4, at = pretty(range(se.preds[,,i])), col = 'blue', col.axis = 'blue' )
		mtext( paste( "Standard deviation in the sample [", ifelse( i==2 | i==5, "CHF", "%"), "]", sep = '' ), side = 4, line = 2.5, col= 'blue', cex = par("cex") )
	}
}
# X: Cross-validation, O: Booststrap, []: 632, V: Apparent error

# General patterns in aggregate errors:
# Bias smaller for mean than for median; se smaller for DoL than for Loss
# Absolute error smaller for DoL than for Loss and for median than for mean; se smaller for DoL than for Loss
# Relative bias and relative absolute error smaller for median than for mean and smaller for Loss than for DoL; se same pattern
# RMSPE smaller for DoL than for Loss; no se available




# Ranks of the models in terms of different prediction errors
#Bias: 11, 12, 10, 7, 9, 8, 1, 5, 3, 6, 2, 4
#Relative bias: 1, 3, 1, 4, 6, 5, 7, 9, 7, 10, 12, 11
#MAE: 7, 9, 9 , 2, 1, 3, 8, 11, 12, 5, 3, 6
#Relative AE: 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11
#RMSE: 8, 11, 10, 3, 2, 3, 7, 9, 12, 6, 1, 5

err.ranks <- rbind( bias = c( 11, 12, 10, 7, 9, 8, 1, 5, 3, 6, 2, 4 ), rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), mae = c( 7, 9.5, 9.5 , 2, 1, 3.5, 8, 11, 12, 5, 3.5, 6 ), rel.ae = c( 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11 ), rmse = c( 8, 11, 10, 3.5, 2, 3.5, 7, 9, 12, 6, 1, 5 ) )

apply( err.ranks, 2, sum ) # 4
apply( err.ranks == apply( err.ranks, 1, min ), 2, sum ) # 1
apply( err.ranks == apply( err.ranks, 1, max ), 2, sum ) # not 9, 11, 2
apply( err.ranks, 2, max ) # 4
apply( err.ranks, 2, min ) # 1, 3, 5, 7, 11

# Ranks of the models in terms of the std deviation for different prediction errors
#Bias: 7.5, 11, 10, 4, 2, 3, 7.5, 9, 12, 6, 1, 5
#Relative bias: 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11
#MAE: 9, 11.5, 10, 3.5, 2, 3.5, 7, 8, 11.5, 6, 1, 5
#Relative AE: 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11

se.ranks <- rbind( se.bias = c( 7.5, 11, 10, 4, 2, 3, 7.5, 9, 12, 6, 1, 5 ), se.rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), se.mae = c( 9, 11.5, 10, 3.5, 2, 3.5, 7, 8, 11.5, 6, 1, 5 ), se.rel.ae = c( 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ) )

apply( se.ranks, 2, sum ) # 4(, 5)
apply( se.ranks == apply( se.ranks, 1, min ), 2, sum ) # 1, 3, 11
apply( se.ranks == apply( se.ranks, 1, max ), 2, sum ) # not 9, 11, 2
apply( se.ranks, 2, max ) # 4
apply( se.ranks, 2, min ) # 1, 3, 11, not 7:9

floor( err.ranks[1:4,] ) + floor( se.ranks )

comb.ranks <- rbind( bias = c( 10, 12, 11, 5, 5, 5, 2, 8, 9, 7, 1, 3 ), rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), mae = c( 8, 11, 9.5, 3, 1, 4, 7, 9.5, 12, 5.5, 2, 5.5 ), rel.ae = c( 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11 ) )
apply( comb.ranks, 2, sum ) # 4 followed by 5, 6, 1

apply( comb.ranks, 2, sum ) # 4(, 5)
apply( comb.ranks == apply( comb.ranks, 1, min ), 2, sum ) # 1, 3, 5, 11
apply( comb.ranks == apply( comb.ranks, 1, max ), 2, sum ) # not 11, 9, 2
apply( comb.ranks, 2, max ) # 4, 6
apply( comb.ranks, 2, min ) # 1, 3, 5, 11, not 7:8


# TBS Loss mean  (best in terms of bias for Loss, has smallest mean bias but larger se than other models) [#7 overall]
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.Loss[,4]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.Loss[,4]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.Loss[,4]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,4]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "TBS model for absolute loss, mean prediction", outer = TRUE, line = 0.8 )

# PTBS DoL mean (overall best in terms of bias when accounting also for se, best in terms of RMSPE, co-best in terms of absolute error) [#11 overall]
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.DoL[,5]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.DoL[,5]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.DoL[,5]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,5]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "PTBS model for relative loss, mean prediction", outer = TRUE, line = 0.8 )

# PTBS DoL median (overall best in terms of absolute error, practically together with PTBS DoL mean) [#5 overall]
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.DoL[,2]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.DoL[,2]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.DoL[,2]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,2]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "PTBS model for relative loss, median prediction", outer = TRUE, line = 0.8 )

# Best in terms of relative bias and relative aboslute error: TBS median Loss and PTBS.seplam median Loss
# TBS Loss median [#1 overall]
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.Loss[,1]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.Loss[,1]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.Loss[,1]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,1]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "TBS model for absolute loss, median prediction", outer = TRUE, line = 0.8 )
# PTBS.seplam Loss median [#3 overall]
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.Loss[,3]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.Loss[,3]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.Loss[,3]-y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.Loss[,3]-y.Loss)/y.Loss)[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "PTBS model with separate lambda for absolute loss, median prediction", outer = TRUE, line = 0.8 )

# TBS median DoL [#4 overall] (best in terms of overall ranks)
par( mfrow = c(2,2), mar = c(4,4,1,1), oma = c(0,0,2,0) )
plot( 1:n.obs, (pred.fullCV.DoL[,1]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Bias [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, 100*((pred.fullCV.DoL[,1]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative bias [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
plot( 1:n.obs, abs(pred.fullCV.DoL[,1]-y.DoL*y.InSum)[ranking], xlab = "Rank of target loss", ylab = "Mean absolute error [CHF]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
plot( 1:n.obs, abs(100*(pred.fullCV.DoL[,1]-y.DoL*y.InSum)/(y.DoL*y.InSum))[ranking], xlab = "Rank of target loss", ylab = "Relative absolute error [%]", mgp = c(2.5,0.8,0), pch = 20, cex = 0.8 )
abline( h = 0, col = 'blue' )
abline( h = 100, col = 'blue', lty = 2 )
title( main = "TBS model for relative loss, median prediction", outer = TRUE, line = 0.8 )
##


##=== Take median of errors as aggregate error --> not good idea
library(abind)
mad.med.bias <- cbind( mod = apply( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss, 2, mad ), resamp = apply( apply( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss, 2:3, mad ), 1, mean ), CV = apply( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss, 2, mad ) )
#apply( mad.med.bias, 1, diff ) # no pattern in the signs

mad.med.rel.bias <- cbind( mod = apply( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss, 2, mad ), resamp = apply( apply( 100*( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss )/y.Loss, 2:3, mad ), 1, mean ), CV = apply( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss, 2, mad ) )
#apply( mad.med.rel.bias, 1, diff ) # no pattern in the signs

mad.med.ae <- cbind( mod = apply( abs( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss ), 2, mad ), resamp = apply( apply( abs( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss ), 2:3, mad ), 1, mean ), CV = apply( abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss ), 2, mad ) )
#apply( mad.med.ae, 1, diff ) # no pattern in the signs

mad.med.rel.ae <- cbind( mod = apply( abs( 100*( cbind( pred.mod.Loss[,1:3], pred.mod.DoL[,1:3], pred.mod.Loss[,4:6], pred.mod.DoL[,4:6] ) - y.Loss )/y.Loss ), 2, mad ), resamp = apply( apply( abs( 100*( abind( pred.resamp.Loss[,1:3,], pred.resamp.DoL[,1:3,], pred.resamp.Loss[,4:6,], pred.resamp.DoL[,4:6,], along = 2 ) - y.Loss )/y.Loss ), 2:3, mad ), 1, mean ), CV = apply( abs( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss ), 2, mad ) )
#apply( mad.med.rel.ae, 1, diff ) # no pattern in the signs
#--> huge variability across the sample

mad.preds <- abind( mad.med.bias = mad.med.bias, mad.med.rel.bias = mad.med.rel.bias, mad.med.ae = mad.med.ae, mad.med.rel.ae = mad.med.rel.ae, along = 3 )

err.med.name <- c( "Median bias [CHF]", "Median relative bias [%]", "Median absolute error [CHF]", "Median relative absolute error [%]" )
# Summary plots for prediction on original scale:
par( mfrow = c(2,2), mar = c( 9, 4, 3.5, 4 ) )
for( i in 1:4 ) {
	plot( rep( 1:12, 4 ) + rep( c(-0.3,-0.1,0.1,0.3), each = 12 ), c( Delta.CV.med.Loss[1:3,i], Delta.CV.med.DoL[1:3,i], Delta.CV.med.Loss[4:6,i], Delta.CV.med.DoL[4:6,i], Delta.B.med.Loss[1:3,i], Delta.B.med.DoL[1:3,i], Delta.B.med.Loss[4:6,i], Delta.B.med.DoL[4:6,i], Delta.632.med.Loss[1:3,i], Delta.632.med.DoL[1:3,i], Delta.632.med.Loss[4:6,i], Delta.632.med.DoL[4:6,i], Delta.app.med.Loss[1:3,i], Delta.app.med.DoL[1:3,i], Delta.app.med.Loss[4:6,i], Delta.app.med.DoL[4:6,i] ), xlab = "", ylab = err.med.name[i], main = "Median prediction error and variability", pch = rep( c(4,1,0,6), each = 12 ), xaxt = 'n', mgp = c(2.5,1,0) )
	axis( 1, at = 1:12, labels = c( rownames(Delta.CV.Loss)[1:3], rownames(Delta.CV.DoL)[1:3], rownames(Delta.CV.Loss)[4:6], rownames(Delta.CV.DoL)[4:6] ), las = 2 )
	#axis( 1, at = 1:12, labels = c( sapply( strsplit( rownames(Delta.CV.Loss)[1:6], ".me" ), '[', 1 ), sapply( strsplit( rownames(Delta.CV.DoL)[1:6], ".me" ), '[', 1 ) ), las = 2 )
	abline( h = 0, col = 'red' )
	abline( v = seq(1.5,12.5,1), lty = 3 )
	par(new = TRUE)
	plot( rep( 1:12, each = 4 ) + rep( c(-0.3,0,0.3,NA), 12 ), c( t(cbind(mad.preds[,,i],NA)) ), axes = FALSE, xlab = "", ylab = "", type = 'l', lty = 6, col = 'blue' )
	axis( 4, at = pretty(range(mad.preds[,,i])), col = 'blue', col.axis = 'blue' )
	mtext( paste( "Standard deviation in the sample [", ifelse( i==2 | i==5, "CHF", "%"), "]", sep = '' ), side = 4, line = 2.5, col= 'blue', cex = par("cex") )
}
# X: Cross-validation, O: Booststrap, []: 632, V: Apparent error

# General patterns:
# Bias: median prediction for Loss performs best.
# Relative bias: Median prediction much better than mean prediction.
c( 4, 2, 3, 5.5, 5.5, 1 ) + c( 3.5, 6, 5, 2, 1, 3.5 )
# Median AE: Median prediction for Loss best, plus maybe TBS.mean.Loss.
# Median relative AE: Median prediction for DoL performs best.

# Bias: 1 or 3
# AE: 1, followed by 3
# Relative bias: 6, followed by 5
# Relative AE: 5, followed by 4 and 6

# After all this seems not a very wise criterion because it "allows" huge errors for single predictions while keeping a small to moderate median error.


# Boxplots
boxplot( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss, ylab = "Bias [CHF]", las = 3 )
abline( h = 0, col = 'blue' )
lines( apply( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss, 2, mean ), col = 'red', lwd = 1.8, lty = 5 )
boxplot( 100*( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss, ylab = "Relative bias [%]" )
abline( h = 0, col = 'blue' )
boxplot( abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss ), ylab = "Mean absolute error [CHF]" )
abline( h = 0, col = 'blue' )
boxplot( 100*abs( cbind( pred.fullCV.Loss[,1:3], pred.fullCV.DoL[,1:3], pred.fullCV.Loss[,4:6], pred.fullCV.DoL[,4:6] ) - y.Loss )/y.Loss, ylab = "Relative absolute error [%]", las = 3 )
abline( h = 100, lty = 2 )
lines( 1:12, c( Delta.CV.Loss[1:3,i], Delta.CV.DoL[1:3,i], Delta.CV.Loss[4:6,i], Delta.CV.DoL[4:6,i] ), col = 'red', lty = 5, lwd = 1.8 )

