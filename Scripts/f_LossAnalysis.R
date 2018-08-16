#setwd("../MOBILAB/Statistik-Hilfe GIUB/Markus_VulnerabilitätFahrhabe")
#load("Loss/Loss_Data.RData")
# contains 384 by 5 dataframes "Contents" and "Structure"

#all( Contents$L_ID == 1:n.obs )
#all( Structure$L_ID == 1:n.obs )
# L_ID is just enumeration

#all( Contents$Loss/Contents$InSum == Contents$DoL )
#all( Structure$Loss/Structure$InSum == Structure$DoL )
# DoL = Loss/InSum


## Functions

BC.transform <- function( lambda, x ) {
	if( abs(lambda) <= 1.e-6 ) log(x) else ((x^lambda) - 1)/lambda
}

BC.backtransform <- function( lambda, x ) {
	if( abs(lambda) <= 1.e-6 ) exp(x) else (1 + lambda*x)^(1/lambda)
}

# -(n/2) [ log{ [ y^{(lambda)}]^T [I - H(lambda)] y^{(lambda)} } - (lambda - 1) log{ ( Prod_{i=1:n} y_i )^{2/n} } ]

PTBS.profllkhd.lambda <- function( lambda, y, x, intercept = TRUE, dummy.zeros = FALSE ) {
	n <- length(y)
	
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	
	if( length(lambda) == 1 ) 
		lambda.y <- lambda.x <- lambda
	else {
		lambda.y <- lambda[1]
		lambda.x <- lambda[2]
	}
	
	y.lambda <- BC.transform( lambda.y, y )
	x.to.transf <- if(is.numeric(intercept)) x[,-intercept,drop=F] else x
	x.lambda <- if(dummy.zeros) apply( x.to.transf, 1:2, function(x.ij){ ifelse( x.ij == 0, x.ij, BC.transform( lambda.x, x.ij ) ) } ) else BC.transform( lambda.x, x.to.transf )
	X.lambda <- if( is.logical(intercept) & intercept[1] ) cbind( 1, x.lambda ) else if( is.logical(intercept) & !intercept[1] ) x.lambda else cbind( x[,intercept,drop=F], x.lambda )

	beta.lambda <- solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda
	
	res.lambda <- c( y.lambda - X.lambda %*% beta.lambda )
	sigma2.lambda <- sum(res.lambda^2)/n
	
	#-1/2 * ( n*log(sigma2.lambda) + n*sigma2.lambda/sigma2.lambda ) + (lambda - 1)*sum(log(y)) 
	-(n/2)*log(sigma2.lambda) - n/2 + (lambda.y - 1)*sum(log(y))
}
PTBS.profllkhd.lamsepopt <- function( lambda.opt, lambda.fix, idx.opt = 1, y, x, intercept, dummy.zeros ) {
	if( idx.opt == 1 )
		lambda <- c(lambda.opt,lambda.fix)
	else
		lambda <- c(lambda.fix,lambda.opt)
	PTBS.profllkhd.lambda( lambda, y, x, intercept = intercept, dummy.zeros = dummy.zeros )
}

prof.shifted <- function( lambda.prof, y, x, f0, idx.lamopt = NULL, intercept, dummy.zeros, intval ) {
	if(is.null(idx.lamopt)) {
		PTBS.profllkhd.lambda( lambda.prof, y, x, intercept, dummy.zeros ) - f0
	} else {
		lambda.opt <- optimize( f = PTBS.profllkhd.lamsepopt, interval = intval, lambda.fix = lambda.prof, idx.opt = idx.lamopt, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros, maximum = TRUE )$maximum
		
		lambda <- if(idx.lamopt==1) { c(lambda.opt,lambda.prof) } else { c(lambda.prof,lambda.opt) }
		PTBS.profllkhd.lambda( lambda, y, x, intercept, dummy.zeros ) - f0
	}
}

PTBS.CI.lambda <- function( y, x, sep.lam = FALSE, intercept = TRUE, dummy.zeros = FALSE, level = 0.95, interval = c(-3,3), lambda.init = c(0.5,0.5) ) {
	if(!sep.lam) {
		prof.opt <- optimize( f = PTBS.profllkhd.lambda, interval = interval, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros, maximum = TRUE )
		lower <- uniroot( f = prof.shifted, interval = c( max(-2,interval[1]), prof.opt$maximum ), y = y, x = x, f0 = prof.opt$objective - qchisq( level, df = 1 )/2, intercept = intercept, dummy.zeros = dummy.zeros )$root
		upper <- uniroot( f = prof.shifted, interval = c( prof.opt$maximum, min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$objective - qchisq( level, df = 1 )/2, intercept = intercept, dummy.zeros = dummy.zeros )$root
		valname <- 'objective'
		maxname <- 'maximum'
	} else {
		prof.opt <- optim( lambda.init, fn = PTBS.profllkhd.lambda, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros, control = list( fnscale = -1 ) )
		print(prof.opt$convergence)
		lower <- c( uniroot( f = prof.shifted, interval = c( max(-2,interval[1]), prof.opt$par[1] ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 2, intval = interval, intercept = intercept, dummy.zeros = dummy.zeros )$root, uniroot( f = prof.shifted, interval = c( max(-2,interval[1]), prof.opt$par[2] ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 1, intval = interval, intercept = intercept, dummy.zeros = dummy.zeros )$root )
		upper <- c( uniroot( f = prof.shifted, interval = c( prof.opt$par[1], min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 2, intval = interval, intercept = intercept, dummy.zeros = dummy.zeros )$root, uniroot( f = prof.shifted, interval = c( prof.opt$par[2], min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 1, intval = interval, intercept = intercept, dummy.zeros = dummy.zeros )$root )
		valname <- 'value'
		maxname <- 'par'
	}
	list( lambda = prof.opt[[maxname]], lower = lower, upper = upper, loglik.opt = prof.opt[[valname]] )
}

PTBS.mle <- function( y, x, lambda = NULL, sep.lam = FALSE, interval = c(-3,3), intercept = TRUE, lambda.init = c(0.5,0.5), dummy.zeros = FALSE ) {  # dummy.zeros only works correctly if the covariate(s) with dummy zeros never take the value 0 themselves!
	n <- length(y)
	
	valname <- 'objective'
	
	if( is.null(lambda) & !sep.lam ) {
		ell.hat <- optimize( f = PTBS.profllkhd.lambda, interval = interval, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros, maximum = TRUE )
		lambda <- ell.hat$maximum
	} else if( is.null(lambda) & sep.lam ) {
		ell.hat <- optim( lambda.init, fn = PTBS.profllkhd.lambda, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros, control = list( fnscale = -1 ) )
		if( ell.hat$conv > 0 ) print(ell.hat$convergence)
		valname <- 'value'
		lambda <- ell.hat$par
	} else {
		ell.hat <- list( objective = PTBS.profllkhd.lambda( lambda, y = y, x = x, intercept = intercept, dummy.zeros = dummy.zeros ) )
	}

	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( length(lambda) == 1 ) 
		lambda.y <- lambda.x <- lambda
	else {
		lambda.y <- lambda[1]
		lambda.x <- lambda[2]
	}

	y.lambda <- BC.transform( lambda.y, y )
	x.to.transf <- if(is.numeric(intercept)) x[,-intercept,drop=F] else x
	x.lambda <- if(dummy.zeros) apply( x.to.transf, 1:2, function(x.ij){ ifelse( x.ij == 0, x.ij, BC.transform( lambda.x, x.ij ) ) } ) else BC.transform( lambda.x, x.to.transf )
	X.lambda <- if( is.logical(intercept) & intercept[1] ) cbind( 1, x.lambda ) else if( is.logical(intercept) & !intercept[1] ) x.lambda else cbind( x[,intercept,drop=F], x.lambda )			

	beta.lambda <- c( solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda )
	
	res.lambda <- c( y.lambda - X.lambda %*% beta.lambda )
	sigma2.lambda <- sum(res.lambda^2)/n  # MLE (biased)
		
	list( lambda.hat = lambda, beta.hat = beta.lambda, sigma2.hat = sigma2.lambda, ell.opt = ell.hat[[valname]] )
}

PTBS.llkhd <- function( theta, y, x, intercept = TRUE, sep.lam = FALSE, dummy.zeros = FALSE ) {
	n <- length(y)
	p <- length(theta)-2-sep.lam
	
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	
	y.lambda <- BC.transform( theta[p+2], y )
	x.to.transf <- if(is.numeric(intercept)) x[,-intercept,drop=F] else x
	x.lambda <- if(dummy.zeros) apply( x.to.transf, 1:2, function(x.ij){ ifelse( x.ij == 0, x.ij, BC.transform( theta[p+2+sep.lam], x.ij ) ) } ) else BC.transform( theta[p+2+sep.lam], x.to.transf )
	X.lambda <- if( is.logical(intercept) & intercept[1] ) cbind( 1, x.lambda ) else if( is.logical(intercept) & !intercept[1] ) x.lambda else cbind( x[,intercept,drop=F], x.lambda )
	
	# beta.lambda <- solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda
	
	res.lambda <- c( y.lambda - X.lambda %*% theta[1:p] )
	#sigma2.lambda <- sum(res.lambda^2)/n
	
	#-1/2 * ( n*log(sigma2.lambda) + n*sigma2.lambda/sigma2.lambda ) + (lambda - 1)*sum(log(y))
	# -1/2 * ( n*log(sigma2.lambda) + sum(res.lambda^2)/sigma2.lambda ) + (lambda - 1)*sum(log(y))
	#-(n/2)*log(theta[3]) + (theta[4] - 1)*sum(log(y))
	-( n*log(theta[p+1]) + sum(res.lambda^2)/theta[p+1] )/2 + (theta[p+2] - 1)*sum(log(y))
}

loglik.other <- function( theta, y, x, intercept = TRUE ) {
	#n <- length(y)
	p <- length(theta)-2

	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	
	X <- if(intercept) cbind( 1, x ) else x
		
	sum( dnorm( ( y^theta[p+2] - c( X %*% theta[1:p] )^theta[p+2] )/theta[p+2], sd = sqrt(theta[p+1]), log = TRUE ) ) + ( theta[p+2] - 1 )*sum(log(y))
}


# Gradient of x.lambda.hat %*% beta.hat as a fct of theta
grad.fitted <- function( theta, x0 ) { # x0 = c(1,x.star) on transformed scale, thus x0 with intercept
	p <- length(theta) - 2
	x.star <- x0[-1]
	
	d.beta <- x0
	d.sigma2 <- 0
	d.lambda <- sum( ( log( 1 + theta[p+2] * x.star )/theta[p+2] * ( x.star + 1/theta[p+2] ) - x.star/theta[p+2] )*theta[2:p] )
	
	c(d.beta, d.sigma2, d.lambda)
}

# Gradient of BC.backtransform( lambda.hat, x.lambda.hat ) as a fct of theta
grad.yhat <- function( theta, x.lam, x.beta = NULL, sep.lam = FALSE, pr = 0.5 ) { # x.lam on transformed scale, 
# = c(1,x.star) thus x.star with intercept
	p <- length(theta) - 2 - sep.lam
	sigma <- sqrt(theta[p+1])
	x.star.orig <- BC.backtransform( theta[p+2+sep.lam], x.lam[-1] )
	if(is.null(x.beta))
		x.beta <- c( x.lam %*% theta[1:p] )
		
	d.beta <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * x.lam  # x.lam = x0(lambda)/lambda
	
	d.sigma2 <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * qnorm(pr)/(2*sigma)
	
	#d.inner <- c( ( x.beta + sigma*qnorm(pr) )/theta[p+2], sum( ( log(x.star.orig) * x.star.orig^theta[p+2+sep.lam] - x.lam[-1] ) * theta[2:p] )/theta[p+2+sep.lam] )
	if(!sep.lam) { # single lambda: A*(B1 + B2) - C
		d.inner <- ( sum( c( 1, x.star.orig^theta[p+2] * log(x.star.orig) ) * theta[1:p] ) + sigma*qnorm(pr) )/theta[p+2]
		#d.inner <- sum(d.inner)
		mult.term2 <- 1
		#d.inner <- ( x.beta + sigma*qnorm(pr) + sum( ( log(x.star.orig) * x.star.orig^theta[p+2] - x.lam[-1] ) * theta[2:p] ) )/theta[p+2]
		# d.lambda <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * ( sum( c( 1, x.star.orig^theta[p+2] * log(x.star.orig) ) * theta[1:p] ) + sigma*qnorm(pr) )/theta[p+2] - ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2]) * log( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )/theta[p+2]^2
	} else { # separate lambda: c( A*B1 - C, A*B2 ) = A*c( B1, B2 ) - c( C, 0 )
		d.inner <- c( ( x.beta + sigma*qnorm(pr) )/theta[p+2], sum( ( log(x.star.orig) * x.star.orig^theta[p+2+sep.lam] - x.lam[-1] ) * theta[2:p] )/theta[p+2+sep.lam] )
		mult.term2 <- c(1,0)
		# d.lambda1 <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * ( x.beta + sigma*qnorm(pr) )/theta[p+2] - ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2]) * log( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )/theta[p+2]^2
		# d.lambda2 <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * sum( ( log(x.star.orig) * x.star.orig^theta[p+2+sep.lam] - x.lam[-1] ) * theta[2:p] )/theta[p+2+sep.lam]
		# d.lambda <- c(d.lambda1, d.lambda2)
	}
	d.lambda <- ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2] - 1) * d.inner - mult.term2 * ( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )^(1/theta[p+2]) * log( 1 + theta[p+2] * x.beta + theta[p+2] * sigma*qnorm(pr) )/theta[p+2]^2	
	c(d.beta, d.sigma2, d.lambda)
}

condmean.orig <- function( theta, x, x.beta = NULL, sep.lam = FALSE ) { # theta = ( beta, sigma2, lambda ), x on transformed scale and with intercept
# lam2 is implicitly contained in x.
	p <- length(theta) - 2 - sep.lam
	if(is.null(x.beta))
		x.beta <- c( x %*% theta[1:p] )
		
	BC.backtransform( theta[p+2], x.beta ) * ( 1 + ( theta[p+1]*(1-theta[p+2]) )/(2*( 1 + theta[p+2] * x.beta )^2) )
}

grad.condmean.orig <- function( theta, x, sep.lam = FALSE ) { # x vector on transformed scale, with intercept
	p <- length(theta) - 2 - sep.lam
	x.beta <- sum( x*theta[1:p] )
	
	x.star.orig <- BC.backtransform( theta[p+2+sep.lam], x[-1] )  # untransformed scale, w/o intercept

	d.beta <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 1) * ( 1 + ( theta[p+1]*(1-theta[p+2])*(1 - 2*theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) * x  # OK
	
	d.sigma2 <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 2) * (1 - theta[p+2])/2 # OK
	
	if(!sep.lam) {
		d.inner <- sum( c( 1, x.star.orig^theta[p+2] * log(x.star.orig) ) * theta[1:p] )/theta[p+2]
		mult.term2 <- 1
		# d.lambda <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 1) * ( 1 + ( theta[p+1]*(1-theta[p+2])*(1 - 2*theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) * sum( c( 1, x.star.orig^theta[p+2] * log(x.star.orig) ) * theta[1:p] )/theta[p+2] - ( 1 + theta[p+2] * x.beta )^(1/theta[p+2]) * ( ( 1 + ( theta[p+1]*(1-theta[p+2]) )/(2*( 1 + theta[p+2] * x.beta )^2) ) * log( 1 + theta[p+2] * x.beta )/theta[p+2]^2 + theta[p+1]/(2*( 1 + theta[p+2] * x.beta )^2) )
		### dependence of x on lambda included
		#( 1 + theta[p+2] * x.beta )^(1/theta[p+2]) * ( x.beta/(theta[p+2]*( 1 + theta[p+2] * x.beta )) - log( 1 + theta[p+2] * x.beta )/theta[p+2]^2 ) * ( 1 + ( theta[p+1]*(1-theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) - theta[p+1]/2 * ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 3) * ( 1 + 2*(1 - theta[p+2])*x.beta )	
	 } else {
		d.inner <- c( sum( x * theta[1:p] )/theta[p+2], sum( ( log(x.star.orig) * x.star.orig^theta[p+2+sep.lam] - x[-1] ) * theta[2:p] )/theta[p+2+sep.lam] )
		mult.term2 <- c(1,0)
	 	# d.lambda1 <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 1) * ( 1 + ( theta[p+1]*(1-theta[p+2])*(1 - 2*theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) * sum( x * theta[1:p] )/theta[p+2] - ( 1 + theta[p+2] * x.beta )^(1/theta[p+2]) * ( ( 1 + ( theta[p+1]*(1-theta[p+2]) )/(2*( 1 + theta[p+2] * x.beta )^2) ) * log( 1 + theta[p+2] * x.beta )/theta[p+2]^2 + theta[p+1]/(2*( 1 + theta[p+2] * x.beta )^2) )
		# d.lambda2 <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 1) * ( 1 + ( theta[p+1]*(1-theta[p+2])*(1 - 2*theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) * sum( ( log(x.star.orig) * x.star.orig^theta[p+2+sep.lam] - x[-1] ) * theta[2:p] )/theta[p+2+sep.lam]
		# d.lambda <- c( d.lambda1, d.lambda2 )
	 }
	 d.lambda <- ( 1 + theta[p+2] * x.beta )^(1/theta[p+2] - 1) * ( 1 + ( theta[p+1]*(1-theta[p+2])*(1 - 2*theta[p+2]) )/(2*(1 + theta[p+2] * x.beta )^2) ) * d.inner - mult.term2 * ( 1 + theta[p+2] * x.beta )^(1/theta[p+2]) * ( ( 1 + ( theta[p+1]*(1-theta[p+2]) )/(2*( 1 + theta[p+2] * x.beta )^2) ) * log( 1 + theta[p+2] * x.beta )/theta[p+2]^2 + theta[p+1]/(2*( 1 + theta[p+2] * x.beta )^2) )
	 
	c(d.beta, d.sigma2, d.lambda)
}

###
regr.llkhd <- function( theta, y, X, sigma = FALSE ) {
	n <- length(y)
	p <- length(theta)-1
		
	# beta.lambda <- solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda
	
	res <- c( y - X %*% theta[1:p] )
	#sigma2.lambda <- sum(res.lambda^2)/n
	
	#-1/2 * ( n*log(sigma2.lambda) + n*sigma2.lambda/sigma2.lambda ) + (lambda - 1)*sum(log(y)) 
	#-(n/2)*log(theta[3]) + (theta[4] - 1)*sum(log(y))
	if(sigma) -( 2*n*log(theta[p+1]) + sum(res^2)/theta[p+1]^2 )/2
	else -( n*log(theta[p+1]) + sum(res^2)/theta[p+1] )/2
}
###

## lambda transforming to a symmetric distribution (p.135 Carroll and Ruppert 1984a)
# More robust version of skewness:
T.skew <- function( lambda, x, y, intercept = TRUE, robust = TRUE ) {  # ! not sure that optimal quantiles !
# !! Only works if x contains no dummy covariates !!
	n <- length(y)
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	y.lambda <- BC.transform( lambda, y )
	x.lambda <- BC.transform( lambda, x )

	X.lambda <- if(intercept) cbind( 1, x.lambda ) else x.lambda

	H.lambda <- X.lambda %*% solve( crossprod(X.lambda) ) %*% t(X.lambda)
	
	res.lambda <- c( ( diag(n) - H.lambda ) %*% y.lambda )

	if(robust) {
		quarts <- quantile( res.lambda, probs = c(0.25,0.5,0.75), type = 8 )
		(quarts[3] - quarts[2])/(quarts[2] - quarts[1]) - 1
	} else {
		mean(res.lambda^3)/(sum(res.lambda^2)/(n-2))^(3/2)
	}
}

## lambda transforming to a homoskedastic distribution (p.135 Carroll and Ruppert 1984a)
# More robust version (Spearman's correlation instead of Pearson)
H.hesk <- function( lambda, x, y, intercept = TRUE, robust = TRUE ) {
# !! Only works if x contains no dummy covariates !!
	n <- length(y)
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	y.lambda <- BC.transform( lambda, y )
	x.lambda <- BC.transform( lambda, x )

	X.lambda <- if(intercept) cbind( 1, x.lambda ) else x.lambda
	
	H.lambda <- X.lambda %*% solve( crossprod(X.lambda) ) %*% t(X.lambda)
	
	res.lambda <- c( ( diag(n) - H.lambda ) %*% y.lambda )

	if(robust)
		cor( abs(res.lambda), c(H.lambda %*% y.lambda), method = "spearman" )
	else
		cor( res.lambda^2, log(c(H.lambda %*% y.lambda)), method = "pearson" )
}

PTBS.llkhd.knot <- function( theta.toOpt, theta = NULL, y, x, intercept = TRUE, var.split = FALSE ) {
	if( is.null(theta) ) {
		theta <- theta.toOpt
	} else {
		theta[is.na(theta)] <- theta.toOpt
	}
	PTBS.llkhd.knot.full( theta, y, x, intercept, var.split )
}

## Constraint that the two lines match at separation point knot.x:
# alpha0 + alpha1 * knot.x = beta0 + beta1 * knot.x
# alpha0 + (alpha1 - beta1) * knot.x = beta0
# For x < knot.x: alpha0 + alpha1 x
# For x >= knot.x: alpha0 + (alpha1 - beta1) * knot.x + beta1 x = alpha0 + alpha1 knot.x + beta1 (x - knot.x)
PTBS.llkhd.knot.full <- function( theta, y, x, intercept = TRUE, var.split = FALSE ) { # theta = c(alpha0, alpha1, beta1, sigma2, lambda, knot.x), knot.x on transformed scale, only if knot.fixed = FALSE
# !! Only works if x contains no dummy covariates !!
	if(!is.logical(var.split) )
		stop("!! var.split must be logical.")
	p <- length(theta) - 3 - var.split
	knot.x <- tail(theta,1)

	n <- length(y)
	
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	
	y.lambda <- BC.transform( theta[p+2+var.split], y )
	x.lambda <- BC.transform( theta[p+2+var.split], x )
	X.lambda <- if(intercept) cbind( 1, pmin( x.lambda, knot.x ), pmax( x.lambda - knot.x, 0 ) ) else cbind( pmin( x.lambda, knot.x ), pmax( x.lambda - knot.x, 0 ) )
		
	res.lambda <- c( y.lambda - X.lambda %*% theta[1:p] )
	
	if(!var.split)
		-( n*log(theta[p+1]) + sum(res.lambda^2)/theta[p+1] )/2 + (theta[p+2] - 1)*sum(log(y))
	else {
		below <- (x.lambda < knot.x)
		n1 <- sum(below)
		n2 <- n - n1
		-( n1*log(theta[p+1]) + n2*log(theta[p+2]) + sum(res.lambda[below]^2)/theta[p+1] + sum(res.lambda[!below]^2)/theta[p+2] )/2 + (theta[p+3] - 1)*sum(log(y))
	}
	# -1/2 * ( n1*log(sigma2.1.lambda) + n2*log(sigma2.2.lambda) + sum(res.lambda[idx.1]^2)/sigma2.1.lambda + sum(res.lambda[idx.2]^2)/sigma2.2.lambda ) + (lambda - 1)*sum(log(y))	
}

PTBS.profllkhd.knot.lambda <- function( lambda, y, x, knot.x, intercept = TRUE ) {
# !! Only works if x contains no dummy covariates !!
	n <- length(y)
		
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	
	y.lambda <- BC.transform( lambda, y )
	x.lambda <- BC.transform( lambda, x )
	X.lambda <- if( (knot.x > min(x.lambda)) & (knot.x < max(x.lambda)) ) cbind( ifelse( intercept, 1, NULL ), pmin( x.lambda, knot.x ), pmax( x.lambda - knot.x, 0 ) ) else cbind( ifelse( intercept, 1, NULL ), x.lambda )

	beta.lambda <- solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda
	
	res.lambda <- c( y.lambda - X.lambda %*% beta.lambda )
	sigma2.lambda <- sum(res.lambda^2)/n
	
	#-1/2 * ( n*log(sigma2.lambda) + n*sigma2.lambda/sigma2.lambda ) + (lambda - 1)*sum(log(y)) 
	-(n/2)*log(sigma2.lambda) - n/2 + (lambda - 1)*sum(log(y))
}

PTBS.knot.CI.lambda <- function( y, x, knot, intercept = TRUE, level = 0.95, interval = c(-3,3) ) {
	prof.opt <- optimize( f = PTBS.profllkhd.knot.lambda, interval = interval, y = y, x = x, knot.x = knot, maximum = TRUE )
	prof.shifted <- function( lambda, y, x, knot.x, f0, intercept = TRUE ) {
		PTBS.profllkhd.knot.lambda( lambda, y, x, knot.x, intercept ) - f0
	}
	lower <- uniroot( f = prof.shifted, interval = c( max(-2,interval[1]), prof.opt$maximum ), y = y, x = x, knot.x = knot, f0 = prof.opt$objective - qchisq( level, df = 1 )/2 )$root
	upper <- uniroot( f = prof.shifted, interval = c( prof.opt$maximum, min(2,interval[2]) ), y = y, x = x, knot.x = knot, f0 = prof.opt$objective - qchisq( level, df = 1 )/2 )$root
	list( lambda = prof.opt$maximum, lower = lower, upper = upper, loglik.opt = prof.opt$objective )
}

PTBS.profllkhd.knot.kn <- function( kn, y, x, lambda = NULL, intercept = TRUE, intval.lam = c(-3,3) ) {
	if(is.null(lambda))
		optimize( f = PTBS.profllkhd.knot.lambda, interval = intval.lam, y = y, x = x, knot.x = kn, intercept = intercept, maximum = TRUE )$objective
	else
		PTBS.profllkhd.knot.lambda( lambda = lambda, y = y, x = x, knot.x = kn, intercept = intercept )
}

PTBS.knot.mle <- function( y, x, lambda = NULL, knot, interval = c(-3,3), intercept = TRUE ) {
# !! Only works if x contains no dummy covariates !!
	n <- length(y)
	
	if( is.null(lambda) ) {
		ell.hat <- optimize( f = PTBS.profllkhd.knot.lambda, interval = interval, y = y, x = x, knot.x = knot, intercept = intercept, maximum = TRUE )
		lambda <- ell.hat$maximum
	} else {
		ell.hat <- list( objective = PTBS.profllkhd.knot.lambda( lambda, y = y, x = x, knot.x = knot, intercept = intercept ) )
	}

	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	y.lambda <- BC.transform( lambda, y )
	x.lambda <- BC.transform( lambda, x ) # !!! only works if x.lambda is a vector !!!
	X.lambda <- if( (knot > min(x.lambda)) & (knot < max(x.lambda)) ) cbind( ifelse( intercept, 1, NULL ), pmin( x.lambda, knot ), pmax( x.lambda - knot, 0 ) ) else cbind( ifelse( intercept, 1, NULL ), x.lambda )

	beta.lambda <- c( solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda )
	
	res.lambda <- c( y.lambda - X.lambda %*% beta.lambda )
	sigma2.lambda <- sum(res.lambda^2)/n  # MLE (biased)
	
	list( lambda.hat = lambda, beta.hat = beta.lambda, sigma2.hat = sigma2.lambda, ell.opt = ell.hat$objective )
}


TBS.RSS.beta <- function( beta.coef, lambda, y, X ) { # y and X on original scale
	if( length(lambda) == 1 ) 
		lambda.y <- lambda.x <- lambda
	else {
		lambda.y <- lambda[1]
		lambda.x <- lambda[2]
	}
	
	y.lambda <- BC.transform( lambda.y, y )
	Xbeta.lambda <- BC.transform( lambda.x, c( X %*% beta.coef ) )
	sum( (y.lambda - Xbeta.lambda)^2 )
}

TBS.profllkhd.lambda <- function( lambda, y, x, intercept = TRUE, beta.init = NULL ) {
	n <- length(y)
	
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	if( dim(x)[1] != n )
		stop( "x should have length(y) entries or rows.")
	
	if( length(lambda) == 1 ) 
		lambda.y <- lambda.x <- lambda
	else {
		lambda.y <- lambda[1]
		lambda.x <- lambda[2]
	}
	
	X <- if(intercept) cbind( 1, x ) else x
		
	if(is.null(beta.init))
		beta.init <- my.lm( y, X )$coef

	beta.lambda <- optim( beta.init, fn = TBS.RSS.beta, lambda = lambda, y = y, X = X )$par #solve( crossprod(X.lambda) ) %*% t(X.lambda) %*% y.lambda
	
	y.lambda <- BC.transform( lambda.y, y )
	Xbeta.lambda <- BC.transform( lambda.x, c( X %*% beta.lambda ) )

	res.lambda <- y.lambda - Xbeta.lambda
	sigma2.lambda <- sum(res.lambda^2)/n
	
	#-1/2 * ( n*log(sigma2.lambda) + n*sigma2.lambda/sigma2.lambda ) + (lambda - 1)*sum(log(y)) 
	-(n/2)*log(sigma2.lambda) - n/2 + (lambda.y - 1)*sum(log(y))
}

TBS.CI.lambda <- function( y, x, sep.lam = FALSE, intercept = TRUE, level = 0.95, interval = c(-3,3), lambda.init = c(0.5,0.5), beta.init = NULL ) { # not sure that this works for sep.lam = TRUE
	if(!sep.lam) { # single lambda
		prof.opt <- optimize( f = TBS.profllkhd.lambda, interval = interval, y = y, x = x, intercept = intercept, beta.init = beta.init, maximum = TRUE )
		lower <- uniroot( f = TBS.prof.shifted, interval = c( max(-2,interval[1]), prof.opt$maximum ), y = y, x = x, f0 = prof.opt$objective - qchisq( level, df = 1 )/2, intercept = intercept, beta.init = beta.init )$root
		upper <- uniroot( f = TBS.prof.shifted, interval = c( prof.opt$maximum, min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$objective - qchisq( level, df = 1 )/2, intercept = intercept, beta.init = beta.init )$root
		valname <- 'objective'
		maxname <- 'maximum'
	} else { # separate lambda
		prof.opt <- optim( lambda.init, fn = TBS.profllkhd.lambda, y = y, x = x, intercept = intercept, beta.init = beta.init, control = list( fnscale = -1 ) )
		print(prof.opt$convergence)
		lower <- c( uniroot( f = TBS.prof.shifted, interval = c( max(-2,interval[1]), prof.opt$par[1] ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 2, intercept = intercept, intval = interval, beta.init = beta.init )$root, uniroot( f = TBS.prof.shifted, interval = c( max(-2,interval[1]), prof.opt$par[2] ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 1, intercept = intercept, intval = interval, beta.init = beta.init )$root )
		upper <- c( uniroot( f = TBS.prof.shifted, interval = c( prof.opt$par[1], min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 2, intercept = intercept, intval = interval, beta.init = beta.init )$root, uniroot( f = TBS.prof.shifted, interval = c( prof.opt$par[2], min(2,interval[2]) ), y = y, x = x, f0 = prof.opt$val - qchisq( level, df = 1 )/2, idx.lamopt = 1, intercept = intercept, intval = interval, beta.init = beta.init )$root )
		valname <- 'value'
		maxname <- 'par'
	}
	list( lambda = prof.opt[[maxname]], lower = lower, upper = upper, loglik.opt = prof.opt[[valname]] )
}

# PTBS.llkhd.knot <- function( theta.toOpt, theta = NULL, y, x, intercept = TRUE, var.split = FALSE ) {
	# if( is.null(theta) ) {
		# theta <- theta.toOpt
	# } else {
		# theta[is.na(theta)] <- theta.toOpt
	# }
	# PTBS.llkhd.knot.full( theta, y, x, intercept, var.split )
# }

TBS.profllkhd.lamsepopt <- function( lambda.opt, lambda.fix, idx.opt = 1, y, x, intercept, beta.init ) {
	if( idx.opt == 1 )
		lambda <- c(lambda.opt,lambda.fix)
	else
		lambda <- c(lambda.fix,lambda.opt)
	TBS.profllkhd.lambda( lambda, y, x, intercept = intercept, beta.init = beta.init )
}

TBS.prof.shifted <- function( lambda.prof, y, x, f0, idx.lamopt = NULL, intercept, intval = c(-3,3), beta.init ) {
	if(is.null(idx.lamopt)) {
		TBS.profllkhd.lambda( lambda.prof, y, x, intercept, beta.init ) - f0
	} else {
		lambda.opt <- optimize( f = TBS.profllkhd.lamsepopt, interval = intval, lambda.fix = lambda.prof, idx.opt = idx.lamopt, y = y, x = x, intercept = intercept, beta.init = beta.init, maximum = TRUE )$maximum
		
		lambda <- if(idx.lamopt==1) { c(lambda.opt,lambda.prof) } else { c(lambda.prof,lambda.opt) }
		TBS.profllkhd.lambda( lambda, y, x, intercept ) - f0
	}
}

TBS.mle <- function( y, x, lambda = NULL, sep.lam = FALSE, interval = c(-3,3), intercept = TRUE, lambda.init = c(0.5,0.5), beta.init = NULL ) {
	n <- length(y)
	
	valname <- 'objective'
	
	if( is.null(dim(x)) ) {
		x <- matrix( x, nc = 1 )
	}
	X <- if(intercept) cbind( 1, x ) else x

	if(is.null(beta.init))
		beta.init <- my.lm( y, X )$coef

	if( is.null(lambda) & !sep.lam ) {
		ell.hat <- optimize( f = TBS.profllkhd.lambda, interval = interval, y = y, x = x, intercept = intercept, beta.init = beta.init, maximum = TRUE )
		lambda <- ell.hat$maximum
	} else if( is.null(lambda) & sep.lam ) {
		ell.hat <- optim( lambda.init, fn = TBS.profllkhd.lambda, y = y, x = x, intercept = intercept, beta.init = beta.init, control = list( fnscale = -1 ) )
		if( ell.hat$convergence > 0 ) print(ell.hat$convergence)
		valname <- 'value'
		lambda <- ell.hat$par
	} else {
		ell.hat <- list( objective = TBS.profllkhd.lambda( lambda, y = y, x = x, intercept = intercept, beta.init = beta.init ) )
	}
	opt.beta <- optim( beta.init, fn = TBS.RSS.beta, lambda = lambda, y = y, X = X, control = list( maxit = 5000 ) )
	if( opt.beta$conv > 0 ) print(opt.beta$conv)
	beta.lambda <- opt.beta$par
	
	if( length(lambda) == 1 ) 
		lambda.y <- lambda.x <- lambda
	else {
		lambda.y <- lambda[1]
		lambda.x <- lambda[2]
	}

	y.lambda <- BC.transform( lambda.y, y )
	Xbeta.lambda <- BC.transform( lambda.x, c(X %*% beta.lambda) )
	
	res.lambda <- y.lambda - Xbeta.lambda
	sigma2.lambda <- sum(res.lambda^2)/n  # MLE (biased)
		
	list( lambda.hat = lambda, beta.hat = beta.lambda, sigma2.hat = sigma2.lambda, ell.opt = ell.hat[[valname]] )
}

TBS.condmean.orig <- function( theta, x, x.beta = NULL ) { # x on original scale
	p <- length(theta) - 2
	if(is.null(x.beta))
		x.beta <- c( x %*% theta[1:p] )
	
	x.beta * ( 1 + ( theta[p+1]*(1-theta[p+2]) )/( 2* x.beta^(2*theta[p+2]) ) )
}

my.lm <- function( y, X ) {
	n <- length(y)
	p <- ncol(X)
	
	XtX.inv <- solve( crossprod(X) )
	beta.hat <- c( XtX.inv %*% t(X) %*% y )
	H <- X %*% XtX.inv %*% t(X)
	s2.hat <- c( y %*% ( diag(n) - H ) %*% y )/(n-p)
	e <- c( ( diag(n) - H ) %*% y )
	covmat <- s2.hat * XtX.inv 
	list( coefs = beta.hat, s2 = s2.hat, e = e, n = n, p = p, se.coefs = sqrt( diag(covmat) ), covbeta = covmat, levs = diag(H) )
}


# Region matrix for Anova
regionnames <- levels(factor(Contents$Region))
Region.mat <- t( sapply( Contents$Region, function(x){ as.numeric( x == regionnames ) } ) ) # for full data set
colnames(Region.mat) <- regionnames
rownames(Region.mat) <- Contents$L_ID


# Some non-random variables
lam.seq.DoL <- seq(-2,2,0.01)
theta.init.DoL <- c( 0, 1, 0.6, 0.5 )
xlam.seq.DoL <- seq(-4.5,0.2,0.025)

lam.seq.Loss <- seq(-2,1,0.01)
theta.init.Loss <- c( 2, 1, 4, 0.1 )
xlam.seq.Loss <- seq(5,180,0.5)

h.seq <- seq(0.001,0.045,0.001)


## Bootstrap indices
# R <- 200
# set.seed(395042)
# boot.idx <- matrix( nr = R, nc = n.obs )
# set.seed(58283)
# seeds <- sample( 1:6000000, 1000 )
# for( r in 1:R ) {
	# set.seed(seeds[r])
	# boot.idx[r,] <- sample( 1:n.obs, n.obs, replace = TRUE )
# }