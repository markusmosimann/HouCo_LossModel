########################################################################################################################
####### Predictions for cantons (non-random CV, transferability assessment according to Wenger and Olden (2012)) #######
########################################################################################################################
# a*(m+1) + b*m = N
# a + b = K --> b = K-a
# a*(m+1) + (K-a)*m = a*m + a + (K-a)*m = K*m + a = N
# a = N - K*m
# b = K - N + K*m = K*(m+1) - N
# Thus:
# N - K*floor(N/K) groups of size m+1 = ceiling(N/K)
# K*ceiling(N/K) - N groups of size m = floor(N/K)

# Group-size function:
  group.sizes <- function( N, k ) {
  	if( N%%k == 0 ) {
  	  rep( N/k, k ) 
  	} else { 
  	  c( rep( ceiling(N/k), N - k*floor(N/k) ), rep( floor(N/k), k*ceiling(N/k) - N ) ) 
    }
  }

# === General definitions ==============================================================================================
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Observations to omit for canton-based analysis ("aov")
  which( Contents$Region %in% c("GE","AI") )
  idx.omit <- c( 145, which( Contents$Region %in% c("GE","AI") ) )

# Vector of cantons
  region.aov <- Contents$Region[-idx.omit]

# Dummy matrix
  Region.mat.aov <- Region.mat[,-which(regionnames %in% c("GE","AI"))][-idx.omit,]

  n.obs.aov <- nrow(Region.mat.aov)
  n.obs.aov.cant <- as.vector( table( region.aov ) )

# insurance sum (for predictions based on DoL)
  y.aov.InSum <- Contents$InSum[-idx.omit]

# List of logical variables indicating the region (for technical reasons)
  rownames(Region.mat.aov) <- NULL
  Region.list.aov <- lapply( seq_len(ncol(Region.mat.aov)), function(i){ apply( Region.mat.aov, 2, as.logical )[,i] } )

tapply( x.aov.DoL, region.aov, range )

  # 2 folds (large separation):
  # {OW,TI}: 192
  # {SZ,UR,VS}: 167 + 19 = 186
  
  # 2 folds (better mixing):
  # {OW,SZ}: 187
  # {TI,UR,VS}: 191
  
# 10 folds (fold size ~40 except VS):
  # OW: 3 groups (size 37, 37, 36)
  # SZ: 2 groups (size 39, 38)
  # TI: 2 groups (size 41, 41)
  # UR: 2 groups (size 45, 45) # 3 groups if 11 folds
  # VS: 1 group (size 19) # exclude 20 more obs from training set (or does this introduce bias?)

# 5 folds (one per canton):
  # OW: 1 group (size 110)
  # SZ: 1 group (size 77)
  # TI: 1 group (size 82)
  # UR: 1 group (size 90)
  # VS: 1 group (size 19) # exclude ~70 more obs from training set (or does this introduce bias?)

# 20 folds (fold size ~20):
  # OW: 6 groups (size 19, 19, 18, 18, 18, 18)
  # SZ:	4 groups (size 20, 19, 19, 19)
  # TI:	4 groups (size 21, 21, 20, 20)
  # UR: 5 groups (size 18 each)
  # VS: 1 group (size 19)

# QUESTION: Is it correct to omit the AI and GE observations here???
# Is the relative total bias/ae now more relevant because the folds are non-random???
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ==================================================================================================================== #

# K = 2 large separation
  gr2.ls <- split( 1:n.obs.aov, apply( Region.mat.aov[,c(2,4:5)], 1, sum ) )
names(gr2.ls) <- c("OW,TI","SZ,UR,VS")

# K = 2 better mixed
  gr2.bm <- split( 1:n.obs.aov, apply( Region.mat.aov[,3:5], 1, sum ) )
  names(gr2.bm) <- c("OW,SZ","TI,UR,VS")
  
# K = 5 (by canton)
  gr5 <- apply( Region.mat.aov, 2, function(x){ which( x == 1 ) } )
  
# K = 10 (olny for DoL)
  gr10 <- unlist( mapply( function(n.gr,c.tf){ lapply( 1:n.gr, function(i){ 
    sort( (which(c.tf)[order( x.aov.DoL[c.tf] )])[seq(i,sum(c.tf),n.gr)] ) } ) }, n.gr = c(3,2,2,2,1), 
    c.tf = Region.list.aov ), recursive = FALSE )
  names(gr10) <- c( paste( rep( colnames(Region.mat.aov)[1:4], c(3,2,2,2) ), 
                           unlist( sapply( c(3,2,2,2), function(ng){ 1:ng } ) ), sep = '' ), 
                    colnames(Region.mat.aov)[5] )

# K = 20 (only for DoL)
gr20 <- unlist( mapply( function(n.gr,c.tf){ lapply( 1:n.gr, function(i){
  sort( (which(c.tf)[order( x.aov.DoL[c.tf] )])[seq(i,sum(c.tf),n.gr)] ) } ) }, n.gr = c(6,4,4,5,1),
  c.tf = Region.list.aov ), recursive = FALSE )
names(gr20) <- c( paste( rep( colnames(Region.mat.aov)[1:4], c(6,4,4,5) ), unlist( sapply( c(6,4,4,5), function(ng){
  1:ng } ) ), sep = '' ), colnames(Region.mat.aov)[5] )


##======================================================================##
##==== Create "joint" strata for both types of data ====================##
##======================================================================##

# rank(vec)[i] = rank of vec[i] among all elements of vec; same ordering as vec
# order(vec)[i] = index of the i-th smallest element of vec; ordering by increasing value of vec, i.e., as in sort(vec)

gr.DoL <- unlist( mapply( function(n.gr,c.tf){ lapply( 1:n.gr, function(i){ 
  sort( (which(c.tf)[order( x.aov.DoL[c.tf] )])[seq(i,sum(c.tf),n.gr)] ) } ) }, n.gr = c(3,2,2,2,1), 
  c.tf = Region.list.aov ), recursive = FALSE )
gr.Loss <- unlist( mapply( function(n.gr,c.tf){ lapply( 1:n.gr, function(i){
  sort( (which(c.tf)[order( x.aov.Loss[c.tf] )])[seq(i,sum(c.tf),n.gr)] ) } ) }, n.gr = c(3,2,2,2,1), 
  c.tf = Region.list.aov ), recursive = FALSE )

# How many common elements in each stratum under stratification into 2:20 strata?
# OW:
  sapply( sapply( 2:20, function(K){ gr.sizes <- group.sizes( n.obs.aov.cant[1], K ); 
  apply( rbind( c(0,cumsum(gr.sizes)[-K]), cumsum(gr.sizes) ), 2, function(gs){ 
    length( intersect( order(x.aov.DoL[Region.list.aov[[1]]])[(1+gs[1]):gs[2]], 
                       order(x.aov.Loss[Region.list.aov[[1]]])[(1+gs[1]):gs[2]] ) ) } )/gr.sizes } ), min )
# OW: Take 5 strata
# SZ:
  sapply( sapply( 2:20, function(K){ gr.sizes <- group.sizes( n.obs.aov.cant[2], K );
  apply( rbind( c(0,cumsum(gr.sizes)[-K]), cumsum(gr.sizes) ), 2, function(gs){ 
    length( intersect( order(x.aov.DoL[Region.list.aov[[2]]])[(1+gs[1]):gs[2]], 
                       order(x.aov.Loss[Region.list.aov[[2]]])[(1+gs[1]):gs[2]] ) ) } )/gr.sizes } ), min )
# SZ: Take 4 strata
# TI:
  sapply( sapply( 2:20, function(K){ gr.sizes <- group.sizes( n.obs.aov.cant[3], K ); 
  apply( rbind( c(0,cumsum(gr.sizes)[-K]), cumsum(gr.sizes) ), 2, function(gs){
    length( intersect( order(x.aov.DoL[Region.list.aov[[3]]])[(1+gs[1]):gs[2]], 
                       order(x.aov.Loss[Region.list.aov[[3]]])[(1+gs[1]):gs[2]] ) ) } )/gr.sizes } ), min )
# TI: Take 4 strata
# UR:
  sapply( sapply( 2:20, function(K){ gr.sizes <- group.sizes( n.obs.aov.cant[4], K ); 
  apply( rbind( c(0,cumsum(gr.sizes)[-K]), cumsum(gr.sizes) ), 2, function(gs){ 
    length( intersect( order(x.aov.DoL[Region.list.aov[[4]]])[(1+gs[1]):gs[2]], 
                       order(x.aov.Loss[Region.list.aov[[4]]])[(1+gs[1]):gs[2]] ) ) } )/gr.sizes } ), min )
# UR: Take 4 strata
  # In each stratum: put common observations, fill rest with obs having ranks closest to stratum
  nb.strata <- c(5,4,4,5)

i <- 4  # let this vary over 1:4
sapply( sapply( 2:20, function(K){ gr.sizes <- group.sizes( n.obs.aov.cant[i], K ); 
apply( rbind( c(0,cumsum(gr.sizes)[-K]), cumsum(gr.sizes) ), 2, function(gs){ 
  length( intersect( order(x.aov.DoL[Region.list.aov[[i]]])[(1+gs[1]):gs[2]], 
                     order(x.aov.Loss[Region.list.aov[[i]]])[(1+gs[1]):gs[2]] ) ) } )/gr.sizes } ), min )

strata <- list()
K <- nb.strata[i] # OW: 5, SZ: 4, TI: 4, UR: 5
gr.size <- group.sizes( n.obs.aov.cant[i], K )
gr.idx <- cumsum(gr.size)

ord.DoL <- order(x.aov.DoL[Region.list.aov[[i]]])
ord.Loss <- order(x.aov.Loss[Region.list.aov[[i]]])

common.elems <- apply( rbind( c(0,gr.idx[-K]), gr.idx ), 2, function(gs){ 
  intersect( ord.DoL[(1+gs[1]):gs[2]], ord.Loss[(1+gs[1]):gs[2]] ) } )  # in 1:n.obs.aov.cant
sapply( common.elems, length )
gr.size - sapply( common.elems, length )

# cbind( DoL = ord.DoL, Loss = ord.Loss )

ord.part.DoL <- ord.DoL[1:gr.idx[1]]
ord.part.Loss <- ord.Loss[1:gr.idx[1]]
common.new <- intersect( ord.part.DoL, ord.part.Loss )

#--- Look at stratum 1:
  gr.size[1] - length(common.new)
  (only.DoL <- setdiff( ord.part.DoL, ord.part.Loss ))
  (only.Loss <- setdiff( ord.part.Loss, ord.part.DoL ))
  sapply( only.DoL, function(x){ which( ord.Loss == x ) } )
  sapply( only.Loss, function(x){ which( ord.DoL == x ) } )
  gr.idx[1]

# OW: Pull 38 into stratum 1 and push 53.
  pull <- 38
  
# SZ: Pull c( 26, 58, 30 ) into stratum 1 and push c( 31, 1, 67 )
  pull <- c( 26, 58, 30 )
  
# TI: Pull c( 40, 9, 48, 28, 44 ) into stratum 1 and push c( 37, 51, 80, 1, 13 )
  pull <- c( 40, 9, 48, 28, 44 )

# UR: Pull 52 into stratum 1 and push 16.
  pull <- 52

strata[[1]] <- c( common.new, pull )
(push.DoL <- setdiff( only.DoL, pull ))
(push.Loss <- setdiff( only.Loss, pull ))

ord.part.DoL <- c( push.DoL, setdiff( ord.DoL[(gr.idx[1]+1):gr.idx[2]], pull ) )
ord.part.Loss <- c( push.Loss, setdiff( ord.Loss[(gr.idx[1]+1):gr.idx[2]], pull ) )
common.new <- intersect( ord.part.DoL, ord.part.Loss )

#---Look at stratum 2:
  gr.size[2] - length(common.new)
  (only.DoL <- setdiff( ord.part.DoL, ord.part.Loss ))
  (only.Loss <- setdiff( ord.part.Loss, ord.part.DoL ))
  sapply( only.DoL, function(x){ which( ord.Loss == x ) } )
  sapply( only.Loss, function(x){ which( ord.DoL == x ) } )
  gr.idx[2]

  # OW: Pull c( 30, 69, 91, 76 ) into stratum 2 and push c( 92, 47, 77, 4 )
    pull <- c( 30, 69, 91, 76 )
  # SZ: Pull c( 63, 65 ) into stratum 2 and push c( 13, 21 )
    pull <- c( 63, 65 )
  # TI: Pull c( 80, 79, 1, 13, 29 ) into stratum 2 and push c( 78, 21, 43, 82, 74 )
    pull <- c( 80, 79, 1, 13, 29 )
  # UR: Pull c( 7, 64 ) into stratum 2 (or ev. c(7,8)??) and push c( 58, 8 )
    pull <- c( 7, 64 )

  strata[[2]] <- c( common.new, pull )
  (push.DoL <- setdiff( only.DoL, pull ))
  (push.Loss <- setdiff( only.Loss, pull ))
  
  ord.part.DoL <- c( push.DoL, setdiff( ord.DoL[(gr.idx[2]+1):gr.idx[3]], pull ) )
  ord.part.Loss <- c( push.Loss, setdiff( ord.Loss[(gr.idx[2]+1):gr.idx[3]], pull ) )
  common.new <- intersect( ord.part.DoL, ord.part.Loss )

#---Look at stratum 3:
  gr.size[3] - length(common.new)
  (only.DoL <- setdiff( ord.part.DoL, ord.part.Loss ))
  (only.Loss <- setdiff( ord.part.Loss, ord.part.DoL ))
  sapply( only.DoL, function(x){ which( ord.Loss == x ) } )
  sapply( only.Loss, function(x){ which( ord.DoL == x ) } )
  gr.idx[3]

  # OW: Pull c( 41, 101 ) into stratum 3 and push c( 10, 42 )
   pull <- c( 41, 101 )
  # SZ: Pull c( 13, 21, 18 ) into stratum 3 and push c( 23, 15, 5 )
    pull <- c( 13, 21, 18 )
  # TI: Pull c( 43, 23, 34, 56, 24, 25 ) into stratum 3 and push c( 33, 70, 82, 39, 2, 14  )
    pull <- c( 43, 23, 34, 56, 24, 25 )
  # UR: Pull c( 33, 60, 49, 27 ) into stratum 3 and push c( 12, 42, 47, 83 )
    Spull <- c( 33, 60, 49, 27 )
  
  strata[[3]] <- c( common.new, pull )
  (push.DoL <- setdiff( only.DoL, pull ))
  (push.Loss <- setdiff( only.Loss, pull ))
  
  ord.part.DoL <- c( push.DoL, setdiff( ord.DoL[(gr.idx[3]+1):gr.idx[4]], pull ) )
  ord.part.Loss <- c( push.Loss, setdiff( ord.Loss[(gr.idx[3]+1):gr.idx[4]], pull ) )
  common.new <- intersect( ord.part.DoL, ord.part.Loss )

#---Look at stratum 4:
  gr.size[4] - length(common.new)
  	strata[[4]] <- common.new
  (only.DoL <- setdiff( ord.part.DoL, ord.part.Loss ))
  (only.Loss <- setdiff( ord.part.Loss, ord.part.DoL ))
  sapply( only.DoL, function(x){ which( ord.Loss == x ) } )
  sapply( only.Loss, function(x){ which( ord.DoL == x ) } )
  gr.idx[4]
  # OW: Pull c( 65, 93, 99, 24, 43, 94, 108 ) into stratum 4 and push c( 20, 110, 79, 45, 46, 5, 15 )
    pull <- c( 65, 93, 99, 24, 43, 94, 108 )
  # UR: Pull c( 72, 70, 18, 47 ) into stratum 4 and push c( 80, 74, 35, 28 )
    pull <- c( 72, 70, 18, 47 )

  strata[[4]] <- c( common.new, pull )
  (push.DoL <- setdiff( only.DoL, pull ))
  (push.Loss <- setdiff( only.Loss, pull ))

  ord.part.DoL <- c( push.DoL, setdiff( ord.DoL[(gr.idx[4]+1):gr.idx[5]], pull ) )
  ord.part.Loss <- c( push.Loss, setdiff( ord.Loss[(gr.idx[4]+1):gr.idx[5]], pull ) )
  common.new <- intersect( ord.part.DoL, ord.part.Loss )

#---Look at stratum 5:
  gr.size[5] - length(common.new)
	strata[[5]] <- common.new
#--------------------#
  sapply( strata, length ) - gr.size
  length(unique(unlist(strata))) - n.obs.aov.cant[i]
  
  strata.cant <- list( OW = strata )
  strata.cant$SZ <- strata
  strata.cant$TI <- strata
  strata.cant$UR <- strata
##########

# Global indices (1:n.obs.aov) for "joint" strata for DoL and Loss
  strata.cant.glob <- mapply( function(strat,reg){ lapply( strat, function(idx){ which(reg)[idx] } ) }, 
                              strat = strata.cant, reg = Region.list.aov[1:4] )
  save( list = c( "strata.cant", "strata.cant.glob" ), file = "JointStrata_nrCV.RData" )

    # #gr10: 
    #   sizes.lists <- list( OW = list(c(8,7,7),c(7,8,7),c(7,7,8),c(8,7,7),c(7,8,7)), SZ = list(c(10,10),c(10,9),c(9,10),c(10,9)),
    #                        TI = list(c(11,10),c(10,11),c(10,10),c(10,10)), UR = list(c(9,9),c(9,9),c(9,9),c(9,9),c(9,9)) )
    #   set.seed(237593)
    #   gr10 <- lapply( unlist( mapply( function(szs,strat){ 
    #     gr.ass <- lapply( szs, function(ss){
    #       sample( x = rep( 1:length(ss), ss ), size = sum(ss) ) 
    #       });
    #     spl <- mapply( split, x = strat, f = gr.ass ); 
    #     lapply( split( spl, f = seq(nrow(spl)) ), unlist ) 
    #   }, szs = sizes.lists, strat = strata.cant.glob ), recursive = FALSE ), sort )
    #   names(gr10) <- sub( pattern = ".", replacement = "", x = names(gr10), fixed = TRUE )
    #   gr10$VS <- which(Region.list.aov[[5]])
    #   
    #   # OW(3): c(8,7,7),c(7,8,7),c(7,7,8),c(8,7,7),c(7,8,7);
    #   # SZ(2): c(10,10),c(10,9),c(9,10),c(10,9);
    #   # TI(2): c(11,10),c(10,11),c(10,10),c(10,10);
    #   # UR(2): c(9,9),c(9,9),c(9,9),c(9,9),c(9,9)
    # 
    # #gr20:
    #   sizes.lists <- list( OW = list(c(4,4,4,4,3,3),c(4,4,4,3,4,3),c(4,4,3,4,3,4),c(4,3,4,3,4,4),c(3,4,3,4,4,4)), 
    #                        SZ = list(rep(5,4),c(5,4,5,5),c(5,5,4,5),c(5,5,5,4)), 
    #                        TI = list(c(6,5,5,5),c(5,6,5,5), rep(5,4),rep(5,4)), 
    #                        UR = list(c(4,4,4,3,3),c(4,4,3,4,3),c(4,3,4,3,4),c(3,4,3,4,4),c(3,3,4,4,4)) )
    #   set.seed(8374628)
    #   gr20 <- lapply( unlist( mapply( function(szs, strat){ gr.ass <- lapply( szs, function(ss){ 
    #     sample( x = rep( 1:length(ss), ss ), size = sum(ss) ) } ); spl <- mapply( split, x = strat, f = gr.ass ); 
    #     lapply( split( spl, f = seq(nrow(spl)) ), unlist ) }, szs = sizes.lists, strat = strata.cant.glob ), 
    #     recursive = FALSE ), sort )
    #    names(gr20) <- sub( pattern = ".", replacement = "", x = names(gr20), fixed = TRUE )
    #   gr20$VS <- which(Region.list.aov[[5]])
    #   
  # OW(6): c(4,4,4,4,3,3),c(4,4,4,3,4,3),c(4,4,3,4,3,4),c(4,3,4,3,4,4),c(3,4,3,4,4,4); 
  # SZ(4): rep(5,4),c(5,4,5,5),c(5,5,4,5),c(5,5,5,4);
  # TI(4): c(6,5,5,5),c(5,6,5,5),rep(5,4),rep(5,4);
  # UR(5): c(4,4,4,3,3),c(4,4,3,4,3),c(4,3,4,3,4),c(3,4,3,4,4),c(3,3,4,4,4);

# Mean ranks of strata
lapply( 1:4, function(i){ sapply( lapply( strata.cant[[i]], function(idx){ 
  ( rank(x.aov.DoL)[Region.list.aov[[i]]][idx] + rank(x.aov.Loss)[Region.list.aov[[i]]][idx] )/2 } ), range ) } )

K.vec <- c( 2, 2, 5, 10, 20 )
nb.K <- length(K.vec)

gr.list <- list( '2.ls' = gr2.ls, '2.bm' = gr2.bm, '5' = gr5, '10' = gr10, '20' = gr20 )

##=======================##
##==== Relative Loss ====##
##=======================##

pred.nrCV.DoL <- array( dim = c( n.obs.aov, 6, nb.K ), 
                        dimnames = list( NULL, c( "DoL.TBS.med", "DoL.PTBS.med", "DoL.PTBS.seplam.med", "DoL.TBS.mean", 
                                                  "DoL.PTBS.mean", "DoL.PTBS.seplam.mean" ), 
                                         names(gr.list) ) ) # n.obs by |{model:pred}| by nb.K

# n.obs by 6 by nb.K by K --> list of arrays (n.obs, 6, K) length nb.K
pred.nrCV.allobs.DoL <- lapply( K.vec, function(K){ array( dim = c( n.obs.aov, 6, K ), 
         dimnames = list( NULL, c( "DoL.TBS.med", "DoL.PTBS.med", "DoL.PTBS.seplam.med", 
                                   "DoL.TBS.mean", "DoL.PTBS.mean", "DoL.PTBS.seplam.mean" ),
                                                      NULL ) ) } )
		
# No adjustment for larger training set for VS!
for( it in 1:nb.K ) {
	K <- K.vec[it]
	groups <- gr.list[[it]]
	
	if( K == 5 ) {
		set.seed(382948)
		suppl.out.VS <- mapply( function( gr, n.out ){ sample( x = gr, size = n.out ) }, gr = groups[1:4], 
		                        n.out = c(21,15,16,18) )
	} else if( K == 10 ) {
		set.seed(46237)
		suppl.out.VS <- mapply( function( gr, n.out ){ sample( x = gr, size = n.out ) }, gr = groups[1:9], 
		                        n.out = c(2,2,2,2,2,3,2,3,2) )
	} else {
		suppl.out.VS <- NULL
	}
	
	for( k in 1:K ) {
	#=======#
		if(k < K){
		  withhold.k <- groups[[k]]
		} else if(k == K){
		  withhold.k <- c( groups[[k]], unlist(suppl.out.VS, use.names = FALSE) )
		}

		# Fit model to y[-groups[[k]]], giving estimated parameters theta.fit.
		# TBS model
		  mle.k <- TBS.mle( y = y.aov.DoL[-withhold.k], x = x.aov.DoL[-withhold.k], intercept = TRUE, beta.init = c(0.1,1) ) 
		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
		  fitted.k <- c( cbind( 1, x.aov.DoL ) %*% mle.k$beta.hat ) # [groups[[k]]]
  		TBS.med <- fitted.k * y.aov.InSum #[groups[[k]]] # predicted median
  		TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) * 
  		  y.aov.InSum #[groups[[k]]] # predicted mean
  		names(TBS.mean) <- NULL
		# PTBS model
  		# Same lambda
  		  mle.k <- PTBS.mle( y = y.aov.DoL[-withhold.k], x = x.aov.DoL[-withhold.k], intercept = TRUE, sep.lam = FALSE )
    		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
    		PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], 
    		                              c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.aov.DoL ) ) %*%
    		                                   mle.k$beta.hat ) ) * y.aov.InSum # [groups[[k]]] [groups[[k]]] # median prediction
    		PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform(
    		  mle.k$lambda.hat[1+0], x.aov.DoL ) ), sep.lam = FALSE ) * y.aov.InSum # [groups[[k]]] [groups[[k]]] # mean prediction
    		names(PTBS.mean) <- NULL
	    # Separate lambda
    		mle.k <- PTBS.mle( y = y.aov.DoL[-withhold.k], x = x.aov.DoL[-withhold.k], intercept = TRUE, sep.lam = TRUE )
    		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
    		PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform(
    		  mle.k$lambda.hat[1+1], x.aov.DoL ) ) %*% mle.k$beta.hat ) ) * y.aov.InSum # [groups[[k]]] [groups[[k]]] # median prediction
    		PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( 
    		  mle.k$lambda.hat[1+1], x.aov.DoL ) ), sep.lam = TRUE ) * y.aov.InSum # [groups[[k]]] [groups[[k]]] # mean prediction
    		names(PTBS.seplam.mean) <- NULL
	#=======#
		pred.nrCV.DoL[groups[[k]],,it] <- cbind( DoL.TBS.med = TBS.med[groups[[k]]], DoL.PTBS.med = PTBS.med[groups[[k]]], 
		                                         DoL.PTBS.seplam.med = PTBS.seplam.med[groups[[k]]],
		                                         DoL.TBS.mean = TBS.mean[groups[[k]]], DoL.PTBS.mean = PTBS.mean[groups[[k]]],
		                                         DoL.PTBS.seplam.mean = PTBS.seplam.mean[groups[[k]]] )
		
		# n.obs by |{model:pred}| by nb.K by K --> list of nb.K arrays (n.obs, |{model:pred}|, K)
		pred.nrCV.allobs.DoL[[it]][,,k] <- cbind( DoL.TBS.med = TBS.med, DoL.PTBS.med = PTBS.med, 
		                                          DoL.PTBS.seplam.med = PTBS.seplam.med, DoL.TBS.mean = TBS.mean, 
		                                          DoL.PTBS.mean = PTBS.mean, DoL.PTBS.seplam.mean = PTBS.seplam.mean )
	}
}
# Prediction error for whole data set
  Delta.nrCV.K.DoL <- aperm( apply( pred.nrCV.DoL, 2:3, err.measures, target = y.aov.DoL * y.aov.InSum ), 
                             c(2,1,3) ) # |{model:pred}| by |{error measures}| by nb.K


## Adjusted K-fold cross-validation
# Prediction errors for whole data set under the models fitted with one fold left out
  pred.err.nrCV.allobs.DoL <- lapply( pred.nrCV.allobs.DoL, function(arr){ 
    aperm( apply( arr, 2:3, err.measures, target = y.aov.DoL * y.aov.InSum ), c(2,1,3) ) } ) 
      # list of nb.K arrays of size |{model:pred}| by |{error measures}| by K

# Weight the previous errors over all folds for each K
  wght.pred.err.nrCV.allobs.DoL <- vapply( mapply( function(m.k, D.k){ 
    apply( m.k/n.obs.aov * aperm( D.k, c(3,1,2) ), 2:3, sum ) }, m.k = lapply( gr.list, sapply, length ), 
    D.k = pred.err.nrCV.allobs.DoL, SIMPLIFY = FALSE ), '[', FUN.VALUE = array( 0.1, dim = c(6,7) ) ) 
      # array |{model:pred}| by |{error measures}| by nb.K

# Adjusted K-fold prediction error
# Delta.AnrCV.K.DoL <- Delta.nrCV.K.DoL + array( Delta.app.DoL, dim = c(6,7,nb.K) ) - wght.pred.err.nrCV.allobs.DoL

# Delta.nrCV.K.DoL - Delta.AnrCV.K.DoL


##=======================##
##===== Absolute loss====##
##=======================##

pred.nrCV.Loss <- array( dim = c( n.obs.aov, 6, nb.K ), 
                         dimnames = list( NULL, c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", 
                                                   "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" ), 
                                          names(gr.list) ) ) # n.obs by |{model:pred}| by nb.K

# n.obs by 6 by nb.K by K --> list of arrays (n.obs, 6, K) length nb.K
pred.nrCV.allobs.Loss <- lapply( K.vec, function(K){ 
  array( dim = c( n.obs.aov, 6, K ), dimnames = list( NULL, c( "Loss.TBS.med", "Loss.PTBS.med", "Loss.PTBS.seplam.med", 
                                                               "Loss.TBS.mean", "Loss.PTBS.mean", "Loss.PTBS.seplam.mean" ),
                                                      NULL ) ) } )
		
# No adjustment for larger training set for VS!
for( it in 1:nb.K ) {
  K <- K.vec[it]
  groups <- gr.list[[it]]
  
  if( K == 5 ) {
    set.seed(382948)
    suppl.out.VS <- mapply( function( gr, n.out ){ sample( x = gr, size = n.out ) }, gr = groups[1:4], n.out = c(21,15,16,18) )
  } else if( K == 10 ) {
    set.seed(46237)
    suppl.out.VS <- mapply( function( gr, n.out ){ sample( x = gr, size = n.out ) }, gr = groups[1:9], n.out = c(2,2,2,2,2,3,2,3,2) )
  } else {
    suppl.out.VS <- NULL
  }
	
	for( k in 1:K ) {
	#=======#
		if(k < K)
			withhold.k <- groups[[k]]
		else if(k == K)
			withhold.k <- c( groups[[k]], unlist(suppl.out.VS, use.names = FALSE) )
		# Fit model to y[-groups[[k]]], giving estimated parameters theta.fit.
		### TBS model
		mle.k <- TBS.mle( y = y.aov.Loss[-withhold.k], x = x.aov.Loss[-withhold.k], intercept = TRUE, 
		                  beta.init = c(1000,0.5), interval = c(-3,1) ) 
		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
		fitted.k <- c( cbind( 1, x.aov.Loss ) %*% mle.k$beta.hat ) # [groups[[k]]]
		TBS.med <- fitted.k #[groups[[k]]] # predicted median
		TBS.mean <- TBS.condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x.beta = fitted.k ) #[groups[[k]]] # predicted mean
		names(TBS.mean) <- NULL
		### PTBS model
		mle.k <- PTBS.mle( y = y.aov.Loss[-withhold.k], x = x.aov.Loss[-withhold.k], intercept = TRUE, interval = c(-3,1), 
		                   sep.lam = FALSE )
		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
		PTBS.med <- BC.backtransform( mle.k$lambda.hat[1], c( cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.aov.Loss ) ) %*% 
		                                                        mle.k$beta.hat ) ) # [groups[[k]]] [groups[[k]]] # median prediction
		PTBS.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), x = cbind( 1, BC.transform( mle.k$lambda.hat[1+0], x.aov.Loss ) ),
		                            sep.lam = FALSE ) # [groups[[k]]] [groups[[k]]] # mean prediction
		names(PTBS.mean) <- NULL
		## Separate lambda
		mle.k <- PTBS.mle( y = y.aov.Loss[-withhold.k], x = x.aov.Loss[-withhold.k], intercept = TRUE, interval = c(-3,1), sep.lam = TRUE )
		# Predict y[groups[[k]]] by yhat( theta.fit, x[groups[[k]]] ).
		PTBS.seplam.med <- BC.backtransform( mle.k$lambda.hat[1], 
		                                     c( cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.aov.Loss ) ) %*% 
		                                          mle.k$beta.hat ) ) # [groups[[k]]] [groups[[k]]] # median prediction
		PTBS.seplam.mean <- condmean.orig( theta = unlist(mle.k[c(2,3,1)]), 
		                                   x = cbind( 1, BC.transform( mle.k$lambda.hat[1+1], x.aov.Loss ) ), 
		                                   sep.lam = TRUE ) # [groups[[k]]] [groups[[k]]] # mean prediction
		names(PTBS.seplam.mean) <- NULL
		
	#=======#
		pred.nrCV.Loss[groups[[k]],,it] <- cbind( Loss.TBS.med = TBS.med[groups[[k]]], Loss.PTBS.med = PTBS.med[groups[[k]]], 
		                                          Loss.PTBS.seplam.med = PTBS.seplam.med[groups[[k]]], 
		                                          Loss.TBS.mean = TBS.mean[groups[[k]]], Loss.PTBS.mean = PTBS.mean[groups[[k]]],
		                                          Loss.PTBS.seplam.mean = PTBS.seplam.mean[groups[[k]]] )
		
		# n.obs by |{model:pred}| by nb.K by K --> list of nb.K arrays (n.obs, |{model:pred}|, K)
		pred.nrCV.allobs.Loss[[it]][,,k] <- cbind( Loss.TBS.med = TBS.med, Loss.PTBS.med = PTBS.med, 
		                                           Loss.PTBS.seplam.med = PTBS.seplam.med, Loss.TBS.mean = TBS.mean, 
		                                           Loss.PTBS.mean = PTBS.mean, Loss.PTBS.seplam.mean = PTBS.seplam.mean )
	}
}
# Prediction error for whole data set
Delta.nrCV.K.Loss <- aperm( apply( pred.nrCV.Loss, 2:3, err.measures, target = y.aov.Loss ), c(2,1,3) ) 
  # |{model:pred}| by |{error measures}| by nb.K


# Adjusted K-fold cross-validation
  # Prediction errors for whole data set under the models fitted with one fold left out
  pred.err.nrCV.allobs.Loss <- lapply( pred.nrCV.allobs.Loss, function(arr){
    aperm( apply( arr, 2:3, err.measures, target = y.aov.Loss ), c(2,1,3) ) } ) 
      # list of nb.K arrays of size |{model:pred}| by |{error measures}| by K

# Weight the previous errors over all folds for each K
wght.pred.err.nrCV.allobs.Loss <- vapply( mapply( function(m.k, D.k){
  apply( m.k/n.obs.aov * aperm( D.k, c(3,1,2) ), 2:3, sum ) }, m.k = lapply( gr.list, sapply, length ), 
  D.k = pred.err.nrCV.allobs.Loss, SIMPLIFY = FALSE ), '[', FUN.VALUE = array( 0.1, dim = c(6,7) ) ) 
    # array |{model:pred}| by |{error measures}| by nb.K

# Adjusted K-fold prediction error
# Delta.AnrCV.K.Loss <- Delta.nrCV.K.Loss + array( Delta.app.Loss, dim = c(6,7,nb.K) ) - wght.pred.err.nrCV.allobs.Loss
# 
# Delta.nrCV.K.Loss - Delta.AnrCV.K.Loss

# Two transferability measures per model (K-fold & adj.K-fold CV)
# Wenger and Olden 2012: 1 number per fold

##=======================##
##======= Summary =======##
##=======================##

library(abind)
se.bias.transfer <- apply( abind( pred.nrCV.Loss[,1:3,], pred.nrCV.DoL[,1:3,], pred.nrCV.Loss[,4:6,], 
                                  pred.nrCV.DoL[,4:6,], along = 2 ) - y.aov.Loss, 2:3, sd )

se.rel.bias.transfer <- apply( 100*( abind( pred.nrCV.Loss[,1:3,], pred.nrCV.DoL[,1:3,], pred.nrCV.Loss[,4:6,], 
                                            pred.nrCV.DoL[,4:6,], along = 2 ) - y.aov.Loss )/y.aov.Loss, 2:3, sd )

se.mae.transfer <- apply( abs( abind( pred.nrCV.Loss[,1:3,], pred.nrCV.DoL[,1:3,], pred.nrCV.Loss[,4:6,], 
                                      pred.nrCV.DoL[,4:6,], along = 2 ) - y.aov.Loss ), 2:3, sd )

se.rel.ae.transfer <- apply( abs( 100*( abind( pred.nrCV.Loss[,1:3,], pred.nrCV.DoL[,1:3,], pred.nrCV.Loss[,4:6,], 
                                               pred.nrCV.DoL[,4:6,], along = 2 ) - y.aov.Loss )/y.aov.Loss ), 2:3, sd )

se.preds.transfer <- abind( se.rmspe = matrix( NA, 12, nb.K ), se.bias = se.bias.transfer, 
                            se.rel.tot.bias = matrix( NA, 12, nb.K ), se.rel.bias = se.rel.bias.transfer, 
                            se.mae = se.mae.transfer, se.rel.tot.ae = matrix( NA, 12, nb.K ), 
                            se.rel.ae = se.rel.ae.transfer, along = 3 )


# Summary plots for prediction on original scale:
# A: 2 large sep, V: 2 better mix, X: 5, O: 10, []: 20
  x11(width = 12, height = 9)
  par( mfrow = c(2,2), mar = c( 9, 4, 3.5, 4 ) )
  for( i in c(2,4,5,7) ) {
  	plot( rep( 1:12, 5 ) + rep( c(-0.36,-0.18,0,0.18,0.36), each = 12 ), 
  	      c( rbind( Delta.nrCV.K.Loss[1:3,i,], Delta.nrCV.K.DoL[1:3,i,], Delta.nrCV.K.Loss[4:6,i,], Delta.nrCV.K.DoL[4:6,i,] ) ),
  	      xlab = "", ylab = err.name[i], main = "Aggregate prediction error and variability", pch = rep( c(2,6,4,1,0), each = 12 ),
  	      xaxt = 'n', mgp = c(2.5,1,0), ylim = range( c( Delta.nrCV.K.Loss[1:3,i,], Delta.nrCV.K.DoL[1:3,i,], 
  	                                                     Delta.nrCV.K.Loss[4:6,i,], Delta.nrCV.K.DoL[4:6,i,]
  	                                                     # , Delta.AnrCV.K.Loss[1:3,i,], Delta.AnrCV.K.DoL[1:3,i,],
  	                                                     # Delta.AnrCV.K.Loss[4:6,i,], Delta.AnrCV.K.DoL[4:6,i,] 
  	                                                     ) ) )
    
  	# points( rep( 1:12, 5 ) + rep( c(-0.36,-0.18,0,0.18,0.36), each = 12 ), 
  	#         c( rbind( 
  	#           Delta.AnrCV.K.Loss[1:3,i,], Delta.AnrCV.K.DoL[1:3,i,], 
  	#           Delta.AnrCV.K.Loss[4:6,i,], Delta.AnrCV.K.DoL[4:6,i,] 
  	#           ) ), pch = rep( c(2,6,4,1,0), each = 12 ), col = 'orange' )
  	axis( 1, at = 1:12, labels = c( rownames(Delta.CV.Loss)[1:3], rownames(Delta.CV.DoL)[1:3], rownames(Delta.CV.Loss)[4:6],
  	                                rownames(Delta.CV.DoL)[4:6] ), las = 2 )
  	# axis( 1, at = 1:12, labels = c( sapply( strsplit( rownames(Delta.CV.Loss)[1:6], ".me" ), '[', 1 ), 
  	#   sapply( strsplit( rownames(Delta.CV.DoL)[1:6], ".me" ), '[', 1 ) ), las = 2 )
  	abline( h = 0, col = 'red' )
  	abline( v = seq(0.5,12.5,1), lty = 3 )
  	if(i!=1) {
  		par(new = TRUE, cex.axis=0.9, mar = c( 9.5, 4, 3.5, 4 ))
  		plot( rep( 1:12, each = 6 ) + rep( c(-0.36,-0.18,0,0.18,0.36,NA), 12 ), c( t(cbind(se.preds.transfer[,,i],NA)) ),
  		      axes = FALSE, xlab = "", ylab = "", type = 'l', lty = 6, col = 'blue' )
  		axis( 4, at = pretty(range(se.preds.transfer[,,i])), col = 'blue', col.axis = 'blue' )
  		mtext( paste( "Standard deviation in the sample [", ifelse( i==2 | i==5, "CHF", "%"), "]", sep = '' ), side = 4, 
  		       line = 2.5, col= 'blue', cex = par("cex") )
  	}
  }
  savePlot( paste( "Figures/nrCV_", modelname, "_ModMetrics", sep = '' ), type = "pdf" )

# A: 2 large sep, V: 2 better mix, X: 5, O: 10, []: 20

# General patterns in aggregate errors for transferability assessment:
# Bias smaller for mean than for median; se smaller for DoL than for Loss
# Absolute error smaller for DoL than for Loss and slightly for median than for mean; se smaller for DoL than for Loss
# Relative bias and relative absolute error smaller for median than for mean and (slightly) smaller for Loss than for DoL; se same pattern
# RMSPE smaller for DoL than for Loss; no se available


# Ranks of the models in terms of different prediction errors
# Bias: 11, 12, 10, 7, 9, 8, 1, 6, 3, 5, 2, 4
# Relative bias: 1, 3, 1, 4, 6, 5, 7, 9, 7, 10, 12, 11
# MAE: 7, 9, 9, 1, 4, 2, 8, 11, 12, 5, 3, 6
# Relative AE: 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11
# RMSE: 8, 11, 10, 2, 4, 2, 7, 9, 12, 6, 1, 4

err.ranks.transfer <- rbind( bias = c( 11, 12, 10, 7, 9, 8, 1, 6, 3, 5, 2, 4 ), 
                             rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ),
                             mae = c( 7, 9.5, 9.5, 1, 4, 2, 8, 11, 12, 5, 3, 6 ), 
                             rel.ae = c( 1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11 ),
                             rmse = c( 8, 11, 10, 2.5, 4, 2.5, 7, 9, 12, 6, 1, 4 ) )

apply( err.ranks.transfer, 2, sum ) # 4
apply( err.ranks.transfer == apply( err.ranks.transfer, 1, min ), 2, sum ) # 1
apply( err.ranks.transfer == apply( err.ranks.transfer, 1, max ), 2, sum ) # not 9, 11, 2
apply( err.ranks.transfer, 2, max ) # 4
apply( err.ranks.transfer, 2, min ) # 1, 3, 4, 7, 11

# Ranks of the models in terms of the std deviation for different prediction errors
# Bias: 8, 10.5, 10.5, 2.5, 4, 2.5, 7, 9, 12, 6, 1, 5
# Relative bias: 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11
# MAE: 8, 10.5, 10.5, 2, 5, 3.5, 7, 9, 12, 6, 1, 3.5
# Relative AE: 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11

se.ranks.transfer <- rbind( se.bias = c( 8, 10.5, 10.5, 2.5, 4, 2.5, 7, 9, 12, 6, 1, 5 ), 
                            se.rel.bias = c( 1.5, 3, 1.5, 4, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ), 
                            se.mae = c( 8, 10.5, 10.5, 2, 5, 3.5, 7, 9, 12, 6, 1, 3.5 ), 
                            se.rel.ae = c( 1.5, 3.5, 1.5, 3.5, 6, 5, 7.5, 9, 7.5, 10, 12, 11 ) )

apply( se.ranks.transfer, 2, sum ) # 4
apply( se.ranks.transfer == apply( se.ranks.transfer, 1, min ), 2, sum ) # 1, 3, 11
apply( se.ranks.transfer == apply( se.ranks.transfer, 1, max ), 2, sum ) # not 9, 11
apply( se.ranks.transfer, 2, max ) # 4
apply( se.ranks.transfer, 2, min ) # 1, 3, 11, not 7:9

