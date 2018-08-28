# R-Script Vulnerability of household contents =========================================================================
# Date:       05.02.2018
# Author:     Markus Mosimann, markus.mosimann@giub.unibe.ch
# Data model: Same for contents and structure (loss on contents + loss on structure = building loss)
#             Variables:
#               ID [INT]: 
#                 > Number or expression to link content data with structure. The data of contents and structure
#                   must be ordered to this number!
#               Date [Date, format YYYY-MM-DD]:
#                 > Date of the flood event, e.g. 2005-08-22
#               Region [Factor]:
#                 > Factor which levels define different regions (non-random cross-validation)
#               Value [INT, > 0]:
#                 > Monetary value of assets
#               s_tot [INT, > 0]:
#                 > Damage, including sum paid out by insurance AND excess. Same unit as "value"
#               DoL [Double, > 0 (and < 1)]:
#                 > Degree of loss (ratio of s_tot on value). DoL gets higher than 1, as soon as clean up work is
#                   included within the loss
#
#=======================================================================================================================
library(lmtest)
library(RColorBrewer)
# 1. Set-up ====================================================================================================
# Inherit scientific notification of numbers
options("scipen"=100)

# Tabelle mit wichtigsten Werten:
paramTable <- as.data.frame(array(data = NA, dim = c(2, 12)))
colnames(paramTable) <- c("Model", "Sp.rho / K.tau", "lambda CI", "lambda Est", "Sigma", "Beta_0", "Beta_0_CI", 
                          "Beta_1", "Beta_1_CI", "R2", "SW_p-value", "BP_p-value")
paramTable$Model <- c("DoL", "Loss")

# load data
load("Data/Loss_Data.RData")

# Load model functions
source("Scripts//f_LossAnalysis.R")


# 2. Data distribution =================================================================================================
# Define the regionnames with enough Data to produce meaningful boxplots:
RegionsToPlot  <- levels(as.factor(c("OW", "TI", "UR", "SZ", "VS")))
sapply(RegionsToPlot, FUN = function(Reg){
  ind <- which(Contents$Region==Reg)
  print(paste(Reg, "Fraction of Content loss to total loss:", round(sum(Contents$Loss[ind]) /
                                                                      sum(Contents$Loss[ind] + Structure$Loss[ind]),4)))
  print(paste(Reg, "Mean Fraction:", round(mean(Contents$Loss[ind]/(Contents$Loss[ind] + Structure$Loss[ind])), 4)))
  print(paste(Reg, "Median Fraction:", round(median(Contents$Loss[ind]/(Contents$Loss[ind] + Structure$Loss[ind])), 4)))
})

sapply(RegionsToPlot, FUN = function(Reg){
  ind <- which(Contents$Region==Reg)
  print(paste(Reg, "Mean ratio of Content:Structure vulnerability:", round(mean(Contents$DoL[ind]/Structure$DoL[ind]), 4)))
  print(paste(Reg, "Median ratio of Content:Structure vulnerability:", round(median(Contents$DoL[ind]/Structure$DoL[ind]), 4)))
})


# Data distribution figures: 
source("Scripts/S_DataDist.R")

# 3. Data analysis =====================================================================================================

# Relative loss model:
source("Scripts/S_RelLossAnalysis.R")

# absolute loss model loss model:
source("Scripts/S_AbsLossAnalysis.R")

# (Leave-one-out-) cross-validation:
source("Scripts/S_CrossVal.R")

# Analysis of variance:
source("Scripts/S_anova.R")

# non-random cross-validation:
source("Scripts/S_nonrandCrossVal.R")

# 4. PTBS plots ========================================================================================================
# CI/PI method 2 was used (check S_RelLossAnalysis.R/S_AbsLossAnalysis.R)
  
  pdf( file = paste( "Figures/Res_", modelname, "_overview.pdf", sep = '' ), width = 11, height = 12)
  # x11( width = 11, height = 12)
  par(cex.axis=1.5, cex.lab=1.6, cex.main=1.6, cex.sub=1.5, mfrow=c(2, 2),  mar = c(8, 4.2, 3.5, 0)+.1)
  
  ##===========##
  ##=== DoL ===##
  ##===========##
  
  ##===== Transformed scale =====##
    # y.lambda vs. X.lambda*beta
      plot( x.lambda.DoL, y.lambda.DoL, xlab = "Structure", ylab = "Contents", pch=20, cex=.5)
      title(main = "a) transformed scale", line = 0.5, cex.main=1.2 )
  
  
    # 95% Prediction interval (mle, lambda known)
      pred95.lwr <- fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL )
      pred95.upr <- fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.DoL )
      polygon(c(xlam.seq.DoL, rev(xlam.seq.DoL)), c(pred95.lwr, rev(pred95.upr)), 
              col="lightblue", border = NA, density = 12, angle = c(-45, 45))
      
    # 50% Prediction interval (mle, lambda known)
      pred95.lwr <- fitted.mle.seq.DoL - qt( 0.75, df = n.obs - p )*sqrt( var.predict.seq.DoL )
      pred95.upr <- fitted.mle.seq.DoL + qt( 0.75, df = n.obs - p )*sqrt( var.predict.seq.DoL )
      polygon(c(xlam.seq.DoL, rev(xlam.seq.DoL)), c(pred95.lwr, rev(pred95.upr)), 
              col="skyblue", border = NA, density = 24, angle = c(-45, 45))
      
    # 95% Confidence interval (mle, lambda known)
      conf.lwr <- fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL )
      conf.upr <- fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.DoL )
      polygon(c(xlam.seq.DoL, rev(xlam.seq.DoL)), c(conf.lwr, rev(conf.upr)), col="lightblue", border = NA)

    # Median fit:
      abline( mle.DoL$beta.hat, col="darkblue", lty=2, lwd=2) 
      
    # if appropriate:
      points( BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], Structure$DoL[idx.outl] ),
              BC.transform( mle.DoL$lambda.hat[1], Contents$DoL[idx.outl] ), pch = 4)
    # X-Y-Line:
      abline(a = 0, b = 1, col="darkgrey", lty=4, lwd=2)
      
    # Out-range
      polygon( x = BC.transform( mle.DoL$lambda.hat[1+sep.lam.global], c( 1.5, 1, 1, 0, 0, 1.5 ) ), 
               y = BC.transform( mle.DoL$lambda.hat[1], c( 0, 0, 1, 1, 1.5, 1.5 ) ), density = 5, angle = -45, 
               col = "lightgrey")
    
    # add legend:
      inset.tmp <- 0.1
      format(adj_R2.DoL, digits=2) # change in adj. R^2 if different!
      legend("topleft", bty="n", bg = "white", inset = c(-0.05, inset.tmp),
             legend=c(expression(paste("adj. R"^2,": 0.67" ))),
             cex=1.2
      ); 
      points(x.lambda.DoL, y.lambda.DoL, pch=20, cex=.5)
      
  
  ##===== Backtransformed scale =====##
    par(mar = c(8, 2.2, 3.5, 2)+.1)
    # y.lambda vs. X.lambda*beta
      plot( x.DoL, y.DoL, xlab = "Structure", ylab = "Contents", pch=20, cex=.5)
      title(main = "b) original scale", line = 0.5, cex.main=1.2 )
    
      
    # 95% Prediction interval (mle, lambda known)
      pred95.lwr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) )
      pred95.upr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) )
      polygon(c(x.seq.DoL, rev(x.seq.DoL)), c(pred95.lwr, rev(pred95.upr)), 
              col="lightblue", border = NA, density = 12, angle = c(-45, 45))
      
    # 50% Prediction interval (mle, lambda known)
      pred50.lwr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.75, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) )
      pred50.upr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.75, df = n.obs - p )*
                                        sqrt( var.predict.seq.DoL ) )
      polygon(c(x.seq.DoL, rev(x.seq.DoL)), c(pred50.lwr, rev(pred50.upr)), 
              col="skyblue", border = NA, density = 24, angle = c(-45, 45))
      
    # 95% Confidence interval (mle, lambda known)
      conf.lwr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL - qt( 0.975, df = n.obs - p )*
                                      sqrt( var.fitted.seq.DoL ) )
      conf.upr <- BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL + qt( 0.975, df = n.obs - p )*
                                      sqrt( var.fitted.seq.DoL ) )
      polygon(c(x.seq.DoL, rev(x.seq.DoL)), c(conf.lwr, rev(conf.upr)), col="lightblue", border = NA)
      
    # Median fit:
      lines( x.seq.DoL, BC.backtransform( mle.DoL$lambda.hat[1], fitted.mle.seq.DoL ), col="darkblue", lty=2, lwd=2 )
      
    # if appropriate:
      points( Structure$DoL[idx.outl], Contents$DoL[idx.outl], pch = 4)
      
    # X-Y-Line:
      abline(a = 0, b = 1, col="darkgrey", lty=4, lwd=2)
      
    # Out-range
      polygon(x = c( 1.5, 1, 1, -1, -1, 1.5 ), y =  c( -1, -1, 1, 1, 1.5, 1.5 ), density = 5, angle = -45, 
              col = "lightgrey")

      points(x.DoL, y.DoL, pch=20, cex=.5)
      
      # Conditional mean on original scale:
      lines( x.seq.DoL, condmean.seq.DoL, col = 'green3', lwd = 1 )
      lines( x.seq.DoL, condmean.seq.DoL - qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
      lines( x.seq.DoL, condmean.seq.DoL + qnorm(0.95)*sqrt(var.condmean.DoL), col = 'green3', lty = 2 )
   
    
  ##============##
  ##=== LOSS ===##
  ##============##
      
  par(mar = c(5, 4.2, 6.5, 0)+.1)
  ##===== Transformed scale =====##
    # y.lambda vs. X.lambda*beta
      plot(x.lambda.Loss, y.lambda.Loss, xlab = "Structure", ylab = "Contents", pch=20, cex=.5)
      title(main = "c) transformed scale", line = 0.5, cex.main=1.2 )
      
    
    # 95% Prediction interval (mle, lambda known)
      pred95.lwr <- fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss )
      pred95.upr <- fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.predict.seq.Loss )
      polygon(c(xlam.seq.Loss, rev(xlam.seq.Loss)), c(pred95.lwr, rev(pred95.upr)), 
              col="lightblue", border = NA, density = 12, angle = c(-45, 45))
    
    # 50% Prediction interval (mle, lambda known)
      pred50.lwr <- fitted.mle.seq.Loss - qt( 0.75, df = n.obs - p )*sqrt( var.predict.seq.Loss )
      pred50.upr <- fitted.mle.seq.Loss + qt( 0.75, df = n.obs - p )*sqrt( var.predict.seq.Loss )
      polygon(c(xlam.seq.Loss, rev(xlam.seq.Loss)), c(pred50.lwr, rev(pred50.upr)), 
              col="skyblue", border = NA, density = 24, angle = c(-45, 45))
    
    # 95% Confidence interval (mle, lambda known)
      conf.lwr <- round(fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), 4)
      conf.upr <- round(fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*sqrt( var.fitted.seq.Loss ), 4)
      polygon(c(xlam.seq.Loss, rev(xlam.seq.Loss)), c(conf.lwr, rev(conf.upr)), col="lightblue", border = NA)
    
    # Median fit:
     abline( mle.Loss$beta.hat, col="darkblue", lty=2, lwd=2) 
    
    # if appropriate:
      points( BC.transform( mle.Loss$lambda.hat[1+sep.lam.global], Structure$Loss[idx.outl] ),
              BC.transform( mle.Loss$lambda.hat[1], Contents$Loss[idx.outl] ), pch = 4)
    # X-Y-Line:
     abline(a = 0, b = 1, col="darkgrey", lty=4, lwd=2)
    
    
    # add legend:
      inset.tmp <- 0.1
      format(adj_R2.Loss, digits=2) # change in adj. R^2 if different!
      legend("topleft", bty="n", bg = "white", inset = c(-0.05, inset.tmp),
             legend=c(expression(paste("adj. R"^2,": 0.62" ))),
             cex=1.2
      ); 
      points(x.lambda.Loss, y.lambda.Loss, pch=20, cex=.5)
      
    
  ##===== Backtransformed scale =====##
    par(mar = c(5, 2.2, 6.5, 2)+.1)
    # y.lambda vs. X.lambda*beta
      plot( x.Loss, y.Loss, xlab = "Structure", ylab = "Contents", pch=20, cex=.5,  xaxt="n", yaxt="n")
      axis(1, at = axTicks(1)[-1], labels = format(axTicks(1)[-1],big.mark = " "))
      axis(2, at = axTicks(2)[-1], labels = format(axTicks(2)[-1],big.mark = " "), mgp=c(3, .7, 0))
      axis(1, at = axTicks(1)[1], labels = axTicks(1)[1])
      axis(2, at = axTicks(2)[1], labels = axTicks(2)[1])
      title(main = "d) original scale", line = 0.5, cex.main=1.2 )
      
    # 95% Prediction interval (mle, lambda known)
      pred95.lwr <- BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) )
      pred95.upr <- BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) )
      polygon(round(c(x.seq.Loss, rev(x.seq.Loss)),0), round(c(pred95.lwr, rev(pred95.upr)),0), 
              col="lightblue", border = NA, density = 12, angle = c(-45, 45))
      
    # 50% Prediction interval (mle, lambda known)
      pred50.lwr <- round(BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.75, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) ), 0)
      pred50.upr <- round(BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.75, df = n.obs - p )*
                                        sqrt( var.predict.seq.Loss ) ), 0)
      polygon(c(x.seq.Loss, rev(x.seq.Loss)), c(pred50.lwr, rev(pred50.upr)), col="skyblue", border = NA, density = 24, angle = c(-45, 45))
      
    # 95% Confidence interval (mle, lambda known)
      conf.lwr <- round(BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss - qt( 0.975, df = n.obs - p )*
                                      sqrt( var.fitted.seq.Loss ) ), 0)
      conf.upr <- round(BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss + qt( 0.975, df = n.obs - p )*
                                      sqrt( var.fitted.seq.Loss ) ), 0)
      polygon(c(x.seq.Loss, rev(x.seq.Loss)), c(conf.lwr, rev(conf.upr)), col="lightblue", border = NA)
    
    # Median fit:
      lines( x.seq.Loss, BC.backtransform( mle.Loss$lambda.hat[1], fitted.mle.seq.Loss ), col="darkblue", lty=2, lwd=2 )
    
    # if appropriate:
      points( Structure$Loss[idx.outl], Contents$Loss[idx.outl], pch = 4)
    
    # X-Y-Line:
      abline(a = 0, b = 1, col="darkgrey", lty=4, lwd=2)
    
      points(x.Loss, y.Loss, pch=20, cex=.5)
    
    # Conditional mean on original scale:
      lines( x.seq.Loss, condmean.seq.Loss, col = 'green3', lwd = 1 )
      lines( x.seq.Loss, condmean.seq.Loss - qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )
      lines( x.seq.Loss, condmean.seq.Loss + qnorm(0.95)*sqrt(var.condmean.Loss), col = 'green3', lty = 2 )
    
    
  ##===================##
  ##=== Main Legend ===##
  ##===================##
    
    par(fig=c(0,1,0,1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=T)
    plot(0, 0, type = "n", bty="n", xaxt="n", yaxt="n")
    title(main = "Degree of loss [-]",line = -2.7)
    title(main = "Monetary loss [CHF]", outer = T, line = -42)
    
    # lines:
      legend(-0.9, 0, ncol = 3, cex=1.5, bg = NA, border = NA, bty = "n",
             legend=c("X = Y"), lty = 4, col="darkgrey")
      legend(-0.2, 0, ncol = 3, cex=1.5, bg = NA, border = NA, bty = "n",
             legend=c("Median"), lty = 2, col="darkblue")
      legend(0.3, 0, ncol = 3, cex=1.5, bg = NA, border = NA, bty = "n",
             legend=c("Cond. mean (+ 95% CI)"), lty = 1, col="green3")
    # areas:
      legend(-0.85, .1, ncol = 3, cex=1.5, bg = NA, bty = "n", border = NA, col = "white", x.intersp = -1,
             density = NA, angle = NA, legend = "95% CI", fill = "lightblue", lty = 2)
      legend(-.15, .1, ncol = 3, cex=1.5, bg = NA, bty = "n", border = NA, col = "white", x.intersp = -1,
             density = 24, angle = -45, legend= "IQR", fill = "lightblue", lty = 4)
      legend(.35, .1, ncol = 3, cex=1.5, bg = NA, bty = "n", border = NA, col = "white", x.intersp = -1,
             density = 12, angle = -45, legend = "95% PI", fill = "skyblue", lty = 1)
    
  dev.off()
    