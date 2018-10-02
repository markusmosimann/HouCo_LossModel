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
library(dunn.test)
library(SuppDists)
library(partitions)
library(abind)
# 1. Set-up ====================================================================================================
# Inherit scientific notification of numbers
options("scipen"=100)

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

# PTBS plots:
source("Scripts/S_PlotMainFig.R")
    
# Stats:
  # CI betas:
    par.DoL[1:2] + qnorm( 0.975 ) * sqrt(diag(covML.DoL)[1:2]) %*% t( c(-1,1) ) 
    par.Loss[1:2] + qnorm( 0.975 ) * sqrt(diag(covML.Loss)[1:2]) %*% t(c(-1,1) )


  # Correlations (non-parametric):
    cor(y.lambda.DoL, x.lambda.DoL, method = "s"); cor(y.lambda.DoL, x.lambda.DoL, method = "k")
    cor(y.lambda.DoL, x.lambda.Loss, method = "s"); cor(y.lambda.DoL, x.lambda.Loss, method = "k")
    
  # Confidence interval lambda:
    CI.lambda.DoL
    CI.lambda.Loss
    
  # MLE of lambda and beta parameters:
    mle.DoL
    mle.Loss
  
  # Sigma:
    sqrt(par.DoL[3])
    sqrt(par.Loss[3])

  # adjusted R2 of regression function:
    adj_R2.Loss
    adj_R2.DoL
    
  # Shapiro-test for normality assumption:
    shapiro.test(x = r.DoL)
    shapiro.test(x = r.Loss)  
    
  # linear model for bptest
    lm.DoL <- lm(y.lambda.DoL ~ x.lambda.DoL)
    lm.Loss <- lm(y.lambda.Loss ~ x.lambda.Loss)
    # Breusch-Pagan test for heteroscedasticity:
      bptest(lm.DoL)
      bptest(lm.Loss)






  