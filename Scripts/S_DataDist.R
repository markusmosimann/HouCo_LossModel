############################################## Chapter 2: Data & Methods ###############################################
# 2.1 Data: ============================================================================================================
pdf(file = "Figures/Figure_2_RatiosAll.pdf", width = 7, height = 6)
par(mfrow=c(1,2))
# ~ Data Distribution --------------------------------------------------------------------------------------------------
  # Share of content damage on total building damage (content + building structure):
  par(mgp = c(3, 3,0))
  boxplot(Contents$Loss/(Structure$Loss+Contents$Loss), ylab = "Share of household content on total", main = "a)",
          names = c(expression(frac("L"["hc"], "L"["tot"])), expression(frac("V"["hc"], "V"["tot"]))), 
          Contents$InSum/(Structure$InSum+Contents$InSum), ylim=c(0, 1), notch = T, yaxt="n")
  par(mgp = c(3, 1, 0))
  axis(2)
  abline(h = seq(0,20,0.1), lty=3, col="lightgrey")
  points(c(mean(Contents$Loss/(Structure$Loss+Contents$Loss)),  
           mean(Contents$InSum/(Structure$InSum+Contents$InSum))), col = "red", pch = 9, cex = 2)
  # legend("topright", legend = "Mean", pch = 9, col = "red")
  
  # Show Ratio of DoL measured on contents relative to the DoL of structure:
  boxplot(Contents$DoL/Structure$DoL, ylim=c(0, 10), notch = T, ylab="Ratio contents to structure degrees of loss",
          xlab = c(expression(frac("DoL"["hc"], "DoL"["bs"]))), main = "b)")
  abline(h = seq(0,20,1), lty=3, col="lightgrey")
  points(mean(Contents$DoL/Structure$DoL), col="red", pch=9, cex=2)
  dev.off()

# ~ Distribution over Regions ------------------------------------------------------------------------------------------
  # Boxplots
  EntriesPerReg <- table(Contents$Region)
  RegionLabels   <- sapply(RegionsToPlot, FUN = function(x){paste(x, "\n n:", EntriesPerReg[x])})
  
  pdf(file = "Figures/Figure_1_Regions.pdf", width = 14, height = 7)
  par(cex.axis=1.5, cex.lab=1.8, mfrow=c(1,2), mar=c(5, 4.5, 4, 1)+.1)
  lg.tmp <- Contents$Region %in% RegionsToPlot
  
  boxplot(Contents$DoL[lg.tmp] ~ Contents$Region[lg.tmp], at = c(1:length(RegionsToPlot)-0.2), varwidth = T, xaxt="n",
          log ="y", col = "#f1a340", ylim = c(0.0001, 1), boxwex = 0.3, ylab=c("Degree of loss [-]"))
  boxplot(Structure$DoL[lg.tmp] ~ Structure$Region[lg.tmp], at = c(1:length(RegionsToPlot)+0.2), varwidth = T, xaxt="n",
          log ="y", col = "#998ec3", boxwex = 0.3 ,add=T)
  axis(1, at = c(1:5), tick = T, labels = RegionLabels, mgp=c(3, 2.7, 0))
  grid(nx = NA, ny = NULL)
  legend("bottomright", fill= c("#f1a340", "#998ec3"), legend= c("Contents", "Structure"), cex = 1.8)
  
  boxplot(Contents$Loss[lg.tmp] ~ Contents$Region[lg.tmp], at = c(1:length(RegionsToPlot)-0.2), varwidth = T, xaxt="n",
          log ="y",  col = "#f1a340", ylim = c(100, 1000000), boxwex = 0.3, ylab = "Monetary loss [CHF]")
  boxplot(Structure$Loss[lg.tmp] ~ Structure$Region[lg.tmp], at = c(1:length(RegionsToPlot)+0.2), varwidth = T,
          log ="y", xaxt="n", col = "#998ec3", ylim = c(100, 1000000), boxwex = 0.3 ,add=T)
  axis(1, at = c(1:5), tick = T, labels = RegionLabels, mgp=c(3, 2.7, 0))
  grid(nx = NA, ny = NULL)
  print("Wilcoxon test for significant differences between contents and structure:")
  # Untersuche, ob signifikante Unterschiede zwischen Fahrhabe und Structure (Wilcoxon, gepaart):
  tmp <- sapply(RegionsToPlot, function(x){
    ind.tmp <- which(Contents$Region==x)
    rbind(round(wilcox.test(Contents$Loss[ind.tmp], Structure$Loss[ind.tmp])$p.value, 5), 
          round(wilcox.test(Contents$DoL[ind.tmp], Structure$DoL[ind.tmp])$p.value, 5))
  }); print(data.frame(tmp, row.names = c("Loss", "DoL")))
  dev.off()
  
  boxplot(Contents$InSum[lg.tmp] ~ Contents$Region[lg.tmp], at = c(1:length(RegionsToPlot)-0.2), varwidth = T, xaxt="n",
          log ="y", col = "#f1a340", boxwex = 0.3, ylab=c("Degree of loss [-]"), ylim = c(20000, 2000000), notch=T)
  boxplot(Structure$InSum[lg.tmp] ~ Structure$Region[lg.tmp], at = c(1:length(RegionsToPlot)+0.2), varwidth = T, xaxt="n",
          log ="y", col = "#998ec3", boxwex = 0.3 ,add=T, notch=T)
  axis(1, at = c(1:5), tick = T, labels = RegionLabels, mgp=c(3, 2.5, 0))
  grid(nx = NA, ny = NULL)
  legend("bottomright", fill= c("#f1a340", "#998ec3"), legend= c("Contents", "Structure"), cex = 1.1)

# Scatter plot
  if(length(levels(as.factor(Structure$Region)))<=12){
    colors <- brewer.pal(length(levels(as.factor(Structure$Region))), "Dark2")
  } else {colors <- rainbow(length(levels(as.factor(Structure$Region))))}
  RegionLabels   <- sapply(levels(as.factor(Structure$Region)), 
                           FUN = function(x){paste(x, " (", sum(Contents$Region==x), ")", sep="")})
  
  pdf(file = "Figures//Figure_3_Scatter.pdf", width = 6, height = 9)
  par(mfcol=c(2,1), mar = c(4.5, 4, 1, 2) + 0.1)
  plot(Structure$DoL, Contents$DoL, ylab="Degree of loss household contents [-]", xlab="Degree of loss on structure [-]",
       col=colors[as.factor(Structure$Region)], pch = 3)
  abline(lm(Contents$DoL ~ Structure$DoL), col="red")
  rect(-1, 1, 2, 2, density = 5, col = "grey")
  rect(1, -1, 2, 1, density = 5, col = "grey")
  
  abline(a = 0, b = 1, col="lightgrey", lty=2, lwd=2)
  grid()
  legend("topleft", legend = RegionLabels, col = colors, pch = 3, bg="white")
  
  plot(Structure$Loss, Contents$Loss, ylab="Monetary loss household contents [CHF]", xlab="Monetary loss on structure [CHF]",
       col=colors[as.factor(Structure$Region)], pch = 3)
  abline(lm(Contents$Loss ~ Structure$Loss), col="red")
  abline(a = 0, b = 1, col="lightgrey", lty=2, lwd=2)
  grid()
  dev.off()
  
  