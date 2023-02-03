##########################################################################################
# A script to run either a Continuous or Discrete Wavelet Transform across a bacterial pop
# Author: Jared Streich
# Date Created: May 2020
# Version 0.1.0
# email: ju0@ornl.gov, if not at ORNL streich.jared@gmail.com
##########################################################################################


##########################################################################################
####################################### Libraries ########################################
##########################################################################################

install.packages("WaveletComp")
library(WaveletComp)


##########################################################################################
###################################### Start Script ######################################
##########################################################################################

##### Set Cols
cols <- colorRampPalette(c("royalblue3", "white", "firebrick3"))(100)

##### Set Directory
setwd("~/Downloads/lines/")

##### List all files to run
list.files(pattern = "*.txt")

fls.lst <- list.files(pattern = "*.txt")


##### Set Loop
i <- 2
for(i in 1:length(fls.lst)){
  sig <- read.delim(fls.lst[i], header = F)
  wvlt <- analyze.wavelet(sig)

  nm <- substr(fls.lst[i], 0, nchar(fls.lst[i])-5)
  
  png(paste(nm, "RawHeatMap.png", sep = ""), width = 1500, height = 800)
    image(1-t(as.matrix(wvlt$Power.pval)), col = cols, 
    main = paste("Wavelet Decomposition Raw Heatmap", fls.lst[i]))
  dev.off()

  png(paste(nm, "Y-Axis_PowerPval.png", sep = ""), width = 1500, height = 800)
  plot(rowSums(1-as.matrix(wvlt$Power.pval)), type = "l", col = "firebrick3", lwd = 2, 
       main = paste("Wavelet Decomposition Power Values", fls.lst[i]))
  dev.off()

  p.power.pval.which.max <- which.max(rowSums(1-as.matrix(wvlt$Power.pval)))
  p.power.pval.max <- max(rowSums(1-as.matrix(wvlt$Power.pval)))
  
  png(paste(nm, "Y-Axis_PowerAve.png", sep = ""), width = 1500, height = 800)
  plot(wvlt$Power.avg, type = "l", col = "firebrick3", lwd = 2, 
       main = paste("Wavelet Decomposition Power Average Values", fls.lst[i]))
  dev.off()
  y.which.max <- which.max(wvlt$Power.avg)
  y.max <- max(wvlt$Power.avg)
  
  png(paste(nm, ".png", sep = "_"), width = 1500, height = 800)
  wt.image(wvlt, color.key = "quantile", n.levels = 360, legend.params = list(lab = "Wavelet Power Levels, mar = 3"), 
    main = paste("Wavelet Decomposition", fls.lst[i]))
  dev.off()

  p <- c(fls.lst[i], p.power.pval.max, p.power.pval.which.max, y.which.max, y.max)
  
  if(i == 1){
    p.prev <- p
  }
  else{
    p.prev <- rbind(p.prev, p)
  }
    
}

##### Set outputs
write.table(p.prev, file = "BactWaveletProfiles_2022-03-15.txt",  quote = F, row.names = F, col.names = F, sep = " ")
