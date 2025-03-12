install.packages(setdiff(c("viridis"), rownames(installed.packages())))

library(viridis)

# Import datasets

# South Galapagos 1

SGalapagos1 <- read.csv("data/SouthGalapagos1.csv", header=TRUE, stringsAsFactors=FALSE)

# Required color scales
mycol = viridis(8)

#setwd("Figures/")

pdf(file = "Figure 2.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 2.png", width = 6000, height = 8000, res = 600)

par(mar=c(5,5,0,0))

plot(SGalapagos1[,2],SGalapagos1[,1], type = "l", xlim = c(-100, 200), ylim = c(3700, 370), col = mycol[5], axes = FALSE, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "", cex.lab = 1.25)
axis(1, at = c(0,50,100,150,200), cex.axis = 1.5, las = 1)
axis(2, line = -1, at = c(370,500,1000,1500,2000,2500,3000,3630), cex.axis = 1.5, las = 1)
mtext("South Galapagos-1 Depth (meters)", side = 2, line = 3, cex = 1.5)
mtext("NGR (gAPI)", side = 1, line = 2.5, at = 100, cex = 1.5)
rect(-80,370,-60,1800, col = rgb(255/255,230/255,25/255),lwd = 1.2)
rect(-80,1800,-60,2117.5, col = rgb(255/255,154/255,82/255),lwd = 1.2)
rect(-80,2117.5,-60,3094.5, col = rgb(127/255,198/255,78/255),lwd = 1.2)
rect(-80,3094.5,-60,3638, col = rgb(52/255,178/255,201/255),lwd = 1.2)
text(-70,1075,"Neogene-Quaternary",srt = 90, cex = 1.5)
text(-70,1960,"Paleogene",srt = 90, cex = 1.3)
text(-70,2600,"Cretaceous",srt = 90, cex = 1.5)
text(-70,3380,"Jurassic",srt = 90, cex = 1.5)
dev.off()
