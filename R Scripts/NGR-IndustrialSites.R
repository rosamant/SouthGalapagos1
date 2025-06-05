install.packages(setdiff(c("RColorBrewer", "viridis"), rownames(installed.packages())))

library(RColorBrewer)
library(viridis)


Picard1 = read.csv("data/Picard 1_ageNGR.csv")
Picard1[,1] = Picard1[,1]/1000
plot(Picard1, type = 'l')


Minilya1 = read.csv("data/Minilya 1_ageNGR.csv")
Minilya1[,1] = Minilya1[,1]/1000
plot(Minilya1, type = 'l')


U1464 = read.csv("data/U1464_ageNGR.csv")
U1464[,1] = U1464[,1]/1000
plot(U1464, type = 'l')


# Import South Galapagos 1 NGR
SG1_agemodel2_i = read.csv("data/SG1_AgeModel2.csv")

############################
# Make Figure 12
###########################

setwd("Figures/")

pdf(file = "Figure 12.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 12.png", width = 5000, height = 6000, res = 600)

par(mar=c(0,5,1,5))

layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.65))

plot(SG1_agemodel2_i, xlim = c(2.49,10.0), xaxs = "i", type = "l", ylim = c(20,65), col = mycol[4], lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(60, seq(20,60,20)), cex.axis = 1.75, col = mycol[4], col.axis = mycol[4])
mtext("NGR (gAPI)", at = 40, side = 4, line = 2.5, cex = 1.0, col = mycol[4])
arrows(x0 = 6.7, y0 = 48, x1 = 5.6, y1 = 55, length = 0.1, lwd = 2)
arrows(x0 = 3.58, y0 = 64, x1 = 3.54, y1 = 43, length = 0.1, lwd = 2)
box()
text(4.85,64,"(a) South Galapagos-1 15°54'S Browse Basin", cex = 1.75, col = mycol[4])

par(mar=c(0,5,0,5))

plot(U1464, type = "l", ylim = c(15,60), xlim = c(2.49,10.0), xaxs = "i", col = mycol[3], lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(2, at = c(60, seq(20,60,20)), cex.axis = 1.75, col = mycol[3], col.axis = mycol[3])
mtext("NGR (gAPI)", at = 40, side = 2, line = 2.5, cex = 1.0, col = mycol[3])
text(4.4,57,"(b) U1464 18°4'S Roebuck Basin", cex = 1.75, col = mycol[3])
arrows(x0 = 6.52, y0 = 27, x1 = 6.45, y1 = 59, length = 0.1, lwd = 2)
arrows(x0 = 3.68, y0 = 49, x1 = 3.53, y1 = 29, length = 0.1, lwd = 2)
box()

plot(Minilya1, type = "l", ylim = c(10,50), xlim = c(2.49,10.0), xaxs = "i", col = mycol[5], lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(4, at = c(45, seq(15,45,15)), cex.axis = 1.75, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 30, side = 4, line = 2.5, cex = 1.0, col = mycol[5])
text(4.5,48,"(c) Minilya-1 18°19'S Roebuck Basin", cex = 1.75, col = mycol[5])
arrows(x0 = 6.5, y0 = 20, x1 = 6.4, y1 = 45, length = 0.1, lwd = 2)
arrows(x0 = 3.72, y0 = 48, x1 = 3.62, y1 = 21, length = 0.1, lwd = 2)
box()

par(mar=c(5,5,0,5))

plot(Picard1, type = "l", ylim = c(0,50), xlim = c(2.49,10.0), xaxs = "i", col = mycol[6], lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(1, at = c(2.58,3.0,3.6,4.0,5.0,5.33,6.0,7.24,8.0,9.0,10.00), labels = c("2.58", "3.0", "3.6", "4.0", "5.0", "5.33", "6.0", "7.24", "8.0", "9.0", "10.00"), cex.axis = 1.75)
axis(2, at = c(50, seq(10,50,20)), cex.axis = 1.75, col = mycol[6], col.axis = mycol[6])
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.25)
mtext("NGR (gAPI)", at = 30, side = 2, line = 2.5, cex = 1.0, col = mycol[6])
text(4.8,48,"(d) Picard-1 18°57'S North Carnarvon Basin", cex = 1.75, col = mycol[6])
arrows(x0 = 6.5, y0 = 27, x1 = 5.6, y1 = 47, length = 0.1, lwd = 2)
arrows(x0 = 3.62, y0 = 37, x1 = 3.55, y1 = 24, length = 0.1, lwd = 2)
box()
rect(2.49,4,2.58,-1.9, col = rgb(255/255,242/255,174/255))
rect(2.58,4,5.33,-1.9, col = rgb(255/255,255/255,153/255))
rect(5.33,4,10.0,-1.9, col = rgb(255/255,255/255,0/255))

text(2.54, 1, "P.", cex = 1.5)
text(3.95, 1, "Pliocene", cex = 1.75)
text(7.63, 1, "Miocene", cex = 1.75)

rect(2.49,4,2.58,10, col = rgb(255/255,237/255,179/255))
rect(2.58,4,3.6,10, col = rgb(255/255,255/255,191/255))
rect(3.6,4,5.33,10, col = rgb(255/255,255/255,179/255))
rect(5.33,4,7.246,10, col = rgb(255/255,255/255,115/255))
rect(7.246,4,10.0,10, col = rgb(255/255,255/255,102/255))


text(2.54, 7, "G.", cex = 1.25)
text(3.08, 7, "Piacenzian", cex = 1.5)
text(4.45, 7, "Zanclean", cex = 1.5)
text(6.28, 7, "Messinian", cex = 1.5)
text(8.63, 7, "Tortonian", cex = 1.5)

dev.off()