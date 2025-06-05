install.packages(setdiff(c("astrochron", "RColorBrewer", "viridis"), rownames(installed.packages())))

#setwd("Figures/")

library(RColorBrewer)
library(viridis)
library(astrochron)


# Import South Galapagos 1 NGR
SG1_agemodel2_i = read.csv("data/SG1_AgeModel2.csv")

############################
# IODP U1482C 
###########################

# Import IODP U1482 NGR 
U1482_SGR = read.delim("data/U1482C_K_U_Th.txt")
U1482_K = data.frame(U1482_SGR[,5], U1482_SGR[,6])

# U1482 Age Model

AgeModelU1482 <-read.csv("data/U1482-AgeModel_NF.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelU1482, type="l")

# Fit a 2nd order polynomial model
fitU1482 <- lm(Age..Ma. ~ poly(Depth..m., 2, raw = TRUE), data = AgeModelU1482)

# Generate predictions
depth_seqU1482 <- seq(1.08, 534.5, 0.1)
age_predU1482 <- predict(fitU1482, newdata = data.frame(Depth..m. = depth_seqU1482))
plot(depth_seqU1482, age_predU1482)

U1482_polynomialAgeModel = data.frame(depth_seqU1482,age_predU1482)

# Tuning the age model data to U1482-NGR
U1482K_on_U1482AGE = tune(U1482_K, U1482_polynomialAgeModel, extrapolate = F)
dev.off()

U1482K_on_U1482AGE = linterp(U1482K_on_U1482AGE, dt = 0.001, start = 2.49) # Interpolation to every 1 kyr data 
colnames(U1482K_on_U1482AGE) <- c("Age (Ma)","K")
dev.off()
plot(U1482K_on_U1482AGE, type="l", xlim = c(2.49, 11), ylim = c(0, 2))

############################
# IODP U1501C 
###########################

U1501_SGR = read.delim("data/U1501C_K_U_Th.txt")

U1501_K = data.frame(U1501_SGR[,5], U1501_SGR[,6])
U1501_K = U1501_K[c(1:2680),]
colnames(U1501_K) <- c("Depth","K")
U1501_K = linterp(U1501_K, dt = 0.10)
dev.off()

# U1501 Age Model

AgeModelU1501 = read.csv("data/Age-Model_PM_U1501.csv")
plot(AgeModelU1501, type="l")

U1501K_on_U1501AGE = tune(U1501_K, AgeModelU1501, extrapolate = F)
colnames(U1501K_on_U1501AGE) <- c("Age (Ma)","K")
dev.off()

U1501K_on_U1501AGE = linterp(U1501K_on_U1501AGE, dt = 0.001, start = 2.49) # Interpolation to every 1 kyr data 
dev.off()
plot(U1501K_on_U1501AGE, type="l", xlim = c(0, 23), ylim = c(0, 3))

############################
# IODP U1447A 
###########################

U1447_SGR = read.delim("data/U1447A_K_U_Th.txt")

U1447_K = data.frame(U1447_SGR[,5], U1447_SGR[,6])
colnames(U1447_K) <- c("Depth","K")
U1447_K = linterp(U1447_K, dt = 0.10)
dev.off()


# U1447A Age Model

AgeModelU1447 = read.csv("data/U1447A_NF_Age_Model.csv")
plot(AgeModelU1447, type="l")

U1447K_on_U1447AGE = tune(U1447_K, AgeModelU1447, extrapolate = F)
colnames(U1447K_on_U1447AGE) <- c("Age (Ma)","K")
dev.off()

U1447K_on_U1447AGE = linterp(U1447K_on_U1447AGE, dt = 0.001, start = 2.49) # Interpolation to every 1 kyr data 
U1447K_on_U1447AGE = mwStats(U1447K_on_U1447AGE, win = 0.015, ends = T, verbose = F, genplot = F) # Moving average window
U1447K_on_U1447AGE = data.frame(U1447K_on_U1447AGE$Center_win, U1447K_on_U1447AGE$Average)

colnames(U1447K_on_U1447AGE) <- c("Age (Ma)","K")
dev.off()

clean_data <- function(data, y_col, threshold) {
    data_clean <- data
    data_clean[[y_col]][data_clean[[y_col]] > threshold] <- NA
    return(data_clean)
}

U1447K_on_U1447AGE_clean <- clean_data(U1447K_on_U1447AGE, "K", threshold = 2.5)

############################
# Make Figure 13
###########################


#setwd("Figures/")

pdf(file = "Fig 13.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 13.png", width = 5000, height = 6000, res = 600)

par(mar=c(0,5,1,5))

layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(1,1,1,1.75))

plot(U1447K_on_U1447AGE_clean$`Age (Ma)`, U1447K_on_U1447AGE_clean$K, xlim = c(2.49,10.7), ylim = c(0.8,2.5), type = "n", xaxs = "i", axes = F, xlab = "", ylab = "")
lines(U1447K_on_U1447AGE_clean$`Age (Ma)`[1:500], U1447K_on_U1447AGE_clean$K[1:500], col = mycol[5], lwd = 2)
lines(U1447K_on_U1447AGE_clean$`Age (Ma)`[661:3811], U1447K_on_U1447AGE_clean$K[661:3811], col = mycol[5], lwd = 2)
lines(U1447K_on_U1447AGE_clean$`Age (Ma)`[4071:7040], U1447K_on_U1447AGE_clean$K[4071:7040], col = mycol[5], lwd = 2)

axis(4, at = c(2.5, seq(1.0,2.5,0.5)), cex.axis = 1.75, col = mycol[5], col.axis = mycol[5])
mtext("K (wt. %)", at = 1.75, side = 4, line = 2.5, cex = 1.2, col = mycol[5])
text(3.1,0.9,"(b) U1447A", cex = 1.75, col = mycol[5])
box()

par(mar=c(0,5,0,5))

plot(U1501K_on_U1501AGE, xlim = c(2.49,10.7), xaxs = "i", type = "l", ylim = c(0.5,2.0), col = mycol[6], lwd = 2, axes = F, xlab = "", ylab = "")
axis(2, at = c(2.0, seq(0.5,2.0,0.5)), cex.axis = 1.75, col = mycol[6], col.axis = mycol[6])
mtext("K (wt. %)", at = 1.5, side = 2, line = 2.5, cex = 1.2, col = mycol[6])
text(3.1,0.6,"(c) U1501C", cex = 1.75, col = mycol[6])
box()

plot(U1482K_on_U1482AGE$`Age (Ma)`, U1482K_on_U1482AGE$K, xlim = c(2.49,10.7), ylim = c(0,2.0), type = "n", xaxs = "i", axes = F, xlab = "", ylab = "")

lines(U1482K_on_U1482AGE$`Age (Ma)`[1:506], U1482K_on_U1482AGE$K[1:506], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[581:1010], U1482K_on_U1482AGE$K[581:1010], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[1181:2230], U1482K_on_U1482AGE$K[1181:2230], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[2360:2481], U1482K_on_U1482AGE$K[2360:2481], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[2516:4961], U1482K_on_U1482AGE$K[2516:4961], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[5030:5500], U1482K_on_U1482AGE$K[5030:5500], col = mycol[3], lwd = 2)
lines(U1482K_on_U1482AGE$`Age (Ma)`[5581:8283], U1482K_on_U1482AGE$K[5581:8283], col = mycol[3], lwd = 2)

axis(4, at = c(0,1,2), cex.axis = 1.75, col = mycol[3], col.axis = mycol[3])
mtext("K (wt. %)", side = 4, line = 2.5, cex = 1.2, col = mycol[3])
text(3.1,0.1,"(d) U1482C", cex = 1.75, col = mycol[3])
box()

par(mar=c(5,5,0,5))

plot(SG1_agemodel2_i, xlim = c(2.49,10.7), xaxs = "i", type = "l", ylim = c(10,65), col = mycol[4], lwd = 2, axes = F, xlab = "", ylab = "")
axis(1, at = c(2.58,3.0,3.6,4.0,5.0,5.33,6.0,7.24,8.0,9.0,10.0,10.70), labels = c("2.58", "3.0", "3.6", "4.0", "5.0", "5.33", "6.0", "7.24", "8.0", "9.0", "10.0","10.7"), cex.axis = 1.75)
axis(2, at = c(65, seq(20,65,15)), cex.axis = 1.75, col = mycol[4], col.axis = mycol[4])
mtext("NGR (gAPI)", side = 2, line = 2.5, cex = 1.2, col = mycol[4])
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.25)
text(3.65,23,"(e) South Galapagos-1", cex = 1.75, col = mycol[4])
box()

rect(2.49,14,2.58,8.0, col = rgb(255/255,242/255,174/255))
rect(2.58,14,5.33,8.0, col = rgb(255/255,255/255,153/255))
rect(5.33,14,10.70,8.0, col = rgb(255/255,255/255,0/255))

text(2.54, 11, "P.", cex = 1.5)
text(3.95, 11, "Pliocene", cex = 1.75)
text(7.98, 11, "Miocene", cex = 1.75)

rect(2.49,20,2.58,14, col = rgb(255/255,237/255,179/255))
rect(2.58,20,3.6,14, col = rgb(255/255,255/255,191/255))
rect(3.6,20,5.33,14, col = rgb(255/255,255/255,179/255))
rect(5.33,20,7.246,14, col = rgb(255/255,255/255,115/255))
rect(7.246,20,10.70,14, col = rgb(255/255,255/255,102/255))


text(2.54, 17, "G.", cex = 1.25)
text(3.08, 17, "Piacenzian", cex = 1.5)
text(4.45, 17, "Zanclean", cex = 1.5)
text(6.28, 17, "Messinian", cex = 1.5)
text(8.98, 17, "Tortonian", cex = 1.5)

dev.off()