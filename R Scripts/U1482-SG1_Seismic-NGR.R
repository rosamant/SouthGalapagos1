install.packages(setdiff(c("dtw", "DescTools", "astrochron"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Required color scales
mycol = viridis(8)

# Import South Galapagos 1 and U1482 datasets

SouthGalapagos1 <- read.csv("data/SouthGalapagos1.csv", header=TRUE, stringsAsFactors=FALSE)
SouthGalapagos1=SouthGalapagos1[c(1:2073),] # Oligocene-Miocene
head(SouthGalapagos1)
plot(SouthGalapagos1, type="l", xlim = c(470, 900), ylim = c(0, 70))

U1482 <- read.csv("data/U1482.csv", header=TRUE, stringsAsFactors=FALSE)
head(U1482)
plot(U1482, type="l", xlim = c(0, 550), ylim = c(0, 40))


# Age Model

AgeModelU1482 <-read.csv("data/U1482-AgeModel_NF.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelU1482, type="l")

# Fit a 2nd order polynomial model
fitU1482 <- lm(Age..Ma. ~ poly(Depth..m., 2, raw = TRUE), data = AgeModelU1482)

# Generate predictions
depth_seqU1482 <- seq(1.08, 534.5, 0.1)
age_predU1482 <- predict(fitU1482, newdata = data.frame(Depth..m. = depth_seqU1482))
plot(depth_seqU1482, age_predU1482)

U1482_polynomialAgeModel = data.frame(depth_seqU1482,age_predU1482)

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets

SouthGalapagos1_interpolated <- linterp(SouthGalapagos1, dt = 0.1, genplot = F)
U1482_interpolated <- linterp(U1482, dt = 0.1, genplot = F)

# Scaling the data
Smean = Gmean(SouthGalapagos1_interpolated$GR)
Sstd = Gsd(SouthGalapagos1_interpolated$GR)
SouthGalapagos1_scaled = (SouthGalapagos1_interpolated$GR - Smean)/Sstd
SouthGalapagos1_rescaled = data.frame(SouthGalapagos1_interpolated$DEPT, SouthGalapagos1_scaled)

Umean = Gmean(U1482_interpolated$GR)
Ustd = Gsd(U1482_interpolated$GR)
U1482_scaled = (U1482_interpolated$GR - Umean)/Ustd
U1482_rescaled = data.frame(U1482_interpolated$DEPT, U1482_scaled)

# Resampling the data using moving window statistics
SouthGalapagos1_scaled = mwStats(SouthGalapagos1_rescaled, cols = 2, win=3, ends = T)
SouthGalapagos1_standardized = data.frame(SouthGalapagos1_scaled$Center_win, SouthGalapagos1_scaled$Average)

U1482_scaled = mwStats(U1482_rescaled, cols = 2, win=3, ends = T)
U1482_standardized = data.frame(U1482_scaled$Center_win, U1482_scaled$Average)

# Plotting the rescaled and resampled data
plot(SouthGalapagos1_standardized, type="l", xlim = c(470, 900), ylim = c(-20, 20), xlab = "South Galapagos1 Resampled Depth", ylab = "Normalized GR")
plot(U1482_standardized, type="l", xlim = c(0, 550), ylim = c(-10, 15), xlab = "U1482 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_U1482_sg1_ap1 <- dtw(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Average, U1482_standardized$U1482_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_U1482_sg1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
SouthGalapagos1_on_U1482_depth = tune(SouthGalapagos1_standardized, cbind(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_U1482_sg1_ap1$index1s], U1482_standardized$U1482_scaled.Center_win[al_U1482_sg1_ap1$index2s]), extrapolate = F)

dev.off()

plot(U1482_standardized, type = "l", ylim = c(-15,15), xlab = "U1482 Depth", ylab = ("South Galapagos-1"))
lines(SouthGalapagos1_on_U1482_depth, type = "l", col = "red")

# Plotting the figure 

#setwd("Figures/")

pdf(file = "Figure 3.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 3.png", width = 6000, height = 4000, res = 600)


plot(U1482_standardized, type = "l", xlim = c(0,540), ylim = c(-15,15), axes = F, xaxs = "i", xlab = "", ylab = "")
lines(SouthGalapagos1_on_U1482_depth, type = "l", col = mycol[5], lwd = 2)
axis(1, at = c(0,100,200,300,400,500,540), cex.axis = 1.25)
axis(2, at = c(15,seq(-15,15,15)), cex.axis = 1.25)
mtext("U1482C Depth (m)", side = 1, line = 2.5, cex = 1.3)
mtext("NGR (gAPI)", side = 2, line = 2.5, cex = 1.3)
legend(x= 350, y = 15, legend = c("U1482C","South Galapagos-1"), col = c("black",mycol[5]),lwd = 2, cex = 1.25, bty = "n")

dev.off()


# Retrieve the corresponding depth values
SG1_DTW_Depth <- SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_U1482_sg1_ap1$index1]
U1482_DTW_Depth <- U1482_standardized$U1482_scaled.Center_win[al_U1482_sg1_ap1$index2]

# Create a data frame to show the depths of South Galapagos 1 and U1482 side by side
SG1_U1482_Depth <- data.frame(SouthGalapagos1_Depth = SG1_DTW_Depth,  U1482_Depth = U1482_DTW_Depth)

