install.packages(setdiff(c("DescTools", "astrochron", "biwavelet", "RColorBrewer", "viridis"), rownames(installed.packages())))

# Set working directory

#setwd("Figures/")

# Importing required libraries
library(DescTools)
library(astrochron)
library(biwavelet)
library(RColorBrewer)
library(viridis)

# Required color scales
mycol = viridis(8)
col = rev(brewer.pal(9,"YlGnBu"))
cols = c(rep("#081d58",15), col)
ccols = c(rep("#081d58",17), col)

# Import South Galapagos 1 GR data
SG1 = read.csv("data/SouthGalapagos1.csv")
SG1 = SG1[c(1:3710),]
SG1 = linterp(SG1, dt = 0.1) # Interpolation to every 10 cm data

dev.off()
plot(SG1, type = "l")

##########
# Step 1: Polynomial Age Model based on seismic and biomarker information
#############

# Importing age-depth data based on seismic horizons, DTW correlation and biomarkers of Site U1482
SG1_DepthAge <- read.csv("data/SG1-DepthAge.csv")
plot(SG1_DepthAge)

# Fit a 2nd order polynomial model
fit <- lm(Age..Ma. ~ poly(Depth..m., 2, raw = TRUE), data = SG1_DepthAge)

# Generate predictions
depth_seq <- seq(475, 1200, 0.1)
age_pred <- predict(fit, newdata = data.frame(Depth..m. = depth_seq))
plot(depth_seq, age_pred)

SG1_polynomialAgeModel = data.frame(depth_seq,age_pred)

# Tuning the NGR data to polynomial age model
SG1_AgeModel = tune(SG1, SG1_polynomialAgeModel, extrapolate = F)

dev.off()
plot(SG1_AgeModel, type = "l")

SG1_interp = linterp(SG1_AgeModel, dt = 0.001, start = 2.54) # Interpolation to every 1 kyr data
SG1_smooth = mwStats(SG1_interp, win = 0.0175, ends = T, verbose = F, genplot = F) # Moving average window

SG1_agemodel1 = data.frame(SG1_smooth[,1],SG1_smooth[,2])
SG1_agemodel1 = linterp(SG1_agemodel1, dt = 0.001, start = 2.54) # Interpolation to every 1 kyr data
dev.off()
plot(SG1_agemodel1, type = "l")

# Writing the polynomial age model to correlate with CENOGRID data in QAnalySeries
write.csv(SG1_agemodel1, file = "data/SG1_AgeModel.csv", row.names = F)

# Detrending the data 
SG1_trend = lowpass(SG1_interp, fcut = 0.25)
SG1_detrend = SG1_interp
SG1_detrend[,2] = SG1_interp[,2] - SG1_trend[,2]

dev.off()
plot(SG1_detrend, type = "l")

# Spectral analysis
SG1_mtm_age_model1 = mtm(SG1_detrend , ntap = 5, padfac = 5, demean = T, detrend = T, siglevel = 0.95, genplot = F, output = 1)
dev.off()
plot(SG1_mtm_age_model1$Frequency, SG1_mtm_age_model1$Power, xlim = c(0,70), ylim = c(1e-05,1), type = "l", log = 'y')
lines(SG1_mtm_age_model1$Frequency, SG1_mtm_age_model1$AR1_95_power, col = "red")

# 95% confidence level 
SG1_2 = linterp(SG1, dt = 0.25)
dev.off()
plot(SG1_2, type = "l")
ar_model <- ar(SG1_2[,2], aic = FALSE, order.max = 1)
rho = ar_model$ar
rho = 0.7
nsim = 1000


a1 = ar1(npts = 7251, dt = 0.1, mean = mean(SG1[,2]), sdev = sd(SG1[,2]), rho = rho, nsim = nsim)

a_agemodel1 = matrix(0, nrow = length(SG1_AgeModel[,1]), ncol = nsim)

a_agemodel1_i = matrix(0, nrow = length(SG1_agemodel1[,1]), ncol = nsim)

MC_spectra1 = matrix(0, nrow = length(SG1_mtm_age_model1[,1]), ncol = nsim)

for (i in 1:nsim) {
  a_trend = lowpass(cbind(SG1[,1],a1[,i]), fcut = 0.005, verbose = F, genplot = F)
  a_detrend = cbind(SG1[,1],a1[,i])
  a_detrend[,2] = a1[,i] - a_trend[,2]
  b = tune(cbind(SG1[,1],a_detrend[,2]), SG1_polynomialAgeModel, extrapolate = F, verbose = F, genplot = F)
  a_agemodel1[,i] = b[,2]
  d = linterp(cbind(SG1_AgeModel[,1],a_agemodel1[,i]), dt = 0.001, start = 2.54, verbose = F, genplot = F)
  a_agemodel1_i[, i] = d[1:nrow(a_agemodel1_i), 2]
  temp1 = mtm(cbind(SG1_agemodel1[,1],a_agemodel1_i[,i]), ntap = 5, padfac = 5, 
              demean = T, detrend = T, siglevel = 0.95, verbose = F, 
              genplot = F, output = 1)
  MC_spectra1[,i] = temp1$Power
  print(i)
}

quantile1_5th = apply(MC_spectra1, 1, quantile, 0.95)


plot(SG1_mtm_age_model1$Frequency, SG1_mtm_age_model1$Power, xlim = c(0,70), ylim = c(1e-05,0.5), type = "l", log = 'y')
lines(temp1$Frequency,quantile1_5th, col = "red")


############################
# Make Figure 6
###########################

#setwd("Figures/")

pdf(file = "Figure 6.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 6.png", width = 6000, height = 5000, res = 600)

par(mar=c(5,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

plot(SG1_agemodel1, type = "l", xlim = c(2.54,21.4), xaxs = "i", col = mycol[5], ylim = c(80,15), axes = F, xlab = "", ylab = "", lwd = 2)
axis(1, at = c(2, seq(2,22,1)), cex.axis = 1.25)
axis(2, at = c(60, seq(20,60,20)), cex.axis = 1.2, col = mycol[5], col.axis = mycol[5])
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.25)
mtext("NGR (gAPI)", at = 40, side = 2, line = 2.5, cex = 1.15, col = mycol[5])
text(3.5,20,"(a)", cex = 1.2)
rect(2.54,82.5,2.58,72, col = rgb(255/255,239/255,175/255))
rect(2.58,82.5,5.33,72, col = rgb(255/255,255/255,153/255))
rect(5.33,82.5,21.669,72, col = rgb(255/255,255/255,0/255))
text(4.01, 77, "Pliocene", cex = 1.6)
text(13.8, 77, "Miocene", cex = 1.6)
box()

par(mar=c(5,5,0,1))

plot(SG1_mtm_age_model1$Frequency, SG1_mtm_age_model1$Power, type = "l", lwd = 2, log = 'y', xaxs = "i", yaxs = "i", 
     xlim = c(0,70), ylim = c(1e-5,1),  cex.lab = 1.25, xlab = "Frequency (cycles/Myr)", ylab= "", cex.axis = 1.25, yaxt = "n")
lines(temp1$Frequency,quantile1_5th, col = "red")
mtext("Spectral Power",side = 2, line = 1, cex = 1.25)
text(3,0.60, "(b)", cex = 1.3)
text(7,0.30, "95%CL", cex = 1, col= "red")
text(10,0.12, "Ecc.", cex = 1.25)
text(24.39,0.070, "Obl.", cex = 1.25)
text(50, 0.025, "Pre.", cex = 1.25)
dev.off()

# Wavelet Analysis
SG1_wavelet = wt(SG1_agemodel1)

plot.biwavelet(SG1_wavelet, ncol = 26, fill.cols = ccols, plot.sig = F, ylim = c(0.015,1.0), xlab = "Age (Ma)", ylab = "Period") 
abline(h=-5.573,lty = 2)
abline(h=-4.605,lty = 2)
abline(h=-3.321,lty = 2)
abline(h=-1.304,lty = 2)


##########
# Step 2: Manually tuned SG1 record to CENOGRID d18O record in QAnalySeries
#############

# Import tie-points from the manual tuning
age_model2 = read.csv("data/QAnaly_CENO2025_TunedAge.csv")

# Tuning the age model to the new age-depth data 
SG1_agemodel2 = tune(SG1_agemodel1, age_model2, extrapolate = F, genplot = F)
SG1_agemodel2_i = linterp(SG1_agemodel2, dt = 0.001, start = 2.49) # Interpolation to every 1 kyr data 
colnames(SG1_agemodel2_i) <- c("Time_Ma", "GR")
write.csv(SG1_agemodel2_i, file = "data/SG1_AgeModel2.csv", row.names = F)

# Wavelet Analysis
SG1_age_model2_wavelet = wt(SG1_agemodel2_i)

# Spectral analysis

SG1_AgeModel2 = tune(SG1_interp, age_model2, extrapolate = F, genplot = F)
SG1_AgeModel2_i = linterp(SG1_AgeModel2, dt = 0.001, start = 2.49) # Interpolation to every 1 kyr data 

# Detrending the data
SG1_agemodel2_trend = lowpass(SG1_AgeModel2_i, fcut = 0.25)

SG1_agemodel2_detrend = SG1_AgeModel2_i
SG1_agemodel2_detrend[,2] = SG1_AgeModel2_i[,2] - SG1_agemodel2_trend[,2]

dev.off()
plot(SG1_agemodel2_detrend, type = "l")

SG1_mtm_agemodel2 = mtm(SG1_agemodel2_detrend , ntap = 5, padfac = 5, demean = T, detrend = T, siglevel = 0.95, genplot = F, output = 1)
dev.off()
plot(SG1_mtm_agemodel2$Frequency, SG1_mtm_agemodel2$Power, xlim = c(0,70), ylim = c(1e-05,1), type = "l", log = 'y')
lines(SG1_mtm_agemodel2$Frequency, SG1_mtm_agemodel2$AR1_95_power, col = "red")

# 95% confidence level 
SG1_agemodel1_2 = linterp(SG1_agemodel1, dt = 0.007)
dev.off()
plot(SG1_agemodel1_2, type = "l")
ar_model2 <- ar(SG1_agemodel1_2[,2], aic = FALSE, order.max = 1)
rho2 = ar_model2$ar
rho2 = 0.85
nsim2 = 1000


a2 = ar1(npts = 18862, dt = 0.1, mean = mean(SG1_agemodel1[,2]), sdev = sd(SG1_agemodel1[,2]), rho = rho2, nsim = nsim2)

a_agemodel2 = matrix(0, nrow = length(SG1_AgeModel2[,1]), ncol = nsim2)

a_agemodel2_i = matrix(0, nrow = length(SG1_AgeModel2_i[,1]), ncol = nsim2)

MC_spectra2 = matrix(0, nrow = length(SG1_mtm_agemodel2[,1]), ncol = nsim2)

for (i in 1:nsim2) {
  a_trend = lowpass(cbind(SG1_agemodel1[,1],a2[,i]), fcut = 0.25, verbose = F, genplot = F)
  a_detrend = cbind(SG1_agemodel1[,1],a2[,i])
  a_detrend[,2] = a2[,i] - a_trend[,2]
  b = tune(cbind(SG1_agemodel1[,1],a_detrend[,2]), age_model2, extrapolate = F, verbose = F, genplot = F)
  a_agemodel2[,i] = b[,2]
  d = linterp(cbind(SG1_AgeModel2[,1],a_agemodel2[,i]), dt = 0.001, start = 2.49, verbose = F, genplot = F)
  a_agemodel2_i[, i] = d[,2]
  temp2 = mtm(cbind(SG1_AgeModel2_i[,1],a_agemodel2_i[,i]), ntap = 5, padfac = 5, 
              demean = T, detrend = T, siglevel = 0.95, verbose = F, 
              genplot = F, output = 1)
  MC_spectra2[,i] = temp2$Power
  print(i)
}

quantile2_5th = apply(MC_spectra2, 1, quantile, 0.95)


plot(SG1_mtm_agemodel2$Frequency, SG1_mtm_agemodel2$Power, xlim = c(0,70), ylim = c(1e-05,0.5), type = "l", log = 'y')
lines(temp2$Frequency,quantile2_5th, col = "red")

############################## ----
# Import ETP and CENOGRID data for wavelet analysis
#####################-----

ETP = etp(tmin = 2490, tmax = 21528, dt = 1)
ETP[,1] = ETP[,1]/1000
ETP_La04_wavelet = wt(ETP)


CENOGRID = read.csv("data/CENOGRID_revised2025.csv")

CENO = data.frame(CENOGRID[,1],CENOGRID[,3])
CENO = CENO[c(2458:11756),]
CENO = linterp(CENO, dt = 0.001, start = 2.49)

CENO_d13C = data.frame(CENOGRID[,1],CENOGRID[,2])
CENO_d13C = CENO_d13C[c(2458:11756),]
CENO_d13C = linterp(CENO_d13C, dt = 0.001, start = 2.49)
dev.off()

CENO_wavelet = wt(CENO)

colnames(CENO) <- c("Time_Ma", "d18O")
colnames(CENO_d13C) <- c("Time_Ma", "d13C")
colnames(SG1_agemodel2_i) <- c("Time_Ma", "GR")


############################
# Make Figure 7
###########################

#setwd("Figures/")

pdf(file = "Figure 7.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 7.png", width = 6000, height = 5000, res = 600)

par(mar=c(5,5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

plot(SG1_agemodel2_i, type = "l", xlim = c(2.49,21.545), xaxs = "i", col = mycol[5], ylim = c(80,15), axes = F, xlab = "", ylab = "", lwd = 2)
axis(1, at = c(2, seq(2,22,1)), cex.axis = 1.25)
axis(2, at = c(60, seq(20,60,20)), cex.axis = 1.2, col = mycol[5], col.axis = mycol[5])
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.25)
mtext("NGR (gAPI)", at = 40, side = 2, line = 2.5, cex = 1.15, col = mycol[5])
text(3.5,20,"(a)", cex = 1.2)
rect(2.49,82.5,2.58,72, col = rgb(255/255,239/255,175/255))
rect(2.58,82.5,5.33,72, col = rgb(255/255,255/255,153/255))
rect(5.33,82.5,21.669,72, col = rgb(255/255,255/255,0/255))
text(4.01, 77, "Pliocene", cex = 1.6)
text(13.8, 77, "Miocene", cex = 1.6)
box()

par(mar=c(5,5,0,1))

plot(SG1_mtm_agemodel2$Frequency, SG1_mtm_agemodel2$Power, type = "l", lwd = 2, log = 'y', xaxs = "i", yaxs = "i", 
     xlim = c(0,70), ylim = c(1e-5,1),  cex.lab = 1.25, xlab = "Frequency (cycles/Myr)", ylab= "", cex.axis = 1.25, yaxt = "n")
lines(temp2$Frequency,quantile2_5th, col = "red")
mtext("Spectral Power",side = 2, line = 1, cex = 1.25)
text(3,0.60, "(b)", cex = 1.3)
text(7,0.30, "95%CL", cex = 1, col= "red")
text(10,0.1, "Ecc.", cex = 1.25)
text(24.39,0.060, "Obl.", cex = 1.25)
text(50, 0.015, "Pre.", cex = 1.25)
dev.off()


############################
# Make Figure 4
###########################

setwd("Figures/")

pdf(file = "Figure 4.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 4.png", width = 6000, height = 4000, res = 600)

par(mar=c(0,8,1,1))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1,1,1,1), heights=c(0.75,0.75,0.80,0.95))

plot.biwavelet(ETP_La04_wavelet, ncol = 24, fill.cols = cols, plot.sig = F, ylim = c(0.015,1.0)
               ,xaxt = "n", xlab = "", yaxt = "n", ylab = "") 
abline(h=-5.573,lty = 2)
abline(h=-4.605,lty = 2)
abline(h=-3.321,lty = 2)
abline(h=-1.304,lty = 2)

par(mar=c(0,8,0,1))

plot.biwavelet(CENO_wavelet, ncol = 24, fill.cols = cols, plot.sig = F, ylim = c(0.015,1.0)
               ,xaxt = "n", xlab = "", yaxt = "n", ylab = "") 
abline(h=-5.573,lty = 2)
abline(h=-4.605,lty = 2)
abline(h=-3.321,lty = 2)
abline(h=-1.304,lty = 2)

plot.biwavelet(SG1_age_model2_wavelet, ncol = 26, fill.cols = ccols, plot.sig = F, ylim = c(0.015,1.0)
               ,xaxt = "n", xlab = "", yaxt = "n", ylab = "")  
abline(h=-5.573,lty = 2)
abline(h=-4.605,lty = 2)
abline(h=-3.321,lty = 2)
abline(h=-1.304,lty = 2)


par(mar=c(5,8,0,1))

plot(SG1_agemodel2_i, type = "l", xlim = c(2.49,21.528), xaxs = "i", col = mycol[5], ylim = c(80,15), axes = F, xlab = "", ylab = "", lwd = 2)
axis(1, at = c(2, seq(2,22,1)), cex.axis = 1.5)
axis(2, at = c(60, seq(20,60,20)), cex.axis = 1.5, col = mycol[5], col.axis = mycol[5])
mtext("Age (Ma)", side = 1, line = 2.5)
mtext("NGR (gAPI)", at = 40, side = 2, line = 2.5, cex = 0.75, col = mycol[5])
box()

rect(2.49,82.5,2.58,72, col = rgb(255/255,239/255,175/255))
rect(2.58,82.5,5.33,72, col = rgb(255/255,255/255,153/255))
rect(5.33,82.5,21.528,72, col = rgb(255/255,255/255,0/255))
text(3.91, 77, "Pliocene", cex = 1.5)
text(13.55, 77, "Miocene", cex = 1.5)

dev.off()


############################
# Make Figure S3
###########################

age2_depth_age <- tune(age_model2, cbind(SG1_polynomialAgeModel[,2], SG1_polynomialAgeModel[,1]), extrapolate = T)

dev.off()
colnames(age2_depth_age) <- c("depth", "age (Ma)")

age1_depth_age = SG1_DepthAge

Events <- c("T Discoaster pentaradiatus","T Discoaster surculus", "T Discoaster tamalis", "T Sphenolithus spp.", "T Reticulofenestra pseudoumbilicus",
            "T Amaurolithus tricorniculatus (/delicatus)", "B Discoaster brouweri", "T Ceratolithus armathus", "B Ceratolithus cristatus",
            "T Orthostylus rugosus", "B Ceratolithus armatus", "T Discoaster quinqueramus", "T Nicklithus amplificus",
            "B Nicklithus amplificus", "B Amaurolithus spp.", "B Discoaster surculus", "B Discoaster berggrenii", 
            "T Catinaster coalitus", "SH", "SH", "SH")

age1_depth_age <- cbind(age1_depth_age, Events)

B_events <- grepl("^B", age1_depth_age$Events)  # Events starting with 'B'
T_events <- grepl("^T", age1_depth_age$Events)  # Events starting with 'T'
SH_events <- age1_depth_age$Events == "SH"      # Events equal to 'SH'

SouthGalapagos1 = SG1_polynomialAgeModel
SouthGalapagos1[,2]=SouthGalapagos1[,1]  
SouthGalapagos1_DepthAge=tune(SouthGalapagos1, age2_depth_age, extrapolate = F, genplot = F)
colnames(SouthGalapagos1_DepthAge) <- c("Time_Ma", "Depth")
write.csv(SouthGalapagos1_DepthAge, file = "data/SouthGalapagos1_DepthAge.csv", row.names = F)

#setwd("Figures/")

pdf(file = "Fig S3.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig S3.png", width = 6000, height = 4000, res = 600)

par(mar=c(5,4,1,0))
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1,1.05), heights=c(1,1))

plot(age1_depth_age[,2], age1_depth_age[,1], type = "n", ylim = c(1300,400), xlim = c(0, 23.03), xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
points(age1_depth_age$Age[T_events], age1_depth_age$Depth[T_events], pch = 25, cex = 2, bg = mycol[8])                       # Color of points
points(age1_depth_age$Age[B_events], age1_depth_age$Depth[B_events], pch = 24, cex = 2, bg = mycol[7])                       # Color of points
points(age1_depth_age$Age[SH_events], age1_depth_age$Depth[SH_events], pch = 15, cex = 2, col = mycol[5])                       # Color of points
lines(SG1_polynomialAgeModel$age_pred, SG1_polynomialAgeModel$depth_seq, col = mycol[3], lwd = 2)

axis(1, at = c(0,2.58,5.33,11.63,15.98,23.03), cex.axis = 0.9)
axis(2, at = c(400, seq(400,1200,200)), cex.axis = 0.9)
box()

rect(0,1300,2.58,1250, col = rgb(255/255,242/255,174/255))
rect(2.58,1300,5.33,1250, col = rgb(255/255,255/255,153/255))
rect(5.33,1300,23.03,1250, col = rgb(255/255,255/255,0/255))

text(1.4, 1275, "Pleist.", cex = 1.05)
text(3.95, 1275, "Plio.", cex = 1.05)
text(14.25, 1275, "Miocene", cex = 1.05)

mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.1)
mtext("Depth (m)", side = 2, line = 2.5, cex = 1.1)
legend(x=8.1,y=450, legend = c("Polynomial Age Model"), col = c(mycol[3]), lwd = 2, cex = 0.9, bty = "n")
legend(x=9,y=480, legend = c("   U1482C Base Biomarker", "   U1482C Top Biomarker", "   Seismic Horizons"), col = c(mycol[7], mycol[8], mycol[5]), 
       pt.bg = c(mycol[7], mycol[8], mycol[5], mycol[1]), pch = c(24,25,15,19), pt.cex = 1.75, cex = 0.9, bty = "n")


par(mar=c(5,4,1,1))

plot(age1_depth_age[,2], age1_depth_age[,1], type = "n", ylim = c(1300,400), xlim = c(0, 23.03), xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
points(age2_depth_age[,2], age2_depth_age[,1], col = mycol[1], pch = 19, cex = 1.5)
lines(SG1_polynomialAgeModel$age_pred, SG1_polynomialAgeModel$depth_seq, col = mycol[3], lwd = 2)

axis(1, at = c(0,2.58,5.33,11.63,15.98,23.03), cex.axis = 0.9)
axis(2, at = c(400, seq(400,1200,200)), cex.axis = 0.9)
box()

rect(0,1300,2.58,1250, col = rgb(255/255,242/255,174/255))
rect(2.58,1300,5.33,1250, col = rgb(255/255,255/255,153/255))
rect(5.33,1300,23.03,1250, col = rgb(255/255,255/255,0/255))

text(1.4, 1275, "Pleist.", cex = 1.05)
text(3.95, 1275, "Plio.", cex = 1.05)
text(14.25, 1275, "Miocene", cex = 1.05)


mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.1)
mtext("Depth (m)", side = 2, line = 2.5, cex = 1.1)
legend(x=5.3,y=450, legend = c("Polynomial Age Model"), col = c(mycol[3]), lwd = 2, cex = 0.9, bty = "n")
legend(x=6,y=480, legend = c("   CENOGRID Tuned Age Model (n=31)"), col = c(mycol[1]), 
       pt.bg = c(mycol[1]), pch = c(19), pt.cex = 1.75, cex = 0.9, bty = "n")

dev.off()


############################
# Make Figure 3
###########################

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 3.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Figure 3.png", width = 6000, height = 6000, res = 600)

par(mar=c(5,4,2,2))

plot(age1_depth_age[,2], age1_depth_age[,1], type = "n", ylim = c(1300,400), xlim = c(0, 23.03), xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
points(age2_depth_age[,2], age2_depth_age[,1], col = mycol[1], pch = 19, cex = 1.15)
lines(SouthGalapagos1_DepthAge$Time_Ma, SouthGalapagos1_DepthAge$Depth, col = mycol[1], lwd = 2)
points(age1_depth_age$Age[T_events], age1_depth_age$Depth[T_events], pch = 25, cex = 1.15, bg = mycol[8])                       # Color of points
points(age1_depth_age$Age[B_events], age1_depth_age$Depth[B_events], pch = 24, cex = 1.15, bg = mycol[7])                       # Color of points
points(age1_depth_age$Age[SH_events], age1_depth_age$Depth[SH_events], pch = 15, cex = 1.15, col = mycol[5])                       # Color of points

axis(1, at = c(0,2.58,5.33,11.63,15.98,23.03), cex.axis = 1.25)
axis(2, at = c(400, seq(400,1200,200)), cex.axis = 1.25)
box()

rect(0,1300,2.58,1250, col = rgb(255/255,242/255,174/255))
rect(2.58,1300,5.33,1250, col = rgb(255/255,255/255,153/255))
rect(5.33,1300,23.03,1250, col = rgb(255/255,255/255,0/255))

text(1.4, 1275, "Pleist.", cex = 1.25)
text(3.95, 1275, "Plio.", cex = 1.25)
text(14.25, 1275, "Miocene", cex = 1.25)

mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.2)
mtext("Depth (m)", side = 2, line = 2.5, cex = 1.2)
legend(x=12,y=480, legend = c("   U1482C Base Bioevent", "   U1482C Top Bioevent", "   Seismic Horizons", "   CENOGRID Tuned Age Model (n=31)"), col = c(mycol[7], mycol[8], mycol[5], mycol[1]), 
       pt.bg = c(mycol[7], mycol[8], mycol[5], mycol[1]), pch = c(24,25,15,19), pt.cex = 1.75, cex = 1.2, bty = "n")
legend(x=11.5,y=580, legend = c("Tuned Age Model"), col = c(mycol[1]), lwd = 2, cex = 1.2, bty = "n")

dev.off()


############################
# Make Figure 5 - d18O & d13C - Burdigalian & Aquitanian
###########################

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 5.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 5.png", width = 6000, height = 4000, res = 600)


par(mar=c(5,5,1,5))

plot(CENO$Time_Ma, CENO$d18O, xlim = c(15.5,23.03), xaxs = "i", type = "l", ylim = c(5.5,1.0), axes = F, xlab = "", ylab = "")

rect(21,5,21.3,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(20.15,5,20.25,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(19.35,5,19.45,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(17.35,5,17.65,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(15.95,5,16.1,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(15.63,5,15.5,0.8, col = rgb(254/255,217/255,154/255), border = NA, )

lines(CENO$Time_Ma, CENO$d18O, lwd = 2)
axis(1, at = c(15.5,17.0,18.0,19.0,20.0,21.0,22.0,23.03), cex.axis = 1.25)
axis(2, at = c(5.0,4.0,3.0,2.0,1.0), cex.axis = 1.25)
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.3)
mtext(expression("Benthic"~delta^18*"O (‰)"), at = 3, side = 2, line = 2.5, cex = 1.3)

rect(15.5,5.35,15.98,5.68, col = rgb(255/255,255/255,0/255))
rect(15.98,5.35,23.03,5.68, col = rgb(255/255,255/255,0/255))
text(19.5, 5.515, "Early Miocene", cex = 1.25)
rect(15.5,5,15.98,5.35, col = rgb(255/255,255/255,77/255))
rect(15.98,5.0,20.44,5.35, col = rgb(255/255,255/255,65/255))
rect(20.44,5.0,23.03,5.35, col = rgb(255/255,255/255,51/255))
text(15.75, 5.175, "La.", cex = 1.25)
text(18.23, 5.175, "Burdigalian", cex = 1.25)
text(21.72, 5.175, "Aquitanian", cex = 1.25)

par(new = T)
plot(SG1_agemodel2_i, type = "l", xlim = c(15.5,23.03), xaxs = "i", col = mycol[5], ylim = c(80,15), lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(65, seq(20,65,15)), cex.axis = 1.25, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 40, side = 4, line = 2.5, cex = 1.3, col = mycol[5])
box()


text(21.15, 67, "Mi-1a", srt = 90)
text(20.2, 66.5, "Mi-1aa", srt = 90)
text(19.4, 66.5, "Mi-1ab", srt = 90)
text(17.5, 67, "Mi-1b", srt = 90)
text(16.02, 68, "Mi2", srt = 90)
text(15.56, 63, "Peak MMCO", srt = 90)
dev.off()


###########################
# Make Figure 6 - d18O & d13C - Serravallian & Langhian
###########################


setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 6.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 6.png", width = 6000, height = 4000, res = 600)

par(mar=c(5,5,1,5))

plot(CENO$Time_Ma, CENO$d18O, xlim = c(11.0,16.5), xaxs = "i", type = "l", ylim = c(5.5,1.0), axes = F, xlab = "", ylab = "")

rect(16.1,5,15.95,0.8, col = rgb(217/255,241/255,247/255), border = NA, )
rect(15.63,5,15.5,0.8, col = rgb(254/255,217/255,154/255), border = NA, )
rect(14.95,5,15.05,0.8, col = rgb(217/255,241/255,247/255), border = NA, )
rect(14.75,5,14.85,0.8, col = rgb(217/255,241/255,247/255), border = NA, )
rect(13.75,5,13.85,0.8, col = rgb(217/255,241/255,247/255), border = NA, )
rect(12.75,5,13,0.8, col = rgb(217/255,241/255,247/255), border = NA, )

lines(CENO$Time_Ma, CENO$d18O, lwd = 2)
axis(1, at = c(11.0,12.0,13.0,14.0,15.0,16.0,16.5), cex.axis = 1.25)
axis(2, at = c(5.0,4.0,3.0,2.0,1.0), cex.axis = 1.25)
mtext(expression("Benthic"~delta^18*"O (‰)"), at = 3, side = 2, line = 2.5, cex = 1.3)
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.3)
text(12.32, 1.5, "Anti-phase", cex = 1.25, col = "red")
text(13.05, 1.5, "In-phase", cex = 1.25)
text(12.6, 4.1, "Phase shift", cex = 1.25, srt = 90)

text(12.0, 3.9, "see Fig. 9", cex = 1.05, col = "red")
text(12.0, 4.15, "for opposite", cex = 1.05, col = "red")
text(12.0, 4.4, "phase-relationship", cex = 1.05, col = "red")

rect(11.0,5.35,11.63,5.68, col = rgb(255/255,255/255,0/255))
rect(11.63,5.35,15.98,5.68, col = rgb(255/255,255/255,0/255))
rect(15.98,5.35,16.5,5.68, col = rgb(255/255,255/255,0/255))
text(13.80, 5.51, "Middle Miocene", cex = 1.25)
rect(11.0,5,11.63,5.35, col = rgb(255/255,255/255,102/255))
rect(11.63,5,13.82,5.35, col = rgb(255/255,255/255,89/255))
rect(13.82,5,15.98,5.35, col = rgb(255/255,255/255,77/255))
rect(15.98,5,16.5,5.35, col = rgb(255/255,255/255,65/255))
text(16.24, 5.175, "Bu.", cex = 1.25)
text(14.9, 5.175, "Langhian", cex = 1.25)
text(12.72, 5.175, "Serravallian", cex = 1.25)
text(11.31, 5.175, "To.", cex = 1.25)

par(new = T)
plot(SG1_agemodel2_i, type = "l", xlim = c(11.0,16.5), xaxs = "i", col = mycol[5], ylim = c(80,15), lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(65, seq(20,65,15)), cex.axis = 1.25, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 40, side = 4, line = 2.5, cex = 1.3, col = mycol[5])
rect(12.7,72.7,12.7,9, col = rgb(217/255,241/255,247/255))
box()

arrows(x0 = 14.1, y0 = 15, x1 = 14.8, y1 = 15, length = 0.1, col = "black", lwd = 2)
arrows(x0 = 13.5, y0 = 15, x1 = 12.8, y1 = 15, length = 0.1, col = "black", lwd = 2)
text(13.8, 15, "MMCT")

text(16.02, 67, "Mi2", srt = 90)
text(15.56, 63, "Peak MMCO", srt = 90)
text(15.0, 66.5, "Mi2a", srt = 90)
text(14.8, 66.5, "Mi3a", srt = 90)
text(13.8, 67, "Mi3", srt = 90)
text(12.87, 67, "Mi4", srt = 90)

dev.off()


############################
# Make Figure 9 - d18O & d13C - Tortonian & Serravallian
###########################

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 9.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 9.png", width = 6000, height = 4000, res = 600)

par(mar=c(5,5,1,5))

plot(CENO$Time_Ma, CENO$d18O, xlim = c(7.246,13.82), xaxs = "i", type = "l", ylim = c(4.5,1.5), axes = F, xlab = "", ylab = "")

rect(12.75,5,13,0.85, col = rgb(217/255,241/255,247/255), border = NA, )
rect(10.7,5,10.83,-2, col = rgb(254/255,217/255,154/255), border = NA, )
rect(10.65,5,10.55,-2, col = rgb(217/255,241/255,247/255), border = NA, )
rect(9.85,5,9.75,-2, col = rgb(217/255,241/255,247/255), border = NA, )

lines(CENO$Time_Ma, CENO$d18O, lwd = 2)
axis(1, at = c(7.25,8.0,9.0,10.0,11.0,12.0,13.0,13.82), cex.axis = 1.25)
axis(2, at = c(3.5,2.5,1.5), cex.axis = 1.25)
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.3)
mtext(expression("Benthic"~delta^18*"O (‰)"), at = 3, side = 2, line = 2.5, cex = 1.3)
text(12.2, 1.7, "Anti-phase", cex = 1.25, col = "red")
text(13.1, 1.7, "In-phase", cex = 1.25)
text(12.6, 3.7, "Phase shift", cex = 1.25, srt = 90)

text(13.37, 3.4, "see Fig. 6", cex = 1.05, col = "red")
text(13.37, 3.55, "for opposite", cex = 1.05, col = "red")
text(13.37, 3.7, "phase", cex = 1.05, col = "red")
text(13.37, 3.85, "relationship", cex = 1.05, col = "red")

rect(7.246,4.36,11.63,4.62, col = rgb(255/255,255/255,0/255))
rect(11.63,4.36,13.82,4.62, col = rgb(255/255,255/255,0/255))
text(12.72, 4.49, "Middle Miocene", cex = 1.25)
text(9.43, 4.49, "Late Miocene", cex = 1.25)
rect(7.246,4.1,11.63,4.36, col = rgb(255/255,255/255,102/255))
rect(11.63,4.1,13.82,4.36, col = rgb(255/255,255/255,89/255))
text(9.44, 4.23, "Tortonian", cex = 1.25)
text(12.71, 4.23, "Serravallian", cex = 1.25)

par(new = T)
plot(SG1_agemodel2_i, type = "l", xlim = c(7.246,13.82), xaxs = "i", col = mycol[5], ylim = c(0,65), lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(65, seq(20,65,15)), cex.axis = 1.25, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 45, side = 4, line = 2.5, cex = 1.3, col = mycol[5])
rect(12.7,70,12.7,9, col = rgb(217/255,241/255,247/255))
box()

#arrows(x0 = 13.85, y0 = 60, x1 = 13.85, y1 = 15, length = 0.1, col = mycol[5], lwd = 2)
#arrows(x0 = 13.85, y0 = 15, x1 = 13.85, y1 = 60, length = 0.1, col = mycol[5], lwd = 2)
#text(13.85, 63, "High", col = mycol[5], cex = 1.2)
#text(13.85, 12, "Low", col = mycol[5], cex = 1.2)

text(9.8, 12, "Mi7", srt = 90)
text(10.6, 12, "Mi6", srt = 90)
text(12.87, 12, "Mi4", srt = 90)

dev.off()


############################
# Make Figure 10 - d18O & d13C - - Messinian & Tortonian
###########################

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 10.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 10.png", width = 6000, height = 4000, res = 600)


par(mar=c(5,5,1,5))

plot(CENO$Time_Ma, CENO$d18O, xlim = c(4.5,8), xaxs = "i", type = "l", ylim = c(4.5,2.0), axes = F, xlab = "", ylab = "")

rect(5.73,5,4.9,-3, col = rgb(254/255,230/255,170/255), border = NA, )
rect(5.81,5,5.73,-3, col = rgb(217/255,241/255,247/255), border = NA, )
rect(5.23,5,5.17,-3, col = rgb(217/255,241/255,247/255), border = NA, )
rect(5.03,5,4.97,-3, col = rgb(217/255,241/255,247/255), border = NA, )
rect(4.78,5.2,4.9,1.87, col = rgb(217/255,241/255,247/255), border = NA, )

lines(CENO$Time_Ma, CENO$d18O, lwd = 2)
axis(1, at = c(4.5,5.0,6.0,7.0,8.0), cex.axis = 1.25)
axis(2, at = c(4.0,3.0,2.0), cex.axis = 1.25)
mtext(expression("Benthic"~delta^18*"O (‰)"), at = 3, side = 2, line = 2.5, cex = 1.3)
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.3)

rect(4.5,4.4,5.33,4.6, col = rgb(255/255,255/255,153/255))
text(4.91, 4.5, "Pliocene", cex = 1.25)
rect(5.33,4.4,8.0,4.6, col = rgb(255/255,255/255,0/255))
text(6.66, 4.5, "Late Miocene", cex = 1.25)

rect(4.5,4.2,5.33,4.4, col = rgb(255/255,255/255,179/255))
rect(5.33,4.2,7.246,4.4, col = rgb(255/255,255/255,115/255))
rect(7.246,4.2,8.0,4.4, col = rgb(255/255,255/255,102/255))
text(4.9, 4.3, "Zanclean", cex = 1.25)
text(6.28, 4.3, "Messinian", cex = 1.25)
text(7.62, 4.3, "Tortonian", cex = 1.25)

par(new = T)
plot(SG1_agemodel2_i, type = "l", xlim = c(4.5,8), xaxs = "i", col = mycol[5], ylim = c(0,65), lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(65, seq(20,65,15)), cex.axis = 1.25, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 45, side = 4, line = 2.5, cex = 1.3, col = mycol[5])
box()

#arrows(x0 = 8.8, y0 = 60, x1 = 8.8, y1 = 15, length = 0.1, col = mycol[5], lwd = 2)
#arrows(x0 = 8.8, y0 = 15, x1 = 8.8, y1 = 60, length = 0.1, col = mycol[5], lwd = 2)
#text(8.8, 63, "High", col = mycol[5], cex = 1.2)
#text(8.8, 12, "Low", col = mycol[5], cex = 1.2)

text(5.77, 12, "TG22", srt = 90)
text(5.2, 11.5, "TG4", srt = 90)
text(5.0, 11, "T2", srt = 90)
text(4.8, 12.5, "Si2", srt = 90)
text(4.87, 12.5, "Si4", srt = 90)

dev.off()


############################
# Make Figure 11 - d18O & d13C - Piacenzian & Zanclean
###########################

setwd("C:/Users/Rohit/OneDrive - Universität Münster/Paper Drafts/Paper 2/Figures/")

pdf(file = "Figure 11.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 11.png", width = 6000, height = 4000, res = 600)

par(mar=c(5,5,1,5))

plot(CENO$Time_Ma, CENO$d18O, xlim = c(2.49,5.333), xaxs = "i", type = "l", ylim = c(5.5,2.0), axes = F, xlab = "", ylab = "")

rect(3.26,5.2,2.78,-3, col = rgb(254/255,230/255,170/255), border = NA, )
rect(2.7,5.2,2.78,1.87, col = rgb(217/255,241/255,247/255), border = NA, )
rect(3.26,5.2,3.34,1.87, col = rgb(217/255,241/255,247/255), border = NA, )
rect(4.33,5.2,4.41,1.87, col = rgb(217/255,241/255,247/255), border = NA, )
rect(4.78,5.2,4.9,1.87, col = rgb(217/255,241/255,247/255), border = NA, )

rect(5.73,5.2,4.9,-3, col = rgb(254/255,230/255,170/255), border = NA, )
rect(5.23,5.2,5.17,-3, col = rgb(217/255,241/255,247/255), border = NA, )
rect(5.03,5.2,4.97,-3, col = rgb(217/255,241/255,247/255), border = NA, )

lines(CENO$Time_Ma, CENO$d18O, lwd = 2)
axis(1, at = c(2.49,3.0,3.5,4.0,4.5,5.0,5.33), cex.axis = 1.25)
axis(2, at = c(5.0,4.0,3.0,2.0), cex.axis = 1.25)
mtext("Age (Ma)", side = 1, line = 2.5, cex = 1.3)
mtext(expression("Benthic"~delta^18*"O (‰)"), at = 3.5, side = 2, line = 2.5, cex = 1.3)

rect(2.49,5.35,2.58,5.64, col = rgb(255/255,239/255,175/255))
text(2.54, 5.5, "Pl.", cex = 1.25)
rect(2.58,5.35,5.33,5.64, col = rgb(255/255,255/255,153/255))
text(3.85, 5.5, "Pliocene", cex = 1.25)

rect(2.49,5.05,2.58,5.35, col = rgb(255/255,237/255,179/255))
rect(2.58,5.05,3.6,5.35, col = rgb(255/255,255/255,191/255))
rect(3.6,5.05,5.33,5.35, col = rgb(255/255,255/255,179/255))
text(2.54, 5.2, "G.", cex = 1.25)
text(3.08, 5.2, "Piacenzian", cex = 1.25)
text(4.45, 5.2, "Zanclean", cex = 1.25)

par(new = T)
plot(SG1_agemodel2_i, type = "l", xlim = c(2.49,5.333), xaxs = "i", col = mycol[5], ylim = c(0,65), lwd = 2, axes = F, xlab = "", ylab = "")
axis(4, at = c(65, seq(20,65,15)), cex.axis = 1.25, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", at = 45, side = 4, line = 2.5, cex = 1.3, col = mycol[5])
box()

#arrows(x0 = 5.1, y0 = 60, x1 = 5.1, y1 = 20, length = 0.1, col = mycol[5], lwd = 2)
#arrows(x0 = 5.1, y0 = 20, x1 = 5.1, y1 = 60, length = 0.1, col = mycol[5], lwd = 2)
#text(5.1, 63, "High", col = mycol[5], cex = 1.2)
#text(5.1, 17, "Low", col = mycol[5], cex = 1.2)

text(3.02, 13, "mPWP")
text(2.74, 11.5, "G6", srt = 90)
text(3.3, 11.5, "M2", srt = 90)
text(4.371, 12.5, "CN4", srt = 90)
text(4.8, 12.5, "Si2", srt = 90)
text(4.87, 12.5, "Si4", srt = 90)
text(5.2, 12.5, "TG4", srt = 90)
text(5.0, 11.5, "T2", srt = 90)

dev.off()