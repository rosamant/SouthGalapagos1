install.packages(setdiff(c("astrochron", "quantmod", "IRISSeismic", "RColorBrewer", "viridis"), rownames(installed.packages())))

#setwd("Figures/")

library(astrochron)
library(quantmod)
library(IRISSeismic)
library(RColorBrewer)
library(viridis)

# Required color scales
mycol = viridis(8)

#### Here, we define the frequency ranges for the phase-analysis (in kyr-1).
#f_obl_low=1/0.060
#f_obl_high=1/0.025
f_ecc_low=1/0.135
f_ecc_high=1/0.085
#f_400_low=1/0.375
#f_400_high=1/0.435

# SG1 NGR and CENOGRID
SG1_agemodel2_i = read.csv("data/SG1_AgeModel2.csv")
SG1_agemodel2_i = SG1_agemodel2_i[c(1:19039),]
plot(SG1_agemodel2_i, type = "l")

CENOGRID = read.csv("data/CENOGRID_revised2025.csv")

CENO1 = data.frame(CENOGRID[,1],CENOGRID[,3])
CENO1 = CENO1[c(2458:11756),]
CENO1 = linterp(CENO1, dt = 0.001, start = 2.49)
dev.off()


SG1_GO=cbind(SG1_agemodel2_i,CENO1[,2]) #Combining NGR and d18O data in one matrix

end=length(SG1_GO[,1])

T1 = round(SG1_GO[1,1], digits = 2) # T1 is the younger limit of the FIRST analysis window (sliding window approach)
if (T1 - SG1_GO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - SG1_GO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(SG1_GO[end,1], digits = 2) # T2 is the older limit of the LAST analysis window (sliding window approach)
if (T2 - SG1_GO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - SG1_GO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2) # T1_all is an array that contains the younger limits of ALL analysis windows that will be applied to this site
T2_all = seq(T1+1.2, T2, by = 0.2) # T2_all is an array that contains the older limits of ALL analysis windows that will be applied to this site

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2) # This array contains the midpoints (in Ma) of all analyses windows. 
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1)) # These matrixes will be filled with phase and coherence output in the "Loop" that follows
#phase_obl=phase_ecc
coh_ecc=phase_ecc
#coh_obl = phase_ecc
phase_400=phase_ecc
coh_400=phase_ecc

# Loop 
for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(SG1_GO[,1]<T2 & SG1_GO[,1]>T1)
  NGR_win1=astrochron::detrend(SG1_GO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(SG1_GO[idx1,c(1,3)], genplot = F)
  
  res_NGR=length(NGR_win1[,1])/(T2-T1)
  dt = signif(1/res_NGR,1)
  time = seq(from = T1, to = T2, by = dt)
  
  NGR_win1=approx(SG1_GO[,1], SG1_GO[,2], xout = time)$y
  d18O_win1=approx(SG1_GO[,1], SG1_GO[,3], xout = time)$y
  
  NGR_win1=ts(NGR_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(NGR_win1,d18O_win1) # Making time series out of the dataset in question
  
  DF <- crossSpectrum(win1, spans=c(3,5)) # Calculating the cross spectrum
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low) # Finding the frequency window that correspond to 100-kyr eccentricity
  coh_ecc[i,2]=max(DF$coh[idx_ecc]) # Finding the frequency with maximum coherence within that frequency window
  idx_coh=which(DF$coh == coh_ecc[i,2]) 
  phase_ecc[i,2]=DF$phase[idx_coh] # Getting the phase result for the frequency with maximum coherence
  
  # Obliquity ----
  #idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  #coh_obl[i,2]=max(DF$coh[idx_obl])
  #idx_coh=which(DF$coh == coh_obl[i,2])
  #phase_obl[i,2]=DF$phase[idx_coh]
  
  # # 405 kyr eccentricity ----
  #idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  #coh_400[i,2]=max(DF$coh[idx_400])
  #idx_coh=which(DF$coh == coh_400[i,2])
  #phase_400[i,2]=DF$phase[idx_coh]
}

phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
 if (coh_ecc[i,2]<0.25) {phase_ecc[i,3]=1}
 else if (coh_ecc[i,2]<0.5) {phase_ecc[i,3]=2} 
 else {phase_ecc[i,3]=3}
}  


#phase_400=cbind(phase_400, phase_400[,1])
#for (i in 1:length(coh_400[,1])){
#  if (coh_400[i,2]<0.25) {phase_400[i,3]=1}
#  else if (coh_400[i,2]<0.5) {phase_400[i,3]=2} 
#  else {phase_400[i,3]=3}
#}  

dev.off()
plot(coh_ecc, ylim = c(0, 1), type="l", xlab = "coh_ecc")
plot(phase_ecc, ylim = c(-3.15, 3.15), xlab = "phase_ecc")

#plot(coh_obl, ylim = c(0, 1), type="l", xlab = "coh_ecc")
#plot(phase_obl, ylim = c(-3.15, 3.15), xlab = "phase_obl")

#plot(coh_400, ylim = c(0, 1), type="l", xlab = "coh_ecc")
#plot(phase_400, ylim = c(-3.15, 3.15), xlab = "phase_obl")
abline(h = 0)

xpi <- seq(-3.2, 0, length.out = 100)
xpi2 <- seq(0, 3.2, length.out = 100)

gradient_colors <- colorRampPalette(c("lightblue", "white", "gray"))(100)
rev_gradient_colors <- rev(gradient_colors)  

############################
# Make Figure 11
###########################

#setwd("/Figures/")

pdf(file = "Figure 11.pdf", width = 8.5, height = 11, paper = "a4")

png(filename = "Fig 11.png", width = 6000, height = 5000, res = 600)

par(mar=c(5,1,5,0))

layout(matrix(c(1, 2, 3, 4, 4, 4), 2, 3, byrow = TRUE), widths = c(1, 1, 1.25), heights = c(1, 0.1))

plot(CENO1$CENOGRID...3.,CENO1$CENOGRID...1., type = "l", xlim = c(0.5,5), ylim = c(21.53,2.49), xaxs = "i", yaxs = "i", col = mycol[3], lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(1, at = c(5, seq(1,5,1)), cex.axis = 1.75, col = mycol[3], col.axis = mycol[3])
mtext(expression("CENOGRID"~delta^18*"O (‰)"), side = 1, line = 2.5, cex = 1.1, col = mycol[3])

par(mar=c(5,0,5,0))
plot(SG1_agemodel2_i$GR, SG1_agemodel2_i$Time_Ma, type = "l", xlim = c(10, 70), ylim = c(21.53,2.49), xaxs = "i", yaxs = "i", col = mycol[5], lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(3, at = c(60, seq(15,60,15)), cex.axis = 1.75, col = mycol[5], col.axis = mycol[5])
mtext("NGR (gAPI)", side = 3, line = 2.5, cex = 1.1, col = mycol[5])

par(mar=c(5,0,5,5))

plot(phase_ecc[,2], phase_ecc[,1], type = "n", xlim = c(-3.2, 4), ylim = c(21.53,2.49), xaxs = "i", yaxs = "i", col = mycol[3], lwd = 2, pch = 16, cex = phase_ecc[,3]/1.5, 
     axes = F, xaxt = "n", xlab = "", ylab = "")

for (i in 1:(length(xpi) - 1)) {
  rect(xpi[i], min(CENO1[,1]), xpi[i+1], max(CENO1[,1]),
       col = gradient_colors[i], border = NA)
}

for (i in 1:(length(xpi2) - 1)) {
  rect(xpi2[i], min(CENO1[,1]), xpi2[i+1], max(CENO1[,1]),
       col = rev_gradient_colors[i], border = NA)
}

points(phase_ecc[,2], phase_ecc[,1], xlim = c(-3.2, 4), ylim = c(21.53,2.49), lwd = 2, pch = 16, cex = phase_ecc[,3]/1.5, xlab = "", ylab = "")

axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), labels = c("-pi", "-pi/2", 0, "pi/2", "pi"), cex.axis = 1.75, col = mycol[3], col.axis = mycol[3])
axis(4, at = c(2.49,5.33,7,9,11,13,15,17,19,21.53), cex.axis = 1.6)
mtext("Phase", side = 1, line = 2.5, at = 0, col = mycol[3])
abline(v=-1.57)
abline(v=0)
abline(v=1.57)
abline(h=11.7)
abline(h=12.7)

rect(3.2,2.49,4,2.58, col = rgb(255/255,239/255,175/255))
rect(3.2,2.58,4,5.33, col = rgb(255/255,255/255,153/255))
rect(3.2,5.33,4,21.53, col = rgb(255/255,255/255,0/255))
text(3.57,3.95, "Pliocene", srt = 90, cex = 2)
text(3.57,13.28, "Miocene", srt= 90, cex = 2)
mtext("Age (Ma)", side = 4, line = 2.5, cex = 1.3)
box()

mtext(text = "Anti-phase", side = 3, line = 1, at = -2.4, cex = 1.05, col = mycol[3])
mtext(text = expression(~delta^18*"O (‰) & NGR"), side = 3, line = 2, at = 0, srt = 45, cex = 1.05, col = mycol[3])
mtext(text = "in phase", side = 3, line = 0.75, at = 0, srt = 45, cex = 1.05, col = mycol[3])
mtext(text = "Anti-phase", side = 3, line = 1, at = 2.4, srt = 45, cex = 1.05, col = mycol[3])

par(mar = c(1, 0, 0, 0))  
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')

coherence_sizes <- c(1, 2, 3)  
legend_x <- seq(0.63, 0.93, length.out = 3)  
legend_y <- rep(0.9, 3) 

points(legend_x, legend_y, pch = 16, cex = coherence_sizes, col = "black")
text(legend_x + 0.01, legend_y, labels = c("< 0.25", "0.25-0.5", "≥ 0.5"), pos = 4, cex = 1.5)
text(0.8, 0.5,"Coherence", cex = 1.5)

dev.off()
