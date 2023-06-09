require(read.table)
require(tidyverse)
require(utils)
library(foreign)
library(grid)
library(KernSmooth)
library(boot)
library(class)
library(cluster)
library(codetools)
library(compiler)
library(lattice)
library(MASS)
library(Matrix)
library(mgcv)
library(nnet)
library(parallel)
library(rpart)
library(spatial)
library(splines)
library(stats4)
library(survival)
library(tcltk)
library(tools)
library(translations)
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
install.packages(c("class", "foreign", "lattice", "nlme", "nnet", "survival"))
library(boot)
library(class)
library(cluster)
library(codetools)
library(compiler)
library(foreign)
library(grid)
library(KernSmooth)
library(lattice)
library(MASS)
library(Matrix)
library(mgcv)
library(nnet)
library(parallel)
library(rpart)
library(spatial)
library(splines)
library(stats4)
library(survival)
library(tcltk)
library(tools)
library(translations)
require(utils)

#***********NFLv1 vs. WML small ROIs****************
# you can set the file name here
fp1 <- file.path('D://16GBBACKUPUSB//BACKUP_USB_SEPTEMBER2014//May Baydoun_folder//HANDLS_PAPER51_NFLSMRI//OUTPUT//R_VOLCANO_PLOT//NFLW1_WMLV.CSV')

# read the file using the readr function read_csv() and assign it to 'bp'
NFLv1_WML<- read.csv(fp1)
head(NFLv1_WML)
str(NFLv1_WML)
table(NFLv1_WML$estimate, exclude=NULL)
# Make a basic volcano plot
with(NFLv1_WML, plot(estimate, -log10(p), pch=20, main="WML vs. NFLv1", ylim=c(0.0,7), xlim=c(-0.5,0.5)))
# Add colored points: red if padj<0.05, orange of estimate>0.20, green if both)
with(subset(NFLv1_WML, p<.05 ), points(estimate, -log10(p), pch=20, col="red"))
with(subset(NFLv1_WML, abs(estimate)>0.20), points(estimate, -log10(p), pch=20, col="orange"))
with(subset(NFLv1_WML, p<.05 & abs(estimate)>0.20), points(estimate, -log10(p), pch=20, col="green"))

  #***********NFL at v2 vs. WML small ROIs****************
  # you can set the file name here
  fp1 <- file.path('D://16GBBACKUPUSB//BACKUP_USB_SEPTEMBER2014//May Baydoun_folder//HANDLS_PAPER51_NFLSMRI//OUTPUT//R_VOLCANO_PLOT//NFLW3_WMLV.CSV')

# read the file using the readr function read_csv() and assign it to 'bp'
NFLv2_WML<- read.csv(fp1)
head(NFLv2_WML)
str(NFLv2_WML)
table(NFLv2_WML$estimate, exclude=NULL)
# Make a basic volcano plot
with(NFLv2_WML, plot(estimate, -log10(p), pch=20, main="WML vs. NFLv2", ylim=c(0.0,7), xlim=c(-0.5,0.5)))
# Add colored points: red if padj<0.05, orange of estimate>0.20, green if both)
with(subset(NFLv2_WML, p<.05 ), points(estimate, -log10(p), pch=20, col="red"))
with(subset(NFLv2_WML, abs(estimate)>0.20), points(estimate, -log10(p), pch=20, col="orange"))
with(subset(NFLv2_WML, p<.05 & abs(estimate)>0.20), points(estimate, -log10(p), pch=20, col="green"))
##