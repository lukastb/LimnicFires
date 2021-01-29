library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)
library(vioplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/LimnicFires") # set working directory to where the folder"LimnicFires" containing the folders "data" and "code" is

#####fit parafac comps#####
# organize data
data <- read.table("data/pf_4comp.txt", header=T)


# create treat, flume and time col
info <- as.vector(data[,1])
comps <- data[,2:5]
treat <- substr(info, start=0,stop=1)
flume <- substr(info, start=3,stop=3)
flume <- paste(flume, treat, sep = "_") # flume names should be unique
time <- substr(info, start=4,stop=4)
info2 <- cbind(treat,flume,time)

data1<-cbind.data.frame(info2,comps)
colnames(data1)[4:7]<-c("comp1","comp2","comp3","comp4")
rownames(data1) <- NULL


# melt data
comp1_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp1")
comp2_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp2")
comp3_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp3")
comp4_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp4")


# make data wide
comp1_wide = dcast(comp1_tall, treat + flume ~ time, value.var="value") 
comp2_wide = dcast(comp2_tall, treat + flume ~ time, value.var="value") 
comp3_wide = dcast(comp3_tall, treat + flume ~ time, value.var="value") 
comp4_wide = dcast(comp4_tall, treat + flume ~ time, value.var="value") 


# NAs in the data -> remove rows with NAs 
comp1_wide_red <- comp1_wide[-c(which(comp1_wide$`1` %in% NA), which(comp1_wide$`8` %in% NA)),]
comp2_wide_red <- comp2_wide[-c(which(comp2_wide$`1` %in% NA), which(comp2_wide$`8` %in% NA)),]
comp3_wide_red <- comp3_wide[-c(which(comp3_wide$`1` %in% NA), which(comp3_wide$`8` %in% NA)),]
comp4_wide_red <- comp4_wide[-c(which(comp4_wide$`1` %in% NA), which(comp4_wide$`8` %in% NA)),]


# create stan data
stan_dat_comp1 = list(
  n_flumes = nrow(comp1_wide_red),
  n_time = ncol(comp1_wide_red) - 2,
  times = as.integer(colnames(comp1_wide_red)[-c(1,2)]),
  treatment = ifelse(grepl("T", comp1_wide_red$flume), 1, 0),
  y = as.matrix(comp1_wide_red[,-c(1,2)])
)

stan_dat_comp2 = list(
  n_flumes = nrow(comp2_wide_red),
  n_time = ncol(comp2_wide_red) - 2,
  times = as.integer(colnames(comp2_wide_red)[-c(1,2)]),
  treatment = ifelse(grepl("T", comp2_wide_red$flume), 1, 0),
  y = as.matrix(comp2_wide_red[,-c(1,2)])
)

stan_dat_comp3 = list(
  n_flumes = nrow(comp3_wide_red),
  n_time = ncol(comp3_wide_red) - 2,
  times = as.integer(colnames(comp3_wide_red)[-c(1,2)]),
  treatment = ifelse(grepl("T", comp3_wide_red$flume), 1, 0),
  y = as.matrix(comp3_wide_red[,-c(1,2)])
)

stan_dat_comp4 = list(
  n_flumes = nrow(comp4_wide_red),
  n_time = ncol(comp4_wide_red) - 2,
  times = as.integer(colnames(comp4_wide_red)[-c(1,2)]),
  treatment = ifelse(grepl("T", comp4_wide_red$flume), 1, 0),
  y = as.matrix(comp4_wide_red[,-c(1,2)])
)


# autoregressive model
fit_gp1 = stan("code/ar1_gp_comp1.stan", data = stan_dat_comp1, iter=10000, control=list(adapt_delta = 0.95)) 
fit_gp2 = stan("code/ar1_gp_comp2.stan", data = stan_dat_comp2, iter=10000, control=list(adapt_delta = 0.95))
fit_gp3 = stan("code/ar1_gp_comp3.stan", data = stan_dat_comp3, iter=10000, control=list(adapt_delta = 0.95))
fit_gp4 = stan("code/ar1_gp_comp4.stan", data = stan_dat_comp4, iter=10000, control=list(adapt_delta = 0.95))


# comp1
parnames1 = c("rho", "a_gp", "sigma", "alpha", "beta")
samps1 = as.array(fit_gp1, pars = parnames1)

# some diagnostics
# mcmc_trace(samps1)
# mcmc_dens_overlay(samps1)

s1 = rstan::extract(fit_gp1, c("alpha", "beta"))

median(s1$alpha) 
quantile(s1$alpha, c(0.05, 0.95))
median(s1$beta) 
quantile(s1$beta, c(0.05, 0.95))



# comp2
parnames2 = c("rho", "a_gp", "sigma", "alpha", "beta")
samps2 = as.array(fit_gp2, pars = parnames2)

# some diagnostics
# mcmc_trace(samps2)
# mcmc_dens_overlay(samps2)

s2 = rstan::extract(fit_gp2, c("alpha", "beta"))

median(s2$alpha) 
quantile(s2$alpha, c(0.05, 0.95))
median(s2$beta) 
quantile(s2$beta, c(0.05, 0.95))



# comp3
parnames3 = c("rho", "a_gp", "sigma", "alpha", "beta")
samps3 = as.array(fit_gp3, pars = parnames3)

# some diagnostics
# mcmc_trace(samps3)
# mcmc_dens_overlay(samps3)

s3 = rstan::extract(fit_gp3, c("alpha", "beta"))

median(s3$alpha) 
quantile(s3$alpha, c(0.05, 0.95))
median(s3$beta) 
quantile(s3$beta, c(0.05, 0.95))





#####plot comp1 2 panel plot#####


# plot the total treatment effect
params1 = as.matrix(fit_gp1, pars = c('alpha', 'beta'))
plotting1 = data.frame(type = c('treatment', 'control'))
params1[,2] = params1[,1] + params1[,2]
plotting1 = data.frame(type = c('control', 'treatment'), median = apply(params1, 2, median), lower = apply(params1, 2, quantile, 0.05), upper = apply(params1, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 20, height = 10)
# create layout for one big and one small plot 
layout(matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE))

# make sorted df for boxplot with C and T of one timmepoint together
comp1_tall_s <- comp1_tall[order(comp1_tall$time),]


# 1st plot
par(mar=c(6,10,3,0)) # sets the bottom, left, top and right margins 
vioplot(comp1_tall_s$value ~ comp1_tall_s$treat + comp1_tall$time, pchMed=19, col=rep(c("gray90", "gray50"), 6), ylim=c(min(comp1_tall$value),max(comp1_tall$value)), xlab="", ylab="", yaxt="n",
        cex=2, cex.axis=2, cex.lab=2, cex.names=2, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5)) # if C and T in plot command use this: names=rep(c("C", "T"), 6)
axis(side=2, at=c(0.46, 0.48, 0.50, 0.52, 0.54), labels = T, cex.axis=3) #y-axis
axis(side=1, at=c(1.5, 4, 6.5, 9, 11.5, 14), labels = c("1", "2", "3", "4", "6", "8"), cex.axis=3, line=0, padj = 0.75) #x-axis
mtext(c("Raman Units", "Component1"), side=2, line=c(5, 7.5), cex=c(1.7, 2))
mtext(side=2, at=max(comp1_tall$value)*1.01, "(a)", cex=2, las=1, adj=2.5)
mtext("Time (h)", side=1, line=4.5, cex=2)
#mtext(c("1", "2", "3", "4", "6", "8"), side=1, line=2.5, cex=1.3, at=c(1.5, 4, 6.5, 9, 11.5, 14))

legend(0.2,0.4725, legend=c("Control", "Treatment"), pch=c(22, 22), pt.bg=c("gray90", "gray50"), bty="n", cex=3)


# 2nd plot
par(mar=c(6,7,4,0)) # sets the bottom, left, top and right margins 
plot(1:2, plotting1$median, xlim=c(0,3), ylim=c(min(plotting1$lower),max(plotting1$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey90", "grey50"), cex=3)
segments(c(1,2), plotting1$lower, y1 = plotting1$upper, lwd=c(1.3, 1.3))
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=3)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=3)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting1$upper)*1.075, "(b)", cex=2, las=1, adj=0.5)



#####plot comp2 2 panel plot#####


# plot the total treatment effect
params2 = as.matrix(fit_gp2, pars = c('alpha', 'beta'))
plotting2 = data.frame(type = c('treatment', 'control'))
params2[,2] = params2[,1] + params2[,2]
plotting2 = data.frame(type = c('control', 'treatment'), median = apply(params2, 2, median), lower = apply(params2, 2, quantile, 0.05), upper = apply(params2, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 20, height = 10)
# create layout for one big and one small plot 
layout(matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE))

# make sorted df for boxplot with C and T of one timmepoint together
comp2_tall_s <- comp2_tall[order(comp2_tall$time),]


# 1st plot
par(mar=c(6,10,3,0)) # sets the bottom, left, top and right margins 
vioplot(comp2_tall_s$value ~ comp2_tall_s$treat + comp2_tall$time, pchMed=19, col=rep(c("gray90", "gray50"), 6), ylim=c(min(comp2_tall$value),max(comp2_tall$value)), xlab="", ylab="", yaxt="n",
        cex=2, cex.axis=2, cex.lab=2, cex.names=2, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5)) # if C and T in plot command use this: names=rep(c("C", "T"), 6)
axis(side=2, at=c(0.25, 0.26, 0.27, 0.28, 0.29, 0.30), labels = T, cex.axis=3) #y-axis
axis(side=1, at=c(1.5, 4, 6.5, 9, 11.5, 14), labels = c("1", "2", "3", "4", "6", "8"), cex.axis=3, line=0, padj = 0.75) #x-axis
mtext(c("Raman Units", "Component2"), side=2, line=c(5, 7.5), cex=c(1.7, 2))
mtext(side=2, at=max(comp2_tall$value)*1.01, "(c)", cex=2, las=1, adj=2.5)
mtext("Time (h)", side=1, line=4.5, cex=2)
#mtext(c("1", "2", "3", "4", "6", "8"), side=1, line=2.5, cex=1.3, at=c(1.5, 4, 6.5, 9, 11.5, 14))

legend(0.2,0.26, legend=c("Control", "Treatment"), pch=c(22, 22), pt.bg=c("gray90", "gray50"), bty="n", cex=3)


# 2nd plot
par(mar=c(6,7,4,0)) # sets the bottom, left, top and right margins 
plot(1:2, plotting2$median, xlim=c(0,3), ylim=c(min(plotting2$lower),max(plotting2$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey90", "grey50"), cex=3)
segments(c(1,2), plotting2$lower, y1 = plotting2$upper, lwd=c(1.3, 1.3))
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=3)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=3)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting2$upper)*1.075, "(d)", cex=2, las=1, adj=0.5)



#####plot comp3 2 panel plot#####



## Try to make the figure nice! maybe put Time (h) on the very letf of the plot
## Move on on the big screen to see how it looks like


# plot the total treatment effect
params3 = as.matrix(fit_gp3, pars = c('alpha', 'beta'))
plotting3 = data.frame(type = c('treatment', 'control'))
params3[,2] = params3[,1] + params3[,2]
plotting3 = data.frame(type = c('control', 'treatment'), median = apply(params3, 2, median), lower = apply(params3, 2, quantile, 0.05), upper = apply(params3, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 20, height = 10)
# create layout for one big and one small plot 
layout(matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE))

# make sorted df for boxplot with C and T of one timmepoint together
comp3_tall_s <- comp3_tall[order(comp3_tall$time),]


# 1st plot
par(mar=c(6,10,3,0)) # sets the bottom, left, top and right margins 
vioplot(comp3_tall_s$value ~ comp3_tall_s$treat + comp3_tall$time, pchMed=19, col=rep(c("gray90", "gray50"), 6), ylim=c(min(comp3_tall$value),max(comp3_tall$value)), xlab="", ylab="", yaxt="n",
        cex=2, cex.axis=2, cex.lab=2, cex.names=2, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5)) # if C and T in plot command use this: names=rep(c("C", "T"), 6)
axis(side=2, at=c(0.19, 0.20, 0.21), labels = T, cex.axis=3) #y-axis
axis(side=1, at=c(1.5, 4, 6.5, 9, 11.5, 14), labels = c("1", "2", "3", "4", "6", "8"), cex.axis=3, line=0, padj = 0.75) #x-axis
mtext(c("Raman Units", "Component3"), side=2, line=c(5, 7.5), cex=c(1.7, 2))
mtext(side=2, at=max(comp3_tall$value)*1.01, "(a)", cex=2, las=1, adj=2.5)
mtext("Time (h)", side=1, line=4.5, cex=2)
#mtext(c("1", "2", "3", "4", "6", "8"), side=1, line=2.5, cex=1.3, at=c(1.5, 4, 6.5, 9, 11.5, 14))

legend(0.2,0.1875, legend=c("Control", "Treatment"), pch=c(22, 22), pt.bg=c("gray90", "gray50"), bty="n", cex=3)


# 2nd plot
par(mar=c(6,7,4,0)) # sets the bottom, left, top and right margins 
plot(1:2, plotting3$median, xlim=c(0,3), ylim=c(-0.05,max(plotting3$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey90", "grey50"), cex=3)
segments(c(1,2), plotting3$lower, y1 = plotting3$upper, lwd=c(1.3, 1.3))
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=3)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=3)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting3$upper)*1.1, "(b)", cex=2, las=1, adj=0.5)



#####plot comp4 2 panel plot#####



## Try to make the figure nice! maybe put Time (h) on the very letf of the plot
## Move on on the big screen to see how it looks like


# plot the total treatment effect
params4 = as.matrix(fit_gp4, pars = c('alpha', 'beta'))
plotting4 = data.frame(type = c('treatment', 'control'))
params4[,2] = params4[,1] + params4[,2]
plotting4 = data.frame(type = c('control', 'treatment'), median = apply(params4, 2, median), lower = apply(params4, 2, quantile, 0.05), upper = apply(params4, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 20, height = 10)
# create layout for one big and one small plot 
layout(matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE))

# make sorted df for boxplot with C and T of one timmepoint together
comp4_tall_s <- comp4_tall[order(comp4_tall$time),]


# 1st plot
par(mar=c(6,10,3,0)) # sets the bottom, left, top and right margins 
vioplot(comp4_tall_s$value ~ comp4_tall_s$treat + comp4_tall$time, pchMed=19, col=rep(c("gray90", "gray50"), 6), ylim=c(min(comp4_tall$value),max(comp4_tall$value)), xlab="", ylab="", yaxt="n",
        cex=2, cex.axis=2, cex.lab=2, cex.names=2, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5)) # if C and T in plot command use this: names=rep(c("C", "T"), 6)
axis(side=2, at=c(0.10, 0.12, 0.14), labels = T, cex.axis=3) #y-axis
axis(side=1, at=c(1.5, 4, 6.5, 9, 11.5, 14), labels = c("1", "2", "3", "4", "6", "8"), cex.axis=3, line=0, padj = 0.75) #x-axis
mtext(c("Raman Units", "Component4"), side=2, line=c(5, 7.5), cex=c(1.7, 2))
mtext(side=2, at=max(comp4_tall$value)*1.0225, "(c)", cex=2, las=1, adj=2.5)
mtext("Time (h)", side=1, line=4.5, cex=2)
#mtext(c("1", "2", "3", "4", "6", "8"), side=1, line=2.5, cex=1.3, at=c(1.5, 4, 6.5, 9, 11.5, 14))

legend(0.2,0.098, legend=c("Control", "Treatment"), pch=c(22, 22), pt.bg=c("gray90", "gray50"), bty="n", cex=3)


# 2nd plot
par(mar=c(6,7,4,0)) # sets the bottom, left, top and right margins 
plot(1:2, plotting4$median, xlim=c(0,3), ylim=c(-0.05,max(plotting4$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey90", "grey50"), cex=3)
segments(c(1,2), plotting4$lower, y1 = plotting4$upper, lwd=c(1.3, 1.3))
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=3)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=3)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting4$upper)*1.1, "(d)", cex=2, las=1, adj=0.5)


