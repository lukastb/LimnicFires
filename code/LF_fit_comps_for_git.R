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
data <- read.table("data/pf_4comp_56samples.txt", header=T)


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


# scale data for better performance of model
comp1_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp1")
comp1_scaled <- scale(comp1_tall$value, scale=T, center=T)
comp1_tall$value <- comp1_scaled

comp2_tall <- melt(data1, id.vars = c("treat", "flume", "time"), measure.vars = "comp2")
comp2_scaled <- scale(comp2_tall$value, scale=T, center=T)
comp2_tall$value <- comp2_scaled


# make data wide
comp1_wide = dcast(comp1_tall, treat + flume ~ time, value.var="value") 
comp2_wide = dcast(comp2_tall, treat + flume ~ time, value.var="value") 


# NAs in the data -> remove rows with NAs 
comp1_wide_red <- comp1_wide[-c(which(comp1_wide$`1` %in% NA), which(comp1_wide$`8` %in% NA)),]
comp2_wide_red <- comp2_wide[-c(which(comp1_wide$`1` %in% NA), which(comp1_wide$`8` %in% NA)),]


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


# autoregressive model
fit_gp1 = stan("code/ar1_gp_comp1.stan", data = stan_dat_comp1, iter=10000, control=list(adapt_delta = 0.95)) 
fit_gp2 = stan("code/ar1_gp_comp2.stan", data = stan_dat_comp2, iter=10000, control=list(adapt_delta = 0.95))


# comp1
parnames1 = c("rho", "a_gp", "sigma", "alpha", "beta")
samps1 = as.array(fit_gp1, pars = parnames1)

# some diagnostics
mcmc_trace(samps1)
mcmc_dens_overlay(samps1)

s1 = rstan::extract(fit_gp1, c("alpha", "beta"))

median(s1$alpha) 
quantile(s1$alpha, c(0.05, 0.95))
median(s1$beta) 
quantile(s1$beta, c(0.05, 0.95))



# comp2
parnames2 = c("rho", "a_gp", "sigma", "alpha", "beta")
samps2 = as.array(fit_gp2, pars = parnames2)

# some diagnostics
mcmc_trace(samps2)
mcmc_dens_overlay(samps2)

s2 = rstan::extract(fit_gp2, c("alpha", "beta"))

median(s2$alpha) 
quantile(s2$alpha, c(0.05, 0.95))
median(s2$beta) 
quantile(s2$beta, c(0.05, 0.95))






#####plot comp1#####

# plot the total treatment effect
params = as.matrix(fit_gp1, pars = c('alpha', 'beta'))
plotting = data.frame(type = c('treatment', 'control'))
params[,2] = params[,1] + params[,2]
plotting = data.frame(type = c('control', 'treatment'), median = apply(params, 2, median), lower = apply(params, 2, quantile, 0.05), upper = apply(params, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 15, height = 10)
# create layout for 3 plots in one window
layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2, byrow = TRUE)) 

# make singel df's for contr and treat -> makes plotting with vioplot easier
comp1_tall_c <- comp1_tall[1:28,]
comp1_tall_t <- comp1_tall[29:56,]

# 1st plot
par(mar=c(5,8,2,1)) # sets the bottom, left, top and right margins 
vioplot(comp1_tall_c$value ~ comp1_tall_c$time, pchMed=19, col=rep(c("gray80"), 6), ylim=c(min(comp1_tall$value),max(comp1_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Control", side=2, line=4, cex=2)
mtext(side=2, at=max(comp1_tall$value)*1.01, "(a)", cex=2, las=1, adj=3)

# 2nd plot
par(mar=c(6,8,2,1)) # sets the bottom, left, top and right margins
vioplot(comp1_tall_t$value ~ comp1_tall_t$time, pchMed=19, col=rep(c("gray20"), 6), at=c(1:6), ylim=c(min(comp1_tall$value),max(comp1_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Treatment", side=2, line=4, cex=2)
mtext(side=2, at=max(comp1_tall$value)*1.01, "(b)", cex=2, las=1, adj=3)
mtext("Time (h)", side=1, line=4, cex=2)

# 3rd plot
par(mar=c(5,8,3,1))
plot(1:2, plotting$median, xlim=c(0,3), ylim=c(-0.2,max(plotting$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey80", "grey20"), cex=2)
segments(c(1,2), plotting$lower, y1 = plotting$upper)
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=2)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=2)
mtext("Component 1", side=2, line=1, cex=2)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting$upper)*1.05, "(c)", cex=2, las=1, adj=1)



#####plot comp2#####

# first make the new params and plotting object for fit_gp2
# plot the total treatment effect
params = as.matrix(fit_gp2, pars = c('alpha', 'beta'))
plotting = data.frame(type = c('treatment', 'control'))
params[,2] = params[,1] + params[,2]
plotting = data.frame(type = c('control', 'treatment'), median = apply(params, 2, median), lower = apply(params, 2, quantile, 0.05), upper = apply(params, 2, quantile, 0.95))



### plot raw data and model results
x11(width = 15, height = 10)
# create layout for 3 plots in one window
layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2, byrow = TRUE)) 

# make singel df's for contr and treat -> makes plotting with vioplot easier
comp2_tall_c <- comp2_tall[1:28,]
comp2_tall_t <- comp2_tall[29:56,]

# 1st plot
par(mar=c(5,8,2,1)) # sets the bottom, left, top and right margins 
vioplot(comp2_tall_c$value ~ comp2_tall_c$time, pchMed=19, col=rep(c("gray80"), 6), ylim=c(min(comp2_tall$value),max(comp2_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Control", side=2, line=4, cex=2)
mtext(side=2, at=max(comp2_tall$value)*1.01, "(a)", cex=2, las=1, adj=3)


# 2nd plot
par(mar=c(6,8,2,1)) # sets the bottom, left, top and right margins
vioplot(comp2_tall_t$value ~ comp2_tall_t$time, pchMed=19, col=rep(c("gray20"), 6), at=c(1:6), ylim=c(min(comp2_tall$value),max(comp2_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Treatment", side=2, line=4, cex=2)
mtext(side=2, at=max(comp2_tall$value)*1.01, "(b)", cex=2, las=1, adj=3)
mtext("Time (h)", side=1, line=4, cex=2)

# 3rd plot
par(mar=c(5,8,3,1))
plot(1:2, plotting$median, xlim=c(0,3), ylim=c(-0.8,max(plotting$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey80", "grey20"), cex=2)
segments(c(1,2), plotting$lower, y1 = plotting$upper)
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=2)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=2)
mtext("Component 2", side=2, line=1, cex=2)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting$upper)*1.05, "(c)", cex=2, las=1, adj=1)



