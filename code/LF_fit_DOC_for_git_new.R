library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)
library(vioplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/LimnicFires") # set working directory to where the folder"LimnicFires" containing the folders "data" and "code" is

doc_dat = fread("data/DOC_molten.txt")

## drop unneeded row numbers
doc_dat[["V1"]] = NULL

# flume names should be unique
doc_dat[, flume := factor(paste(flume, treat, sep="_"))]

doc_wide = as(pivot_wider(doc_dat, id_cols=c("flume", "time"), names_from="time"), "data.table")
stan_dat = list(
  n_flumes = nrow(doc_wide),
  n_time = ncol(doc_wide) - 1,
  times = as.integer(colnames(doc_wide)[-1]),
  treatment = ifelse(grepl("T", doc_wide$flume), 1, 0),
  y = as.matrix(doc_wide[,-1])
  
)



# autoregressive model
fit_gp = stan("code/ar1_gp.stan", data = stan_dat, iter=10000, control=list(adapt_delta = 0.95))
parnames = c("rho", "a_gp", "sigma", "alpha", "beta")
samps = as.array(fit_gp, pars = parnames)


# some diagnostics
mcmc_trace(samps)
mcmc_dens_overlay(samps)


# extract parameters
s1 = rstan::extract(fit_gp, c("alpha", "beta"))


# get mean DOC for contr and treat as well as effect size
dat = as.matrix(fit_gp, pars = c('alpha', 'beta'))
dat1 = data.frame(type = c('treatment', 'control'))
dat[,2] = dat[,1] + dat[,2]
dat_all = data.frame(type = c('control', 'treatment'), median = apply(dat, 2, median), lower = apply(dat, 2, quantile, 0.05), upper = apply(dat, 2, quantile, 0.95))
dat_all


ci_disp = function(x, digits = 3, lower = 0.05, upper = 0.95) {
  paste0(round(median(x), digits), " (",
         round(quantile(x, lower), digits), ", ",
         round(quantile(x, upper), digits), ")"
  )
}


# show results from gp 
matrix(sapply(s1, ci_disp, digits = 3), nrow = 2,
       dimnames = list(c("alpha", "beta"), "GP")) # the beta here is the average effect size!




#####plot DOC in 2 panel plot#####


# plot the total treatment effect
params1 = as.matrix(fit_gp, pars = c('alpha', 'beta'))
plotting1 = data.frame(type = c('treatment', 'control'))
params1[,2] = params1[,1] + params1[,2]
plotting1 = data.frame(type = c('control', 'treatment'), median = apply(params1, 2, median), lower = apply(params1, 2, quantile, 0.05), upper = apply(params1, 2, quantile, 0.95))


# plot raw data and model results
x11(width = 20, height = 10)
# create layout for one big and one small plot 
layout(matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE))


# 1st plot
par(mar=c(6,10,3,0)) # sets the bottom, left, top and right margins 
vioplot(doc_dat$value ~ doc_dat$treat + doc_dat$time, pchMed=19, col=rep(c("gray90", "gray50"), 6), ylim=c(min(doc_dat$value),max(doc_dat$value)), xlab="", ylab="", yaxt="n",
        cex=2, cex.axis=2, cex.lab=2, cex.names=2, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5)) # if C and T in plot command use this: names=rep(c("C", "T"), 6)
axis(side=2, at=c(3, 4, 5, 6, 7), labels = T, cex.axis=3) #y-axis
axis(side=1, at=c(1.5, 4, 6.5, 9, 11.5, 14), labels = c("1", "2", "3", "4", "6", "8"), cex.axis=3, line=0, padj = 0.75) #x-axis
mtext(c(expression(DOC~(mg~l^-1))), side=2, line=5, cex=2)
mtext(side=2, at=max(doc_dat$value)*1.03, "(a)", cex=2, las=1, adj=2.5)
mtext("Time (h)", side=1, line=4.5, cex=2)
#mtext(c("1", "2", "3", "4", "6", "8"), side=1, line=2.5, cex=1.3, at=c(1.5, 4, 6.5, 9, 11.5, 14))

legend(0.2,7, legend=c("Control", "Treatment"), pch=c(22, 22), pt.bg=c("gray90", "gray50"), bty="n", cex=3)


# 2nd plot
par(mar=c(6,7,4,0)) # sets the bottom, left, top and right margins 
plot(1:2, plotting1$median, xlim=c(0,3), ylim=c(round(min(plotting1$lower)), round(max(plotting1$upper))), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey90", "grey50"), cex=3)
segments(c(1,2), plotting1$lower, y1 = plotting1$upper, lwd=c(1.3, 1.3))
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=3)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=3)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting1$upper)*1.025, "(b)", cex=2, las=1, adj=1)





