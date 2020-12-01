library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)
library(vioplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("") # set working directory to where the folder"LimnicFires" containing the folders "data" and "code" is

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




### plot the results
x11(width = 15, height = 10)
# create layout for 3 plots in one window
layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2, byrow = TRUE)) 


# make singel df's for contr and treat -> makes plotting with vioplot easier
doc_tall_c <- doc_dat[which(doc_dat$treat=="C")]
doc_tall_t <- doc_dat[which(doc_dat$treat=="T")]


# 1st plot
par(mar=c(5,8,2,1))
vioplot(doc_tall_c$value ~ doc_tall_c$time, pchMed=19, col=rep(c("gray80"), 6), ylim=c(min(doc_dat$value), max(doc_dat$value)), xlab="", ylab="", cex=2, cex.axis=2)
mtext("Control", side=2, line=4, cex=2)
mtext(side=2, at=max(doc_dat$value)*1.01, "(a)", cex=2, las=1, adj=3)


# 2nd plot
par(mar=c(6,8,2,1)) # sets the bottom, left, top and right margins
vioplot(doc_tall_t$value ~ doc_tall_t$time, pchMed=19, col=rep(c("gray20"), 6), ylim=c(min(doc_dat$value), max(doc_dat$value)), xlab="", ylab="", cex=2, cex.axis=2)
mtext("Treatment", side=2, line=4, cex=2)
mtext(side=2, at=max(doc_dat$value)*1.01, "(b)", cex=2, las=1, adj=3)
mtext("Time (h)", side=1, line=4, cex=2)


# 3rd plot
par(mar=c(5,8,3,1))
plot(1:2, plotting$median, ylim=c(3,4), xlim=c(0,3), xaxt='n', yaxt='n', xlab="", bty='n', ylab="", pch=c(21,21), bg=c("grey80", "grey20"), cex=2, cex.lab=2, cex.axis=2)
segments(c(1,2), plotting$lower, y1 = plotting$upper)
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=2)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=2)
mtext(expression(DOC~(mg~l^-1)), side=2, line=1, cex=2)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting$upper)*1.01, "(c)", cex=2, las=1, adj=1)


