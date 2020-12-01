library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)
library(vioplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/LimnicFires") # set working directory to where the folder"LimnicFires" containing the folders "data" and "code" is


#####fit SUVA254#####
# organize data
dom_dat <- fread("data/DOM.txt")
dom_dat[["V1"]] = NULL # drop unneeded row numbers


# flume names should be unique
dom_dat[, flume := factor(paste(flume, treat, sep="_"))]


# extract SUVA data only and make data wide
suva_wide = dcast(dom_dat, treat + flume ~ time, value.var="SUVA254")


# create stan data
stan_dat_suva = list(
  n_flumes = nrow(suva_wide),
  n_time = ncol(suva_wide) - 2,
  times = as.integer(colnames(suva_wide)[-c(1,2)]),
  treatment = ifelse(grepl("T", suva_wide$flume), 1, 0),
  y = as.matrix(suva_wide[,-c(1,2)])
  
)


# autoregressive model
fit_gp = stan("code/ar1_gp.stan", data = stan_dat_suva, iter=10000, control=list(adapt_delta = 0.95))
parnames = c("rho", "a_gp", "sigma", "alpha", "beta")
samps = as.array(fit_gp, pars = parnames)


# some diagnostics
mcmc_trace(samps)
mcmc_dens_overlay(samps)


# show results from fit_gp
ci_disp = function(x, digits = 3, lower = 0.05, upper = 0.95) {
  paste0(round(median(x), digits), " (",
         round(quantile(x, lower), digits), ", ",
         round(quantile(x, upper), digits), ")"
  )
}

s1 = rstan::extract(fit_gp, c("alpha", "beta"))

matrix(sapply(s1, ci_disp, digits = 3), nrow = 2,
       dimnames = list(c("alpha", "beta"), "GP")) # the beta here is the average effect size!


#####plot#####

# plot raw data and model results
params = as.matrix(fit_gp, pars = c('alpha', 'beta'))
plotting = data.frame(type = c('treatment', 'control'))
params[,2] = params[,1] + params[,2]
plotting = data.frame(type = c('control', 'treatment'), median = apply(params, 2, median), lower = apply(params, 2, quantile, 0.05), upper = apply(params, 2, quantile, 0.95))


suva_tall <- melt(dom_dat, id.vars = c("treat", "flume", "time"), measure.vars = "SUVA254")


x11(width = 15, height = 10)
# create layout for 3 plots in one window
layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2, byrow = TRUE)) 

# make singel df's for contr and treat -> makes plotting with vioplot easier
suva_tall_c <- suva_tall[1:28,]
suva_tall_t <- suva_tall[29:56,]

# 1st plot
par(mar=c(5,8,2,1)) # sets the bottom, left, top and right margins 
vioplot(suva_tall_c$value ~ suva_tall_c$time, pchMed=19, col=rep(c("gray80"), 6), ylim=c(min(suva_tall$value),max(suva_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Control", side=2, line=4, cex=2)
mtext(side=2, at=max(suva_tall$value)*1.01, "(a)", cex=2, las=1, adj=3)


# 2nd plot
par(mar=c(6,8,2,1)) # sets the bottom, left, top and right margins
vioplot(suva_tall_t$value ~ suva_tall_t$time, pchMed=19, col=rep(c("gray20"), 6), at=c(1:6), ylim=c(min(suva_tall$value),max(suva_tall$value)), xlab="", ylab="", cex=2, cex.axis=2, cex.lab=2)
mtext("Treatment", side=2, line=4, cex=2)
mtext(side=2, at=max(suva_tall$value)*1.01, "(b)", cex=2, las=1, adj=3)
mtext("Time (h)", side=1, line=4, cex=2)


# 3rd plot
par(mar=c(5,8,3,1))
plot(1:2, plotting$median, xlim=c(0,3), ylim=c(3.2,max(plotting$upper)), xaxt='n', yaxt='n', xlab='', bty='n', ylab="", pch=c(21,21), bg=c("grey80", "grey20"), cex=2)
segments(c(1,2), plotting$lower, y1 = plotting$upper)
axis(1, at=c(1,2), labels=c('', ''), cex=2, cex.lab=2, cex.axis=2)
axis(2, line=-5, cex=2, cex.lab=2, cex.axis=2)
mtext("SUVA 254", side=2, line=1, cex=2)
mtext(c("Control", "Treatment"), at=c(1,2), side=1, line=2.5, cex=2)
mtext(side=2, at=max(plotting$upper)*1.001, "(c)", cex=2, las=1, adj=1)

