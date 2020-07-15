library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

dom_dat = fread("data/DOM.txt")
doc_dat = fread("data/DOC_molten.txt")

## drop unneeded row numbers
dom_dat[["V1"]] = NULL
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

fit = stan("code/ar1_simple.stan", data = stan_dat, iter=5000)
mcmc_dens_overlay(as.array(fit))

# compare with a latent variable model, which is a bit more correct
# some fitting complications required a bit of tuning, hence diagnostic plots and changing
# adapt_delta
fit_gp = stan("code/ar1_gp.stan", data = stan_dat, iter=10000, control=list(adapt_delta = 0.95))
parnames = c("rho", "a_gp", "sigma", "alpha", "beta")
samps = as.array(fit_gp, pars = parnames)
mcmc_trace(samps)
mcmc_dens_overlay(samps)


# compare credible intervals for the two params we care about
s1 = rstan::extract(fit, c("alpha", "beta"))
s2 = rstan::extract(fit_gp, c("alpha", "beta"))


ci_disp = function(x, digits = 3, lower = 0.05, upper = 0.95) {
	paste0(round(median(x), digits), " (",
		round(quantile(x, lower), digits), ", ",
		round(quantile(x, upper), digits), ")"
	)
}

matrix(c(sapply(s1, ci_disp, digits = 3), sapply(s2, ci_disp, digits = 3)), nrow = 2,
	dimnames = list(c("alpha", "beta"), c("AR1", "GP")))

## take home message from this
## we see similar effect size (around 0.4) with both models
## the AR1 model is a very stupid way of thinking about time dependence, the treatment
## 		basically evolves in a constant way over time, to infinity
##		As a result, the time dependence is weak, because this is a stupid way to think about
##		time and the model does best if it doesn't worry about time too much 
##		(see the gamma param from fit, it is small)
## The GP model is the opposite; it can model fukkin anything, and any slightest bit of covariance
##		in time can be captured, because it is completely non-linear and quite wiggly (the smaller the rho parameter)
##		the more wiggly the model, FYI
##		Thus, a TON of variation can be captured by the GP term, and this tends to increase uncertainty and 
##		widen the confidence intervals for the regression parameters (which are also fairly inflexible)
## Recommendations going forward
##		1. Investigate how evolution in time actually looks; we can alter the AR1 model to be smarter
##		2. Consider alternative likelihood; is y really normal? It can't be negative, so I guess not
##		3. Compare these fits to the (really dumb and obvoiously wrong) parameters from a GLM with no time dependence
##		4. If no progress on the above, consider whether a HMM is worth the trouble (intermediate in wigglyness)
##		5. If no, decide whether formal correctness but overestimated variance (the GP) is preferable
##				to stupid time model with probably underestimated variance (the AR1)
##		6. Do all of the above for each response variable of interest, the answer isn't always the same
