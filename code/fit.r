library(data.table)
library(tidyr)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("C:/Users/thuile/PhD/GitHub/LimnicFires")

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



#####going forward#####

# I think I'll go with the AR model as it is simpler and looking at the effect sizes it does the same job as the GP model! However, the GP model will also get some more attention down the line...

##    1. Evolution over time, i.e. is there an effect of time on DOC -> DOC slightly increases over time, more so in the treatment than in the control!
par(mar=c(5,5,2,2))
plot(as.character(doc_dat$time[doc_dat$treat=="C"]), doc_dat$value[doc_dat$treat=="C"], type="p", ylim=c(2.5,7), pch=21, bg="grey80", cex=2, xlab="Time (h)", ylab= "DOC (mg/l)", cex.lab=2, cex.axis=2)
points(as.character(doc_dat$time[doc_dat$treat=="T"]), doc_dat$value[doc_dat$treat=="T"], type="p", ylim=c(2,7), pch=21, bg="grey20", cex=2)

mod_contr <- lm(doc_dat$value[doc_dat$treat=="C"] ~ as.numeric(as.character(doc_dat$time[doc_dat$treat=="C"])))
abline(mod_contr, col="grey80", lwd=2)
mod_treat <- lm(doc_dat$value[doc_dat$treat=="T"] ~ as.numeric(as.character(doc_dat$time[doc_dat$treat=="T"])))
abline(mod_treat, col="grey20", lwd=2)

legend("topleft", legend=c("Treatment", "Control"), pch=c(21,21), pt.bg=c("grey20", "grey80"), cex=2.1, y.intersp = 0.75, bty="n")



##    2. Investigate likelihood and distribution
# check distribution of DOC -> right skewed
hist(doc_dat$value, freq = FALSE)
lines(density(doc_dat$value))

# more "normal looking", but still a bit off (too much?)
hist(log(doc_dat$value), freq = FALSE)
lines(density(log(doc_dat$value)))



##    3. compare fits with GLM 
mod <- glm(value ~ time + treat, family = gaussian, data = doc_dat)
summary(mod) 

# visualize the model
par(mar=c(5,5,2,2))
plot(as.character(doc_dat$time[doc_dat$treat=="C"]), doc_dat$value[doc_dat$treat=="C"], type="p",  pch=21, bg="grey80", cex=2, xlab="Time (h)", ylab= "DOC (mg/l)", cex.lab=2, cex.axis=2)
points(as.character(doc_dat$time[doc_dat$treat=="T"]), doc_dat$value[doc_dat$treat=="T"], type="p",  pch=21, bg="grey20", cex=2)

# maybe these polygons and lines are complete nonsense? If we look at the plot before treat and contr have different slopes!
x_line_treat = data.frame(time=seq(min(doc_dat$time), max(doc_dat$time), length.out = 1000), treat= rep("T", 1000))
y_line_treat = predict(mod, newdata = x_line_treat, type = "link", se.fit = TRUE)

polygon(c(x_line_treat$time, rev(x_line_treat$time)),(c(y_line_treat$fit + y_line_treat$se.fit, rev(y_line_treat$fit - y_line_treat$se.fit))), border=NA, col="#66666680")
lines(x_line_treat$time, y_line_treat$fit, col="black")


x_line_contr = data.frame(time=seq(min(doc_dat$time), max(doc_dat$time), length.out = 1000), treat= rep("C", 1000))
y_line_contr = predict(mod, newdata = x_line_contr, type = "link", se.fit = TRUE)

polygon(c(x_line_contr$time, rev(x_line_contr$time)),(c(y_line_contr$fit + y_line_contr$se.fit, rev(y_line_contr$fit - y_line_contr$se.fit))), border=NA, col="#66666640")
lines(x_line_contr$time, y_line_contr$fit, col="black")


# Just as a reminder: These are the params of the two rstan models:
#       AR1                    GP
# alpha "2.719 (2.559, 2.886)" "3.321 (3.094, 3.556)"
# beta  "0.433 (0.317, 0.548)" "0.405 (0.081, 0.732)"
# 
# My interpretation of this (tell me if I'm wrong):
# 1) the intercept (alpha) of the glm solution is more similar to the GP model than the AR1
# 2) the effect size (beta) is with 0.391 very similar to the ones from the AR1 and the GP model, though closer again to the GP model
# 
# What my brain does not understand:
# 1) Why is the intercept so "low" in the AR1 model in comparison to the GLM and GP model?
# 
# What my brain thinks to know:
# 1) Over all it seems that out of 3 models two are more alike (GP and GLM)
# 2) Given this one could say going with the GP model would make sense; However, the AR1 is way less complex and the overall messagse is still the same
# 3) Hence I'd stick wsith the AR1 and maybe try to tweak it a bit into getting loser to the GP regarding it's parameters
# 
# What do you think???



##    5. Generally, I'd still stick with the AR1. Mostly because I think that in this case it makes more sense to stick with the easier model. Also, blown up variacne with such few data points seems a bit wrong...

##    5a) Before moving to other varaibels I'll try to see if log transforming DOC or changing gamma makes a difference


### reducing gamma
fit_red_gamma = stan("code/ar1_simple_gamma_red.stan", data = stan_dat, iter=5000)
mcmc_dens_overlay(as.array(fit_red_gamma))
mcmc_dens_overlay(as.array(fit))


# compare credible intervals for the two params we care about; including the model with reduced gamma
s1 = rstan::extract(fit, c("alpha", "beta"))
s2 = rstan::extract(fit_gp, c("alpha", "beta"))
s3 = rstan::extract(fit_red_gamma, c("alpha", "beta"))


matrix(c(sapply(s1, ci_disp, digits = 3), sapply(s2, ci_disp, digits = 3), sapply(s3, ci_disp, digits = 3)), nrow = 2,
       dimnames = list(c("alpha", "beta"), c("AR1", "GP", "AR1_red_gamma")))
# -> reducing gamma to half does not really make a difference!!!



### log transform DOC
stan_dat_log = list(
  n_flumes = nrow(doc_wide),
  n_time = ncol(doc_wide) - 1,
  times = as.integer(colnames(doc_wide)[-1]),
  treatment = ifelse(grepl("T", doc_wide$flume), 1, 0),
  y = as.matrix(log(doc_wide[,-1]))
  
)

# reducing gamma to half does not really make a difference!!!
fit_log = stan("code/ar1_simple.stan", data = stan_dat_log, iter=5000)
mcmc_dens_overlay(as.array(fit_log))
mcmc_dens_overlay(as.array(fit))


# compare credible intervals for the two params we care about; including the model with reduced gamma and the model with log transformed DOC
s1 = rstan::extract(fit, c("alpha", "beta"))
s2 = rstan::extract(fit_gp, c("alpha", "beta"))
s3 = rstan::extract(fit_red_gamma, c("alpha", "beta"))
s4 = rstan::extract(fit_log, c("alpha", "beta"))


matrix(c(sapply(s1, ci_disp, digits = 3), sapply(s2, ci_disp, digits = 3), sapply(s3, ci_disp, digits = 3), sapply(s4, ci_disp, digits = 3)), nrow = 2,
       dimnames = list(c("alpha", "beta"), c("AR1", "GP", "AR1_red_gamma", "AR1_log")))
# -> log transforming does not seem to make a big difference; scales and absolute numbers change, but the message is the same



#### FINAL CONCLUSION: AR1 seems to do just fine
#### Possible way forward: adapt some priors: 
#### alpha ~ normal(0, 10) seems reasonanle 
#### beta ~ normal(0, 10) -> would it make sense to decrease this parameter? Is this not saying that we expect that an effect size of 10 is possible? 
#### If that's the case, we could easiliy go to 5 and that would still be a high number. What are your thoughts on this?
#### Reducing beta to 5 in the chunk of code below!!!



### reducing beta
fit_red_beta = stan("code/ar1_simple_beta_red.stan", data = stan_dat, iter=5000)
mcmc_dens_overlay(as.array(fit_red_beta))
mcmc_dens_overlay(as.array(fit))


s1 = rstan::extract(fit, c("alpha", "beta"))
s2 = rstan::extract(fit_gp, c("alpha", "beta"))
s3 = rstan::extract(fit_red_gamma, c("alpha", "beta"))
s5 = rstan::extract(fit_red_beta, c("alpha", "beta"))


matrix(c(sapply(s1, ci_disp, digits = 3), sapply(s2, ci_disp, digits = 3), sapply(s3, ci_disp, digits = 3), sapply(s5, ci_disp, digits = 3)), nrow = 2,
       dimnames = list(c("alpha", "beta"), c("AR1", "GP", "AR1_red_gamma", "AR1_red_beta")))
# -> seem that also playing around with beta is not changing a lot... Hence AR1 as you built it is the one to go with. Or do you have second thoughts?





