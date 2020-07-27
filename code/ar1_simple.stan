// a simple lag-1 autocorrelation model
// because time steps are not homogeneous, we compute the delta-t in hours

data {
	int<lower=1> n_flumes; // number of flumes
	int<lower=1> n_time; // number of time steps

	int<lower=0> times [n_time]; // in hours
	int<lower=0, upper=1> treatment [n_flumes]; // 0 if control, 1 if treatment, one per flume
	matrix[n_flumes, n_time] y; // response variable
}
parameters {
	real alpha; // global intercept
	real beta; // treatment effect
	real gamma; // strength of time dependence
	real<lower=0> sigma;
}
transformed parameters {
	matrix[n_flumes, n_time] Ey;
	for (f in 1:n_flumes) {
		Ey[f,1] = alpha + beta * treatment[f];
		for(t in 2:n_time) {
			int dt = times[t] - times[t-1];
			Ey[f,t] = alpha + beta * treatment[f] + gamma * y[t-1] * dt;
		}
	}

}
model {
	for (f in 1:n_flumes) {
		y[f,1] ~ normal(Ey[f,1], sigma);
		for(t in 2:n_time) {
			y[f,t] ~ normal(Ey[f,t], sigma);
		}
	}

	// these are quite uninformative priors
	alpha ~ normal(0, 10);
	beta ~ normal(0, 10);
	gamma ~ normal(0, 5); // can reduce the scale here to favor stationarity (gamma should be between -1 and 1 generally) 
	sigma ~ exponential(0.2);
}
