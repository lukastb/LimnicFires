// model time dependence using a gaussian process
// covariance between points is modelled as a function of the time between observations
// the mean function is a function of the treatment and thus we have one value per flume
// we assume NO covariance among flumes
// thus we have a constant (among flumes) VCV matrix 

data {
	int<lower=1> n_flumes; // number of flumes
	int<lower=1> n_time; // number of time steps

	real<lower=0> times [n_time]; // in hours; this is the "distance" variable for the GP

	int<lower=0, upper=1> treatment [n_flumes]; // 0 if control, 1 if treatment, one per flume; predictor for mean
	matrix[n_flumes, n_time] y; // response variable/outcome
}
parameters {
	// GP hyperparameters
	real<lower=0> rho;
	real<lower=0> a_gp;
	real<lower=0> sigma;

	// regression parameters
	real alpha; // intercept
	real beta; // treatment effect
}
transformed parameters {
	matrix[n_flumes, n_time] mu;
	for(f in 1:n_flumes) {
		for(t in 1:n_time) {
			mu[f,t] = alpha + beta * treatment[f];
		}
	}
	
}
model {
	matrix[n_time, n_time] L_K;
	matrix[n_time, n_time] K = cov_exp_quad(times, a_gp, rho);
	real sq_sigma = square(sigma);

	// diagonal elements
	for (n in 1:n_time)
		K[n, n] = K[n, n] + sq_sigma;

	L_K = cholesky_decompose(K);


	// likelihood
	for(f in 1:n_flumes)
		y[f,] ~ multi_normal_cholesky(mu[f,], L_K);

	// GP hyperpriors
	rho ~ inv_gamma(5, 5);
	a_gp ~ std_normal();
	sigma ~ std_normal();

	// regression priors
	alpha ~ normal(0, 0.01);
	beta ~ normal(0, 1);

}
