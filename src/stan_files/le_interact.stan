functions {
#include /chunks/le_functions.stan
}
data {
#include /chunks/le_data.stan
}
parameters {

  matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age;
  vector[nk_incrate_time] coef_incrate_time_young;
  vector[nk_natmx_time] coef_natmx_time;
  vector[nk_natmx_age-1] param_natmx_age;
  vector<upper=0>[STEPS_time-artstart_tIDX] dt_log_artrr;
  real<lower=0> sigma_incrate_time_age;
  real<lower=0> sigma_incrate_time;
  real<lower=0> sigma_incrate_age;
  real<lower=0> sigma_natmx_time;
  real<lower=0> sigma_natmx_age;
  real<lower=0> sigma_art;
  real<lower=0>  hivsurv_shape;
  real hivsurv_scale_b0_centered;
  real hivsurv_scale_b1_centered;
}
transformed parameters{

  vector[nk_incrate_time] coef_incrate_time;
  row_vector[nk_incrate_age] coef_incrate_age;
  vector[nk_natmx_age] coef_natmx_age;
  real hivsurv_scale_b0; // intercept for log-scale parameter at age 30
  real hivsurv_scale_b1; // slope for log-scale parameter per 10 years age change

  coef_incrate_time = row_means(coef_incrate_time_age);
  coef_incrate_age = col_means(coef_incrate_time_age);
    
  for(i in 1:nk_natmx_age)
    if (i < fixcoef_natmx_age){
      coef_natmx_age[i] = param_natmx_age[i];
    } else if (i == fixcoef_natmx_age) {
      coef_natmx_age[i] = -sum(param_natmx_age);
    } else {
      coef_natmx_age[i] = param_natmx_age[i-1];
    }

  // Informative prior derived based on Todd et al. (2007) Weibull regression.
  hivsurv_scale_b0 = 0.25*hivsurv_scale_b0_centered + 2.55; // prior ~normal(2.55, 0.25)
  hivsurv_scale_b1 = 0.05*hivsurv_scale_b1_centered + -0.2; // prior ~normal(-0.2, 0.05)
}
model {

  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  sigma_incrate_time_age ~ cauchy(0, var_incrate_time_age);
  sigma_incrate_time ~ cauchy(0, var_incrate_time);
  sigma_incrate_age ~ cauchy(0, var_incrate_age);
  sigma_natmx_time ~ cauchy(0, var_natmx_time);
  sigma_natmx_age ~ cauchy(0, var_natmx_age);
  sigma_art ~ cauchy(0, var_art);

  //////////////////////
  //  Spline penalty  //
  //////////////////////

  {
    matrix[nk_incrate_time, nk_incrate_age] resid_incrate_time_age;
    vector[nk_incrate_time*nk_incrate_age] vec_resid_incrate_time_age;

    for(j in 1:nk_incrate_age)
      for(i in 1:nk_incrate_time)
	resid_incrate_time_age[i,j] = coef_incrate_time_age[i,j] - coef_incrate_time[i] - coef_incrate_age[j];

    vec_resid_incrate_time_age = to_vector(resid_incrate_time_age);
    target += -(nk_incrate_time-1)*(nk_incrate_age-1)*log(sigma_incrate_time_age) -
		       1/(2*sigma_incrate_time_age*sigma_incrate_time_age) * (vec_resid_incrate_time_age' * Pcar_prec_incrate * vec_resid_incrate_time_age);

    D_incrate_time * coef_incrate_time ~ normal(0, sigma_incrate_time);
    D_incrate_age * to_vector(coef_incrate_age) ~ normal(0, sigma_incrate_age);
    
    D_incrate_time * coef_incrate_time_young ~ normal(0, sigma_incrate_time); // DIFFERENT SIGMA FOR YOUNGS??
    
    D_natmx_time * coef_natmx_time ~ normal(0, sigma_natmx_time);
    D_natmx_age * coef_natmx_age ~ normal(0, sigma_natmx_age);
    D_art * dt_log_artrr ~ normal(0, sigma_art);
  }

  ////////////////////////////////////
  // Prior on survival distribution //
  ////////////////////////////////////

  // Informative prior in transformed parameters
  hivsurv_scale_b0_centered ~ normal(0, 1); 
  hivsurv_scale_b1_centered ~ normal(0, 1);  

  // Informative prior on Weibull shape parameter -- mean 2, 95% quantiles (1.0, 3.3)
  hivsurv_shape ~ gamma(12.0, 6.0);

#include /chunks/le_likelihood.stan
}
