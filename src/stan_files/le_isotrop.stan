functions {
#include "chunks/le_functions.stan"
}
data {
#include "chunks/le_data.stan"
}
parameters {

  matrix[nk_incrate_time, nk_incrate_age] coef_incrate_time_age[NSITES];
  vector[nk_incrate_time] coef_incrate_time_young[NSITES];
  vector[nk_natmx_time] coef_natmx_time[NSITES];
  vector[nk_natmx_age-1] param_natmx_age[NSITES];
  vector<upper=0>[STEPS_time*NSITES-sum(artstart_tIDX)] dt_log_artrr;
  real<lower=0> sigma_incrate_time_age;
  real<lower=0> sigma_natmx_time;
  real<lower=0> sigma_natmx_age;
  real<lower=0> sigma_art;
  real<lower=0>  hivsurv_shape[NSITES];
  real hivsurv_scale_b0_centered[NSITES];
  real hivsurv_scale_b1_centered[NSITES];
  
  real<lower=0> sigma_incrate_time;
}
transformed parameters{

  vector[nk_natmx_age] coef_natmx_age[NSITES];
  real hivsurv_scale_b0[NSITES]; // intercept for log-scale parameter at age 30
  real hivsurv_scale_b1[NSITES]; // slope for log-scale parameter per 10 years age change

  for(k in 1:NSITES) {
    for(i in 1:nk_natmx_age)
      if (i < fixcoef_natmx_age){
        coef_natmx_age[k,i] = param_natmx_age[k,i];
      } else if (i == fixcoef_natmx_age) {
        coef_natmx_age[k,i] = -sum(param_natmx_age[k]);
      } else {
        coef_natmx_age[k,i] = param_natmx_age[k,i-1];
      }
    
    // Informative prior derived based on Todd et al. (2007) Weibull regression.
    hivsurv_scale_b0[k] = 0.25*hivsurv_scale_b0_centered[k] + 2.55; // prior ~normal(2.55, 0.25)
    hivsurv_scale_b1[k] = 0.05*hivsurv_scale_b1_centered[k] + -0.2; // prior ~normal(-0.2, 0.05)
  }
}
model {


  ////////////////////////////////
  //  Priors on variance terms  //
  ////////////////////////////////

  sigma_incrate_time_age ~ cauchy(0, var_incrate_time_age);
  sigma_natmx_time ~ cauchy(0, var_natmx_time);
  sigma_natmx_age ~ cauchy(0, var_natmx_age);
  sigma_art ~ cauchy(0, var_art);
  sigma_incrate_time ~ cauchy(0, var_incrate_time);

  //////////////////////
  //  Spline penalty  //
  //////////////////////

  {
    vector[nk_incrate_time*nk_incrate_age] vec_coef_incrate_time_age[NSITES];
    int ind_art_x_s;
    int ind_art_x_e;
    int ind_art_y_s;
    int ind_art_y_e;

    for (k in 1:NSITES) {
      vec_coef_incrate_time_age[k] = to_vector(coef_incrate_time_age[k]);
      target += -nk_incrate_time*nk_incrate_age*log(sigma_incrate_time_age) -
		       1/(2*sigma_incrate_time_age*sigma_incrate_time_age) * (vec_coef_incrate_time_age[k]' * Pcar_prec_incrate * vec_coef_incrate_time_age[k]);
		  
		  D_incrate_time * coef_incrate_time_young[k] ~ normal(0, sigma_incrate_time); // DIFFERENT SIGMA FOR YOUNGS??
		  D_natmx_time * coef_natmx_time[k] ~ normal(0, sigma_natmx_time);
      D_natmx_age * coef_natmx_age[k] ~ normal(0, sigma_natmx_age);
            // sparse matrices instead??
      if(k==1) {
        ind_art_x_s = 1 + STEPS_time*(k-1) - pen_ord_art*(k-1);
      } else {
        ind_art_x_s = 1 + STEPS_time*(k-1) - sum(artstart_tIDX[1:(k-1)]) - pen_ord_art*(k-1);
      }
      ind_art_x_e = STEPS_time*k - sum(artstart_tIDX[1:k]) - pen_ord_art*k;
      if (k==1) {
        ind_art_y_s = 1 + STEPS_time*(k-1);
      } else {
        ind_art_y_s = 1 + STEPS_time*(k-1) - sum(artstart_tIDX[1:(k-1)]);
      }
      ind_art_y_e = STEPS_time*k - sum(artstart_tIDX[1:k]);
      D_art[ind_art_x_s:ind_art_x_e,ind_art_y_s:ind_art_y_e] * dt_log_artrr[ind_art_y_s:ind_art_y_e] ~ normal(0, sigma_art);
    }
  }
  
  ////////////////////////////////////
  // Prior on survival distribution //
  ////////////////////////////////////

  for (k in 1:NSITES){
    // Informative prior in transformed parameters
    hivsurv_scale_b0_centered[k] ~ normal(0, 1); 
    hivsurv_scale_b1_centered[k] ~ normal(0, 1);
    
    // Informative prior on Weibull shape parameter -- mean 2, 95% quantiles (1.0, 3.3)
    hivsurv_shape[k] ~ gamma(12.0, 6.0);
  }

#include /chunks/le_likelihood.stan
}
