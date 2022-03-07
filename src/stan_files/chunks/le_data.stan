  // STATE SPACE PARAMETERS
  int<lower=1> NSITES;
  
  real<lower=0> dt;
  int<lower=1> STEPS_time;
  int<lower=1> STEPS_age;
  int<lower=1, upper=STEPS_age> STEPS_age_test;
  int<lower=1, upper=STEPS_time> artstart_tIDX[NSITES];
  int<lower=1, upper=STEPS_age> test_aIDX;

  vector[STEPS_time] x_time;
  vector[STEPS_age] x_age;
  
  // COHORT DATA
  int<lower=1> NCOH[NSITES];

  int coh_cIDX[sum(NCOH)]; 
  int<lower=1, upper=STEPS_time> coh_minexpose_tIDX[sum(NCOH)];
  int<lower=1, upper=STEPS_time-1> coh_maxexpose_tIDX[sum(NCOH)];
  int<lower=1> coh_nexit[sum(NCOH)];

  // EXIT DATA
  int<lower=min(NCOH)> NEXIT[NSITES];

  int<lower=1, upper=STEPS_time> exdat_tIDX[sum(NEXIT)];
  int<lower=1, upper=STEPS_time> exdat_minexpose_tIDX[sum(NEXIT)];
  int<lower=1, upper=STEPS_time> exdat_maxexpose_tIDX[sum(NEXIT)];
  int<lower=1> exdat_ndat[sum(NEXIT)];
  
  // AGGREGATE INDIVIDUAL DATA
  int<lower=min(NEXIT)> NAGGR[NSITES];

  int<lower=1, upper=STEPS_time> aggr_exposestart_tIDX[sum(NAGGR)];
  int<lower=1, upper=STEPS_time> aggr_exposeend_tIDX[sum(NAGGR)];
  int<lower=0, upper=1> aggr_hivpos[sum(NAGGR)];
  int<lower=0, upper=1> aggr_death[sum(NAGGR)];
  int<lower=0, upper=1> aggr_deathinterv[sum(NAGGR)];
  int<lower=0, upper=STEPS_time> aggr_deathinterv_DUR[sum(NAGGR)];
  vector[sum(NAGGR)] aggr_nrepl;


  // MODEL PARAMETERS

  // incidence model
  int<lower=5> nk_incrate_time;
  int<lower=5> nk_incrate_age;

  matrix[STEPS_time, nk_incrate_time] X_incrate_time;
  matrix[STEPS_time-1, nk_incrate_time] Xmid_incrate_time;
  matrix[STEPS_age_test, nk_incrate_age] X_incrate_age;
  matrix[STEPS_age_test-1, nk_incrate_age] Xmid_incrate_age;

  int<lower=0> pen_ord_incrate;

  matrix[nk_incrate_time-pen_ord_incrate, nk_incrate_time] D_incrate_time;
  matrix[nk_incrate_age-pen_ord_incrate, nk_incrate_age] D_incrate_age;
  matrix[nk_incrate_time*nk_incrate_age, nk_incrate_time*nk_incrate_age] Pcar_prec_incrate;


  // non-HIV mortality model
  int<lower=5> nk_natmx_time;
  int<lower=5> nk_natmx_age;

  matrix[STEPS_time, nk_natmx_time] X_natmx_time;
  matrix[STEPS_time-1, nk_natmx_time] Xmid_natmx_time;
  matrix[STEPS_age, nk_natmx_age] X_natmx_age;
  matrix[STEPS_age-1, nk_natmx_age] Xmid_natmx_age;

  int<lower=0> pen_ord_natmx_time;
  int<lower=0> pen_ord_natmx_age;

  matrix[nk_natmx_time-pen_ord_natmx_time, nk_natmx_time] D_natmx_time;
  matrix[nk_natmx_age-pen_ord_natmx_age, nk_natmx_age] D_natmx_age;

  int<lower=1, upper=nk_natmx_time> fixcoef_natmx_time;
  int<lower=1, upper=nk_natmx_age> fixcoef_natmx_age;
  int<lower=1, upper=nk_incrate_age> fixcoef_incrate_age;



  // ART model
  
  // matrix[STEPS_time, nk_art] X_art;
  // matrix[STEPS_time-1, nk_art] Xmid_art;

  int<lower=0> pen_ord_art;
  // matrix[STEPS_time-artstart_tIDX-pen_ord_art, STEPS_time-artstart_tIDX] D_art[NSITES];
  matrix[STEPS_time*NSITES-sum(artstart_tIDX)-pen_ord_art*NSITES, STEPS_time*NSITES-sum(artstart_tIDX)] D_art; // this will be weird!

  
  // HIV survival model
  // matrix[STEPS_time-1, STEPS_age-1] hivmx_dur_a0;     // sequenced [1:DUR, 1:STEPS_age]
  // matrix[STEPS_time-1, STEPS_age-1] hivsurv_dur_a0;   // sequenced [1:DUR, 1:STEPS_age]
  // matrix[STEPS_time-1, STEPS_age-1] hivmxMID_dur_a0;   // sequenced [1:(STEPS_time-1), 1:STEPS_age]
  // real hivsurv_shape;
  // vector[STEPS_age-1] hivsurv_scale_a0;

  // Prior variance
  real<lower=0> var_incrate_time_age;
  real<lower=0> var_incrate_time;
  real<lower=0> var_incrate_age;
  real<lower=0> var_natmx_time;
  real<lower=0> var_natmx_age;
  real<lower=0> var_art;

}
transformed data {

  vector[STEPS_age-1] X_hivsurv_age;
  X_hivsurv_age = ((x_age[2:] - dt/2) - 30) / 10; // centered on age 30 at seroconversion, per 10 years of age
