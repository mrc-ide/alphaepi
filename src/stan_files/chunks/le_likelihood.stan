  ///////////////////////////////////////////////////////
  //  Construct incidence rate and mortality matrices  //
  ///////////////////////////////////////////////////////

  {
    matrix[STEPS_time-1, STEPS_age-1] incrateMID_time_age[NSITES];
    matrix[STEPS_time, STEPS_age] cumavoid_time_age[NSITES];
    matrix[STEPS_time-1, STEPS_age-1] cumavoidMID_time_age[NSITES];
    matrix[STEPS_time, STEPS_age] natmx_time_age[NSITES];
    matrix[STEPS_time, STEPS_age] natsurv_time_age[NSITES];
    vector[STEPS_time] artrr[NSITES]; // THESE WILL BE WEIRD
    vector[STEPS_time-1] artrr_MID[NSITES]; // THESE WILL BE WEIRD

    // HIV survival model
    vector[STEPS_age-1] hivsurv_scale_a0[NSITES];
    
    matrix[STEPS_time-1, STEPS_age-1] hivmx_dur_a0[NSITES];     // sequenced [1:DUR, 1:STEPS_age]
    matrix[STEPS_time-1, STEPS_age-1] hivsurv_dur_a0[NSITES];   // sequenced [1:DUR, 1:STEPS_age]
    matrix[STEPS_time-1, STEPS_age-1] hivmxMID_dur_a0[NSITES];   // sequenced [1:(STEPS_time-1), 1:STEPS_age]
    
    for (i in 1:NSITES) {
        incrateMID_time_age[i] = rep_matrix(1e-100, STEPS_time-1, STEPS_age-1);
        incrateMID_time_age[i, , test_aIDX:(STEPS_age-1)] = exp(Xmid_incrate_time * coef_incrate_time_age[i] * Xmid_incrate_age');
        incrateMID_time_age[i, , 1] = exp(Xmid_incrate_time * coef_incrate_time_young[i]);
        cumavoid_time_age[i] = exp(-dt*diagCumSum(incrateMID_time_age[i]));
        cumavoidMID_time_age[i] = block(cumavoid_time_age[i], 1, 1, STEPS_time-1, STEPS_age-1) .* exp(-dt/2*incrateMID_time_age[i]);
        
        natmx_time_age[i] = exp(X_natmx_time * coef_natmx_time[i]) * exp(X_natmx_age * coef_natmx_age[i])';
        natsurv_time_age[i] = exp(-dt*diagCumSum(exp(Xmid_natmx_time * coef_natmx_time[i]) * exp(Xmid_natmx_age * coef_natmx_age[i])'));
        
        hivsurv_scale_a0[i] = exp(hivsurv_scale_b0[i] + X_hivsurv_age * hivsurv_scale_b1[i]);
        
        hivmx_dur_a0[i] = create_hivmx_dur_a0(hivsurv_shape[i], hivsurv_scale_a0[i], STEPS_time-1, dt);
        hivsurv_dur_a0[i] = create_log_hivsurv_dur_a0(hivsurv_shape[i], hivsurv_scale_a0[i], STEPS_time-1, dt);
        hivmxMID_dur_a0[i] = diff_hivmxMID_dur_a0(hivsurv_dur_a0[i], dt);
        hivsurv_dur_a0[i] = exp(hivsurv_dur_a0[i]);
    }
    
    {
      // COME BACK TO THIS - GOING TO BE WEIRD!!
      vector[STEPS_time*NSITES-sum(artstart_tIDX)] log_artrr; // this is weird - since artrr is length STEPS_time, but left for time being
      int art_ind_s;
      int art_ind_e;
      for(k in 1:NSITES){
          if (k==1) {
            art_ind_s = 1 + STEPS_time*(k-1);
          } else {
            art_ind_s = 1 + STEPS_time*(k-1) - sum(artstart_tIDX[1:(k-1)]);
          }
          art_ind_e = STEPS_time*k - sum(artstart_tIDX[1:k]);
          log_artrr[art_ind_s:art_ind_e] = dt*cumulative_sum(dt_log_artrr[art_ind_s:art_ind_e]);
          for(i in 1:(STEPS_time-artstart_tIDX[k])){
              artrr[k,artstart_tIDX[k]+i] =  exp(log_artrr[i+art_ind_s-1]);
              artrr_MID[k,artstart_tIDX[k]+i-1] = exp(log_artrr[i+art_ind_s-1] - dt/2*dt_log_artrr[i+art_ind_s-1]);
              }
      }
      
      
    }
    
    ///////////////////////////////////////
    //  Calculate individual likelihood  //
    ///////////////////////////////////////

    // does this need to be in its own little {}??
    {
      int cinds; // cohort start ind
      int cinde; // cohort end ind
    
      int einds; // exit start ind
      int einde; // exit end ind

      int ainds; // aggregate start ind
      int ainde; // aggregate end ind

      for(i in 1:NSITES) {
        cinds = 1 + sum(NCOH[1:(i-1)]);
        cinde = sum(NCOH[1:i]);
        
        einds = 1 + sum(NEXIT[1:(i-1)]);
        einde = sum(NEXIT[1:i]);
        
        ainds = 1 + sum(NAGGR[1:(i-1)]);
        ainde = sum(NAGGR[1:i]);
        
        //art start needs to get fixed
        
        target += calc_ll_cohexit(coh_cIDX[cinds:cinde], coh_minexpose_tIDX[cinds:cinde], coh_maxexpose_tIDX[cinds:cinde], 
                                coh_nexit[cinds:cinde], exdat_tIDX[einds:einde], exdat_minexpose_tIDX[einds:einde], 
                                exdat_maxexpose_tIDX[einds:einde], exdat_ndat[einds:einde], aggr_exposestart_tIDX[ainds:ainde], 
                                aggr_exposeend_tIDX[ainds:ainde], aggr_death[ainds:ainde], aggr_deathinterv[ainds:ainde], 
                                aggr_deathinterv_DUR[ainds:ainde], aggr_hivpos[ainds:ainde], aggr_nrepl[ainds:ainde],
                                cumavoid_time_age[i], cumavoidMID_time_age[i], incrateMID_time_age[i],
                                hivsurv_dur_a0[i], hivmx_dur_a0[i], hivmxMID_dur_a0[i], artrr[i], 
                                artrr_MID[i], artstart_tIDX[i], natsurv_time_age[i], natmx_time_age[i], dt);
      }
    }
  }
