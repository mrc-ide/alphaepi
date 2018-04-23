  vector log_inv_logit_vec(vector x){

    vector[rows(x)] val;

    for(i in 1:rows(x))
      val[i] = log_inv_logit(x[i]);

    return val;
  }

  vector inv_logit_vec(vector x){

    vector[rows(x)] val;

    for(i in 1:rows(x))
      val[i] = inv_logit(x[i]);

    return val;
  }

  matrix diagCumSum(matrix x){

    matrix[rows(x)+1, cols(x)+1] val;

    for(i in 1:rows(val))
      val[i,1] = 0;
    for(j in 2:cols(val)){
      val[1,j] = 0;
      for(i in 2:rows(val))
        val[i,j] = val[i-1,j-1] + x[i-1,j-1];
    }

    return val;
  }

  row_vector col_means(matrix X) {
    row_vector[cols(X)] val;
    val = rep_row_vector(1.0/rows(X), rows(X)) * X;
    return val;
  }

  vector row_means(matrix X) {
    vector[rows(X)] val;
    val = X * rep_vector(1.0/cols(X), cols(X));
    return val;
  }

  real weibull_hazard(real y, real alpha, real sigma){
    return (alpha/sigma)*(y/sigma)^(alpha-1);
  }

  
  ////////////////////////////////////////////////////////////////////////////////
  // create_hivmx_dur_a0: functions to create survival distribution and hazard  //
  //                      by duration of infection and age at seroconversion.   //
  ////////////////////////////////////////////////////////////////////////////////

  matrix create_hivmx_dur_a0(real shape, vector scale_a0, int steps_dur, real dt){
    matrix[steps_dur, rows(scale_a0)] y;
    for(j in 1:rows(scale_a0)){
      for(i in 1:steps_dur)
	y[i,j] = weibull_hazard(i*dt-dt/2, shape, scale_a0[j]);
    }
    return y;
  }
  
  matrix create_log_hivsurv_dur_a0(real shape, vector scale_a0, int steps_dur, real dt){
    matrix[steps_dur, rows(scale_a0)] y;
    for(j in 1:rows(scale_a0)){
      for(i in 1:steps_dur)
	y[i,j] = weibull_ccdf_log(i*dt-dt/2, shape, scale_a0[j]);
    }
    return y;
  }
  
  matrix diff_hivmxMID_dur_a0(matrix x, real dt){
    matrix[rows(x), cols(x)] y;
    y[1,] = x[1,];
    for(j in 1:cols(x))
      for(i in 2:rows(x))
        y[i,j] = x[i,j] - x[i-1,j];
    y = -y/dt;
    return y;
  }

  //////////////////////////////////////////////////////////////////////////////
  // calc_phivsurv: probability of surviving until exit_tIDX given infection  //
  //                at midpoint of (tIDX, aIDX) + 0:(exposeDUR-1).            //
  //////////////////////////////////////////////////////////////////////////////

  vector calc_phivsurv(int tIDX, int aIDX, int exposeDUR, int exit_tIDX,
                       matrix hivsurv_dur_a0, matrix hivmxMID_dur_a0, vector artrr_MID, int artstart_tIDX, real dt){
    vector[exposeDUR] phivsurv;

    if(exit_tIDX > artstart_tIDX){

      int i_tIDX;
      int i_aIDX;

      for(ii in 1:exposeDUR){
        i_tIDX = tIDX + ii - 1;
        i_aIDX = aIDX + ii - 1;

        if(i_tIDX < artstart_tIDX)
          phivsurv[ii] = hivsurv_dur_a0[artstart_tIDX-i_tIDX, i_aIDX] *
            exp(-dt*dot_product(sub_col(hivmxMID_dur_a0, artstart_tIDX-i_tIDX+1, i_aIDX, exit_tIDX - artstart_tIDX),
                                segment(artrr_MID, artstart_tIDX, exit_tIDX - artstart_tIDX)));
        else
          phivsurv[ii] = exp(-dt*dot_product(sub_col(hivmxMID_dur_a0, 1, i_aIDX, exit_tIDX - i_tIDX), segment(artrr_MID, i_tIDX, exit_tIDX - i_tIDX)));
      }
    } else{
      int taoIDX;
      taoIDX = exit_tIDX - tIDX + 1;
      for(ii in 1:exposeDUR)
    	phivsurv[ii] = hivsurv_dur_a0[taoIDX-ii, aIDX+ii-1];
    }

    return(phivsurv);
  }
  
  vector calc_hivmx(int tIDX, int aIDX, int exposeDUR, int exit_tIDX,
                    matrix hivmx_dur_a0, vector artrr, int artstart_tIDX){
    int taoIDX;
    vector[exposeDUR] hivmx;

    taoIDX = exit_tIDX - tIDX + 1;
    
    for(ii in 1:exposeDUR)
      hivmx[ii] = hivmx_dur_a0[taoIDX-ii, aIDX+ii-1];
    
    if(exit_tIDX > artstart_tIDX)
      hivmx = hivmx * artrr[exit_tIDX];
    
    return(hivmx);
  }


  //////////////////////////////////////////////////////////////////////////
  // calc_log_psurventry: calculate log probability of survival to entry  //
  //                      into the cohort (account for left truncation).  //
  //////////////////////////////////////////////////////////////////////////

  // Note: psurventry = phivp + phivn with hivpos = 0 and death = 0.
  //       Coded as separate function for efficiency of omitting unneeded calculations.

  real calc_log_psurventry(int entry_tIDX, int entry_aIDX,
                           matrix cumavoid_time_age, matrix cumavoidMID_time_age, matrix incrateMID_time_age,
                           matrix hivsurv_dur_a0, matrix hivmxMID_dur_a0, vector artrr_MID, int artstart_tIDX,
                           matrix natsurv_time_age, real dt){

    real log_psurventry;
    int expose_DUR;
    int expose_tIDX;
    int expose_aIDX;

    expose_DUR = min(entry_tIDX, entry_aIDX) - 1;
    expose_tIDX = entry_tIDX - expose_DUR;
    expose_aIDX = entry_aIDX - expose_DUR;

    if(expose_DUR > 0){
      vector[expose_DUR] phivsurv;
      vector[expose_DUR] log_psurventry_i;

      phivsurv = calc_phivsurv(expose_tIDX, expose_aIDX, expose_DUR, entry_tIDX,
    				     hivsurv_dur_a0, hivmxMID_dur_a0, artrr_MID, artstart_tIDX, dt);
      for(ii in 1:expose_DUR)
    	log_psurventry_i[ii] = cumavoidMID_time_age[expose_tIDX+ii-1, expose_aIDX+ii-1] *
    	                       incrateMID_time_age[expose_tIDX+ii-1, expose_aIDX+ii-1] *
    	                       phivsurv[ii];

      log_psurventry = log((dt*sum(log_psurventry_i) + cumavoid_time_age[entry_tIDX, entry_aIDX])*natsurv_time_age[entry_tIDX, entry_aIDX]);
    } else
      log_psurventry = 0.0;

    return(log_psurventry);
  }


  /////////////////////////////////////////////////////////////////////////////
  //  calc_phivn: probability of remaining HIV- and survival status at exit  //
  /////////////////////////////////////////////////////////////////////////////

  real calc_phivn(int exit_tIDX, int exit_aIDX, int death, int hivpos,
                  matrix cumavoid_time_age, matrix natsurv_time_age, matrix natmx_time_age){
    real phivn;

    if(hivpos)
      phivn = 0;
    else{
      phivn = cumavoid_time_age[exit_tIDX, exit_aIDX] * natsurv_time_age[exit_tIDX, exit_aIDX];
      if(death)
        phivn = phivn*natmx_time_age[exit_tIDX, exit_aIDX];
    }

    return(phivn);
  }


  ////////////////////////////////////////////////////////////////////////////
  //  calc_phivp: probability of becoming HIV+ and survival status at exit  //
  ////////////////////////////////////////////////////////////////////////////

  real calc_phivp(int exit_tIDX, int exit_aIDX, int expose_tIDX, int expose_aIDX, int expose_DUR, int death,
                  matrix cumavoidMID_time_age, matrix incrateMID_time_age,
                  matrix hivsurv_dur_a0, matrix hivmx_dur_a0, matrix hivmxMID_dur_a0,
                  vector artrr, vector artrr_MID, int artstart_tIDX,
                  matrix natsurv_time_age, matrix natmx_time_age,
                  real dt){

    real phivp;

    if(expose_DUR > 0){
      vector[expose_DUR] phivsurv;
      vector[expose_DUR] integrand;
      
      phivsurv = calc_phivsurv(expose_tIDX, expose_aIDX, expose_DUR, exit_tIDX,
       				hivsurv_dur_a0, hivmxMID_dur_a0, artrr_MID, artstart_tIDX, dt);
	
      for(ii in 1:expose_DUR)
    	integrand[ii] = cumavoidMID_time_age[expose_tIDX+ii-1, expose_aIDX+ii-1] *
    	                 incrateMID_time_age[expose_tIDX+ii-1, expose_aIDX+ii-1] *
    	                 phivsurv[ii];

      if(death)
        integrand = integrand .* (calc_hivmx(expose_tIDX, expose_aIDX, expose_DUR, exit_tIDX,
                                              hivmx_dur_a0, artrr, artstart_tIDX) +
                                   natmx_time_age[exit_tIDX, exit_aIDX]);

      phivp = dt * sum(integrand) * natsurv_time_age[exit_tIDX, exit_aIDX];
    } else
      phivp = 0.0;

    return phivp;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //  create_phivp_mat: full matrix for probability of becoming HIV+ and surviving  //
  ////////////////////////////////////////////////////////////////////////////////////

  matrix create_phivp_mat(matrix cumavoidMID_time_age, matrix incrateMID_time_age,
			  matrix hivsurv_dur_a0, matrix hivmxMID_dur_a0,
			  vector artrr_MID, int artstart_tIDX,
			  matrix natsurv_time_age, real dt){

    matrix[rows(natsurv_time_age), cols(natsurv_time_age)] phivp;
    int steps_time;
    int steps_age;

    steps_time = rows(natsurv_time_age);
    steps_age = cols(natsurv_time_age);
    

    for(cidx in (1-steps_age):(steps_time-1)){

      int min_tidx;
      int max_tidx;
      vector[steps_time-1] incdens_coh;
      vector[steps_time-1] phivsurv_coh;

      min_tidx = max(1, 1+cidx);
      max_tidx = min(steps_time, steps_age+cidx);

      for(tidx in min_tidx:max_tidx)
	if(tidx > min_tidx){
	  int tmin1;
	  int exp_DUR;
	  exp_DUR = tidx - min_tidx;

	  tmin1 = tidx-1;
	  incdens_coh[tmin1] = cumavoidMID_time_age[tmin1, tmin1-cidx] * incrateMID_time_age[tmin1, tmin1-cidx];
	  if(tidx <= artstart_tIDX)
	    for(t0idx in min_tidx:(tidx-1))
	      phivsurv_coh[t0idx] = hivsurv_dur_a0[tidx-t0idx, t0idx-cidx];
	  else{
	    phivsurv_coh[tmin1] = 1.0;
	    for(t0idx in min_tidx:(tidx-1))
	      phivsurv_coh[t0idx] = phivsurv_coh[t0idx] * exp(-dt*hivmxMID_dur_a0[tidx-t0idx, t0idx-cidx] * artrr_MID[tmin1]);
	  }
	    
	  phivp[tidx, tidx-cidx] = dt * sum(segment(incdens_coh, min_tidx, exp_DUR) .* segment(phivsurv_coh, min_tidx, exp_DUR)) * natsurv_time_age[tidx, tidx-cidx];
	} else
	  phivp[tidx, tidx-cidx] = 0.0;
    }

    return phivp;
  }
                  

  
  /////////////////////////////////////////////////////
  //  calc_ll_ind: log-likelihood for an individual  //
  /////////////////////////////////////////////////////

  real calc_ll_ind(int entry_tIDX, int entry_aIDX, int exit_tIDX, int exit_aIDX,
                   int expose_tIDX, int expose_aIDX, int expose_DUR,
		   int death, int deathinterv, int deathinterv_DUR, int hivpos,
                   matrix cumavoid_time_age, matrix cumavoidMID_time_age, matrix incrateMID_time_age,
                   matrix hivsurv_dur_a0, matrix hivmx_dur_a0, matrix hivmxMID_dur_a0,
                   vector artrr, vector artrr_MID, int artstart_tIDX,
                   matrix natsurv_time_age, matrix natmx_time_age,
                   real dt){

    real log_psurventry;
    real phivn;
    real phivp;
    real phivn_d;
    real phivp_d;


    log_psurventry = calc_log_psurventry(entry_tIDX, entry_aIDX,
					 cumavoid_time_age, cumavoidMID_time_age, incrateMID_time_age,
					 hivsurv_dur_a0, hivmxMID_dur_a0, artrr_MID, artstart_tIDX,
					 natsurv_time_age, dt);

    phivn = calc_phivn(exit_tIDX, exit_aIDX, death, hivpos,
                       cumavoid_time_age, natsurv_time_age, natmx_time_age);

    phivp = calc_phivp(exit_tIDX, exit_aIDX, expose_tIDX, expose_aIDX, expose_DUR, death,
		       cumavoidMID_time_age, incrateMID_time_age,
		       hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
		       artrr, artrr_MID, artstart_tIDX,
		       natsurv_time_age, natmx_time_age,
		       dt);

    if(deathinterv){

      // Note: if person is HIV+, exposure interval ends at HIV test, otherwise it extends to exit,
      //       and needs to be extended to the end of the death interval.
      int expose_dinterv_DUR;
      if(hivpos)
	expose_dinterv_DUR = expose_DUR;
      else
	expose_dinterv_DUR = expose_DUR + deathinterv_DUR;
      
      phivn_d = calc_phivn(exit_tIDX+deathinterv_DUR, exit_aIDX++deathinterv_DUR, death, hivpos,
                        cumavoid_time_age, natsurv_time_age, natmx_time_age);

      phivp_d = calc_phivp(exit_tIDX+deathinterv_DUR, exit_aIDX+deathinterv_DUR, expose_tIDX, expose_aIDX, expose_dinterv_DUR, death,
			    cumavoidMID_time_age, incrateMID_time_age,
			    hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
			    artrr, artrr_MID, artstart_tIDX,
			    natsurv_time_age, natmx_time_age,
			    dt);
    } else {
      phivn_d = 0;
      phivp_d = 0;
    }

    return(log(phivn+phivp-phivn_d-phivp_d) - log_psurventry);
  }


  /////////////////////////////////////////////////////////////////
  //  calc_ll_coh: calculate the log likelihood for cohort data  //
  /////////////////////////////////////////////////////////////////

  vector calc_ll_coh(int[] entry_tIDX, int[] entry_aIDX, int[] exit_tIDX, int[] exit_aIDX,
                     int[] expose_tIDX, int[] expose_aIDX, int[] expose_DUR,
		     int[] death, int[] deathinterv, int[] deathinterv_DUR, int[] hivpos,
                     matrix cumavoid_time_age, matrix cumavoidMID_time_age, matrix incrateMID_time_age,
                     matrix hivsurv_dur_a0, matrix hivmx_dur_a0, matrix hivmxMID_dur_a0,
                     vector artrr, vector artrr_MID, int artstart_tIDX,
                     matrix natsurv_time_age, matrix natmx_time_age,
                     real dt){

    vector[size(entry_tIDX)] ll_coh;

    for(i in 1:size(entry_tIDX))
      ll_coh[i] = calc_ll_ind(entry_tIDX[i], entry_aIDX[i], exit_tIDX[i], exit_aIDX[i],
			      expose_tIDX[i], expose_aIDX[i], expose_DUR[i],
			      death[i], deathinterv[i], deathinterv_DUR[i], hivpos[i],
			      cumavoid_time_age, cumavoidMID_time_age, incrateMID_time_age,
			      hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
                               artrr, artrr_MID, artstart_tIDX,
			      natsurv_time_age, natmx_time_age,
			      dt);

    return ll_coh;
  }
  
  real calc_ll_coh_aggr(int[] cohIDX, int[] exit_tIDX, int[] exposestart_tIDX, int[] exposeend_tIDX,
			int[] death, int[] deathinterv, int[] deathinterv_DUR, int[] hivpos, vector nrepl,
			matrix cumavoid_time_age, matrix cumavoidMID_time_age, matrix incrateMID_time_age,
			matrix hivsurv_dur_a0, matrix hivmx_dur_a0, matrix hivmxMID_dur_a0,
			vector artrr, vector artrr_MID, int artstart_tIDX,
			matrix natsurv_time_age, matrix natmx_time_age,
			real dt){

    vector[rows(nrepl)] ll_one;
    
    for(i in 1:rows(ll_one)){
      real phivp;
      real phivn;
      real phivp_d;
      real phivn_d;

      phivp = calc_phivp(exit_tIDX[i], exit_tIDX[i]-cohIDX[i], exposestart_tIDX[i], exposestart_tIDX[i]-cohIDX[i], exposeend_tIDX[i]-exposestart_tIDX[i]+1, death[i],
			 cumavoidMID_time_age, incrateMID_time_age,
			 hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
			  artrr, artrr_MID, artstart_tIDX,
			 natsurv_time_age, natmx_time_age,
			 dt);
      phivn = calc_phivn(exit_tIDX[i], exit_tIDX[i]-cohIDX[i], death[i], hivpos[i],
			 cumavoid_time_age, natsurv_time_age, natmx_time_age);

      if(deathinterv[i]){
	  // If person is HIV+, exposure interval ends at HIV test, otherwise it extends to exit,
	  // and needs to be extended to the end of the death interval.
	  int expose_dinterv_DUR;
	  if(hivpos[i])
	    expose_dinterv_DUR = exposeend_tIDX[i]-exposestart_tIDX[i]+1;
	  else
	    expose_dinterv_DUR = exposeend_tIDX[i]-exposestart_tIDX[i]+1 + deathinterv_DUR[i];

	  phivp_d = calc_phivp(exit_tIDX[i]+deathinterv_DUR[i], exit_tIDX[i]+deathinterv_DUR[i]-cohIDX[i],
			       exposestart_tIDX[i], exposestart_tIDX[i]-cohIDX[i], expose_dinterv_DUR, death[i],
			       cumavoidMID_time_age, incrateMID_time_age,
			       hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
			       artrr, artrr_MID, artstart_tIDX,
			       natsurv_time_age, natmx_time_age,
			       dt);
	  phivn_d = calc_phivn(exit_tIDX[i]+deathinterv_DUR[i], exit_tIDX[i]+deathinterv_DUR[i]-cohIDX[i], death[i], hivpos[i],
			       cumavoid_time_age, natsurv_time_age, natmx_time_age);
      } else {
	phivp_d = 0;
	phivn_d = 0;
      }
      
      ll_one[i] = log(phivp+phivn-phivp_d-phivn_d);
    }
    
    return dot_product(ll_one, nrepl);
  }


  // calculate likelihood by taking each cohort in order and then summing along exit times sequentially
  real calc_ll_cohexit(int[] coh_cIDX, int[] coh_minexpose_tIDX, int[] coh_maxexpose_tIDX, int[] coh_nexit,
		       int[] exdat_tIDX, int[] exdat_minexpose_tIDX, int[] exdat_maxexpose_tIDX, int[] exdat_ndat,
		       int[] aggr_exposestart_tIDX, int[] aggr_exposeend_tIDX,
		       int[] aggr_death, int[] aggr_deathinterv, int[] aggr_deathinterv_DUR, int[] aggr_hivpos, vector aggr_nrepl,
		       matrix cumavoid_time_age, matrix cumavoidMID_time_age, matrix incrateMID_time_age,
		       matrix hivsurv_dur_a0, matrix hivmx_dur_a0, matrix hivmxMID_dur_a0,
		       vector artrr, vector artrr_MID, int artstart_tIDX,
		       matrix natsurv_time_age, matrix natmx_time_age,
		       real dt){
    
    vector[rows(aggr_nrepl)] ll_one;
    int aggr_i;
    int exdat_i;
    int max_aggr_i;
    int max_exdat_i;
    aggr_i = 0;
    exdat_i = 0;
    max_aggr_i = 0;
    max_exdat_i = 0;

    for(coh_i in 1:size(coh_cIDX)){
      
      int cidx;   // cohort id
      vector[rows(cumavoidMID_time_age)] incdens_coh;
      vector[rows(cumavoidMID_time_age)] phivsurv_coh;
      int eidx_prev; // previous exit time (for summing survival)
      
      cidx = coh_cIDX[coh_i];
      eidx_prev = 0; // initialize to 0

      for(tidx in coh_minexpose_tIDX[coh_i]:coh_maxexpose_tIDX[coh_i])
	incdens_coh[tidx] = cumavoidMID_time_age[tidx, tidx-cidx] * incrateMID_time_age[tidx, tidx-cidx];

      for(ii in artstart_tIDX:coh_maxexpose_tIDX[coh_i])
	phivsurv_coh[ii] = 1.0;
      
      max_exdat_i = max_exdat_i + coh_nexit[coh_i];
      while(exdat_i < max_exdat_i){

	int eidx;
	real phivn_exit;

	exdat_i = exdat_i + 1;
	eidx = exdat_tIDX[exdat_i];

	phivn_exit = cumavoid_time_age[eidx, eidx-cidx] * natsurv_time_age[eidx, eidx-cidx];

	// calculate phivsurv_coh: probability of surviving to eidx for someone in cidx conditional
	//                         on time of infection.
	if(eidx <= artstart_tIDX){ // lookup phivsurv
	  for(t0idx in exdat_minexpose_tIDX[exdat_i]:exdat_maxexpose_tIDX[exdat_i])
	    phivsurv_coh[t0idx] = hivsurv_dur_a0[eidx-t0idx, t0idx-cidx];

	} else {

	  if(eidx_prev < artstart_tIDX){
	    for(t0idx in exdat_minexpose_tIDX[exdat_i]:(artstart_tIDX-1))
	      phivsurv_coh[t0idx] = hivsurv_dur_a0[artstart_tIDX - t0idx, t0idx-cidx];
	    eidx_prev = artstart_tIDX;
	  }

	  // calculate probability of surviving from eidx_prev to eidx for each t0idx.
	  for(t0idx in exdat_minexpose_tIDX[exdat_i]:exdat_maxexpose_tIDX[exdat_i]){

	      int a0idx;
	      int t0idxminus1;
	      real cumhaz;
	      
	      a0idx = t0idx-cidx;
	      t0idxminus1 = t0idx-1;
	      cumhaz = 0;

	      for(jj in max(t0idx, eidx_prev):(eidx-1))
		cumhaz = cumhaz + hivmxMID_dur_a0[jj-t0idxminus1, a0idx] * artrr_MID[jj];

	      phivsurv_coh[t0idx] = phivsurv_coh[t0idx] * exp(-dt*cumhaz);
	  }
	}

	eidx_prev = eidx;
	
	max_aggr_i = max_aggr_i + exdat_ndat[exdat_i];
	while(aggr_i < max_aggr_i){
	  
	  real phivp;
	  real phivn;
	  int i_expose_tIDX;
	  int i_expose_DUR;

	  aggr_i = aggr_i+1;
	  
	  // calc phivn
	  if(aggr_hivpos[aggr_i])
	    phivn = 0;
	  else{
	    phivn = phivn_exit;
	    if(aggr_death[aggr_i])
	      phivn = phivn*natmx_time_age[eidx, eidx-cidx];
	  }
	  
	  // calc phivp
	  i_expose_tIDX = aggr_exposestart_tIDX[aggr_i];
	  i_expose_DUR = aggr_exposeend_tIDX[aggr_i] - i_expose_tIDX + 1;
	  if(i_expose_DUR > 0){
	    
	    vector[i_expose_DUR] integrand;
	    	    
	    integrand = segment(incdens_coh, i_expose_tIDX, i_expose_DUR) .* segment(phivsurv_coh, i_expose_tIDX, i_expose_DUR);
	    
	    if(aggr_death[aggr_i]){
	      vector[i_expose_DUR] hivmx;
	      for(t0idx in i_expose_tIDX:aggr_exposeend_tIDX[aggr_i])
		hivmx[t0idx-i_expose_tIDX+1] = hivmx_dur_a0[eidx-t0idx, t0idx-cidx];
	      if(eidx > artstart_tIDX)
		hivmx = hivmx * artrr[eidx];
	      
	      integrand = integrand .* (hivmx + natmx_time_age[eidx, eidx-cidx]);
	    }
	    
	    phivp = dt * sum(integrand) * natsurv_time_age[eidx, eidx-cidx];
	    
	  } else
	    phivp = 0.0;
	  
	  if(aggr_deathinterv[aggr_i]){

	      // Note: Challenge with interval censored deaths is that it uses an exit date in the future.
	      // Number of interval censored deaths should be few relative to overall likelihood calculation,
	      // so just use these functions for calculating. But it's probably possible to revise this algorithm,
	      // to handle interval censored deaths more efficiently (e.g. sum forward HIV survival from the
	      // current exit.

	      int dur_i;
	      real phivn_d;
	      real phivp_d;

	      dur_i = aggr_deathinterv_DUR[aggr_i];
	      phivn_d = calc_phivn(eidx+dur_i, eidx+dur_i-cidx, 0, aggr_hivpos[aggr_i], // death = 0, hivpos = 0
				    cumavoid_time_age, natsurv_time_age, natmx_time_age);


	      // If person is HIV+, exposure interval ends at HIV test, otherwise it extends to exit,
	      // and needs to be extended to the end of the death interval.
	      if(!aggr_hivpos[aggr_i])
		i_expose_DUR = i_expose_DUR + dur_i;

	      phivp_d = calc_phivp(eidx+dur_i, eidx+dur_i-cidx, i_expose_tIDX, i_expose_tIDX-cidx, i_expose_DUR, 0, // death = 0
	      			cumavoidMID_time_age, incrateMID_time_age,
	      			hivsurv_dur_a0, hivmx_dur_a0, hivmxMID_dur_a0,
	      			artrr, artrr_MID, artstart_tIDX,
	      			natsurv_time_age, natmx_time_age,
	      			dt);

	      phivn = phivn - phivn_d;
	      phivp = phivp - phivp_d;
	  }
	  
	  ll_one[aggr_i] = log(phivn+phivp);
	}
      }
    }

    return dot_product(ll_one, aggr_nrepl);
  }
