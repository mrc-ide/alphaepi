
#####################################
####  HIV survival distribution  ####
#####################################

hivsurv <- function(tao, a0, param){
  ## tao: time (years) since infection
  ## a0:  age (years) at infection

  ## Note: from Bellan et al. Lancet 2013, based on CASCADE Lancet 2000.
  shape <- 2.3
  scale <- 166/(shape*a0^0.53)
  return(exp(-(tao/scale)^shape))
}

hivsurvhaz <- function(tao, a0, param){
  shape <- 2.3
  scale <- 166/(shape*a0^0.53)
  return(shape/scale*(tao/scale)^(shape-1))
}


########################
####  Prepare data  ####
########################
                                   

prepare.stan.data <- function(sites = NULL, sexes = NULL, dat = NULL, dt = 0.1,
                              min.age = 15.0, max.age = 60.0,
                              min.time = 1980.5, max.time = 2011.5,
                              natmxstart.time, artstart.time,
                              cohortstart.time=min.time, cohortend.time=max.time,
                              ## nk.time = 7, nk.age = 7, nk.natmx = 5, nk.art = 5,
                              k.dt = 5, nk.art=5,
                              hivsurv.shape=2.0,
                              time.pen=TRUE, age.pen=TRUE, cohort.pen=FALSE,
                              pen.ord.incrate=1L, pen.ord.natmx.time=1L, pen.ord.natmx.age=1L, pen.ord.art=1L,
                              nsamp=NULL, hivonly=FALSE, hivelig=FALSE,
                              var_incrate_time_age=2.5, var_incrate_time=2.5, var_incrate_age=2.5, 
                              var_natmx_time=2.5, var_natmx_age=2.5, var_art=2.5){

  ## hivonly: if TRUE, only use HIV test data, don't use any residency episode data (no mortality).
  ## hivelig: indicates individuals only included if they have some HIV status information, so left truncate
  ##          at first HIV status information (inclusion conditional on survival to that point).

  if(is.null(dat)){
    dat <- do.call(rbind, mapply(prepare.interval.data, sites, sexes, min.age, max.age, cohortstart.time, cohortend.time, hivonly, hivelig, SIMPLIFY=FALSE))
  }

  ## Select sub-sample 
  if(is.null(nsamp))
    samp <- 1:nrow(dat)  # all data
  else if (length(nsamp) == 1)
    samp <- sample(nrow(dat), nsamp)
  else
    samp <- nsamp  # specific indices specified

  dat <- dat[samp,]


  ## Discretise the dataset
  dat <- discretise.cohort.data(dat, dt=0.2, min.age, max.age, min.time, max.time)
  
  ##  Create aggregated cohort data for likelihood
  aggrdat <- aggregate.cohort.data(dat)
  aggr <- aggrdat$aggr
  exitdat <- aggrdat$exitdat
  cohdat <- aggrdat$cohdat



  min.timeTS <- round(min.time / dt)
  max.timeTS <- round(max.time / dt)
  min.ageTS <- round(min.age / dt)
  max.ageTS <- round(max.age / dt)


  artstart.timeTS <- round(artstart.time / dt)
  natmxstart.timeTS <- round(natmxstart.time / dt)
  
  artstart.tIDX <- as.integer(artstart.timeTS - min.timeTS) + 1L
  natmxstart.tIDX <- as.integer(natmxstart.timeTS - min.timeTS) + 1L


  ## ###################### ##
  ##  Prepare spline model  ##
  ## ###################### ##

  ## Model based on equally spaced knots through domain at interval k.dt

  x.time <- min.timeTS:max.timeTS*dt
  x.age <- min.ageTS:max.ageTS*dt
  STEPS_time=length(x.time)
  STEPS_age=length(x.age)

  ## incidence model
  k.incrate.time <- k.dt*(floor(min.time / k.dt) - 3L):(ceiling(max.time / k.dt) + 3L)
  k.incrate.age <- k.dt*(floor(min.age / k.dt) - 3L):(ceiling(max.age / k.dt) + 3L)

  nk_incrate_time <- length(k.incrate.time)-4L
  nk_incrate_age <- length(k.incrate.age)-4L
  
  x.alltime <- (min.timeTS*2):(max.timeTS*2)*(dt/2)
  x.allage <- (min.ageTS*2):(max.ageTS*2)*(dt/2)
  data_length <- max(length(x.alltime),length(x.allage))
  ytemp <- rep(1,data_length)
  inc_spline <- mgcv::jagam(ytemp ~ s(x.time, k=nk_incrate_time, bs="ps", m=c(2,pen.ord.incrate)) +
                              s(x.age, k=nk_incrate_age, bs="ps", m=c(2,pen.ord.incrate)),
                            data=data.frame(x.time = rep(x.alltime,length.out=data_length),
                                            x.age = rep(x.allage,length.out=data_length), ytemp = ytemp),
                            file=tempfile(), knots = list(x.time = k.incrate.time, x.age = k.incrate.age), 
                            centred=FALSE)[[1]]
  
  X_incrate_time <- inc_spline$X[which(1:length(x.alltime) %% 2!=0),2:(1+nk_incrate_time)]
  Xmid_incrate_time <- inc_spline$X[which(1:length(x.alltime) %% 2==0),2:(1+nk_incrate_time)]
  X_incrate_age <- inc_spline$X[which(1:length(x.allage) %% 2!=0), (2+nk_incrate_time):(1+nk_incrate_time+nk_incrate_age)]
  Xmid_incrate_age <- inc_spline$X[which(1:length(x.allage) %% 2==0),(2+nk_incrate_time):(1+nk_incrate_time+nk_incrate_age)]
  
  # X_incrate_time <- splines::splineDesign(k.incrate.time, x.time, outer.ok=TRUE)
  # Xmid_incrate_time <- splines::splineDesign(k.incrate.time, x.time[-1]-dt/2, outer.ok=TRUE)
  # X_incrate_age <- splines::splineDesign(k.incrate.age, x.age, outer.ok=TRUE)
  # Xmid_incrate_age <- splines::splineDesign(k.incrate.age, x.age[-1]-dt/2)

  D_incrate_time <- inc_spline$smooth[[1]]$D
  D_incrate_age <- inc_spline$smooth[[2]]$D
  
  # D_incrate_time <- diff(diag(nk_incrate_time), differences=pen.ord.incrate)
  # D_incrate_age <- diff(diag(nk_incrate_age), differences=pen.ord.incrate)

  ## Create precision matrix for bivariate incrate smoothing
  Pcar_prec_incrate <- matrix(0, nk_incrate_time*nk_incrate_age, nk_incrate_time*nk_incrate_age)

  ## construct lower triangle with appropriate edges for time (vertical), age (horizontal),
  ## or cohort (diagonal) edges [column-major order].
  if(time.pen)
    diag(Pcar_prec_incrate[-1,]) <- rep(rep(c(0, -1), c(1, nk_incrate_time-1)), nk_incrate_age)[-1]
  if(age.pen)
    diag(Pcar_prec_incrate[-(1:nk_incrate_time),]) <- -1
  if(cohort.pen)
    diag(Pcar_prec_incrate[-(1:(nk_incrate_time+1)),]) <- rep(rep(c(0, -1), c(1, nk_incrate_time-1)), nk_incrate_age-1)[-1]
  
  Pcar_prec_incrate <- Pcar_prec_incrate + t(Pcar_prec_incrate)
  diag(Pcar_prec_incrate) <- -rowSums(Pcar_prec_incrate)

  Pcar_prec_incrate <- expm::"%^%"(Pcar_prec_incrate, pen.ord.incrate)


  ## non-HIV mortality model

  x.natmx <- natmxstart.timeTS:max.timeTS*dt

  k.natmx.time <- k.dt*(floor(natmxstart.time / k.dt) - 3L):(ceiling(max.time / k.dt) + 3L)
  k.natmx.age <- k.incrate.age

  nk_natmx_time <- length(k.natmx.time)-4L
  nk_natmx_age <- length(k.natmx.age)-4L
  
  x.allnat <- (natmxstart.timeTS*2):(max.timeTS*2)*(dt/2)
  x.time.natspline <- c(rep(x.allnat[1], natmxstart.tIDX-1L), x.allnat)
  data_length <- max(length(x.time.natspline),length(x.allage))
  ytemp <- rep(1,data_length)
  natmx_spline <- mgcv::jagam(ytemp ~ s(x.time,k=nk_natmx_time,bs="ps", m=c(2,pen.ord.natmx.time)) +
                                s(x.age,k=nk_natmx_age,bs="ps", m=c(2,pen.ord.natmx.age)),
                              data=data.frame(x.time=rep(x.time.natspline,length.out=data_length), 
                                              x.age = rep(x.allage,length.out=data_length),
                                              ytemp = ytemp),file=tempfile(),
                              knots = list(x.time = k.natmx.time, x.age = k.natmx.age), centred=FALSE)[[1]]

  # X_natmx_time <- splines::splineDesign(k.natmx.time, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx), outer.ok=TRUE)
  # Xmid_natmx_time <- splines::splineDesign(k.natmx.time, c(rep(x.natmx[1], natmxstart.tIDX-1L), x.natmx[-1]-dt/2))
  # X_natmx_age <- splines::splineDesign(k.natmx.age, x.age, outer.ok=TRUE)
  # Xmid_natmx_age <- splines::splineDesign(k.natmx.age, x.age[-1]-dt/2)
  
  X_natmx_time <- natmx_spline$X[which(c(rep(TRUE,natmxstart.tIDX-1L),1:length(x.allnat) %% 2!=0)), 2:(nk_natmx_time+1)]
  Xmid_natmx_time <- natmx_spline$X[which(c(rep(TRUE,natmxstart.tIDX-1L),1:length(x.allnat) %% 2==0)), 2:(nk_natmx_time+1)]
  X_natmx_age <- natmx_spline$X[which(1:length(x.allage) %% 2!=0), (nk_natmx_time+2):(nk_natmx_time+nk_natmx_age+1)]
  Xmid_natmx_age <- natmx_spline$X[which(1:length(x.allage) %% 2==0), (nk_natmx_time+2):(nk_natmx_time+nk_natmx_age+1)]
  
  # D_natmx_time<- diff(diag(nk_natmx_time), differences=pen.ord.natmx.time)
  # D_natmx_age <- diff(diag(nk_natmx_age), differences=pen.ord.natmx.age)

  D_natmx_time <- natmx_spline$smooth[[1]]$D
  D_natmx_age <- natmx_spline$smooth[[2]]$D
  
  ## ART model

  ## x.art <- artstart.timeTS:max.timeTS*dt
  ## X.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splines::splineDesign(k.art, x.art, outer.ok=TRUE))
  ## Xmid.art <- rbind(matrix(0, artstart.tIDX-1L, nk.art), splines::splineDesign(k.art, x.art[-1]-dt/2))
  ## P.art <- diff(diag(nk.art), differences=1)
  
  D_art <- diff(diag(STEPS_time - artstart.tIDX), differences=pen.ord.art)

  ## ##################################### ##
  ##  Calculate HIV survival lookup table  ##
  ## ##################################### ##

  ## log_hivmx_dur_a0 <- log(outer(1:(max.timeTS-min.timeTS)*dt - dt/2, min.ageTS:(max.ageTS-1L)*dt + dt/2, hivsurvhaz, param=NULL))
  ## log_hivsurv_dur_a0 <- log(outer(1:(max.timeTS-min.timeTS)*dt - dt/2, min.ageTS:(max.ageTS-1L)*dt + dt/2, hivsurv, param=NULL))
  ## log_hivmxMID_dur_a0 <- log(-diff(rbind(0, log_hivsurv_dur_a0))/dt)


  ## hivmx_dur_a0 <- exp(log_hivmx_dur_a0)
  ## hivsurv_dur_a0 <- exp(log_hivsurv_dur_a0)
  ## hivmxMID_dur_a0 <- exp(log_hivmxMID_dur_a0)

    
  ## ######################## ##
  ##  Create Stan input data  ##
  ## ######################## ##

  stan.data <- list(dt                    = dt,
                    STEPS_time            = STEPS_time,
                    STEPS_age             = STEPS_age,
                    artstart_tIDX         = artstart.tIDX,
                    ## INDIVIDUAL DATA  ##
                    N                     = nrow(dat),
                    id                    = dat$id,
                    entry_tIDX            = dat$entry.tIDX,
                    exit_tIDX             = dat$exit.tIDX,
                    exposestart_tIDX      = dat$exposestart.tIDX,
                    entry_aIDX            = dat$entry.aIDX,
                    exit_aIDX             = dat$exit.aIDX,
                    exposestart_aIDX      = dat$exposestart.aIDX,
                    expose_DUR            = dat$expose.DUR,
                    hivpos                = dat$hivpos,
                    death                 = dat$death,
                    deathinterv           = dat$deathinterv,
                    deathinterv_DUR       = dat$deathinterv_DUR,
                    ## COHORT DATA ##
                    NCOH                  = nrow(cohdat),
                    coh_cIDX              = cohdat$coh_cIDX,
                    coh_minexpose_tIDX    = cohdat$coh_minexpose_tIDX,
                    coh_maxexpose_tIDX    = cohdat$coh_maxexpose_tIDX,
                    coh_nexit             = cohdat$coh_nexit,
                    coh_ndat              = cohdat$coh_ndat,
                    ## EXIT DATA ##
                    NEXIT                 = nrow(exitdat),
                    exdat_cIDX            = exitdat$exdat_cIDX,
                    exdat_tIDX            = exitdat$exdat_tIDX,
                    exdat_minexpose_tIDX  = exitdat$exdat_minexpose_tIDX,
                    exdat_maxexpose_tIDX  = exitdat$exdat_maxexpose_tIDX,
                    exdat_ndat            = exitdat$exdat_ndat,
                    ## AGGR DATA ##
                    NAGGR                 = nrow(aggr),
                    aggr_cIDX             = aggr$cIDX,
                    aggr_exit_tIDX        = aggr$exit.tIDX,
                    aggr_exposestart_tIDX = aggr$exposestart.tIDX,
                    aggr_exposeend_tIDX   = aggr$exposeend.tIDX,
                    aggr_death            = aggr$death,
                    aggr_deathinterv      = aggr$deathinterv,
                    aggr_deathinterv_DUR  = aggr$deathinterv_DUR,
                    aggr_hivpos           = aggr$hivpos,
                    aggr_nrepl            = aggr$nrepl,
                    ## MODEL PARAMETERS ##
                    x_time                = x.time,
                    x_age                 = x.age,
                    ## x_art                 = x.art,
                    x_natmx               = x.natmx,
                    ## incidence model
                    nk_incrate_time       = nk_incrate_time,
                    nk_incrate_age        = nk_incrate_age,
                    k_incrate_time        = k.incrate.time,
                    k_incrate_age         = k.incrate.age,
                    X_incrate_time        = X_incrate_time,
                    Xmid_incrate_time     = Xmid_incrate_time,
                    X_incrate_age         = X_incrate_age,
                    Xmid_incrate_age      = Xmid_incrate_age,
                    D_incrate_time        = D_incrate_time,
                    D_incrate_age         = D_incrate_age,
                    pen_ord_incrate       = pen.ord.incrate,
                    Pcar_prec_incrate     = Pcar_prec_incrate,
                    ## non-HIV mortality model
                    nk_natmx_time         = nk_natmx_time,
                    nk_natmx_age          = nk_natmx_age,
                    k_natmx_time          = k.natmx.time,
                    k_natmx_age           = k.natmx.age,
                    X_natmx_time          = X_natmx_time,
                    Xmid_natmx_time       = Xmid_natmx_time,
                    X_natmx_age           = X_natmx_age,
                    Xmid_natmx_age        = Xmid_natmx_age,
                    pen_ord_natmx_time    = pen.ord.natmx.time,
                    pen_ord_natmx_age     = pen.ord.natmx.age,
                    D_natmx_time          = D_natmx_time,
                    D_natmx_age           = D_natmx_age,
                    fixcoef_natmx_time    = as.integer(nk_natmx_time/2),
                    fixcoef_natmx_age     = as.integer(nk_natmx_age/2),
                    ## ART model
                    ## nk_art                = nk.art,
                    ## k_art                 = k.art,
                    ## X_art                 = X.art,
                    ## Xmid_art              = Xmid.art,
                    pen_ord_art           = pen.ord.art,
                    D_art                 = D_art,
                    ## Variance of priors
                    var_incrate_time_age  = var_incrate_time_age,
                    var_incrate_time      = var_incrate_time,
                    var_incrate_age       = var_incrate_age,
                    var_natmx_time        = var_natmx_time,
                    var_natmx_age         = var_natmx_age,
                    var_art               = var_art,
                    ## HIV survival model
                    hivsurv_shape         = hivsurv.shape)
                    ## log_hivmx_dur_a0      = log_hivmx_dur_a0,
                    ## log_hivsurv_dur_a0    = log_hivsurv_dur_a0,
                    ## log_hivmxMID_dur_a0   = log_hivmxMID_dur_a0,
                    ## hivmx_dur_a0          = hivmx_dur_a0,
                    ## hivsurv_dur_a0        = hivsurv_dur_a0,
                    ## hivmxMID_dur_a0       = hivmxMID_dur_a0)
}


#### Calculate likelihood

Rcalc_ll_coh <- function(stand, modpred){
  sum(calc_ll_coh(stand$entry_tIDX, stand$entry_aIDX, stand$exit_tIDX, stand$exit_aIDX,
                  stand$exposestart_tIDX, stand$exposestart_aIDX, stand$expose_DUR,
                  stand$death, stand$deathinterv, stand$deathinterv_DUR, stand$hivpos,
                  modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                  modpred$hivsurv_dur_a0, modpred$hivmx_dur_a0, modpred$hivmxMID_dur_a0,
                  modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                  modpred$natsurv_time_age, modpred$natmx_time_age,
                  stand$dt))
}

Rcalc_ll_coh_aggr <- function(stand, modpred){
  calc_ll_coh_aggr(stand$aggr_cIDX, stand$aggr_exit_tIDX, stand$aggr_exposestart_tIDX, stand$aggr_exposeend_tIDX,
                   stand$aggr_death, stand$aggr_deathinterv, stand$aggr_deathinterv_DUR, stand$aggr_hivpos, stand$aggr_nrepl,
                   modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                   modpred$hivsurv_dur_a0, modpred$hivmx_dur_a0, modpred$hivmxMID_dur_a0,
                   modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                   modpred$natsurv_time_age, modpred$natmx_time_age,
                   stand$dt)
}

Rcalc_ll_cohexit <- function(stand, modpred){
  calc_ll_cohexit(stand$coh_cIDX, stand$coh_minexpose_tIDX, stand$coh_maxexpose_tIDX, stand$coh_nexit,
                  stand$exdat_tIDX, stand$exdat_minexpose_tIDX, stand$exdat_maxexpose_tIDX, stand$exdat_ndat,
                  stand$aggr_exposestart_tIDX, stand$aggr_exposeend_tIDX,
                  stand$aggr_death, stand$aggr_deathinterv, stand$aggr_deathinterv_DUR, stand$aggr_hivpos, stand$aggr_nrepl,
                  modpred$cumavoid_time_age, modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                  modpred$hivsurv_dur_a0, modpred$hivmx_dur_a0, modpred$hivmxMID_dur_a0,
                  modpred$artrr, modpred$artrr_MID, stand$artstart_tIDX,
                  modpred$natsurv_time_age, modpred$natmx_time_age,
                  stand$dt)
}
