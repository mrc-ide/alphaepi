##################################
####                          ####
####    Analysis functions    ####
####                          ####
##################################

invlogit <- function(x) exp(x)/(1+exp(x))

convert.stan.params <- function(par, stand){

  nk_time <- stand$nk_incrate_time
  nk_age <- stand$nk_incrate_age
  nk_natmx_time <- stand$nk_natmx_time
  nk_art <- stand$nk_art
  
  coef_incrate_time_age.idx <- 1:(nk_time*nk_age)
  coef_natmx_time.idx <- nk_time*nk_age + 1:nk_natmx
  coef_natmx_age.idx <- nk_time*nk_age + nk_natmx + 1:(nk_age-1L)
  coef_art.idx <- nk_time*nk_age + nk_natmx + (nk_age-1L) + 1:nk_art

  coef_incrate_time_age <- matrix(par[coef_incrate_time_age.idx], nk_time, nk_age)

  coef_natmx_time <- par[coef_natmx_time.idx]
  
  coef_natmx_age <- rep(0, nk_age)
  coef_natmx_age[-stand$fixcoef_age_idx] <- par[coef_natmx_age.idx]
    
  coef_natmx_time_age <- outer(coef_natmx_time, coef_natmx_age, "+")
  
  coef_art <- par[coef_art.idx]

  return(list(coef_incrate_time_age = coef_incrate_time_age,
              coef_natmx_time_age = coef_natmx_time_age,
              coef_natmx_time = coef_natmx_time,
              coef_natmx_age = coef_natmx_age,
              coef_art = coef_art))
}

create.param.list <- function(stanfit){
  param <- rstan::extract(stanfit)
  param <- lapply(seq_along(param$lp__), function(ii) list(coef_incrate_time_age   = param$coef_incrate_time_age[ii,,],
                                                           coef_natmx_time         = param$coef_natmx_time[ii,],
                                                           coef_natmx_age          = param$coef_natmx_age[ii,],
                                                           coef_natmx_time_age     = outer(param$coef_natmx_time[ii,], param$coef_natmx_age[ii,], "+"),
                                                           dt_log_artrr            = param$dt_log_artrr[ii,],
                                                           sigma2_incrate_time_age = param$sigma2_incrate_time_age[ii],
                                                           sigma2_natmx_time       = param$sigma2_natmx_time[ii],
                                                           sigma2_natmx_age        = param$sigma2_natmx_age[ii],
                                                           sigma2_art              = param$sigma2_art[ii],
                                                           hivsurv_shape           = param$hivsurv_shape[ii],
                                                           hivsurv_scale_b0        = param$hivsurv_scale_b0[ii],
                                                           hivsurv_scale_b1        = param$hivsurv_scale_b1[ii]))
  return(param)
}

create.modpred <- function(param, stand){

  incrateMID_time_age <- exp(stand$Xmid_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age))
  cumavoid_time_age <- exp(-stand$dt*diagCumSum(incrateMID_time_age))
  cumavoidMID_time_age <- cumavoid_time_age[1:(stand$STEPS_time-1), 1:(stand$STEPS_age-1)]*exp(-stand$dt/2*incrateMID_time_age)

  natmx_time_age <- exp(stand$X_natmx_time %*% param$coef_natmx_time_age %*% t(stand$X_natmx_age))
  natsurv_time_age <- exp(-stand$dt*diagCumSum(exp(stand$Xmid_natmx_time %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age))))

  log_artrr <- stand$dt*cumsum(param$dt_log_artrr);
  artrr <- rep(1.0, stand$STEPS_time)
  artrr_MID <- rep(1.0, stand$STEPS_time-1L)
  artrr[(stand$artstart_tIDX+1L):stand$STEPS_time] <- exp(log_artrr)
  artrr_MID[stand$artstart_tIDX:(stand$STEPS_time-1L)] <- exp(log_artrr - stand$dt/2*param$dt_log_artrr)

  hivsurv_scale_a0 <- exp(param$hivsurv_scale_b0 + param$hivsurv_scale_b1 * ((stand$x_age[-1] - stand$dt/2) - 30) / 10) # !!! HARD-CODED DESIGN MATRIX
  hivmx_dur_a0 <- create_hivmx_dur_a0(param$hivsurv_shape, hivsurv_scale_a0, stand$STEPS_time-1L, stand$dt)
  log_hivsurv_dur_a0 <- create_log_hivsurv_dur_a0(param$hivsurv_shape, hivsurv_scale_a0, stand$STEPS_time-1L, stand$dt)
  hivmxMID_dur_a0 <- diff_hivmxMID_dur_a0(log_hivsurv_dur_a0, stand$dt)
  hivsurv_dur_a0 <- exp(log_hivsurv_dur_a0)

  return(list(incrateMID_time_age = incrateMID_time_age,
              cumavoid_time_age = cumavoid_time_age,
              cumavoidMID_time_age = cumavoidMID_time_age,
              natmx_time_age = natmx_time_age,
              natsurv_time_age = natsurv_time_age,
              artrr = artrr,
              artrr_MID = artrr_MID,
              hivmx_dur_a0 = hivmx_dur_a0,
              hivsurv_dur_a0 = hivsurv_dur_a0,
              hivmxMID_dur_a0 = hivmxMID_dur_a0))
}


cumincid.period <- function(param, stand){
  log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age)

  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  return(setNames(cumincid.period, stand$x_time))
}

cumnatmort.period <- function(param, stand){
  log_natmx_time_age <- stand$X_natmx_time %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age)

  cumnatmort.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_natmx_time_age)))
  return(setNames(cumnatmort.period, stand$x_time))
}

prev <- function(tidx, aidx, modpred, stand){
  exposeDUR <- min(tidx, aidx)-1L
  phivp <- calc_phivp(tidx, aidx, tidx-exposeDUR, aidx-exposeDUR, exposeDUR, 0,
                      modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                      modpred$hivsurv_dur_a0, modpred$hivmx_dur_a0, modpred$hivmxMID_dur_a0, modpred$artrr, modpred$artrr_MID,
                      stand$artstart_tIDX, modpred$natsurv_time_age, modpred$natmx_time_age, stand$dt)
  phivn <- calc_phivn(tidx, aidx, 0, 0, modpred$cumavoid_time_age, modpred$natsurv_time_age, modpred$natmx_time_age)
  return(phivp/(phivp+phivn))
}


Rcreate_phivp_mat <- function(modpred, stand, art=TRUE){
  if(!art)
    stand$artstart_tIDX <- stand$STEPS_time
  create_phivp_mat(modpred$cumavoidMID_time_age, modpred$incrateMID_time_age, modpred$hivsurv_dur_a0, modpred$hivmxMID_dur_a0, modpred$artrr_MID, stand$artstart_tIDX, modpred$natsurv_time_age, stand$dt)
}
Rcreate_phivn_mat <- function(modpred){ return(modpred$cumavoid_time_age * modpred$natsurv_time_age)}

calc.psurv <- function(param, stand){
  modpred <- create.modpred(param, stand)
  phivn <- Rcreate_phivn_mat(modpred)
  phivp <- Rcreate_phivp_mat(modpred, stand)
  phivp.noart <- Rcreate_phivp_mat(modpred, stand, art=FALSE)
  psurv <- phivn+phivp
  psurv.noart <- phivn+phivp.noart
  psurv.nohiv <- modpred$natsurv_time_age
  list(psurv=psurv, psurv.noart=psurv.noart, psurv.nohiv=psurv.nohiv, prev=phivp/psurv, prev.noart=phivp.noart/psurv.noart)
}

calc.le <- function(psurv, dt){
  psurv.last <- psurv[-nrow(psurv), -ncol(psurv)]
  psurv.curr <- psurv[-1,-1]
  dt*colSums(apply(1.0 - (psurv.last-psurv.curr)/psurv.last, 1, cumprod))
}

calc.cumincid <- function(param, stand, years=diff(range(mod$stand$x_age))){
  aidx <- seq_len(years/stand$dt)
  log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  setNames(cumincid.period, stand$x_time)
}

calc.incid <- function(param, stand){ exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age)) }

calc.prev <- function(param, stand){
  modpred <- create.modpred(param, stand)
  phivn <- Rcreate_phivn_mat(modpred)
  phivp <- Rcreate_phivp_mat(modpred, stand)
  return(phivp/(phivn+phivp))
}



###############################
####  Aggregate functions  ####
###############################

library(parallel)

calc.45q15 <- function(psurv, dt){
  aidx <- seq_len(45/dt) # assumes psurv starts at age 15
  psurv.last <- psurv[-nrow(psurv), -ncol(psurv)]
  psurv.curr <- psurv[-1,-1]
  qx <- (psurv.last[,aidx]-psurv.curr[,aidx])/psurv.last[,aidx]
  1.0 - exp(rowSums(log(1.0-qx)))
}

add.le <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(le       = sapply(lapply(psurvobj, "[[", "psurv"), calc.le, mod$stand$dt),
       le.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.le, mod$stand$dt),
       le.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.le, mod$stand$dt))
}

add.45q15 <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(q4515       = sapply(lapply(psurvobj, "[[", "psurv"), calc.45q15, mod$stand$dt),
       q4515.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.45q15, mod$stand$dt),
       q4515.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.45q15, mod$stand$dt))
}

add.mx <- function(mod){
  param <- create.param.list(mod$fit)
  psurvobj <- mclapply(param, calc.psurv, mod$stand)
  list(le       = sapply(lapply(psurvobj, "[[", "psurv"), calc.le, mod$stand$dt),
       le.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.le, mod$stand$dt),
       le.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.le, mod$stand$dt),
       q4515       = sapply(lapply(psurvobj, "[[", "psurv"), calc.45q15, mod$stand$dt),
       q4515.noart = sapply(lapply(psurvobj, "[[", "psurv.noart"), calc.45q15, mod$stand$dt),
       q4515.nohiv = sapply(lapply(psurvobj, "[[", "psurv.nohiv"), calc.45q15, mod$stand$dt))
}


who.standard.pop <- setNames(c(179177, 178022, 177011, 176133, 175366, 174700, 174144, 173666, 173277, 172933, 172633, 172333, 171999, 171600, 171155, 170699, 170188, 169522, 168666, 167655, 166566, 165466, 164344, 163233, 162122, 160944, 159711, 158511, 157355, 156211, 155044, 153766, 152333, 150700, 148900, 147022, 145100, 143077, 140911, 138655, 136299, 133911, 131633, 129499, 127455, 125388, 123200, 120911, 118455, 115855, 113200, 110488, 107577, 104433, 101133, 9774, 9435, 9095, 8757, 8423, 8087, 7750, 7425, 7113, 6812, 6517, 6221, 5923, 5618, 5311, 5005, 4706, 4412, 4125, 3844, 3569, 3300, 3035, 2773, 2518, 2271, 2033, 1807, 1594, 1203, 1025, 863, 718, 589, 475, 374, 287, 212, 152, 106, 76, 62, 66, 90, 50), 0:99)



add.incprev <- function(mod){
  param <- create.param.list(mod$fit)

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  system.time(prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param))))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  mod$prev <- prev

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  mod$incid <- incid

  ## cumincid
  mod$cumincid <- sapply(param, calc.cumincid, mod$stand, years=45) # incidence between age 15 and 60

  ## WHO prev
  aidx <- 1:(35/mod$stand$dt+1L)
  age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidx])$y  # prevalence age 15 to 50
  age.dist <- age.dist/sum(age.dist)
  mod$who15to49prev <- apply(sweep(mod$prev[, aidx,], 2, age.dist, "*"), c(1,3), sum)

  return(mod)
}


##########################
####  Plot functions  ####
##########################

cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), border=NA, ...)

transp <- function(col, alpha=0.5)
  return(apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha)))

plot.timebysex <- function(m.post, f.post, main, ylim=c(0,1), xlim=c(1990, 2015), ylab=""){
  par(tcl=-0.25, las=1, mgp=c(2, 0.5, 0), mar=c(2.1, 3.1, 2.1, 1.1))
  xx <- as.numeric(rownames(m.post))
  matplot(xx, cbind(apply(m.post, 1, median), apply(f.post, 1, median)),
        type="l", ylim=ylim, xlim=xlim, col=c("royalblue", "deeppink"), lty=1, lwd=1.5,
        main=main, xlab="", ylab=ylab)
  cred.region(xx, apply(m.post, 1, quantile, c(0.025, 0.975)), col=transp("royalblue"))
  cred.region(xx, apply(f.post, 1, quantile, c(0.025, 0.975)), col=transp("deeppink"))
  legend("topleft", legend=c("men", "women"), col=c("royalblue", "deeppink"), lty=1.5, bg="white", cex=0.8)
}

plot.ci <- function(x, y, se=NULL, lower=NULL, upper=NULL, col=1, pch=20, cex=0.8, lty=1, lwd=1){
  points(x, y, pch=pch, col=col, cex=cex)
  segments(x0=x, y0=if(is.null(lower)) y-qnorm(0.975)*se else lower, y1=if(is.null(upper)) y+qnorm(0.975)*se else upper, col=col, lty=lty, lwd=lwd)
}


plot.le <- function(le, stand, main, ylim=c(35, 65), xlim=c(1980, 2015)){
  par(las=1, tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(2.1, 3.1, 2.1, 0.5))
  ##
  xx <- stand$x_time[-1]
  idx.le <- 2:stand$STEPS_time - 1
  idx.noart <- stand$artstart_tIDX:stand$STEPS_time-1
  idx.nohiv <- which(xx %in% stand$x_natmx)
  ##
  plot(xx[idx.le], rowMeans(le$le[idx.le,]), type='n', lty=1, ylim=ylim, xlim=xlim,
       ylab=expression(e[15]), xlab="", main=main)
  ##
  cred.region(xx[idx.noart], apply(le$le.noart[idx.noart,], 1, quantile, c(0.025, 0.975)), col=transp(2))
  cred.region(xx[idx.nohiv], apply(le$le.nohiv[idx.nohiv,], 1, quantile, c(0.025, 0.975)), col=transp(3))
  cred.region(xx[idx.le], apply(le$le[idx.le,], 1, quantile, c(0.025, 0.975)), col=transp(1))
  ##
  lines(xx[idx.noart], rowMeans(le$le.noart[idx.noart,]), col=2, lwd=1.5)
  lines(xx[idx.nohiv], rowMeans(le$le.nohiv[idx.nohiv,]), col=3, lwd=1.5)
  lines(xx[idx.le], rowMeans(le$le[idx.le,]), col=1, lwd=1.5)
  ##
  return(invisible(NULL))
}

plot.45q15 <- function(q4515, stand, main, ylim=c(0, 0.8), xlim=c(1980, 2015)){
  par(las=1, tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(2.1, 3.1, 2.1, 0.5))
  ##
  xx <- stand$x_time[-1]
  idx.est <- 2:stand$STEPS_time - 1
  idx.noart <- stand$artstart_tIDX:stand$STEPS_time-1
  idx.nohiv <- which(xx %in% stand$x_natmx)
  ##
  plot(xx[idx.est], rowMeans(q4515$q4515[idx.est,]), type='n', lty=1, ylim=ylim, xlim=xlim,
       ylab=expression(""[45]~q[15]), xlab="", main=main)
  ##
  cred.region(xx[idx.noart], apply(q4515$q4515.noart[idx.noart,], 1, quantile, c(0.025, 0.975)), col=transp(2))
  cred.region(xx[idx.nohiv], apply(q4515$q4515.nohiv[idx.nohiv,], 1, quantile, c(0.025, 0.975)), col=transp(3))
  cred.region(xx[idx.est], apply(q4515$q4515[idx.est,], 1, quantile, c(0.025, 0.975)), col=transp(1))
  ##
  lines(xx[idx.noart], rowMeans(q4515$q4515.noart[idx.noart,]), col=2, lwd=1.5)
  lines(xx[idx.nohiv], rowMeans(q4515$q4515.nohiv[idx.nohiv,]), col=3, lwd=1.5)
  lines(xx[idx.est], rowMeans(q4515$q4515[idx.est,]), col=1, lwd=1.5)
  ##
  return(invisible(NULL))
}

