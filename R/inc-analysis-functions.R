##################################
####                          ####
####    Analysis functions    ####
####                          ####
##################################


prev <- function(tidx, aidx, modpred, stand){
  exposeDUR <- min(tidx, aidx)-1L
  phivp <- calc_phivp(tidx, aidx, tidx-exposeDUR, aidx-exposeDUR, exposeDUR, 0,
                      modpred$cumavoidMID_time_age, modpred$incrateMID_time_age,
                      modpred$hivsurv_dur_a0, modpred$hivmx_dur_a0, modpred$hivmxMID_dur_a0, modpred$artrr, modpred$artrr_MID,
                      stand$artstart_tIDX, modpred$natsurv_time_age, modpred$natmx_time_age, stand$dt)
  phivn <- calc_phivn(tidx, aidx, 0, 0, modpred$cumavoid_time_age, modpred$natsurv_time_age, modpred$natmx_time_age)
  return(phivp/(phivp+phivn))
}


calc.cumincid <- function(param, stand, aidx = 1:stand$STEPS_age,incid){
  # aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  # log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  log_incrate_time_age <- log(incid[,aidx])
  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  setNames(cumincid.period, stand$x_time)
}

calc.period.cumincid.proj <- function(incid, stand, ages){
  log_incrate_time_age <- log(incid[,colnames(incid)>=ages[1] & colnames(incid)<ages[2]])
  cumincid.period <- 1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age)))
  setNames(cumincid.period, rownames(incid))
}

# cumulative incidence at each age for cohorts
calc.cumincid.coh.mult <- function(param, stand, years=diff(range(mod$stand$x_age))){
  # aidx <- seq_len(years/stand$dt)
  aidx <- which(stand$x_age==years[1]):which(stand$x_age==years[2])
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  exp_incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  exp_incrate_time_age <- exp_incrate_time_age[,aidx]
  # log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  # exp_incrate_time_age <- exp(log_incrate_time_age)
  rates_for_calc <- sapply(split(exp_incrate_time_age, col(exp_incrate_time_age) -
                                   row(exp_incrate_time_age)), cumsum)
  cumincid.coh <- lapply(rates_for_calc,function(x) {1.0 - exp(-stand$dt * x)} )
  yearbirth <- c((stand$x_time-15), stand$x_time[1] - stand$x_age[aidx[-1]])
  yearbirth <- sort(yearbirth,decreasing=TRUE)
  setNames(cumincid.coh, yearbirth)
}

calc.incid.proj <- function(param, stand, years=diff(range(mod$stand$x_age)),karman=FALSE) {
  # aidx <- seq_len(years/stand$dt)
  if(karman==TRUE) {
    tidx <- which(stand$x_time<=2012)
  } else {
    tidx <- 1:length(stand$x_time)
  }
  nreps <- length(seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  newtime <- c(stand$x_time[tidx],seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  aidx <- which(stand$x_age==years[1]):which(stand$x_age==years[2])
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  exp_incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  exp_incrate_time_age <- exp_incrate_time_age[tidx,aidx]
  exp_incrate_time_age <- rbind(exp_incrate_time_age,
                                matrix(rep(exp_incrate_time_age[nrow(exp_incrate_time_age),],nreps),nrow=nreps,byrow=TRUE))
  return(exp_incrate_time_age)
}

calc.cumincid.coh.mult.proj <- function(param, stand, years=diff(range(mod$stand$x_age)),karman=FALSE){
  # aidx <- seq_len(years/stand$dt)
  if(karman==TRUE) {
    tidx <- which(stand$x_time<=2012)
  } else {
    tidx <- 1:length(stand$x_time)
  }
  nreps <- length(seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  newtime <- c(stand$x_time[tidx],seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  aidx <- which(stand$x_age==years[1]):which(stand$x_age==years[2])
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  exp_incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  exp_incrate_time_age <- exp_incrate_time_age[tidx,aidx]
  exp_incrate_time_age <- rbind(exp_incrate_time_age,
                                matrix(rep(exp_incrate_time_age[nrow(exp_incrate_time_age),],nreps),nrow=nreps,byrow=TRUE))
  # log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  # exp_incrate_time_age <- exp(log_incrate_time_age)
  rates_for_calc <- sapply(split(exp_incrate_time_age, col(exp_incrate_time_age) -
                                   row(exp_incrate_time_age)), cumsum)
  cumincid.coh <- lapply(rates_for_calc,function(x) {1.0 - exp(-stand$dt * x)} )
  yearbirth <- c((newtime-15),newtime[1] - stand$x_age[aidx[-1]])
  yearbirth <- sort(yearbirth,decreasing=TRUE)
  setNames(cumincid.coh, yearbirth)
}

calc.incid.proj2 <- function(param, stand, years=diff(range(mod$stand$x_age)),karman=FALSE,reldec,
                             meetyear) {
  # aidx <- seq_len(years/stand$dt)
  if(karman==TRUE) {
    tidx <- which(stand$x_time<=2012)
  } else {
    tidx <- 1:length(stand$x_time)
  }
  nreps <- length(seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  newtime <- c(stand$x_time[tidx],seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  timediff <- seq((max(stand$x_time[tidx]) + 0.2),max(newtime),0.2) - max(stand$x_time[tidx])
  if(length(reldec)==1) {
    alldec <- reldec^timediff
  } else {
    steps <- length(seq(max(stand$x_time[tidx])+0.2,meetyear,0.2))
    converge_dec <- exp(seq(log(reldec[1]),(log(reldec[2])+log(reldec[1]))/2,length=steps))
    steps_same <- length(seq(meetyear+0.2,1995+years[2],0.2))
    reldec_all <- c(converge_dec,rep(exp((log(reldec[2])+log(reldec[1]))/2),steps_same))
    alldec <- cumprod(reldec_all^.2)
  }
  aidx <- which(stand$x_age==years[1]):which(stand$x_age==years[2])
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  exp_incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  exp_incrate_time_age <- exp_incrate_time_age[tidx,aidx]
  extra_decreasing <- sweep(matrix(rep(exp_incrate_time_age[nrow(exp_incrate_time_age),],nreps),nrow=nreps,byrow=TRUE),
                            MARGIN=1,alldec,`*`)
  exp_incrate_time_age <- rbind(exp_incrate_time_age,
                                extra_decreasing)
  return(exp_incrate_time_age)
}

calc.cumincid.coh.mult.proj2 <- function(param, stand, years=diff(range(mod$stand$x_age)),karman=FALSE,reldec,
                                         meetyear){
  # aidx <- seq_len(years/stand$dt)
  if(karman==TRUE) {
    tidx <- which(stand$x_time<=2012)
  } else {
    tidx <- 1:length(stand$x_time)
  }
  nreps <- length(seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  newtime <- c(stand$x_time[tidx],seq(max(stand$x_time[tidx])+.2,1995+years[2],.2))
  timediff <- seq((max(stand$x_time[tidx]) + 0.2),max(newtime),0.2) - max(stand$x_time[tidx])
  if(length(reldec)==1) {
    alldec <- reldec^timediff
  } else {
    steps <- length(seq(max(stand$x_time[tidx])+0.2,meetyear,0.2))
    converge_dec <- exp(seq(log(reldec[1]),(log(reldec[2])+log(reldec[1]))/2,length=steps))
    steps_same <- length(seq(meetyear+0.2,1995+years[2],0.2))
    reldec_all <- c(converge_dec,rep(exp((log(reldec[2])+log(reldec[1]))/2),steps_same))
    alldec <- cumprod(reldec_all^.2)
  }
  aidx <- which(stand$x_age==years[1]):which(stand$x_age==years[2])
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  exp_incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  exp_incrate_time_age <- exp_incrate_time_age[tidx,aidx]
  extra_decreasing <- sweep(matrix(rep(exp_incrate_time_age[nrow(exp_incrate_time_age),],nreps),nrow=nreps,byrow=TRUE),
                            MARGIN=1,alldec,`*`)
  exp_incrate_time_age <- rbind(exp_incrate_time_age,
                                extra_decreasing)
  # log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  # exp_incrate_time_age <- exp(log_incrate_time_age)
  rates_for_calc <- sapply(split(exp_incrate_time_age, col(exp_incrate_time_age) -
                                   row(exp_incrate_time_age)), cumsum)
  cumincid.coh <- lapply(rates_for_calc,function(x) {1.0 - exp(-stand$dt * x)} )
  yearbirth <- c((newtime-15),newtime[1] - stand$x_age[aidx[-1]])
  yearbirth <- sort(yearbirth,decreasing=TRUE)
  setNames(cumincid.coh, yearbirth)
}

calc.avgageinf <- function(param, stand, mod, aidx, years = range(mod$stand$x_time), w, incid) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  midw <- (w[,-ncol(w)] + w[,-1])/2

  # log_incrate_time_age <- stand$X_incrate_time[tidx,] %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  #pinfect <- exp(log_incrate_time_age) * exp(-stand$dt * exp(log_incrate_time_age))
  # pinfect <- exp(log_incrate_time_age)
  pinfect <- incid[tidx,aidx]
  avgageinf <- rowSums((pinfect * midw) %*% stand$x_age[aidx])/rowSums(pinfect * midw)
  setNames(avgageinf, stand$x_time[tidx])
}

calc.numnewinf <- function(param, stand, mod, aidx,
                           years = range(mod$stand$x_time), susc, incid) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  midsusc <- (susc[,-ncol(susc)] + susc[,-1])/2

  # log_incrate_time_age <- stand$X_incrate_time[tidx,] %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  # pinfect <- exp(log_incrate_time_age)
  pinfect <- incid[tidx,aidx]
  numnewinf <- pinfect * midsusc
  return(numnewinf)
}

calc.props_by_agrgr <- function(param, stand, mod, ages = range(mod$stand$x_age),
                                years = range(mod$stand$x_time), w, byyears = 10) {
  startind <- which(mod$stand$x_age==ages[1])
  endind <- which(mod$stand$x_age==ages[2])
  aidx <- startind:(endind-1)
  # aidx <- seq_len((ages[2]-ages[1])/stand$dt)
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  midw <- (w[,-ncol(w)] + w[,-1])/2

  log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age)
  extra <- matrix(log(1e-15),nrow=nrow(incrate_time_age),ncol=24)
  youngs <- stand$X_incrate_time %*% param$coef_incrate_time_young
  log_incrate_time_age <- cbind(youngs,extra,log_incrate_time_age)
  log_incrate_time_age <- log_incrate_time_age[tidx,aidx]
  # log_incrate_time_age <- stand$X_incrate_time[tidx,] %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  pinfect <- exp(log_incrate_time_age) * exp(-stand$dt * exp(log_incrate_time_age))

  pinfect_w <- pinfect * midw / rowSums(midw)
  prop_infec <- pinfect_w
  # prop_infec <- pinfect_w / rowSums(pinfect_w) # should i do this now or when plotting??
  prop_infec_gr <- sapply(seq(1,ncol(prop_infec),byyears/stand$dt), function(i) {
    rowSums(prop_infec[,c(i:(i+byyears/stand$dt-1))], na.rm=T)})
  agenames <- c(paste0(seq(ages[1],ages[2]-byyears,byyears),"-",seq(ages[1]+byyears-1,ages[2]-1,byyears)))
  colnames(prop_infec_gr) <- agenames
  return(prop_infec_gr)
}

calc.ageinfunw <- function(param, stand, mod, aidx, years = range(mod$stand$x_time), incid) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  pinfect <- incid[tidx,aidx]
  # log_incrate_time_age <- stand$X_incrate_time[tidx,] %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age[aidx,])
  # pinfect <- exp(log_incrate_time_age)
  # pinfect <- exp(log_incrate_time_age) * exp(-stand$dt * exp(log_incrate_time_age))
  avgageinf <- rowSums((pinfect) %*% stand$x_age[aidx])/rowSums(pinfect)
  setNames(avgageinf, stand$x_time[tidx])
}

calc.incid <- function(param, stand){
  # exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  incrate_time_age <- exp(stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$X_incrate_age))
  extra <- matrix(1e-15,nrow=nrow(incrate_time_age),ncol=24)
  youngs <- exp(stand$X_incrate_time %*% param$coef_incrate_time_young)
  incrate_time_age <- cbind(youngs,extra,incrate_time_age)
  return(incrate_time_age)
}

calc.prev <- function(param, stand){
  modpred <- create.modpred(param, stand)
  phivn <- Rcreate_phivn_mat(modpred)
  phivp <- Rcreate_phivp_mat(modpred, stand)
  return(phivp/(phivn+phivp))
}


###############################
####  Aggregate functions  ####
###############################

# library(parallel)

add.medinc <- function(mod, ages=c(15,60)) {

  param <- create.param.list(mod$fit)
  stand <- mod$stand
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  medinc <- apply(incid,1:2,median)
  out$incid <- medinc

  return(out)
}


add.incprev <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$incid <- incid

  ## cumincid
  incidlist <- lapply(seq(dim(incid)[3]), function(x) incid[ , , x])
  out$cumincid <- mapply(function(x,y) calc.cumincid(x,mod$stand,aidx=aidx,incid=y),
                         param, incidlist)
  rownames(out$cumincid) <- mod$stand$x_time

  ## WHO prev
  aidxwho <- which(mod$stand$x_age==15):which(mod$stand$x_age==50)
  age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidxwho])$y  # prevalence age 15 to 50
  age.dist <- age.dist/sum(age.dist)
  out$who15to49prev <- apply(sweep(out$prev[, aidxwho,], 2, age.dist, "*"), c(1,3), sum)
  rownames(out$who15to49prev) <- mod$stand$x_time

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  out$pop <- pop
  # plus pop full ages
  popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
                                                      years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  out$susc_prop <- susc_prop

  ## Prevalence weighted by UNPop national
  prev_pop <- apply(prev, 3, function(x) add_prev_pop(mod, x, pop_prop, aidx, years))
  out$prev_pop <- prev_pop

  ## Incidence weighted by UNpop national susceptibles
  if(years[2]>max(mod$stand$x_time)){
    topyear <- max(mod$stand$x_time)
  } else { topyear<- years[2]}
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==topyear)
  incid_w <- matrix(nrow=length(times),ncol = dim(prev)[[3]])
  for(i in 1:1000) incid_w[,i] <- rowSums(susc_prop[,,i]*incid[tidx,aidx,i])/rowSums(susc_prop[,,i])
  rownames(incid_w) <- times
  out$incid_w <- incid_w

  ## Mean of period age-specific incidence
  ageinfunw <- mapply(function(x,y) calc.ageinfunw(x,stand,mod,aidx,years,y), param, incidlist)
  out$ageinfunw <- ageinfunw

  ## Average age of infection (among susceptible population)
  susclist <- lapply(seq(dim(susc_prop)[3]), function(x) susc_prop[ , , x])# make a list to go in mapply - clean this!
  avgageinf <- mapply(function(x,y,z){calc.avgageinf(x,mod$stand,mod,aidx[-length(aidx)],years,y,z)}, param, susclist, incidlist)
  out$avgageinf <- avgageinf

  ## Number of new infections by age/time
  nsusc2 <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),years))
  dim(nsusc2) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(nsusc2) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  nsusc2list <- lapply(seq(dim(nsusc2)[3]), function(x) nsusc2[ , , x])# make a list to go in mapply - clean this!
  numnewinf <- mapply(function(x,y,z){calc.numnewinf(x,mod$stand,mod,1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),years,y,z)},
                      param, nsusc2list, incidlist)
  dim(numnewinf) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  dimnames(numnewinf) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  out$numnewinf <- numnewinf

  ## Proportion of new infections by age
  prop_newinf <- apply(numnewinf[,aidx[-length(aidx)],],3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(prop_newinf) <- c(length(times),length(aidx)-1,dim(prev)[3])
  dimnames(prop_newinf) <- list(times,mod$stand$x_age[aidx[-length(aidx)]],NULL)
  out$prop_newinf <- prop_newinf

  prop_newinf5 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,5/stand$dt), function(i) {
    rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/5,dim(prev)[3])
  dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  out$prop_newinf5 <- prop_newinf5

  # Add 10 yr age groups
  # prop_newinf10 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,10/stand$dt), function(i) {
  #   rowSums(x[,c(i:(i+10/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/10,dim(prev)[3])
  # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)

  prop_newinfplus <- apply(numnewinf,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(prop_newinfplus) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  dimnames(prop_newinfplus) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  out$prop_newinfplus <- prop_newinfplus

  prop_newinfplus5 <- apply(prop_newinfplus,3,function(x) sapply(seq(1,((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt-1)
                                                                     ,5/stand$dt), function(i) {
                                                                       rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  dim(prop_newinfplus5) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)*.2/5,dim(prev)[3])
  dimnames(prop_newinfplus5) <- list(times,c(paste0(seq(min(stand$x_age),min(max(mod$stand$x_age),65)-5,5),"-",
                                                    seq(min(stand$x_age)+5-1,min(max(mod$stand$x_age),65)-1,5))),NULL)
  out$prop_newinfplus5 <- prop_newinfplus5

  return(out)
}

add.refgroup <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  # out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  # out$incid <- incid

  ## cumincid
  incidlist <- lapply(seq(dim(incid)[3]), function(x) incid[ , , x])
  # out$cumincid <- mapply(function(x,y) calc.cumincid(x,mod$stand,aidx=aidx,incid=y),
  #                        param, incidlist)
  # rownames(out$cumincid) <- mod$stand$x_time

  ## WHO prev
  # aidxwho <- which(mod$stand$x_age==15):which(mod$stand$x_age==50)
  # age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidxwho])$y  # prevalence age 15 to 50
  # age.dist <- age.dist/sum(age.dist)
  # out$who15to49prev <- apply(sweep(out$prev[, aidxwho,], 2, age.dist, "*"), c(1,3), sum)
  # rownames(out$who15to49prev) <- mod$stand$x_time

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  # out$pop <- pop
  # plus pop full ages
  popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  # out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
                                                      years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  # out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  # out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  # out$susc_prop <- susc_prop

  ## Prevalence weighted by UNPop national
  prev_pop <- apply(prev, 3, function(x) add_prev_pop(mod, x, pop_prop, aidx, years))
  # out$prev_pop <- prev_pop

  ## Incidence weighted by UNpop national susceptibles
  if(years[2]>max(mod$stand$x_time)){
    topyear <- max(mod$stand$x_time)
  } else { topyear<- years[2]}
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==topyear)
  incid_w <- matrix(nrow=length(times),ncol = dim(prev)[[3]])
  for(i in 1:1000) incid_w[,i] <- rowSums(susc_prop[,,i]*incid[tidx,aidx,i])/rowSums(susc_prop[,,i])
  rownames(incid_w) <- times
  out$incid_w <- incid_w

  ## Mean of period age-specific incidence
  # ageinfunw <- mapply(function(x,y) calc.ageinfunw(x,stand,mod,aidx,years,y), param, incidlist)
  # out$ageinfunw <- ageinfunw

  ## Average age of infection (among susceptible population)
  susclist <- lapply(seq(dim(susc_prop)[3]), function(x) susc_prop[ , , x])# make a list to go in mapply - clean this!
  avgageinf <- mapply(function(x,y,z){calc.avgageinf(x,mod$stand,mod,aidx[-length(aidx)],years,y,z)}, param, susclist, incidlist)
  out$avgageinf <- avgageinf

  ## Number of new infections by age/time
  nsusc2 <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),years))
  dim(nsusc2) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(nsusc2) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  nsusc2list <- lapply(seq(dim(nsusc2)[3]), function(x) nsusc2[ , , x])# make a list to go in mapply - clean this!
  numnewinf <- mapply(function(x,y,z){calc.numnewinf(x,mod$stand,mod,1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),years,y,z)},
                      param, nsusc2list, incidlist)
  dim(numnewinf) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  dimnames(numnewinf) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  # out$numnewinf <- numnewinf

  ## Proportion of new infections by age
  prop_newinf <- apply(numnewinf[,aidx[-length(aidx)],],3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(prop_newinf) <- c(length(times),length(aidx)-1,dim(prev)[3])
  dimnames(prop_newinf) <- list(times,mod$stand$x_age[aidx[-length(aidx)]],NULL)
  # out$prop_newinf <- prop_newinf

  prop_newinf5 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,5/stand$dt), function(i) {
    rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/5,dim(prev)[3])
  dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  out$prop_newinf5 <- prop_newinf5

  # Add 10 yr age groups
  # prop_newinf10 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,10/stand$dt), function(i) {
  #   rowSums(x[,c(i:(i+10/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/10,dim(prev)[3])
  # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)

  # prop_newinfplus <- apply(numnewinf,3, function(x) sweep(x,1,rowSums(x),'/'))
  # dim(prop_newinfplus) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  # dimnames(prop_newinfplus) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  # out$prop_newinfplus <- prop_newinfplus
  #
  # prop_newinfplus5 <- apply(prop_newinfplus,3,function(x) sapply(seq(1,((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt-1)
  #                                                                    ,5/stand$dt), function(i) {
  #                                                                      rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinfplus5) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)*.2/5,dim(prev)[3])
  # dimnames(prop_newinfplus5) <- list(times,c(paste0(seq(min(stand$x_age),min(max(mod$stand$x_age),65)-5,5),"-",
  #                                                   seq(min(stand$x_age)+5-1,min(max(mod$stand$x_age),65)-1,5))),NULL)
  # out$prop_newinfplus5 <- prop_newinfplus5

  return(out)
}

add.refgroup2 <- function(mod, ages=c(15,60), years=c(1990,2017),
                         sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                         pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  # out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  # out$incid <- incid

  ## cumincid
  incidlist <- lapply(seq(dim(incid)[3]), function(x) incid[ , , x])
  # out$cumincid <- mapply(function(x,y) calc.cumincid(x,mod$stand,aidx=aidx,incid=y),
  #                        param, incidlist)
  # rownames(out$cumincid) <- mod$stand$x_time

  ## WHO prev
  # aidxwho <- which(mod$stand$x_age==15):which(mod$stand$x_age==50)
  # age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidxwho])$y  # prevalence age 15 to 50
  # age.dist <- age.dist/sum(age.dist)
  # out$who15to49prev <- apply(sweep(out$prev[, aidxwho,], 2, age.dist, "*"), c(1,3), sum)
  # rownames(out$who15to49prev) <- mod$stand$x_time

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  # out$pop <- pop
  # plus pop full ages
  popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  # out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
                                                      years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  # out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  # out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  # out$susc_prop <- susc_prop

  ## Prevalence weighted by UNPop national
  prev_pop <- apply(prev, 3, function(x) add_prev_pop(mod, x, pop_prop, aidx, years))
  # out$prev_pop <- prev_pop

  ## Incidence weighted by UNpop national susceptibles
  if(years[2]>max(mod$stand$x_time)){
    topyear <- max(mod$stand$x_time)
  } else { topyear<- years[2]}
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==topyear)
  incid_w <- matrix(nrow=length(times),ncol = dim(prev)[[3]])
  for(i in 1:1000) incid_w[,i] <- rowSums(susc_prop[,,i]*incid[tidx,aidx,i])/rowSums(susc_prop[,,i])
  rownames(incid_w) <- times
  # out$incid_w <- incid_w

  ## Mean of period age-specific incidence
  # ageinfunw <- mapply(function(x,y) calc.ageinfunw(x,stand,mod,aidx,years,y), param, incidlist)
  # out$ageinfunw <- ageinfunw

  ## Average age of infection (among susceptible population)
  susclist <- lapply(seq(dim(susc_prop)[3]), function(x) susc_prop[ , , x])# make a list to go in mapply - clean this!
  avgageinf <- mapply(function(x,y,z){calc.avgageinf(x,mod$stand,mod,aidx[-length(aidx)],years,y,z)}, param, susclist, incidlist)
  out$avgageinf <- avgageinf

  ## Number of new infections by age/time
  nsusc2 <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),years))
  dim(nsusc2) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(nsusc2) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  nsusc2list <- lapply(seq(dim(nsusc2)[3]), function(x) nsusc2[ , , x])# make a list to go in mapply - clean this!
  numnewinf <- mapply(function(x,y,z){calc.numnewinf(x,mod$stand,mod,1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),years,y,z)},
                      param, nsusc2list, incidlist)
  dim(numnewinf) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  dimnames(numnewinf) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)

  ncases_ages <- numnewinf[,colnames(numnewinf)>=ages[1] & colnames(numnewinf)<ages[2],]
  ncases_ages <- apply(ncases_ages,c(1,3),sum)

  # out$numnewinf <- numnewinf
  #
  # ## Proportion of new infections by age
  # prop_newinf <- apply(numnewinf[,aidx[-length(aidx)],],3, function(x) sweep(x,1,rowSums(x),'/'))
  # dim(prop_newinf) <- c(length(times),length(aidx)-1,dim(prev)[3])
  # dimnames(prop_newinf) <- list(times,mod$stand$x_age[aidx[-length(aidx)]],NULL)
  # # out$prop_newinf <- prop_newinf
  #
  # prop_newinf5 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,5/stand$dt), function(i) {
  #   rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/5,dim(prev)[3])
  # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  # out$prop_newinf5 <- prop_newinf5
  #
  # # Add 10 yr age groups
  # # prop_newinf10 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,10/stand$dt), function(i) {
  # #   rowSums(x[,c(i:(i+10/stand$dt-1))], na.rm=T)}))
  # # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/10,dim(prev)[3])
  # # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  #
  # # prop_newinfplus <- apply(numnewinf,3, function(x) sweep(x,1,rowSums(x),'/'))
  # # dim(prop_newinfplus) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  # # dimnames(prop_newinfplus) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  # # out$prop_newinfplus <- prop_newinfplus
  # #
  # # prop_newinfplus5 <- apply(prop_newinfplus,3,function(x) sapply(seq(1,((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt-1)
  # #                                                                    ,5/stand$dt), function(i) {
  # #                                                                      rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  # # dim(prop_newinfplus5) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)*.2/5,dim(prev)[3])
  # # dimnames(prop_newinfplus5) <- list(times,c(paste0(seq(min(stand$x_age),min(max(mod$stand$x_age),65)-5,5),"-",
  # #                                                   seq(min(stand$x_age)+5-1,min(max(mod$stand$x_age),65)-1,5))),NULL)
  # # out$prop_newinfplus5 <- prop_newinfplus5

  return(ncases_ages)
}

add.incages <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$incid <- incid

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  out$pop <- pop
  # plus pop full ages
  popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
                                                      years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  out$susc_prop <- susc_prop

  ## Incidence weighted by UNpop national susceptibles
  if(years[2]>max(mod$stand$x_time)){
    topyear <- max(mod$stand$x_time)
  } else { topyear<- years[2]}
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==topyear)
  incid_w <- matrix(nrow=length(times),ncol = dim(incid)[[3]])
  for(i in 1:1000) incid_w[,i] <- rowSums(susc_prop[,,i]*incid[tidx,aidx,i])/rowSums(susc_prop[,,i])
  rownames(incid_w) <- times
  out$incid_w <- incid_w

  return(incid_w)
}

add.prevages <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  # system.time(incid <- lapply(param, calc.incid, mod$stand))
  # incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  # dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  # out$incid <- incid

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  out$pop <- pop
  # plus pop full ages
  popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
                                                      years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  out$susc_prop <- susc_prop

  age10 <- seq(ages[1],ages[2],10)

  prev_pop <- list()
  for(i in 1:(length(age10)-1)) {
    aidx_new <- which(stand$x_age==age10[i]):which(stand$x_age==age10[i+1])
    ## Prevalence weighted by UNPop national
    prev_pop[[i]] <- apply(prev, 3, function(x)
      add_prev_pop(mod, x, pop_prop[,colnames(pop_prop)>=age10[i] & colnames(pop_prop)<=age10[i+1]],
                                                           aidx_new, years))
  }
  names(prev_pop) <- paste0("prev",age10[-length(age10)],(age10[-1]-1))

  return(prev_pop)
}

add.inc_isoadd <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)

  return(incid)
}

add.incid_w <- function(mod, ages=c(15,60), years=c(1990,2017),
                        sex = "Female", country = "Malawi", nonUNPOP = FALSE,
                        pop_data = pop_data, ind, tot=2){
  param <- create.param.list(mod$fit, ind, tot)
  stand <- mod$stand
  stand$artstart_tIDX <- stand$artstart_tIDX[ind]
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]
  if(stand$x_age[1]!=10){
    stop("ERROR: fix prev function for non-age 10 start!")
  }
  aidx <- which(stand$x_age==ages[1]):which(stand$x_age==ages[2])
  out <- list()

  ## prev
  system.time(prev <- mclapply(param, calc.prev, mod$stand))
  prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$prev <- prev

  # not sure this will work - likely need to update so it's "x_age_inc" to not mess up
  # where x_age gets used for natural mortality
  # stand$x_age <- stand$x_age[stand$test_aIDX:stand$STEPS_age]
  # mod$stand$x_age <- mod$stand$x_age[mod$stand$test_aIDX:mod$stand$STEPS_age]

  ## incid
  system.time(incid <- lapply(param, calc.incid, mod$stand))
  incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
  dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
  out$incid <- incid

  # ## cumincid
  # out$cumincid <- sapply(param, calc.cumincid, mod$stand, aidx=aidx) # incidence between age 15 and 60
  # rownames(out$cumincid) <- mod$stand$x_time
  #
  # ## WHO prev
  # aidxwho <- which(mod$stand$x_age==15):which(mod$stand$x_age==50)
  # age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidxwho])$y  # prevalence age 15 to 50
  # age.dist <- age.dist/sum(age.dist)
  # out$who15to49prev <- apply(sweep(out$prev[, aidxwho,], 2, age.dist, "*"), c(1,3), sum)
  # rownames(out$who15to49prev) <- mod$stand$x_time

  # Nest the next few in a if (nonUNPOP==FALSE){}

  ## UNPop
  if(!"10-14" %in% filter(pop_data,Location==country & Sex==sex)$Age) {
    extrarow <- filter(pop_data,Location==country & Sex==sex & Age=="15-19")
    extrarow$Age <- "10-14"
    extrarow$age_low <- 10
    extrarow$age_high <- 15
    pop_data <- rbind(pop_data,extrarow)
  }
  if(!"55-59" %in% filter(pop_data,Location==country & Sex==sex)$Age) {
    extrarow <- filter(pop_data,Location==country & Sex==sex & Age=="50-54")
    extrarow$Age <- "55-59"
    extrarow$age_low <- 55
    extrarow$age_high <- 60
    pop_data <- rbind(pop_data,extrarow)
  }

  pop <- add_pop(mod, pop_data, ages, country, sex, years)
  out$pop <- pop
  # plus pop full ages
  # popplus <- add_pop(mod, pop_data, ages = c(min(mod$stand$x_age),min(max(mod$stand$x_age),65)), country, sex, years)
  # out$popplus <- popplus
  pop_prop <- sweep(pop,1,rowSums(pop),`/`)
  out$pop_prop <- pop_prop

  ## Number of HIV cases in UNpop
  # ncases <- apply(prev, 3, function(x) add_ncases_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),
  #                                                     years))
  times <- seq(years[1],min(max(mod$stand$x_time),years[2]),by=mod$stand$dt)
  # dim(ncases) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  # dimnames(ncases) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  # out$ncases <- ncases

  ## Number susceptible
  nsusc <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, pop, aidx,years))
  dim(nsusc) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(nsusc) <- list(times,mod$stand$x_age[aidx],NULL)
  out$nsusc <- nsusc
  ## Susceptible proportion
  susc_prop <- apply(nsusc,3, function(x) sweep(x,1,rowSums(x),'/'))
  dim(susc_prop) <- c(length(times),length(aidx),dim(prev)[3])
  dimnames(susc_prop) <- list(times,mod$stand$x_age[aidx],NULL)
  out$susc_prop <- susc_prop

  ## Prevalence weighted by UNPop national
  prev_pop <- apply(prev, 3, function(x) add_prev_pop(mod, x, pop_prop, aidx, years))
  out$prev_pop <- prev_pop

  ## Incidence weighted by UNpop national susceptibles
  if(years[2]>max(mod$stand$x_time)){
    topyear <- max(mod$stand$x_time)
  } else { topyear<- years[2]}
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==topyear)
  incid_w <- matrix(nrow=length(times),ncol = dim(incid)[[3]])
  for(i in 1:1000) incid_w[,i] <- rowSums(susc_prop[,,i]*incid[tidx,aidx,i])/rowSums(susc_prop[,,i])
  rownames(incid_w) <- times
  out$incid_w <- incid_w

  # ## Mean of period age-specific incidence
  # ageinfunw <- sapply(param, calc.ageinfunw, stand, mod, aidx, years)
  # out$ageinfunw <- ageinfunw
  #
  # ## Average age of infection (among susceptible population)
  # susclist <- lapply(seq(dim(susc_prop)[3]), function(x) susc_prop[ , , x])# make a list to go in mapply - clean this!
  # avgageinf <- mapply(function(x,y){calc.avgageinf(x,mod$stand,mod,aidx[-length(aidx)],years,y)}, param, susclist)
  # out$avgageinf <- avgageinf
  #
  # ## Number of new infections by age/time
  # nsusc2 <- apply(prev, 3, function(x) add_nsusc_pop(mod, x, popplus, 1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),years))
  # dim(nsusc2) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1),dim(prev)[3])
  # dimnames(nsusc2) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt+1)],NULL)
  # nsusc2list <- lapply(seq(dim(nsusc2)[3]), function(x) nsusc2[ , , x])# make a list to go in mapply - clean this!
  # numnewinf <- mapply(function(x,y){calc.numnewinf(x,mod$stand,mod,1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),years,y)}, param, nsusc2list)
  # dim(numnewinf) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  # dimnames(numnewinf) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  # out$numnewinf <- numnewinf
  #
  # ## Proportion of new infections by age
  # prop_newinf <- apply(numnewinf[,aidx[-length(aidx)],],3, function(x) sweep(x,1,rowSums(x),'/'))
  # dim(prop_newinf) <- c(length(times),length(aidx)-1,dim(prev)[3])
  # dimnames(prop_newinf) <- list(times,mod$stand$x_age[aidx[-length(aidx)]],NULL)
  # out$prop_newinf <- prop_newinf
  #
  # prop_newinf5 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,5/stand$dt), function(i) {
  #   rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/5,dim(prev)[3])
  # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  # out$prop_newinf5 <- prop_newinf5
  #
  # # Add 10 yr age groups
  # # prop_newinf10 <- apply(prop_newinf,3,function(x) sapply(seq(1,length(aidx)-1,10/stand$dt), function(i) {
  # #   rowSums(x[,c(i:(i+10/stand$dt-1))], na.rm=T)}))
  # # dim(prop_newinf5) <- c(length(times),(length(aidx)-1)*stand$dt/10,dim(prev)[3])
  # # dimnames(prop_newinf5) <- list(times, c(paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+5-1,ages[2]-1,5))), NULL)
  #
  # prop_newinfplus <- apply(numnewinf,3, function(x) sweep(x,1,rowSums(x),'/'))
  # dim(prop_newinfplus) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt),dim(prev)[3])
  # dimnames(prop_newinfplus) <- list(times,mod$stand$x_age[1:((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)],NULL)
  # out$prop_newinfplus <- prop_newinfplus
  #
  # prop_newinfplus5 <- apply(prop_newinfplus,3,function(x) sapply(seq(1,((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt-1)
  #                                                                    ,5/stand$dt), function(i) {
  #                                                                      rowSums(x[,c(i:(i+5/stand$dt-1))], na.rm=T)}))
  # dim(prop_newinfplus5) <- c(length(times),((min(max(mod$stand$x_age),65)-min(mod$stand$x_age))/stand$dt)*.2/5,dim(prev)[3])
  # dimnames(prop_newinfplus5) <- list(times,c(paste0(seq(min(stand$x_age),min(max(mod$stand$x_age),65)-5,5),"-",
  #                                                   seq(min(stand$x_age)+5-1,min(max(mod$stand$x_age),65)-1,5))),NULL)
  # out$prop_newinfplus5 <- prop_newinfplus5

  out <- out$incid_w

  return(out)
}

# library('dplyr')

add.propinfage <- function(mod,prev,pop_data=pop_data,ages=c(15,55),
                           country = "Malawi",sex = "Female", years=c(1990,2015),
                           byyears=10){
  param <- create.param.list(mod$fit)
  # prev <- mclapply(param, calc.prev, mod$stand)

  standard_pop <- as.data.frame(subset(pop_data,Sex==sex &
                                         Location==country))
  rownames(standard_pop) <- standard_pop$age_low + (standard_pop$age_high - standard_pop$age_low)/2
  standard_pop <- dplyr::select(standard_pop, - c("Location","Sex","Age","age_low","age_high"))
  standard_pop <- standard_pop[, names(standard_pop)>=years[1] & names(standard_pop)<=years[2]]
  startind <- which(mod$stand$x_age==ages[1])
  endind <- which(mod$stand$x_age==ages[2])
  aidx <- startind:endind
  age_dist <- apply(standard_pop,2,function(x) {approx(as.numeric(rownames(standard_pop)), x, mod$stand$x_age[aidx])$y})
  age_dist <- apply(age_dist,1,function(x) {approx(as.numeric(names(standard_pop)), x, seq(years[1],years[2],by=.2))$y})
  age_dist <- sweep(age_dist,1,rowSums(age_dist),`/`)

  weights <- lapply(prev, function(x) {gen_uninfec_weights(age_dist,x,mod,ages,years)})
  props_by_agegr <- Map(function(x,y){calc.props_by_agrgr(x,mod$stand,mod,ages,years,y,byyears=byyears)}, param, weights)
  return(props_by_agegr)
}

add_pop <- function(mod,pop_data = pop_data,ages=c(15,60),country="Malawi",
                    sex="Female",years=c(1990,2017)) {
  standard_pop <- as.data.frame(subset(pop_data,Sex==sex &
                                         Location==country))
  rownames(standard_pop) <- standard_pop$age_low + (standard_pop$age_high - standard_pop$age_low)/2
  standard_pop <- dplyr::select(standard_pop, - c("Location","Sex","Age","age_low","age_high"))
  standard_pop <- standard_pop[, names(standard_pop)>=years[1] & names(standard_pop)<=years[2]]
  aidx <- which(mod$stand$x_age==ages[1]):which(mod$stand$x_age==ages[2])
  # aidx <-   1:((ages[2]-ages[1])/mod$stand$dt+1L)
  age_dist <- apply(standard_pop,2,function(x) {approx(as.numeric(rownames(standard_pop)), x, mod$stand$x_age[aidx])$y})
  age_dist <- age_dist / ((1/mod$stand$dt)*5)
  age_dist <- apply(age_dist,1,function(x) {approx(as.numeric(names(standard_pop)), x, seq(years[1],years[2],by=.2))$y})
  dimnames(age_dist) <- list(seq(years[1],years[2],mod$stand$dt), seq(ages[1],ages[2],mod$stand$dt))
  return(age_dist)
}

add_ncases_pop <- function(mod, prev, pop, aidx,years=c(1990,2017)) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  pop <- pop[1:length(tidx),]
  pop_prop <- pop[1:length(tidx),]
  ncases <- prev[tidx,aidx]*pop
  return(ncases)
}

add_nsusc_pop <- function(mod, prev, pop, aidx,years=c(1990,2017)) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  pop <- pop[1:length(tidx),]
  pop_prop <- pop[1:length(tidx),]
  nsusc <- pop - prev[tidx,aidx]*pop
  return(nsusc)
}

add_prev_pop <- function(mod, prev, pop_prop, aidx,years=c(1990,2017)) {
  if(years[2]>max(mod$stand$x_time)){
    years[2] <- max(mod$stand$x_time)
  }
  tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
  pop_prop <- pop_prop[1:length(tidx),]
  prev_pop <- rowSums(prev[tidx,aidx]*pop_prop)/rowSums(pop_prop)
  return(prev_pop)
}

# add_prev_pop <- function(mod,pop_data,prev,ages=c(15,60),country="Malawi",
#                          sex="Female",years=c(1990,2017)){
#   param <- create.param.list(mod$fit)
#   # prev <- mclapply(param, calc.prev, mod$stand)
#
#   standard_pop <- as.data.frame(subset(pop_data,Sex==sex &
#                                          Location==country))
#   rownames(standard_pop) <- standard_pop$age_low + (standard_pop$age_high - standard_pop$age_low)/2
#   standard_pop <- select(standard_pop, - c("Location","Sex","Age","age_low","age_high"))
#   standard_pop <- standard_pop[, names(standard_pop)>=years[1] & names(standard_pop)<=years[2]]
#   aidx <-   1:((ages[2]-ages[1])/mod$stand$dt+1L)
#   age_dist <- apply(standard_pop,2,function(x) {approx(as.numeric(rownames(standard_pop)), x, mod$stand$x_age[aidx])$y})
#   age_dist <- apply(age_dist,1,function(x) {approx(as.numeric(names(standard_pop)), x, seq(years[1],years[2],by=.2))$y})
#   age_dist <- sweep(age_dist,1,rowSums(age_dist),`/`)
#   aidx <-   1:((ages[2]-ages[1])/mod$stand$dt+1L)
#   if(years[2]>max(mod$stand$x_time)){
#     years[2] <- max(mod$stand$x_time)
#   }
#   tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
#   age_dist <- age_dist[1:length(tidx),]
#   prev_pop <- sapply(prev,function(x) rowSums(x[tidx,aidx]*age_dist)/rowSums(age_dist))
#   rownames(prev_pop) <- seq(years[1],years[2],mod$stand$dt)
#   return(prev_pop)
# }

# add.incpop <- function(mod,pop_data=pop_data,prev,ages=c(15,55),
#                        country="Malawi",sex="Female",years=c(1990,2015)){
#   param <- create.param.list(mod$fit)
#   # prev <- mclapply(param, calc.prev, mod$stand)
#
#   standard_pop <- as.data.frame(subset(pop_data,Sex==sex &
#                                          Location==country))
#   rownames(standard_pop) <- standard_pop$age_low + (standard_pop$age_high - standard_pop$age_low)/2
#   standard_pop <- select(standard_pop, - c("Location","Sex","Age","age_low","age_high"))
#   standard_pop <- standard_pop[, names(standard_pop)>=years[1] & names(standard_pop)<=years[2]]
#   aidx <-   1:((ages[2]-ages[1])/mod$stand$dt+1L)
#   age_dist <- apply(standard_pop,2,function(x) {approx(as.numeric(rownames(standard_pop)), x, mod$stand$x_age[aidx])$y})
#   age_dist <- apply(age_dist,1,function(x) {approx(as.numeric(names(standard_pop)), x, seq(years[1],years[2],by=.2))$y})
#   age_dist <- sweep(age_dist,1,rowSums(age_dist),`/`)
#   weights <- lapply(prev, function(x) {gen_uninfec_weights(age_dist,x,mod,ages,years)})
#   if(years[2]>max(mod$stand$x_time)){
#     years[2] <- max(mod$stand$x_time)
#   }
#   tidx <- which(mod$stand$x_time==years[1]):which(mod$stand$x_time==years[2])
#   incid <- lapply(param, calc.incid, mod$stand)
#   incid_w <- mapply(function(x,y) rowSums(x[tidx,aidx]*y)/rowSums(y),incid,weights)
#   rownames(incid_w) <- seq(years[1],years[2],mod$stand$dt)
#   return(incid_w)
# }

# add.ageinfecunw <- function(mod,ages=c(15,60), years=c(1990,2015)){
#   param <- create.param.list(mod$fit)
#   #avgageinf <- mapply(function(x,y){calc.avgageinf(x,mod$stand,mod,ages,years,y)}, param, weights)
#   avgageunw <- sapply(param, calc.ageinfunw, mod$stand, mod, ages = ages, years=years) # avg age of infection between 15 and 60
#   return(avgageunw)
# }

# # function for average age of infection
# add.ageinfec <- function(mod,pop_data=pop_data,ages=c(15,60),
#                          country = "Malawi",sex = "Female", years=c(1990,2015)){
#   param <- create.param.list(mod$fit)
#   # prev <- mclapply(param, calc.prev, mod$stand)
#
#   standard_pop <- as.data.frame(subset(pop_data,Sex==sex &
#                                          Location==country))
#   rownames(standard_pop) <- standard_pop$age_low + (standard_pop$age_high - standard_pop$age_low)/2
#   standard_pop <- select(standard_pop, - c("Location","Sex","Age","age_low","age_high"))
#   standard_pop <- standard_pop[, names(standard_pop)>=years[1] & names(standard_pop)<=years[2]]
#   aidx <-   1:((ages[2]-ages[1])/mod$stand$dt+1L)
#   age_dist <- apply(standard_pop,2,function(x) {approx(as.numeric(rownames(standard_pop)), x, mod$stand$x_age[aidx])$y})
#   age_dist <- apply(age_dist,1,function(x) {approx(as.numeric(names(standard_pop)), x, seq(years[1],years[2],by=.2))$y})
#   age_dist <- sweep(age_dist,1,rowSums(age_dist),`/`)
#
#   weights <- lapply(prev, function(x) {gen_uninfec_weights(age_dist,x,mod,ages,years)})
#   avgageinf <- mapply(function(x,y){calc.avgageinf(x,mod$stand,mod,ages,years,y)}, param, weights)
#   # avgageinf <- sapply(param, calc.avgageinf, mod$stand, years=45) # avg age of infection between 15 and 60
#   return(avgageinf)
# }

add.prev <-  function(mod){
  param <- create.param.list(mod$fit)
  prev <- mclapply(param, calc.prev, mod$stand)
  return(prev)
}

add_cumincid_coh <- function(mod, years=c(15,45), ind, tot=2) {
  param <- create.param.list(mod$fit,ind,tot)

  cumincid_coh <- lapply(param,calc.cumincid.coh.mult, mod$stand, years=years) # cumulative incidence
  return(cumincid_coh)
}

add_incid_proj <- function(mod, years=c(15,55), ind, tot=2, karman=FALSE) {
  param <- create.param.list(mod$fit,ind,tot)

  incid <- lapply(param,calc.incid.proj, mod$stand, years=years,karman) # incidence
  if(karman==TRUE) {
    tidx <- which(mod$stand$x_time<=2012)
  } else {
    tidx <- 1:length(mod$stand$x_time)
  }
  newtime <- c(mod$stand$x_time[tidx],seq(max(mod$stand$x_time[tidx])+.2,1995+years[2],.2))
  newage <- seq(years[1],years[2],0.2)
  incid <- array(unlist(incid), c(length(newtime), length(newage), length(param)))
  dimnames(incid) <- list(newtime, newage, NULL)
  return(incid)
}

add_incid_proj2 <- function(mod, years=c(15,55), ind, tot=2, karman=FALSE, reldec,
                            meetyear=2022) {
  param <- create.param.list(mod$fit,ind,tot)

  incid <- lapply(param,calc.incid.proj2, mod$stand, years=years,karman,reldec,meetyear) # incidence
  if(karman==TRUE) {
    tidx <- which(mod$stand$x_time<=2012)
  } else {
    tidx <- 1:length(mod$stand$x_time)
  }
  newtime <- c(mod$stand$x_time[tidx],seq(max(mod$stand$x_time[tidx])+.2,1995+years[2],.2))
  newage <- seq(years[1],years[2],0.2)
  incid <- array(unlist(incid), c(length(newtime), length(newage), length(param)))
  dimnames(incid) <- list(newtime, newage, NULL)
  return(incid)
}

add_period_cumincid_proj <- function(incid,mod,ages=c(15,55)) {
  cumincid <- apply(incid,3,calc.period.cumincid.proj,stand = mod$stand,ages = ages)
  return(cumincid)
}

add_cumincid_coh_proj <- function(mod, years=c(15,55), ind, tot=2, karman=FALSE) {
  param <- create.param.list(mod$fit,ind,tot)

  cumincid_coh <- lapply(param,calc.cumincid.coh.mult.proj, mod$stand, years=years,karman) # cumulative incidence
  return(cumincid_coh)
}

add_cumincid_coh_proj2 <- function(mod, years=c(15,55), ind, tot=2, karman=FALSE, reldec,
                                   meetyear=2022) {
  param <- create.param.list(mod$fit,ind,tot)

  cumincid_coh <- lapply(param,calc.cumincid.coh.mult.proj2, mod$stand, years=years,karman,reldec,
                         meetyear) # cumulative incidence
  return(cumincid_coh)
}

find_max_int <- function(newinf,window_len=10,ages = c(15,55)) {
  aidx <- which(as.numeric(colnames(newinf))==ages[1]):which(as.numeric(colnames(newinf))==ages[2])
  newinf <- newinf[,aidx,]
  newinf <- apply(newinf,c(1,2),median)
  winds <- sapply(1:(ncol(newinf)-window_len/.2+1),
                  function(x) rowSums(newinf[,x:(x+(window_len/.2)-1)]))
  ageinds <- apply(winds,1,function(x) which(x==max(x)))
  prop <- apply(winds,1,function(x) x[which(x==max(x))])/rowSums(newinf)
  lowage <- as.numeric(colnames(newinf)[ageinds])
  highage <- lowage + window_len
  years <- as.numeric(rownames(newinf))
  return(data.frame(lowage = lowage, highage=highage, years = years, prop=prop))
}

pops_study <- function(res,studyname,ages,years) {
  study_res <- subset(res,site==studyname)
  study_res <- filter(study_res,!is.na(sex))
  if(min(study_res$entry)>years[1]) {
    years[1] <- ceiling(min(study_res$entry))
  }
  if(max(study_res$exit)<years[2]) {
    years[2] <- floor(max(study_res$exit))
  }
  n <- length(years[1]:years[2])
  study_res <- lapply(years[1]:years[2], function(x) {
    dat <- filter(study_res,entry<=x & exit>=x);
    dat <- mutate(dat, age = cut(x - dob, seq(ages[1],ages[2],by=5),right=FALSE,
                                 labels=paste0(seq(ages[1],ages[2]-5,by=5),"-",seq(ages[1]+4,ages[2]-1,by=5))));
    summarise(group_by(dat, sex, age),
              count=n())})
  study_pop_tot <- merge(study_res[[1]],study_res[[2]],by=c("sex","age"),all=TRUE)
  for(i in 3:n) {
    study_pop_tot <- merge(study_pop_tot,study_res[[i]],by=c("sex","age"),all=TRUE)
  }
  names(study_pop_tot)[3:(n+2)] <- seq(floor(years[1]),floor(years[2]))
  study_pop <- filter(study_pop_tot,age %in%
                        paste0(seq(ages[1],ages[2]-5,5),"-",seq(ages[1]+4,ages[2],5)))
  study_pop$Location <- studyname
  names(study_pop)[1:2] <- c("Sex","Age")
  return(study_pop)
}

pops_study_inc <- function(res,studyname,ages,years) {
  test <- do.call("rbind",res)
  test2 <- summarise(group_by(test,year,age,sex), pyrs=mean(pyears))
  test2 <- mutate(test2, agecat = cut(age, seq(15,55,5),right=FALSE,
                                      labels = paste0(seq(15,50,5),"-",seq(19,54,5))))
  test2 <- summarise(group_by(test2,year,agecat,sex),pyrs=sum(pyrs))
  out <- tidyr::pivot_wider(test2, id_cols=c(agecat,sex),
                     names_from=year, values_from=pyrs)
  yearsadd <- as.character(c(1985:2019))[!1985:2019 %in% names(out)]
  for(i in yearsadd) {
    out <- mutate(out, !!i := NA_real_)
  }
  out$Location <- studyname
  out$agecat <- as.character(out$agecat)
  out$sex <- as.character(out$sex)
  out <- mutate(out,
                           age_low = as.numeric(substr(agecat,1,regexpr("-",agecat)-1)),
                           age_high = as.numeric(substr(agecat,
                                                        regexpr("-",agecat)+1,nchar(agecat))) + 1)
  out <- rename(out, Age = agecat, Sex = sex)
  out <- select(out,Location,Sex,Age,as.character(1985:2019),age_low,age_high)
  out <- ungroup(out)
  return(out)
}

clean_for_pop_pyramid <- function(data) {
  site <- data$Location
  data <- dplyr::select(data,-"Location")

  data$Age <- as.factor(data$Age)
  data <- reshape(data, varying = list(3:(length(data))), v.names = "pop",
                  idvar = c("Sex","Age"), direction="long",
                  times = as.numeric(names(data)[3:length(data)]))
  data$site <- site[1]

  data <- mutate(group_by(data,Sex,time),
                 prop = pop/sum(pop))
  return(data)
}


##########################
####  Plot functions  ####
##########################

pop_pyr_plot <- function(dat) {
  # library('gganimate')
  ggplot(data=dat, aes(x = Age, y = pop, fill=Sex)) +
    geom_col(data = subset(dat, Sex=="Female"),position=position_dodge(0.9)) +
    geom_col(data = subset(dat, Sex=="Male"), aes(y=-pop),position=position_dodge(0.9)) +
    geom_hline(aes(yintercept=0)) +
    annotate("text",x="50-54",y=-1*max(dat$pop[dat$Sex=="Male"])+1000,label="Men",size=8,vjust=0) +
    annotate("text",x="50-54",y=max(dat$pop[dat$Sex=="Female"])-1000,label="Women",size=8,vjust=0) +
    coord_flip() +
    # scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    # facet_wrap(~Location,scales="free_x")+
    labs(subtitle = 'Year: {floor(frame_time)}', x = 'Age', y = 'Population') +
    transition_time(time) +
    ease_aes('linear') +
    # scale_y_continuous(breaks=seq(-.3,.3,.1),
    #                    labels=c(as.character(c(seq(.3,0,-.1), seq(.1,.3,.1))))) +
    scale_y_continuous(breaks = seq(-floor(max(dat$pop[dat$Sex=="Male"])/1000)*1000,
                                    floor(max(dat$pop[dat$Sex=="Female"])/1000)*1000,
                                    1000),
                       labels = as.character(c(seq(floor(max(dat$pop[dat$Sex=="Male"])/1000)*1000,0,-1000),
                                               seq(1000,floor(max(dat$pop[dat$Sex=="Female"])/1000)*1000,1000)))) +
    scale_fill_discrete("") +
    theme(legend.position = "none")+
    ggtitle(dat$site[1])
}

pop_pyr_plot_incprev <- function(datf,datm,main,avgagef,avgagem) {
  datage <- rbind(mutate(avgagef,sex="Female"), mutate(avgagem, sex="Male"))
  dat <- rbind(mutate(datf,sex="Female"), mutate(datm,sex="Male"))
  dat$type <- as.factor(dat$type)
  dat$year <- as.numeric(dat$year)
  dat <- merge(dat,datage)
  dat <- filter(dat,floor(year)==year)
  dat$type <- factor(dat$type,levels=c("numnewinf","ncases","nsusc"))
  dat$avgagetrans <- (dat$avgage-17.5)/5+1
  ggplot(data=dat, aes(x = age, y = values, fill=type)) +
    geom_col(data = subset(dat, sex=="Female")) +
    geom_col(data = subset(dat, sex=="Male"), aes(y=-values)) +
    geom_hline(aes(yintercept=0)) +
    geom_segment(data = subset(dat, sex=="Female"),
                 aes(x=avgagetrans,y=0,yend=max(dat$values[dat$sex=="Female"]),xend=avgagetrans,lty="Avg age of infection"))+
    geom_segment(data = subset(dat, sex=="Male"),
                 aes(x=avgagetrans,y=0,yend=-max(dat$values[dat$sex=="Male"]),xend=avgagetrans,lty="Avg age of infection"))+
    annotate("text",x="50-54",y=-.75*max(dat$values[dat$sex=="Male"]),label="Men",size=3,vjust=0) +
    annotate("text",x="50-54",y=.75*max(dat$values[dat$sex=="Female"]),label="Women",size=3,vjust=0) +
    geom_text(data = subset(dat, sex=="Female"),size=2.5,
              aes(x=avgagetrans+.5,y=.9*max(dat$values[dat$sex=="Female"]),label=round(avgage,1)))+
    geom_text(data = subset(dat, sex=="Male"),size=2.5,
              aes(x=avgagetrans+.5,y=-.9*max(dat$values[dat$sex=="Male"]),label=round(avgage,1)))+
    coord_flip() +
    # scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    facet_wrap(~year,ncol=4)+
    scale_linetype_manual("",values="dashed") +
    # scale_y_continuous(breaks=seq(-.3,.3,.1),
    #                    labels=c(as.character(c(seq(.3,0,-.1), seq(.1,.3,.1))))) +
    scale_y_continuous("Population",labels=abs) +
    scale_fill_discrete("",labels=c("New infections","Prevalent cases","HIV-negative"),
                        l=55) +
    theme(legend.position = "bottom")+
    ggtitle(main)+
    xlab("Age")
}

pop_pyr_plot_nonanim <- function(dat,yearrange=c(2000,2017),main) {
  dat <- do.call(rbind,dat)
  dat <- subset(dat,time>=yearrange[1] & time<=yearrange[2])
  if(main=="uMkhanyakude/South Africa") {
    dat$site <- factor(dat$site,levels=c("uMkhanyakude","South Africa"))
  }
  if(length(unique(dat$site))==1) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 55, c = 100)[1:n]
    }
    colplot <- gg_color_hue(2)[2]
    out <- ggplot(data=dat, aes(x = Age, y = prop,fill=site)) +
      geom_col(data = subset(dat, Sex=="Female"),position=position_dodge(0.9),col=colplot) +
      geom_col(data = subset(dat, Sex=="Male"), aes(y=-prop),position=position_dodge(0.9),col=colplot) +
      geom_hline(aes(yintercept=0)) +
      annotate("text",x="45-49",y=-.17,label="Men",size=3,vjust=0) +
      annotate("text",x="45-49",y=.17,label="Women",size=3,vjust=0) +
      coord_flip() +
      # scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      facet_wrap(~time,ncol=4)+
      # scale_y_continuous(breaks=seq(-.3,.3,.1),
      #                    labels=c(as.character(c(seq(.3,0,-.1), seq(.1,.3,.1))))) +
      scale_y_continuous(breaks = seq(-.3,.3,.1),
                         labels = as.character(c(seq(.3,0,-.1), seq(.1,.3,.1)))) +
      scale_fill_manual("",values = colplot) +
      theme(legend.position = "bottom")+
      ggtitle(main) +
      ylab("Proportion")
  } else {
    out <- ggplot(data=dat, aes(x = Age, y = prop, fill=site)) +
      geom_col(data = subset(dat, Sex=="Female"),position=position_dodge(0.9)) +
      geom_col(data = subset(dat, Sex=="Male"), aes(y=-prop),position=position_dodge(0.9)) +
      geom_hline(aes(yintercept=0)) +
      annotate("text",x="45-49",y=-.17,label="Men",size=3,vjust=0) +
      annotate("text",x="45-49",y=.17,label="Women",size=3,vjust=0) +
      coord_flip() +
      # scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      facet_wrap(~time,ncol=4)+
      # scale_y_continuous(breaks=seq(-.3,.3,.1),
      #                    labels=c(as.character(c(seq(.3,0,-.1), seq(.1,.3,.1))))) +
      scale_y_continuous(breaks = seq(-.3,.3,.1),
                         labels = as.character(c(seq(.3,0,-.1), seq(.1,.3,.1)))) +
      scale_fill_discrete("",l=55) +
      theme(legend.position = "bottom")+
      ggtitle(main) +
      ylab("Proportion")
  }
  return(out)
}

pop_pyr_prop_plot <- function(dat) {
  # library('gganimate')
  ggplot(data=dat, aes(x = Age, y = prop, fill=Sex)) +
    geom_col(data = subset(dat, Sex=="Female"),position=position_dodge(0.9)) +
    geom_col(data = subset(dat, Sex=="Male"), aes(y=-prop),position=position_dodge(0.9)) +
    geom_hline(aes(yintercept=0)) +
    annotate("text",x="50-54",y=-.2,label="Men",size=8,vjust=0) +
    annotate("text",x="50-54",y=.2,label="Women",size=8,vjust=0) +
    coord_flip() +
    # scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    # facet_wrap(~Location,scales="free_x")+
    labs(subtitle = 'Year: {floor(frame_time)}', x = 'Age', y = 'Population') +
    transition_time(time) +
    ease_aes('linear') +
    # scale_y_continuous(breaks=seq(-.3,.3,.1),
    #                    labels=c(as.character(c(seq(.3,0,-.1), seq(.1,.3,.1))))) +
    scale_y_continuous(breaks = seq(-.3,.3,.1),
                       labels = as.character(c(seq(.3,0,-.1), seq(.1,.3,.1)))) +
    scale_fill_discrete("") +
    theme(legend.position = "none")+
    ggtitle(dat$site[1])
}

plot_inc_by_age <- function(dat,xlim=c(2000,2017),agelim,ymax=NA,main,
                            legendonly=FALSE) {
  out <- ggplot(data=filter(dat,year<=xlim[2] & year>=xlim[1] & age %in% agelim),
                aes(x=year,col=age)) +
    geom_line(aes(y=median*100,lty="Median"),size=1) +
    geom_line(aes(y=mean*100,lty="Mean"),size=1) +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    scale_y_continuous(expand=c(0,0),limits=c(0,ymax))+
    theme_bw()+
    labs(x="Year",
         y="Incidence (per 100 py)",
         color="Age",
         title=main,
         linetype="")+
    theme(legend.position="bottom")
  if(legendonly==FALSE) {
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

plot.cuminc <- function(cuminc, stand, main, ylim=c(0, 1), xlim=c(1985, 2015)){
  par(las=1, tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(2.1, 3.1, 2.1, 0.5))
  ##
  xx <- stand$x_time[-1]
  idx.est <- 2:stand$STEPS_time - 1
  ##
  ###EDITING HERE - q4515.obs instead of q4515
  plot(xx[idx.est], rowMeans(cuminc$cumincid[idx.est,"i45_15",]), type='n', lty=1, ylim=ylim, xlim=xlim,
       ylab="", xlab="", main=main)
  ##
  cred.region(xx[idx.est], apply(cuminc$cumincid[idx.est,"i45_15",], 1, quantile, c(0.025, 0.975)), col=transp(2))
  lines(xx[idx.est],  rowMeans(cuminc$cumincid[idx.est,"i45_15",]), col=2, lwd=1.5)
  ##
  return(invisible(NULL))
}

plot.cuminc_update <- function(cuminc, stand, main, ylim=c(0, 1), xlim=c(1985, 2015)){
  ##
  xx <- stand$x_time[-1]
  idx.est <- 2:stand$STEPS_time - 1
  data <- data.frame(x = xx[idx.est],
                     y = rowMeans(cuminc$cumincid[idx.est,"i45_15",]),
                     ymin = apply(cuminc$cumincid[idx.est,"i45_15",],1,quantile,0.025),
                     ymax = apply(cuminc$cumincid[idx.est,"i45_15",],1,quantile,0.975))
  out <- ggplot(data = data, aes(y=y,x=x)) +
    geom_line(col="seagreen",size=1.2)+
    geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),alpha=.2,col=NA,fill="seagreen")+
    theme_bw()+
    xlab("") +
    ylab("Cumulative Incidence (15-60y/o)")+
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=ylim)
  out
}

plot_cuminc_empir <- function(cuminc, stand, main, ylim=c(0, 1), xlim=c(1985, 2015),
                              emp){
  ##
  data <- data.frame(x = as.numeric(rownames(cuminc)),
                     y = apply(cuminc,1,median),
                     ymin = apply(cuminc,1,quantile,0.025),
                     ymax = apply(cuminc,1,quantile,0.975))
  emp <- emp[is.finite(emp$se),]
  out <- ggplot() +
    geom_line(data = data, aes(y=y,x=x),col="seagreen",size=1.2)+
    geom_ribbon(data=data,aes(x=x,ymin=ymin,ymax=ymax),alpha=.2,col=NA,fill="seagreen")+
    theme_bw()+
    geom_point(data=emp,aes(x=year,y=1-coef)) +
    geom_linerange(data=emp,aes(x=year,ymin=1-coef-1.96*se,
                                ymax=1-coef+1.96*se)) +
    xlab("") +
    ylab("Cumulative Incidence (15-54y/o)")+
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim=ylim)
  out
}

plot.prev_update <- function(prev, stand, main, ylim=c(0, 1), xlim=c(1985, 2015)){
  ##
  data <- data.frame(x = as.numeric(rownames(prev)),
                     y = apply(prev,1,median),
                     ymin = apply(prev,1,quantile,0.025),
                     ymax = apply(prev,1,quantile,0.975))
  out <- ggplot(data = data, aes(y=y,x=x)) +
    geom_line(col="seagreen",size=1.2)+
    geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),alpha=.2,col=NA,fill="seagreen")+
    theme_bw()+
    xlab("") +
    ylab("Prevalence (15-55y/o)")+
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=ylim)
  out
}

plot.prev_two_emilio <- function(prevm, prevf, main, ylim=c(0, 1), xlim=c(1985, 2015),
                          artyear,legendonly=FALSE){
  ##
  data <- data.frame(x = as.numeric(rownames(prevm)),
                     ym = apply(prevm,1,median),
                     yminm = apply(prevm,1,quantile,0.025),
                     ymaxm = apply(prevm,1,quantile,0.975),
                     yf = apply(prevf,1,median),
                     yminf = apply(prevf,1,quantile,0.025),
                     ymaxf = apply(prevf,1,quantile,0.975))
  out <- ggplot(mapping=aes(x=x)) +
    geom_line(data = data,aes(y=ym,col="Men",lty="Men"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminm,ymax=ymaxm, fill="Men"),alpha=.2,col=NA)+
    geom_line(data = data,aes(y=yf,col="Women",lty="Men"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminf,ymax=ymaxf,fill="Women"),alpha=.2,col=NA)+
    geom_vline(aes(xintercept=artyear,lty="ART start",color="ART start",fill="ART start"),size=1.2,
               show.legend = FALSE) +
    theme_bw()+
    xlab("") +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=ylim) +
    theme(legend.position="bottom")+
    ylab("Incidence (15-54y/o)")
  if(legendonly==FALSE){
    out <- out +
      scale_fill_manual("",values = c("#f1a340", "#998ec3","#FFFFFFFF"), limits=c("Men","Women","ART start")) +
      scale_color_manual("",values = c("#f1a340", "#998ec3","grey30"), limits=c("Men","Women","ART start")) +
      scale_linetype_manual("",values = c("solid","solid","twodash"),limits=c("Men","Women","ART start")) +
      theme(legend.position="none")
  } else {
    out <- out +
      scale_fill_manual("",values = c("#f1a340", "#998ec3","#FFFFFFFF"), limits=c("Men","Women","ART start")) +
      scale_color_manual("",values = c("#f1a340", "#998ec3","black"), limits=c("Men","Women","ART start")) +
      scale_linetype_manual("",values = c("solid","solid","twodash"),limits=c("Men","Women","ART start")) +
      theme(legend.position="bottom",legend.title=element_text(size=20),
                       legend.text=element_text(size=15))
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

# Relative cumulative incidence plot

f_test <- Vectorize(function(x) {
  if(x>0) return(exp(x*3.4-1.7))
  x + 0
})

plot_cuminc_ratio <- function(datm,datf,main,xlim = c(2000,2017),ylim=c(0,1),
                              legendonly=FALSE) {
  dat <- datf/datm
  meanr <- apply(dat,1,function(x) exp(mean(log(x))))
  lowr <- apply(dat,1,quantile,0.025)
  highr <- apply(dat,1,quantile,0.975)
  data <- data.frame(x = as.numeric(rownames(datm)),
                     ym = apply(datm,1,median),
                     yminm = apply(datm,1,quantile,0.025),
                     ymaxm = apply(datm,1,quantile,0.975),
                     yf = apply(datf,1,median),
                     yminf = apply(datf,1,quantile,0.025),
                     ymaxf = apply(datf,1,quantile,0.975),
                     rrm = meanr,
                     rrmin = lowr,
                     rrmax = highr)
  out <- ggplot(mapping=aes(x=x)) +
    geom_hline(yintercept=0.5,size=1,col="grey30",alpha=0.7,lty="dotted") +
    geom_line(data = data,aes(y=ym,col="Men",lty="Men"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminm,ymax=ymaxm, fill="Men"),alpha=.2,col=NA)+
    geom_line(data = data,aes(y=yf,col="Women",lty="Men"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminf,ymax=ymaxf,fill="Women"),alpha=.2,col=NA)+
    geom_line(data=data,aes(y=(log(rrm)+1.7)/3.4,col="W/M ratio",lty="W/M ratio"),size=1.2) +
    geom_ribbon(data=data,aes(ymin=(log(rrmin)+1.7)/3.4,ymax=(log(rrmax)+1.7)/3.4,fill="W/M ratio"),alpha=0.2,col=NA) +
    theme_bw()+
    xlab("") +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=ylim,
                       sec.axis=sec_axis(~f_test(.), name="Sex ratio (Women/Men)", breaks=c(0.2,0.5,1,2,4))) +
    theme(legend.position="bottom")+
    ylab("Cumulative Incidence (15-54y/o)") +
    scale_fill_manual("",values = c("#f1a340", "#998ec3","#1b9e77"), limits=c("Men","Women","W/M ratio")) +
    scale_color_manual("",values = c("#f1a340", "#998ec3","#1b9e77"), limits=c("Men","Women","W/M ratio")) +
    scale_linetype_manual("",values = c("solid","solid","twodash"),limits=c("Men","Women","W/M ratio"))
  if(legendonly==FALSE) {
    out <- out + theme(legend.position="none")
  } else {
    out <- out +
      theme(legend.position="bottom",legend.title=element_text(size=20),
            legend.text=element_text(size=15))
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  return(out)
}

plot_cuminc_ratio2 <- function(datm,datf,main,xlim = c(2000,2017),ylim=c(0.2,5),
                              legendonly=FALSE) {
  dat <- datf/datm
  meanr <- apply(dat,1,function(x) exp(mean(log(x))))
  lowr <- apply(dat,1,quantile,0.025)
  highr <- apply(dat,1,quantile,0.975)
  data <- data.frame(x = as.numeric(rownames(datm)),
                     ym = apply(datm,1,median),
                     yminm = apply(datm,1,quantile,0.025),
                     ymaxm = apply(datm,1,quantile,0.975),
                     yf = apply(datf,1,median),
                     yminf = apply(datf,1,quantile,0.025),
                     ymaxf = apply(datf,1,quantile,0.975),
                     rrm = meanr,
                     rrmin = lowr,
                     rrmax = highr)
  out <- ggplot(mapping=aes(x=x)) +
    geom_hline(yintercept=1,size=1,col="grey30",alpha=0.7,lty="dotted") +
    geom_line(data=data,aes(y=rrm),size=1.2, col="seagreen") +
    geom_ribbon(data=data,aes(ymin=rrmin,ymax=rrmax),alpha=0.2,col=NA,fill="seagreen") +
    theme_bw()+
    xlab("") +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_log10(expand=c(0,0),limits=ylim) +
    theme(legend.position="bottom")+
    ylab("Sex Ratio of Cumulative Incidence (W/M)")
  return(out)
}

plot_prop_emilio2 <- function(dat,main,xlim = c(2000,2017),ages = "15-24") {
  test <- dat[,"15-24",]
  meanp <- apply(test,1,mean)
  lowp <- apply(test,1,quantile,0.025)
  highp <- apply(test,1,quantile,0.975)
  if(ages=="15-24") {
    plotdat <- data.frame(x = as.numeric(names(meanp)),
                          rrm = meanp,
                          rrmin = lowp,
                          rrmax = highp)
  } else {
    plotdat <- data.frame(x = as.numeric(names(meanp)),
                          rrm = 1- meanp,
                          rrmin = 1- lowp,
                          rrmax = 1- highp)
  }
  out <- ggplot(data=plotdat,mapping=aes(x=x)) +
    geom_line(aes(y=rrm),size=1.2, col="seagreen") +
    geom_ribbon(aes(ymin=rrmin,ymax=rrmax),alpha=0.2,col=NA,fill="seagreen") +
    theme_bw()+
    xlab("") +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=c(0,1)) +
    theme(legend.position="bottom")
  if(ages=="15-24") {
    out <- out + ylab("Prop of inf in 15-24 y/o")
  } else {
    out <- out + ylab("Prop of inf in 25-54 y/o")
  }
  return(out)
}

plot.prev_two <- function(prevm, prevf, main, ylim=c(0, 1), xlim=c(1985, 2015),legendonly=FALSE,
                          incprev = "prev"){
  ##
  data <- data.frame(x = as.numeric(rownames(prevm)),
                     ym = apply(prevm,1,median),
                     yminm = apply(prevm,1,quantile,0.025),
                     ymaxm = apply(prevm,1,quantile,0.975),
                     yf = apply(prevf,1,median),
                     yminf = apply(prevf,1,quantile,0.025),
                     ymaxf = apply(prevf,1,quantile,0.975))
  out <- ggplot(mapping=aes(x=x)) +
    geom_line(data = data,aes(y=ym,col="Men"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminm,ymax=ymaxm, fill="Men"),alpha=.2,col=NA)+
    geom_line(data = data,aes(y=yf,col="Women"),size=1.2)+
    geom_ribbon(data = data,aes(ymin=yminf,ymax=ymaxf,fill="Women"),alpha=.2,col=NA)+
    theme_bw()+
    xlab("") +
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    scale_y_continuous(expand=c(0,0),limits=ylim) +
    scale_color_manual("",values = c("#f1a340", "#998ec3"), limits=c("Men","Women")) +
    scale_fill_manual("",values = c("#f1a340", "#998ec3"), limits=c("Men","Women")) +
    theme(legend.position="bottom")
  if(incprev=="prev"){
    out <- out + ylab("Prevalence (15-54y/o)")
  }
  if(incprev == "inc") {
    out <- out + ylab("Incidence (15-54y/o)")
  }
  if(incprev == "cuminc") {
    out <- out + ylab("Cumulative Incidence (15-54y/o)")
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

plot.inc_update <- function(prev, stand, main, ylim=c(0, 1), xlim=c(1985, 2015)){
  ##
  data <- data.frame(x = as.numeric(rownames(prev)),
                     y = apply(prev,1,median),
                     ymin = apply(prev,1,quantile,0.025),
                     ymax = apply(prev,1,quantile,0.975))
  if(main=="Karonga men"){
    data <- subset(data,x<2015)
  }
  out <- ggplot(data = data, aes(y=y,x=x)) +
    geom_line(col="seagreen",size=1.2)+
    geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),alpha=.2,col=NA,fill="seagreen")+
    theme_bw()+
    xlab("") +
    ylab("Incidence (15-54y/o)")+
    # scale_x_continuous(expand=c(0,0)) +
    ggtitle(main) +
    # scale_y_continuous(expand=c(0,0))
    coord_cartesian(xlim=xlim,ylim=ylim,expand=FALSE)
  out
}

# library('ggplot2')

agespec_inc_dat <- function(agespec,xlim=c(15,55),years=c(2000,2015),offset){
  meds <- apply(agespec,1:2,mean)
  y <- c(meds)*100
  year <- rep(rownames(meds),length(colnames(meds)))
  age <- as.numeric(sapply(colnames(meds),rep,length(rownames(meds))))
  ymin <- c(apply(agespec,1:2,quantile,0.025))*100
  ymax <- c(apply(agespec,1:2,quantile,0.975))*100
  data <- data.frame(y = y, age = age, year = year, ymin = ymin, ymax = ymax)
  data <- data[as.numeric(as.character(data$year))<=max(years) &
                 as.numeric(as.character(data$year))>=min(years) &
                 data$age>=min(xlim) & data$age<=max(xlim),]
  data <- data[(as.numeric(as.character(data$year))-offset) %% 5==0,]

  return(data)
}

plot_agespec_inc <- function(agespec, main, ylim=NA, #xlim=c(15,55),
                             years = c(2000,2015),rib=TRUE, legendonly=FALSE,
                             offset=0){
  n <- (years[2]-years[1])/5+1
  hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
  names(hues) <- as.character(seq(years[1],years[2],5))
  data <- agespec_inc_dat(agespec,years=years,offset=offset)
  if(substr(main,1,16)=="uMkhanyakude men" | substr(main,1,18)=="uMkhanyakude women") {
    data <- subset(data,year!="2000")
  }
  if(substr(main,1,11)=="Karonga men" | substr(main,1,13)=="Karonga women") {
    data <- subset(data,year!="2000")
    data <- subset(data,year!="2015")
  }
  out <- ggplot(data=data) +
    geom_line(aes(x=age,y=y,col=year,linetype=year),size=1.2)+
    theme_bw()+
    xlab("Age") +
    ylab("Incidence (per 100 py)")+
    scale_x_continuous(expand=c(0,0)) +
    scale_color_manual("Year",values=hues)+ scale_fill_manual("Year",values=hues)+
    ggtitle(main) +
    theme(legend.position="bottom")+
    scale_linetype_manual("Year",values = c("solid","solid","solid","solid"))
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=age,ymin=ymin,ymax=ymax,fill=year),alpha=.1,col=NA)
  }
  # if(substr(main,1,10)=="Masaka men" | substr(main,1,12)=="Masaka women" |
  #    substr(main,1,11)=="Rakai women") {
  #   out <- out +
  #     scale_linetype_manual(values = c("solid","solid","dashed","solid"))
  # }
  # if(substr(main,1,18)=="uMkhanyakude women") {
  #   out <- out +
  #     scale_linetype_manual(values=c("solid","dashed","solid","solid"))
  # }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

agespec_inc_dat_overtime <- function(agespec,years=c(2000,2017),offset){
  meds <- c(apply(agespec$inc1524,1,mean),apply(agespec$inc2534,1,mean),
            apply(agespec$inc3544,1,mean),apply(agespec$inc4554,1,mean))
  y <- c(meds)*100
  year <- as.numeric(names(y))
  age <- c(rep("15-24",nrow(agespec$inc1524)),rep("25-34",nrow(agespec$inc2534)),
           rep("35-44",nrow(agespec$inc3544)),rep("45-54",nrow(agespec$inc4554)))
  ymin <- c(apply(agespec$inc1524,1,quantile,0.025),apply(agespec$inc2534,1,quantile,0.025),
            apply(agespec$inc3544,1,quantile,0.025),apply(agespec$inc4554,1,quantile,0.025))*100
  ymax <- c(apply(agespec$inc1524,1,quantile,0.975),apply(agespec$inc2534,1,quantile,0.975),
            apply(agespec$inc3544,1,quantile,0.975),apply(agespec$inc4554,1,quantile,0.975))*100
  data <- data.frame(y = y, age = age, year = year, ymin = ymin, ymax = ymax)
  data <- filter(data,year>=years[1] & year<=years[2])
  return(data)
}

plot_agespec_inc_overtime <- function(agespec, main, ylim=NA, #xlim=c(15,55),
                             years = c(2000,2017),rib=TRUE, legendonly=FALSE,
                             offset=0){
  n <- 4
  hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
  names(hues) <- c("15-24","25-34","35-44","45-54")
  data <- agespec_inc_dat_overtime(agespec,years=years,offset=offset)

  if(substr(main,1,11)=="Karonga men" | substr(main,1,13)=="Karonga women" |
     substr(main,1,14)=="Manicaland men" | substr(main,1,16)=="Manicaland women") {
    data <- subset(data,year<=2012)
  }
  if(substr(main,1,9)=="Rakai men" | substr(main,1,11)=="Rakai women") {
    data <- subset(data,age!="45-54")
  }
  out <- ggplot(data=data) +
    geom_line(aes(x=year,y=y,col=age,linetype=age),size=1.2)+
    theme_bw()+
    xlab("Year") +
    ylab("Incidence (per 100 py)")+
    scale_x_continuous(expand=c(0,0),limits=years) +
    scale_color_manual("Age",values=hues)+ scale_fill_manual("Age",values=hues)+
    ggtitle(main) +
    theme(legend.position="bottom") +
    scale_linetype_manual("Age",values = c("solid","solid","solid","solid"))
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax,fill=age),alpha=.15,col=NA)
  }
  # if(substr(main,1,13)=="Karonga women" | substr(main,1,12)=="Kisesa women") {
  #   out <- out +
  #     scale_linetype_manual(values = c("solid","solid","dashed","solid"))
  # }
  # if(substr(main,1,10)=="Kisesa men") {
  #   out <- out +
  #     scale_linetype_manual(values = c("solid","solid","dashed","dashed"))
  # }
  # if(substr(main,1,10)=="Masaka men") {
  #   out <- out +
  #     scale_linetype_manual(values = c("solid","dashed","solid","solid"))
  # }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

agespec_prev_dat_overtime <- function(agespec,years=c(2000,2017),offset){
  meds <- c(apply(agespec$prev1524,1,mean),apply(agespec$prev2534,1,mean),
            apply(agespec$prev3544,1,mean),apply(agespec$prev4554,1,mean))
  y <- c(meds)
  year <- as.numeric(names(y))
  age <- c(rep("15-24",nrow(agespec$prev1524)),rep("25-34",nrow(agespec$prev2534)),
           rep("35-44",nrow(agespec$prev3544)),rep("45-54",nrow(agespec$prev4554)))
  ymin <- c(apply(agespec$prev1524,1,quantile,0.025),apply(agespec$prev2534,1,quantile,0.025),
            apply(agespec$prev3544,1,quantile,0.025),apply(agespec$prev4554,1,quantile,0.025))
  ymax <- c(apply(agespec$prev1524,1,quantile,0.975),apply(agespec$prev2534,1,quantile,0.975),
            apply(agespec$prev3544,1,quantile,0.975),apply(agespec$prev4554,1,quantile,0.975))
  data <- data.frame(y = y, age = age, year = year, ymin = ymin, ymax = ymax)
  data <- filter(data,year>=years[1] & year<=years[2])
  return(data)
}

plot_agespec_prev_overtime <- function(agespec, main, ylim=NA, #xlim=c(15,55),
                                      years = c(2000,2017),rib=TRUE, legendonly=FALSE,
                                      offset=0){
  n <- 4
  hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
  names(hues) <- c("15-24","25-34","35-44","45-54")
  data <- agespec_prev_dat_overtime(agespec,years=years,offset=offset)

  if(main=="Karonga men" | main=="Karonga women" | main=="Manicaland men" | main=="Manicaland women") {
    data <- subset(data,year<=2012)
  }
  if(main=="Rakai men" | main=="Rakai women") {
    data <- subset(data,age!="45-54")
  }
  out <- ggplot(data=data) +
    geom_line(aes(x=year,y=y,col=age),size=1.2)+
    theme_bw()+
    xlab("Year") +
    ylab("Prevalence")+
    scale_x_continuous(expand=c(0,0),limits=years) +
    scale_color_manual("Age",values=hues)+ scale_fill_manual("Age",values=hues)+
    ggtitle(main) +
    theme(legend.position="bottom")
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax,fill=age),alpha=.15,col=NA)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

plot_avgageinfect <- function(avgageinf, main, ylim=NA, xlim=c(1985,2017),
                              rib=TRUE){
  y <- apply(avgageinf,1,median)
  year <- as.numeric(names(y))
  ymin <- apply(avgageinf,1,quantile,0.025)
  ymax <- apply(avgageinf,1,quantile,0.975)
  data <- data.frame(y = y, year = year, ymin = ymin, ymax = ymax)
  data <- data[data$year<=max(xlim) &
                 data$year>=min(xlim),]
  out <- ggplot(data=data) +
    geom_line(aes(x=year,y=y),col="seagreen",size=1.2)+
    theme_bw()+
    xlab("Year") +
    ylab("Average age of infection")+
    scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
    ggtitle(main) +
    theme(legend.position="bottom")
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax),alpha=.5,col=NA,fill="seagreen")
  }
  out
}

plot_avgageinfect_both <- function(avgageinf,avgageinfunw, main, ylim=NA, xlim=c(1985,2017),
                                   rib=TRUE,legendonly=FALSE,whichtype = "sex"){
  y <- c(apply(avgageinf,1,median),apply(avgageinfunw,1,median))
  year <- as.numeric(names(y))
  ymin <- c(apply(avgageinf,1,quantile,0.025),apply(avgageinfunw,1,quantile,0.025))
  ymax <- c(apply(avgageinf,1,quantile,0.975),apply(avgageinfunw,1,quantile,0.975))
  if(whichtype=="sex") {
    type <- c(rep("Men",nrow(avgageinf)),rep("Women",nrow(avgageinfunw)))
  } else {
    type <- c(rep("weighted",nrow(avgageinf)),rep("unweighted",nrow(avgageinfunw)))
  }
  data <- data.frame(y = y, year = year, ymin = ymin, ymax = ymax,type=type)
  data <- data[data$year<=max(xlim) &
                 data$year>=min(xlim),]
  out <- ggplot(data=data) +
    geom_line(aes(x=year,y=y,col=type),size=1.2)+
    theme_bw()+
    xlab("Year") +
    ylab("Average age of infection")+
    scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
    ggtitle(main) +
    theme(legend.position="bottom")+
    scale_color_brewer("",palette="Set1",direction=-1) + scale_fill_brewer("",palette="Set1",direction=-1)
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax,fill=type),alpha=.2,col=NA)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

plot_avgageinfect_four <- function(avgageinf,avgageinfunw,avgageinf2,avgageinfunw2, main, ylim=NA, xlim=c(1985,2017),
                                   rib=TRUE,legendonly=FALSE,whichtype = "sex"){
  library(RColorBrewer)
  y <- c(apply(avgageinf,1,median),apply(avgageinfunw,1,median),
         apply(avgageinf2,1,median),apply(avgageinfunw2,1,median))
  year <- as.numeric(names(y))
  ymin <- c(apply(avgageinf,1,quantile,0.025),apply(avgageinfunw,1,quantile,0.025),
            apply(avgageinf2,1,quantile,0.025),apply(avgageinfunw2,1,quantile,0.025))
  ymax <- c(apply(avgageinf,1,quantile,0.975),apply(avgageinfunw,1,quantile,0.975),
            apply(avgageinf2,1,quantile,0.975),apply(avgageinfunw2,1,quantile,0.975))
  if(whichtype=="sex") {
    type <- c(rep("Men",nrow(avgageinf)),rep("Women",nrow(avgageinfunw)),
              rep("Men",nrow(avgageinf2)),rep("Women",nrow(avgageinfunw2)))
    esttype <- c(rep("Natl",2*nrow(avgageinf)),rep("Study",2*nrow(avgageinf2)))
  } else {
    type <- c(rep("weighted",nrow(avgageinf)),rep("unweighted",nrow(avgageinfunw)))
  }
  data <- data.frame(y = y, year = year, ymin = ymin, ymax = ymax,type=type,esttype=esttype)
  data$ymin <- ifelse(data$esttype=="Study",NA,data$ymin)
  data$ymax <- ifelse(data$esttype=="Study",NA,data$ymax)
  data <- data[data$year<=max(xlim) &
                 data$year>=min(xlim),]
  if(legendonly==FALSE){
    data <- mutate(data, sex = case_when(type=="Men" & esttype=="Natl" ~ "Men",
                                         type=="Women" & esttype == "Natl" ~ "Women",
                                            type=="Men" & esttype=="Study" ~ "Men2",
                                            type=="Women" & esttype == "Study" ~ "Women2",
                                            TRUE ~ NA_character_))
    data$sex <- factor(data$sex)
    cols <- brewer.pal(2, "Set1")[1:2]
    out <- ggplot(data=data) +
      geom_line(aes(x=year,y=y,col=sex,lty=esttype),size=1.2)+
      theme_bw()+
      xlab("Year") +
      ylab("Average age of infection")+
      scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
      ggtitle(main) +
      theme(legend.position="bottom")+
      scale_color_manual("",values = c(cols[2],cols[2],cols[1],cols[1])) +
      scale_fill_manual("",values = c(cols[2],cols[2],cols[1],cols[1])) +
      scale_linetype_discrete("")
    if(is.na(ylim[1])) {
      out <- out + scale_y_continuous(expand=c(0,0))
    } else if(!is.na(ylim[1])) {
      out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
    }
    if(rib==TRUE) {
      out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax,fill=sex),alpha=.2,col=NA)
    }
    out <- out + theme(legend.position="none")
  } else {
    out <- ggplot(data=data) +
      geom_line(aes(x=year,y=y,col=type,lty=esttype),size=1.2)+
      theme_bw()+
      xlab("Year") +
      ylab("Average age of infection")+
      scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
      ggtitle(main) +
      theme(legend.position="bottom")+
      scale_color_brewer("",palette="Set1",direction=-1) + scale_fill_brewer("",palette="Set1",direction=-1) +
      scale_linetype_discrete("")
    if(is.na(ylim[1])) {
      out <- out + scale_y_continuous(expand=c(0,0))
    } else if(!is.na(ylim[1])) {
      out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
    }
    if(rib==TRUE) {
      out <- out + geom_ribbon(aes(x=year,ymin=ymin,ymax=ymax,fill=type),alpha=.2,col=NA)
    }
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

irr_fx <- function(incid,refgr,xlim=c(15,55)){
  incid <- incid[rownames(incid)>=refgr,colnames(incid)>=xlim[1] & colnames(incid)<=xlim[2],]
  refyr <- incid[rownames(incid)==as.character(refgr),,,drop=FALSE]
  irr <- array(dim=dim(incid),dimnames=dimnames(incid))
  for(i in 1:nrow(incid)) {
    irr[i,,] <- incid[i,,,drop=FALSE]/refyr
  }
  irr <- irr[rownames(irr)>refgr,,,drop=FALSE]
  meds <- apply(irr,1:2,median)
  y <- c(meds)
  year <- rep(rownames(meds),length(colnames(meds)))
  age <- as.numeric(sapply(colnames(meds),rep,length(rownames(meds))))
  ymin <- c(apply(irr,1:2,quantile,0.025))
  ymax <- c(apply(irr,1:2,quantile,0.975))
  data <- data.frame(y = y, age = age, year = year, ymin = ymin, ymax = ymax)
  return(data)
}

plot_irr <- function(incid, main, ylim=NA, xlim=c(15,55),
                     refgr=2000,rib=TRUE, legendonly=FALSE){
  data <- irr_fx(incid,refgr,xlim)
  n <- 3
  hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
  names(hues) <- c("2005","2010","2015")
  out <- ggplot(data=data) +
    geom_line(aes(x=age,y=y,col=year),size=1.5)+
    theme_bw()+
    geom_hline(yintercept=1)+
    xlab("Age") +
    ylab(paste0("IRR (relative to ",refgr,")"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_color_manual("Year",values=hues)+ scale_fill_manual("Year",values=hues)+
    ggtitle(main) +
    theme(legend.position="bottom")
  if(is.na(ylim[1])) {
    out <- out + scale_y_log10(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_log10(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=age,ymin=ymin,ymax=ymax,fill=year),alpha=.15,col=NA)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

irr_fx_group <- function(incid,refgr,xlim=c(15,55),yrdiff){
  incid <- incid[rownames(incid)>=refgr,colnames(incid)>=xlim[1] & colnames(incid)<xlim[2],]
  dt <- as.numeric(colnames(incid)[2]) - as.numeric(colnames(incid)[1])
  output <- apply(incid,3,function(x) {
    sapply(seq(1,ncol(x),yrdiff/dt), function(i) {
      rowMeans(x[,c(i:(i+5/dt-1))], na.rm=T)})
  })
  updated_output <- array(output,c(nrow(incid),(xlim[2]-xlim[1])/yrdiff,length(incid[1,1,])))
  colnames(updated_output) <- paste0(seq(xlim[1],xlim[2]-yrdiff+1,yrdiff),"-",seq(xlim[1]+yrdiff-1,xlim[2]-1,yrdiff))
  rownames(updated_output) <- rownames(incid)

  refyr <- updated_output[rownames(updated_output)==as.character(refgr),,,drop=FALSE]
  irr <- array(dim=dim(updated_output),dimnames=dimnames(updated_output))
  for(i in 1:nrow(updated_output)) {
    irr[i,,] <- updated_output[i,,,drop=FALSE]/refyr
  }
  irr <- irr[rownames(irr)>refgr,,,drop=FALSE]
  #irr <- t(irr)
  meds <- apply(irr,1:2,median)
  y <- c(meds)
  year <- rep(rownames(meds),length(colnames(meds)))
  age <- c(sapply(colnames(meds),rep,length(rownames(meds))))
  ymin <- c(apply(irr,1:2,quantile,0.025))
  ymax <- c(apply(irr,1:2,quantile,0.975))
  data <- data.frame(y = y, age = age, year = year, ymin = ymin, ymax = ymax)
  return(data)
}

plot_irr_grouped <- function(incid, main, ylim=NA, xlim=c(15,55),
                             refgr=2000, legendonly=FALSE,yrdiff=5){
  data <- irr_fx_group(incid,refgr,xlim,yrdiff)
  if(refgr==2000) {
    n <- 3
    hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
    names(hues) <- c("2005","2010","2015")
  }
  if(refgr==2005){
    n <- 2
    hues <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
    names(hues) <- c("2010","2015")
  }

  out <- ggplot(data=data,aes(x=age,y=y,col=year)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=ymin,ymax=ymax),position=position_dodge(width=0.5),
                  width=0.5) +
    theme_bw()+
    geom_hline(yintercept=1) +
    scale_color_manual("Year",values=hues)+ scale_fill_manual("Year",values=hues)+
    xlab("Age") +
    ylab(paste0("IRR (relative to ",refgr,")"))+
    ggtitle(main) +
    theme(legend.position="bottom")
  if(is.na(ylim[1])) {
    out <- out + scale_y_log10(expand=c(0,0.01))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_log10(expand=c(0,0.01),limits=ylim)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

# data_for_infect_by_age <- function(input,yearbase = 1990) {
#   fulldat <- array(unlist(input), c(nrow(input[[1]]), ncol(input[[1]]), length(input)))
#   meds <- apply(fulldat,1:2,median)
#   meds <- meds/rowSums(meds)
#   y <- c(meds)
#   yearstot <- seq(yearbase,yearbase + (nrow(fulldat)-1)/5,by=.2)
#   year <- rep(yearstot,ncol(meds))
#   age <- c(sapply(colnames(input[[1]]),rep,nrow(meds)))
#   data <- data.frame(y = y, age = age, year = year)
#   return(data)
# }

data_for_infect_by_age <- function(fulldat) {
  meds <- apply(fulldat,1:2,median)
  meds <- meds/rowSums(meds)
  y <- c(meds)
  yearstot <- as.numeric(rownames(fulldat))
  year <- rep(yearstot,ncol(meds))
  age <- c(sapply(colnames(fulldat),rep,nrow(meds)))
  data <- data.frame(y = y, age = age, year = year)
  return(data)
}


plot_infect_by_age <- function(input,yearbase = 1990, main,
                               legendonly=FALSE, xlim=c(1990,2015)) {
  data <- data_for_infect_by_age(input)
  if(n_distinct(data$age)==4){
    data$age <- data$age %>% factor(levels= c("45-54","35-44","25-34","15-24"))
  }
  if(n_distinct(data$age)==8) {
    data$age <- data$age %>% factor(levels= c("50-54","45-49","40-44","35-39","30-34","25-29","20-24","15-19"))
  }
  if(n_distinct(data$age)==11) {
    data$age <- data$age %>% factor(levels= c("60-64","55-59","50-54","45-49","40-44","35-39","30-34","25-29","20-24","15-19","10-14"))
  }
  if(main=="Manicaland men" | main=="Manicaland women") {
    data <- subset(data,year<2012)
  }
  out <- ggplot(data, aes(x = year, y = y, fill = age)) +
    geom_area(position = 'stack',alpha=.8) +
    theme_bw()+
    xlab("Year") +
    ylab("Proportion") +
    scale_fill_brewer("Age",palette="Paired") +
    ggtitle(main) +
    theme(legend.position="bottom") +
    scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(reverse=T))
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

plot_infect_by_age_emilio <- function(input,yearbase = 1990, main,
                               legendonly=FALSE, xlim=c(1990,2015)) {
  data <- data_for_infect_by_age(input)
  if(nlevels(data$age)==4){
    data$age <- data$age %>% factor(levels= c("45-54","35-44","25-34","15-24"))
  }
  if(nlevels(data$age)==8) {
    data$age <- data$age %>% factor(levels= c("50-54","45-49","40-44","35-39","30-34","25-29","20-24","15-19"))
  }
  if(nlevels(data$age)==11) {
    data$age <- data$age %>% factor(levels= c("60-64","55-59","50-54","45-49","40-44","35-39","30-34","25-29","20-24","15-19","10-14"))
  }
  if(main=="Manicaland men" | main=="Manicaland women") {
    data <- subset(data,year<2012)
  }
  out <- ggplot(data, aes(x = year, y = y, fill = age)) +
    geom_area(position = 'stack',alpha=.8) +
    theme_bw()+
    xlab("Year") +
    ylab("Proportion") +
    scale_fill_manual("Age",values = c(brewer.pal(7,"Blues")[2:7],
                                       brewer.pal(8,"Paired")[7:8])) +
    ggtitle(main) +
    theme(legend.position="bottom") +
    scale_x_continuous(expand=c(0,0),limits=xlim, breaks=seq((xlim[1]),(xlim[2]),by=5)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(reverse=T))
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}


numnewinf_plot <- function(data,main,ages,years,legendonly=FALSE,ylim){
  infs <- apply(data,c(1,2),median)
  infs <- infs[,colnames(infs)>=ages[1] & colnames(infs)<ages[2]]
  numinfs <- c(infs)
  age <- as.numeric(unlist(lapply(colnames(infs),rep,nrow(infs))))
  year <- as.numeric(rep(rownames(infs),ncol(infs)))
  dat <- data.frame(age = age, year = year, numinfs = numinfs)
  dat$agecat <- cut(dat$age,seq(ages[1],ages[2],by=5),right=FALSE,
                    labels=c(paste0(seq(ages[1],ages[2]-5,by=5),"-",seq(ages[1]+4,ages[2]-1,5))))
  dat$agecat <- factor(dat$agecat,levels=rev(levels(dat$agecat)))
  dat <- summarise(group_by(dat,year,agecat), tot = sum(numinfs))
  out <- ggplot(data=dat,aes(y=tot,x=year,fill=agecat)) +
    geom_area() +
    scale_x_continuous(limits=c(years[1],years[2]),expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,ylim)) +
    theme_bw() +
    ylab("Number of new infections") +
    xlab("") +
    scale_fill_discrete("Ages",breaks=rev(levels(dat$agecat))) +
    ggtitle(main) +
    theme(legend.position="bottom")
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

ptiles_plot <- function(data,main,yearlim,agelim,legendonly=FALSE) {
  data$ptiles <- as.factor(data$ptiles)
  data$ptiles <- factor(data$ptiles,levels=c(rev(levels(data$ptiles))[9],rev(levels(data$ptiles))[-9]))
  out <- ggplot(data,aes(y=ages,x=years,fill=ptiles)) +
    geom_area(position=position_identity()) +
    geom_line(data = filter(data,ptiles=="50%"),aes(y=ages,x=years,linetype="Median")) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),limits=c(yearlim[1],yearlim[2])) +
    coord_cartesian(ylim=c(agelim[1],agelim[2])) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(main) +
    scale_fill_manual("Percentiles",values = c("10%" = viridis(5)[1], "20%" = viridis(5)[2],
                                               "30%" = viridis(5)[3], "40%" = viridis(5)[4],
                                               "50%" = viridis(5)[5], "60%" = viridis(5)[5],
                                               "70%" = viridis(5)[4], "80%" = viridis(5)[3],
                                               "90%" = viridis(5)[2], "100%" = viridis(5)[1]),
                      breaks = c("10%","20%","30%","40%","50%"),
                      labels = c("0-100%","10-90%","20-80%","30-70%","40-60%")) +
    scale_linetype_discrete("") +
    ylab("Age") +
    xlab("Year") +
    theme(legend.position="bottom")
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

data_for_agebands <- function(x) {
  years <- rep(as.numeric(rownames(x)),10)
  ages <- c(x[,!colnames(x) %in% c("0%","5%","25%","75%","95%")])
  ptiles <- c(sapply(colnames(x[,!colnames(x) %in% c("0%","5%","25%","75%","95%")]),rep,length(years)/10))
  data <- data.frame(years = years, ages = ages, ptiles=ptiles)
  return(data)
}

data2_for_infect_by_age <- function(input,sexname,nnas) {
  # input <- lapply(input,function(x) x/rowSums(x))
  # fulldat <- array(unlist(input), c(nrow(input[[1]]), ncol(input[[1]]), length(input)))
  fulldat <- input
  lastyr <- fulldat[nrow(fulldat[,,1])-nnas,,,drop=FALSE]
  ## is this sensible??
  meds <- apply(lastyr,1:2,median)
  y <- c(meds)
  ymin <- c(apply(lastyr,1:2,quantile,0.025))
  ymax <- c(apply(lastyr,1:2,quantile,0.975))
  age <- colnames(input)
  sex <- rep(sexname,length(y))
  data <- data.frame(y = y, ymin=ymin, ymax = ymax, age = age, sex = sex)
  return(data)
}

plot_infect_prop_mr <- function(inputf,inputm,yearbase = 1990, main,
                                legendonly=FALSE,nnas=2,ylim=c(0,1)) {
  # if(main=="Manicaland"){
  #   inputf <- lapply(inputf,function(x) x[1:((2012-yearbase)*5),])
  #   inputm <- lapply(inputm,function(x) x[1:((2012-yearbase)*5),])
  #   nnas <- 0
  # }
  dataf <- data2_for_infect_by_age(inputf,"Female",nnas)
  datam <- data2_for_infect_by_age(inputm,"Male",nnas)
  data <- rbind(dataf,datam)
  year <- round(yearbase + (nrow(inputf)-1-nnas)*.2)
  out <- ggplot(data, aes(x=age,y=y,fill=sex)) +
    geom_col(aes(col=sex),position=position_dodge(width=.92)) +
    geom_errorbar(aes(ymin=ymin,ymax=ymax),position=position_dodge(width=.92),width=.1,col="darkgrey") +
    ggtitle(paste0(main," (",year,")")) +
    theme_bw()+
    scale_y_continuous("Proportion",limits=ylim,expand=c(0,0))  +
    scale_x_discrete("Age")+
    # scale_colour_discrete("Sex") +
    # scale_fill_discrete("Sex")+
    theme(legend.position="bottom") +
    scale_color_brewer("Sex",palette="Set1") +
    scale_fill_brewer("Sex",palette="Set1")
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}


ptiles_plot2 <- function(data,main,yearlim,agelim,legendonly=FALSE) {
  data$ptiles <- factor(data$quant, levels = c(rev(levels(factor(data$quant)))[7],
                                               rev(levels(factor(data$quant)))[-7]))
  data <- dplyr::filter(data,ptiles!="0")
  data$ages <- data$val
  data$years <- data$year
  out <- ggplot(data,aes(y=ages,x=years,fill=ptiles)) +
    geom_area(position=position_identity()) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),limits=c(yearlim[1],yearlim[2])) +
    coord_cartesian(ylim=c(agelim[1],agelim[2])) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(main) +
    scale_fill_manual("Percentiles",values = c("12.5" = viridis(4)[1], "25" = viridis(4)[2],
                                               "37.5" = viridis(4)[3], "62.5" = viridis(4)[4],
                                               "75" = viridis(4)[3], "87.5" = viridis(4)[2],
                                               "100" = viridis(4)[1]),
                      breaks = c("12.5","25","37.5","62.5"),
                      labels = c("All","75%","50%","25%")) +
    scale_linetype_discrete("") +
    ylab("Age") +
    xlab("Year") +
    theme(legend.position="bottom")
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}




# Cohort cumulative incidence plot
cumincid_coh_plot <- function(cuminc,stand,main,years = as.character(seq(1970,1995,by=5)),
                              rib = TRUE, ylim = c(0,1),legendonly=FALSE,xlim) {
  data <- data.frame(y = numeric(), ymin = numeric(), ymax = numeric(), age = numeric(), coh = character())
  for(i in 1:length(years)){
    curyear <- years[i]
    subdat <- array(unlist(lapply(cuminc,"[",curyear)),
                    dim=c(length(unlist(cuminc[[1]][curyear])),length(cuminc))) #each row is new age
    # if(main=="Karonga men") {
    #   subdat <- subdat[-c(nrow(subdat),nrow(subdat)-1,nrow(subdat)-2,
    #                       nrow(subdat)-3,nrow(subdat)-4),]
    # }
    y <- apply(subdat,1,median)
    ymin <- apply(subdat,1,quantile,0.025)
    ymax <- apply(subdat,1,quantile,0.975)
    age <- seq(15.1,length(y)*stand$dt+15.1 - stand$dt,by=stand$dt)
    year <- rep(curyear,length(y))
    data <- rbind(data,data.frame(y=y,ymin=ymin,ymax=ymax,age=age,coh=year))
  }
  out <- ggplot(data=data,aes(y=y,x=age,col=coh)) +
    geom_line(size=1.2) +
    theme_bw()+
    xlab("Age") +
    ylab("Cumulative Incidence")+
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    theme(legend.position="bottom") +
    scale_colour_discrete("Cohort Birth Year") +
    scale_fill_discrete("Cohort Birth Year")
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=age,ymin=ymin,ymax=ymax, fill=coh),alpha=.2,col=NA)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

# Cohort cumulative incidence plot Projection
cumincid_coh_plot_proj <- function(cuminc,stand,main,years = as.character(seq(1970,1995,by=5)),
                              rib = TRUE, ylim = c(0,1),legendonly=FALSE,xlim) {
  data <- data.frame(y = numeric(), ymin = numeric(), ymax = numeric(), age = numeric(), coh = character())
  for(i in 1:length(years)){
    curyear <- years[i]
    subdat <- array(unlist(lapply(cuminc,"[",curyear)),
                    dim=c(length(unlist(cuminc[[1]][curyear])),length(cuminc))) #each row is new age
    # if(main=="Karonga men") {
    #   subdat <- subdat[-c(nrow(subdat),nrow(subdat)-1,nrow(subdat)-2,
    #                       nrow(subdat)-3,nrow(subdat)-4),]
    # }
    y <- apply(subdat,1,median)
    ymin <- apply(subdat,1,quantile,0.025)
    ymax <- apply(subdat,1,quantile,0.975)
    age <- seq(15.1,length(y)*stand$dt+15.1 - stand$dt,by=stand$dt)
    year <- rep(curyear,length(y))
    data <- rbind(data,data.frame(y=y,ymin=ymin,ymax=ymax,age=age,coh=year))
  }
  yearchange <- max(stand$x_time)
  if(main %in% c("Manicaland men","Manicaland women","Karonga men","Karonga women")) {
    yearchange <- 2012
  }
  data$year <- data$age + as.numeric(as.character(data$coh))
  data$proj <- factor(ifelse(data$year>yearchange,"Projected","Observed"))
  out <- ggplot(data=data,aes(y=y,x=age,col=coh)) +
    geom_line(aes(lty=proj),size=1.2) +
    theme_bw()+
    xlab("Age") +
    ylab("Cumulative Incidence")+
    scale_x_continuous(expand=c(0,0),limits=xlim) +
    ggtitle(main) +
    theme(legend.position="bottom") +
    scale_colour_discrete("Cohort Birth Year") +
    scale_fill_discrete("Cohort Birth Year") +
    scale_linetype_discrete("")
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  if(rib==TRUE) {
    out <- out + geom_ribbon(aes(x=age,ymin=ymin,ymax=ymax, fill=coh),alpha=.2,col=NA)
  }
  if(legendonly==FALSE){
    out <- out + theme(legend.position="none")
  } else {
    tmp <- ggplot_gtable(ggplot_build(out))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  out
}

# # Cohort cumulative incidence plot Projection with a reduction
# cumincid_coh_plot_proj <- function(cuminc,stand,main,years = as.character(seq(1970,1995,by=5)),
#                                    rib = TRUE, ylim = c(0,1),legendonly=FALSE,xlim,incid_w) {
#   data <- data.frame(y = numeric(), ymin = numeric(), ymax = numeric(), age = numeric(), coh = character())
#   for(i in 1:length(years)){
#     curyear <- years[i]
#     subdat <- array(unlist(lapply(cuminc,"[",curyear)),
#                     dim=c(length(unlist(cuminc[[1]][curyear])),length(cuminc))) #each row is new age
#     # if(main=="Karonga men") {
#     #   subdat <- subdat[-c(nrow(subdat),nrow(subdat)-1,nrow(subdat)-2,
#     #                       nrow(subdat)-3,nrow(subdat)-4),]
#     # }
#     y <- apply(subdat,1,median)
#     ymin <- apply(subdat,1,quantile,0.025)
#     ymax <- apply(subdat,1,quantile,0.975)
#     age <- seq(15.1,length(y)*stand$dt+15.1 - stand$dt,by=stand$dt)
#     year <- rep(curyear,length(y))
#     data <- rbind(data,data.frame(y=y,ymin=ymin,ymax=ymax,age=age,coh=year))
#   }
#   yearchange <- max(stand$x_time)
#   if(main %in% c("Manicaland men","Manicaland women","Karonga men","Karonga women")) {
#     yearchange <- 2012
#   }
#   data$year <- data$age + as.numeric(as.character(data$coh))
#   data$proj <- factor(ifelse(data$year>yearchange,"Projected","Observed"))
#   out <- ggplot(data=data,aes(y=y,x=age,col=coh)) +
#     geom_line(aes(lty=proj),size=1.2) +
#     theme_bw()+
#     xlab("Age") +
#     ylab("Cumulative Incidence")+
#     scale_x_continuous(expand=c(0,0),limits=xlim) +
#     ggtitle(main) +
#     theme(legend.position="bottom") +
#     scale_colour_discrete("Cohort Birth Year") +
#     scale_fill_discrete("Cohort Birth Year") +
#     scale_linetype_discrete("")
#   if(is.na(ylim[1])) {
#     out <- out + scale_y_continuous(expand=c(0,0))
#   } else if(!is.na(ylim[1])) {
#     out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
#   }
#   if(rib==TRUE) {
#     out <- out + geom_ribbon(aes(x=age,ymin=ymin,ymax=ymax, fill=coh),alpha=.2,col=NA)
#   }
#   if(legendonly==FALSE){
#     out <- out + theme(legend.position="none")
#   } else {
#     tmp <- ggplot_gtable(ggplot_build(out))
#     leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
#     out <- cowplot::ggdraw() +
#       cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
#   }
#   out
# }



cumincid_proj_plot_2 <- function(datm,datf,xlim,ylim,main,legendonly=FALSE) {
  dat25 <- data.frame(med = c(apply(datm[,"cuminc25",],1,median),apply(datf[,"cuminc25",],1,median)),
                      lb = c(apply(datm[,"cuminc25",],1,quantile,.025),
                             apply(datf[,"cuminc25",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc25",],1,quantile,0.975),
                             apply(datf[,"cuminc25",],1,quantile,0.975)),
                      yob = c(datm[,"yob",1],datf[,"yob",1]),
                      proj = factor(c(datm[,"proj25",1],datf[,"proj25",1]),levels=c(0,1),
                                    labels=c("Observed","Projected")),
                      sex = c(rep("Male",nrow(datm)),rep("Female",nrow(datf))))
  dat35 <- data.frame(med = c(apply(datm[,"cuminc35",],1,median),apply(datf[,"cuminc35",],1,median)),
                      lb = c(apply(datm[,"cuminc35",],1,quantile,.025),
                             apply(datf[,"cuminc35",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc35",],1,quantile,0.975),
                             apply(datf[,"cuminc35",],1,quantile,0.975)),
                      yob = c(datm[,"yob",1],datf[,"yob",1]),
                      proj = factor(c(datm[,"proj35",1],datf[,"proj35",1]), levels=c(0,1),
                                    labels=c("Observed","Projected")),
                      sex = c(rep("Male",nrow(datm)),rep("Female",nrow(datf))))
  dat50 <- data.frame(med = c(apply(datm[,"cuminc50",],1,median),apply(datf[,"cuminc50",],1,median)),
                      lb = c(apply(datm[,"cuminc50",],1,quantile,.025),
                             apply(datf[,"cuminc50",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc50",],1,quantile,0.975),
                             apply(datf[,"cuminc50",],1,quantile,0.975)),
                      yob = c(datm[,"yob",1],datf[,"yob",1]),
                      proj = factor(c(datm[,"proj50",1],datf[,"proj50",1]),levels=c(0,1),
                                    labels=c("Observed","Projected")),
                      sex = c(rep("Male",nrow(datm)),rep("Female",nrow(datf))))
  fmts <- rep("atop(atop(textstyle(%s), textstyle(bold(%s))),NA)",6) # two rows
  labs2 <- do.call(sprintf,list(fmts,c("Y25","2000","2005","2010","2015","2020"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'")) )
  ex2 <- parse(text=labs2)
  outlab <- test <- parse(text="atop('Y25/35/50 = Year turned 25, 35 or 50',bold('YOB = Year of birth'))")
  p25 <- ggplot(data=dat25,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(alpha=0.2,col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(y = "Cum incidence age 25"
         # ,
         # title = "Cumulative incidence at age 25"
         )+
    scale_linetype_discrete("",drop=FALSE) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
    scale_color_brewer(outlab,palette="Set1") +
    scale_fill_brewer(outlab,palette="Set1") +
    guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))
  labs2 <- do.call(sprintf,list(fmts,c("Y35","2010","2015","2020","2025","2030"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p35 <- ggplot(data=dat35,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(alpha=0.2,col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
         y = "Cum incidence age 35"
         # ,title = "Cumulative incidence at age 35"
         )+
    scale_linetype_discrete(drop=FALSE) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_brewer("",palette="Set1") + scale_fill_brewer("",palette="Set1")
  labs2 <- do.call(sprintf,list(fmts,c("Y50","2025","2030","2035","2040","2045"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p50 <- ggplot(data=dat50,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(alpha=0.2,col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
         y = "Cum incidence age 50"
         # ,
         # title = "Cumulative incidence at age 50"
         ) +
    scale_linetype_discrete(drop=FALSE)+
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_brewer("",palette="Set1") + scale_fill_brewer("",palette="Set1")
  if(legendonly==FALSE) {
    out <- grid.arrange(p25,p35,p50,nrow=1,top=grid::textGrob(main, x = .05, hjust = 0, gp = grid::gpar(fontsize = 16)))
  } else {
    p25 <- p25 + theme(legend.position="bottom",legend.title=element_text(size=12),
                       legend.text=element_text(size=12))
    tmp <- ggplot_gtable(ggplot_build(p25))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  return(out)
}

cumincid_proj_plot_3 <- function(datm,datf,datm2,datf2,xlim,ylim,main,legendonly=FALSE) {
  proj2m25 <- datm2[datm2[,"proj25",1]==1,,]
  proj2f25 <- datf2[datf2[,"proj25",1]==1,,]
  dat25 <- data.frame(med = c(apply(datm[,"cuminc25",],1,median),apply(proj2m25[,"cuminc25",],1,median),
                              apply(datf[,"cuminc25",],1,median),apply(proj2f25[,"cuminc25",],1,median)),
                      lb = c(apply(datm[,"cuminc25",],1,quantile,.025),apply(proj2m25[,"cuminc25",],1,quantile,.025),
                             apply(datf[,"cuminc25",],1,quantile,.025),apply(proj2f25[,"cuminc25",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc25",],1,quantile,0.975),apply(proj2m25[,"cuminc25",],1,quantile,.975),
                             apply(datf[,"cuminc25",],1,quantile,0.975),apply(proj2f25[,"cuminc25",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m25[,"yob",1],datf[,"yob",1],proj2f25[,"yob",1]),
                      proj = factor(c(datm[,"proj25",1],proj2m25[,"proj25",1]*2,datf[,"proj25",1],
                                      proj2f25[,"proj25",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m25)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f25))))
  proj2m35 <- datm2[datm2[,"proj35",1]==1,,]
  proj2f35 <- datf2[datf2[,"proj35",1]==1,,]
  dat35 <- data.frame(med = c(apply(datm[,"cuminc35",],1,median),apply(proj2m35[,"cuminc35",],1,median),
                              apply(datf[,"cuminc35",],1,median),apply(proj2f35[,"cuminc35",],1,median)),
                      lb = c(apply(datm[,"cuminc35",],1,quantile,.025),apply(proj2m35[,"cuminc35",],1,quantile,.025),
                             apply(datf[,"cuminc35",],1,quantile,.025),apply(proj2f35[,"cuminc35",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc35",],1,quantile,0.975),apply(proj2m35[,"cuminc35",],1,quantile,.975),
                             apply(datf[,"cuminc35",],1,quantile,0.975),apply(proj2f35[,"cuminc35",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m35[,"yob",1],datf[,"yob",1],proj2f35[,"yob",1]),
                      proj = factor(c(datm[,"proj35",1],proj2m35[,"proj35",1]*2,datf[,"proj35",1],
                                      proj2f35[,"proj35",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m35)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f35))))
  proj2m50 <- datm2[datm2[,"proj50",1]==1,,]
  proj2f50 <- datf2[datf2[,"proj50",1]==1,,]
  dat50 <- data.frame(med = c(apply(datm[,"cuminc50",],1,median),apply(proj2m50[,"cuminc50",],1,median),
                              apply(datf[,"cuminc50",],1,median),apply(proj2f50[,"cuminc50",],1,median)),
                      lb = c(apply(datm[,"cuminc50",],1,quantile,.025),apply(proj2m50[,"cuminc50",],1,quantile,.025),
                             apply(datf[,"cuminc50",],1,quantile,.025),apply(proj2f50[,"cuminc50",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc50",],1,quantile,0.975),apply(proj2m50[,"cuminc50",],1,quantile,.975),
                             apply(datf[,"cuminc50",],1,quantile,0.975),apply(proj2f50[,"cuminc50",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m50[,"yob",1],datf[,"yob",1],proj2f50[,"yob",1]),
                      proj = factor(c(datm[,"proj50",1],proj2m50[,"proj50",1]*2,datf[,"proj50",1],
                                      proj2f50[,"proj50",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m50)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f50))))
  fmts <- rep("atop(atop(textstyle(%s), textstyle(bold(%s))),NA)",6) # two rows
  labs2 <- do.call(sprintf,list(fmts,c("Y25","2000","2005","2010","2015","2020"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'")) )
  ex2 <- parse(text=labs2)
  outlab <- test <- parse(text="atop('Y25/35/50 = Year turned 25, 35 or 50',bold('YOB = Year of birth'))")
  p25 <- ggplot(data=dat25,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(y = "Cum incidence age 25"
         # ,
         # title = "Cumulative incidence at age 25"
    )+
    scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))+
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y35","2010","2015","2020","2025","2030"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p35 <- ggplot(data=dat35,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 35"
      # ,title = "Cumulative incidence at age 35"
    )+
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y50","2025","2030","2035","2040","2045"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p50 <- ggplot(data=dat50,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 50"
      # ,
      # title = "Cumulative incidence at age 50"
    ) +
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted"))+
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  if(legendonly==FALSE) {
    out <- arrangeGrob(p25,p35,p50,nrow=1,top=grid::textGrob(main, x = .05, hjust = 0, gp = grid::gpar(fontsize = 16)))
  } else {
    dat25forleg <- dat25
    dat25forleg$sex <- ifelse(dat25forleg$sex=="Female" | dat25forleg$sex=="Female2","Women","Men")
    dat25forleg$sex <- factor(dat25forleg$sex,levels=c("Women","Men"))
    p25forleg <- ggplot(data=dat25forleg,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
      geom_line(aes(lty=proj),size=1.2) +
      geom_ribbon(alpha=0.2,col=NA) +
      theme_bw() +
      scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
      scale_y_continuous(expand=c(0,0),lim=ylim) +
      labs(y = "Cum incidence age 25"
           # ,
           # title = "Cumulative incidence at age 25"
      )+
      scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
      theme(legend.position="none",axis.title.x=element_blank(),
            plot.margin = unit(c(5.5,5.5,0,5.5), "pt"),
            legend.key.size = unit(0.5,"in")) +
      scale_color_brewer(outlab,palette="Set1") +
      scale_fill_brewer(outlab,palette="Set1") +
      guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))
    p25 <- p25forleg + theme(legend.position="bottom",legend.title=element_text(size=12),
                       legend.text=element_text(size=12)) +
      guides(lty=guide_legend(nrow=2,byrow=TRUE))
    tmp <- ggplot_gtable(ggplot_build(p25))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  return(out)
}

cumincid_proj_plot_obsonly <- function(datm,datf,datm2,datf2,xlim=c(1985,2015),ylim,main,legendonly=FALSE) {
  proj2m25 <- datm2[datm2[,"proj25",1]==1,,]
  proj2f25 <- datf2[datf2[,"proj25",1]==1,,]
  dat25 <- data.frame(med = c(apply(datm[,"cuminc25",],1,median),apply(proj2m25[,"cuminc25",],1,median),
                              apply(datf[,"cuminc25",],1,median),apply(proj2f25[,"cuminc25",],1,median)),
                      lb = c(apply(datm[,"cuminc25",],1,quantile,.025),apply(proj2m25[,"cuminc25",],1,quantile,.025),
                             apply(datf[,"cuminc25",],1,quantile,.025),apply(proj2f25[,"cuminc25",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc25",],1,quantile,0.975),apply(proj2m25[,"cuminc25",],1,quantile,.975),
                             apply(datf[,"cuminc25",],1,quantile,0.975),apply(proj2f25[,"cuminc25",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m25[,"yob",1],datf[,"yob",1],proj2f25[,"yob",1]),
                      proj = factor(c(datm[,"proj25",1],proj2m25[,"proj25",1]*2,datf[,"proj25",1],
                                      proj2f25[,"proj25",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m25)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f25))))
  proj2m35 <- datm2[datm2[,"proj35",1]==1,,]
  proj2f35 <- datf2[datf2[,"proj35",1]==1,,]
  dat35 <- data.frame(med = c(apply(datm[,"cuminc35",],1,median),apply(proj2m35[,"cuminc35",],1,median),
                              apply(datf[,"cuminc35",],1,median),apply(proj2f35[,"cuminc35",],1,median)),
                      lb = c(apply(datm[,"cuminc35",],1,quantile,.025),apply(proj2m35[,"cuminc35",],1,quantile,.025),
                             apply(datf[,"cuminc35",],1,quantile,.025),apply(proj2f35[,"cuminc35",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc35",],1,quantile,0.975),apply(proj2m35[,"cuminc35",],1,quantile,.975),
                             apply(datf[,"cuminc35",],1,quantile,0.975),apply(proj2f35[,"cuminc35",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m35[,"yob",1],datf[,"yob",1],proj2f35[,"yob",1]),
                      proj = factor(c(datm[,"proj35",1],proj2m35[,"proj35",1]*2,datf[,"proj35",1],
                                      proj2f35[,"proj35",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m35)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f35))))
  proj2m50 <- datm2[datm2[,"proj50",1]==1,,]
  proj2f50 <- datf2[datf2[,"proj50",1]==1,,]
  dat50 <- data.frame(med = c(apply(datm[,"cuminc50",],1,median),apply(proj2m50[,"cuminc50",],1,median),
                              apply(datf[,"cuminc50",],1,median),apply(proj2f50[,"cuminc50",],1,median)),
                      lb = c(apply(datm[,"cuminc50",],1,quantile,.025),apply(proj2m50[,"cuminc50",],1,quantile,.025),
                             apply(datf[,"cuminc50",],1,quantile,.025),apply(proj2f50[,"cuminc50",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc50",],1,quantile,0.975),apply(proj2m50[,"cuminc50",],1,quantile,.975),
                             apply(datf[,"cuminc50",],1,quantile,0.975),apply(proj2f50[,"cuminc50",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m50[,"yob",1],datf[,"yob",1],proj2f50[,"yob",1]),
                      proj = factor(c(datm[,"proj50",1],proj2m50[,"proj50",1]*2,datf[,"proj50",1],
                                      proj2f50[,"proj50",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m50)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f50))))
  fmts <- rep("atop(atop(textstyle(%s), textstyle(bold(%s))),NA)",6) # two rows
  labs2 <- do.call(sprintf,list(fmts,c("Y25","2000","2005","2010","2015","2020"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'")) )
  ex2 <- parse(text=labs2)
  outlab <- test <- parse(text="atop('Y25/35/50 = Year turned 25, 35 or 50',bold('YOB = Year of birth'))")
  p25 <- ggplot(data=filter(dat25,proj=="Observed"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(y = "Cum incidence age 25"
         # ,
         # title = "Cumulative incidence at age 25"
    )+
    scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))+
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y35","2010","2015","2020","2025","2030"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p35 <- ggplot(data=filter(dat35,proj=="Observed"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 35"
      # ,title = "Cumulative incidence at age 35"
    )+
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y50","2025","2030","2035","2040","2045"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p50 <- ggplot(data=filter(dat50,proj=="Observed"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 50"
      # ,
      # title = "Cumulative incidence at age 50"
    ) +
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted"))+
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  if(legendonly==FALSE) {
    out <- grid.arrange(p25,p35,p50,nrow=1,top=grid::textGrob(main, x = .05, hjust = 0, gp = grid::gpar(fontsize = 16)))
  } else {
    dat25forleg <- dat25
    dat25forleg$sex <- ifelse(dat25forleg$sex=="Female" | dat25forleg$sex=="Female2","Women","Men")
    p25forleg <- ggplot(data=dat25forleg,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
      geom_line(aes(lty=proj),size=1.2) +
      geom_ribbon(alpha=0.2,col=NA) +
      theme_bw() +
      scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
      scale_y_continuous(expand=c(0,0),lim=ylim) +
      labs(y = "Cum incidence age 25"
           # ,
           # title = "Cumulative incidence at age 25"
      )+
      scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
      theme(legend.position="none",axis.title.x=element_blank(),
            plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
      scale_color_brewer(outlab,palette="Set1") +
      scale_fill_brewer(outlab,palette="Set1") +
      guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))
    p25 <- p25forleg + theme(legend.position="bottom",legend.title=element_text(size=12),
                             legend.text=element_text(size=12)) +
      guides(lty=guide_legend(nrow=2,byrow=TRUE))
    tmp <- ggplot_gtable(ggplot_build(p25))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  return(out)
}

cumincid_proj_plot_obscons <- function(datm,datf,datm2,datf2,xlim=c(1985,2015),ylim,main,legendonly=FALSE) {
  proj2m25 <- datm2[datm2[,"proj25",1]==1,,]
  proj2f25 <- datf2[datf2[,"proj25",1]==1,,]
  dat25 <- data.frame(med = c(apply(datm[,"cuminc25",],1,median),apply(proj2m25[,"cuminc25",],1,median),
                              apply(datf[,"cuminc25",],1,median),apply(proj2f25[,"cuminc25",],1,median)),
                      lb = c(apply(datm[,"cuminc25",],1,quantile,.025),apply(proj2m25[,"cuminc25",],1,quantile,.025),
                             apply(datf[,"cuminc25",],1,quantile,.025),apply(proj2f25[,"cuminc25",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc25",],1,quantile,0.975),apply(proj2m25[,"cuminc25",],1,quantile,.975),
                             apply(datf[,"cuminc25",],1,quantile,0.975),apply(proj2f25[,"cuminc25",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m25[,"yob",1],datf[,"yob",1],proj2f25[,"yob",1]),
                      proj = factor(c(datm[,"proj25",1],proj2m25[,"proj25",1]*2,datf[,"proj25",1],
                                      proj2f25[,"proj25",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m25)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f25))))
  proj2m35 <- datm2[datm2[,"proj35",1]==1,,]
  proj2f35 <- datf2[datf2[,"proj35",1]==1,,]
  dat35 <- data.frame(med = c(apply(datm[,"cuminc35",],1,median),apply(proj2m35[,"cuminc35",],1,median),
                              apply(datf[,"cuminc35",],1,median),apply(proj2f35[,"cuminc35",],1,median)),
                      lb = c(apply(datm[,"cuminc35",],1,quantile,.025),apply(proj2m35[,"cuminc35",],1,quantile,.025),
                             apply(datf[,"cuminc35",],1,quantile,.025),apply(proj2f35[,"cuminc35",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc35",],1,quantile,0.975),apply(proj2m35[,"cuminc35",],1,quantile,.975),
                             apply(datf[,"cuminc35",],1,quantile,0.975),apply(proj2f35[,"cuminc35",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m35[,"yob",1],datf[,"yob",1],proj2f35[,"yob",1]),
                      proj = factor(c(datm[,"proj35",1],proj2m35[,"proj35",1]*2,datf[,"proj35",1],
                                      proj2f35[,"proj35",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m35)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f35))))
  proj2m50 <- datm2[datm2[,"proj50",1]==1,,]
  proj2f50 <- datf2[datf2[,"proj50",1]==1,,]
  dat50 <- data.frame(med = c(apply(datm[,"cuminc50",],1,median),apply(proj2m50[,"cuminc50",],1,median),
                              apply(datf[,"cuminc50",],1,median),apply(proj2f50[,"cuminc50",],1,median)),
                      lb = c(apply(datm[,"cuminc50",],1,quantile,.025),apply(proj2m50[,"cuminc50",],1,quantile,.025),
                             apply(datf[,"cuminc50",],1,quantile,.025),apply(proj2f50[,"cuminc50",],1,quantile,.025)),
                      ub = c(apply(datm[,"cuminc50",],1,quantile,0.975),apply(proj2m50[,"cuminc50",],1,quantile,.975),
                             apply(datf[,"cuminc50",],1,quantile,0.975),apply(proj2f50[,"cuminc50",],1,quantile,.975)),
                      yob = c(datm[,"yob",1],proj2m50[,"yob",1],datf[,"yob",1],proj2f50[,"yob",1]),
                      proj = factor(c(datm[,"proj50",1],proj2m50[,"proj50",1]*2,datf[,"proj50",1],
                                      proj2f50[,"proj50",1]*2),levels=c(0,1,2),
                                    labels=c("Observed","Proj, Curr","Proj, Cont Red")),
                      sex = c(rep("Male",nrow(datm)),rep("Male2",nrow(proj2m50)),
                              rep("Female",nrow(datf)),rep("Female2",nrow(proj2f50))))
  fmts <- rep("atop(atop(textstyle(%s), textstyle(bold(%s))),NA)",6) # two rows
  labs2 <- do.call(sprintf,list(fmts,c("Y25","2000","2005","2010","2015","2020"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'")) )
  ex2 <- parse(text=labs2)
  outlab <- test <- parse(text="atop('Y25/35/50 = Year turned 25, 35 or 50',bold('YOB = Year of birth'))")
  p25 <- ggplot(data=filter(dat25,proj=="Observed" | proj=="Proj, Curr"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(y = "Cum incidence age 25"
         # ,
         # title = "Cumulative incidence at age 25"
    )+
    scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))+
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y35","2010","2015","2020","2025","2030"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p35 <- ggplot(data=filter(dat35,proj=="Observed" | proj=="Proj, Curr"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 35"
      # ,title = "Cumulative incidence at age 35"
    )+
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted")) +
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  labs2 <- do.call(sprintf,list(fmts,c("Y50","2025","2030","2035","2040","2045"),
                                c("YOB","'1975'","'1980'","'1985'","'1990'","'1995'") ) )
  ex2 <- parse(text=labs2)
  p50 <- ggplot(data=filter(dat50,proj=="Observed" | proj=="Proj, Curr"),aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
    geom_line(aes(lty=proj),size=1.2) +
    geom_ribbon(aes(alpha=sex),col=NA) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2,lim=xlim) +
    scale_y_continuous(expand=c(0,0),lim=ylim) +
    labs(
      y = "Cum incidence age 50"
      # ,
      # title = "Cumulative incidence at age 50"
    ) +
    scale_linetype_manual(drop=FALSE,values = c("solid","dashed","dotted"))+
    theme(legend.position="none",axis.title.x=element_blank(),
          plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
    scale_color_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_fill_manual(outlab,values = c("Female"="#E41A1C","Female2"="#E41A1C","Male"="#377EB8","Male2"="#377EB8"),
                      breaks=c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men")) +
    scale_alpha_manual(values = c("Female"=0.1,"Female2"=0.1,"Male"=0.1,"Male2"=0.1),
                       breaks = c("Female","Female2","Male","Male2"),labels=c("Women","Women","Men","Men"))
  if(legendonly==FALSE) {
    out <- grid.arrange(p25,p35,p50,nrow=1,top=grid::textGrob(main, x = .05, hjust = 0, gp = grid::gpar(fontsize = 16)))
  } else {
    dat25forleg <- dat25
    dat25forleg$sex <- ifelse(dat25forleg$sex=="Female" | dat25forleg$sex=="Female2","Women","Men")
    p25forleg <- ggplot(data=dat25forleg,aes(x=yob+15,y=med,ymin=lb,ymax=ub,col=sex,fill=sex)) +
      geom_line(aes(lty=proj),size=1.2) +
      geom_ribbon(alpha=0.2,col=NA) +
      theme_bw() +
      scale_x_continuous(expand=c(0,0),breaks=c(1985,1990,1995,2000,2005,2010),labels=ex2) +
      scale_y_continuous(expand=c(0,0),lim=ylim) +
      labs(y = "Cum incidence age 25"
           # ,
           # title = "Cumulative incidence at age 25"
      )+
      scale_linetype_manual("",drop=FALSE,values = c("solid","dashed","dotted")) +
      theme(legend.position="none",axis.title.x=element_blank(),
            plot.margin = unit(c(5.5,5.5,0,5.5), "pt")) +
      scale_color_brewer(outlab,palette="Set1") +
      scale_fill_brewer(outlab,palette="Set1") +
      guides(lty = guide_legend(order = 2),col = guide_legend(order = 1), fill=guide_legend(order=1))
    p25 <- p25forleg + theme(legend.position="bottom",legend.title=element_text(size=12),
                             legend.text=element_text(size=12)) +
      guides(lty=guide_legend(nrow=2,byrow=TRUE))
    tmp <- ggplot_gtable(ggplot_build(p25))
    leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
    out <- cowplot::ggdraw() +
      cowplot::draw_grob(grid::grobTree(tmp$grobs[[leg]]))
  }
  return(out)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cumincid_coh_plot_oneage <- function(cuminc,stand,main,years = as.character(seq(1970,1995,by=5)),
                                     ylim = c(0,.65),legendonly=FALSE,xlim,ncolors=6) {
  data <- data.frame(y = numeric(), ymin = numeric(), ymax = numeric(), age = numeric(), coh = character())
  for(i in 1:length(years)){
    curyear <- years[i]
    subdat <- array(unlist(lapply(cuminc,"[",curyear)),
                    dim=c(length(unlist(cuminc[[1]][curyear])),length(cuminc))) #each row is new age
    # if(main=="Karonga men") {
    #   subdat <- subdat[-c(nrow(subdat),nrow(subdat)-1,nrow(subdat)-2,
    #                       nrow(subdat)-3,nrow(subdat)-4),]
    # }
    y <- apply(subdat,1,median)
    ymin <- apply(subdat,1,quantile,0.025)
    ymax <- apply(subdat,1,quantile,0.975)
    age <- seq(15.1,length(y)*stand$dt+15.1 - stand$dt,by=stand$dt)
    year <- rep(curyear,length(y))
    data <- rbind(data,data.frame(y=y,ymin=ymin,ymax=ymax,age=age,coh=year))
    data <- subset(data,age==24.9)
  }
  if(main=="Manicaland men" | main=="Manicaland women" | main=="Karonga men" | main=="Karonga women") {
    data <- rbind(data,data.frame(y=NA,ymin=NA,ymax=NA,age=24.9,coh="1990"))
  }
  colors <- gg_color_hue(ncolors)
  colscheme <- c(colors[1:(ncolors-1)])
  names(colscheme) <- levels(droplevels(data$coh))
  out <- ggplot(data=data,aes(y=y,x=coh)) +
    geom_col(aes(col=coh,fill=coh),position = position_dodge(width=1)) +
    geom_errorbar(aes(ymin=ymin,ymax=ymax),position = position_dodge(width=1),width=.3,
                  col="black")+
    theme_bw()+
    xlab("Cohort Year of Birth") +
    ylab("Cumulative Incidence (age 25)")+
    ggtitle(main) +
    theme(legend.position="none") +
    scale_color_manual(values=colscheme) +
    scale_fill_manual(values=colscheme)
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  out
}

cumincid_coh_plot_50 <- function(cuminc,stand,main,years = as.character(seq(1945,1965,by=5)),
                                     ylim = c(0,1),legendonly=FALSE,ncolors=6) {
  data <- data.frame(y = numeric(), ymin = numeric(), ymax = numeric(), age = numeric(), coh = character())
  for(i in 1:length(years)){
    curyear <- years[i]
    subdat <- array(unlist(lapply(cuminc,"[",curyear)),
                    dim=c(length(unlist(cuminc[[1]][curyear])),length(cuminc))) #each row is new age
    # if(main=="Karonga men") {
    #   subdat <- subdat[-c(nrow(subdat),nrow(subdat)-1,nrow(subdat)-2,
    #                       nrow(subdat)-3,nrow(subdat)-4),]
    # }
    y <- apply(subdat,1,median)
    ymin <- apply(subdat,1,quantile,0.025)
    ymax <- apply(subdat,1,quantile,0.975)
    if((stand$x_time[1] - 15)>=as.numeric(curyear)) {
      bottom <- stand$x_time[1] - as.numeric(curyear) + 0.1
    } else { bottom <- 15.1}
    age <- seq(bottom,length(y)*stand$dt+bottom - stand$dt,by=stand$dt)
    year <- rep(curyear,length(y))
    data <- rbind(data,data.frame(y=y,ymin=ymin,ymax=ymax,age=age,coh=year))
    data <- subset(data,round(age,1)==49.9)
  }
  if(main=="Manicaland men" | main=="Manicaland women" | main=="Karonga men" | main=="Karonga women") {
    data <- rbind(data,data.frame(y=NA,ymin=NA,ymax=NA,age=24.9,coh="1965"))
  }
  colors <- gg_color_hue(ncolors)
  colscheme <- c(colors[1:(ncolors-1)])
  names(colscheme) <- levels(droplevels(data$coh))
  out <- ggplot(data=data,aes(y=y,x=coh)) +
    geom_col(aes(col=coh,fill=coh),position = position_dodge(width=1)) +
    geom_errorbar(aes(ymin=ymin,ymax=ymax),position = position_dodge(width=1),width=.3,
                  col="black")+
    theme_bw()+
    xlab("Cohort Year of Birth") +
    ylab("Cumulative Incidence (age 50)")+
    ggtitle(main) +
    theme(legend.position="none") +
    scale_color_manual(values=colscheme) +
    scale_fill_manual(values=colscheme)
  if(is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0))
  } else if(!is.na(ylim[1])) {
    out <- out + scale_y_continuous(expand=c(0,0),limits=ylim)
  }
  out
}


prev_post_plot <- function(dat,titlename) {
  ggplot(data = filter(dat,!is.na(surveyround)),aes(x=agecat)) +
    geom_crossbar(aes(ymin=prevlow80,ymax=prevhigh80,y=prevmean),fill="seagreen",colour=NA,alpha=.7)+
    geom_crossbar(aes(ymin=prevlow95,ymax=prevhigh95,y=prevmean),fill="seagreen",colour=NA,alpha=.4)+
    geom_crossbar(aes(y=prevmean,ymin=prevmean,ymax=prevmean),colour="seagreen",fatten=1,width=0.78) +
    geom_point(aes(y=prevalence)) +
    facet_wrap(~(.5*round(year/.5)),ncol=3) +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    scale_y_continuous(expand=c(0,0),lim=c(0,NA)) +
    scale_x_discrete(labels=c("15-19","20-24","25-29","30-34","35-39",
                              "40-44","45-49","50-54")) +
    labs(x="Age",
         y="HIV Prevalence",
         title=titlename)
}

