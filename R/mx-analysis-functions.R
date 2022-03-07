calc.e15 <- function(sx, dt, stand){
  aidx <- which(colnames(sx)==as.character(15 + dt)):(which(colnames(sx)==max(stand$x_age)))
  dt*colSums(apply(sx[,aidx], 1, cumprod))
}

calc.45q15 <- function(sx, dt, stand){
  # aidx <- seq_len(45/dt) # assumes psurv starts at age 15
  aidx <- which(colnames(sx)==as.character(15 + dt)):(which(colnames(sx)==as.character(60)))
  1.0 - exp(rowSums(log(sx[,aidx])))
}


calc.sx <- function(param, stand, base_tIDX=which(stand$x_time == 2005),ind){

  psurvobj <- calc.psurv(param, stand, ind)
  psurv <- psurvobj$psurv
  psurv.nohiv <- psurvobj$psurv.nohiv
  psurv.noart <- psurvobj$psurv.noart

  psurv.last <- psurv[-nrow(psurv), -ncol(psurv)]
  psurv.curr <- psurv[-1,-1]
  sx <- 1.0 - (psurv.last-psurv.curr)/psurv.last
  
  psurv.nohiv.last <- psurv.nohiv[-nrow(psurv.nohiv), -ncol(psurv.nohiv)]
  psurv.nohiv.curr <- psurv.nohiv[-1,-1]
  sx.nohiv <- 1.0 - (psurv.nohiv.last-psurv.nohiv.curr)/psurv.nohiv.last
  
  psurv.noart.last <- psurv.noart[-nrow(psurv.noart), -ncol(psurv.noart)]
  psurv.noart.curr <- psurv.noart[-1,-1]
  sx.noart <- 1.0 - (psurv.noart.last-psurv.noart.curr)/psurv.noart.last

  dimnames(sx) <- dimnames(sx.nohiv) <-   dimnames(sx.noart) <- list(stand$x_time[-1], stand$x_age[-1])
  
  sx.hiv <- sx / sx.nohiv              # probability of dying from HIV
  sx.noart.hiv <- sx.noart / sx.nohiv  # probability fo dying from HIV, no ART counterfactual

  ## Calculate sx for HIV positive persons
  sx.hivp <- sx.nohiv * (1 - (1-sx.hiv) / psurvobj$prev[-1,-1])
  sx.hivp.noart <- sx.nohiv * (1 - (1-sx.noart.hiv) / psurvobj$prev.noart[-1,-1])
  
  ## Calculate no-background mortality counterfactual as probability of dying from HIV
  ## times the non-HIV death probability in the counterfactual year

  ## Constant counterfactual mortality from base_tIDX (default 2005)
  omitIDX <- -(1:(base_tIDX-2))
  sx.natconst05 <- sweep(sx.hiv[omitIDX,], 2, sx.nohiv[base_tIDX-1,], "*")              # total prob of dying, constant natural mortality
  sx.noart.natconst05 <- sweep(sx.noart.hiv[omitIDX,], 2, sx.nohiv[base_tIDX-1,], "*")  # total prob of dying
  sx.hivconst05 <- sweep(sx.nohiv[omitIDX,], 2, sx.hiv[base_tIDX-1,], "*")

  ## Constant counterfactual mortality from the time of ART start
  artstart_tIDX <- stand$artstart_tIDX
  omitIDX <- -(1:(artstart_tIDX-2))
  sx.natconstART <- sweep(sx.hiv[omitIDX,], 2, sx.nohiv[artstart_tIDX-1,], "*")              # total prob of dying, constant natural mortality
  sx.noart.natconstART <- sweep(sx.noart.hiv[omitIDX,], 2, sx.nohiv[artstart_tIDX-1,], "*")  # total prob of dying
  sx.hivconstART <- sweep(sx.nohiv[omitIDX,], 2, sx.hiv[artstart_tIDX-1,], "*")

  ## Restrict intervals for no-ART counterfactual
  sx.noart <- sx.noart[-(1:(artstart_tIDX-2)),]
  if(length(stand$x_natmx) < length(stand$x_time))
    sx.nohiv <- sx.nohiv[-(1:(length(stand$x_time) - length(stand$x_natmx) - 1)),]
                       
  sx.list <- list(obs        = sx,
                  noart      = sx.noart,
                  nohiv      = sx.nohiv,
                  hivpos     = sx.hivp,
                  hivpos.noart = sx.hivp.noart,
                  natconst05 = sx.natconst05,
                  hivconst05 = sx.hivconst05,
                  noart.natconst05 = sx.noart.natconst05,
                  natconstART = sx.natconstART,
                  hivconstART = sx.hivconstART,
                  noart.natconstART = sx.noart.natconstART)
  return(sx.list)
}

calc.mx.output <- function(param, stand, base_tIDX, sx.times, ind){
  
  sxlist <- calc.sx(param, stand, base_tIDX, ind) ### Does this need to be updated b/c starting model @ 10?  I think no...

  e15 <- lapply(sxlist, calc.e15, stand$dt, stand)
  q4515 <- lapply(sxlist, calc.45q15, stand$dt, stand)
  names(e15) <- paste("e15", names(e15), sep=".")
  names(q4515) <- paste("q4515", names(q4515), sep=".")

  extract.sx <- function(sx, times, dt){
    sx_tIDX <- as.integer(round((times - as.numeric(rownames(sx)[1])) / dt)) + 1L
    sx_tIDX <- sx_tIDX[sx_tIDX >=1 & sx_tIDX <= nrow(sx)]
    sx[sx_tIDX,]
  }

  sxlist <- lapply(sxlist[c("obs", "noart", "nohiv")], extract.sx, sx.times, stand$dt)
  names(sxlist) <- paste("sx", names(sxlist), sep=".")

  return(c(e15, q4515, sxlist))
}
  
add.mx <- function(mod, base.time = 2005, sx.times=seq(1990, 2015, 5), ind, tot=2){
  
  print(paste("starting add.mx()", Sys.time()))
  
  mod$stand$artstart_tIDX <- mod$stand$artstart_tIDX[ind]

  base_tIDX <- which.min(abs(mod$stand$x_time - base.time))
  
  param <- create.param.list(mod$fit,ind,tot)
  system.time(output <- lapply(param, calc.mx.output, mod$stand, base_tIDX, sx.times, ind)) # mclapply causing problems so lapply for now
  output <- do.call(mapply, c(list(FUN=cbind, SIMPLIFY=FALSE), output))  # reorder by outcome

  return(output)
}


calc.cumincid <- function(param, stand, ages=min(stand$xage), nyears=diff(range(stand$x_age))){

  log_incrate_time_age <- stand$X_incrate_time %*% param$coef_incrate_time_age %*% t(stand$Xmid_incrate_age)
  min.age <- min(stand$x_age)
  
  extract.cumincid <- function(a, n){
    aidx <- (a - min.age) / stand$dt + seq_len(n / stand$dt)
    1.0 - exp(-stand$dt * rowSums(exp(log_incrate_time_age[,aidx])))
  }
  cumincid.period <- mapply(extract.cumincid, ages, nyears)
  dimnames(cumincid.period) <- list(stand$x_time, paste("i", nyears, "_", ages, sep=""))
  return(cumincid.period)
}


# add.incprev <- function(mod, incprev.times=seq(1970, 2015, 5),
#                         cumincid.ages=c(15, 15, 15, 25, 35, 45, 55, 15, 20),
#                         cumincid.nyears=c(85, 45, 10, 10, 10, 10, 10, 5, 5)){
#   stand <- mod$stand
#   param <- create.param.list(mod$fit)
# 
# 
#   ## prev
#   system.time(prev <- mclapply(param, calc.prev, mod$stand))
#   system.time(prev <- array(unlist(prev), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param))))
#   dimnames(prev) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
# 
#   ## incid
#   system.time(incid <- lapply(param, calc.incid, mod$stand))
#   incid <- array(unlist(incid), c(mod$stand$STEPS_time, mod$stand$STEPS_age, length(param)))
#   dimnames(incid) <- list(mod$stand$x_time, mod$stand$x_age, NULL)
# 
#   ## cumincid
#   cumincid <- array(sapply(param, calc.cumincid, mod$stand, cumincid.ages, cumincid.nyears),
#                     c(mod$stand$STEPS_time, length(cumincid.ages), length(param)))
#   dimnames(cumincid) <- list(stand$x_time, paste("i", cumincid.nyears, "_", cumincid.ages, sep=""), NULL)
# 
#   ## WHO prev
#   aidx <- 1:(35/mod$stand$dt+1L)
#   age.dist <- approx(0:99+0.5, who.standard.pop, mod$stand$x_age[aidx])$y  # prevalence age 15 to 50
#   age.dist <- age.dist/sum(age.dist)
#   who15to49prev <- apply(sweep(prev[, aidx,], 2, age.dist, "*"), c(1,3), sum)
# 
#   incprev_tIDX <- as.integer(round((incprev.times - stand$x_time[1]) / stand$dt)) + 1L
#   incprev_tIDX <- incprev_tIDX[incprev_tIDX >=1 & incprev_tIDX <= length(stand$x_time)]
# 
#   list(prev=prev[incprev_tIDX,,],
#        incid=incid[incprev_tIDX,,],
#        cumincid=cumincid,
#        who15to49prev=who15to49prev)
# }



#############################
####  Non-HIV mortality  ####
#############################

calc.nat45q15 <- function(param, stand){
  aidx <- seq_len(45/stand$dt)
  log_natmx_time_age <- stand$X_natmx_time[stand$x_time %in% stand$x_natmx,] %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age[aidx,])
  nat45q15 <- 1.0 - exp(-stand$dt * rowSums(exp(log_natmx_time_age)))
  setNames(nat45q15, stand$x_natmx)
}

calc.nat.e15 <- function(param, stand){
  lognatmx <- stand$X_natmx_time[stand$x_time %in% stand$x_natmx,] %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age)
  nate15 <- stand$dt*colSums(exp(-stand$dt*apply(exp(lognatmx), 1, cumsum)))
  setNames(nate15, stand$x_natmx)
}


calc.lognatmx2010 <- function(param, stand, years=70){
  aidx <- seq_len(years/stand$dt+1)
  lognatmx <- stand$X_natmx_time[stand$x_time == 2010,] %*% param$coef_natmx_time_age %*% t(stand$X_natmx_age[aidx,])
  setNames(lognatmx, stand$x_age[aidx])
}

calc.natsx2010 <- function(param, stand, years=70){
  aidx <- seq_len(years/stand$dt)
  lognatmx <- stand$X_natmx_time[stand$x_time == 2010,] %*% param$coef_natmx_time_age %*% t(stand$Xmid_natmx_age[aidx,])
  sx <- exp(-stand$dt*tapply(exp(lognatmx), floor(stand$x_age[aidx]), sum))
  return(sx)
}

add.nonhivmx <- function(mod){
  param <- create.param.list(mod$fit)
  mod$nat45q15 <- sapply(param, calc.nat45q15, mod$stand)
  mod$nat.e15 <- sapply(param, calc.nat.e15, mod$stand)
  mod$lognatmx2010 <- sapply(param, calc.lognatmx2010, mod$stand)
  mod$natsx2010 <- sapply(param, calc.natsx2010, mod$stand)
  return(mod)
}




sum.sx <- function(mod){
  x <- mod$natsx2010
  data.frame(mean=rowMeans(x),
             cil=apply(x, 1, quantile, 0.025),
             ciu=apply(x, 1, quantile, 0.975))
}

estci <- function(x)
  data.frame(mean=rowMeans(x),
             cil=apply(x, 1, quantile, 0.025),
             ciu=apply(x, 1, quantile, 0.975))

