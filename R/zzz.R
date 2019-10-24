.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
.onAttach <- function(libname, pkgname){
##  rstan::expose_stan_functions(stanmodels$functions, envir = as.environment("package:alphaepi"))
}
