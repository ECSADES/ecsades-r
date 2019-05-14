#' Drawing a sample from a fitted joint distribution
#'
#' @description 
#' This function draws a random sample from an existing fitted Weibull-log-normal
#' (\code{wln}) or Heffernan-Tawn (\code{ht}) model for a given duration.
#' 
#' @param jdistr a \code{wln} or \code{ht} joint distribution object, as the output from function \code{\link{fit_wln}}
#' and \code{\link{fit_ht}} respectively.
#' 
#' @param sim_year the period of the simulate data in equivalent number of years.
#' 
#' @param perturbed_ht_residuals whether or not to the simulation from the Heffernan-Tawn (\code{ht}) model should
#' add a perturbation when sampling from the empirical residual distribution (default value is \code{TRUE}).  This
#' argument is ignored when the input object is a Weibull-log-normal (\code{wln}) distribution.
#'
#' @return The function returns \code{data.table} object with two columns - wave
#' height \code{hs} and wave period and \code{tp}.
#'
#' @examples
#' 
#' # Load data
#' data(noaa_ww3)
#' 
#' # Fit Heffernan-Tawn model with 0.95 marginal and dependence thresholds
#' noaa_ht = fit_ht(data = noaa_ww3, npy = nrow(noaa_ww3)/10, margin_thresh_count = 100, dep_thresh_count = 100)
#' 
#' # Draw the equivalent of 100-year of data from the fitted Heffernan-Tawn model with or without perturbed residuals
#' sim_ht1 = sample_jdistr(jdistr = noaa_ht, sim_year = 100, perturbed_ht_residuals = FALSE)
#' sim_ht2 = sample_jdistr(jdistr = noaa_ht, sim_year = 100)
#' par(mfrow=c(1,2))
#' plot(sim_ht1, pch=20, cex=.5, main="HT without perturbation")
#' plot(sim_ht2, pch=20, cex=.5, main="HT with perturbation")
#' 
#' # Fit Weibull-log-normal model
#' noaa_wln = fit_wln(data = noaa_ww3, npy = nrow(noaa_ww3)/10)
#'
#' # Draw the equivalent of 100-year of data from the fitted Weibull-log-normal model
#' sim_wln = sample_jdistr(jdistr = noaa_wln, sim_year = 100)
#' plot(sim_wln, pch=20, cex=.5, main="WLN")
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{fit_wln}}
#' 
#' @export
sample_jdistr = function(jdistr, sim_year, perturbed_ht_residuals = TRUE){
  
  if(class(jdistr)=="ht"){
    res = .sample_ht(jdistr, sim_year, perturbed_ht_residuals)
  }else if(class(jdistr)=="wln"){
    res = .sample_wln(jdistr, sim_year)
  }else{
    stop("Input class of distribution not supported.")
  }
  return(res)
}


# Heffernan-Tawn ----------------------------------------------------------
.sample_ht = function(ht, sim_year, perturbed_ht_residuals){
  
  set.seed(.seed_sampling)
  n_sim = sim_year*ht$npy

  # Sample overall freq
  res = ht$dep$dep_data[sample.int(.N, size = n_sim, replace = T), .(u_hs, u_tp)]
  rot_wb = res[, .N^(-1/6)]
  res[, u_hs:=pnorm(qnorm(u_hs)+rnorm(.N, 0, rot_wb))]
  res[, u_tp:=pnorm(qnorm(u_tp)+rnorm(.N, 0, rot_wb))]
  n_hs = res[u_hs>ht$dep$p_dep_thresh & u_hs>u_tp, .N]
  n_tp = res[u_tp>ht$dep$p_dep_thresh & u_tp>=u_hs, .N]
  
  # Sample hs & tp extreme in unif
  dep_data_hs = .sample_ht1(n_hs, ht$dep$hs$par, ht$dep$hs$resid, ht$dep$p_dep_thresh, perturbed_ht_residuals)
  res[u_hs>ht$dep$p_dep_thresh & u_hs>u_tp, c("u_hs", "u_tp"):=dep_data_hs[, .(u_cond, u_dep)]]
  dep_data_tp = .sample_ht1(n_tp, ht$dep$tp$par, ht$dep$tp$resid, ht$dep$p_dep_thresh, perturbed_ht_residuals)
  res[u_tp>ht$dep$p_dep_thresh & u_tp>=u_hs, c("u_hs", "u_tp"):=dep_data_tp[, .(u_dep, u_cond)]]

  # Convert to original scale
  n = length(ht$margin$hs$emp)
  res[, hs:=quantile(ht$margin$hs$emp, u_hs)]
  res[u_hs>ht$margin$p_margin_thresh, u_gpd:=(u_hs-ht$margin$p_margin_thresh)/(1-ht$margin$p_margin_thresh)]
  res[u_hs>ht$margin$p_margin_thresh, hs:=evd::qgpd(
    p = u_gpd, loc = ht$margin$hs$par[1],
    scale = ht$margin$hs$par[2], shape = ht$margin$hs$par[3])]
  res[, u_gpd:=NULL]
  res[, tp:=quantile(ht$margin$tp$emp, u_tp)]
  res[u_tp>ht$margin$p_margin_thresh, u_gpd:=(u_tp-ht$margin$p_margin_thresh)/(1-ht$margin$p_margin_thresh)]
  res[u_tp>ht$margin$p_margin_thresh, tp:=evd::qgpd(
    p = u_gpd, loc = ht$margin$tp$par[1],
    scale = ht$margin$tp$par[2], shape = ht$margin$tp$par[3])]
  
  return(res[, .(hs, tp)])
}

.sample_ht1 = function(n_sim, par, resid, p_dep_thresh, perturbed_ht_residuals){
  set.seed(.seed_sampling)
  a = par[1]
  b = par[2]
  sim_data = data.table(u_cond=runif(n_sim, p_dep_thresh, 1))
  sim_data[u_cond<0.5, l_cond:=log(2*u_cond)]
  sim_data[u_cond>=0.5, l_cond:=-log(2-2*u_cond)]
  resid_max = sim_data[, (1-a)*l_cond^(1-b)]
  if(perturbed_ht_residuals){
    resid_bw = bw.SJ(resid)
    extended_resid = resid+rnorm(length(resid)*100, 0, resid_bw)
  }else{
    extended_resid = resid
  }
  resid_sample = sapply(resid_max, function(x)sample(extended_resid[extended_resid<x],1))
  sim_data[, l_dep:=resid_sample*l_cond^b+a*l_cond]
  sim_data[l_dep<0, u_dep:=.5*exp(l_dep)]
  sim_data[l_dep>=0, u_dep:=1-.5*exp(-l_dep)]
  return(sim_data[, .(u_cond, u_dep)])
}



# Weibull log-normal ------------------------------------------------------
.sample_wln = function(wln, sim_year){
  set.seed(.seed_sampling)
  n_sim = round(wln$npy*sim_year)
  
  # Sample hs
  hs = rweibull(n_sim, shape = wln$hs$par["shape"], scale = wln$hs$par["scale"]) + wln$hs$par["loc"]
  
  # sample tp cond on hs
  tp_norm_mean = wln$tp$par[1] + wln$tp$par[2] * (hs ^ wln$tp$par[3])
  tp_norm_sd = wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs)
  tp = exp(rnorm(n_sim, tp_norm_mean, tp_norm_sd))
  
  res = data.table(hs, tp)
  return(res)
}

