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
#' data(ww3_pk)
#' 
#' # Fit Heffernan-Tawn model with 0.95 marginal and dependence thresholds
#' ht = fit_ht(data = ww3_pk, npy = nrow(ww3_pk)/10, margin_thresh_count = 100, dep_thresh_count = 100)
#' 
#' # Draw the equivalent of 100-year of data from the fitted Heffernan-Tawn model with or without perturbed residuals
#' sim_ht1 = sample_jdistr(jdistr = ht, sim_year = 100, perturbed_ht_residuals = FALSE)
#' sim_ht2 = sample_jdistr(jdistr = ht, sim_year = 100)
#' par(mfrow=c(1,2))
#' plot(sim_ht1, pch=20, cex=.5, main="HT without perturbation")
#' plot(sim_ht2, pch=20, cex=.5, main="HT with perturbation")
#' 
#' # Fit Weibull-log-normal model
#' wln = fit_wln(data = ww3_pk, npy = nrow(ww3_pk)/10)
#'
#' # Draw the equivalent of 100-year of data from the fitted Weibull-log-normal model
#' sim_wln = sample_jdistr(jdistr = wln, sim_year = 100)
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
  res[, hs:=.convert_unif_to_origin(
    unif = u_hs, p_thresh = ht$margin$p_margin_thresh, gpd_par = ht$margin$hs$par, emp = ht$margin$hs$emp)]
  res[, tp:=.convert_unif_to_origin(
    unif = u_tp, p_thresh = ht$margin$p_margin_thresh, gpd_par = ht$margin$tp$par, emp = ht$margin$tp$emp)]
  
  return(res[, .(hs, tp)])
}

.sample_ht1 = function(n_sim, par, resid, p_dep_thresh, perturbed_ht_residuals){
  set.seed(.seed_sampling)
  a = par[1]
  b = par[2]
  sim_data = data.table(u_cond=runif(n_sim, p_dep_thresh, 1))
  sim_data[, l_cond:=.convert_unif_to_lap(u_cond)]
  resid_max = sim_data[, (1-a)*l_cond^(1-b)]
  if(perturbed_ht_residuals){
    resid_bw = bw.SJ(resid)
    extended_resid = resid+rnorm(length(resid)*100, 0, resid_bw)
  }else{
    extended_resid = resid
  }
  resid_sample = sapply(resid_max, function(x)sample(extended_resid[extended_resid<x],1))
  sim_data[, l_dep:=resid_sample*l_cond^b+a*l_cond]
  sim_data[, u_dep:=.convert_lap_to_unif(l_dep)]
  return(sim_data[, .(u_cond, u_dep)])
}

.sample_ht_is = function(ht, target_rp){
  
  # Common
  n_sim = round(ht$npy*target_rp*.target_rp_ub)
  tail_prob = 1/(target_rp*.target_rp_lb*ht$npy)
  min_r = qnorm(1-tail_prob)
  n_tail = n_sim*exp(-min_r^2/2)

  # Generate importance samples outside the Rosenblatt transformed contour of tail_rtrp
  set.seed(.seed_sampling)
  calc = data.table(ur2 = runif(n_tail, min=pchisq(q = min_r^2,df = 2)))
  calc[, r:=sqrt(qchisq(p=ur2, df = 2))]
  calc[, dir:=runif(.N, min=0, max=2*pi)]
  
  # Generate by Tp|Hs model
  calc_hs = calc[, .(u_hs=pnorm(cos(dir)*r), u_tp_hs=pnorm(sin(dir)*r))]
  calc_hs[, lap_hs:=.convert_unif_to_lap(u_hs)]
  calc_hs[u_hs>=ht$dep$p_dep_thresh,lap_tp:=
            lap_hs*ht$dep$hs$par[["a"]]+lap_hs^ht$dep$hs$par[["b"]]*quantile(ht$dep$hs$resid, u_tp_hs)]
  calc_hs[u_hs>=ht$dep$p_dep_thresh, u_tp:=.convert_lap_to_unif(lap_tp)]
  calc_hs = calc_hs[u_hs>ht$dep$p_dep_thresh][u_hs>u_tp]
  calc_hs[, dir:=atan2(qnorm(u_hs), qnorm(u_tp))]
  
  # Generate by Tp|Hs model
  calc_tp = calc[, .(u_hs_tp=pnorm(cos(dir)*r), u_tp=pnorm(sin(dir)*r))]
  calc_tp[, lap_tp:=.convert_unif_to_lap(u_tp)]
  calc_tp[u_tp>=ht$dep$p_dep_thresh,lap_hs:=
            lap_tp*ht$dep$tp$par[["a"]]+lap_tp^ht$dep$tp$par[["b"]]*quantile(ht$dep$tp$resid, u_hs_tp)]
  calc_tp[u_tp>=ht$dep$p_dep_thresh, u_hs:=.convert_lap_to_unif(lap_hs)]
  calc_tp = calc_tp[u_tp>ht$dep$p_dep_thresh][u_tp>u_hs]
  calc_tp[, dir:=atan2(qnorm(u_hs), qnorm(u_tp))]
  
  # Generate low/low corner
  n_ll = calc[, .N]-calc_hs[, .N]-calc_tp[, .N]
  set.seed(.seed_sampling)
  res = ht$dep$dep_data[sample.int(.N, size = calc[, .N]/.target_rp_lb, replace = T), .(u_hs, u_tp)]
  # print(res[, .N])
  rot_wb = res[, .N^(-1/6)]
  res[, u_hs:=pnorm(qnorm(u_hs)+rnorm(.N, 0, rot_wb))]
  res[, u_tp:=pnorm(qnorm(u_tp)+rnorm(.N, 0, rot_wb))]
  res[, r:=sqrt(qnorm(u_hs)^2+qnorm(u_tp)^2)]
  res[, dir:=atan2(qnorm(u_hs), qnorm(u_tp))]
  res_ll = res[dir>max(calc_hs$dir) | dir<min(calc_tp$dir)]
  res_ll[, dir_grp:=round(dir, 3)]
  calc_ll = res_ll[, .SD[r>quantile(r, 1-n_ll/res_ll[,.N])], .(dir_grp)]
  
  # Merge
  out = rbind(calc_hs[, .(u_hs, u_tp)], calc_tp[, .(u_hs, u_tp)], calc_ll[, .(u_hs, u_tp)])
  out[, hs:=.convert_unif_to_origin(
    unif = u_hs, p_thresh = ht$margin$p_margin_thresh,
    gpd_par = ht$margin$hs$par, emp = ht$margin$hs$emp)]
  out[, tp:=.convert_unif_to_origin(
    unif = u_tp, p_thresh = ht$margin$p_margin_thresh,
    gpd_par = ht$margin$tp$par, emp = ht$margin$tp$emp)]
  
  return(out[, .(hs, tp)])
}


# Weibull log-normal ------------------------------------------------------
.sample_wln = function(wln, sim_year){
  set.seed(.seed_sampling)
  n_sim = round(wln$npy*sim_year)
  
  # Sample hs
  hs = rweibull(n_sim, shape = wln$hs$par["shape"], scale = wln$hs$par["scale"]) + wln$hs$par["loc"]
  
  # sample tp cond on hs
  tp_norm_mean = wln$tp$par[1] + wln$tp$par[2] * (hs ^ wln$tp$par[3])
  # tp_norm_sd = wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs)
  tp_norm_sd = pmax(.limit_zero, wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs))
  tp = exp(rnorm(n_sim, tp_norm_mean, tp_norm_sd))
  
  res = data.table(hs, tp)
  return(res)
}

.sample_wln_is = function(wln, target_rp){
  
  n_sim = round(wln$npy*target_rp*.target_rp_ub)
  tail_prob = 1/(target_rp*.target_rp_lb*wln$npy)
  min_r = qnorm(1-tail_prob)
  n_tail = n_sim*exp(-min_r^2/2)
  
  # Generate importance samples outside the Rosenblatt transformed contour of tail_rtrp
  set.seed(.seed_sampling)
  calc = data.table(ur2 = runif(n_tail, min=pchisq(q = min_r^2,df = 2)))
  calc[, r:=sqrt(qchisq(p=ur2, df = 2))]
  calc[, dir:=runif(.N, min=0, max=2*pi)]
  calc[, u_hs:=pnorm(cos(dir)*r)]
  calc[, u_tp:=pnorm(sin(dir)*r)]
  
  # Sample hs
  calc[, hs:= qweibull(u_hs, shape = wln$hs$par["shape"], scale = wln$hs$par["scale"]) + wln$hs$par["loc"]]
  
  # sample tp cond on hs
  calc[, tp_norm_mean:= wln$tp$par[1] + wln$tp$par[2] * (hs ^ wln$tp$par[3])]
  # calc[, tp_norm_sd:= wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs)]
  calc[, tp_norm_sd:= pmax(.limit_zero, wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs))]
  calc[, tp:= exp(qnorm(u_tp, tp_norm_mean, tp_norm_sd))]
  
  # Return
  return(calc[, .(hs, tp)])
}
