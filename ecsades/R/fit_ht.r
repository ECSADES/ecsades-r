#' Fitting Heffernan-Tawn model to wave data

#' @description 
#' This function fits a Heffernan-Tawn (\code{ht}) model to the given wave data. The marginal distributions are
#' a mixture model with an empirical distribution below a specified threshold, and the generalised Pareto distribution
#' (GPD) above the threshold. The dependence structure is the conditional extreme value model proposed in
#' Heffernan and Tawn (2004).
#' 
#' @param data the wave data in the form a \code{data.table} with wave height \code{hs} and wave period
#' \code{tp} as columns
#' 
#' @param npy the number of data points per year, usually estimated by the number of rows in the
#' supplied wave data divided by the total period of data coverage (in years)
#' 
#' @param margin_thresh_count an integer indicating the number of exceedances strcitly above the threshold
#' used for fitting the Generalised Pareto Distribution (GPD) for the upper tail of the marginal distribution
#' for \code{hs} and \code{tp}
#' 
#' @param dep_thresh_count an integer indicating the number of exceedances (of the conditioning variable)
#' strcitly above the threshold used for fitting the conditional multivariate extreme value model of
#' Heffernan-Tawn (2004)
#' 
#' @details
#' The input \code{data} must be a \code{data.table} object with \code{hs} and \code{tp} columns. This can be
#' generated by reading a CSV file using function \code{\link[data.table]{fread}}.
#' 
#' By the asymptotic properties of the GPD and the HT model, the model bias is only negligible when the 
#' thresholds are sufficiently high. This function prompts the user with the exceedance probability implied by
#' their choie of the arguments \code{margin_thresh_count} and \code{dep_thresh_count} to remind them that
#' the exceedance probability needs to be quite small for the fit of GPD and HT to be valid.
#' 
#' @return An joint distribution object of class \code{ht} containing the key information of a fitted
#' Heffernan-Tawn model, including the marginal distribution parameters and the dependence model parameters.
#'
#' @examples
#' # Load data
#' data(ww3_pk)
#' 
#' # Fit Heffernan-Tawn model
#' ht = fit_ht(data = ww3_pk, npy = nrow(ww3_pk)/10, margin_thresh_count = 100, dep_thresh_count = 100)
#' 
#' @references Heffernan, Janet & A. Tawn, Jonathan. (2004). A Conditional Approach for Multivariate Extreme Values.
#' Journal of the Royal Statistical Society Series B. 66. 497-546. 10.1111/j.1467-9868.2004.02050.x. 
#' 
#' @seealso \code{\link{fit_wln}}, \code{\link{sample_jdistr}}
#' 
#' @export
fit_ht = function(data, npy, margin_thresh_count, dep_thresh_count){
  
  data = as.data.table(data)
  res = list()
  class(res) = "ht"
  res$npy = npy
  
  # Marginal models
  hs_thresh = sort(data$hs, decreasing = T)[margin_thresh_count+1]
  tp_thresh = sort(data$tp, decreasing = T)[margin_thresh_count+1]
  p_margin_thresh = 1-margin_thresh_count/(data[,.N])
  
  res$margin = list(
    p_margin_thresh = p_margin_thresh,
    hs = .fit_gpd_emp(data$hs, thresh = hs_thresh),
    tp = .fit_gpd_emp(data$tp, thresh = tp_thresh))
  
  # Dep models
  dep_data = copy(data)
  dep_data[, u_hs:=rank(hs, ties.method = "random")/(.N+1)]
  dep_data[, u_tp:=rank(tp, ties.method = "random")/(.N+1)]
  dep_data[u_hs<0.5, l_hs:=log(2*u_hs)]
  dep_data[u_hs>=0.5, l_hs:=-log(2-2*u_hs)]
  dep_data[u_tp<0.5, l_tp:=log(2*u_tp)]
  dep_data[u_tp>=0.5, l_tp:=-log(2-2*u_tp)]

  op_hs = nlminb(
    start = c(1,0), objective = .nll_ht_pair,
    lower = c(-1, -100), upper = c(1, 1-.limit_zero),
    cond_var = dep_data[frank(u_hs)>.N-dep_thresh_count, l_hs],
    dep_var = dep_data[frank(u_hs)>.N-dep_thresh_count, l_tp])
  
  op_tp = nlminb(
    start = c(1,0), objective = .nll_ht_pair,
    lower = c(-1, -100), upper = c(1, 1-.limit_zero),
    cond_var = dep_data[frank(u_tp)>.N-dep_thresh_count, l_tp],
    dep_var = dep_data[frank(u_tp)>.N-dep_thresh_count, l_hs])
  
  p_dep_thresh = 1-dep_thresh_count/(data[,.N])
  res$dep = list(
    p_dep_thresh = p_dep_thresh,
    hs = list(
      par = c(a = op_hs$par[1], b = op_hs$par[2]),
      resid = dep_data[frank(u_hs)>.N-dep_thresh_count, (l_tp-op_hs$par[1]*l_hs)/(l_hs^op_hs$par[2])]), 
    tp = list(
      par = c(a = op_tp$par[1], b = op_tp$par[2]),
      resid = dep_data[frank(u_tp)>.N-dep_thresh_count, (l_hs-op_tp$par[1]*l_tp)/(l_tp^op_tp$par[2])]),
    dep_data = dep_data)

  # Return
  message(paste0(
    "The input margin_thresh_count corresponds to an exceedance probability of ",
    signif(1-p_margin_thresh, 2),"."))
  message(paste0(
    "The input dep_thresh_count corresponds to an exceedance probability of ",
    signif(1-p_dep_thresh, 2),"."))
  return(res)
}


.nll_ht_pair = function(theta, cond_var, dep_var){
  a = theta[1]
  b = theta[2]
  resid = (dep_var-a*cond_var)/(cond_var^b)
  mean_resid = mean(resid)
  sd_resid = sd(resid)
  
  nll = -dnorm(
    x = dep_var, log = TRUE,
    mean = a*cond_var+(cond_var^b)*mean_resid,
    sd = (cond_var^b)*sd_resid)
  
  return(sum(nll))
}

.fit_gpd_emp = function(data, thresh){
  gpd = evd::fpot(data, threshold = thresh, cmax = F, std.err = F)
  
  res = list()
  class(res) = "gpd"
  res$par = c(loc = thresh, gpd$param[1], gpd$param[2])
  res$emp = data
  res$conv = gpd$convergence
  
  return(res)
}

.convert_unif_to_lap = function(unif){
  calc = data.table(unif=unif, lap=NA_real_)
  calc[unif>0.5, lap:=-log(2-2*unif)]
  calc[unif<=0.5, lap:=log(2*unif)]
  return(calc$lap)
}

.convert_lap_to_unif = function(lap){
  calc = data.table(lap=lap, unif=NA_real_)
  calc[lap>0, unif:=1-.5*exp(-lap)]
  calc[lap<=0, unif:=.5*exp(lap)]
  return(calc$unif)
}






