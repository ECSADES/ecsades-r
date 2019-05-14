#' Estimating the isodensity contours
#'
#' @description
#' This function estimates the isodensity environmental contours for a 
#' given sample data or a fitted joint distribution of class \code{ht} or \code{wln} object.
#' 
#' @param jdistr the (optional) fitted joint distribution as generated by function \code{\link{fit_ht}}
#' or \code{\link{fit_wln}}.  See details.
#' 
#' @param sample_data the (optional) existing sample wave data containing \code{hs} and \code{tp} as columns. See details.
#' 
#' @param sample_data_npy the (optional) number of data points per year of the sample data. See details.
#' 
#' @param output_rp the required return periods (in years) for the estimated contours. Note that return periods
#' that are substantially larger than the data coverage may lead to a long processing time
#' 
#' @param n_grid the resolution used for the numerical integration of joint probability density. 
#' The default value is 100. See details.
#' 
#' @param range_multiplier the range in \code{hs} and \code{tp} used for the numerical integration
#' of the joint probability density, as a multiplier of the actual range of the sample data.  See details.
#' 
#' @details
#' The estimated contour for a given return period \code{rp} is a set of points with constant probability
#' density, such that the probability of a point lying in the interior of the contour is \code{1-1/rp}.
#' 
#' The probability is estimated using numerical integration over a region defined by \code{range_multiplier}
#' at a resolution specified by \code{n_grid}. For exmaple, setting \code{range_multiplier = 1.5} and 
#' \code{n_grid = 100} means the joint probability is estimated over a mesh grid of \eqn{100 * 100}
#' from \eqn{min(hs)/1.5} to \eqn{max(hs)*1.5} and similarly for \code{tp}. Higher \code{range_multiplier}
#' and \code{n_grid} mean less estimation error, but longer processign time.
#' 
#' The input \code{sample_data} must be a \code{data.table} object with \code{hs} and \code{tp} columns. This can be
#' generated by reading a CSV file using function \code{\link[data.table]{fread}}.
#' 
#' This function can be applied to either a fitted joint distribution of class \code{ht} or \code{wln}. This
#' happens when the user supplies a valid \code{ht} or \code{wln} object, and ignores arguments \code{sample_data}
#' and \code{sample_data_npy}.
#' 
#' This function can also be applied to an existing sample data. To do this, the user needs
#' to supply both the sample data \code{sample_data} and the corresponding number of data points per year
#' \code{sample_data_npy}. The sample data can be generated using function \code{\link{sample_jdistr}} or otherwise, 
#' but must contain \code{hs} and \code{tp} as columns.
#' 
#' The current implementation of the estimation method cannot efficiently deal with the lower-left quadrant
#' of a fitted \code{ht} object where both \code{hs} and \code{tp} are below the dependence treshold.
#' The resulting output contours do not include that quadrant to avoid excessive computation.
#' 
#' Providing \code{sample_data} and \code{sample_data_npy} will make the input \code{jdistr} redundant.
#' 
#' @return A set of estimated environmental contours with the specified return periods in the format
#' of a \code{data.table} with \code{rp}, \code{hs}, and \code{tp} as columns.
#'
#' @examples
#' # Load data
#' data(noaa_ww3)
#' 
#' # Fit the Weibull-log-normal distribution to hs and tp
#' wln = fit_wln(data = noaa_ww3, npy = nrow(noaa_ww3)/10)
#' 
#' # Estimate the GJE contours using the Hs return levels as 
#' iso = estimate_iso(jdistr = wln,  output_rp = c(1,10,100))
#'   
#' # Plot output
#' plot_ec(iso, noaa_ww3)
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{fit_wln}}, \code{\link{sample_jdistr}}, \code{\link{plot_ec}}
#' 
#' @export
estimate_iso = function(jdistr, output_rp, n_grid=100, range_multiplier=1.5){
  
  if(class(jdistr)=="ht"){
    res = .estimate_iso_from_ht(ht = jdistr, output_rp = output_rp, n_grid, range_multiplier)
  }else if(class(jdistr)=="wln"){
    res = .estimate_iso_from_wln(wln = jdistr, output_rp = output_rp, n_grid, range_multiplier)
  }else{
    stop("Input class of distribution not supported.")
  }
  return(res)
}


.estimate_iso_from_ht = function(ht, output_rp, n_grid, range_multiplier){
  
  # Set up grid on Laplace scale
  max_u = 1-1/max(output_rp)/ht$npy
  max_lap = -log(2-2*max_u)*range_multiplier
  thresh_lap = -log(2-2*ht$dep$p_dep_thresh)
  grid_lap = seq(-max_lap, max_lap, length.out = n_grid*2+1)
  calc = data.table(l_hs = grid_lap, l_tp = rep(grid_lap, each=n_grid*2+1))

  # Upper 
  bw_resid_hs = bw.SJ(ht$dep$hs$resid)
  bw_resid_tp = bw.SJ(ht$dep$tp$resid)
  calc[l_hs>=l_tp & l_hs>=thresh_lap, resid_hs:= (l_tp-l_hs*ht$dep$hs$par["a"])/(l_hs^ht$dep$hs$par["b"])]
  calc[
    !is.na(resid_hs),
    f:= 1/l_hs^ht$dep$hs$par["b"] *
      mean(dnorm(resid_hs-ht$dep$hs$resid, 0, bw_resid_hs))*.5*exp(-abs(l_hs)),
    .(l_hs, l_tp)]
  
  calc[l_hs<l_tp & l_tp>=thresh_lap, resid_tp:= (l_hs-l_tp*ht$dep$tp$par["a"])/(l_tp^ht$dep$tp$par["b"])]
  calc[
    !is.na(resid_tp),
    f:= 1/l_tp^ht$dep$tp$par["b"] *
      mean(dnorm(resid_tp-ht$dep$tp$resid, 0, bw_resid_tp))*.5*exp(-abs(l_tp)),
    .(l_hs, l_tp)]
  calc[is.na(f), f:=0]
  
  # Estimate contous on Laplace scale
  sorted_f = sort(calc$f, decreasing = T)
  ex_prob = ht$dep$dep_data[, mean(l_hs>=thresh_lap | l_tp>=thresh_lap)]
  prob = cumsum(sorted_f)/sum(sorted_f)*ex_prob+(1-ex_prob)
  dens_level = sorted_f[sapply(output_rp, function(rp)sum(prob<=(1-1/rp/ht$npy)))]
  contour_lines = contourLines(x=grid_lap, y=grid_lap, z=matrix(calc$f, n_grid*2+1, n_grid*2+1), levels=dens_level)
  contour_lap = rbindlist(lapply(
    X = contour_lines,
    FUN = function(dat)data.table(level=dat$level, l_hs=dat$x, l_tp=dat$y)))
  
  # Convert to original scale
  contour_lap[, u_hs:=(1+sign(l_hs))/2-(sign(l_hs))/2*exp(-abs(l_hs))]
  contour_lap[, u_tp:=(1+sign(l_tp))/2-(sign(l_tp))/2*exp(-abs(l_tp))]
  n = length(ht$margin$hs$emp)
  contour_lap[, hs:=quantile(ht$margin$hs$emp, u_hs)]
  contour_lap[u_hs>ht$margin$p_margin_thresh, u_gpd:=(u_hs-ht$margin$p_margin_thresh)/(1-ht$margin$p_margin_thresh)]
  contour_lap[u_hs>ht$margin$p_margin_thresh, hs:=evd::qgpd(
    p = u_gpd, loc = ht$margin$hs$par[1],
    scale = ht$margin$hs$par[2], shape = ht$margin$hs$par[3])]
  contour_lap[, u_gpd:=NULL]
  
  contour_lap[, tp:=quantile(ht$margin$tp$emp, u_tp)]
  contour_lap[u_tp>ht$margin$p_margin_thresh, u_gpd:=(u_tp-ht$margin$p_margin_thresh)/(1-ht$margin$p_margin_thresh)]
  contour_lap[u_tp>ht$margin$p_margin_thresh, tp:=evd::qgpd(
    p = u_gpd, loc = ht$margin$tp$par[1],
    scale = ht$margin$tp$par[2], shape = ht$margin$tp$par[3])]
  
  # Return
  res = merge(
    contour_lap,
    data.table(level=dens_level, rp=output_rp))[l_hs>=thresh_lap | l_tp>=thresh_lap, .(rp, hs, tp)]
  return(res)
}


.estimate_iso_from_wln = function(wln, output_rp, n_grid, range_multiplier){
  
  # Generate sample data
  sample_data = .sample_wln(wln, sim_year = max(output_rp)*.rp_multiplier) 
  
  # Calculate density
  dens = .dwln(sample_data$hs, sample_data$tp, hs_par = wln$hs$par, tp_par = wln$tp$par)
  dens_level = quantile(dens, 1/(output_rp*wln$npy), names = F)
  
  # Contour lines
  grid_hs = sample_data[, seq(min(hs)/range_multiplier, max(hs)*range_multiplier, length.out = n_grid)]
  grid_tp = sample_data[, seq(min(tp)/range_multiplier, max(tp)*range_multiplier, length.out = n_grid)]
  grid_dens = .dwln(
    hs = matrix(grid_hs, n_grid, n_grid),
    tp = t(matrix(grid_tp, n_grid, n_grid)),
    hs_par = wln$hs$par, tp_par = wln$tp$par)
  
  contour_lines = contourLines(x=grid_hs, y=grid_tp, z=grid_dens, levels=dens_level)
  
  ## Return
  calc = rbindlist(lapply(
    X = contour_lines,
    FUN = function(dat)data.table(level=dat$level, hs=dat$x, tp=dat$y)))
  res = merge(calc, data.table(level=dens_level, rp=output_rp))[, .(rp, hs, tp)]
  return(res)
}

.dwln = function(hs, tp, hs_par, tp_par){
  f_hs = dweibull(hs-hs_par["loc"], shape=hs_par["shape"], scale=hs_par["scale"])
  norm_mean = tp_par[1] + tp_par[2] * log(hs + tp_par[3])
  norm_var = tp_par[4] + tp_par[5] * exp(tp_par[6] * (hs ^ tp_par[7]))
  f_tp_hs = dlnorm(tp, norm_mean, sqrt(norm_var))
  return(f_hs*f_tp_hs)
}
