#' Estimating the IFORM contours
#'
#' @description
#' This function estimates the Invserse FORM environmental contours for a 
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
#' @param n_point the number of points to output around each contour. The default value is 100.
#' 
#' @param alpha0 the omission factor for the IFORM contours, as proposed in Winterstein et al. (1993).
#' The default value is 0, meaning no inflation to the estimated contours.
#' 
#' @details
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
#' iform = estimate_iform(
#'   jdistr = wln,
#'   output_rp = c(1,10,100),
#'   alpha0 = 0.01)
#'   
#' # Plot output
#' plot_ec(iform, noaa_ww3)
#' 
#' @references
#' Winterstein, S. R., Ude, T. C., Cornell, C. A., Bjerager, P., Haver, S., 1993. Environmental parameters
#' for extreme variable: Inverse Form with omission factors. In: Proc. 6th Int. Conf. on
#' Structural Safety and Reliability, Innsbruck, Austria.
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{fit_wln}}, \code{\link{sample_jdistr}}, \code{\link{plot_ec}}
#' 
#' @export
estimate_iform = function(
  jdistr, sample_data = NULL, sample_data_npy = NULL,
  output_rp, n_point=100, alpha0=0){

  if(is.null(sample_data)){
    
    if(class(jdistr)=="wln"){
      
      res = .estimate_iform_from_wln(wln = jdistr, output_rp = output_rp, n_point = n_point, alpha0 = alpha0)
      
    }else if(class(jdistr)=="ht"){
      
      sim_year = max(output_rp*.rp_multiplier)
      res = .estimate_iform_from_sample_data(
        sample_data = .sample_ht(jdistr, sim_year),
        sample_data_npy = jdistr$npy, output_rp = output_rp, n_point = n_point, alpha0 = alpha0)
      
    }else{
      
      stop("Input class of distribution not")
      
    }
  }else{
    
    res = .estimate_iform_from_sample_data(
      sample_data = sample_data, sample_data_npy = sample_data_npy,
      output_rp = output_rp, n_point = n_point, alpha0 = alpha0)
  }

  return(res)
}


.estimate_iform_from_sample_data = function(sample_data, sample_data_npy, output_rp, n_point, alpha0){

  
  # Calculate beta and inflate
  prob = 1 - 1/(sample_data_npy*output_rp)
  beta = qnorm(prob)
  beta_star = beta / (sqrt(1-(alpha0^2)))
  
  # Define points
  angle = seq(0, 2*pi, length.out = n_point+1)
  calc = data.table(rp = rep(output_rp, each=n_point+1), angle = angle)
  calc[, u1:=rep(beta_star, each=n_point+1)*cos(angle)]
  calc[, u2:=rep(beta_star, each=n_point+1)*sin(angle)]
  calc[, p1:=pnorm(u1)]
  calc[, p2:=pnorm(u2)]
  calc[, hs:=quantile(sample_data$hs, p1, names = F)]
  
  for(i in 1:calc[, .N]){
    cond_distr_bw = calc[rp==calc[i]$rp, bw.SJ(hs)]
    td_data = sample_data[abs(hs-calc$hs[i])<=cond_distr_bw/2]
    if(td_data[,.N]<.knn){
      td_data = sample_data[sort.list(abs(hs-calc$hs[i]))[1:.knn]]
    }
    calc[i, tp:=quantile(td_data$tp, p2, names = F)]
  }
  
  # Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}


.estimate_iform_from_wln = function(wln, output_rp, n_point, alpha0){
  
  # Calculate beta and inflate
  prob = 1 - 1/(wln$npy*output_rp)
  beta = qnorm(prob)
  beta_star = beta / (sqrt(1-(alpha0^2)))
  
  # Define points
  angle = seq(0, 2*pi, length.out = n_point+1)
  calc = data.table(rp = rep(output_rp, each=n_point+1), angle = angle)
  calc[, u1:=rep(beta_star, each=n_point+1)*cos(angle)]
  calc[, u2:=rep(beta_star, each=n_point+1)*sin(angle)]
  calc[, p1:=pnorm(u1)]
  calc[, p2:=pnorm(u2)]
  calc[, hs:=qweibull(p1, shape=wln$hs$par["shape"], scale=wln$hs$par["scale"])+wln$hs$par["loc"]]
  calc[, m:=wln$tp$par[1] + wln$tp$par[2] * (hs ^ wln$tp$par[3])]
  calc[, s:=wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * hs)]
  calc[, tp:=qlnorm(p2, m, s)]
  
  # Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}

