#' Estimating the direct sampling contours
#'
#' @description
#' This function estimates the direct sampling contours, as proposed in Huesby et al. (2015), for a given
#' sample data or a fitted joint distribution of class \code{ht} or \code{wln} object.
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
#' @param n_point the number of points to output around each contour. The default value is 100
#' 
#' @param standardize whether or not to apply the recommended standardization as detailed in Husbey et al. (2015).
#' The default value is \code{TRUE}.
#' 
#' @details
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
#' # Estimate the DSC based on fitted Heffernan-Tawn model
#' ht = fit_ht(data = noaa_ww3, npy = nrow(noaa_ww3)/10, p_margin_thresh = 0.95, p_dep_thresh = 0.95)
#' dsc_ht = estimate_dsc(jdistr = ht, output_rp = c(1,10,100))
#' 
#' # Estimate the DSC based on sample data
#' wln = fit_wln(data = noaa_ww3, npy = nrow(noaa_ww3)/10)
#' sample_data = sample_jdistr(jdistr = wln, sim_year = 1000)
#' dsc_data = estimate_dsc(sample_data = sample_data, sample_data_npy = nrow(sample_data)/1000, output_rp = c(1,10,100))
#' 
#' @references
#' Huseby, A., Vanem, E., Natvig, B., 2015. Alternative environmental contours for structural reliability
#' analysis. Struct. Saf. 54, 32-45.
#' 
#' @export
estimate_dsc = function(
  jdistr, sample_data = NULL, sample_data_npy = NULL,
  output_rp, n_point=100, standardize = TRUE){
  
  ## Generate sample data
  if(is.null(sample_data)){
    sample_data = sample_jdistr(jdistr = jdistr, sim_year = max(output_rp)*.rp_multiplier)  
    npy = jdistr$npy
  }else{
    npy = sample_data_npy
  }

  
  ## Standardization
  ec_data = copy(sample_data)
  if(standardize){
    st_mu = colMeans(sample_data)
    st_sigma = apply(sample_data, 2, sd)
    st_rho = cor(ec_data)[1,2]
    ec_data[, hs:= (sample_data$hs-st_mu["hs"])/st_sigma["hs"]]
    ec_data[, tp:= ((sample_data$tp-st_mu["tp"])/st_sigma["tp"]-st_rho*hs)/sqrt(1-st_rho^2)]
  }
    
  ## C_theta calculation
  ex_prob = 1/(npy*output_rp)
  theta = seq(0, 2*pi, length.out = n_point+1)
  d_theta = 2*pi/(n_point)
  calc = data.table(rp = rep(output_rp, each=n_point+1), theta = theta, ex_prob = rep(ex_prob, each=n_point+1))
  calc[, c_theta:=quantile(cos(theta[1])*ec_data$tp+sin(theta[1])*ec_data$hs, 1-ex_prob), .(theta)]
  calc[, c_theta_left:=c_theta[c(.N-1, 1:(.N-2), .N-1)], .(rp)]
  calc[, c_theta_right:=c_theta[c(2:.N, 2)], .(rp)]
  calc[, c_theta_prime:=(c_theta_right-c_theta_left)/(2*d_theta)]
  calc[, tp:=(c_theta*cos(theta)-c_theta_prime*sin(theta))]
  calc[, hs:=(c_theta_prime*cos(theta)+c_theta*sin(theta))]

  ## Rev-standardization
  if(standardize){
    st_data = calc[, .(hs, tp)]
    calc[, hs:=st_sigma["hs"]*st_data$hs+st_mu["hs"]]
    calc[, tp:=st_sigma["tp"]*(st_rho*st_data$hs+st_data$tp*sqrt(1-st_rho^2))+st_mu["tp"]]
  }
  
  ## Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}


