#' Estimating the generalised joint exceedance (GJE) environmental contours
#'
#' @description
#' This function estimates the generalised joint exceedance (GJE) environmental contours, as proposed in Jonathan et al. (2014), for a
#' given sample data set or a fitted joint distribution of class \code{ht} or \code{wln}.
#' 
#' @param object either a fitted joint distribution of class \code{wln} or \code{ht}, or an existing sample data
#' set containing \code{hs} and \code{tp} as columns. See details.
#' 
#' @param output_rp the required return periods (in years) for the estimated contours.
#' 
#' @param type the type of joint exceedance definition. The default value is 2. See details.
#' 
#' @param npy the (optional) number of data points per year. This argument is only needed if the input
#' \code{object} is a sample data set.
#' 
#' @param n_point the number of points to output around each contour. The default value is 100.
#' 
#' @details
#' The joint exceedance contour is a collection of points corresponding to a constant joint exceedance proability.
#' The definition is one-sided or two-sided depending on the input argument \code{type}:
#' 
#' \itemize{
#' 
#'   \item If \code{type = 1}, the \eqn{T}-year contour passes through any point \eqn{(hs, tp)} such that the 
#'   joint exceedance probability \eqn{P(Hs > hs & Tp > tp)} is one in \eqn{T} years.
#'   
#'   \item If \code{type = 2}, the \eqn{T}-year contour passes through an anchor point \eqn{(hs0, tp0)}
#'   where \eqn{hs0} is the (\eqn{T/2})-year wave heights, and \eqn{tp0} is the conditional median of the
#'   wave period given the wave height is \eqn{hs0}. The joint exceedance probability is defined as
#'   \eqn{P(Hs > hs & Tp > tp)} for all \eqn{tp > tp0} and \eqn{P(Hs > hs & Tp < tp)} for all \eqn{tp < tp0}.
#' 
#' }
#' 
#' The function can be applied to an existing sample set, which can be generated by function \code{\link{sample_jdistr}}
#' or otherwise by importing a CSV file using function \code{\link[data.table]{fread}}. The input sample data
#' object must contain \code{hs} and \code{tp} as columns. The user must also specify the number of points per year
#' \code{npy}.
#' 
#' Alternatively this function can be applied to a fitted joint distribution object
#' generated by function \code{\link{fit_ht}} or \code{\link{fit_wln}}.  When applied to a
#' \code{wln} or \code{ht} object, the contour estimation makes use of the importance
#' sampling technique proposed in Huseby et. al. (2014) to improve the efficiency of the function.
#' 
#' 
#' @return
#' A set of estimated GJE contours with the specified return periods in the format
#' of a \code{data.table} with \code{rp}, \code{hs}, and \code{tp} as columns.
#'
#' @examples
#' # Estimating type-1 GJE contours based on a fitted model
#' data(ww3_pk)
#' wln = fit_wln(data = ww3_pk, npy = nrow(ww3_pk)/10)
#' ec_wln = estimate_gje(object = wln, output_rp = c(10,100,1000,10000), type = 1)
#' plot_ec(ec = ec_wln, raw_data = ww3_pk)
#' 
#' # Estimating type-2 GJE contours based on a sample data set
#' data(ww3_ts)
#' ec_data = estimate_gje(object = ww3_ts, output_rp = c(.5, 1, 2), type = 2, npy = ww3_ts[, .N/10])
#' plot_ec(ec = ec_data, raw_data = ww3_ts)
#' 
#' @references
#' Jonathan, P., Ewans, K., Flynn, J., 2014. On the estimation of ocean engineering design contours.
#' ASME J. Offshore Mech. Arct. Eng. 136:041101.
#' 
#' Huseby, A., Vanem, E., Natvig, B., 2014. A new Monte Carlo method for environmental contour estimation.
#' Conference Proceedings. DOI:10.1201/b17399-286.
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{fit_wln}}, \code{\link{sample_jdistr}}, \code{\link{plot_ec}}
#' 
#' @export
estimate_gje = function(
  object, output_rp, type = 2, npy = NULL, n_point = 100){
  
  if("wln" %in% class(object)){
    
    res = .estimate_gje_from_wln(wln = object, output_rp = output_rp, type = type, n_point = n_point)
    
    if(!is.null(npy)){
      warning("Argument npy will be imported from the provided joint distribution object.  The supplied npy is ignored.")
    }
    
  }else if("ht" %in% class(object)){
    
    res = .estimate_gje_from_ht(ht = object, output_rp = output_rp, type = type, n_point = n_point)
    
    if(!is.null(npy)){
      warning("Argument npy will be imported from the provided joint distribution object.  The supplied npy is ignored.")
    }
    
  }else if("data.table" %in% class(object)){
    
    if(is.null(npy)){
      stop("Argument npy must be provided for estimating contours from sample data.")
    }
    
    res_list=list()
    for(this_rp in output_rp){
      this_ap_hs = quantile(object$hs, 1-type/npy/this_rp)
      this_contour = .estimate_gje_from_data(
        sample_data = object, n_jex = object[,.N/npy]/this_rp,
        ap_hs = this_ap_hs, type=type, n_point = n_point)  
      res_list[[as.character(this_rp)]] = cbind(rp = this_rp, this_contour)
    }
    res = rbindlist(res_list)
    
  }else{
    stop("The input object must be of class ht, wln or data.table.")
  }
  
  return(res)
}


.estimate_gje_from_data = function(sample_data, n_jex, ap_hs, type, n_point){
  if(!type %in% c(1,2)){
    stop("Only type 1 or 2 is supported.")
  }
  if(type==2){
    calc = data.table(hs1 = seq(ap_hs, min(sample_data$hs), length.out = ceiling((n_point+1)/2)))
    calc[, n_hs_tail:=sample_data[hs>=hs1, .N], .(hs1)]
    calc[, p_tp_upper:=1-n_jex/n_hs_tail]
    calc[, p_tp_lower:=n_jex/n_hs_tail]
    calc[, tp_upper:=sample_data[hs>=hs1, quantile(tp, p_tp_upper)], .(hs1)]
    calc[, tp_lower:=sample_data[hs>=hs1, quantile(tp, p_tp_lower)], .(hs1)]
    out = rbind(calc[.N:2, .(hs=hs1, tp=tp_upper)], calc[, .(hs=hs1, tp=tp_lower)])
  }else{
    # browser()
    calc = data.table(hs1 = seq(ap_hs, min(sample_data$hs), length.out = n_point))
    calc[, tp:=sample_data[hs>=hs1, sort(tp, decreasing = T)[n_jex]], .(hs1)]
    out = calc[, .(hs=hs1, tp)]  
  }
  return(out)
}


.estimate_gje_from_wln = function(wln, output_rp, type, n_point){
  
  res_list=list()
  for(this_rp in output_rp){
    this_sample_data = .sample_wln_is(wln, target_rp=this_rp)
    this_out = .estimate_gje_from_data(
      sample_data = this_sample_data, n_jex = .target_rp_ub,
      ap_hs = this_sample_data[, sort(hs, decreasing = T)[type*.target_rp_ub]],
      type = type, n_point = n_point)
    res_list[[as.character(this_rp)]] = cbind(rp=this_rp, this_out)
  }
  res = rbindlist(res_list)
  return(res)
}

.estimate_gje_from_ht = function(ht, output_rp, type, n_point){
  
  res_list=list()
  for(this_rp in output_rp){
    this_sample_data = .sample_ht_is(ht, target_rp=this_rp)
    this_out = .estimate_gje_from_data(
      sample_data = this_sample_data, n_jex = .target_rp_ub,
      ap_hs = this_sample_data[, sort(hs, decreasing = T)[type*.target_rp_ub]],
      type = type, n_point = n_point)
    res_list[[as.character(this_rp)]] = cbind(rp=this_rp, this_out)
  }
  res = rbindlist(res_list)
}
