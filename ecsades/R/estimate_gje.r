#' Estimating the generalised joint exceedance contours
#'
#' @description
#' This function estimates the generalised joint exceedance contours, as proposed in Jonathan et al. (2014), for a
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
#' that are substantially larger than the data coverage may lead to a long processing time.
#' 
#' @param n_point the number of points to output around each contour. The default value is 100.
#' 
#' @param ref_tp the wave period for the reference point. The default value is 0. See details.
#' 
#' @param ref_hs the wave height for the reference point. The default value is 0. See details.
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
#' For each return period of the output contours, the corresponding reference point (\code{hs}, \code{tp}) divides
#' the space into (up to) four quadrants.  The joint exceedance region of a given point
#' is then the corner formed by the point as the vertex, two semi-straight lines parallel to the \code{hs} and \code{tp}
#' axes respectively, and facing away from the reference point.
#' 
#' Please refer to Jonathan et al. (2014) for the mathematical formulation of the definition.
#' 
#' @return A set of estimated environmental contours with the specified return periods in the format
#' of a \code{data.table} with \code{rp}, \code{hs}, and \code{tp} as columns.
#'
#' @examples
#' # Load data
#' data(noaa_ww3)
#' 
#' # Generate a 2000-year sample data based on fitted HT
#' ht = fit_ht(data = noaa_ww3, npy = nrow(noaa_ww3)/10, p_margin_thresh = 0.95, p_dep_thresh = 0.95)
#' n_sim_year = 2000
#' sample_data = sample_jdistr(jdistr = ht, sim_year = n_sim_year)
#' 
#' # Estimate the 10, 100, 1000-year return levels of Hs based on the simulated data
#' out_rp = c(1, 10, 100)
#' rl_hs = sample_data[sort.list(hs, decreasing = T)[n_sim_year/out_rp], hs]
#'
#' # Estimate the median of the Tp given Hs is equal to these return levels
#' ref_tp = c()
#' for(i in 1:3){
#'   ref_tp[i] = sample_data[abs(hs-rl_hs[i])<=1, median(tp)]
#' }
#' 
#' # Estimate the GJE contours using the Hs return levels as 
#' gje_data = estimate_gje(
#'   sample_data = sample_data,
#'   sample_data_npy = nrow(sample_data)/n_sim_year,
#'   output_rp = out_rp,
#'   ref_tp = ref_tp)
#'   
#' # Plot output
#' plot_ec(gje_data, noaa_ww3)
#' 
#' @references
#' Jonathan, P., Ewans, K., Flynn, J., 2014. On the estimation of ocean engineering design contours.
#' ASME J. Offshore Mech. Arct. Eng. 136:041101.
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{fit_wln}}, \code{\link{sample_jdistr}}, \code{\link{plot_ec}}
#' 
#' @export
estimate_gje = function(
  jdistr, sample_data = NULL, sample_data_npy = NULL,
  output_rp,
  n_point = 100,
  ref_tp = rep(0, length(output_rp)),
  ref_hs = rep(0, length(output_rp))){
  
  ## Generate sample data
  if(is.null(sample_data)){
    sample_data = sample_jdistr(jdistr = jdistr, sim_year = max(output_rp)*.rp_multiplier)  
    npy = jdistr$npy
  }else{
    sample_data = copy(sample_data)
    npy = sample_data_npy
  }
  sample_data[, q:=NA_character_]
  
  ## Main loop
  res = list()
  for(i in 1:length(output_rp)){
    this_rp = output_rp[i]
    
    ## Define quadrants
    sample_data[, ":="(c("qhs", "qtp"), list(abs(hs-ref_hs[i]), abs(ref_tp[i]-tp)))]
    sample_data[hs>ref_hs[i] & tp>ref_tp[i], q:="hh"]
    sample_data[hs>ref_hs[i] & tp<=ref_tp[i], q:="hl"]
    sample_data[hs<=ref_hs[i] & tp>ref_tp[i], q:="lh"]
    sample_data[hs<=ref_hs[i] & tp<=ref_tp[i], q:="ll"]
    
    ## Define contours per quadrant
    target_prob = 1/(this_rp*npy)
    target_n = target_prob*sample_data[, .N]
    list_out = list()
    n_q = sample_data[, uniqueN(q)]
    
    for(this_q in unique(sample_data$q)){
      qdata = sample_data[q==this_q]
      res_hs = lapply(
        X = qdata[, seq(min(qhs), quantile(qhs, 1-target_prob), length.out=n_point/n_q)],
        FUN = .find_gje_point_hs, target_n = target_n, qdata = qdata)
      res_tp = lapply(
        X = qdata[, seq(min(qhs), quantile(qtp, 1-target_prob), length.out=n_point/n_q)],
        FUN = .find_gje_point_tp, target_n = target_n, qdata = qdata)
      list_out[[this_q]] = cbind(rbind(rbindlist(res_hs), rbindlist(res_tp)), q=this_q)
    }
    this_out = cbind(rp=this_rp, rbindlist(list_out))
    
    ## Convert (qhs, qtp) to original hs, tp
    this_out[q=="hh", ":="(c("hs", "tp"), list(qhs+ref_hs[i],qtp+ref_tp[i]))]
    this_out[q=="hl", ":="(c("hs", "tp"), list(qhs+ref_hs[i],-qtp+ref_tp[i]))]
    this_out[q=="lh", ":="(c("hs", "tp"), list(-qhs+ref_hs[i],qtp+ref_tp[i]))]
    this_out[q=="ll", ":="(c("hs", "tp"), list(-qhs+ref_hs[i],-qtp+ref_tp[i]))]
    this_out = this_out[sort.list(atan2(x=hs-ref_hs[i], y=tp-ref_tp[i]))]
    res[[i]] = this_out
  }
  
  ## Return
  res = rbindlist(res)[,.(rp,hs,tp)]
  return(res)
}

.find_gje_point_hs = function(qhs0, target_n, qdata){
  ub_n_qdata = qdata[, sum(qhs>qhs0)]
  idx_valid = target_n<=ub_n_qdata
  
  qtp0 = qdata[qhs>qhs0, quantile(qtp, 1-target_n[idx_valid]/ub_n_qdata)]
  res = data.table(qhs=rep(qhs0, sum(idx_valid)), qtp=qtp0)
  return(res)
}

.find_gje_point_tp = function(qtp0, target_n, qdata){
  ub_n_qdata = qdata[, sum(qtp>qtp0)]
  idx_valid = target_n<=ub_n_qdata
  
  qhs0 = qdata[qtp>qtp0, quantile(qhs, 1-target_n[idx_valid]/ub_n_qdata)]
  res = data.table(qhs = qhs0, qtp=rep(qtp0, sum(idx_valid)))
  return(res)
}


