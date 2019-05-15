#' Fitting Weibull-log-normal model to wave data
#'
#' @description 
#' This function fits a Weibull-log-normal (\code{ht}) model to the given wave data, such that the wave height
#' \code{hs} follows a translated (or 3-parameter) Weibull distribution, and the wave period given the wave
#' follows a conditional log-normal distributuion with the location and scale parameters as functions
#' of the corresponding \code{hs} value.
#' 
#' @param data the wave data in the form a \code{data.table} with wave height \code{hs} and wave period
#' \code{tp} as columns.
#' 
#' @param npy the number of data points per year, usually estimated by the number of rows in the
#' supplied wave data divided by the total period of data coverage (in years).
#' 
#' @param hs_constraint_range a multiplier to the maximum \code{hs} in the data within such that the standard
#' deviation of the log-normal conditional distribution of \code{tp} given \code{hs} is constrained to be strictly
#' positive (see details for more information) over this user-defined range. The default value is 1.5.
#' 
#' @param weighted_tp_fit whether or not the estimation of the conditional model for \code{tp} uses a weighted
#' likelihood (see details). The default option is \code{FALSE}.
#' 
#' @details
#' The input \code{data} must be a \code{data.table} object with \code{hs} and \code{tp} columns. This can be
#' generated by reading a CSV file using function \code{\link[data.table]{fread}}.
#' 
#' The formulation of the conditional distribution is given by
#' \deqn{log(tp | hs=h) ~ N(\mu(h), \sigma(h)^2)}
#' where the mean and the standard deviation are
#' \deqn{\mu(h) = a_0 + a_1 h^a_2 and \sigma(h) = b_0 + b_1 exp(h * b_2)}
#' 
#' The strictly positive constraint for the standard deviation has a large influence on the estimated values for
#' coefficients \eqn{b_0}, \eqn{b_1} and \eqn{b_2}. In general, enforcing the standardard deviation to be strictly
#' positive over a larger range (larger \code{hs_constraint_range}) may result in a poorer fit of the model
#' to the provided observation data.
#' 
#' There is an option to use weighted likelihood when fitting the log-normal distribution to \code{tp|hs}. When this
#' option is used, the observation data are divided to bins of equal sizes (every whole number in \code{hs}) and the
#' weight for each data point is inversely proportional to the total number of points within the same bin.  This
#' option effectively weighs up points in the tail of the distribution, but weighs down the points in the centre
#' of the distribution.
#' 
#' @return A joint distribution object of class \code{wln} containing the key information of a fitted
#' Weibull-log-normal model, including the three parameters of the Weibull distribution for \code{hs} (as a named
#' numeric vector) and the six parameters of the conditional log-normal distribution for \code{tp} given \code{hs}
#' (as an unamed numeric vector, in the order of \eqn{{a_0, a_1, a_2, b_0, b_1, b_2}}).
#'
#' It is possible to replace one or multiple \code{hs} or \code{tp} parameters in an existing \code{wln} object
#' (see the examples provided below).  Note the resulting object may violate the input \code{hs_constraint_range}.
#' The current version of the package does not automatically check the validity of these user-input parameters.
#' Therefore users are advised to perform the check independently.
#' 
#' @examples
#' # Load data
#' data(ww3_pk)
#' 
#' # Fit Weibull-log-normal distribution 
#' wln1 = fit_wln(data = ww3_pk, npy = nrow(ww3_pk)/10)
#' 
#' # Fit Weibull-log-normal distribution with additional options
#' wln2 = fit_wln(data = ww3_pk, npy = nrow(ww3_pk)/10, hs_constraint_range = 3, weighted_tp_fit = TRUE)
#' 
#' # Update the hs parameters for object wln1
#' wln3 = copy(wln1)
#' wln3$hs$par[["loc"]] = 0.66
#' wln3$hs$par[["scale"]] = 2.2
#' wln3$hs$par[["shape"]] = 1.8
#' 
#' # Update the tp parameters for object wln2
#' wln4 = copy(wln2)
#' wln4$tp$par[1] = 2.5
#' wln4$tp$par[2] = -0.015
#' 
#' 
#' @references 
#' Haver, Sverre & Winterstein, Steven. (2009). Environmental Contour Lines: A Method for Estimating Long
#' Term Extremes by a Short Term Analysis. Transactions - Society of Naval Architects and Marine Engineers. 116. 
#' 
#' @seealso \code{\link{fit_ht}}, \code{\link{sample_jdistr}}
#' 
#' @export
fit_wln = function(data, npy, hs_constraint_range = 1.5, weighted_tp_fit = FALSE){
  
  res = list()
  class(res) = "wln"
  res$npy = npy
  res$hs = .fit_weibull(data = data$hs)
  res$tp = .fit_iform_lnorm(
    hs = data$hs, tp = data$tp,
    hs_constraint_range = hs_constraint_range, weighted_tp_fit = weighted_tp_fit)
  return(res)
}


# Conditional log-normal used for IFROM -----------------------------------

.fit_iform_lnorm = function(hs, tp, hs_constraint_range, weighted_tp_fit){
  
  log_tp = log(tp)
  theta0 = c(mean(log_tp), 0, 0, sd(log_tp), 0, 0)
  op = optim(
    par = theta0,
    fn = .nll_iform_lnorm,
    log_tp = log_tp, hs = hs,
    hs_constraint_range = hs_constraint_range,
    weighted_tp_fit = weighted_tp_fit,
    control = list(maxit=2e3))
  
  res = list()
  class(res) = "iform_lnorm"
  res$par = op$par
  res$conv = op$convergence
  
  return(res)
}

.nll_iform_lnorm = function(theta, log_tp, hs, hs_constraint_range, weighted_tp_fit){
  
  # Constraints
  mean_range = theta[1] + theta[2] * c(max(hs)*hs_constraint_range, .limit_zero)^theta[3]
  sd_range = theta[4] + theta[5] * exp(theta[6] * c(max(hs)*hs_constraint_range, .limit_zero))
  
  if(any(is.na(sd_range)) || any(is.na(mean_range)) || min(sd_range)<.limit_zero){
    return(.limit_inf)
  }
  
  # mean & sd
  norm_mean = theta[1] + theta[2]* (hs^theta[3])
  norm_sd = theta[4] + theta[5] * exp(theta[6] * hs)
  if (any(norm_sd <= 0)){
    return(.limit_inf)
  }
  
  # Weights
  if(weighted_tp_fit){
    freq = data.table(hs, rhs =round(hs*2))
    freq[, weight:=1/.N, .(rhs)]
    weight = freq$weight
  }else{
    weight = 1
  }
  
  # return
  nll = -dnorm(log_tp, norm_mean, norm_sd, log=TRUE)*weight - log_tp*weight
  res = min(.limit_inf, sum(nll))
  return(res)
}


# Weibull fit ----------------------------------------------------

.fit_weibull = function(data, method){
  theta0 = c(min(data-min(data))/2, sd(data), 1)
  op = nlminb(
    start = theta0,
    objective = .nll_weibull3, data = data,
    lower = c(.limit_zero, .limit_zero, .limit_zero),
    upper = c(min(data)-.limit_zero, .limit_inf, .limit_inf))
  
  res = list()
  res$par = c(loc=op$par[1], scale=op$par[2], shape=op$par[3])
  res$conv = op$convergence
  class(res) = "weibull"
  return(res)
}

.nll_weibull3 = function(theta, data) {
  if(theta[1]>=min(data) || theta[2]<=0 || theta[3]<=0){
    res = .limit_inf
  }else{
    nll = dweibull(data-theta[1], scale=theta[2], shape=theta[3], log = TRUE)
    res = -sum(nll)
  }
  return(res)
}