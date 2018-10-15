fit_wln = function(data, npy){
  
  res = list()
  class(res) = "wln"
  res$npy = npy
  res$hs = .fit_weibull(data = data$hs)
  res$tp = .fit_iform_lnorm(hs = data$hs, tp = data$tp)
  return(res)
}


# Conditional log-normal used for IFROM -----------------------------------

.fit_iform_lnorm = function(hs, tp){
  
  log_tp = log(tp)
  # theta = c(d, e, f, g, k, m, q)
  theta0 = c(d=mean(log_tp), e=0, f=0, g=sd(log_tp), k=0, m=0, q=0)
  op = optim(
    par = theta0,
    fn = .nll_iform_lnorm,
    log_tp = log_tp, hs = hs, control = list(maxit=2e3))
  
  res = list()
  class(res) = "iform_lnorm"
  res$par = op$par
  res$conv = op$convergence
  
  return(res)
}

.nll_iform_lnorm = function(theta, log_tp, hs){
  
  # Constraints
  if(theta[3] <= 0 | theta[6]>0 | theta[4]<0 | theta[4]+theta[5]<0){
    return(.limit_inf)
  }
  
  # mean & sd
  norm_mean = theta[1] + theta[2] * log(hs + theta[3])
  norm_var = theta[4] + theta[5] * exp(theta[6] * (hs ^ theta[7]))
  if (any(norm_var <= 0)){
    return(.limit_inf)
  }
  
  # return
  nll = -dnorm(log_tp, norm_mean, sqrt(norm_var), log=TRUE)
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