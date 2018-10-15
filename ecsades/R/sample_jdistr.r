# Sampling function -------------------------------------------------------

.sample_jdistr = function(jdistr, sim_year){
  
  if(class(jdistr)=="ht"){
    res = .sample_ht(jdistr, sim_year)
  }else if(class(jdistr)=="wln"){
    res = .sample_wln(jdistr, sim_year)
  }else{
    stop("Input class of distribution not supported.")
  }
  return(res)
}


# Heffernan-Tawn ----------------------------------------------------------

.sample_ht = function(ht, sim_year){
  
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
  dep_data_hs = .sample_ht1(n_hs, ht$dep$hs$par, ht$dep$hs$resid, ht$dep$p_dep_thresh)
  res[u_hs>ht$dep$p_dep_thresh & u_hs>u_tp, c("u_hs", "u_tp"):=dep_data_hs[, .(u_cond, u_dep)]]
  dep_data_tp = .sample_ht1(n_tp, ht$dep$tp$par, ht$dep$tp$resid, ht$dep$p_dep_thresh)
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

.sample_ht1 = function(n_sim, par, resid, p_dep_thresh){
  set.seed(.seed_sampling)
  a = par[1]
  b = par[2]
  sim_data = data.table(u_cond=runif(n_sim, p_dep_thresh, 1))
  sim_data[u_cond<0.5, l_cond:=log(2*u_cond)]
  sim_data[u_cond>=0.5, l_cond:=-log(2-2*u_cond)]
  resid_max = sim_data[, (1-a)*l_cond^(1-b)]
  extended_resid = resid+rnorm(length(resid)*100, 0, bw.SJ(resid))
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
  tp_norm_mean = wln$tp$par[1] + wln$tp$par[2] * log(hs + wln$tp$par[3])
  tp_norm_var = wln$tp$par[4] + wln$tp$par[5] * exp(wln$tp$par[6] * (hs ^ wln$tp$par[7]))
  tp = exp(rnorm(n_sim, tp_norm_mean, sqrt(tp_norm_var)))
  
  res = data.table(hs, tp)
  return(res)
}

