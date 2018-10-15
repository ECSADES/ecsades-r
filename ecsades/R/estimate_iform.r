estimate_iform = function(jdistr, output_rp, n_point, alpha){
  if(class(jdistr)=="ht"){
    res = .estimate_iform_from_ht(ht = jdistr, output_rp = output_rp, n_point = n_point, alpha = alpha)
  }else if(class(jdistr)=="wln"){
    res = .estimate_iform_from_wln(wln = jdistr, output_rp = output_rp, n_point = n_point, alpha = alpha)
  }else{
    stop("Input class of distribution not supported.")
  }
  return(res)
}


.estimate_iform_from_ht = function(ht, output_rp, n_point, alpha){
  
  # Simulate from ht
  sim_year = max(output_rp*.rp_multiplier)
  sample_data = .sample_ht(ht, sim_year)
  
  # Calculate beta and inflate
  prob = 1 - 1/(ht$npy*output_rp)
  beta = qnorm(prob)
  beta_star = beta / (sqrt(1-(alpha^2)))
  
  # Define points
  angle = seq(0, 2*pi, length.out = n_point+1)
  calc = data.table(rp = rep(output_rp, each=n_point+1), angle = angle)
  calc[, u1:=rep(beta_star, each=n_point+1)*cos(angle)]
  calc[, u2:=rep(beta_star, each=n_point+1)*sin(angle)]
  calc[, p1:=pnorm(u1)]
  calc[, p2:=pnorm(u2)]
  calc[, hs:=quantile(sample_data$hs, p1, names = F)]
  
  cond_distr_bw = calc[, diff(range(hs))/(n_point/2)]
  for(i in 1:calc[, .N]){
    tp_hs_distr = sample_data[abs(hs-calc$hs[i])<=cond_distr_bw*2, tp]
    calc[i, tp:=quantile(tp_hs_distr, p2, names=F)]
  }
  
  # Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}


.estimate_iform_from_wln = function(wln, output_rp, n_point, alpha){
  
  # Calculate beta and inflate
  prob = 1 - 1/(wln$npy*output_rp)
  beta = qnorm(prob)
  beta_star = beta / (sqrt(1-(alpha^2)))
  
  
  # Define points
  angle = seq(0, 2*pi, length.out = n_point+1)
  calc = data.table(rp = rep(output_rp, each=n_point+1), angle = angle)
  calc[, u1:=rep(beta_star, each=n_point+1)*cos(angle)]
  calc[, u2:=rep(beta_star, each=n_point+1
                 )*sin(angle)]
  calc[, p1:=pnorm(u1)]
  calc[, p2:=pnorm(u2)]
  calc[, hs:=qweibull(p1, shape=wln$hs$par["shape"], scale=wln$hs$par["scale"])+wln$hs$par["loc"]]
  calc[, m:=wln$tp$par["d"] + wln$tp$par["e"] * log(hs + wln$tp$par["f"])]
  calc[, s:=sqrt(wln$tp$par["g"] + wln$tp$par["k"] * exp(wln$tp$par["m"] * (hs ^ wln$tp$par["q"])))]
  calc[, tp:=qlnorm(p2, m, s)]
  
  # Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}

