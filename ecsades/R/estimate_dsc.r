# DSC ---------------------------------------------------------------------

estimate_dsc = function(
  jdistr, sample_data = NULL, sample_data_npy = NULL,
  output_rp, n_point=100, standardize = TRUE){
  
  ## Generate sample data
  if(is.null(sample_data)){
    sample_data = .sample_jdistr(jdistr = jdistr, sim_year = max(output_rp)*.rp_multiplier)  
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


