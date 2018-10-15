# JEC ---------------------------------------------------------------------

estimate_jec = function(jdistr, output_rp, n_point, ref_tp=0, ref_hs=0){
  
  ## Generate sample data
  sample_data = .sample_jdistr(jdistr = jdistr, sim_year = max(output_rp)*.rp_multiplier)  
  

  ## Return
  res = calc[, .(rp, hs, tp)]
  return(res)
}
