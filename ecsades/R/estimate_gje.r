# GJE ---------------------------------------------------------------------

estimate_gje = function(
  jdistr, output_rp,
  n_point = 100,
  ref_tp = rep(0, length(output_rp)),
  ref_hs = rep(0, length(output_rp))){
  
  ## Generate sample data
  sample_data = .sample_jdistr(jdistr = jdistr, sim_year = max(output_rp)*.rp_multiplier)  
  sample_data[, q:=NA_character_]
  
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
    target_prob = 1/(this_rp*jdistr$npy)
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


