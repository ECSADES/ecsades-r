## 2h
plot_ec = function(ec, raw_data=NULL, hs_x=TRUE, save_to_file = NULL, ...){
  
  # Base layer
  p0 = ggplot()
  
  # Point layer
  if(!is.null(raw_data)){
    if(hs_x){
      p_point = geom_point(aes(x=hs, y=tp), data=raw_data)
    }else{
      p_point = geom_point(aes(x=tp, y=hs), data=raw_data)
    }
  }else{
    p_point = NULL
  }
  
  # Contour layer
  if(hs_x){
    p_contour = geom_path(aes(x=hs, y=tp, group=rp, colour=factor(rp)), data=ec)
  }else{
    p_contour = geom_path(aes(x=tp, y=hs, group=rp, colour=factor(rp)), data=ec)
  }
  
  # Output
  p_out = p0 + p_point + p_contour+
    scale_colour_discrete("Return period (y)")
  if(!is.null(save_to_file)){
    ggsave(save_to_file, p_out, ...)
  }
  return(p_out)
}
