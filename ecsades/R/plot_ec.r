#' Plotting the environmental contours
#'
#' @description
#' This function prints on screen a given set of environmental contours and saves the image
#' to a file with optional customisation.
#' 
#' @param ec the environmental contours in the format of a \code{data.table} with
#' \code{rp}, \code{hs} and \code{tp} as columns.  These can be generated using the contour
#' estimation functions such as \code{\link{estimate_dsc}}, \code{\link{estimate_iform}}, etc.
#' 
#' @param raw_data the (optional) raw data behind the the fitted contours. Providing this data
#' will add the corresponding points to the final plot. The default value is \code{NULL}.
#' 
#' @param hs_x whether or not the x-axis is the wave height \code{hs}. Setting \code{hs_x = TRUE}
#' means the plot is \code{(hs, tp)}; otherwise it would be \code{(tp, hs)}. The default
#' value is \code{TRUE}.
#' 
#' @param save_to_file an (optional) file name for saving the image to. The file type can be 
#' anything that is supported by the \code{\link{ggsave}} function. The default value is
#' \code{NULL} which means the image is not saved.
#' 
#' @param ... additional (optional) arguments passed onto function \code{\link{ggsave}}.
#' 
#' @return A \code{\link{ggplot}} object containing the input environmental contours
#' and other optional specified items.
#'
#' @examples
#' # Load data
#' data(noaa_ww3)
#' 
#' # Fit the Heffernan-Tawn model to the data
#' ht = fit_ht(
#'   data = noaa_ww3,
#'   npy = nrow(noaa_ww3)/10,
#'   margin_thresh_count = 100,
#'   dep_thresh_count = 100)
#' 
#' # Estimate the DSC contours using the Hs return levels as 
#' dsc = estimate_dsc(jdistr = ht,  output_rp = c(1,10,100))
#'   
#' # Plot output
#' plot_ec(
#'   ec = dsc, raw_data = noaa_ww3, hs_x = FALSE, save_to_file = "dsc_ht.png",
#'   width = 6, height = 6, units = "in")
#' 
#' @export
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
