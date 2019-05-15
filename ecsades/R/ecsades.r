#' @details
#' This package contains functions to fit joint probability distributions to bivariate wave data (wave height and 
#' wave period), estimate environmental contours based on the fitted distribution or directly on a given sample 
#' data, output the coordinates of the vertices of the contours, and generate a plot of the contours.
#' 
#' The choices for the joint probability distributions include the Heffernan-Tawn model \code{\link{fit_ht}} and
#' the Weibull-log-normal distribution \code{\link{fit_wln}}. The choices for the environmental contours are the
#' direct sampling contours \code{\link{estimate_dsc}}, IFORM
#' contours \code{\link{estimate_iform}}, generalised joint exceedance contours \code{\link{estimate_gje}}, and the
#' isodensity contours \code{\link{estimate_iso}}.
#' 
#' This package has been developed as part of the
#' \href{https://www.dnvgl.com/technology-innovation/sri/maritime-transport/ecsades-project.html}{ECSADES} project
#' and is free to use. The
#' statistical models and contour estimation methods included in this package are discussed in detail in the project
#' paper Ross et al. (2018).
#' 
#' @examples
#' # Load sample data from NOAA's WaveWatch III project
#' data(ww3_pk)
#' 
#' # Fit the Heffernan-Tawn model to the data
#' ht = fit_ht(
#'   data = ww3_pk,
#'   npy = nrow(ww3_pk)/10,
#'   margin_thresh_count = 100,
#'   dep_thresh_count = 100)
#' 
#' # Estimate the DSC contours using the Hs return levels as 
#' dsc = estimate_dsc(jdistr = ht,  output_rp = c(1,10,100))
#'   
#' # Plot output and save to a file
#' plot_ec(
#'   ec = dsc, raw_data = ww3_pk, hs_x = FALSE, save_to_file = "dsc_ht.png",
#'   width = 6, height = 6, units = "in")
#' 
#' @references 
#' Emma Ross, Ole Christian Astrup, Elzbieta BitnerGregersen, Nigel Bunn, Graham Feld, Ben Gouldby, Arne Huseby,
#' Ye Liu, David Randell, Erik Vanem, Philip Jonathan, 2018. On environmental contours for marine and coastal design.
#' (submitted)
"_PACKAGE"
.limit_inf = 1e7
.limit_zero = 1e-7
.seed_sorting = 170899
.seed_sampling = 121110
.rp_multiplier = 50
.knn = 100
