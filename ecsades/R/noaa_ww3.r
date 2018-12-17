#' Sample wave data by NOAA
#' 
#' @description 
#' This sample dataset contains wave peak data (large Hs and/or Tp) as extracted from
#' the 3-hourly (before March 1st 2010) and 1-hourly (after March 1st 2010) wave time series data 
#' provided by NOAA's WaveWatch III project for location 30.5N, 240E.
#'
#' @docType data
#' 
#' @name noaa_ww3
#' 
#' @format An object of class \code{data.table} with columns \code{time}, \code{hs} and \code{tp};
#' see \code{\link{data.table}}.
#'
#' @references The WAVEWATCH III Development Group., 2016:
#' User manual and system documentation of WAVEWATCH III version 5.16.
#' NOAA / NWS / NCEP / MMAB Technical Note 329, 326 pp. + Appendices.
#'
#' @source \href{http://polar.ncep.noaa.gov/waves/wavewatch/}{NOAA WAVEWATCH III Model}
#'
#' @examples
#' data(noaa_ww3)


NULL