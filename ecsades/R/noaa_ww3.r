#' Sample wave data by NOAA
#'
#' This sample dataset contains wave data (Hs and Tp) provided by the WaveWatch III project
#' by NOAA for location 30.5N, 240E.  It contains 3-hourly data from 1 January 2006 to 
#' 31 December 2016.
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
#' # Loading the 3-hourly time series data
#' data(noaa_ww3_ts)
#' 
#' # Loading the declustered event peak data
#' data(noaa_ww3_peak)

NULL