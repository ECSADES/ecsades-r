#' Sample wave data by NOAA
#' 
#' @description 
#' This sample dataset contains the wave height and wave period time series provided by 
#' NOAA's WaveWatch III project for location 30.5N, 240E. Two variations of the dataset are
#' included:
#' 
#' 1) the hourly time series (\code{ww3_ts}), interpolated down from 3-hourly data before
#' March 2010, and
#' 
#' 2) the extracted wave peak data (\code{ww3_pk}), where the peaks are extracted base
#' on the peak Hs and/or peak Tp.
#'
#' @docType data
#' 
#' @name ww3_data
#' @aliases ww3_pk ww3_ts
#' 
#' @format An object of class \code{data.table} with columns \code{time}, \code{hs} and \code{tp};
#' see \code{\link{data.table}}.
#'
#' @references The WAVEWATCH III Development Group., 2016:
#' User manual and system documentation of WAVEWATCH III version 5.16.
#' NOAA / NWS / NCEP / MMAB Technical Note 329, 326 pp. + Appendices.
#'
#' @source \href{https://polar.ncep.noaa.gov/waves/}{NOAA WAVEWATCH III Model}
#'
#' @examples
#' data(ww3_ts)
#' data(ww3_pk)

NULL