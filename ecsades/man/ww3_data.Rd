% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ww3_data.r
\docType{data}
\name{ww3_data}
\alias{ww3_data}
\alias{ww3_pk}
\alias{ww3_ts}
\title{Sample wave data by NOAA}
\format{An object of class \code{data.table} with columns \code{time}, \code{hs} and \code{tp};
see \code{\link{data.table}}.}
\source{
\href{https://polar.ncep.noaa.gov/waves/}{NOAA WAVEWATCH III Model}
}
\description{
This sample dataset contains the wave height and wave period time series provided by 
NOAA's WaveWatch III project for location 30.5N, 240E. Two variations of the dataset are
included:

1) the hourly time series (\code{ww3_ts}), interpolated down from 3-hourly data before
March 2010, and

2) the extracted wave peak data (\code{ww3_pk}), where the peaks are extracted base
on the peak Hs and/or peak Tp.
}
\examples{
data(ww3_ts)
data(ww3_pk)
}
\references{
The WAVEWATCH III Development Group., 2016:
User manual and system documentation of WAVEWATCH III version 5.16.
NOAA / NWS / NCEP / MMAB Technical Note 329, 326 pp. + Appendices.
}
