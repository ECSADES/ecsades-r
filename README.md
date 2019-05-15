# ECSADES
Environmental contours for safe design of ships and other marine structures using R.

## Introduction
This package contains functions to fit joint probability distributions to bivariate wave data (wave height and wave period), estimate environmental contours based on the fitted distribution or directly on a given sample data, output the coordinates of the vertices of the contours, and generate a plot of the contours.

The choices for the joint probability distributions include the Heffernan-Tawn model and the Weibull-log-normal distribution. The choices for the environmental contours are the direct sampling contours (Huseby et al. 2015), IFORM contours (Winterstein et al. 1993), generalised joint exceedance contours (Jonathan et al. 2014), and isodensity contours.

The statistical models and contour estimation methods included in this package are discussed in detail in the project paper Ross et al. (2018).

## Getting Started
The user will need the ```devtools``` package to install directly from a github repository.
```
install.packages("devtools")
```
The following command will install the latest versioin of the package on user's computer.
```
devtools::install_github(repo="ECSADES/ecsades-r", branch='master', subdir='ecsades', dependencies = TRUE)
```
An example of applying the functions in the package can be found in the package help file.
```
?ecsades::ecsades
```
There are two sample datasets provided in the package, which are used in examples in the help files.
```
data(ww3_pk, package = "ecsades")
data(ww3_ts, package = "ecsades")
```
Alternatively, the user could use function ```fread()``` from the ```data.table``` package to read in their own wave data.

## Author
* [Ye Liu](https://github.com/yeliuhrw) - y.liu@hrwallingford.com

## Acknowledgments
This work was part-funded by the European Union ERA-NET project entitled "Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES)". The author thanks Erik Vanem and Arne Husbey for their help on parts of the code, and the project team for useful discussions.
