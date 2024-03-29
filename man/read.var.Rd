% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catch.R
\name{read.var}
\alias{read.var}
\title{Reading Discard and Catch by number and box}
\usage{
read.var(FG, nc.data, is.C = NULL, grp, by.box = FALSE, time)
}
\arguments{
\item{FG}{Csv with the Functional groups}

\item{nc.data}{NetCDF file with the information of catch}

\item{is.C}{Default TRUE,  if the data needed is Catch (TRUE) or Discard (FALSE)}

\item{grp}{Specific group to get from the csv file}

\item{by.box}{Default FALSE. If the information is needed by box (TRUE) or not (FALSE)}

\item{time}{Time vector}
}
\value{
A list witht the information for all the functional groups
}
\description{
Reading Discard and Catch by number and box
}
\author{
Demiurgo
}
