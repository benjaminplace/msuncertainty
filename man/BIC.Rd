% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_ms.R
\name{BIC}
\alias{BIC}
\title{Generate a base ion chromatogram object}
\usage{
BIC(mzml, MS1 = TRUE, charge = "ALL")
}
\arguments{
\item{mzml}{list object parsed from mzML file using `mzMLtoR` function}

\item{MS1}{logical to include only MS1 (TRUE) or all MS levels (FALSE)}

\item{charge}{character specify the ionization polarity to include ("positive", "negative", or "ALL")}
}
\value{
named list object containing time (`time)`) and base ion intensity (`intensity`)
}
\description{
Generate a base ion chromatogram object
}
