% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_ms.R
\name{write_spectrum}
\alias{write_spectrum}
\title{Save mass spectrum to .msp file}
\usage{
write_spectrum(ms, file, precursor = "", name = "R-exported Spectrum")
}
\arguments{
\item{ms}{data.frame mass spectrum generated from `extract_ms` function}

\item{file}{character string of file name to write mass spectrum to (include extension)}

\item{precursor}{numeric value of the precursor of the saved mass spectrum (not required)}

\item{name}{character string of the title of the saved mass spectrum}
}
\value{
exported text file containing mass spectrum
}
\description{
Save mass spectrum to .msp file
}
