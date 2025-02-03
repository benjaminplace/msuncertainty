#' mzML I/O Functions ----------------------------------------------------------

#' Read mzML file into R environment
#'
#' Based function to read a standard mzML file into an R object. This includes unzipping any gzip-based spectra.
#'
#' @param mzmlfile string path to .mzML file
#'
#' @return list object containing mzML data including metadata
#' @importFrom XML xmlToList
#' @importFrom base64enc base64decode
#' @importFrom magrittr %>%
#'
#' @section note: currently this format does not support the use of lockmass scans, this will be incorporated in future updates
#'
#' @export
#'
#' @examples
#' mzmlfile <- system.file('extdata', 'PFAS_STD.mzML', package = 'msuncertainty')
#' mzml <- mzMLtoR(mzmlfile)
#'
mzMLtoR <- function(mzmlfile = file.choose()) {
  mzml <- XML::xmlToList(mzmlfile)

  compression <- "none"
  if ("zlib compression" %in% unlist(mzml$mzML$run$spectrumList[[1]]$binaryDataArrayList[["binaryDataArray"]])) {
    compression = "gzip" #check for zlib compression for memDecompress
  }

  scans <- which(names(mzml$mzML$run$spectrumList) == "spectrum")

  masses <- lapply(scans, function(x) {
      x.1 <- which(names(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList) == "binaryDataArray")
      x.2 <- sapply(x.1, function(y) "m/z array" %in% unlist(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[y]]))
      x.3 <- which(x.2)
      x.4 <- mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[x.3]]
      unzipms(x.4[["binary"]], compression = compression, size = switch(get_float(x.4), float32bit = 4, float64bit = 8))
  })

  intensities <- lapply(scans, function(x) {
    x.1 <- which(names(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList) == "binaryDataArray")
    x.2 <- sapply(x.1, function(y) "intensity array" %in% unlist(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[y]]))
    x.3 <- which(x.2)
    x.4 <- mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[x.3]]
    unzipms(x.4[["binary"]], compression = compression, size = switch(get_float(x.4), float32bit = 4, float64bit = 8))
  })

  for (x in 1:length(scans)) {
    mzml$mzML$run$spectrumList[[scans[x]]]$masses <- masses[[x]]
    mzml$mzML$run$spectrumList[[scans[x]]]$intensities <- intensities[[x]]
  }

  mzml
}


#' Get encoding float size from single binary data array
#'
#' @param binaryDataArrayList list from binary data array within a single mzml scan
#'
#' @returns character string of `float32bit` or `float64bit`
#' @export
#'
get_float <- function(binaryDataArrayList) {
  names <- sapply(which(names(binaryDataArrayList) == "cvParam"), function(y) binaryDataArrayList[[y]]["name"])
  if ("32-bit float" %in% names) {return("float32bit")}
  if ("64-bit float" %in% names) {return("float64bit")}
}

#' Unzip mzML mass spectra
#'
#' Internal function used by mzMLtoR function to unzip mass spectral binary character strings.
#'
#' @param binary character string representing
#' @param compression character listing the type of compression used ("none" or "gzip")
#' @param size integer float size (4 or 8)
#'
#' @return double vector containing unzipped MS data. If binary is empty, returns NULL.
#' @importFrom base64enc base64decode
#'
#' @export
#'
unzipms <- function(binary, compression = "none", size = 8) {
  if (is.null(binary) | length(binary) == 0) {
    out <- NULL
  }
  if (!is.null(binary) & length(binary) > 0) {
    out <- binary %>%
      base64decode() %>%
      memDecompress(type = compression) %>%
      readBin(what = "double", n = length(.)%/%size, size = size)
  }
  out
}
