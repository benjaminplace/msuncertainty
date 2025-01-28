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
    which(names(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList) == "binaryDataArray") %>%
      sapply(., function(y) "m/z array" %in% unlist(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[y]])) %>%
      which() %>%
      mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[.]] %>%
      .[["binary"]] %>%
      unzipms(compression = compression)
  })

  intensities <- lapply(scans, function(x) {
    which(names(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList) == "binaryDataArray") %>%
      sapply(., function(y) "intensity array" %in% unlist(mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[y]])) %>%
      which() %>%
      mzml$mzML$run$spectrumList[[x]]$binaryDataArrayList[[.]] %>%
      .[["binary"]] %>%
      unzipms(compression = compression)
  })

  for (x in 1:length(scans)) {
    mzml$mzML$run$spectrumList[[scans[x]]]$masses <- masses[[x]]
    mzml$mzML$run$spectrumList[[scans[x]]]$intensities <- intensities[[x]]
  }

  mzml
}

#' Unzip mzML mass spectra
#'
#' Internal function used by mzMLtoR function to unzip mass spectral binary character strings.
#'
#' @param binary character string representing
#' @param compression
#'
#' @return double vector containing unzipped MS data. If binary is empty, returns NULL.
#' @importFrom base64enc base64decode
#'
#' @export
#'
unzipms <- function(binary, compression = "none") {
  if (is.null(binary) | length(binary) == 0) {
    out <- NULL
  }
  if (!is.null(binary) & length(binary) > 0) {
    out <- binary %>%
      base64decode() %>%
      memDecompress(type = compression) %>%
      readBin(what = "double", n = length(.)%/%8, size = 8)
  }
  out
}
