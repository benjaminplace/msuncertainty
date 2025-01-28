#' Base mzML analysis Functions ----------------------------------------------------------



#' Get the intensity of a specific m/z value for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#' @param minmass numeric minimum m/z value
#' @param maxmass numeric maximum m/z value
#'
#' @export
getionint <- function(mzml, i, minmass, maxmass) {
  sum(mzml$mzML$run$spectrumList[[i]]$intensities[which(mzml$mzML$run$spectrumList[[i]]$masses >= minmass & mzml$mzML$run$spectrumList[[i]]$masses <= maxmass)])
}

#' Get the total ion intensity for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getTIC <- function(mzml, i) {
  as.numeric(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]) == "cvParam"), function(x) {if(!"total ion current" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]][[x]][which(names(mzml$mzML$run$spectrumList[[i]][[x]]) == "value")]})))
}

#' Get the base ion intensity for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getBIC <- function(mzml, i) {
  as.numeric(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]) == "cvParam"), function(x) {if(!"base peak intensity" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]][[x]][which(names(mzml$mzML$run$spectrumList[[i]][[x]]) == "value")]})))
}

#' Get the base ion m/z value for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getbaseion <- function(mzml, i) {
  as.numeric(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]) == "cvParam"), function(x) {if(!"base peak m/z" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]][[x]][which(names(mzml$mzML$run$spectrumList[[i]][[x]]) == "value")]})))
}

#' Get the time for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
gettime <- function(mzml, i) {
  as.numeric(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]$scanList$scan) == "cvParam"), function(x) {if(!"scan start time" %in% mzml$mzML$run$spectrumList[[i]]$scanList$scan[[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]]$scanList$scan[[x]][which(names(mzml$mzML$run$spectrumList[[i]]$scanList$scan[[x]]) == "value")]})))
}

#' Get the MS level for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getmslevel <- function(mzml, i) {
  as.integer(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]) == "cvParam"), function(x) {if(!"ms level" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]][[x]][which(names(mzml$mzML$run$spectrumList[[i]][[x]]) == "value")]})))
}

#' Get the precursor ion m/z value for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getprecursor <- function(mzml, i) {
  if (!"precursorList" %in% names(mzml$mzML$run$spectrumList[[i]])) {return(0)}
  as.numeric(do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]$precursorList$precursor$selectedIonList$selectedIon) == "cvParam"), function(x) {if(!"selected ion m/z" %in% mzml$mzML$run$spectrumList[[i]]$precursorList$precursor$selectedIonList$selectedIon[[x]]) {return(NULL)}; mzml$mzML$run$spectrumList[[i]]$precursorList$precursor$selectedIonList$selectedIon[[x]][which(names(mzml$mzML$run$spectrumList[[i]]$precursorList$precursor$selectedIonList$selectedIon[[x]]) == "value")]})))
}

#' Get the ionization polarity for a specific scan
#' Internal function
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param i integer scan number
#'
#' @export
getcharge <- function(mzml, i) {
  do.call(c, lapply(which(names(mzml$mzML$run$spectrumList[[i]]) == "cvParam"), function(x) {o <- NULL; if("positive scan" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {o <- 1}; if("negative scan" %in% mzml$mzML$run$spectrumList[[i]][[x]]) {o <- -1}; o}))
}

#' Generate a total ion chromatogram object
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param MS1 logical to include only MS1 (TRUE) or all MS levels (FALSE)
#' @param charge character specify the ionization polarity to include ("positive", "negative", or "ALL")
#'
#' @return named list object containing time (`time)`) and total ion intensity (`intensity`)
#' @export
#'
TIC <- function(mzml, MS1 = TRUE, charge = "ALL") {
  # creates a total ion chromatogram object with atomic properties of $intensity and $time (used for plot.chrom)
  scans <- which(names(mzml$mzML$run$spectrumList) == "spectrum")
  if (MS1 == TRUE & charge == "ALL") {
    precursors <- sapply(scans, getprecursor, mzml = mzml)
    scans <- scans[which(precursors == 0)]
  }
  if (MS1 == TRUE & charge != "ALL") {
    precursors <- sapply(scans, getprecursor, mzml = mzml)
    charges <- sapply(scans, getcharge, mzml = mzml)
    scans <- scans[which(precursors == 0 & charges == charge)]
  }
  if (MS1 == FALSE & charge != "ALL") {
    charges <- sapply(scans, getcharge, mzml = mzml)
    scans <- scans[which(charges == charge)]
  }
  intensities <- sapply(scans, getTIC, mzml = mzml)
  times <- sapply(scans, gettime, mzml = mzml)
  results <- c()
  results$intensity <- intensities
  results$time <- times
  results
}

#' Generate a base ion chromatogram object
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param MS1 logical to include only MS1 (TRUE) or all MS levels (FALSE)
#' @param charge character specify the ionization polarity to include ("positive", "negative", or "ALL")
#'
#' @return named list object containing time (`time)`) and base ion intensity (`intensity`)
#' @export
#'
BIC <- function(mzml, MS1 = TRUE, charge = "ALL") {
  # creates a base ion pair chromatogram object with atomic properties of $intensity and $time (used for plot.chrom)
  scans <- which(names(mzml$mzML$run$spectrumList) == "spectrum")
  if (MS1 == TRUE & charge == "ALL") {
    precursors <- sapply(scans, getprecursor, mzml = mzml)
    scans <- scans[which(precursors == 0)]
  }
  if (MS1 == TRUE & charge != "ALL") {
    precursors <- sapply(scans, getprecursor, mzml = mzml)
    charges <- sapply(scans, getcharge, mzml = mzml)
    scans <- scans[which(precursors == 0 & charges == charge)]
  }
  if (MS1 == FALSE & charge != "ALL") {
    charges <- sapply(scans, getcharge, mzml = mzml)
    scans <- scans[which(charges == charge)]
  }
  intensities <- sapply(scans, getBIC, mzml = mzml)
  times <- sapply(scans, gettime, mzml = mzml)
  results <- c()
  results$intensity <- intensities
  results$time <- times
  results
}

#' Generate an extracted ion chromatogram object for a specific m/z value
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param mass numeric m/z value for generating the extracted ion chromatogram
#' @param error numeric relative error (in ppm)
#' @param minerror numeric minimum absolute error for extracted ion
#' @param MS1 logical to include only MS1 (TRUE) or all MS levels (FALSE)
#' @param charge character specify the ionization polarity to include ("positive", "negative", or "ALL")
#'
#' @return named list object containing time (`time)`) and extracted ion intensity (`intensity`)
#' @export
#'
EIC <- function(mzml, mass, error, minerror = 0.001, MS1 = TRUE, charge = "ALL") {
  # creates a extracted ion chromatogram object (with error in ppm) with atomic properties of $intensity and $time (used for plot.chrom)
  mass.diff <- max(error*mass/(10^6), minerror)
  minmass <- mass - mass.diff
  maxmass <- mass + mass.diff
  scans <- which(names(mzml$mzML$run$spectrumList) == "spectrum")
  precursors <- sapply(scans, getprecursor, mzml = mzml)
  charges <- sapply(scans, getcharge, mzml = mzml)
  if (MS1 == TRUE & charge == "ALL") {scans <- scans[which(precursors == 0)]}
  if (MS1 == TRUE & charge != "ALL") {scans <- scans[which(precursors == 0 & charges == charge)]}
  if (MS1 == FALSE & charge != "ALL") {scans <- scans[which(charges == charge)]}
  intensities <- sapply(scans, getionint, mzml = mzml, minmass = minmass, maxmass = maxmass)
  times <- sapply(scans, gettime, mzml = mzml)
  results <- c()
  results$intensity <- intensities
  results$time <- times
  results
}


#' Get local minimum (peak start and end) in chromatographic peak
#' Internal Function
#'
#' @param x numeric vector representing a intensities of a chromatographic peak
#' @param peak integer index of local maximum
#' @param width integer number of points across peak to use to calculate slope
#' @param slope numeric slope limit that determines peak start and end
#'
#' @return numeric vector containing the index of the start and end of the peak
#' @export
#'
peak_desc <- function(x, peak, width = 1, slope = 0.00044) {
  # Internal function used to determine the peak width in a single dimension
  rx <- peak
  lx <- peak
  rightwidth <- peak
  leftwidth <- peak
  rightslope = -slope - 1
  leftslope = -slope - 1
  while (rightslope <= -slope & rx <= length(x)) {
    rightslope <- slopeest(x, rx:min(c(rx + width, length(x))))
    if (is.na(rightslope)) {rightslope = -slope}
    rightwidth <- rx
    rx <- rx + 1
  }
  while (leftslope <= -slope & lx > 0) {
    leftslope <- slopeest(x, lx:max(c(rx-width, 1)))
    if (is.na(leftslope)) {leftslope = -slope}
    leftwidth <- lx
    lx <- lx - 1
  }
  c(leftwidth, rightwidth)
}

#' Calculate slope across points
#'
#' @param x numeric vector containing two or more values
#' @param points integer vector containing index of points of `x` to calculate slope
#'
#' @return numeric value of estimated slope.
#' @export
#'
slopeest <- function(x, points) {
  xs <- 1:length(points)
  ys <- x[points]
  sum((xs-mean(xs))*(ys-mean(ys)))/sum((xs-mean(xs))^2)
}

#' Determine peak width at specific level of chromatographic peak
#' Internal Function
#'
#' @param chrom numeric vector representing a intensities of a chromatographic peak
#' @param peak integer index of local maximum
#' @param width integer number of points across peak to use to calculate slope
#' @param slope numeric slope limit that determines peak start and end
#' @param level numeric indicating the relative peak height to determine the width (0.1 means the width is determined at 10 percent of the peak height)
#'
#' @return numeric vector containing the index of the start and end of the peak at the specified peak height
#' @export
#'
peak_width <- function(chrom, peak, width = 100, slope = 0.0001, level = 0.1) {
  widths.init <- peak_desc(chrom, peak, width = width, slope = slope)
  leftmin <- (chrom[peak] - chrom[widths.init[1]])*level
  rightmin <- (chrom[peak] - chrom[widths.init[1]])*level
  leftwidth <- which(chrom[widths.init[1]:peak] > leftmin)
  leftwidth <- leftwidth[1] + widths.init[1] - 1
  rightwidth <- which(chrom[widths.init[2]:peak] > leftmin)
  rightwidth <- widths.init[2] - rightwidth[1] + 1
  widths <- c(leftwidth, rightwidth)
  widths
}


#' Generate mass spectrum of one or more scans
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param scans integer vector of the one or more scan indices to combine into a single average spectrum
#' @param mz.round integer value number of decimals to round the m/z values
#'
#' @importFrom stats ave
#'
#' @return data.frame containing mass spectrum `mz` and `int` values
#' @export
#'
extract_ms <- function(mzml, scans, mz.round = 4) {
  # creates mass spectrum object for list of scan numbers (not times)
  full.ms <- NULL
  for (j in scans) {
    x.ms <- cbind(mzml$mzML$run$spectrumList[[j]]$masses, mzml$mzML$run$spectrumList[[j]]$intensities)
    x.ms[,1] <- round(x.ms[,1], digits = mz.round)
    x.ms <- cbind(x.ms[,1], ave(x.ms[,2], x.ms[,1], FUN = sum))
    x.ms <- x.ms[!duplicated(x.ms),]
    x.ms
    if (j == scans[1]) {
      full.ms <- x.ms
    }
    else {
      int.ms <- merge(full.ms, x.ms, by = 1, all = TRUE)
      full.ms <- cbind(int.ms[,1], rowSums(int.ms[,2:3], na.rm = TRUE))
    }
  }
  full.ms[,1] <- round(full.ms[,1], digits = mz.round)
  full.ms <- cbind(full.ms[,1], ave(full.ms[,2], full.ms[,1], FUN = sum))
  full.ms <- full.ms[!duplicated(full.ms),]
  colnames(full.ms) <- c("mz", "int")
  full.ms
}

#' Save mass spectrum to .msp file
#'
#' @param ms data.frame mass spectrum generated from `extract_ms` function
#' @param file character string of file name to write mass spectrum to (include extension)
#' @param precursor numeric value of the precursor of the saved mass spectrum (not required)
#' @param name character string of the title of the saved mass spectrum
#'
#' @importFrom utils write.table
#'
#' @return exported text file containing mass spectrum
#' @export
#'
write_spectrum <- function(ms, file, precursor = "", name = "R-exported Spectrum") {
  # saves mass spectrum as a NIST MS Search-readable file with precursor ion designated

  peaks <- nrow(ms)
  names <- paste("Name: ", name)
  precursor <- paste("PrecursorMZ: ", precursor)
  peaks <- paste("Num Peaks: ", peaks)
  header <- rbind(names, precursor, peaks)
  write.table(header, file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (!is.na(ms[1])) {
    write.table(ms, file, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  if (is.na(ms[1])) {
    write.table(cbind(0,0), file, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  write("\n\n", file, append = TRUE)
}

#just in case I forget in other functions
peak.desc <- function(...) {
  warning("This function is deprecated, use `peak_desc` for future functions.")
  peak_desc(...)
}

peak.width <- function(...) {
  warning("This function is deprecated, use `peak_width` for future functions.")
  peak_width(...)
}

extract.ms <- function(...) {
  warning("This function is deprecated, use `extract_ms` for future functions.")
  extract_ms(...)
}



