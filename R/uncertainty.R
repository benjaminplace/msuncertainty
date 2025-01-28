#' Spectral Uncertainty Functions ----------------------------------------------------------

#' Get scan numbers for a specific m/z value
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param mass numeric m/z value of interest
#' @param peak_start numeric time (in minutes) of the start of the chromatographic peak
#' @param peak_end  numeric time (in minutes) of the end of the chromatographic peak
#' @param isowidth numeric instrument isolation width (in Da) of the precursor ion, see below note.
#'
#' @section note: for DDA, isowidth is the Q1 isolation width, for SWATH/SONAR the isowidth is half the scan window width, and for DIA/AIF isowidth should be `NA`
#'
#' @return numeric vector containing the scan numbers for the MS1 and relevant MS2 scans for the specific precursor ion
#' @export
#'
get_peak_scans <- function(mzml, mass, peak_start, peak_end, isowidth = 0.7) {
  scans <- which(names(mzml$mzML$run$spectrumList) == "spectrum")
  times <- sapply(scans, gettime, mzml = mzml)
  startscan <- which.min(abs(times - peak_start))
  endscan <- which.min(abs(times - peak_end))
  subscans <- startscan:endscan
  mslevels <- sapply(subscans, getmslevel, mzml = mzml)
  ms1scans <- subscans[which(mslevels == 1)]
  ms2scans <- subscans[which(mslevels == 2)]
  if (!is.na(isowidth)) {
    precursors <- sapply(ms2scans, getprecursor, mzml = mzml)
    ms2scans <- ms2scans[which(precursors - isowidth <= mass & precursors + isowidth >= mass)]
  }
  if (length(ms2scans) < 1) {warning("There is no MS2 data available based on the parameters.")}
  totalscans <- sort(c(ms1scans, ms2scans))
  totalscans
}

#' Create peak list from mzML data
#'
#' The function extracts the relevant information and sorts it into nested lists for
#' use in the uncertainty functions
#'
#' @param mzml list object parsed from mzML file using `mzMLtoR` function
#' @param scans numeric vector containing scans to extract into peak list relevant, preferably produced by `get_peak_scans` function
#'
#' @section note: This function does not check to see if MS2 data is relevant to a specific peak, the `scans` parameter should be calculated through the `get_peak_scans` function
#'
#' @return nested list of all data
#' @export
#'
create_peak_list_from_mzml <- function(mzml, scans) {
  lapply(scans, function(x) list(
    masses =  extract_ms(mzml, x)[,1],
    intensities = extract_ms(mzml, x)[,2],
    totalion = getTIC(mzml, x),
    baseion = getbaseion(mzml, x),
    base_int = getBIC(mzml, x),
    time = gettime(mzml, x),
    mslevel = getmslevel(mzml,x)))
}

#' Create peak list from SQL ms_data table
#'
#' The function extracts the relevant information and sorts it into nested lists for
#' use in the uncertainty functions
#'
#' @param ms_data extraction of the ms_data from the SQL table for a specified peak
#'
#' @return nested list of all data
#' @export
#'
create_peak_list <- function(ms_data) {
  lapply(1:nrow(ms_data), function(x) list(
    masses =  as.numeric(unlist(strsplit(ms_data$measured_mz[x]," "))),
    intensities = as.numeric(unlist(strsplit(ms_data$measured_intensity[x], " "))),
    totalion = sum(as.numeric(unlist(strsplit(ms_data$measured_intensity[x], " ")))),
    baseion = ms_data$baseion[x],
    base_int = ms_data$base_int[x],
    time = ms_data$scantime[x],
    mslevel = ms_data$ms_n[x]))
}


#' Create peak table for MS2 data
#'
#' Takes a nested peak list and creates a peak table for easier determination of
#' uncertainty of the measurement for MS2 data.
#'
#' @param peak list result of the `create_peak_list` function
#' @param mass the exact mass of the compound of interest
#' @param masserror the mass accuracy (in ppm) of the instrument data
#' @param minerror the minimum mass error (in Da) of the instrument data
#' @param int0 the default setting for intensity values for missing m/z values
#'
#' @return nested list of dataframes containing all MS2 data for the peak
#' @export
#'
create_peak_table_ms2 <- function(peak, mass, masserror = 5, minerror = 0.002, int0 = NA) {
  #creates MS2 data tables from a peak object including the mass and intensity (separate tables) into a peak table object
  scans <- 1:length(peak)
  mslevels <- sapply(scans, function(x) peak[[x]]$mslevel)
  #need to include error function to error if no MS2 data
  peak.EIC <- getEIC_from_peaklist(peak, mass, masserror, minerror, mslevel = "MS1")
  ms1scans <- scans[which(mslevels == 1)]
  ms2scans <- scans[which(mslevels == 2)]
  if (length(ms2scans) < 1) {stop("There is no MS2 data in the peak list")}
  ms1scans_peak <- ms2scans - 1
  ms1scans_ind <- which(ms1scans %in% ms1scans_peak)
  peakindex <- which(ms2scans == ms1scans_peak[which.max(peak.EIC$intensity[ms1scans_ind])]+1)
  mslist <- lapply(ms2scans, function(x) cbind(peak[[x]]$masses, peak[[x]]$intensities))
  merged <- mergems(mslist, peakindex, masserror = masserror, minerror = minerror)
  peaktable_int <- matrix(NA, nrow = length(merged), ncol = length(ms2scans))
  peaktable_mass <- matrix(int0, nrow = length(merged), ncol = length(ms2scans))
  for (i in 1:length(merged)) {
    peaktable_int[i,merged[[i]]$scans] <- merged[[i]]$ints
    peaktable_mass[i,merged[[i]]$scans] <- merged[[i]]$masses
  }
  out <- list(peaktable_int = peaktable_int, peaktable_mass = peaktable_mass, EIC = peak.EIC, ms1scans = ms1scans, ms2scans = ms2scans)
  attr(out, "mslevel") <- 2
  out
}

#' Create peak table for MS1 data
#'
#' Takes a nested peak list and creates a peak table for easier determination of
#' uncertainty of the measurement for MS1 data.
#'
#' @param peak list result of the `create_peak_list` function
#' @param mass the exact mass of the compound of interest
#' @param masserror the mass accuracy (in ppm) of the instrument data
#' @param minerror the minimum mass error (in Da) of the instrument data
#' @param int0 the default setting for intensity values for missing m/z values
#'
#' @return nested list of dataframes containing all MS2 data for the peak
#' @export
#'
create_peak_table_ms1 <- function(peak, mass, masserror = 5, minerror = 0.002, int0 = NA) {
  #creates MS1 data tables from a peak object including the mass and intensity (separate tables) into a peak table object
  scans <- 1:length(peak)
  mslevels <- sapply(scans, function(x) peak[[x]]$mslevel)
  peak.EIC <- getEIC_from_peaklist(peak, mass, masserror, minerror, mslevel = "MS1")
  ms1scans <- scans[which(mslevels == 1)]
  ms2scans <- scans[which(mslevels == 2)]
  peakindex <- which.max(peak.EIC$intensity)
  mslist <- lapply(ms1scans, function(x) cbind(peak[[x]]$masses, peak[[x]]$intensities))
  merged <- mergems(mslist, peakindex, masserror = masserror, minerror = minerror)
  peaktable_int <- matrix(NA, nrow = length(merged), ncol = length(ms1scans))
  peaktable_mass <- matrix(int0, nrow = length(merged), ncol = length(ms1scans))
  for (i in 1:length(merged)) {
    peaktable_int[i,merged[[i]]$scans] <- merged[[i]]$ints
    peaktable_mass[i,merged[[i]]$scans] <- merged[[i]]$masses
  }
  out <- list(peaktable_int = peaktable_int, peaktable_mass = peaktable_mass, EIC = peak.EIC, ms1scans = ms1scans, ms2scans = ms2scans)
  attr(out, "mslevel") <- 1
  out
}

#' Get Extracted Ion Chromatogram (EIC) from a peak list object
#' Internal function
#'
#' @param peak list object result of the `create_peak_list` function
#' @param mass numeric m/z value of targeted mass
#' @param masserror numeric relative mass error (in ppm)
#' @param minerror numeric minimum absolute mass error (in Da)
#' @param mslevel character designating an EIC from MS1 ("MS1") or MS2 ("MS2") data
#'
#' @return data.frame object containing the EIC for a specified m/z value
#' @export
#'
getEIC_from_peaklist <- function(peak, mass, masserror = 5, minerror = 0.002, mslevel = "MS1") {
  mass.diff <- max(masserror*mass/(10^6), minerror)
  minmass <- mass - mass.diff
  maxmass <- mass + mass.diff
  time <- sapply(peak, function(x) x$time)
  int <- sapply(peak, function(x) sum(x$intensities[which(x$masses >= minmass & x$masses <= maxmass)]))
  if (mslevel == "MS1") {
    mslevels <- sapply(peak, function(x) x$mslevel)
    int <- int[which(mslevels == 1)]
    time <- time[which(mslevels == 1)]
  }
  if (mslevel == "MS2") {
    mslevels <- sapply(peak, function(x) x$mslevel)
    int <- int[which(mslevels == 2)]
    time <- time[which(mslevels == 2)]
  }
  list(time = time, intensity = int)
}


#' Generate a merged mass spectrum comprised of multiple individual mass spectra
#' Internal function
#'
#' @param mslist list object containing multiple data.frames with individual mass spectra
#' @param peakindex integer representing the initial mass spectrum in the `mslist` to start merging process
#' @param masserror numeric relative mass error (in ppm)
#' @param minerror numeric minimum absolute mass error (in Da)
#'
#' @return list object containing merged masses, intensities, and scans
#' @export
#'
mergems <- function(mslist, peakindex, masserror = 5, minerror = 0.001) {
  ms <- mslist[[peakindex]]
  if (nrow(ms) == 1) {
    list(masses = ms[1,1], ints = ms[1,2], scans = peakindex)
  }
  mergedms <- lapply(1:nrow(ms), function(x) list(masses = ms[x,1], ints = ms[x,2], scans = peakindex))
  if (length(mslist) == 1) {
    return(mergedms)
  }
  ind <- seq_len(length(mslist))
  inds <- order(abs(ind - peakindex))
  inds <- inds[-1]
  for (i in inds) {
    if (nrow(mslist[[i]]) > 0) {
      newms <- lapply(mergedms, function(x) mergemass(x$masses, x$ints, x$scans, mslist[[i]], ind = i, masserror = masserror, minerror = minerror))
      leftin <- do.call(c, lapply(mergedms, function(x) getmergedind(x$masses, mslist[[i]], masserror = masserror, minerror = minerror)))
      leftout <- setdiff(1:nrow(mslist[[i]]), leftin)
      if (length(leftout) > 0) {
        extrams <- lapply(leftout, function(x) list(masses = mslist[[i]][x,1], ints = mslist[[i]][x,2], scans = i))
        mergedms <- append(newms, extrams)
      }
    }
  }
  mergedms
}


#' Merges a set of m/z values of a mass spectra based on a list of m/z values
#' Internal function
#'
#' @param masses numeric vector containing m/z values to be merged together
#' @param ints numeric vector containing int values to be merged together
#' @param scans integer vector
#' @param addms numeric data.frame containing m/z and intensity values to be merged into the masses and ints
#' @param ind integer index number of the scan
#' @param masserror numeric relative mass error (in ppm)
#' @param minerror numeric minimum absolute mass error (in Da)
#'
#' @return list object containing all merged masses
#' @export
#'
mergemass <- function(masses, ints, scans, addms, ind, masserror, minerror) {
  minmass <- min(min(masses) - min(masses)*masserror/(10^6),min(masses) - minerror)
  maxmass <- max(max(masses) + max(masses)*masserror/(10^6),max(masses) + minerror)
  inds <- which(addms[,1] >= minmass & addms[,1] <= maxmass)
  masses <- c(masses, addms[inds,1])
  ints <- c(ints, addms[inds,2])
  scans <- c(scans, rep(ind, length(inds)))
  output <- list(masses = masses, ints = ints, scans = scans)
  output
}


#' Get indexes of m/z values to be merged together
#' Internal Function
#'
#' @param masses numeric vector containing m/z values to be merged together
#' @param addms numeric data.frame containing m/z and intensity values to be merged into the masses and ints
#' @param masserror numeric relative mass error (in ppm)
#' @param minerror numeric minimum absolute mass error (in Da)
#'
#' @return integer vector containing m/z values to be merged together
#' @export
#'
getmergedind <- function(masses, addms, masserror, minerror) {
  minmass <- min(min(masses) - min(masses)*masserror/(10^6),min(masses) - minerror)
  maxmass <- max(max(masses) + max(masses)*masserror/(10^6),max(masses) + minerror)
  inds <- which(addms[,1] >= minmass & addms[,1] <= maxmass)
  inds
}

#' Generate consensus mass spectrum
#'
#' The function calculates the uncertainty mass spectrum for a single peak table based
#' on specific settings described in https://doi.org/10.1021/jasms.0c00423
#'
#' @param peaktable result of the `create_peak_table_ms1` or  `create_peak_table_ms1` function
#' @param correl Minimum correlation coefficient between the target ions and the base ion intensity of the targeted m/z to be included in the mass spectrum
#' @param ph Minimum chromatographic peak height from which to extract MS2 data for the mass spectrum
#' @param freq minimum observational frequency of the target ions to be included in the mass spectrum
#' @param normfn the normalization function typically "mean" or "sum" for normalizing the intensity values
#' @param cormethod the correlation method used for calculating the correlation, see `cor` function for methods
#'
#' @importFrom stats sd
#' @importFrom stats cor
#'
#' @return nested list of dataframes containing all MS1 and MS2 data for the peak
#' @export
#'
get_ums <- function(peaktable, correl = NULL, ph = NULL, freq = NULL, normfn = "sum", cormethod = "pearson") {

  if (attr(peaktable, "mslevel") == 2) {eic_cor <- eic <- peaktable$EIC$intensity[sapply(peaktable$ms2scans, function(x) which.min(abs(x - peaktable$ms1scans)))]}
  if (attr(peaktable, "mslevel") == 1) {eic_cor <- eic <- peaktable$EIC$intensity}
  if (!is.null(correl)) {
    if (!is.na(correl)) {
      suppressWarnings(cors <- apply(peaktable$peaktable_int, 1, function(y) cor(eic, y, use = "complete.obs", method = cormethod)))
      ind <- which(cors >= correl)
      peaktable$peaktable_int <- data.frame(peaktable$peaktable_int[ind,], fix.empty.names = FALSE)
      peaktable$peaktable_mass <- data.frame(peaktable$peaktable_mass[ind,], fix.empty.names = FALSE)
    }
  }
  if (!is.null(ph)) {
    if (!is.na(ph)) {
      ph = ph/100
      minph <- (ph * (max(eic) - min(eic)))+min(eic)
      scans <- which(eic >= minph)
      peaktable$peaktable_int <- data.frame(peaktable$peaktable_int[,scans], fix.empty.names = FALSE)
      peaktable$peaktable_mass <- data.frame(peaktable$peaktable_mass[,scans], fix.empty.names = FALSE)
    }
  }
  if (!is.null(freq)) {
    if (!is.na(freq)) {
      tot <- ncol(peaktable$peaktable_mass)
      freq <- tot*freq/100
      ns <- apply(peaktable$peaktable_mass, 1, function(x) length(which(!is.na(x))))
      ind <- which(ns >= freq)
      peaktable$peaktable_int <- data.frame(peaktable$peaktable_int[ind,], fix.empty.names = FALSE)
      peaktable$peaktable_mass <- data.frame(peaktable$peaktable_mass[ind,], fix.empty.names = FALSE)
    }
  }
  peaktable$peaktable_int <- do.call(cbind, lapply(1:ncol(peaktable$peaktable_int), function(i) peaktable$peaktable_int[,i]/get(normfn)(peaktable$peaktable_int[,i], na.rm = TRUE)))
  mz <- apply(peaktable$peaktable_mass, 1, mean, na.rm = TRUE)
  mz.u <- apply(peaktable$peaktable_mass, 1, sd, na.rm = TRUE)
  int <- apply(peaktable$peaktable_int, 1, mean, na.rm = TRUE)
  int.u <- apply(peaktable$peaktable_int, 1, sd, na.rm = TRUE)
  n <- apply(peaktable$peaktable_mass, 1, function(x) length(which(!is.na(x))))
  x <- which(!is.nan(mz) & !is.nan(int))
  out <- data.frame(mz = mz[x], mz.u = mz.u[x], int = int[x], int.u = int.u[x], n = n[x])

  attr(out, "numscans") <- ncol(peaktable$peaktable_mass)

  out
}

#' Generate consensus mass spectrum plot
#'
#'
#' @param ms table uncertainty mass spectrum generated from `get_ums` function
#' @param xlim numeric vector containing the start and stop values of the x (m/z) axis (NULL if using default)
#' @param ylim numeric vector containing the start and stop values of the y (intensity) axis (NULL if using default)
#' @param main character string Title of the mass spectral plot
#' @param color character string color of the mass spectrum plot lines
#' @param size numeric width of the mass spectrum plot lines
#' @param removal numeric value of intensity to which remove mass spectrum peaks that are below specified intensity
#'
#' @import ggplot2
#' @importFrom utils packageVersion
#'
#' @return ggplot object
#' @export
#'
plot_ms <- function(ms, xlim = NULL, ylim = NULL, main = "Mass Spectrum", color = "black", size = 1, removal = 0) {
  ms <- ms[which(ms$int >= removal),]
  int <- ms$int
  intu <- ms$int.u
  mz <- ms$mz
  mzu <- ms$mz.u
  if (packageVersion("ggplot2") >= '3.4.0') {
    ggplot(data.frame(mz = mz, int = int)) + geom_linerange(aes(x = mz, ymin = 0, ymax = int), color = color, linewidth = size) + geom_pointrange(aes(x = mz, ymin = 0, ymax = int, y= int), shape = 20, size = 0.5) + geom_errorbar(aes(x = mz, ymin = int - intu, ymax = int + intu, width = 0.01), color = "red", na.rm = TRUE, linetype = 5) + geom_errorbarh(aes(y = int, xmin = mz - mzu, xmax = mz + mzu, height = 0.01), color = "red", na.rm = TRUE, linetype = 5) + ggtitle(main) + xlab("m/z") + ylab("Relative Intensity") + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)+ theme_bw()
  } else {
    ggplot(data.frame(mz = mz, int = int)) + geom_linerange(aes(x = mz, ymin = 0, ymax = int), color = color, size = size) + geom_pointrange(aes(x = mz, ymin = 0, ymax = int, y= int), shape = 20, size = 0.5) + geom_errorbar(aes(x = mz, ymin = int - intu, ymax = int + intu, width = 0.01), color = "red", na.rm = TRUE, linetype = 5) + geom_errorbarh(aes(y = int, xmin = mz - mzu, xmax = mz + mzu, height = 0.01), color = "red", na.rm = TRUE, linetype = 5) + ggtitle(main) + xlab("m/z") + ylab("Relative Intensity") + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)+ theme_bw()
  }
}


#' Write an uncertainty mass spectrum to a NIST .msp file
#'
#' @param ums uncertainty mass spectrum object generated by `get_ums` function
#' @param file character string to write the msp file to (include the extension)
#' @param precursor numeric value of the precursor ion m/z
#' @param name character string title to give the written mass spectrum
#' @param include_uncertainty logical to include the uncertainty calculations
#'
#' @importFrom utils write.table
#'
#' @return exports a mass spectrum to the designated file
#' @export
#'
write_ums <- function(ums, file, precursor = "", name = "R-exported Spectrum", include_uncertainty = TRUE) {
  peaks <- nrow(ums)
  names <- paste("Name: ", name)
  precursor <- paste("PrecursorMZ: ", precursor)
  labelorder <- paste0("LabelOrder: m/z, intensity")
  if (include_uncertainty) {  labelorder <- paste0("LabelOrder: m/z, intensity, m/z uncertainty, intensity uncertainty")}
  peaks <- paste("Num Peaks: ", peaks)
  header <- rbind(names, precursor, labelorder, peaks)
  write.table(header, file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (!is.na(ums[1,1])) {
    if (!include_uncertainty) {
      write.table(ums[,c(1,3)], file, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
    }
    if (include_uncertainty) {
      write.table(ums[,c(1,3,2,4)], file, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
    }
  }
  if (is.na(ums[1,1])) {
    write.table(cbind(0,0,0,0), file, append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  write("\n\n", file, append = TRUE)
  print(paste0("Mass spectrum written to file: ", file))
}

getEIC <- function(...) {
  warning("This function is deprecated, use `getEIC_from_peaklist` for future functions.")
  getEIC_from_peaklist(...)
}
