#!/usr/bin/env Rscript

#' wl-14-11-2022, Mon: XCMS Processing. Modified from Tom's R script.
#' wl-03-05-2023, Wed: Revise
#' wl-26-07-2023, Wed: Parallel computation. Only work for Linux and macOS
#' wl-28-07-2023, Fri: Parallel computation for windows (without 'bpstart')
#' wl-02-08-2023, Wed: Save peak detection results in case the parallel
#' computation fails in retention time correction.
#' -------------------install necessary packages
#' Only need to run this section on fresh installation
install.packages('tidyverse')
install.packages('tictoc')
install.package('BiocManager')
BiocManager::install('xcms')


#'---------------------Load packages-------------------------------------
rm(list = ls())
pkgs <- c("tidyverse", "xcms", "MSnbase", "parallel", "BiocParallel", "tictoc")
invisible(lapply(pkgs, library, character.only = TRUE))
# options(help_type = "html")


#' ----------------------- USER INPUT -----------------------------------
#' User Inputs: **User should change in this section**

#' mzML data file directory
in_dir <- "E:/Raw_data/Thapanee/19May2026/PoS/mzML/"
setwd(in_dir)
#' result directory
out_dir <- in_dir

#' Assay type description must be compatible to be used in a vector
#' assigning a directory e.g. no "\" "/" characters
assay_type <- "May26Pos"

#' A string common and unique to all QC sample file names
qc_pattern <- "QC"


#' A string common and unique to all blank sample file names
blank_pattern <- "Blank"

polarity <- 1    # polarity  1: positive; 0 - negative
ms <- 1             # msLevel cannot be NULL. Must be an integer, such as 1

sep <- ","              # change here for: ",", "\t" and " "; "," for csv, "\t" for tsv and " " for txt


#' ----------------------- PARAMETER SETTING ------------------------------
#' parallel computation
No <- parallel::detectCores() / 2 - 2   # default is 2
if (.Platform$OS.type == "windows") {
  #  register(SnowParam(workers = No, progressbar = TRUE))
  register(SnowParam(workers = 8, progressbar = TRUE)) #Specific for Intel CPU with 8 P-cores (i9-14900)
  #' register(SerialParam()) # works
  #' register(SnowParam())   # works
  #' register(bpstart(SnowParam(workers = No, progressbar = TRUE))) # No
} else {
  #  register(bpstart(MulticoreParam(workers = No,  progressbar = TRUE)))
  register(bpstart(MulticoreParam(workers = No,  progressbar = TRUE)))
}
#' check parallel
bpparam()
#' parameters for peak detection using the centWave method
cwp <- CentWaveParam(
    ppm = 25,
    peakwidth = c(5, 20),
    snthresh = 10,
    prefilter = c(3, 100),
    mzCenterFun = "wMean",
    integrate = 1L,
    mzdiff = 0.01,
    fitgauss = FALSE,
    noise = 2,
    verboseColumns = TRUE
  )


#' files to be processed
file_ext <- "\\.(mz[X]{0,1}ML|cdf)"    # no need to change
files <- list.files(in_dir, pattern = file_ext, full.names = TRUE,
                    recursive = FALSE, ignore.case = TRUE)
#files <- files[-25:-29]
# wl-12-05-2023, Fri: construct metadata
sample_name <- str_remove(basename(files), file_ext)
sample_group <- case_when(
  str_detect(sample_name, regex(qc_pattern, ignore_case = TRUE)) ~ "qc",
  str_detect(sample_name, regex(blank_pattern, ignore_case = TRUE)) ~ "blank",
  .default = "sample"
)

pd <- data.frame(sample_name = sample_name, sample_group = sample_group)

#' select middle QC as target for alignment

cqc_pattern <- "CQC"
qc_idx <- which(sample_group == "qc")
qc_name <- sample_name[qc_idx]
qc_group <- case_when(
     str_detect(qc_name, regex(cqc_pattern, ignore_case = TRUE)) ~ "cqc",
     .default = "qc"
)
cqc_idx <- which(qc_group == "cqc")
if (!isEmpty(cqc_idx))
  qc_idx <- qc_idx[-cqc_idx]

qc_num <- str_remove(sample_name[qc_idx],"QC")
qc_num <- str_remove(qc_num,"_Pos")
qc_num <- as.numeric(qc_num)
#centre_idx <- which(sample_name == paste("QC", floor(median(qc_num, na.rm = TRUE)), sep = ""))
#centre_idx <- which(sample_name == paste(paste("QC", floor(median(qc_num, na.rm = TRUE)), sep = ""), sep = ""))
centre_name <- paste(paste("QC", floor(median(qc_num, na.rm = TRUE)), sep = ""), sep = "")
centre_name <- paste(paste(centre_name, "_Pos", sep = ""), sep = "")
centre_idx <- which(sample_name == centre_name)

#' parameters for retention time adjustment based on the obiwarp method
#obiwp <- ObiwarpParam(centerSample = centre_idx, binSize = 1, gapInit = 0.4, gapExtend = 2.4)
obiwp <- ObiwarpParam(centerSample = centre_idx, binSize = 0.6, gapInit = 0.4, gapExtend = 2.4)

#' parameters for the peak grouping based on peak densities
pdp <- PeakDensityParam(sampleGroups = pd$sample_group,
                        bw = 5, minFraction = 0.25)

#' -------------------------- RUN ------------------------------
#' Read data for XCMS
print("Loading mzML data...")
data <- MSnbase::readMSData(files = files,
                            pdata = new("AnnotatedDataFrame", pd),
                            mode = "onDisk")

#' Polarity filtering (1: pos; 0 - neg)

if (!is.null(polarity)) {
  data <- MSnbase::filterPolarity(data, polarity = polarity)
}
data2 <- filterEmptySpectra(data)

#' Peak detection
print("Performing peak detection...")
tic()
xset <- findChromPeaks(data2, param = cwp, msLevel = ms)
cat("** Peak Detection Time **:  ")
toc()

#' Retention time correction/alignment (wl: for old version use 'rector')
print("Performing retention time correction...")
#' wl-02-08-2023, Wed: no parallel computation
register(SerialParam())  
tic()
xset <- adjustRtime(xset, param = obiwp)
cat("** Retention itime correction time **:  ")
toc()

#' Peak grouping(wl: for old version use 'group')
print("Performing peak grouping...")
xset <- groupChromPeaks(xset, param = pdp, msLevel = ms)

#' Extract the feature definitions
print("Preparing data for export...")
feature_definitions <- as.data.frame(featureDefinitions(xset)) %>%
  select(-contains("peakidx"))

#' Produce the intensity matrix
feature_matrix <- as.data.frame(featureValues(xset, value = "into"))

#' Combine feature definitions and intensity values
peak_table <- cbind(feature_definitions, feature_matrix)

#' Write data matrix table

if (sep == ",") {
  ext <- ".csv"
} else if (sep == "\t") {
  ext <- ".tsv"
} else {                 # for default sep = " " or anything else
  ext <- ".txt"
}
#' add rownames into data frame
tab <- peak_table %>% rownames_to_column(., var = "name")
write.table(tab,
            file = paste0(out_dir, assay_type, ext),
            sep = sep, row.names = FALSE, na = "")

# Fill missing peaks
xset_fill <- fillChromPeaks(xset)
print("Preparing data for export...")
feature_definitions <- as.data.frame(featureDefinitions(xset_fill)) %>%
  select(-contains("peakidx"))

#' Produce the intensity matrix
feature_matrix <- as.data.frame(featureValues(xset_fill, value = "into"))

#' Combine feature definitions and intensity values
peak_table <- cbind(feature_definitions, feature_matrix)

#' Write data matrix table
sep <- ","              # change here for: ",", "\t" and " "
if (sep == ",") {
  ext <- ".csv"
} else if (sep == "\t") {
  ext <- ".tsv"
} else {                 # for default sep = " " or anything else
  ext <- ".txt"
}
#' add rownames into data frame
tab <- peak_table %>% rownames_to_column(., var = "name")
write.table(tab,
            file = paste0(out_dir, assay_type, "_fill" ,ext),
            sep = sep, row.names = FALSE, na = "")
save.image("Results_final.RData")

