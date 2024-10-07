# Project:      Routine cold storage leads to hyperacute graft loss in pig-to-primate kidney xenotransplantation; hypothermic machine perfusion may be preferred preservation modality in xenotransplantation
# Description:  Data Pre-Processing (Quality Control)
# Author:       Adam Luo

# Clear environment
rm(list = ls())

# Load packages
library(stats)
library(tidyverse)
library(readxl)
library(writexl)
library(utils)
library(scales)
library(BiocManager)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(reshape2)
library(preprocessCore)
library(umap)
library(GSVA)
library(DESeq2)
library(ggsignif)
library(ggprism)
library(ggpubr)
library(ggrepel)

# Set working directory
setwd(path.expand("~"))

# Load data into GeoMxSet object
datadir <- file.path("./Data/RawFiles")
DCCFiles <- dir(file.path(datadir, "dcc"),
                pattern = ".dcc$",
                full.names = TRUE,
                recursive = TRUE)
PKCFile <- dir(file.path(datadir, "pkc"),
               pattern = ".pkc$",
               full.names = TRUE)
AnnotationFile <- dir(file.path(datadir, "annotation"),
                      pattern = ".xlsx$",
                      full.names= TRUE)

dsp_raw <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                  pkcFiles = PKCFile,
                                  phenoDataFile = AnnotationFile,
                                  phenoDataSheet = "Template",
                                  phenoDataDccColName = "Sample_ID")

# Shift count values of 0 to 1
dsp_one <- shiftCountsOne(dsp_raw, useDALogic = TRUE)

# Add quality control flags for probe-level data
dsp_one <- setBioProbeQCFlags(dsp_one,
                              qcCutoffs = list(minProbeRatio = 0.1,
                                               percentFailGrubbs = 20), 
                              removeLocalOutliers = FALSE)
ProbeQCResults <- fData(dsp_one)[["QCFlags"]]
ProbeQCPassed <- 
  subset(dsp_one, 
         fData(dsp_one)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(dsp_one)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dsp_QC <- ProbeQCPassed

# Collapse probe-level data to target-level data
target_dsp_QC <- aggregateCounts(dsp_QC)

# Remove "no template control" data
target_dsp_QC <- target_dsp_QC[, pData(target_dsp_QC)$'Slide Name' != "No Template Control"]

# Calculate and store LOQ values for each ROI
cutoff <- 2
minLOQ <- 2
LOQ <- data.frame(row.names = colnames(target_dsp_QC))
pkcs <- annotation(dsp_raw)
module <- gsub(".pkc", "", pkcs)
vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
if(all(vars[1:2] %in% colnames(pData(target_dsp_QC)))) {
  LOQ[, module] <-
    pmax(minLOQ,
         pData(target_dsp_QC)[, vars[1]] *
           pData(target_dsp_QC)[, vars[2]] ^ cutoff)}
pData(target_dsp_QC)$LOQ <- LOQ$Hs_R_NGS_WTA_v1.0

# Perform manual segment filtering
target_dsp_QC <- target_dsp_QC[, pData(target_dsp_QC)$Genetics == "sGalKO-Xeno"]

# Generate Boolean matrix of gene-level counts greater than LOQ per segment
LOQ_Mat <- data.frame(row.names = colnames(target_dsp_QC))
LOQ_Mat <- t(LOQ_Mat)
Mat_i <- t(esApply(target_dsp_QC, MARGIN = 1,
                   FUN = function(x) {x > LOQ}))
LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
LOQ_Mat <- LOQ_Mat[fData(target_dsp_QC)$TargetName, ]

# Save segment-level gene detection rate data to phenoData
pData(target_dsp_QC)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_dsp_QC)$GeneDetectionRate <-
  pData(target_dsp_QC)$GenesDetected / nrow(target_dsp_QC)

# Generate target-level segment detection rate data
LOQ_Mat <- LOQ_Mat[, colnames(target_dsp_QC)]
fData(target_dsp_QC)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_dsp_QC)$SegmentDetectionRate <-
  fData(target_dsp_QC)$DetectedSegments / nrow(pData(target_dsp_QC))

# Set filtering threshold & create new object for filtered data
threshold <- c(5) # percent (can modify & add more thresholds)
for (i in 1:length(threshold)) {
  new <- paste0("target_dsp_F", threshold[i])
  assign(new, target_dsp_QC[fData(target_dsp_QC)$SegmentDetectionRate >= (threshold[i] / 100) |
                               fData(target_dsp_QC)$TargetName %in% "NegProbe-WTX", 
                             pData(target_dsp_QC)$GeneDetectionRate >= (threshold[i] / 100)])
}

# Choose the dataset for downstream analysis
target_dsp_F <- target_dsp_F5
