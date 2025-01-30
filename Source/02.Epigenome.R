## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/functions.R")
###########################################################################
## To load various epigenetic tracks from UCSC, ENCODE, etc. databases
## Input: bigWig; bed (broad vs. narrow peak); 
## Species: human (hg38), mouse (mm10), (other genome builds are converted by CrossMap.py)
## K562 Tracks: 
##    - scATAC-seq narrow and broad peaks
##    - DNase narrow and broad peaks 
##    - FAIRE narrow and broad peaks
##    - RNA- and scRNA-seq, exonic and intronic expression level
##    - GRO-seq, exonic and intronic expression, inferred transcripts (regions)
##    - PRO-seq, inferred transcripts (regions)
##    - GRO-cap, inferred TSS (regions)
##    - Methylome: no/low/mid/high/full methylation regions
##    - Broad and narrow histone marks: narrow and broad peaks, respectively
##    - ChIP-seq TFBS locations: narrow peaks
###########################################################################
## 1. K562
## 1.1 single-cell K562 ATAC-seq From Greenleaf et al. 
## GSE65360; peaks called by MACS2
K562ATACSampleMeta <- read.csv("Data/GSE65360/info/SraRunTable.txt", as.is = TRUE, check.names = FALSE)
K562ATACSampleMeta <- subset(K562ATACSampleMeta, source_name == "singles-K562-rep1" | source_name == "singles-K562-rep2" | source_name == "singles-K562-rep3")
K562ATACSampleMeta <- K562ATACSampleMeta[!duplicated(K562ATACSampleMeta[["Sample Name"]]), ]
rownames(K562ATACSampleMeta) <- K562ATACSampleMeta[["Sample Name"]]
K562ATACSampleIDs <- c(paste("GSM", 1596831:1597118, sep = ""), "K562UntreatedRep1", "K562UntreatedRep2",  "K562UntreatedRep3", "K562UntreatedAll")
K562ATACNarrowPeaks <- sapply(K562ATACSampleIDs, function(sampleID) {
    message(sampleID)
    Genome$standardizeSeqInfo(Genome$import.narrowPeak(sprintf("Data/GSE65360/analyzed/MACS2/%s/filtered_peaks.narrowPeak", sampleID)), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human) 
}, simplify = FALSE)
K562ATACBroadPeaks <- sapply(K562ATACSampleIDs, function(sampleID) { 
    message(sampleID) 
    Genome$standardizeSeqInfo(Genome$import.broadPeak(sprintf("Data/GSE65360/analyzed/MACS2/%s/filteredBroad_peaks.broadPeak", sampleID)), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human) 
}, simplify = FALSE)
## Because MACS2 generates multiple peaks spanning same region but with different summits. Here we don't 
## care about these "summit variants", so we just filter out the duplicate peaks.
K562ATACNarrowPeaks <- lapply(K562ATACNarrowPeaks, Genome$removeDupPeaks)
K562ATACSignals <- sapply(K562ATACSampleIDs, function(sampleID) {
    message(sampleID) 
    Genome$standardizeSeqInfo(rtracklayer::import(sprintf("Data/GSE65360/analyzed/bowtie2X2k/%s/filtered.bigWig", sampleID)), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
}, simplify = FALSE)
saveRDS(K562ATACSampleMeta, file = "Data/Epigenome/K562ATACSampleMeta.RDS")
saveRDS(K562ATACSampleIDs, file = "Data/Epigenome/K562ATACSampleIDs.RDS")
saveRDS(K562ATACNarrowPeaks, file = "Data/Epigenome/K562ATACNarrowPeaks.RDS")
dir.create("Report/Epigenome/K562ATACNarrowPeaks", FALSE, TRUE)
rtracklayer::export(K564ATACNarrowPeaks$K562UntreatedAll, con = "Report/Epigenome/K562ATACNarrowPeaks/K562UntreatedAll.bed.gz")
dir.create("Data/Epigenome", FALSE, TRUE)
saveRDS(K562ATACBroadPeaks, file = "Data/Epigenome/K562ATACBroadPeaks.RDS")
saveRDS(K562ATACSignals, file = "Data/Epigenome/K562ATACSignals.RDS")

## 1.2 DNase-seq from ENCODE
K562DNaseFiles <- list(
    Rep1 = c(Signal = "ENCFF868NHV.bigWig", NarrowPeak = "ENCFF941GAQ.bed", BroadPeak = "ENCFF581VXJ.bed", Bam = "ENCFF823AJQ.bam"), 
    Rep2 = c(Signal = "ENCFF227MPB.bigWig", NarrowPeak = "ENCFF196JWI.bed", BroadPeak = "ENCFF416PBA.bed", Bam = "ENCFF851DHZ.bam")
)
K562DNaseNarrowPeaks <- sapply(c("Rep1", "Rep2"), function(ID) { 
    Genome$standardizeSeqInfo(Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/DNase-seq/hg38/K562/%s", K562DNaseFiles[[ID]]["NarrowPeak"])), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human) 
}, simplify = FALSE)
K562DNaseBroadPeaks <- sapply(c("Rep1", "Rep2"), function(ID) { 
    Genome$standardizeSeqInfo(Genome$import.broadPeak(sprintf("Data/Database/ENCODE/DNase-seq/hg38/K562/%s", K562DNaseFiles[[ID]]["BroadPeak"])), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human) 
}, simplify = FALSE)

## Note, broadPeaks have -1 as signalValue and pValue
K562DNaseSignals <- sapply(c("Rep1", "Rep2"), function(ID) { 
    Genome$standardizeSeqInfo(rtracklayer::import(sprintf("Data/Database/ENCODE/DNase-seq/hg38/K562/%s", K562DNaseFiles[[ID]]["Signal"])), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
}, simplify = FALSE)
saveRDS(K562DNaseNarrowPeaks, file = "Data/Epigenome/K562DNaseNarrowPeaks.RDS")
dir.create("Report/Epigenome/K562DNaseNarrowPeaks", FALSE, TRUE)
rtracklayer::export(K562DNaseNarrowPeaks$Rep1, con = "Report/Epigenome/K562DNaseNarrowPeaks/ENCFF941GAQ.bed.gz")
rtracklayer::export(K562DNaseNarrowPeaks$Rep2, con = "Report/Epigenome/K562DNaseNarrowPeaks/ENCFF196JWI.bed.gz")
saveRDS(K562DNaseBroadPeaks, file = "Data/Epigenome/K562DNaseBroadPeaks.RDS")
saveRDS(K562DNaseSignals, file = "Data/Epigenome/K562DNaseSignals.RDS")

## 1.3 FAIRE-seq
## Note, we don't have broadPeaks for FAIRE-seq
K562FAIRENarrowPeaks <- list(standardizeSeqInfo(Genome$import.narrowPeak("Data/Database/ENCODE/FAIRE-seq/hg38/K562/ENCFF000TLT.bed"), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human))
## Note, narrowPeaks have -1 as pValue.
K562FAIRESignals <- list(standardizeSeqInfo(rtracklayer::import(sprintf("Data/Database/ENCODE/FAIRE-seq/hg38/K562/%s", "ENCFF000TLE.bigWig")), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human))
saveRDS(K562FAIRENarrowPeaks, file = "Data/Epigenome/K562FAIRENarrowPeaks.RDS")
dir.create("Report/Epigenome/K562FAIRENarrowPeaks", FALSE, TRUE)
rtracklayer::export(K562FAIRENarrowPeaks[[1]], con = "Report/Epigenome/K562FAIRENarrowPeaks/ENCFF000TLT.bed.gz")
saveRDS(K562FAIRESignals, file = "Data/Epigenome/K562FAIRESignals.RDS")

## 1.4 RRBS CpG methylation
## As of rtrackyaler v1.42.2 we don't need import.bedMethyl any more.
K562RRBSMethy <- list(
    Rep1 = rtracklayer::import("Data/Database/wgEncode/MethylRRBS/hg38/wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed"), 
    Rep2 = rtracklayer::import("Data/Database/wgEncode/MethylRRBS/hg38/wgEncodeHaibMethylRrbsK562HaibSitesRep2.bed")
)
K562RRBSMethy <- lapply(K562RRBSMethy, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
saveRDS(K562RRBSMethy, file = "Data/Epigenome/K562RRBSMethy.RDS")

## Stratify RRBS methylation data into four categories: low, mid-low, mid-high, high
## 1. merge Rep1 and Rep2 first
K562RRBSMethy <- unname(unlist(as(K562RRBSMethy, "GRangesList")))
## 2. filter out loci with less than 5 coverage
mean(score(K562RRBSMethy) < 10)
## [1] 0.3424202
K562RRBSMethy <- K562RRBSMethy[score(K562RRBSMethy) >= 10]
K562RRBSMethy$methylPerc <- as.integer(mcols(K562RRBSMethy)[["blockSizes"]])
methylLevels <- cut(0:100, breaks = c(0, 25, 50, 75, 100), right = TRUE)
methylLevels[1] <- "(0,25]"
names(methylLevels) <- c('(0,25]' = "Low", '(25,50]' = "MidLow", '(50,75]' = "MidHigh", '(75,100)' = "High")[methylLevels]
K562RRBSMethy$methylLevel <- names(methylLevels[K562RRBSMethy$methylPerc + 1])
K562RRBSMethyLevels <- split(K562RRBSMethy, f = K562RRBSMethy$methylLevel)
K562RRBSMethyLevels <- K562RRBSMethyLevels[c("Low", "MidLow", "MidHigh", "High")]
lengths(K562RRBSMethyLevels)
##     Low  MidLow MidHigh    High
## 1036075  114318  106518  315017
saveRDS(K562RRBSMethyLevels, file = "Data/Epigenome/K562RRBSMethyLevels.RDS")
dir.create("Report/Epigenome/K562RRBSMethyLevels", FALSE, TRUE)
for (q in c("Low", "MidLow", "MidHigh", "High")) {
    rtracklayer::export(K562RRBSMethyLevels[[q]], con = sprintf("Report/Epigenome/K562RRBSMethyLevels/%s.bed.gz", q))
}

## 1.5 Replication origin from GSE46189
K562RepOrigPeaks <- list(Genome$standardizeSeqInfo(rtracklayer::import("Data/GSE46189/GEO/hg38/GSE46189_Ori-Peak.bed"), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human))
saveRDS(K562RepOrigPeaks, file = "Data/Epigenome/K562RepOrigPeaks.RDS")
dir.create("Report/Epigenome/K562RepOrigPeaks", FALSE, TRUE)
rtracklayer::export(K562RepOrigPeaks[[1]], con = "Report/Epigenome/K562RepOrigPeaks/GSE46189.bed.gz")

## 1.6 ChIP-seq
## 1.6.1 TFBS 
K562ChIPMetadataTF <- read.csv("Data/Database/ENCODE/ChIP-seq/hg38/K562/TF/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(K562ChIPMetadataTF) <- K562ChIPMetadataTF[, "File accession"]
K562ChIPMetadataTF <- subset(K562ChIPMetadataTF, `File Status` == "released" & Assembly == "GRCh38")
K562ChIPSampleIDsTF <- subset(K562ChIPMetadataTF, `Output type` == "optimal idr thresholded peaks")[, "File accession"]
isDuplicated <- duplicated(K562ChIPMetadataTF[K562ChIPSampleIDsTF, "Experiment target"])
K562ChIPSampleIDsTF <- K562ChIPSampleIDsTF[!isDuplicated]
K562ChIPMetadataTF <- data.frame(SampleID = K562ChIPSampleIDsTF, Target = K562ChIPMetadataTF[K562ChIPSampleIDsTF, "Experiment target"], stringsAsFactors = FALSE)
K562ChIPMetadataTF$Target <- sub("-human", "", K562ChIPMetadataTF$Target)
dim(K562ChIPMetadataTF)
## [1] 322   2
## Further remove eGFP tagged TFs
K562ChIPMetadataTF <- subset(K562ChIPMetadataTF, !grepl("^eGFP", Target))
dim(K562ChIPMetadataTF)
## [1] 258   2
K562ChIPTFs <- with(K562ChIPMetadataTF, sapply(SampleID, function(sampleID) {
    message(sampleID)
    GRs <- Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/ChIP-seq/hg38/K562/TF/%s.bed", sampleID)) 
    Genome$standardizeSeqInfo(GRs, seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
}, simplify = FALSE))
names(K562ChIPTFs) <- K562ChIPMetadataTF[, "Target"]
write.csv(K562ChIPMetadataTF, file = "Report/Epigenome/K562ChIPMetadataTF.csv")
saveRDS(K562ChIPTFs, file = "Data/Epigenome/K562ChIPTFs.RDS")
dir.create("Report/Epigenome/K562ChIPTFs", FALSE, TRUE)
for (x in names(K562ChIPTFs)) {
    message(x)
    rtracklayer::export(K562ChIPTFs[[x]], con = sprintf("Report/Epigenome/K562ChIPTFs/%s.bed.gz", x))
}

## 1.6.2 Broad histone marks
K562ChIPMetadataBHM <- read.csv("Data/Database/ENCODE/ChIP-seq/hg38/K562/BroadHistoneMark/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(K562ChIPMetadataBHM) <- K562ChIPMetadataBHM[, "File accession"]
K562ChIPMetadataBHM <- subset(K562ChIPMetadataBHM, `File Status` == "released" & Assembly == "GRCh38")
K562ChIPSampleIDsBHM <- subset(K562ChIPMetadataBHM, `Output type` == "replicated peaks" | `Output type` == "stable peaks")[, "File accession"]
isDuplicated <- duplicated(K562ChIPMetadataBHM[K562ChIPSampleIDsBHM, "Experiment target"])
K562ChIPSampleIDsBHM <- K562ChIPSampleIDsBHM[!isDuplicated]
K562ChIPMetadataBHM <- data.frame(SampleID = K562ChIPSampleIDsBHM, Target = K562ChIPMetadataBHM[K562ChIPSampleIDsBHM, "Experiment target"], stringsAsFactors = FALSE)
K562ChIPMetadataBHM$Target <- sub("-human", "", K562ChIPMetadataBHM$Target)
dim(K562ChIPMetadataBHM)
## [1] 6 2
K562ChIPBHMs <- with(K562ChIPMetadataBHM, sapply(SampleID, function(sampleID) {
    message(sampleID)
    GRs <- Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/ChIP-seq/hg38/K562/BroadHistoneMark/%s.bed", sampleID)) 
    Genome$standardizeSeqInfo(GRs, seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
}, simplify = FALSE))
names(K562ChIPBHMs) <- K562ChIPMetadataBHM[, "Target"]
write.csv(K562ChIPMetadataBHM, file = "Report/Epigenome/K562ChIPMetadataBHM.csv")
saveRDS(K562ChIPBHMs, file = "Data/Epigenome/K562ChIPBHMs.RDS")
dir.create("Report/Epigenome/K562ChIPBHMs", FALSE, TRUE)
for (x in names(K562ChIPBHMs)) {
    message(x)
    rtracklayer::export(K562ChIPBHMs[[x]], con = sprintf("Report/Epigenome/K562ChIPBHMs/%s.bed.gz", x))
}

## 1.6.3 Narrow histone marks
K562ChIPMetadataNHM <- read.csv("Data/Database/ENCODE/ChIP-seq/hg38/K562/NarrowHistoneMark/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(K562ChIPMetadataNHM) <- K562ChIPMetadataNHM[, "File accession"]
K562ChIPMetadataNHM <- subset(K562ChIPMetadataNHM, `File Status` == "released" & Assembly == "GRCh38")
K562ChIPSampleIDsNHM <- subset(K562ChIPMetadataNHM, `Output type` == "replicated peaks")[, "File accession"]
isDuplicated <- duplicated(K562ChIPMetadataNHM[K562ChIPSampleIDsNHM, "Experiment target"])
K562ChIPSampleIDsNHM <- K562ChIPSampleIDsNHM[!isDuplicated]
K562ChIPMetadataNHM <- data.frame(SampleID = K562ChIPSampleIDsNHM, Target = K562ChIPMetadataNHM[K562ChIPSampleIDsNHM, "Experiment target"], stringsAsFactors = FALSE)
K562ChIPMetadataNHM$Target <- sub("-human", "", K562ChIPMetadataNHM$Target)
dim(K562ChIPMetadataNHM)
## [1] 4 2
K562ChIPNHMs <- with(K562ChIPMetadataNHM, sapply(SampleID, function(sampleID) {
    message(sampleID)
    GRs <- Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/ChIP-seq/hg38/K562/NarrowHistoneMark/%s.bed", sampleID)) 
    Genome$standardizeSeqInfo(GRs, seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
}, simplify = FALSE))
names(K562ChIPNHMs) <- K562ChIPMetadataNHM[, "Target"]

K562ChIPNHMsH3K27ac <- Genome$standardizeSeqInfo(Genome$import.narrowPeak("Data/Database/ENCODE/ChIP-seq/hg38/K562/NarrowHistoneMark/ENCFF045OHM.bed"), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
K562ChIPNHMs <- c(K562ChIPNHMs, list(H3K27ac = K562ChIPNHMsH3K27ac))
write.csv(K562ChIPMetadataNHM, file = "Report/Epigenome/K562ChIPMetadataNHM.csv")
saveRDS(K562ChIPNHMs, file = "Data/Epigenome/K562ChIPNHMs.RDS")
dir.create("Report/Epigenome/K562ChIPNHMs", FALSE, TRUE)
for (x in names(K562ChIPNHMs)) {
    message(x)
    rtracklayer::export(K562ChIPNHMs[[x]], con = sprintf("Report/Epigenome/K562ChIPNHMs/%s.bed.gz", x))
}

## 1.7 GRO-seq nascent RNA run-on sequencing data
K562GROFiles <- list(Plus = "Data/GSE60454/GEO/hg38/GSM1480325_K562_GROseq_plus.bigWig", Minus = "Data/GSE60454/GEO/hg38/GSM1480325_K562_GROseq_minus.bigWig")
K562GROSignals <- lapply(K562GROFiles, function(filename) Genome$standardizeSeqInfo(rtracklayer::import(filename), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human))
K562GROPeaks <- list(Homer = rtracklayer::import("Data/GSE60454/analyzed/Sample_GSM1480325/HOMERnoSplice/TagDir/transcripts.gtf"))
K562GROPeaks <- lapply(K562GROPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)
saveRDS(K562GROPeaks, file = "Data/Epigenome/K562GROPeaks.RDS")
saveRDS(K562GROSignals, file = "Data/Epigenome/K562GROSignals.RDS")
dir.create("Report/Epigenome/K562GROPeaks", FALSE, TRUE)
rtracklayer::export(K562GROPeaks[[1]], con = sprintf("Report/Epigenome/K562GROPeaks/%s.bed.gz", "GSM1480325"))

## 1.8 super-enhancer
K562SEPeaks <- sort(Genome$standardizeSeqInfo(import("Data/Database/dbSUPER/hg38/K562.bed"), seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE))
score(K562SEPeaks) <- 1 ## Original score is peak number index, we need to replace this misleading information.
saveRDS(K562SEPeaks, file = "Data/Epigenome/K562SEPeaks.RDS")
dir.create("Report/Epigenome/K562SEPeaks", FALSE, TRUE)
rtracklayer::export(K562SEPeaks, con = sprintf("Report/Epigenome/K562SEPeaks/%s.bed.gz", "dbSUPER"))

## 1.9 DRIP-seq R-loop DNA
K562RloopPeaks <- list(GSM1720619 = sort(Genome$standardizeSeqInfo(rtracklayer::import("Data/GSE70189/GEO/hg38/GSM1720619_K562_DRIP_peaks.bed"), seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)))
saveRDS(K562RloopPeaks, file = "Data/Epigenome/K562RloopPeaks.RDS")
K562RloopSignals <- list(GSM1720619 = sort(Genome$standardizeSeqInfo(rtracklayer::import("Data/Datasets/GSE70189/GEO/hg38/GSM1720619_K562_DRIP.bigWig"), seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)))
saveRDS(K562RloopSignals, file = "Data/Epigenome/K562RloopSignals.RDS")
dir.create("Report/Epigenome/K562RloopPeaks", FALSE, TRUE)
rtracklayer::export(K562RloopPeaks[["GSM1720619"]], con = "Report/Epigenome/K562RloopPeaks/GSM1720619.bed.gz")

###########################################################################
## 2. Brain cells
###########################################################################
## 2.1 human brain
## 2.1.1 human brain ATAC-seq
## Since we have only one sample, we don't need metadata; just do it manually.
HumanBrainATACNarrowPeaks <- list(CerebellumMaleAdult20years = Genome$import.narrowPeak("Data/Database/ENCODE/ATAC-seq/hg38/Brain/CerebellumMaleAdult20years/MACS2/ENCFF114WRE_peaks.narrowPeak"))
HumanBrainATACNarrowPeaks <- lapply(HumanBrainATACNarrowPeaks,  Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)
HumanBrainATACNarrowPeaks <- lapply(HumanBrainATACNarrowPeaks, Genome$removeDupPeaks)
HumanBrainATACBroadPeaks <- list(CerebellumMaleAdult20years = Genome$import.broadPeak("Data/Database/ENCODE/ATAC-seq/hg38/Brain/CerebellumMaleAdult20years/MACS2/ENCFF114WRE_peaks.broadPeak"))
HumanBrainATACBroadPeaks <- lapply(HumanBrainATACBroadPeaks,  Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)
HumanBrainATACBroadPeaks <- lapply(HumanBrainATACBroadPeaks, Genome$removeDupPeaks)
saveRDS(HumanBrainATACNarrowPeaks, file = "Data/Epigenome/HumanBrainATACNarrowPeaks.RDS")
saveRDS(HumanBrainATACBroadPeaks, file = "Data/Epigenome/HumanBrainATACBroadPeaks.RDS")
HumanBrainATACNarrowPeaks <- readRDS("Data/Epigenome/HumanBrainATACNarrowPeaks.RDS")
dir.create("Report/Epigenome/HumanBrainATACNarrowPeaks", FALSE, TRUE)
rtracklayer::export(HumanBrainATACNarrowPeaks$CerebellumMaleAdult20years, con = "Report/Epigenome/HumanBrainATACNarrowPeaks/CerebellumMaleAdult20years.bed.gz")

## 2.1.2 human brain DNase-seq
HumanBrainDNaseMetadata <- read.csv("Data/Database/ENCODE/DNase-seq/hg38/Brain/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
HumanBrainDNaseMetadata <- subset(HumanBrainDNaseMetadata, `Assembly` == "GRCh38" & `File Status` == "released")
HumanBrainDNaseInDir <- "Data/Database/ENCODE/DNase-seq/hg38/Brain"
HumanBrainDNaseNarrowPeaks <- list(
    CerebellumAstrocyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed narrowPeak" & `Biosample term name` == "astrocyte of the cerebellum"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.narrowPeak(sprintf("%s/%s.narrowPeak", HumanBrainDNaseInDir, x)))), 
    HippocampusAstrocyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed narrowPeak" & `Biosample term name` == "astrocyte of the hippocampus"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.narrowPeak(sprintf("%s/%s.narrowPeak", HumanBrainDNaseInDir, x)))), 
    MicrovascularEndothelialCell = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed narrowPeak" & `Biosample term name` == "brain microvascular endothelial cell"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.narrowPeak(sprintf("%s/%s.narrowPeak", HumanBrainDNaseInDir, x)))), 
    ChoroidPlexusEpithelialCell = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed narrowPeak" & `Biosample term name` == "choroid plexus epithelial cell"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.narrowPeak(sprintf("%s/%s.narrowPeak", HumanBrainDNaseInDir, x)))), 
    Pericyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed narrowPeak" & `Biosample term name` == "brain pericyte"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.narrowPeak(sprintf("%s/%s.narrowPeak", HumanBrainDNaseInDir, x))))
)
HumanBrainDNaseNarrowPeaks <- lapply(HumanBrainDNaseNarrowPeaks, function(X) "names<-"(X, paste0("Rep", names(X))))
HumanBrainDNaseNarrowPeaks <- rapply(HumanBrainDNaseNarrowPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE, how = "replace")
HumanBrainDNaseBroadPeaks <- list(
    CerebellumAstrocyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed broadPeak" & `Biosample term name` == "astrocyte of the cerebellum"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.broadPeak(sprintf("%s/%s.broadPeak", HumanBrainDNaseInDir, x)))), 
    HippocampusAstrocyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed broadPeak" & `Biosample term name` == "astrocyte of the hippocampus"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.broadPeak(sprintf("%s/%s.broadPeak", HumanBrainDNaseInDir, x)))), 
    MicrovascularEndothelialCell = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed broadPeak" & `Biosample term name` == "brain microvascular endothelial cell"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.broadPeak(sprintf("%s/%s.broadPeak", HumanBrainDNaseInDir, x)))), 
    ChoroidPlexusEpithelialCell = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed broadPeak" & `Biosample term name` == "choroid plexus epithelial cell"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.broadPeak(sprintf("%s/%s.broadPeak", HumanBrainDNaseInDir, x)))), 
    Pericyte = with(subset(HumanBrainDNaseMetadata, `File format` == "bigBed broadPeak" & `Biosample term name` == "brain pericyte"), tapply(`File accession`, INDEX = `Biological replicate(s)`, function(x) Genome$import.broadPeak(sprintf("%s/%s.broadPeak", HumanBrainDNaseInDir, x))))
)
HumanBrainDNaseBroadPeaks <- lapply(HumanBrainDNaseBroadPeaks, function(X) "names<-"(X, paste0("Rep", names(X))))
HumanBrainDNaseBroadPeaks <- rapply(HumanBrainDNaseBroadPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE, how = "replace")
write.csv(HumanBrainDNaseMetadata, file = "Report/Epigenome/HumanBrainDNaseMetadata.csv")
saveRDS(HumanBrainDNaseNarrowPeaks, file = "Data/Epigenome/HumanBrainDNaseNarrowPeaks.RDS")
saveRDS(HumanBrainDNaseBroadPeaks, file = "Data/Epigenome/HumanBrainDNaseBroadPeaks.RDS")
HumanBrainDNaseNarrowPeaks <- readRDS(file = "Data/Epigenome/HumanBrainDNaseNarrowPeaks.RDS")
dir.create("Report/Epigenome/HumanBrainDNaseNarrowPeaks", FALSE, TRUE)
for (x in names(HumanBrainDNaseNarrowPeaks)) { 
    rtracklayer::export(Reduce(union, HumanBrainDNaseNarrowPeaks[[x]]), con = sprintf("Report/Epigenome/HumanBrainDNaseNarrowPeaks/%s.bed.gz", x)) 
}

## 2.1.3 human brain FAIRE-seq
HumanBrainFAIRENarrowPeaks <- list(Astrocyte = Genome$import.narrowPeak("Data/Database/ENCODE/FAIRE-seq/hg38/Brain/Astrocycyte/ENCFF000TOD.narrowPeak"))
HumanBrainFAIRENarrowPeaks <- lapply(HumanBrainFAIRENarrowPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE)
saveRDS(HumanBrainFAIRENarrowPeaks, file = "Data/Epigenome/HumanBrainFAIRENarrowPeaks.RDS")
HumanBrainFAIRENarrowPeaks <- readRDS(file = "Data/Epigenome/HumanBrainFAIRENarrowPeaks.RDS")
dir.create("Report/Epigenome/HumanBrainFAIRENarrowPeaks", FALSE, TRUE)
rtracklayer::export(HumanBrainFAIRENarrowPeaks$Astrocyte, con = "Report/Epigenome/HumanBrainFAIRENarrowPeaks/Astrocyte.bed.gz")

## 2.1.4 human brain WGBS
## Only SK-N-SH available, this is neuronblastoma cell line, since cancer cell line has distinct methylome than normal cells, we don't include the data here.

## 2.2 Mouse brain (cerebellum, whole brain)
## 2.2.1 Mouse brain ATAC-seq
MouseBrainATACInDir <- "Data/Database/ENCODE/ATAC-seq/mm10/Brain"
MouseBrainATACTissueTypes <- list.dirs("Data/Database/ENCODE/ATAC-seq/mm10/Brain", recursive = FALSE, full.names = FALSE)
MouseBrainATACMetadata <- read.csv("Data/Database/ENCODE/ATAC-seq/mm10/Brain/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
MouseBrainATACMetadata <- subset(MouseBrainATACMetadata, `Assembly` == "mm10" & `File Status` == "released")
MouseBrainATACAccToRep <- MouseBrainATACMetadata[, "Biological replicate(s)"]
names(MouseBrainATACAccToRep) <- MouseBrainATACMetadata[, "File accession"]
## Here we merge .bam and call peaks from the merged data
MouseBrainATACNarrowPeaks <- sapply(MouseBrainATACTissueTypes, function(tissue) {
    message(tissue) 
    peaks <- sapply("merged", function(acc) { 
    Genome$import.narrowPeak(sprintf("%s/%s/MACS2/%s_peaks.narrowPeak", MouseBrainATACInDir, tissue, acc)) 
    }, simplify = FALSE) 
}, simplify = FALSE)
MouseBrainATACNarrowPeaks <- rapply(MouseBrainATACNarrowPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace") 
MouseBrainATACNarrowPeaks <- sapply(MouseBrainATACTissueTypes, function(tissue) MouseBrainATACNarrowPeaks[[tissue]]$merged, simplify = FALSE)
MouseBrainATACNarrowPeaks <- sapply(MouseBrainATACTissueTypes, function(tissue) { 
    message(tissue) 
    Genome$removeDupPeaks(MouseBrainATACNarrowPeaks[[tissue]]) 
}, simplify = FALSE)
MouseBrainATACBroadPeaks <- sapply(MouseBrainATACTissueTypes, function(tissue) { 
    message(tissue) 
    peaks <- sapply("merged", function(acc) { 
    Genome$import.broadPeak(sprintf("%s/%s/MACS2/%s_peaks.broadPeak", MouseBrainATACInDir, tissue, acc)) 
    }, simplify = FALSE) 
}, simplify = FALSE)
MouseBrainATACBroadPeaks <- rapply(MouseBrainATACBroadPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace") 
MouseBrainATACBroadPeaks <- sapply(MouseBrainATACTissueTypes, function(tissue) MouseBrainATACBroadPeaks[[tissue]]$merged, simplify = FALSE)
write.csv(MouseBrainATACMetadata, file = "Report/Epigenome/MouseBrainATACMetadata.csv")
saveRDS(MouseBrainATACNarrowPeaks, file = "Data/Epigenome/MouseBrainATACNarrowPeaks.RDS")
saveRDS(MouseBrainATACBroadPeaks, file = "Data/Epigenome/MouseBrainATACBroadPeaks.RDS")

MouseBrainATACNarrowPeaks <- readRDS(file = "Data/Epigenome/MouseBrainATACNarrowPeaks.RDS")
dir.create("Report/Epigenome/MouseBrainATACNarrowPeaks", FALSE, TRUE)
for (x in grep("Adult", names(MouseBrainATACNarrowPeaks), value = TRUE)) {
    rtracklayer::export(MouseBrainATACNarrowPeaks[[x]], con = sprintf("Report/Epigenome/MouseBrainATACNarrowPeaks/%s.bed.gz", x))
}

## 2.2.2 Mouse Brain DNase-seq
MouseBrainDNaseMetadata <- read.csv("Data/Database/ENCODE/DNase-seq/mm10/Brain/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
MouseBrainDNaseMetadata <- subset(MouseBrainDNaseMetadata, `Assembly` == "mm10" & `File Status` == "released")
MouseBrainDNaseInDir <- "Data/Database/ENCODE/DNase-seq/mm10/Brain"
MouseBrainDNaseTissueTypes <- list.dirs("Data/Database/ENCODE/DNase-seq/mm10/Brain", recursive = FALSE, full.names = FALSE)
MouseBrainDNaseAccToRep <- MouseBrainDNaseMetadata[, "Biological replicate(s)"]
names(MouseBrainDNaseAccToRep) <- MouseBrainDNaseMetadata[, "File accession"]
MouseBrainDNaseNarrowPeakAccs <- subset(MouseBrainDNaseMetadata, `File format` == "bigBed narrowPeak")[, "File accession"]
MouseBrainDNaseBroadPeakAccs <- subset(MouseBrainDNaseMetadata, `File format` == "bigBed broadPeak")[, "File accession"]
MouseBrainDNaseNarrowPeaks <- sapply(MouseBrainDNaseTissueTypes, function(tissue) { 
    accs <- list.files(sprintf("%s/%s", MouseBrainDNaseInDir, tissue), pattern = ".bigBed", recursive = FALSE, full.names = FALSE)
    accs <- sub("\\.bigBed", "", accs)
    accs <- accs[accs %in% MouseBrainDNaseNarrowPeakAccs]
    repNames <- paste0("Rep", MouseBrainDNaseAccToRep[accs])
    peaks <- sapply(accs, function(acc) { 
        message(tissue, " ", acc)
        rtracklayer::import(sprintf("%s/%s.narrowPeak", MouseBrainDNaseInDir, acc)) 
    }, simplify = FALSE)
    names(peaks) <- repNames
    peaks 
}, simplify = FALSE)
MouseBrainDNaseNarrowPeaks <- lapply(MouseBrainDNaseNarrowPeaks, function(X) {
    repNames <- names(X) 
    ind <- order(repNames) 
    X[ind] 
})
MouseBrainDNaseNarrowPeaks <- rapply(MouseBrainDNaseNarrowPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace")
MouseBrainDNaseNarrowPeaks <- lapply(MouseBrainDNaseNarrowPeaks, function(X) {
    X <- as(object = X, Class = "GRangesList") 
    reduce(Reduce(union, X), ignore.strand = TRUE) 
})
MouseBrainDNaseBroadPeaks <- sapply(MouseBrainDNaseTissueTypes, function(tissue) { 
    accs <- list.files(sprintf("%s/%s", MouseBrainDNaseInDir, tissue), pattern = ".bigBed", recursive = FALSE, full.names = FALSE)
    accs <- sub("\\.bigBed", "", accs)
    accs <- accs[accs %in% MouseBrainDNaseBroadPeakAccs]
    repNames <- paste0("Rep", MouseBrainDNaseAccToRep[accs])
    peaks <- sapply(accs, function(acc) { 
        message(tissue, " ", acc)
        rtracklayer::import(sprintf("%s/%s.broadPeak", MouseBrainDNaseInDir, acc)) 
    }, simplify = FALSE)
    names(peaks) <- repNames
    peaks 
}, simplify = FALSE)
MouseBrainDNaseBroadPeaks <- lapply(MouseBrainDNaseBroadPeaks, function(X) { 
    repNames <- names(X) 
    ind <- order(repNames) 
    X[ind] 
})
MouseBrainDNaseBroadPeaks <- rapply(MouseBrainDNaseBroadPeaks, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace")
MouseBrainDNaseBroadPeaks <- lapply(MouseBrainDNaseBroadPeaks, function(X) { 
    X <- as(object = X, Class = "GRangesList") 
    reduce(Reduce(union, X), ignore.strand = TRUE) 
})
write.csv(MouseBrainDNaseMetadata, file = "Report/Epigenome/MouseBrainDNaseMetadata.csv")
saveRDS(MouseBrainDNaseNarrowPeaks, file = "Data/Epigenome/MouseBrainDNaseNarrowPeaks.RDS")
saveRDS(MouseBrainDNaseBroadPeaks, file = "Data/Epigenome/MouseBrainDNaseBroadPeaks.RDS")

MouseBrainDNaseNarrowPeaks <- readRDS(file = "Data/Epigenome/MouseBrainDNaseNarrowPeaks.RDS")
dir.create("Report/Epigenome/MouseBrainDNaseNarrowPeaks", FALSE, TRUE)
for (x in grep("Adult", names(MouseBrainDNaseNarrowPeaks), value = TRUE)) {
    rtracklayer::export(MouseBrainDNaseNarrowPeaks[[x]], con = sprintf("Report/Epigenome/MouseBrainDNaseNarrowPeaks/%s.bed.gz", x))
}

## 2.2.3 Mouse Brain FAIRE-seq
## ENCODE does not have data.

## 2.2.4 Mouse Brain ChIP-seq
MouseBrainChIPMetadataTF <- read.csv("Data/Database/ENCODE/ChIP-seq/mm10/Brain/TF/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(MouseBrainChIPMetadataTF) <- MouseBrainChIPMetadataTF[, "File accession"]
MouseBrainChIPMetadataTF <- subset(MouseBrainChIPMetadataTF, `File Status` == "released" & Assembly == "mm10")
write.csv(MouseBrainChIPMetadataTF, file = "Data/Database/ENCODE/ChIP-seq/mm10/Brain/TF/metadata_cleaned.csv", row.names= FALSE)
MouseBrainChIPMetadataTF <- subset(MouseBrainChIPMetadataTF,  `File format` == "bigBed narrowPeak" & `Output type` == "optimal idr thresholded peaks")
MouseBrainChIPSampleIDsTF <- MouseBrainChIPMetadataTF[, "File accession"]
MouseBrainChIPTissuesTF <- MouseBrainChIPMetadataTF[, "Biosample term name"]
MouseBrainChIPRepNamesTF <- MouseBrainChIPMetadataTF[, "Biological replicate(s)"]
MouseBrainChIPTargetsTF <- sub("-mouse$", "", MouseBrainChIPMetadataTF[, "Experiment target"])
MouseBrainChIPMetadataTF <- data.frame(SampleID = MouseBrainChIPSampleIDsTF, Tissue = MouseBrainChIPTissuesTF, Target = MouseBrainChIPTargetsTF, stringsAsFactors = FALSE)
MouseBrainChIPMetadataTF
##      SampleID         Tissue Target
## 1 ENCFF113BDJ       midbrain   CTCF
## 2 ENCFF400OMH cortical plate POLR2A
## 3 ENCFF687MGE olfactory bulb   CTCF
## 4 ENCFF029GRM     cerebellum POLR2A
## 5 ENCFF923TYN      forebrain   CTCF
## 6 ENCFF868NVS      hindbrain   CTCF
## 7 ENCFF404ARG     cerebellum   CTCF
## 8 ENCFF548OFL olfactory bulb POLR2A
## 9 ENCFF779CYW cortical plate   CTCF
MouseBrainChIPTFs <- by(MouseBrainChIPMetadataTF, INDICES = MouseBrainChIPMetadataTF$Tissue, function(X) {
    targets <- X[, "Target"]
    accs <- X[, "SampleID"]
    peaks <- sapply(accs, function(acc) { 
        message(acc)
        Genome$standardizeSeqInfo(rtracklayer::import(sprintf("Data/Database/ENCODE/ChIP-seq/mm10/Brain/TF/%s.narrowPeak", acc)), seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE) 
    }, simplify = FALSE) 
    names(peaks) <- targets
    peaks
}, simplify = FALSE)
names(MouseBrainChIPTFs)
## [1] "cerebellum"     "cortical plate" "forebrain"      "hindbrain"     
## [5] "midbrain"       "olfactory bulb"
MouseBrainChIPTFs <- as(MouseBrainChIPTFs, "list")
names(MouseBrainChIPTFs) <- c("Cerebellum", "CorticalPlate", "ForeBrain", "HindBrain", "MidBrain", "OlfactoryBulb")
write.csv(MouseBrainChIPMetadataTF, file = "Report/Epigenome/MouseBrainChIPMetadataTF.csv")
saveRDS(MouseBrainChIPTFs, file = "Data/Epigenome/MouseBrainChIPTFs.RDS")
## We finally decided not to include TFs because of the lack of age information from the metadata sheet. 

MouseBrainChIPExptMetaBHM <- read.csv("Data/Database/ENCODE/ChIP-seq/mm10/Brain/BroadHistoneMark/Experiment Report 2019_1_21.tsv", sep = "\t", as.is = TRUE, check.names = FALSE, skip = 1)
rownames(MouseBrainChIPExptMetaBHM) <- MouseBrainChIPExptMetaBHM[, "Accession"]
MouseBrainChIPMetadataBHM <- read.csv("Data/Database/ENCODE/ChIP-seq/mm10/Brain/BroadHistoneMark/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(MouseBrainChIPMetadataBHM) <- MouseBrainChIPMetadataBHM[, "File accession"]
MouseBrainChIPMetadataBHM <- subset(MouseBrainChIPMetadataBHM, `File Status` == "released" & Assembly == "mm10")
## Because Metadata does not contain mouse age information, we need to refer to the external table ExptMeta, whose key is Accession (experiment ID).
MouseBrainChIPMetadataBHM <- data.frame(MouseBrainChIPMetadataBHM, MouseBrainChIPExptMetaBHM[MouseBrainChIPMetadataBHM[, "Experiment accession"], ], stringsAsFactors = FALSE, check.names = FALSE)
write.csv(MouseBrainChIPMetadataBHM, file = "Data/Database/ENCODE/ChIP-seq/mm10/Brain/BroadHistoneMark/metadata_cleaned.csv", row.names= FALSE)
table(MouseBrainChIPMetadataBHM[, "File format"])
## 
## bigBed narrowPeak            bigWig 
##               396               594
table(MouseBrainChIPMetadataBHM[, "Output type"])
## 
## fold change over control                    peaks         replicated peaks 
##                      297                      297                       99 
##           signal p-value 
##                      297
MouseBrainChIPMetadataBHM <- subset(MouseBrainChIPMetadataBHM, `Output type` == "replicated peaks")
dim(MouseBrainChIPMetadataBHM)
## [1] 99 80
MouseBrainChIPSampleIDsBHM <- MouseBrainChIPMetadataBHM[, "File accession"]
MouseBrainChIPTissuesBHM <- MouseBrainChIPMetadataBHM[, "Biosample term name"]
MouseBrainChIPRepNamesBHM <- MouseBrainChIPMetadataBHM[, "Biological replicate(s)"]
MouseBrainChIPTargetsBHM <- sub("-mouse$", "", MouseBrainChIPMetadataBHM[, "Experiment target"])
MouseBrainChIPBiosampleInfoBHM <- MouseBrainChIPMetadataBHM[, "Biosample summary"]
MouseBrainChIPBiosampleInfoBHM <- sapply(strsplit(gsub("\\(|\\)", "", MouseBrainChIPBiosampleInfoBHM), " "), function(x) { 
                                             sub("Days$", "days", paste0(tools::toTitleCase(x[-1]), collapse = "")) 
})
MouseBrainChIPMetadataBHM <- data.frame(SampleID = MouseBrainChIPSampleIDsBHM, Tissue = MouseBrainChIPTissuesBHM, Target = MouseBrainChIPTargetsBHM, BiosampleInfo = MouseBrainChIPBiosampleInfoBHM, stringsAsFactors = FALSE)
MouseBrainChIPMetadataBHM
##       SampleID         Tissue   Target
##      SampleID         Tissue   Target                BiosampleInfo
## 1 ENCFF597QLV      forebrain H3K27me3      ForebrainEmbryo10.5days
## 2 ENCFF787SFY olfactory bulb  H3K4me1 OlfactoryBulbMaleAdult8Weeks
## 3 ENCFF484RAH       midbrain H3K36me3       MidbrainEmbryo10.5days
## 4 ENCFF042GKY       midbrain  H3K9me3       MidbrainEmbryo12.5days
## ...
## 98 ENCFF910SGA  midbrain H3K27me3  MidbrainEmbryo11.5days
## 99 ENCFF893KVM  midbrain H3K27me3  MidbrainEmbryo12.5days
tapply(MouseBrainChIPMetadataBHM[, 1], INDEX = list(MouseBrainChIPMetadataBHM[, 2], MouseBrainChIPMetadataBHM[, 3]), FUN = length)
##                H3K27me3 H3K36me3 H3K4me1 H3K9me3
## cerebellum            1       NA       1      NA
## cortical plate       NA       NA       1      NA
## forebrain             8        8       8       8
## hindbrain             8        8       8       8
## midbrain              8        7       8       8
## olfactory bulb       NA       NA       1      NA
tapply(MouseBrainChIPMetadataBHM[, 1], INDEX = list(MouseBrainChIPMetadataBHM[, 4], MouseBrainChIPMetadataBHM[, 3]), FUN = length)
##                              H3K27me3 H3K36me3 H3K4me1 H3K9me3
## CerebellumMaleAdult8Weeks           1       NA       1      NA
## CorticalPlateMaleAdult8Weeks       NA       NA       1      NA
## ForebrainEmbryo10.5days             1        1       1       1
## ForebrainEmbryo11.5days             1        1       1       1
## ForebrainEmbryo12.5days             1        1       1       1
## ForebrainEmbryo13.5days             1        1       1       1
## ForebrainEmbryo14.5days             1        1       1       1
## ForebrainEmbryo15.5days             1        1       1       1
## ForebrainEmbryo16.5days             1        1       1       1
## ForebrainPostnatal0days             1        1       1       1
## HindbrainEmbryo10.5days             1        1       1       1
## HindbrainEmbryo11.5days             1        1       1       1
## HindbrainEmbryo12.5days             1        1       1       1
## HindbrainEmbryo13.5days             1        1       1       1
## HindbrainEmbryo14.5days             1        1       1       1
## HindbrainEmbryo15.5days             1        1       1       1
## HindbrainEmbryo16.5days             1        1       1       1
## HindbrainPostnatal0days             1        1       1       1
## MidbrainEmbryo10.5days              1        1       1       1
## MidbrainEmbryo11.5days              1        1       1       1
## MidbrainEmbryo12.5days              1       NA       1       1
## MidbrainEmbryo13.5days              1        1       1       1
## MidbrainEmbryo14.5days              1        1       1       1
## MidbrainEmbryo15.5days              1        1       1       1
## MidbrainEmbryo16.5days              1        1       1       1
## MidbrainPostnatal0days              1        1       1       1
## OlfactoryBulbMaleAdult8Weeks       NA       NA       1      NA
MouseBrainChIPBHMs <- by(MouseBrainChIPMetadataBHM, INDICES = MouseBrainChIPMetadataBHM$BiosampleInfo, function(X) {
    targets <- X[, "Target"]
    accs <- X[, "SampleID"]
    peaks <- sapply(accs, function(acc) { 
    message(acc, "...") 
    Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/ChIP-seq/mm10/Brain/BroadHistoneMark/%s.narrowPeak", acc)) 
    }, simplify = FALSE) 
    names(peaks) <- targets
    peaks
}, simplify = FALSE)
MouseBrainChIPBiosampleLabsBHM <- names(MouseBrainChIPBHMs)
MouseBrainChIPBiosampleLabsBHM
##  [1] "CerebellumMaleAdult8Weeks"    "CorticalPlateMaleAdult8Weeks"
##  [3] "ForebrainEmbryo10.5days"      "ForebrainEmbryo11.5days"     
##  [5] "ForebrainEmbryo12.5days"      "ForebrainEmbryo13.5days"     
##  [7] "ForebrainEmbryo14.5days"      "ForebrainEmbryo15.5days"     
##  [9] "ForebrainEmbryo16.5days"      "ForebrainPostnatal0days"     
## [11] "HindbrainEmbryo10.5days"      "HindbrainEmbryo11.5days"     
## [13] "HindbrainEmbryo12.5days"      "HindbrainEmbryo13.5days"     
## [15] "HindbrainEmbryo14.5days"      "HindbrainEmbryo15.5days"     
## [17] "HindbrainEmbryo16.5days"      "HindbrainPostnatal0days"     
## [19] "MidbrainEmbryo10.5days"       "MidbrainEmbryo11.5days"      
## [21] "MidbrainEmbryo12.5days"       "MidbrainEmbryo13.5days"      
## [23] "MidbrainEmbryo14.5days"       "MidbrainEmbryo15.5days"      
## [25] "MidbrainEmbryo16.5days"       "MidbrainPostnatal0days"      
## [27] "OlfactoryBulbMaleAdult8Weeks"
MouseBrainChIPBHMs <- as(MouseBrainChIPBHMs, "list")
names(MouseBrainChIPBHMs) <- MouseBrainChIPBiosampleLabsBHM
MouseBrainChIPBHMs <- rapply(MouseBrainChIPBHMs, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace")
write.csv(MouseBrainChIPMetadataBHM, file = "Report/Epigenome/MouseBrainChIPMetadataBHM.csv")
saveRDS(MouseBrainChIPBHMs, file = "Data/Epigenome/MouseBrainChIPBHMs.RDS")

MouseBrainChIPBHMs <- readRDS(file = "Data/Epigenome/MouseBrainChIPBHMs.RDS")
dir.create("Report/Epigenome/MouseBrainChIPBHMs", FALSE, TRUE)
for (x in grep("Adult", names(MouseBrainChIPBHMs), value = TRUE)) {
    for (y in names(MouseBrainChIPBHMs[[x]])) {
        rtracklayer::export(MouseBrainChIPBHMs[[x]][[y]], con = sprintf("Report/Epigenome/MouseBrainChIPBHMs/%s_%s.bed.gz", x, y))
    }
}

MouseBrainChIPExptMetaNHM <- read.csv("Data/Database/ENCODE/ChIP-seq/mm10/Brain/NarrowHistoneMark/Experiment Report 2019_1_21.tsv", sep = "\t", as.is = TRUE, check.names = FALSE, skip = 1)
rownames(MouseBrainChIPExptMetaNHM) <- MouseBrainChIPExptMetaNHM[, "Accession"]
MouseBrainChIPMetadataNHM <- read.csv("Data/Database/ENCODE/ChIP-seq/mm10/Brain/NarrowHistoneMark/metadata.tsv", sep = "\t", as.is = TRUE, check.names = FALSE)
rownames(MouseBrainChIPMetadataNHM) <- MouseBrainChIPMetadataNHM[, "File accession"]
MouseBrainChIPMetadataNHM <- subset(MouseBrainChIPMetadataNHM, `File Status` == "released" & Assembly == "mm10")
## Because Metadata does not contain mouse age information, we need to refer to the external table ExptMeta, whose key is Accession (experiment ID)
MouseBrainChIPMetadataNHM <- data.frame(MouseBrainChIPMetadataNHM, MouseBrainChIPExptMetaNHM[MouseBrainChIPMetadataNHM[, "Experiment accession"], ], stringsAsFactors = FALSE, check.names = FALSE)
write.csv(MouseBrainChIPMetadataNHM, file = "Data/Database/ENCODE/ChIP-seq/mm10/Brain/NarrowHistoneMark/metadata_cleaned.csv", row.names= FALSE)
table(MouseBrainChIPMetadataNHM[, "File format"])
## 
## bigBed narrowPeak            bigWig 
##               384               576 
table(MouseBrainChIPMetadataNHM[, "Output type"])
## 
## fold change over control                    peaks         replicated peaks 
##                      288                      288                       96 
##           signal p-value 
##                      288 
MouseBrainChIPMetadataNHM <- subset(MouseBrainChIPMetadataNHM, `Output type` == "replicated peaks")
dim(MouseBrainChIPMetadataNHM)
## [1] 96 80
MouseBrainChIPSampleIDsNHM <- MouseBrainChIPMetadataNHM[, "File accession"]
MouseBrainChIPTissuesNHM <- MouseBrainChIPMetadataNHM[, "Biosample term name"]
MouseBrainChIPRepNamesNHM <- MouseBrainChIPMetadataNHM[, "Biological replicate(s)"]
MouseBrainChIPTargetsNHM <- sub("-mouse$", "", MouseBrainChIPMetadataNHM[, "Experiment target"])
MouseBrainChIPBiosampleInfoNHM <- MouseBrainChIPMetadataNHM[, "Biosample summary"]
MouseBrainChIPBiosampleInfoNHM <- sapply(strsplit(gsub("\\(|\\)", "", MouseBrainChIPBiosampleInfoNHM), " "), function(x) { sub("Days$", "days", paste0(toTitleCase(x[-1]), collapse = "")) })
MouseBrainChIPMetadataNHM <- data.frame(SampleID = MouseBrainChIPSampleIDsNHM, Tissue = MouseBrainChIPTissuesNHM, Target = MouseBrainChIPTargetsNHM, BiosampleInfo = MouseBrainChIPBiosampleInfoNHM, stringsAsFactors = FALSE)
MouseBrainChIPMetadataNHM
##      SampleID         Tissue  Target                BiosampleInfo
## 1 ENCFF831XNV     cerebellum H3K4me3    CerebellumMaleAdult8Weeks
## 2 ENCFF169RNW      hindbrain H3K4me3      HindbrainPostnatal0days
## 3 ENCFF414QQR      hindbrain H3K27ac      HindbrainEmbryo13.5days
## 4 ENCFF887RGY      forebrain H3K27ac      ForebrainEmbryo10.5days
## ...
## 98 ENCFF910SGA  midbrain H3K27me3  MidbrainEmbryo11.5days
## 99 ENCFF893KVM  midbrain H3K27me3  MidbrainEmbryo12.5days
tapply(MouseBrainChIPMetadataNHM[, 1], INDEX = list(MouseBrainChIPMetadataNHM[, 2], MouseBrainChIPMetadataNHM[, 3]), FUN = length)
##                H3K27me3 H3K36me3 H3K4me1 H3K9me3
## cerebellum            1       NA       1      NA
## cortical plate       NA       NA       1      NA
## forebrain             8        8       8       8
## hindbrain             8        8       8       8
## midbrain              8        7       8       8
## olfactory bulb       NA       NA       1      NA
tapply(MouseBrainChIPMetadataNHM[, 1], INDEX = list(MouseBrainChIPMetadataNHM[, 4], MouseBrainChIPMetadataNHM[, 3]), FUN = length)
##                              H3K27ac H3K4me2 H3K4me3 H3K9ac
## CerebellumMaleAdult8Weeks          1      NA       1     NA
## CorticalPlateMaleAdult8Weeks       1      NA       1     NA
## ForebrainEmbryo10.5days            1      NA       1     NA
## ForebrainEmbryo11.5days            1       1       1      1
## ForebrainEmbryo12.5days            1       1       1      1
## ForebrainEmbryo13.5days            1       1       1      1
## ForebrainEmbryo14.5days            1       1       1      1
## ForebrainEmbryo15.5days            1       1       1      1
## ForebrainEmbryo16.5days            1       1       1      1
## ForebrainPostnatal0days            1       1       1      1
## HindbrainEmbryo10.5days            1      NA       1     NA
## HindbrainEmbryo11.5days            1       1       1      1
## HindbrainEmbryo12.5days            1       1       1      1
## HindbrainEmbryo13.5days            1       1       1      1
## HindbrainEmbryo14.5days            1       1       1      1
## HindbrainEmbryo15.5days            1       1       1      1
## HindbrainEmbryo16.5days            1       1       1      1
## HindbrainPostnatal0days            1       1       1      1
## MidbrainEmbryo10.5days             1      NA       1     NA
## MidbrainEmbryo11.5days             1       1       1      1
## MidbrainEmbryo12.5days             1       1       1      1
## MidbrainEmbryo13.5days             1       1       1      1
## MidbrainEmbryo14.5days             1       1       1      1
## MidbrainEmbryo15.5days             1       1       1      1
## MidbrainEmbryo16.5days             1       1       1      1
## MidbrainPostnatal0days             1       1       1      1
## OlfactoryBulbMaleAdult8Weeks       1      NA       1     NA
MouseBrainChIPNHMs <- by(MouseBrainChIPMetadataNHM, INDICES = MouseBrainChIPMetadataNHM$BiosampleInfo, function(X) {
    targets <- X[, "Target"]
    accs <- X[, "SampleID"]
    peaks <- sapply(accs, function(acc) { 
        message(acc, "...") 
        Genome$import.narrowPeak(sprintf("Data/Database/ENCODE/ChIP-seq/mm10/Brain/NarrowHistoneMark/%s.narrowPeak", acc)) 
    }, simplify = FALSE) 
    names(peaks) <- targets
    peaks
}, simplify = FALSE)
MouseBrainChIPBiosampleLabsNHM <- names(MouseBrainChIPNHMs)
MouseBrainChIPBiosampleLabsNHM
##  [1] "CerebellumMaleAdult8Weeks"    "CorticalPlateMaleAdult8Weeks"
##  [3] "ForebrainEmbryo10.5days"      "ForebrainEmbryo11.5days"     
##  [5] "ForebrainEmbryo12.5days"      "ForebrainEmbryo13.5days"     
##  [7] "ForebrainEmbryo14.5days"      "ForebrainEmbryo15.5days"     
##  [9] "ForebrainEmbryo16.5days"      "ForebrainPostnatal0days"     
## [11] "HindbrainEmbryo10.5days"      "HindbrainEmbryo11.5days"     
## [13] "HindbrainEmbryo12.5days"      "HindbrainEmbryo13.5days"     
## [15] "HindbrainEmbryo14.5days"      "HindbrainEmbryo15.5days"     
## [17] "HindbrainEmbryo16.5days"      "HindbrainPostnatal0days"     
## [19] "MidbrainEmbryo10.5days"       "MidbrainEmbryo11.5days"      
## [21] "MidbrainEmbryo12.5days"       "MidbrainEmbryo13.5days"      
## [23] "MidbrainEmbryo14.5days"       "MidbrainEmbryo15.5days"      
## [25] "MidbrainEmbryo16.5days"       "MidbrainPostnatal0days"      
## [27] "OlfactoryBulbMaleAdult8Weeks"
MouseBrainChIPNHMs <- as(MouseBrainChIPNHMs, "list")
names(MouseBrainChIPNHMs) <- MouseBrainChIPBiosampleLabsNHM
MouseBrainChIPNHMs <- rapply(MouseBrainChIPNHMs, Genome$standardizeSeqInfo, seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE, how = "replace")
write.csv(MouseBrainChIPMetadataNHM, file = "Report/Epigenome/MouseBrainChIPMetadataNHM.csv")
saveRDS(MouseBrainChIPNHMs, file = "Data/Epigenome/MouseBrainChIPNHMs.RDS")

MouseBrainChIPNHMs <- readRDS(file = "Data/Epigenome/MouseBrainChIPNHMs.RDS")
dir.create("Report/Epigenome/MouseBrainChIPNHMs", FALSE, TRUE)
for (x in grep("Adult", names(MouseBrainChIPNHMs), value = TRUE)) {
    for (y in names(MouseBrainChIPNHMs[[x]])) {
        rtracklayer::export(MouseBrainChIPNHMs[[x]][[y]], con = sprintf("Report/Epigenome/MouseBrainChIPNHMs/%s_%s.bed.gz", x, y))
    }
}
###########################################################################
## 3. non-B-form DNA database
###########################################################################
HumanNonBDNAclasses <- c("inverted_repeats", "direct_repeats", "g-quadruplex_forming_repeats", "mirror_repeats", "z-dna_motifs", "short_tandem_repeats", "a-phased_repeats")
HumanNonBDNAPeaks <- lapply(seq_along(HumanNonBDNAclasses), function(i) { 
    message(HumanNonBDNAclasses[i])
    Genome$standardizeSeqInfo(rtracklayer::import(sprintf("Data/nonBDB/2.0/hg38/%s.gtf.gz", HumanNonBDNAclasses[i])), seqInfo = Genome$MainSeqInfo$human, seqLevels = Genome$MainSeqLevels$human, prune = TRUE) 
})
names(HumanNonBDNAPeaks) <- c("Inverted", "Direct", "GQuadruplexForming", "Mirror", "ZDNAMotifs", "ShortTandem", "APhased")
HumanNonBDNAPeaks <- lapply(HumanNonBDNAPeaks, function(X) {
    mcols(X)[["sequence"]] <- NULL
    X 
})
HumanNonBDNAPeaks <- lapply(HumanNonBDNAPeaks, function(X) { 
    mcols(X)[["composition"]] <- NULL
    X 
})
saveRDS(HumanNonBDNAPeaks, file = "Data/Epigenome/HumanNonBDNAPeaks.RDS")
HumanNonBDNAPeaks <- lapply(HumanNonBDNAPeaks, function(X) {
    mcols(X)[["score"]] <- 0
    X 
})

dir.create("Report/Epigenome/HumanNonBDNAPeaks", FALSE, TRUE)
for (q in c("Inverted", "Direct", "GQuadruplexForming", "Mirror", "ZDNAMotifs", "ShortTandem", "APhased")) {
    rtracklayer::export(HumanNonBDNAPeaks[[q]], con = sprintf("Report/Epigenome/HumanNonBDNAPeaks/%s.bed.gz", q))
}

MouseNonBDNAclasses <- c("inverted_repeats", "direct_repeats", "g-quadruplex_forming_repeats", "mirror_repeats", "z-dna_motifs", "short_tandem_repeats", "a-phased_repeats")
MouseNonBDNAPeaks <- lapply(seq_along(MouseNonBDNAclasses), function(i) { 
    message(MouseNonBDNAclasses[i])
    Genome$standardizeSeqInfo(rtracklayer::import(sprintf("Data/nonBDB/2.0/mm10/%s.gtf.gz", MouseNonBDNAclasses[i])), seqInfo = Genome$MainSeqInfo$mouse, seqLevels = Genome$MainSeqLevels$mouse, prune = TRUE) 
})
names(MouseNonBDNAPeaks) <- c("Inverted", "Direct", "GQuadruplexForming", "Mirror", "ZDNAMotifs", "ShortTandem", "APhased")
MouseNonBDNAPeaks <- lapply(MouseNonBDNAPeaks, function(X) { 
    mcols(X)[["sequence"]] <- NULL
    X 
})
MouseNonBDNAPeaks <- lapply(MouseNonBDNAPeaks, function(X) {
    mcols(X)[["composition"]] <- NULL
    X 
})
saveRDS(MouseNonBDNAPeaks, file = "Data/Epigenome/MouseNonBDNAPeaks.RDS")
MouseNonBDNAPeaks <- lapply(MouseNonBDNAPeaks, function(X) { 
    mcols(X)[["score"]] <- 0
    X 
})

dir.create("Report/Epigenome/MouseNonBDNAPeaks", FALSE, TRUE)
for (q in c("Inverted", "Direct", "GQuadruplexForming", "Mirror", "ZDNAMotifs", "ShortTandem", "APhased")) {
    rtracklayer::export(MouseNonBDNAPeaks[[q]], con = sprintf("Report/Epigenome/MouseNonBDNAPeaks/%s.bed.gz", q))
}
