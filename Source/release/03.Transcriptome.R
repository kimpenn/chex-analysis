## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/release/functions.R")

.RERUN.K562GRORNA <- FALSE
.RERUN.K562RNARNA <- FALSE
.RERUN.K562scRNA <- FALSE
.RERUN.HumanBrainRNA <- FALSE
.RERUN.MouseBrainRNA <- FALSE
###########################################################################
## K562 transcriptome
###########################################################################
if (.RERUN.K562GRORNA) {
    K562GROExonCnts <- read.csv("Data/GSE60454/analyzed/Sample_GSM1480325/verseNoSplice_s1/Sample_GSM1480325.verse.exon.txt", row.names = 1, as.is = TRUE, sep = "\t")
    K562GROIntronCnts <- read.csv("Data/GSE60454/analyzed/Sample_GSM1480325/verseNoSplice_s1/Sample_GSM1480325.verse.intron-lev1-lev2.txt", row.names = 1, as.is = TRUE, sep = "\t")
    K562GROGeneCnts <- Stats$AddMatRows(K562GROExonCnts, K562GROIntronCnts)
    write.csv(K562GROExonCnts, file = gzfile("Report/release/Transcriptome/K562GROExonCnts.csv.gz"))
    write.csv(K562GROIntronCnts, file = gzfile("Report/release/Transcriptome/K562GROIntronCnts.csv.gz"))
    write.csv(K562GROGeneCnts, file = gzfile("Report/release/Transcriptome/K562GROGeneCnts.csv.gz"))
} else {
    K562GROExonCnts <- read.csv("Report/release/Transcriptome/K562GROExonCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
    K562GROIntronCnts <- read.csv("Report/release/Transcriptome/K562GROIntronCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
    K562GROGeneCnts <- read.csv("Report/release/Transcriptome/K562GROGeneCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
}
K562GROExonExprs <- log2(1 + K562GROExonCnts)
K562GROIntronExprs <- log2(1 + K562GROIntronCnts)
K562GROGeneExprs <- log2(1 + K562GROGeneCnts)
saveRDS(K562GROExonExprs, file = "Data/release/Transcriptome/K562GROExonExprs.RDS")
saveRDS(K562GROIntronExprs, file = "Data/release/Transcriptome/K562GROIntronExprs.RDS")
saveRDS(K562GROGeneExprs, file = "Data/release/Transcriptome/K562GROGeneExprs.RDS")

if (.RERUN.K562RNARNA) {
    K562RNAExonCnts <- read.csv("Data/GSE32213/analyzed/Sample_GSM798057/verse/Sample_GSM798057.verse.exon.cnts.txt", row.names = 1, as.is = TRUE, sep = "\t")
    K562RNAIntronCnts <- read.csv("Data/GSE32213/analyzed/Sample_GSM798057/verse/Sample_GSM798057.verse.intron-lev1-lev2.cnts.txt", row.names = 1, as.is = TRUE, sep = "\t")
    K562RNAGeneCnts <- Stats$AddMatRows(K562RNAExonCnts, K562RNAIntronCnts)
    write.csv(K562RNAExonCnts, file = gzfile("Report/release/Transcriptome/K562RNAExonCnts.csv.gz"))
    write.csv(K562RNAIntronCnts, file = gzfile("Report/release/Transcriptome/K562RNAIntronCnts.csv.gz"))
    write.csv(K562RNAGeneCnts, file = gzfile("Report/release/Transcriptome/K562RNAGeneCnts.csv.gz"))
} else {
    K562RNAExonCnts <- read.csv("Report/release/Transcriptome/K562RNAExonCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
    K562RNAIntronCnts <- read.csv("Report/release/Transcriptome/K562RNAIntronCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
    K562RNAGeneCnts <- read.csv("Report/release/Transcriptome/K562RNAGeneCnts.csv.gz", as.is = TRUE, check.names = FALSE, row.names = 1)
}
K562RNAExonExprs <- log2(1 + K562RNAExonCnts)
K562RNAIntronExprs <- log2(1 + K562RNAIntronCnts)
K562RNAGeneExprs <- log2(1 + K562RNAGeneCnts)
saveRDS(K562RNAExonExprs, file = "Data/release/Transcriptome/K562RNAExonExprs.RDS")
saveRDS(K562RNAIntronExprs, file = "Data/release/Transcriptome/K562RNAIntronExprs.RDS")
saveRDS(K562RNAGeneExprs, file = "Data/release/Transcriptome/K562RNAGeneExprs.RDS")

## K562 single-cell expression from Perturb-seq, but it seems the wild-type cells are not contained in the downloaded raw data (10X BAM file).
## so we have to use the GEO wild-type expression matrix, which contains exonic counts only. 
if (.RERUN.K562scRNA) {
    K562scRNAExonCnts <- read.csv("Data/GSE90063/GEO/GSE90063_k562_umi_wt.txt.gz", sep = "\t", as.is = TRUE, check.names = FALSE)
    dim(K562scRNAExonCnts)
    ## [1] 35635  5410
    rIDs <- K562scRNAExonCnts[, 1]
    geneSymbols <- sapply(strsplit(rIDs, "hg19_"), function(k) k[length(k)])
    sum(duplicated(geneSymbols))
    ## [1] 2697
    K562scRNAExonCnts <- as.matrix(K562scRNAExonCnts[, -1])
    K563scRNAExonCntsSummarized <- Stats$summarizeMatrix(K562scRNAExonCnts, by = geneSymbols, FUN = "mean", ncores = 12)
    K562scRNAExonCntsSummarized <- round(K562scRNAExonCntsSummarized)
    class(K562scRNAExonCntsSummarized) <- "integer"
    write.csv(K562scRNAExonCntsSummarized, file = gzfile("Report/release/Transcriptome/K562scRNAExonCntsSummarized.csv.gz"))
} else {
    K562scRNAExonCntsSummarized <- read.csv("Report/release/Transcriptome/K562scRNAExonCntsSummarized.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
}
K562scRNAGeneCntsSummarized <- K562scRNAExonCntsSummarized
K562scRNAExonExprs <- log2(1 + NGS$RLE2(K562scRNAExonCntsSummarized, threshold = 0.7))
K562scRNAGeneExprs <- K562scRNAExonExprs
saveRDS(K562scRNAExonExprs, file = "Data/release/Transcriptome/K562scRNAExonExprs.RDS")
saveRDS(K562scRNAGeneExprs, file = "Data/release/Transcriptome/K562scRNAGeneExprs.RDS")

###########################################################################
## Human brain in-house scRNA
###########################################################################
if (.RERUN.HumanBrainRNA) {
    HumanLabRNACnts <- read.csv("Data/SCTE/Processed/inc3/counts_exon_human.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
    HumanLabRNASampleIDs <- colnames(HumanLabRNACnts) 
    HumanLabRNAMetadata <- read.csv("Data/SCTE/Processed/inc3/metadata_statsFixed.csv", as.is = TRUE, check.names = FALSE)
    rownames(HumanLabRNAMetadata) <- HumanLabRNAMetadata[, "Sample Name"]
    HumanBrainRNAMetadata <- subset(HumanLabRNAMetadata, `Sample Name` %in% HumanLabRNASampleIDs & Region == "brain" & `Harvest Compound` != "TIVA" & `Harvest Method` == "pipette") #& VERSE.transcriptome == "/lab/repo/resources/verse/hg38.gencode25.gtf")
    HumanBrainRNASampleIDs <- intersect(HumanLabRNASampleIDs, HumanBrainRNAMetadata[, "Sample Name"])
    length(HumanBrainRNASampleIDs)
    ## [1] 372

    ## Hotfix for neuron and astrocyte labels
    HumanJennRNAMetadata <- read.csv("Data/release/Transcriptome/BrainSampleIDsFromJenn.csv", as.is = TRUE, check.names = FALSE)
    rownames(HumanJennRNAMetadata) <- HumanJennRNAMetadata[, "Sample"]
    HumanJennRNAMetadata <- subset(HumanJennRNAMetadata, `cell type` == "Neurons" | `cell type` == "Astrocytes")
    HumanRepoJennSampleIDs <- intersect(HumanBrainRNAMetadata[, "Sample Name"], HumanJennRNAMetadata[, "Sample"])
    HumanBrainRNAMetadata[HumanRepoJennSampleIDs, "Cell Class"] <- unname(c(Astrocytes = "astrocyte", Neurons = "neuron")[HumanJennRNAMetadata[HumanRepoJennSampleIDs, "cell type"]])
    write.csv(HumanBrainRNAMetadata, "Report/release/Transcriptome/HumanBrainRNAMetadata_CellTypeFixed.csv", row.names = FALSE)

    HumanBrainRNASampleIDs_failed_uniqMapPerc <- HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "STAR.Uniq-mapped Perc"] < 70 | is.na(HumanBrainRNAMetadata[, "STAR.Uniq-mapped Perc"]), "Sample Name"] 
    HumanBrainRNASampleIDs_failed_uniqMapReads <- HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] < 1e6 | is.na(HumanBrainRNAMetadata[, "STAR.Uniq-mapped Reads"]), "Sample Name"] 
    HumanBrainRNASampleIDs_failed_mismatchPerc <- HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "STAR.Mismatch-per-base Perc"] > 5 | is.na(HumanBrainRNAMetadata[, "STAR.Mismatch-per-base Perc"]), "Sample Name"] 
    HumanBrainRNASampleIDs_failed_spikeinPerc <- HumanBrainRNAMetadata[HumanBrainRNAMetadata[!is.na(HumanBrainRNAMetadata[, "VERSE.Perc: Spike-In"]), "VERSE.Perc: Spike-In"] > 10, "Sample Name"]
    HumanBrainRNASampleIDs_failed_mitoPerc <- na.omit(HumanBrainRNAMetadata[100 * HumanBrainRNAMetadata[, "VERSE.mito: Reads Counted"] / HumanBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] > 20, "Sample Name"])
    HumanBrainRNASampleIDs_failed_exonPerc <- na.omit(HumanBrainRNAMetadata[100 * (HumanBrainRNAMetadata[, "VERSE.exons Level 1,2: non-SpikeIn Reads Counted"] + HumanBrainRNAMetadata[, "VERSE.exons Level 3: Reads Counted"]) / HumanBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] < 40, "Sample Name"])
    HumanBrainRNASampleIDs_failed_antiRatio <- na.omit(HumanBrainRNAMetadata[HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "Strand Specific"], "VERSE.anti-exons: Reads Counted"] / (HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "Strand Specific"], "VERSE.exons Level 1,2: non-SpikeIn Reads Counted"] + HumanBrainRNAMetadata[HumanBrainRNAMetadata[, "Strand Specific"], "VERSE.exons Level 3: Reads Counted"]) > 2, "Sample Name"])

    HumanBrainRNASampleIDs_failed <- Reduce(union, list(HumanBrainRNASampleIDs_failed_uniqMapPerc, HumanBrainRNASampleIDs_failed_uniqMapReads, HumanBrainRNASampleIDs_failed_mismatchPerc, HumanBrainRNASampleIDs_failed_spikeinPerc, HumanBrainRNASampleIDs_failed_mitoPerc, HumanBrainRNASampleIDs_failed_exonPerc, HumanBrainRNASampleIDs_failed_antiRatio))
    length(HumanBrainRNASampleIDs_failed)
    ## [1] 150
    HumanBrainRNASampleIDs <- setdiff(HumanBrainRNASampleIDs, HumanBrainRNASampleIDs_failed)
    length(HumanBrainRNASampleIDs)
    ## [1] 222
    HumanBrainRNAMetadata <- HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ]

    ## We further separate neurons and astrocytes
    table(HumanBrainRNAMetadata[, "Cell Class"])
    ## 
    ## astrocyte    neuron 
    ##        29       193 
    ## after Jenn correction
    ## astrocyte    neuron 
    ##        92       130

    HumanAstroRNAMetadta <- subset(HumanBrainRNAMetadata, `Cell Class` == "astrocyte")
    HumanNeuronRNAMetadta <- subset(HumanBrainRNAMetadata, `Cell Class` == "neuron")
    HumanAstroRNASampleIDs <- HumanAstroRNAMetadta[, "Sample Name"]
    HumanNeuronRNASampleIDs <- HumanNeuronRNAMetadta[, "Sample Name"]
    grep("[a-z]$", sort(HumanNeuronRNASampleIDs), value = TRUE)
    ##  [1] "BRP114b" "BRP213b" "BRP214b" "BRP298b" "BRP299b" "BRP311b" "BRP312b"
    ##  [8] "BRP316b" "BRP316c" "BRP316d" "BRP317b" "BRP317c" "BRP317d" "BRP318b"
    ## [15] "BRP318c" "BRP318d" "BRP319b" "BRP319c" "BRP319d" "BRP321b" "BRP331b"
    ## [22] "BRP332b" "BRP386b" "BRP387b" "BRP396b" "BRP397b"
    HumanNeuronRNASampleIDs <- setdiff(HumanNeuronRNASampleIDs, grep("[a-z]$", sort(HumanNeuronRNASampleIDs), value = TRUE))
    length(HumanNeuronRNASampleIDs)
    ## [1] 104
    grep("[a-z]$", sort(HumanAstroRNASampleIDs), value = TRUE)
    ## [1] "BRP230b" "BRP231b"
    HumanAstroRNASampleIDs <- setdiff(HumanAstroRNASampleIDs, grep("[a-z]$", sort(HumanAstroRNASampleIDs), value = TRUE))
    length(HumanAstroRNASampleIDs)
    ## [1] 90

    ## We further remove technical replicates
    HumanBrainRNASampleIDs <- setdiff(HumanBrainRNASampleIDs, grep("[a-z]$", sort(HumanBrainRNASampleIDs), value = TRUE))
    length(HumanBrainRNASampleIDs)
    ## [1] 194

    HumanBrainRNAVerse <- sapply(HumanBrainRNASampleIDs, function(suid) { Repo$getVERSEcounts( EID = HumanBrainRNAMetadata[suid, "Experiment Name"], SID = HumanBrainRNAMetadata[suid, "Sample Name"], features = c("exon", "intron-lev1-lev2")) }, simplify = FALSE)

    sum(sapply(HumanBrainRNAVerse, function(x) ifelse(is.null(x[["exon"]]), TRUE, nrow(x[["exon"]]) == 0)))
    ## [1] 0
    HumanBrainRNAVerseGIDs <- list(Exon = sapply(HumanBrainRNAVerse, function(X) rownames(X[["exon"]])), Intron = sapply(HumanBrainRNAVerse, function(X) rownames(X[["intron-lev1-lev2"]])))
    table(sapply(HumanBrainRNAVerseGIDs[["Exon"]], length))
    ## 
    ## 50188 50526 
    ##     6   188 
    table(sapply(HumanBrainRNAVerseGIDs[["Intron"]], length))
    ## 
    ## 34449 34778 
    ##     6   188 
    names(which(sapply(HumanBrainRNAVerseGIDs[["Exon"]], length) == 50188))
    ## [1] "BRP230" "BRP231" "BRP213" "BRP214" "Brp403" "BRP411"
    names(which(sapply(HumanBrainRNAVerseGIDs[["Intron"]], length) == 34449))
    ## [1] "BRP230" "BRP231" "BRP213" "BRP214" "Brp403" "BRP411"

    ## We just exclude those samples whose exon number is 50188 and intron number 34449
    HumanBrainRNASampleIDs <- setdiff(HumanBrainRNASampleIDs, names(which(sapply(HumanBrainRNAVerseGIDs[["Exon"]], length) == 50188)))
    length(HumanBrainRNASampleIDs)
    ## [1] 188
    HumanBrainRNAVerse <- HumanBrainRNAVerse[HumanBrainRNASampleIDs]
    HumanBrainRNAMetadata <- HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ]
    write.csv(HumanBrainRNAMetadata, file = "Report/release/Transcriptome/HumanBrainRNAMetadata.csv")

    HumanAstroRNASampleIDs <- intersect(HumanAstroRNASampleIDs, HumanBrainRNASampleIDs)
    HumanNeuronRNASampleIDs <- intersect(HumanNeuronRNASampleIDs, HumanBrainRNASampleIDs)

    HumanBrainRNAGIDs <- list(Exon = sapply(HumanBrainRNAVerse, function(X) rownames(X[["exon"]])), Intron = sapply(HumanBrainRNAVerse, function(X) rownames(X[["intron-lev1-lev2"]])))
    dim(HumanBrainRNAGIDs[[1]])
    ## [1] 50526   188
    dim(HumanBrainRNAGIDs[[2]])
    ## [1] 34778   188
    HumanBrainRNAGIDs <- lapply(HumanBrainRNAGIDs, function(X) { K <- apply(X, 2, sort); apply(K, 1, unique) })

    HumanBrainRNAExonCnts <- sapply(HumanBrainRNAVerse, function(X) X[["exon"]][HumanBrainRNAGIDs$Exon, 1])
    rownames(HumanBrainRNAExonCnts) <- HumanBrainRNAGIDs[["Exon"]]
    HumanBrainRNAIntronCnts <- sapply(HumanBrainRNAVerse, function(X) X[["intron-lev1-lev2"]][HumanBrainRNAGIDs$Intron, 1])
    rownames(HumanBrainRNAIntronCnts) <- HumanBrainRNAGIDs[["Intron"]]

    write.csv(HumanBrainRNAExonCnts, file = gzfile("Report/release/Transcriptome/HumanBrainRNAExonCnts_withSpikein.csv.gz"))
    write.csv(HumanBrainRNAIntronCnts, file = gzfile("Report/release/Transcriptome/HumanBrainRNAIntronCnts.csv.gz"))
    HumanBrainRNAGIDs$ExonNoSpikein <- HumanBrainRNAGIDs[["Exon"]][!grepl("^spikeIn", HumanBrainRNAGIDs[["Exon"]])]
    lengths(HumanBrainRNAGIDs)
    ##          Exon        Intron ExonNoSpikein 
    ##         50526         34778         50427 
    HumanBrainRNAExonCnts <- HumanBrainRNAExonCnts[HumanBrainRNAGIDs[["ExonNoSpikein"]], ]
    write.csv(HumanBrainRNAExonCnts, file = gzfile("Report/release/Transcriptome/HumanBrainRNAExonCnts_noSpikein.csv.gz"))
} else {
    HumanBrainRNAExonCnts <- read.csv("Report/release/Transcriptome/HumanBrainRNAExonCnts_noSpikein.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
    HumanBrainRNAIntronCnts <- read.csv("Report/release/Transcriptome/HumanBrainRNAIntronCnts.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
}
HumanBrainRNAGeneCnts <- Stats$AddMatRows(HumanBrainRNAExonCnts, HumanBrainRNAIntronCnts)
write.csv(HumanBrainRNAGeneCnts, file = gzfile("Report/release/Transcriptome/HumanBrainRNAGeneCnts.csv.gz"))

if (.RERUN.HumanBrainRNA) {
    HumanBrainRNASampleIDs <- colnames(HumanBrainRNAExonCnts)
    HumanBrainRNAMetadata <- read.csv("Report/release/Transcriptome/HumanBrainRNAMetadata_CellTypeFixed.csv", as.is = TRUE, check.names = FALSE)
    rownames(HumanBrainRNAMetadata) <- HumanBrainRNAMetadata[, "Sample Name"]
    HumanBrainRNAMetadata <- HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ]
    HumanAstroRNASampleIDs <- subset(HumanBrainRNAMetadata, `Cell Class` == "astrocyte")[, "Sample Name"]
    HumanNeuronRNASampleIDs <- subset(HumanBrainRNAMetadata, `Cell Class` == "neuron")[, "Sample Name"]

    HumanBrainRNAMitoVerse <- sapply(HumanBrainRNASampleIDs, function(SID) { Repo$getVERSEcounts(exptName = paste0("E.", HumanBrainRNAMetadata[SID, "Experiment Name"]), sampleID = HumanBrainRNAMetadata[SID, "Sample Name"], features = c("mito")) }, simplify = FALSE)
    HumanBrainRNAMitoGIDs <- rownames(HumanBrainRNAMitoVerse[[1]][["mito"]])
    HumanBrainRNAMitoCnts <- sapply(HumanBrainRNAMitoVerse, function(Verse) Verse[["mito"]][HumanBrainRNAMitoGIDs, ])
    rownames(HumanBrainRNAMitoCnts) <- HumanBrainRNAMitoGIDs
    write.csv(HumanBrainRNAMitoCnts, file = "Report/release/Transcriptome/HumanBrainRNAMitoCnts.csv")
} else {
    HumanBrainRNAMitoCnts <- read.csv(file = "Report/release/Transcriptome/HumanBrainRNAMitoCnts.csv", row.names = 1, check.names = FALSE, as.is = TRUE)
}

###########################################################################
## Mouse brain in-house single-cell RNA-seq data
###########################################################################
## Astrocytes untreated (cortex, whole cell)
## 505	BRP564
## 505	BRP565
## 505	BRP566
## 505	BRP567
## 505	BRP568
## 505	BRP569
## 505	BRP570
## 505	BRP571
## 505	BRP572
## 505	BRP573 (unique mapping rate is 8%, unusable)

## Neurons untreated (cortex, soma)
## 487	BRP525 (very low "VERSE.Perc: Exons Level 1,2" 9.09%)
## 487	BRP526 (outlier suspect)
## 487	BRP527
## 487	BRP528
## 487	BRP529
## 487	BRP534
## 487	BRP535
## 487	BRP536
## 489	BRP537
## 489	BRP538
## 489	BRP539
## 489	BRP540
## 489	BRP541
## 489	BRP542
## 489	BRP543

if (.RERUN.MouseBrainRNA) {
    ## - Untreated astrocytes (2017-10-11)
    MouseAstroRNAExptIDs <- rep(505, 10)
    MouseAstroRNASampleIDs <- paste0("BRP", 564:573)

    ## - Untreated neurons
    MouseNeuronRNAExptIDs <- c(rep(487, 8), rep(489, 7))
    MouseNeuronRNASampleIDs <- paste0("BRP", c(525:529, 534:543))

    ## - Untreated all
    MouseBrainRNAExptIDs <- c(MouseAstroRNAExptIDs, MouseNeuronRNAExptIDs)
    MouseBrainRNASampleIDs <- c(MouseAstroRNASampleIDs, MouseNeuronRNASampleIDs)

    cat(MouseBrainRNASampleIDs, file = "Report/release/Transcriptome/MouseBrainRNASampleIDs.txt", sep = ", ")
    ## Download from GUI the metadata
    MouseBrainRNAMetadata <- read.csv("Report/release/Transcriptome/MouseBrainRNAMetadata.xls", as.is = TRUE, check.names = FALSE, sep = "\t")
    rownames(MouseBrainRNAMetadata) <- MouseBrainRNAMetadata[, "Sample Name"]

    MouseBrainRNASampleIDs_failed_uniqMapPerc <- MouseBrainRNAMetadata[MouseBrainRNAMetadata[, "STAR.Uniq-mapped Perc"] < 70 | is.na(MouseBrainRNAMetadata[, "STAR.Uniq-mapped Perc"]), "Sample Name"] 
    MouseBrainRNASampleIDs_failed_uniqMapReads <- MouseBrainRNAMetadata[MouseBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] < 1e7 | is.na(MouseBrainRNAMetadata[, "STAR.Uniq-mapped Reads"]), "Sample Name"] 
    MouseBrainRNASampleIDs_failed_mismatchPerc <- MouseBrainRNAMetadata[MouseBrainRNAMetadata[, "STAR.Perc Bases Mismatched"] > 5 | is.na(MouseBrainRNAMetadata[, "STAR.Perc Bases Mismatched"]), "Sample Name"] 
    MouseBrainRNASampleIDs_failed_spikeinPerc <- MouseBrainRNAMetadata[MouseBrainRNAMetadata[!is.na(MouseBrainRNAMetadata[, "VERSE.Perc: Spike-In"]), "VERSE.Perc: Spike-In"] > 10, "Sample Name"]
    MouseBrainRNASampleIDs_failed_mitoPerc <- na.omit(MouseBrainRNAMetadata[100 * MouseBrainRNAMetadata[, "VERSE.mito: Reads Counted"] / MouseBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] > 20, "Sample Name"])
    MouseBrainRNASampleIDs_failed_exonPerc <- na.omit(MouseBrainRNAMetadata[100 * MouseBrainRNAMetadata[, "VERSE.exons Level 1,2: non-SpikeIn Reads Counted"] / MouseBrainRNAMetadata[, "STAR.Uniq-mapped Reads"] < 40, "Sample Name"])

    MouseBrainRNASampleIDs_failed <- Reduce(union, list(MouseBrainRNASampleIDs_failed_uniqMapPerc, MouseBrainRNASampleIDs_failed_uniqMapReads, MouseBrainRNASampleIDs_failed_mismatchPerc, MouseBrainRNASampleIDs_failed_spikeinPerc, MouseBrainRNASampleIDs_failed_mitoPerc, MouseBrainRNASampleIDs_failed_exonPerc))
    length(MouseBrainRNASampleIDs_failed)
    ## [1] 5
    MouseBrainRNASampleIDs <- setdiff(MouseBrainRNASampleIDs, MouseBrainRNASampleIDs_failed)
    length(MouseBrainRNASampleIDs)
    ## [1] 20
    MouseBrainRNAMetadata <- MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ]

    ## Read the NGS statistics and GUI metadata
    MouseBrainRNAVerse <- sapply(MouseBrainRNASampleIDs, function(suid) { Repo$getVerseCounts( EID = MouseBrainRNAMetadata[suid, "Experiment Name"], SID = MouseBrainRNAMetadata[suid, "Sample Name"], features = c("exon", "intron")) }, simplify = FALSE)

    sum(sapply(MouseBrainRNAVerse, function(x) ifelse(is.null(x[["exon"]]), TRUE, nrow(x[["exon"]]) == 0)))
    ## [1] 0
    MouseBrainRNAVerseGIDs <- list(Exon = sapply(MouseBrainRNAVerse, function(X) rownames(X[["exon"]])), Intron = sapply(MouseBrainRNAVerse, function(X) rownames(X[["intron"]])))
    table(sapply(MouseBrainRNAVerseGIDs[["Exon"]], length))
    ##      1 
    ## 470900 
    table(sapply(MouseBrainRNAVerseGIDs[["Intron"]], length))
    ## 
    ##      1 
    ## 408440 

    write.csv(MouseBrainRNAMetadata, file = "Report/release/Transcriptome/MouseBrainRNAMetadata.csv")

    MouseAstroRNASampleIDs <- subset(MouseBrainRNAMetadata, `Cell Class` == "astrocyte")[, "Sample Name"]
    MouseNeuronRNASampleIDs <- subset(MouseBrainRNAMetadata, `Cell Class` == "neuron")[, "Sample Name"]
    MouseNeuronRNASampleIDs
    ##  [1] "BRP528" "BRP529" "BRP534" "BRP535" "BRP536" "BRP537" "BRP538" "BRP539"
    ##  [9] "BRP540" "BRP541" "BRP542" "BRP543"
    MouseAstroRNASampleIDs
    ## [1] "BRP564" "BRP566" "BRP567" "BRP568" "BRP569" "BRP570" "BRP571" "BRP572"

    MouseBrainRNAGIDs <- list(Exon = sapply(MouseBrainRNAVerse, function(X) rownames(X[["exon"]])), Intron = sapply(MouseBrainRNAVerse, function(X) rownames(X[["intron"]])))
    dim(MouseBrainRNAGIDs[[1]])
    ## [1] 23545   188
    dim(MouseBrainRNAGIDs[[2]])
    ## [1] 20422   188
    MouseBrainRNAGIDs <- lapply(MouseBrainRNAGIDs, function(X) { K <- apply(X, 2, sort); apply(K, 1, unique) })

    MouseBrainRNAExonCnts <- sapply(MouseBrainRNAVerse, function(X) X[["exon"]][MouseBrainRNAGIDs$Exon, 1])
    rownames(MouseBrainRNAExonCnts) <- MouseBrainRNAGIDs[["Exon"]]
    MouseBrainRNAIntronCnts <- sapply(MouseBrainRNAVerse, function(X) X[["intron"]][MouseBrainRNAGIDs$Intron, 1])
    rownames(MouseBrainRNAIntronCnts) <- MouseBrainRNAGIDs[["Intron"]]

    write.csv(MouseBrainRNAExonCnts, file = gzfile("Report/release/Transcriptome/MouseBrainRNAExonCnts_withSpikein.csv.gz"))
    write.csv(MouseBrainRNAIntronCnts, file = gzfile("Report/release/Transcriptome/MouseBrainRNAIntronCnts.csv.gz"))
    MouseBrainRNAGIDs$ExonNoSpikein <- MouseBrainRNAGIDs[["Exon"]][!grepl("^spikeIn", MouseBrainRNAGIDs[["Exon"]])]
    lengths(MouseBrainRNAGIDs)
    ##          Exon        Intron ExonNoSpikein 
    ##         23545         20422         23453
    MouseBrainRNAExonCnts <- MouseBrainRNAExonCnts[MouseBrainRNAGIDs[["ExonNoSpikein"]], ]
    write.csv(MouseBrainRNAExonCnts, file = gzfile("Report/release/Transcriptome/MouseBrainRNAExonCnts_noSpikein.csv.gz"))
} else {
    MouseBrainRNAExonCnts <- read.csv("Report/release/Transcriptome/MouseBrainRNAExonCnts_noSpikein.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
    MouseBrainRNAIntronCnts <- read.csv("Report/release/Transcriptome/MouseBrainRNAIntronCnts.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
}
MouseBrainRNAGeneCnts <- Stats$AddMatRows(MouseBrainRNAExonCnts, MouseBrainRNAIntronCnts)
write.csv(MouseBrainRNAGeneCnts, file = gzfile("Report/release/Transcriptome/MouseBrainRNAGeneCnts.csv.gz"))

if (.RERUN.MouseBrainRNA) {
    MouseBrainRNASampleIDs <- colnames(MouseBrainRNAExonCnts)
    MouseBrainRNAMetadata <- read.csv("Report/release/Transcriptome/MouseBrainRNAMetadata.csv", as.is = TRUE, check.names = FALSE)
    rownames(MouseBrainRNAMetadata) <- MouseBrainRNAMetadata[, "Sample Name"]
    MouseBrainRNAMetadata <- MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ]
    MouseAstroRNASampleIDs <- subset(MouseBrainRNAMetadata, `Cell Class` == "astrocyte")[, "Sample Name"]
    MouseNeuronRNASampleIDs <- subset(MouseBrainRNAMetadata, `Cell Class` == "neuron")[, "Sample Name"]

    MouseBrainRNAMitoVerse <- sapply(MouseBrainRNASampleIDs, function(SID) { Repo$getVERSEcounts(exptName = paste0("E.", MouseBrainRNAMetadata[SID, "Experiment Name"]), sampleID = MouseBrainRNAMetadata[SID, "Sample Name"], features = c("mito")) }, simplify = FALSE)
    MouseBrainRNAMitoGIDs <- rownames(MouseBrainRNAMitoVerse[[1]][["mito"]])
    MouseBrainRNAMitoCnts <- sapply(MouseBrainRNAMitoVerse, function(Verse) Verse[["mito"]][MouseBrainRNAMitoGIDs, ])
    rownames(MouseBrainRNAMitoCnts) <- MouseBrainRNAMitoGIDs
    write.csv(MouseBrainRNAMitoCnts, file = "Report/release/Transcriptome/MouseBrainRNAMitoCnts.csv")
} else {
    MouseBrainRNAMitoCnts <- read.csv(file = "Report/release/Transcriptome/MouseBrainRNAMitoCnts.csv", row.names = 1, check.names = FALSE, as.is = TRUE)
}


if (.RERUN.MouseInterneuronRNA) {
    ## Download from GUI the metadata
    MouseInterneuronRNAMetadata <- read.csv("Data/release/Transcriptome/MouseInterneuronRNASampleMetadata.xls", as.is = TRUE, check.names = FALSE, sep = "\t")
    rownames(MouseInterneuronRNAMetadata) <- MouseInterneuronRNAMetadata[, "Sample Name"]
    MouseInterneuronRNASampleIDs <- MouseInterneuronRNAMetadata[, "Sample Name"]
    MouseInterneuronRNAExptIDs <- MouseInterneuronRNAMetadata[, "Experiment Name"]
    cat(MouseInterneuronRNASampleIDs, file = "Data/release/Transcriptome/MouseInterneuronRNASampleIDs.txt", sep = ",")
    ## Because the old metadata lacks some important stats (e.g. mismatch rate per base), we need to run our procedure.
    MouseInterneuronRNAMetadata1 <- as.data.frame(t(sapply(MouseInterneuronRNASampleIDs, function(SID) { Repo$getNGSstats(exptID = MouseInterneuronRNAMetadata[SID, "Experiment Name"], sampleID = SID, features = c("exon", "intron", "mito", "intergenic", "lines_sines")) })), check.names = FALSE, stringsAsFactors = FALSE)
    MouseInterneuronRNAMetadata[["STAR.Perc Bases Mismatched"]] <- as.numeric(MouseInterneuronRNAMetadata1[["STAR.Mismatch-per-base Perc"]])

    MouseInterneuronRNASampleIDs_failed_uniqMapPerc <- MouseInterneuronRNAMetadata[MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Perc"] < 20 | is.na(MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Perc"]), "Sample Name"] 
    MouseInterneuronRNASampleIDs_failed_uniqMapReads <- MouseInterneuronRNAMetadata[MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Reads"] < 1e6 | is.na(MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Reads"]), "Sample Name"] 
    MouseInterneuronRNASampleIDs_failed_mismatchPerc <- MouseInterneuronRNAMetadata[MouseInterneuronRNAMetadata[, "STAR.Perc Bases Mismatched"] > 10 | is.na(MouseInterneuronRNAMetadata[, "STAR.Perc Bases Mismatched"]), "Sample Name"] 
    MouseInterneuronRNASampleIDs_failed_spikeinPerc <- MouseInterneuronRNAMetadata[MouseInterneuronRNAMetadata[!is.na(MouseInterneuronRNAMetadata[, "VERSE.Perc: Spike-In"]), "VERSE.Perc: Spike-In"] > 10, "Sample Name"]
    ## Because it is heart muscle cell, we would expect more mito counts; hence we relax the threshold to be 60
    MouseInterneuronRNASampleIDs_failed_mitoPerc <- na.omit(MouseInterneuronRNAMetadata[100 * MouseInterneuronRNAMetadata[, "VERSE.mito: Reads Counted"] / MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Reads"] > 60, "Sample Name"])
    MouseInterneuronRNASampleIDs_failed_exonPerc <- na.omit(MouseInterneuronRNAMetadata[100 * MouseInterneuronRNAMetadata[, "VERSE.exons Level 1,2: non-SpikeIn Reads Counted"] / MouseInterneuronRNAMetadata[, "STAR.Uniq-mapped Reads"] < 10, "Sample Name"])

    MouseInterneuronRNASampleIDs_failed <- Reduce(union, list(MouseInterneuronRNASampleIDs_failed_uniqMapPerc, MouseInterneuronRNASampleIDs_failed_uniqMapReads, MouseInterneuronRNASampleIDs_failed_mismatchPerc, MouseInterneuronRNASampleIDs_failed_spikeinPerc, MouseInterneuronRNASampleIDs_failed_mitoPerc, MouseInterneuronRNASampleIDs_failed_exonPerc))
    length(MouseInterneuronRNASampleIDs_failed)
    ## [1] 15
    cat(MouseInterneuronRNASampleIDs_failed, file = "Report/release/Transcriptome/MouseInterneuronRNASampleIDs_failed.txt", sep = "\n")
    MouseInterneuronRNASampleIDs <- setdiff(MouseInterneuronRNASampleIDs, MouseInterneuronRNASampleIDs_failed)
    length(MouseInterneuronRNASampleIDs)
    ## [1] 8
    MouseInterneuronRNAMetadata <- MouseInterneuronRNAMetadata[MouseInterneuronRNASampleIDs, ]
    write.csv(MouseInterneuronRNAMetadata, file = "Report/release/Transcriptome/MouseInterneuronRNAMetadata.csv")

    ## Read the NGS statistics and GUI metadata
    MouseInterneuronRNAVerse <- sapply(MouseInterneuronRNASampleIDs, function(suid) { Repo$getVERSEcounts(exptName = paste0("E.", MouseInterneuronRNAMetadata[suid, "Experiment Name"]), sampleID = MouseInterneuronRNAMetadata[suid, "Sample Name"], features = c("exon", "intron")) }, simplify = FALSE)

    sum(sapply(MouseInterneuronRNAVerse, function(x) ifelse(is.null(x[["exon"]]), TRUE, nrow(x[["exon"]]) == 0)))
    ## [1] 0
    MouseInterneuronRNAVerseGIDs <- list(Exon = sapply(MouseInterneuronRNAVerse, function(X) rownames(X[["exon"]]), simplify = FALSE), Intron = sapply(MouseInterneuronRNAVerse, function(X) rownames(X[["intron"]]), simplify = FALSE))
    table(sapply(MouseInterneuronRNAVerseGIDs[["Exon"]], length))
    ##  23545
    ##     101
    table(sapply(MouseInterneuronRNAVerseGIDs[["Intron"]], length))
    ## 
    ##  20422 
    ##     10 
    MouseInterneuronRNAVerseGIDs <- lapply(MouseInterneuronRNAVerseGIDs, function(X) do.call(cbind, X))
    MouseInterneuronRNAGIDs <- lapply(MouseInterneuronRNAVerseGIDs, function(X) { K <- apply(X, 2, sort); apply(K, 1, unique) })
    lengths(MouseInterneuronRNAGIDs)
    ##  Exon Intron 
    ## 23545  20422 

    MouseInterneuronRNAExonCnts <- sapply(MouseInterneuronRNAVerse, function(X) X[["exon"]][MouseInterneuronRNAGIDs$Exon, 1])
    rownames(MouseInterneuronRNAExonCnts) <- MouseInterneuronRNAGIDs[["Exon"]]
    MouseInterneuronRNAIntronCnts <- sapply(MouseInterneuronRNAVerse, function(X) X[["intron"]][MouseInterneuronRNAGIDs$Intron, 1])
    rownames(MouseInterneuronRNAIntronCnts) <- MouseInterneuronRNAGIDs[["Intron"]]

    write.csv(MouseInterneuronRNAExonCnts, file = gzfile("Report/release/Transcriptome/MouseInterneuronRNAExonCnts_withSpikein.csv.gz"))
    write.csv(MouseInterneuronRNAIntronCnts, file = gzfile("Report/release/Transcriptome/MouseInterneuronRNAIntronCnts.csv.gz"))
    MouseInterneuronRNAGIDs$ExonNoSpikein <- MouseInterneuronRNAGIDs[["Exon"]][!grepl("^spikeIn", MouseInterneuronRNAGIDs[["Exon"]])]
    lengths(MouseInterneuronRNAGIDs)
    ##          Exon        Intron ExonNoSpikein 
    ##         23545         20422         23453 
    MouseInterneuronRNAExonCnts <- MouseInterneuronRNAExonCnts[MouseInterneuronRNAGIDs[["ExonNoSpikein"]], ]
    write.csv(MouseInterneuronRNAExonCnts, file = gzfile("Report/release/Transcriptome/MouseInterneuronRNAExonCnts_noSpikein.csv.gz"))
} else {
    MouseInterneuronRNAExonCnts <- read.csv("Report/release/Transcriptome/MouseInterneuronRNAExonCnts_noSpikein.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
    MouseInterneuronRNAIntronCnts <- read.csv("Report/release/Transcriptome/MouseInterneuronRNAIntronCnts.csv.gz", row.names = 1, as.is = TRUE, check.names = FALSE)
}
MouseInterneuronRNAGeneCnts <- Stats$AddMatRows(MouseInterneuronRNAExonCnts, MouseInterneuronRNAIntronCnts)
write.csv(MouseInterneuronRNAGeneCnts, file = gzfile("Report/release/Transcriptome/MouseInterneuronRNAGeneCnts.csv.gz"))

if (.RERUN.MouseInterneuronRNA) {
    MouseInterneuronRNASampleIDs <- colnames(MouseInterneuronRNAExonCnts)
    MouseInterneuronRNAMetadata <- read.csv("Report/release/Transcriptome/MouseInterneuronRNAMetadata.csv", as.is = TRUE, check.names = FALSE)
    rownames(MouseInterneuronRNAMetadata) <- MouseInterneuronRNAMetadata[, "Sample Name"]
    MouseInterneuronRNAMetadata <- MouseInterneuronRNAMetadata[MouseInterneuronRNASampleIDs, ]
    table(MouseInterneuronRNAMetadata[["Cell Class"]])
    ## 
    ## neuron
    ##     10
    table(MouseInterneuronRNAMetadata[["Cell Type"]])
    ## fast spiking 
    ##           10

    MouseInterneuronRNAMitoVerse <- sapply(MouseInterneuronRNASampleIDs, function(SID) { Repo$getVERSEcounts(exptName = paste0("E.", MouseInterneuronRNAMetadata[SID, "Experiment Name"]), sampleID = MouseInterneuronRNAMetadata[SID, "Sample Name"], features = c("mito")) }, simplify = FALSE)
    MouseInterneuronRNAMitoGIDs <- rownames(MouseInterneuronRNAMitoVerse[[1]][["mito"]])
    MouseInterneuronRNAMitoCnts <- sapply(MouseInterneuronRNAMitoVerse, function(Verse) Verse[["mito"]][MouseInterneuronRNAMitoGIDs, ])
    rownames(MouseInterneuronRNAMitoCnts) <- MouseInterneuronRNAMitoGIDs
    write.csv(MouseInterneuronRNAMitoCnts, file = gzfile("Report/release/Transcriptome/MouseInterneuronRNAMitoCnts.csv.gz"))
} else {
    MouseInterneuronRNAMitoCnts <- read.csv(file = "Report/release/Transcriptome/MouseInterneuronRNAMitoCnts.csv.gz", row.names = 1, check.names = FALSE, as.is = TRUE)
}

###########################################################################
## in-house Brain RNA counts normalization and log2 transformation
###########################################################################
HumanBrainRNAExonExprs <- log2(1 + NGS$RLE2(HumanBrainRNAExonCnts, threshold = 0.7))
HumanBrainRNAIntronExprs <- log2(1 + NGS$RLE2(HumanBrainRNAIntronCnts, threshold = 0.7))
HumanBrainRNAGeneExprs <- log2(1 + NGS$RLE2(HumanBrainRNAGeneCnts, threshold = 0.7))
MouseBrainRNAExonExprs <- log2(1 + NGS$RLE2(MouseBrainRNAExonCnts, threshold = 0.7))
MouseBrainRNAIntronExprs <- log2(1 + NGS$RLE2(MouseBrainRNAIntronCnts, threshold = 0.7))
MouseBrainRNAGeneExprs <- log2(1 + NGS$RLE2(MouseBrainRNAGeneCnts, threshold = 0.7))

MouseInterneuronRNAExonExprs <- log2(1 + NGS$RLE2(MouseInterneuronRNAExonCnts, threshold = 0.7))
MouseInterneuronRNAIntronExprs <- log2(1 + NGS$RLE2(MouseInterneuronRNAIntronCnts, threshold = 0.7))
MouseInterneuronRNAGeneExprs <- log2(1 + NGS$RLE2(MouseInterneuronRNAGeneCnts, threshold = 0.7))

saveRDS(HumanBrainRNAExonExprs, file = "Data/release/Transcriptome/HumanBrainRNAExonExprs.RDS")
saveRDS(HumanBrainRNAIntronExprs, file = "Data/release/Transcriptome/HumanBrainRNAIntronExprs.RDS")
saveRDS(HumanBrainRNAGeneExprs, file = "Data/release/Transcriptome/HumanBrainRNAGeneExprs.RDS")
saveRDS(MouseBrainRNAExonExprs, file = "Data/release/Transcriptome/MouseBrainRNAExonExprs.RDS")
saveRDS(MouseBrainRNAIntronExprs, file = "Data/release/Transcriptome/MouseBrainRNAIntronExprs.RDS")
saveRDS(MouseBrainRNAGeneExprs, file = "Data/release/Transcriptome/MouseBrainRNAGeneExprs.RDS")

saveRDS(MouseInterneuronRNAExonExprs, file = "Data/release/Transcriptome/MouseInterneuronRNAExonExprs.RDS")
saveRDS(MouseInterneuronRNAIntronExprs, file = "Data/release/Transcriptome/MouseInterneuronRNAIntronExprs.RDS")
saveRDS(MouseInterneuronRNAGeneExprs, file = "Data/release/Transcriptome/MouseInterneuronRNAGeneExprs.RDS")

HumanAstroRNAExonExprs <- HumanBrainRNAExonExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
HumanNeuronRNAExonExprs <- HumanBrainRNAExonExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]
MouseAstroRNAExonExprs <- MouseBrainRNAExonExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
MouseNeuronRNAExonExprs <- MouseBrainRNAExonExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]
HumanAstroRNAIntronExprs <- HumanBrainRNAIntronExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
HumanNeuronRNAIntronExprs <- HumanBrainRNAIntronExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]
MouseAstroRNAIntronExprs <- MouseBrainRNAIntronExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
MouseNeuronRNAIntronExprs <- MouseBrainRNAIntronExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]
HumanAstroRNAGeneExprs <- HumanBrainRNAGeneExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
HumanNeuronRNAGeneExprs <- HumanBrainRNAGeneExprs[,  subset(HumanBrainRNAMetadata[HumanBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]
MouseAstroRNAGeneExprs <- MouseBrainRNAGeneExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "astrocyte")[, "Sample Name"]]
MouseNeuronRNAGeneExprs <- MouseBrainRNAGeneExprs[,  subset(MouseBrainRNAMetadata[MouseBrainRNASampleIDs, ], `Cell Class` == "neuron")[, "Sample Name"]]

saveRDS(HumanAstroRNAExonExprs, file = "Data/release/Transcriptome/HumanAstroRNAExonExprs.RDS")
saveRDS(HumanNeuronRNAExonExprs, file = "Data/release/Transcriptome/HumanNeuronRNAExonExprs.RDS")
saveRDS(MouseAstroRNAExonExprs, file = "Data/release/Transcriptome/MouseAstroRNAExonExprs.RDS")
saveRDS(MouseNeuronRNAExonExprs, file = "Data/release/Transcriptome/MouseNeuronRNAExonExprs.RDS")
saveRDS(HumanAstroRNAIntronExprs, file = "Data/release/Transcriptome/HumanAstroRNAIntronExprs.RDS")
saveRDS(HumanNeuronRNAIntronExprs, file = "Data/release/Transcriptome/HumanNeuronRNAIntronExprs.RDS")
saveRDS(MouseAstroRNAIntronExprs, file = "Data/release/Transcriptome/MouseAstroRNAIntronExprs.RDS")
saveRDS(MouseNeuronRNAIntronExprs, file = "Data/release/Transcriptome/MouseNeuronRNAIntronExprs.RDS")
saveRDS(HumanAstroRNAGeneExprs, file = "Data/release/Transcriptome/HumanAstroRNAGeneExprs.RDS")
saveRDS(HumanNeuronRNAGeneExprs, file = "Data/release/Transcriptome/HumanNeuronRNAGeneExprs.RDS")
saveRDS(MouseAstroRNAGeneExprs, file = "Data/release/Transcriptome/MouseAstroRNAGeneExprs.RDS")
saveRDS(MouseNeuronRNAGeneExprs, file = "Data/release/Transcriptome/MouseNeuronRNAGeneExprs.RDS")
