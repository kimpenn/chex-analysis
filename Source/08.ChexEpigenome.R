## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/functions.R")

SampleInfoFull <- read.csv("Data/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfoFull) <- SampleInfoFull[["SampleID"]]
SampleInfoFull <- subset(SampleInfoFull, IsOut == "N")
sampleIDsFull <- SampleInfoFull[, "SampleID"]
sampleIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual")[, "SampleID"]
sampleIDs <- subset(SampleInfoFull, CompType == "Biol")[, "SampleID"]
SampleInfo <- SampleInfoFull[sampleIDs, ]

bioGroups <- c(
    "K562", "K562TPAnone", "K562TPA15min", "K562TPA1hr", "K562TPA2hr", "K562TPA24hr", 
    "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture",
    "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice"
)
sampleIDsVirtmaxByBioGroup <- sapply(bioGroups, function(bioGroup) tail(subset(SampleInfoFull, (IsNegCtrl == "N" | BioGroup == "NoCell" | BioGroup == "K562MungBean" | BioGroup == "HBR") & ProbeType == "Positive" & BioGroup == bioGroup), 1)[["SampleID"]])

GRs_human <- GRanges(seqnames = seqnames(Genome$MainSeqInfo$human), IRanges(start = 1, end = seqlengths(Genome$MainSeqInfo$human)), seqinfo = Genome$MainSeqInfo$human)
GRs_mouse <- GRanges(seqnames = seqnames(Genome$MainSeqInfo$mouse), IRanges(start = 1, end = seqlengths(Genome$MainSeqInfo$mouse)), seqinfo = Genome$MainSeqInfo$mouse)
GRs_blacklisted_human <- Genome$standardizeSeqInfo(rtracklayer::import("Data/ENCODE/Blacklist/Kundaje/hg38.blacklist.bed"), seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
GRs_blacklisted_mouse <- Genome$standardizeSeqInfo(rtracklayer::import("Data/ENCODE/Blacklist/Kundaje/mm10.blacklist.bed"), seqInfo = Genome$MainSeqInfo$mouse, prune = TRUE, seqLevels = Genome$MainSeqLevels$mouse)

GRsNoBlacklisted_human <- setdiff(GRs_human, GRs_blacklisted_human)
GRsNoBlacklisted_mouse <- setdiff(GRs_mouse, GRs_blacklisted_mouse)
rtracklayer::export(GRsNoBlacklisted_human, con = "Report/AlignmentArtifacts/GRsNoBlacklisted_human.bed")
rtracklayer::export(GRsNoBlacklisted_mouse, con = "Report/AlignmentArtifacts/GRsNoBlacklisted_mouse.bed")

qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
GRsABreadCmate5End <- readRDS(sprintf("Data/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh))
GRsABreadCmate5EndVirtual <- sapply(sampleIDsVirtual, function(sampleID) {
    message(sampleID)
    sourceIDs <- SampleInfoFull[sampleID, "SourceIDsNoOut"]
    SIDs <- strsplit(sourceIDs, ",")[[1]]
    grs <- GRsABreadCmate5End[SIDs]
    Reduce(c, grs)
}, simplify = FALSE)
GRsABreadCmate5EndFull <- c(GRsABreadCmate5End, GRsABreadCmate5EndVirtual)[sampleIDsFull]
GRsABreadCmate5EndExtFull <- lapply(GRsABreadCmate5EndFull, function(GRs) flank(GRs, width = 1000, start = TRUE, both = TRUE))

bioTrtGroupsVirtmax <- c(
    "HumanAstroCulturePositiveAll", "HumanNeuronCulturePositiveMerged", "HumanInterneuronCulturePositiveMerged", 
    "MouseAstroCulturePositiveAll", "MouseNeuronCulturePositiveAll", "MouseNeuronSlicePositiveAll", "MouseInterneuronSlicePositiveMerged"
)
bioTrtGroups <- sub("Positive.*$", "", bioTrtGroupsVirtmax)
names(bioTrtGroupsVirtmax) <- bioTrtGroups

GRsABreadCmate5EndExtFull_K562PositiveAll_reduced <- reduce(GRsABreadCmate5EndExtFull$K562PositiveAll, ignore.strand=TRUE)
mcols(GRsABreadCmate5EndExtFull_K562PositiveAll_reduced)[["name"]] <- NULL
GRsABreadCmate5EndExtFull_K562PositiveAll_reduced <- trim(GRsABreadCmate5EndExtFull_K562PositiveAll_reduced)
start(GRsABreadCmate5EndExtFull_K562PositiveAll_reduced)[start(GRsABreadCmate5EndExtFull_K562PositiveAll_reduced) <= 0] <- 1
rtracklayer::export(GRsABreadCmate5EndExtFull_K562PositiveAll_reduced, con = sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/GRsABreadCmate5EndExtFull_K562PositiveAll_reduced.bed", qualOutInPair, mapqTh))
 
GRsABreadCmate5EndExtFull_virtmax_reduced <- lapply(bioTrtGroupsVirtmax, function(x) reduce(GRsABreadCmate5EndExtFull[[x]], ignore.strand=TRUE))

CvgsABreadCmate5EndFull <- lapply(GRsABreadCmate5EndFull, function(GRs) coverage(GRs))
CvgsABreadCmate5EndExtFull <- lapply(GRsABreadCmate5EndExtFull, function(GRs) coverage(GRs))

###########################################################################
## Load K562 epigenomes
###########################################################################
K562ATACNarrowPeaks <- readRDS("Data/Epigenome/K562ATACNarrowPeaks.RDS")
K562ATACBroadPeaks <- readRDS("Data/Epigenome/K562ATACBroadPeaks.RDS")
K562DNaseNarrowPeaks <- readRDS("Data/Epigenome/K562DNaseNarrowPeaks.RDS")
K562DNaseBroadPeaks <- readRDS("Data/Epigenome/K562DNaseBroadPeaks.RDS")
K562FAIRENarrowPeaks <- readRDS("Data/Epigenome/K562FAIRENarrowPeaks.RDS")
K562RRBSMethyLevels <- readRDS("Data/Epigenome/K562RRBSMethyLevels.RDS")
K562RepOrigPeaks <- readRDS("Data/Epigenome/K562RepOrigPeaks.RDS")
K562ChIPTFs <- readRDS("Data/Epigenome/K562ChIPTFs.RDS")
K562ChIPNHMs <- readRDS("Data/Epigenome/K562ChIPNHMs.RDS")
K562ChIPBHMs <- readRDS("Data/Epigenome/K562ChIPBHMs.RDS")
K562SEPeaks <- readRDS("Data/Epigenome/K562SEPeaks.RDS")
K562GROPeaks <- readRDS("Data/Epigenome/K562GROPeaks.RDS")
K562RloopPeaks <- readRDS("Data/Epigenome/K562RloopPeaks.RDS")
K562NonBDNAPeaks <- readRDS("Data/Epigenome/HumanNonBDNAPeaks.RDS")

K562coreChIPTFs <- c("CTCF", "MYC", "POLR2A", "POLR2B", "GATA1", "GATA2")
K562coreChIPNHMs <- c("H3K27ac", "H3K4me3", "H3K9ac", "H3K4me2", "H2AFZ")
K562coreChIPBHMs <- c("H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3", "H3K79me2", "H3K9me1")

###########################################################################
## Prepare Cvgs for epigenomes
###########################################################################
CvgsK562ATACNarrowPeak <- coverage(K562ATACNarrowPeaks[["K562UntreatedAll"]], weight = 1L)
CvgsK562ATACBroadPeak <- coverage(K562ATACBroadPeaks[["K562UntreatedAll"]], weight = 1L)
CvgsK562DNaseNarrowPeak <- lapply(K562DNaseNarrowPeaks, coverage, weight = 1L)
CvgsK562DNaseBroadPeak <- lapply(K562DNaseBroadPeaks, coverage, weight = 1L)
CvgsK562DNaseNarrowPeak <- CvgsK562DNaseNarrowPeak[["Rep1"]] + CvgsK562DNaseNarrowPeak[["Rep2"]]
CvgsK562DNaseBroadPeak <- CvgsK562DNaseBroadPeak[["Rep1"]] + CvgsK562DNaseBroadPeak[["Rep2"]]
CvgsK562FAIRENarrowPeak <- coverage(K562FAIRENarrowPeaks[[1]], weight = 1L)

CvgsK562RRBSMethyLevels <- lapply(K562RRBSMethyLevels, coverage, weight = 1L)
names(CvgsK562RRBSMethyLevels) <- paste0("Methyl", names(CvgsK562RRBSMethyLevels))
CvgsK562RepOrigPeak <- coverage(K562RepOrigPeaks[[1]], weight = 1L)
CvgsK562ChIPTFs <- lapply(K562ChIPTFs, coverage, weight = 1L)
CvgsK562ChIPNHMs <- lapply(K562ChIPNHMs, coverage, weight = 1L)
CvgsK562ChIPBHMs <- lapply(K562ChIPBHMs, coverage, weight = 1L)
CvgsK562SEPeak <- coverage(K562SEPeaks, weight = 1L)
CvgsK562GROPeak <- coverage(K562GROPeaks[["Homer"]], weight = 1L)
 
CvgsK562EpiCore <- c(
    list(ATAC = CvgsK562ATACNarrowPeak), 
    list(DNase = CvgsK562DNaseNarrowPeak), 
    list(FAIRE = CvgsK562FAIRENarrowPeak),
    list(GRO = CvgsK562GROPeak), 
    list(SE = CvgsK562SEPeak), 
    list(RepOrig = CvgsK562RepOrigPeak), 
    CvgsK562RRBSMethyLevels, 
    CvgsK562ChIPNHMs, 
    CvgsK562ChIPBHMs, 
    CvgsK562ChIPTFs[grep("POL", names(K562ChIPTFs))]
    )

CvgsK562EpiExt <- c(
    list(ATAC = CvgsK562ATACNarrowPeak), 
    list(DNase = CvgsK562DNaseNarrowPeak), 
    list(FAIRE = CvgsK562FAIRENarrowPeak),
    list(GRO = CvgsK562GROPeak), 
    list(SE = CvgsK562SEPeak), 
    list(RepOrig = CvgsK562RepOrigPeak), 
    CvgsK562RRBSMethyLevels, 
    CvgsK562ChIPNHMs, 
    CvgsK562ChIPBHMs, 
    CvgsK562ChIPTFs
    )

CvgsK562ChexEpiCore <- c(list(CHEX = CvgsABreadCmate5EndExtFull[["K562PositiveAll"]]), CvgsK562EpiCore)
CvgsK562ChexEpiExt <- c(list(CHEX = CvgsABreadCmate5EndExtFull[["K562PositiveAll"]]), CvgsK562EpiExt)
assaysCore <- names(CvgsK562ChexEpiCore) 
assaysExt <- names(CvgsK562ChexEpiExt)

###########################################################################
## Use chromosomal bins as features
###########################################################################
binsizes <- c(2e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6)
binsizeLabs <- c("2k", "5k", "10k", "50k", "100k", "500k", "1m", "5m")
names(binsizes) <- binsizeLabs
BinsListHuman <- sapply(binsizes, function(binsize) Genome$getTilesFromSeqLengths(Genome$MainSeqLengths$human, seqInfo = Genome$MainSeqInfo$human, width = binsize), simplify = FALSE)
MainSeqLevelsNoMY <- lapply(Genome$MainSeqLevels, function(seqLevels) seqLevels[!seqLevels %in% c("chrM", "chrY")])
BinsListHumanNoMY <- lapply(BinsListHuman, function(bins) { seqLevels <- MainSeqLevelsNoMY[["human"]]; bins[seqLevels] })
BinLabsListHumanNoMY <- lapply(BinsListHumanNoMY, function(BinsList) lapply(BinsList, function(Bin) as.character(Bin)))

PRsListK562ChexEpiCore <- sapply(binsizeLabs, function(binsizeLab) lapply(CvgsK562ChexEpiCore, function(Cvgs) { viewSums(Views(Cvgs[MainSeqLevelsNoMY$human], BinsListHumanNoMY[[binsizeLab]])) }), simplify = FALSE)
PRsListK562ChexEpiCore <- sapply(binsizeLabs, function(binsizeLab) sapply(assaysCore, function(assay) { structure(unlist(PRsListK562ChexEpiCore[[binsizeLab]][[assay]]), names = unlist(BinLabsListHumanNoMY[[binsizeLab]])) }), simplify = FALSE)
bPRsListK562ChexEpiCore <- lapply(PRsListK562ChexEpiCore, function(X) { X <- X > 0; class(X) <- "integer"; X })

bPRsListK562ChexEpiCoreJac <- sapply(binsizeLabs, function(binsizeLab) sapply(assaysCore, function(j) sapply(assaysCore, function(i) { message(binsizeLab, " ", j, " ", i); x <- bPRsListK562ChexEpiCore[[binsizeLab]][, i]; y <- bPRsListK562ChexEpiCore[[binsizeLab]][, j]; z <- x + y; sum(z == 2)/sum(z == 1 | z == 2) })), simplify = FALSE)

dirname <- sprintf("Report/ChexEpigenome/AssaysOverlap/HclustHeatmap/%s_%s", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)
for (binsizeLab in binsizeLabs) {
    write.csv(bPRsListK562ChexEpiCoreJac[[binsizeLab]], file = sprintf("%s/bPRsListK562ChexEpiCoreJacSim_%s.csv", dirname, binsizeLab))
}

bPRsListK562ChexEpiCoreJac <- sapply(binsizeLabs, function(binsizeLab) read.csv(sprintf("%s/bPRsListK562ChexEpiCoreJacSim_%s.csv", dirname, binsizeLab), as.is = TRUE, check.names = FALSE, row.names = 1), simplify = FALSE)
bPRsListK562ChexEpiCoreJacHclusts <- lapply(bPRsListK562ChexEpiCoreJac, function(sim) { ape::as.phylo(hclust(as.dist(1 - sim))) })

pal <- brewer.pal(4, "Set1")
names(pal) <- c("ATAC", "DNase", "FAIRE", "CHEX")
colorsCore <- pal[assaysCore]
colorsCore[is.na(colorsCore)] <- "black"
names(colorsCore)[is.na(names(colorsCore))] <- setdiff(assaysCore, c("ATAC", "DNase", "FAIRE", "CHEX"))
colorsCore <- colorsCore[assaysCore]

colorsExt <- pal[assaysExt]
colorsExt[is.na(colorsExt)] <- "gray"
names(colorsExt)[is.na(names(colorsExt))] <- setdiff(assaysExt, c("ATAC", "DNase", "FAIRE", "CHEX"))
colorsExt <- colorsExt[assaysExt]
colorsExt[assaysCore] <- "black"
colorsExt[c("ATAC", "DNase", "FAIRE", "CHEX")] <- pal
cexExt <- ifelse(assaysExt %in% assaysCore, 1, 0.6)
names(cexExt) <- assaysExt

dirname <- sprintf("Report/ChexEpigenome/AssaysOverlap/HclustHeatmap/%s_%s", qualOutInPair, mapqTh)
pdf(sprintf("%s/bPRsListK562ChexEpiCoreHclust_Jaccard.pdf", dirname), height = 7, width = 5)
for (binsizeLab in binsizeLabs) {
    hc <- bPRsListK562ChexEpiCoreJacHclusts[[binsizeLab]]
    plot(hc, tip.color = colorsCore, main = paste0("binsize = ", binsizeLab))
    axisPhylo(side = 1)
}
dev.off()

###########################################################################
## Load extended K562 epigenomes
###########################################################################
K562ATACNarrowPeaks <- readRDS("Data/Epigenome/K562ATACNarrowPeaks.RDS")
K562DNaseNarrowPeaks <- readRDS("Data/Epigenome/K562DNaseNarrowPeaks.RDS")
K562FAIRENarrowPeaks <- readRDS("Data/Epigenome/K562FAIRENarrowPeaks.RDS")
K562RRBSMethyLevels <- readRDS("Data/Epigenome/K562RRBSMethyLevels.RDS")
K562RepOrigPeaks <- readRDS("Data/Epigenome/K562RepOrigPeaks.RDS")
K562SEPeaks <- readRDS("Data/Epigenome/K562SEPeaks.RDS")
K562GROPeaks <- readRDS("Data/Epigenome/K562GROPeaks.RDS")
K562RloopPeaks <- readRDS("Data/Epigenome/K562RloopPeaks.RDS")
K562NonBDNAPeaks <- readRDS("Data/Epigenome/K562NonBDNAPeaks.RDS")
K562ChIPNHMs <- readRDS("Data/Epigenome/K562ChIPNHMs.RDS")
K562ChIPBHMs <- readRDS("Data/Epigenome/K562ChIPBHMs.RDS")
K562ChIPTFs <- readRDS("Data/Epigenome/K562ChIPTFs.RDS")
K562Epigenomes <- list(
    ATAC = K562ATACNarrowPeaks$K562UntreatedAll, 
    DNase = K562DNaseNarrowPeaks$Rep1, 
    FAIRE = K562FAIRENarrowPeaks[[1]], 
    MethyLow = K562RRBSMethyLevels[["Low"]], 
    MethyMidLow = K562RRBSMethyLevels[["MidLow"]], 
    MethyMidHigh = K562RRBSMethyLevels[["MidHigh"]], 
    MethyHigh = K562RRBSMethyLevels[["High"]], 
    GRO = K562GROPeaks[[1]],
    SE = K562SEPeaks, 
    RepOrig = K562RepOrigPeaks[[1]], 
    Rloop = K562RloopPeaks[["GSM1720619"]],
    NonBDNAPeaks_Direct = K562NonBDNAPeaks[["Direct"]], 
    NonBDNAPeaks_Inverted = K562NonBDNAPeaks[["Inverted"]], 
    NonBDNAPeaks_Mirror = K562NonBDNAPeaks[["Mirror"]], 
    NonBDNAPeaks_ZDNAMotifs = K562NonBDNAPeaks[["ZDNAMotifs"]], 
    NonBDNAPeaks_ShortTandem = K562NonBDNAPeaks[["ShortTandem"]], 
    NonBDNAPeaks_APhased = K562NonBDNAPeaks[["APhased"]], 
    NonBDNAPeaks_GQuadruplexForming = K562NonBDNAPeaks[["GQuadruplexForming"]]
) 
K562Epigenomes <- c(K562Epigenomes, K562ChIPBHMs, K562ChIPNHMs)
K562Epigenomes <- c(K562Epigenomes, "names<-"(K562ChIPTFs, paste0("TF_", names(K562ChIPTFs))))

ncores <- 10
epigenomes <- names(K562Epigenomes)
K562AssayEpiAssoc <- sapply(c("CHEX", "ATAC", "DNase", "FAIRE"), function(a) {
    if (a == "CHEX") { 
        q <- reduce(flank(GRsFull[["K562PositiveAll"]], width = 1000, start = TRUE, both = TRUE), ignore.strand = TRUE) 
    } else { 
        q <- K562Epigenomes[[a]]
    }
    res <- mclapply(epigenomes, function(x) {
        message(a, " ", x)
        Genome$testGRsOverlap(query = q, subject = K562Epigenomes[[x]], universe = GRsNoBlacklisted_human, ignore.strand = TRUE)
    }, mc.cores = ncores)
    names(res) <- epigenomes
    res
}, simplify = FALSE)
K562AssayEpiOddsRatio <- sapply(K562AssayEpigenomeAssoc, function(X) sapply(X, function(x) unname(x[1])))
dirname <- sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)
write.csv(K562AssayEpiOddsRatio, file = sprintf("%s/K562AssayEpiLog2oddsRatio.csv", dirname))

K562AssayEpiOddsRatio_rank <- K562AssayEpiOddsRatio[epigenomes[-c(1:3)], ]
K562AssayEpiOddsRatio_rank <- apply(K562AssayEpiOddsRatio_rank, 2, rank)
K562AssayEpiOddsRatio_rank <- K562AssayEpiOddsRatio_rank / nrow(K562AssayEpiOddsRatio_rank)
rownames(K562AssayEpiOddsRatio_rank) <- sub("K562ChIP", "", rownames(K562AssayEpiOddsRatio_rank))
rownames(K562AssayEpiOddsRatio_rank) <- sub("K562", "", rownames(K562AssayEpiOddsRatio_rank))
pdf(sprintf("%s/FourAssayFoldRankQuantiles.pdf", dirname), height = 24, width = 4)
pheatmap(K562AssayEpiOddsRatio_rank, fontsize_row = 4)
dev.off()

###########################################################################
## Human brain
###########################################################################
bioGroups_human <- c("HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture")
sampleIDsVirtmax_human <- sampleIDsVirtmaxByBioGroup[bioGroups_human]

HumanBrainATACNarrowPeaks <- readRDS("Data/Epigenome/HumanBrainATACNarrowPeaks.RDS")
HumanBrainDNaseNarrowPeaks <- readRDS("Data/Epigenome/HumanBrainDNaseNarrowPeaks.RDS")
HumanBrainDNaseNarrowPeaks <- lapply(HumanBrainDNaseNarrowPeaks, function(X) Reduce(union, X))
HumanBrainFAIRENarrowPeaks <- readRDS("Data/Epigenome/HumanBrainFAIRENarrowPeaks.RDS")
HumanNonBDNAPeaks <- readRDS(file = "Data/Epigenome/HumanNonBDNAPeaks.RDS")

names(HumanBrainATACNarrowPeaks) <- paste0("HumanBrainATACNarrowPeaks_", names(HumanBrainATACNarrowPeaks))
names(HumanBrainDNaseNarrowPeaks) <- paste0("HumanBrainDNaseNarrowPeaks_", names(HumanBrainDNaseNarrowPeaks))
names(HumanBrainFAIRENarrowPeaks) <- paste0("HumanBrainFAIRENarrowPeaks_", names(HumanBrainFAIRENarrowPeaks))
names(HumanNonBDNAPeaks) <- paste0("HumanNonBDNAPeaks_", names(HumanNonBDNAPeaks))

HumanBrainEpigenomes <- c(HumanBrainATACNarrowPeaks, HumanBrainDNaseNarrowPeaks, HumanBrainFAIRENarrowPeaks, HumanNonBDNAPeaks)

epigenomes_human <- names(HumanBrainEpigenomes)
HumanBrainAssayEpiAssoc <- sapply(sampleIDsVirtmax_human, function(a) {
    q <- reduce(flank(GRsFull[[a]], width = 1000, start = TRUE, both = TRUE), ignore.strand = TRUE) 
    res <- mclapply(epigenomes_human, function(x) {
        message(a, " ", x)
        Genome$testGRsOverlap(query = q, subject = HumanBrainEpigenomes[[x]], universe = GRsNoBlacklisted_human, ignore.strand = TRUE)
    }, mc.cores = 2)
    names(res) <- epigenomes_human
    res
}, simplify = FALSE)

## get log2OR and pval
HumanBrainAssayEpiOddsRatio <- sapply(HumanBrainAssayEpiAssoc, function(X) sapply(X, function(x) unname(x[1])))
HumanBrainAssayEpiPval <- sapply(HumanBrainAssayEpiAssoc, function(X) sapply(X, function(x) unname(x[2])))
write.csv(HumanBrainAssayEpiOddsRatio, file = sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/HumanBrainAssayEpiOddsRatio.csv", qualOutInPair, mapqTh))
write.csv(HumanBrainAssayEpiPval, file = sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/HumanBrainAssayEpiPval.csv", qualOutInPair, mapqTh))

rownames(HumanBrainAssayEpiOddsRatio) <- gsub("HumanBrain|NarrowPeaks|ChIPNHMs_|ChIPBHMs_", "", rownames(HumanBrainAssayEpiOddsRatio))
rownames(HumanBrainAssayEpiOddsRatio) <- gsub("HumanNonBDNAPeaks", "NonBDNA", rownames(HumanBrainAssayEpiOddsRatio))
pdf(sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/HumanBrainAssayEpiOddsRatio_brain-nonbdna_heatmap.pdf", qualOutInPair, mapqTh), width = 5, height = 6)
pheatmap(HumanBrainAssayEpiOddsRatio, col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), breaks = seq(-max(abs(HumanBrainAssayEpiOddsRatio)), max(abs(HumanBrainAssayEpiOddsRatio)), length.out = 255), angle = 45, cellheight = 12, cellwidth = 12, border = NA)
dev.off()

###########################################################################
## Mouse brain
###########################################################################
bioGroups_mouse <- c("MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice")
sampleIDsVirtmax_mouse <- sampleIDsVirtmaxByBioGroup[bioGroups_mouse]

MouseBrainATACNarrowPeaks <- readRDS(file = "Data/Epigenome/MouseBrainATACNarrowPeaks.RDS")
MouseBrainDNaseNarrowPeaks <- readRDS("Data/Epigenome/MouseBrainDNaseNarrowPeaks.RDS")
MouseBrainChIPTFs <- readRDS(file = "Data/Epigenome/MouseBrainChIPTFs.RDS")
MouseBrainChIPNHMs <- readRDS(file = "Data/Epigenome/MouseBrainChIPNHMs.RDS")
MouseBrainChIPBHMs <- readRDS(file = "Data/Epigenome/MouseBrainChIPBHMs.RDS")
MouseBrainSEPeaks <- readRDS(file = "Data/Epigenome/MouseBrainSEPeaks.RDS")
## We finally excluded VISTA as they are too sparse.
MouseNonBDNAPeaks <- readRDS(file = "Data/Epigenome/MouseNonBDNAPeaks.RDS")

names(MouseBrainATACNarrowPeaks) <- paste0("MouseBrainATACNarrowPeaks_", names(MouseBrainATACNarrowPeaks))
names(MouseBrainDNaseNarrowPeaks) <- paste0("MouseBrainDNaseNarrowPeaks_", names(MouseBrainDNaseNarrowPeaks))
MouseBrainChIPTFs_flat <- list()
for (x in names(MouseBrainChIPTFs)) {
    for (y in names(MouseBrainChIPTFs[[x]])) {
        MouseBrainChIPTFs_flat[[paste0("MouseBrainChIPTFs_", x, "_", y)]] <- MouseBrainChIPTFs[[x]][[y]]
    }
}
MouseBrainChIPNHMs_flat <- list()
for (x in names(MouseBrainChIPNHMs)) {
    for (y in names(MouseBrainChIPNHMs[[x]])) {
        MouseBrainChIPNHMs_flat[[paste0("MouseBrainChIPNHMs_", x, "_", y)]] <- MouseBrainChIPNHMs[[x]][[y]]
    }
}
MouseBrainChIPBHMs_flat <- list()
for (x in names(MouseBrainChIPBHMs)) {
    for (y in names(MouseBrainChIPBHMs[[x]])) {
        MouseBrainChIPBHMs_flat[[paste0("MouseBrainChIPBHMs_", x, "_", y)]] <- MouseBrainChIPBHMs[[x]][[y]]
    }
}
names(MouseBrainSEPeaks) <- paste0("MouseBrainSEPeaks_", names(MouseBrainSEPeaks))
names(MouseNonBDNAPeaks) <- paste0("MouseNonBDNAPeaks_", names(MouseNonBDNAPeaks))

MouseBrainEpigenomes <- c(MouseBrainATACNarrowPeaks, MouseBrainDNaseNarrowPeaks, MouseBrainChIPTFs_flat, MouseBrainChIPNHMs_flat, MouseBrainChIPBHMs_flat, MouseBrainSEPeaks, MouseNonBDNAPeaks)

epigenomes_mouse <- names(MouseBrainEpigenomes)
MouseBrainAssayEpiAssoc <- sapply(sampleIDsVirtmax_mouse, function(a) {
    q <- reduce(flank(GRsFull[[a]], width = 1000, start = TRUE, both = TRUE), ignore.strand = TRUE) 
    res <- mclapply(epigenomes_mouse, function(x) {
        message(a, " ", x)
        Genome$testGRsOverlap(query = q, subject = MouseBrainEpigenomes[[x]], universe = GRsNoBlacklisted_mouse, ignore.strand = TRUE)
    }, mc.cores = 2) # more than 64GB is needed if mc.cores > 3
    names(res) <- epigenomes_mouse
    res
}, simplify = FALSE)

## get log2OR and pval
MouseBrainAssayEpiOddsRatio <- sapply(MouseBrainAssayEpiAssoc, function(X) sapply(X, function(x) unname(x[1])))
MouseBrainAssayEpiPval <- sapply(MouseBrainAssayEpiAssoc, function(X) sapply(X, function(x) unname(x[2])))
write.csv(MouseBrainAssayEpiOddsRatio, file = sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/MouseBrainAssayEpiOddsRatio.csv", qualOutInPair, mapqTh))
write.csv(MouseBrainAssayEpiPval, file = sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/MouseBrainAssayEpiPval.csv", qualOutInPair, mapqTh))

MouseBrainAssayEpiOddsRatio <- MouseBrainAssayEpiOddsRatio[grepl("NonBDNA|Adult", rownames(MouseBrainAssayEpiOddsRatio)), ]
rownames(MouseBrainAssayEpiOddsRatio) <- gsub("MouseBrain|NarrowPeaks|ChIPNHMs_|ChIPBHMs_", "", rownames(MouseBrainAssayEpiOddsRatio))
rownames(MouseBrainAssayEpiOddsRatio) <- gsub("MouseNonBDNAPeaks", "NonBDNA", rownames(MouseBrainAssayEpiOddsRatio))
pdf(sprintf("Report/ChexEpigenome/GenomicAssocTest/%s_%s/MouseBrainAssayEpiOddsRatio_adultchip-nonbdna_heatmap.pdf", qualOutInPair, mapqTh), width = 7, height = 6)
pheatmap(MouseBrainAssayEpiOddsRatio, col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), breaks = seq(-max(abs(MouseBrainAssayEpiOddsRatio)), max(abs(MouseBrainAssayEpiOddsRatio)), length.out = 255), angle = 45, cellheight = 12, cellwidth = 12, border = NA)
dev.off()
