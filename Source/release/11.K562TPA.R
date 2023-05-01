## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/release/functions.R")

SampleInfoFull <- read.csv("Data/release/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[, "SampleID"]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
SampleInfoVirtual <- subset(SampleInfoFull, CompType == "Virtual")
sampleIDs <- SampleInfo[, "SampleID"]
sampleIDsVirtual <- SampleInfoVirtual[, "SampleID"]

bioGroups <- c(
    "K562", "K562TPAnone", "K562TPA15min", "K562TPA1hr", "K562TPA2hr", "K562TPA24hr", 
    "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture",
    "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice", 
    "MouseLungEpithelialCulture", "MouseLungEndothelialCulture",
    "MouseKidneyEpithelialCulture",
    "MouseCardiomyoCulture", "MouseCardiomyoSlice",
    "RatCardiomyoCulture", 
    "NoCell", 
    "K562MungBean", 
    "HBR"
)
bioGroupCols <- structure(if (length(bioGroups) <= 9) { brewer.pal(n = length(bioGroups), "Set1") } else { colorRampPalette(brewer.pal(n = 9, "Set1"))(length(bioGroups)) }, names = bioGroups)
sampleIDsByBioGroup <- sapply(bioGroups, function(bioGroup) subset(SampleInfoFull, CompType == "Biol" & BioGroup == bioGroup)[["SampleID"]], simplify = FALSE)
sampleIDsVirtmaxByBioGroup <- sapply(bioGroups, function(bioGroup) tail(subset(SampleInfoFull, BioGroup == bioGroup), 1)[["SampleID"]], simplify = FALSE)
bioGroup2sampleIDVirtmax <- unlist(sampleIDsVirtmaxByBioGroup)
bioGroup2Species <- sapply(bioGroups, function(bioGroup) subset(SampleInfoFull, BioGroup == bioGroup)[1, "Species"])

qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
GRs <- readRDS(filename)
GRsVirtual <- sapply(sampleIDsVirtual, function(sampleID) {
    message(sampleID)
    sourceIDs <- SampleInfoFull[sampleID, "SourceIDsNoOut"]
    SIDs <- strsplit(sourceIDs, ",")[[1]]
    grs <- GRs[SIDs]
    Reduce(c, grs)
}, simplify = FALSE)
GRsFull <- c(GRs, GRsVirtual)[sampleIDsFull]

EnsFeatures <- readRDS("Data/release/GenomicFeatures/EnsFeatures.RDS")
FeatureIDsMainNoMY <- readRDS("Data/release/GenomicFeatures/FeatureIDsMainNoMY.RDS")
EnsGenesMainNoMY_human <- EnsFeatures[["Gene"]][["human"]][FeatureIDsMainNoMY[["human"]][["Gene"]]]

EnsTSSdn5kMainNoMY_human <- flank(EnsGenesMainNoMY_human, width = -5000, start = TRUE)
EnsTSSup5kMainNoMY_human <- flank(EnsGenesMainNoMY_human, width = 5000, start = TRUE)
EnsTSSFlank5kMainNoMY_human <- flank(EnsGenesMainNoMY_human, width = 5000, both = TRUE)

bioGroupsK562TPA <- c("K562TPAnone", "K562TPA15min", "K562TPA1hr", "K562TPA2hr", "K562TPA24hr")
sampleIDsVirtmaxK562TPA <- bioGroup2sampleIDVirtmax[bioGroupsK562TPA]

GRsVirtmatK562TPA <- GRsFull[sampleIDsVirtmaxK562TPA]

dnOverlapsK562TPA <- lapply(sampleIDsVirtmaxK562TPA, function(SID) {
    GRs <- GRsFull[[SID]]
    findOverlaps(subject = EnsTSSdn5kMainNoMY_human, query = GRs, ignore.strand = TRUE)
})
upOverlapsK562TPA <- lapply(sampleIDsVirtmaxK562TPA, function(SID) {
    GRs <- GRsFull[[SID]]
    findOverlaps(subject = EnsTSSup5kMainNoMY_human, query = GRs, ignore.strand = TRUE)
})

mappedTSSdn5kMainNoMY_human <- lapply(dnOverlapsK562TPA, function(Ov) {
    idx <- unique(subjectHits(Ov))
    EnsTSSdn5kMainNoMY_human[idx]     
})
mappedTSSup5kMainNoMY_human <- lapply(upOverlapsK562TPA, function(Ov) {
    idx <- unique(subjectHits(Ov))
    EnsTSSup5kMainNoMY_human[idx]     
})

mappedStrandsTSSdn5kVirtmaxK562TPA <- lapply(mappedTSSdn5kMainNoMY_human, strand)
mappedStrandsTSSup5kVirtmaxK562TPA <- lapply(mappedTSSup5kMainNoMY_human, strand)

mappedGRsTSSdn5kVirtmaxK562TPA <- sapply(bioGroupsK562TPA, function(bioGroup) {
    SID <- bioGroup2sampleIDVirtmax[bioGroup]
    GRs <- GRsFull[[SID]]
    Ov <- dnOverlapsK562TPA[[bioGroup]]
    idx <- unique(queryHits(Ov))
    GRs[idx]
}, simplify = FALSE)

mappedGRsTSSup5kVirtmaxK562TPA <- sapply(bioGroupsK562TPA, function(bioGroup) {
    SID <- bioGroup2sampleIDVirtmax[bioGroup]
    GRs <- GRsFull[[SID]]
    Ov <- upOverlapsK562TPA[[bioGroup]]
    idx <- unique(queryHits(Ov))
    GRs[idx]
}, simplify = FALSE)

mappedCvgsTSSdn5kVirtmaxK562TPA <- lapply(mappedGRsTSSdn5kVirtmaxK562TPA, coverage)
mappedCvgsTSSup5kVirtmaxK562TPA <- lapply(mappedGRsTSSup5kVirtmaxK562TPA, coverage)

## TSS +5kb
CntMatsTSSdn5kMainNoMYVirtmaxK562TPA <- sapply(bioGroupsK562TPA, function(bioGroup) {
    strands <- as.character(mappedStrandsTSSdn5kVirtmaxK562TPA[[bioGroup]])
    cvgs <- mappedCvgsTSSdn5kVirtmaxK562TPA[[bioGroup]]
    grs <- mappedTSSdn5kMainNoMY_human[[bioGroup]]
    symbols <- mcols(grs)[["symbol"]]
    coords <- as.character(grs)
    positions <- paste(coords, symbols, sep = " ")
    cnts <- cvgs[grs]
    cnts <- mapply(FUN = function(X, Y) { if (Y == "-") { rev(as.integer(X)) } else  { as.integer(X) } }, X = cnts, Y = strands, SIMPLIFY = FALSE)
    cnts <- do.call(rbind, cnts)
    rownames(cnts) <- positions
    cnts
}, simplify = FALSE)
CntSumsTSSdn5kMainNoMYVirtmaxK562TPA <- sapply(CntMatsTSSdn5kMainNoMYVirtmaxK562TPA, colSums)
CntSumCulsumsTSSdn5kMainNoMYVirtmaxK562TPA <- apply(CntSumsTSSdn5kMainNoMYVirtmaxK562TPA, 2, function(x) cumsum(x))
TotCntsTSSdn5kMainNoMYVirtmaxK562TPA <- colSums(CntSumsTSSdn5kMainNoMYVirtmaxK562TPA)
CntSumCulprobsTSSdn5kMainNoMYVirtmaxK562TPA <- t(t(CntSumCulsumsTSSdn5kMainNoMYVirtmaxK562TPA) / TotCntsTSSdn5kMainNoMYVirtmaxK562TPA)

## TSS -5kb
CntMatsTSSup5kMainNoMYVirtmaxK562TPA <- sapply(bioGroupsK562TPA, function(bioGroup) {
    strands <- as.character(mappedStrandsTSSup5kVirtmaxK562TPA[[bioGroup]])
    cvgs <- mappedCvgsTSSup5kVirtmaxK562TPA[[bioGroup]]
    grs <- mappedTSSup5kMainNoMY_human[[bioGroup]]
    symbols <- mcols(grs)[["symbol"]]
    coords <- as.character(grs)
    positions <- paste(coords, symbols, sep = " ")
    cnts <- cvgs[grs]
    cnts <- mapply(FUN = function(X, Y) { if (Y == "-") { rev(as.integer(X)) } else  { as.integer(X) } }, X = cnts, Y = strands, SIMPLIFY = FALSE)
    cnts <- do.call(rbind, cnts)
    rownames(cnts) <- positions
    cnts
}, simplify = FALSE)
CntSumsTSSup5kMainNoMYVirtmaxK562TPA <- sapply(CntMatsTSSup5kMainNoMYVirtmaxK562TPA, colSums)
CntSumCulsumsTSSup5kMainNoMYVirtmaxK562TPA <- apply(CntSumsTSSup5kMainNoMYVirtmaxK562TPA, 2, function(x) cumsum(x))
TotCntsTSSup5kMainNoMYVirtmaxK562TPA <- colSums(CntSumsTSSup5kMainNoMYVirtmaxK562TPA)
CntSumCulprobsTSSup5kMainNoMYVirtmaxK562TPA <- t(t(CntSumCulsumsTSSup5kMainNoMYVirtmaxK562TPA) / TotCntsTSSup5kMainNoMYVirtmaxK562TPA)

## TSS -5 to 5kb
CntSumsTSSflank5kMainNoMYVirtmaxK562TPA <- rbind(CntSumsTSSup5kMainNoMYVirtmaxK562TPA, CntSumsTSSdn5kMainNoMYVirtmaxK562TPA)
CntSumCulsumsTSSflank5kMainNoMYVirtmaxK562TPA <- apply(CntSumsTSSflank5kMainNoMYVirtmaxK562TPA, 2, function(x) cumsum(x))
TotCntsTSSflank5kMainNoMYVirtmaxK562TPA <- colSums(CntSumsTSSflank5kMainNoMYVirtmaxK562TPA)
CntSumCulprobsTSSflank5kMainNoMYVirtmaxK562TPA <- t(t(CntSumCulsumsTSSflank5kMainNoMYVirtmaxK562TPA) / TotCntsTSSflank5kMainNoMYVirtmaxK562TPA)

dirname <- "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict"
dir.create(dirname, FALSE, TRUE)
filename <- file.path(dirname, "CntSumCulprobsTSSflank5kMainNoMYVirtmaxK562TPA.pdf")
pdf(filename, width = 4.5, height = 4.5)
par(ps = 11, lend = 2, ljoin = 1, bty = "L", mar = c(2.5, 2.5, 0, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
matplot(x = seq(-4999, 5000, by = 1), y = CntSumCulprobsTSSflank5kMainNoMYVirtmaxK562TPA, type = "l", col = bioGroupCols[bioGroupsK562TPA], xlab = "Distance to TSS (bp)", ylab = "Cumulative distribution", main = "", lwd = 2)
legend("topleft", legend = sub("^K562TPA", "", bioGroupsK562TPA), lty = seq(bioGroupsK562TPA), lwd = 2, col = bioGroupCols[bioGroupsK562TPA], box.col = NA)
abline(b = 1/10000, a = 0.5, lty = 2, col = "red")
dev.off()

## AUC-type of stats
FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA <- 1 / ((colSums(CntSumCulprobsTSSflank5kMainNoMYVirtmaxK562TPA[1:5000, ]) + colSums(1 - CntSumCulprobsTSSflank5kMainNoMYVirtmaxK562TPA[5001:10000, ])) / 2500)
filename <- file.path(dirname, "FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA.pdf")
pdf(filename, width = 5, height = 3)
par(ps = 11, lend = 2, ljoin = 1, bty = "L", mar = c(0, 5, 3, 1), oma = c(0, 1, 1.5, 0), mgp = c(1.5, 0.5, 0))
barplot(rev(FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA), horiz = TRUE, las = 1, xaxt = "n", xlim = c(0, 1.5), border = FALSE)
axis(side = 3, at = seq(0, 1.5, by = 0.5))
title(main = "Fold of TSS enrichment", outer = TRUE)
dev.off()
dirname <- "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict"
dir.create(dirname, FALSE, TRUE)
write.csv(FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA, file = "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA.csv")

###########################################################################
## Diffrential analysis
###########################################################################
Species <- c("human", "mouse")
featureTypesByGene <- c("GeneExt")
qualOutInPairs <- c("ABreadCmate5End", "ABread5End", "Aread5End", "Bread5End", "Cmate5End", "Deither5End")[1]
mapqThs <- c("ge20_le0.1_strict", "ge30_le0.1", "ge20_le0.1", "ge10_le0.1")[1]

PFsFiltered <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        filename <- sprintf("Data/release/PrimingRateGene/PFs%sFiltered_%s.RDS", qualOutInPair, th)
        message(filename)
        readRDS(filename)
    }, simplify = FALSE)
}, simplify = FALSE)

PFsFilteredVirtual <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        sapply(Species, function(species) { 
            SIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual" & Species == species)[, "SampleID"]
            sapply(featureTypesByGene, function(featureTypeByGene) {
                sapply(SIDsVirtual, function(SIDVirtual) {
                    sourceIDs <- SampleInfoFull[SIDVirtual, "SourceIDsNoOut"]
                    SIDs <- strsplit(sourceIDs, ",")[[1]]
                    message(th, " ", qualOutInPair, " ", featureTypeByGene, " ", species, " ", SIDVirtual)
                    PF <- PFsFiltered[[th]][[qualOutInPair]][[species]][[featureTypeByGene]]
                    rowSums(PF[, SIDs])
                })
            }, simplify = FALSE)
        }, simplify = FALSE)
    }, simplify = FALSE)
}, simplify = FALSE)

PFsFilteredFull <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        sapply(Species, function(species) { 
            sapply(featureTypesByGene, function(featureTypeByGene) {
                message(th, " ", qualOutInPair, " ", species, " ", featureTypeByGene)
                    SIDs <- subset(SampleInfoFull, CompType == "Biol" & Species == species)[, "SampleID"]
                    SIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual" & Species == species)[, "SampleID"]
                    PF <- PFsFiltered[[th]][[qualOutInPair]][[species]][[featureTypeByGene]]
                    PFVirtual <- PFsFilteredVirtual[[th]][[qualOutInPair]][[species]][[featureTypeByGene]]
                    cbind(PF, PFVirtual)
            }, simplify = FALSE)
        }, simplify = FALSE)
    }, simplify = FALSE)
}, simplify = FALSE)

###########################################################################
## Acute-response genes
###########################################################################
sampleIDs_K562TPAnone <- subset(SampleInfo, BioGroup == "K562TPAnone")[, "SampleID"]
sampleIDs_K562TPA15min <- subset(SampleInfo, BioGroup == "K562TPA15min")[, "SampleID"]
PFsFiltered_K562TPA_GeneExt <- PFsFiltered$ge20_le0.1_strict$ABreadCmate5End$human$GeneExt[, c(sampleIDs_K562TPAnone, sampleIDs_K562TPA15min)]

th <- 0.6
sum(rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPAnone] > 0) >= th & rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPA15min] > 0) == 0)
## [1] 61
sum(rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPAnone] > 0) == 0 & rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPA15min] > 0) >= th)
## [1] 6
## in total 67 genes
genes <- rownames(PFsFiltered_K562TPA_GeneExt)
diffgenes <- genes[
    (rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPAnone] > 0) >= th & rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPA15min] > 0) == 0) | 
    (rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPAnone] > 0) == 0 & rowMeans(PFsFiltered_K562TPA_GeneExt[, sampleIDs_K562TPA15min] > 0) >= th)
]
length(diffgenes)
## [1] 67
cat(diffgenes, file = "Report/release/K562TPA/DiffPriming/ABreadCmate5EndFiltered_ge20_le0.1_strict/DiffPrimStats_Binary_GeneExt_K562TPA15min_vs_K562TPAnone_symbol.txt", sep = "\n")

diffgenes <- readLines("Report/release/K562TPA/DiffPriming/ABreadCmate5EndFiltered_ge20_le0.1_strict/DiffPrimStats_Binary_GeneExt_K562TPA15min_vs_K562TPAnone_symbol.txt")
pdf("Report/release/K562TPA/DiffPriming/ABreadCmate5EndFiltered_ge20_le0.1_strict/DiffPrimStats_Binary_GeneExt_K562TPA15min_vs_K562TPAnone.pdf", width = 10, height = 3)
pheatmap(t(log2(1 + PFsFiltered_K562TPA_GeneExt[diffgenes, ])), border_color = NA, annotation_row = SampleInfo[c(sampleIDs_K562TPAnone, sampleIDs_K562TPA15min), "BioGroup", drop = FALSE], show_rownames = FALSE, angle = 90, cellheight = 8, cellwidth = 8, treeheight_row = 20, treeheight_col = 20)
dev.off()

###########################################################################
## Time-dependent gene-body changes corresponding to the TSS changes
###########################################################################
sampleIDs_K562TPAnone <- subset(SampleInfo, BioGroup == "K562TPAnone")[, "SampleID"]
sampleIDs_K562TPA15min <- subset(SampleInfo, BioGroup == "K562TPA15min")[, "SampleID"]
sampleIDs_K562TPA1hr <- subset(SampleInfo, BioGroup == "K562TPA1hr")[, "SampleID"]
sampleIDs_K562TPA2hr <- subset(SampleInfo, BioGroup == "K562TPA2hr")[, "SampleID"]
sampleIDs_K562TPA24hr <- subset(SampleInfo, BioGroup == "K562TPA24hr")[, "SampleID"]
PFsFiltered_K562TPA_GeneExt_series <- PFsFiltered$ge20_le0.1_strict$ABreadCmate5End$human$GeneExt[, c(sampleIDs_K562TPAnone, sampleIDs_K562TPA15min, sampleIDs_K562TPA1hr, sampleIDs_K562TPA2hr, sampleIDs_K562TPA24hr)]
PFsLibNorm_K562TPA_GeneExt_series <- cbind(
    rowSums(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPAnone])/sum(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPAnone]), 
    rowSums(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA15min])/sum(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA15min]), 
    rowSums(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA1hr])/sum(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA1hr]), 
    rowSums(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA2hr])/sum(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA2hr]), 
    rowSums(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA24hr])/sum(PFsFiltered_K562TPA_GeneExt_series[, sampleIDs_K562TPA24hr])
)
colnames(PFsLibNorm_K562TPA_GeneExt_series) <- c("K562TPAnone", "K562TPA15min", "K562TPA1hr", "K562TPA2hr", "K562TPA24hr")
FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman <- apply(PFsLibNorm_K562TPA_GeneExt_series, 1, function(y) cor(y, FoldEnrichTSSflank5kMainNoMYVirtmaxK562TPA, method = "spearman"))
write.csv(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman, file = "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman.csv")
pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(256)
breaks <- seq(0, max(log2(1 + 1e6 * PFsLibNorm_K562TPA_GeneExt_series[names(na.omit(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman[abs(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman) >= 0.9])), ])), length.out = 255)
pdf("Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_+0.9_pheatmap.pdf", width = 5, height = 5)
pheatmap(log2(1 + 1e6 * PFsLibNorm_K562TPA_GeneExt_series[names(na.omit(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman[FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman >= 0.9])), ]), cluster_row = TRUE, cluster_col = FALSE, border = FALSE, annotation_col = SampleInfo[colnames(PFsFiltered_K562TPA_GeneExt_series), "BioGroup", drop = FALSE], show_colnames = FALSE, annotation_colors = list(BioGroup = structure(bioGroupCols[bioGroupsK562TPA], names = colnames(PFsLibNorm_K562TPA_GeneExt_series))), fontsize_row = 6, cellheight = 6, color = pal, breaks = breaks)
dev.off()
pdf("Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_-0.9_pheatmap.pdf", width = 5, height = 5)
pheatmap(log2(1 + 1e6 * PFsLibNorm_K562TPA_GeneExt_series[names(na.omit(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman[FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman <= -0.9])), ]), cluster_row = TRUE, cluster_col = FALSE, border = FALSE, annotation_col = SampleInfo[colnames(PFsFiltered_K562TPA_GeneExt_series), "BioGroup", drop = FALSE], show_colnames = FALSE, annotation_colors = list(BioGroup = structure(bioGroupCols[bioGroupsK562TPA], names = colnames(PFsLibNorm_K562TPA_GeneExt_series))), fontsize_row = 6, cellheight = 6, color = pal, breaks = breaks)
dev.off()

###########################################################################
## Functional enrichment analysis 
###########################################################################
OrgDbGOSymbolList <- readRDS(file = "Data/release/GenomicFeatures/Annotations/OrgDbGOSymbolList.RDS")
OrgDbGONameList <- readRDS(file = "Data/release/GenomicFeatures/Annotations/OrgDbGONameList.RDS")
GOonts <- c("MF", "BP", "CC")

pos <- enricher(gene = names(na.omit(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman[FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman >= 0.9])), pvalueCutoff = Inf, qvalueCutoff = Inf, universe = names(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman), pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, TERM2GENE = OrgDbGOSymbolList$human$MF, TERM2NAME = OrgDbGONameList$human$MF)
neg <- enricher(gene = names(na.omit(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman[FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman <= -0.9])), pvalueCutoff = Inf, qvalueCutoff = Inf, universe = names(FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman), pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, TERM2GENE = OrgDbGOSymbolList$human$MF, TERM2NAME = OrgDbGONameList$human$MF)

go_mf_pos <- as.data.frame(pos)
go_mf_neg <- as.data.frame(neg)
write.csv(go_mf_pos, file = "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_+0.9_enrich_GO_MF.csv", row.names = FALSE)
write.csv(go_mf_neg, file = "Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_-0.9_enrich_GO_MF.csv", row.names = FALSE)

go_mf_pos$OddsRatio <- sapply(strsplit(go_mf_pos[, "GeneRatio"], "/"), function(x) as.integer(x[1])/as.integer(x[2])) / sapply(strsplit(go_mf_pos[, "BgRatio"], "/"), function(x) as.integer(x[1])/as.integer(x[2]))
go_mf_neg$OddsRatio <- sapply(strsplit(go_mf_neg[, "GeneRatio"], "/"), function(x) as.integer(x[1])/as.integer(x[2])) / sapply(strsplit(go_mf_neg[, "BgRatio"], "/"), function(x) as.integer(x[1])/as.integer(x[2]))

padj_th <- 0.1
go_mf_pos_sig <- go_mf_pos[go_mf_pos[, "p.adjust"] < padj_th, ]
go_mf_neg_sig <- go_mf_neg[go_mf_neg[, "p.adjust"] < padj_th, ]
pdf("Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_+0.9_enrich_GO_MF.pdf", height = 3.5, width = 16)
ggplot(go_mf_pos_sig, aes(y = Description, x = p.adjust, fill = log2(OddsRatio))) + geom_bar(stat = "identity") + scale_y_discrete(limits = go_mf_pos_sig[rev(order(go_mf_pos_sig[, "p.adjust"])), "Description"]) + xlab("BH p-value") + ylab("") + theme_classic(16) + scale_fill_gradientn(colours=brewer.pal(9,"Reds"))
dev.off()

pdf("Report/release/K562TPA/CumulCurves/ABreadCmate5EndFiltered_ge20_le0.1_strict/FoldEnrichTSSflank5kLibNormCnts_GeneExt_spearman_-0.9_enrich_GO_MF.pdf", height = 8, width = 16)
ggplot(go_mf_neg_sig, aes(y = Description, x = p.adjust, fill = log2(OddsRatio))) + geom_bar(stat = "identity") + scale_y_discrete(limits = go_mf_neg_sig[rev(order(go_mf_neg_sig[, "p.adjust"])), "Description"]) + xlab("BH p-value") + ylab("") + theme_classic(16) + scale_fill_gradientn(colours = brewer.pal(9,"Blues"))
dev.off()
