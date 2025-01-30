## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/functions.R")

EnsFeatures <- readRDS("Data/GenomicFeatures/EnsFeatures.RDS")
EnsFeatureIDsChrM <- readRDS("Data/GenomicFeatures/EnsFeatureIDsChrM.RDS")
EnsGeneIDsChrM <- list(human = EnsFeatureIDsChrM$human$Gene, mouse = EnsFeatureIDsChrM$mouse$Gene, rat = EnsFeatureIDsChrM$rat$Gene, none = EnsFeatureIDsChrM$mouse$Gene)
EnsSymbolsChrM <- list(human = mcols(EnsFeatures$Gene$human[EnsGeneIDsChrM$human])[["symbol"]], 
                       mouse = mcols(EnsFeatures$Gene$mouse[EnsGeneIDsChrM$mouse])[["symbol"]], 
                       rat = mcols(EnsFeatures$Gene$rat[EnsGeneIDsChrM$rat])[["symbol"]], 
                       none = mcols(EnsFeatures$Gene$mouse[EnsGeneIDsChrM$none])[["symbol"]]
)

EnsFeaturesGene_human_mito_df <- as.data.frame(EnsFeatures$Gene$human[EnsFeatureIDsChrM$human$Gene])
rownames(EnsFeaturesGene_human_mito_df) <- EnsFeaturesGene_human_mito_df$symbol
EnsFeaturesGene_human_mito_df$Color <- rep(c("grey", "green", "grey", "yellow", "grey", "orange", "grey"), c(1, 5, 3, 19, 3, 4, 2))
EnsFeaturesGene_human_mito_df$Block <- c("grey" = "Gap", "green" = "Region1", "yellow" = "Region2", "orange" = "Region3")[as.character(EnsFeaturesGene_human_mito_df$Color)]

EnsFeaturesGene_mouse_mito_df <- as.data.frame(EnsFeatures$Gene$mouse[EnsFeatureIDsChrM$mouse$Gene])
rownames(EnsFeaturesGene_mouse_mito_df) <- EnsFeaturesGene_mouse_mito_df$symbol
EnsFeaturesGene_mouse_mito_df$Color <- rep(c("green", "grey", "yellow", "grey", "orange", "grey"), c(16, 2, 4, 2, 11, 2))
EnsFeaturesGene_mouse_mito_df$Block <- c("grey" = "Gap", "green" = "Region1", "yellow" = "Region2", "orange" = "Region3")[as.character(EnsFeaturesGene_mouse_mito_df$Color)]

SampleInfoFull <- read.csv("Data/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
SampleInfoFull <- subset(SampleInfoFull, IsOut == "N")
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[["SampleID"]]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
sampleIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual")[, "SampleID"]
sampleIDs <- SampleInfo$SampleID

bioGroups_human <- c("K562", "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture")
bioGroups_mouse <- c("MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice")
bioGroups <- c(bioGroups_human, bioGroups_mouse)

sampleIDsVirtmaxByBioGroup <- sapply(bioGroups, function(bioGroup) tail(subset(SampleInfoFull, IsNegCtrl == "N" & ProbeType == "Positive" & BioGroup == bioGroup), 1)[["SampleID"]])
sampleIDsFullByBioGroup <- sapply(bioGroups, function(bioGroup) subset(SampleInfoFull, IsNegCtrl == "N" & ProbeType == "Positive" & BioGroup == bioGroup)[["SampleID"]])


qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
dirname <- sprintf("Report/Mito/%s_%s", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)

filename <- sprintf("Data/PrimingRateGene/PFs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
PFs <- readRDS(filename)

###########################################################################
## Number of priming sites per mt-gene
## Number of genes primed per sample
###########################################################################
## Mt-genes are small, so we should use the gene-body intervals without extension. 
PFs_mito_human_Gene <- PFs$human$Gene[EnsSymbolsChrM$human, subset(SampleInfoFull, BioGroup %in% bioGroups_human & Species == "human" & IsOut == "N" & IsNegCtrl == "N" & CompType == "Biol")[["SampleID"]]]
PFs_mito_human_Gene <- PFs_mito_human_Gene[, colSums(PFs_mito_human_Gene > 0) > 0]
sampleIDs_human <- colnames(PFs_mito_human_Gene)
SampleInfo_human <- SampleInfoFull[sampleIDs_human, ]
SampleInfo_human$BioGroup <- factor(SampleInfo_human$BioGroup, levels = bioGroups_human)
SampleInfo_human$Composition <- factor(SampleInfo_human$Composition, levels = c("Single", "Pooled"))
SampleInfo_human_K562 <- subset(SampleInfo_human, BioGroup %in% bioGroups_human_K562)
sampleIDs_human_K562 <- SampleInfo_human_K562[, "SampleID"]
sampleIDs_human_sorted <- sampleIDs_human[order(SampleInfo_human$Composition, SampleInfo_human$BioGroup)]
sampleIDs_human_K562_sorted <- sampleIDs_human_K562[order(SampleInfo_human_K562$Composition, SampleInfo_human_K562$BioGroup)]
PFs_mito_human_Gene <- PFs_mito_human_Gene[, sampleIDs_human_sorted]

colnames(PFs_mito_human_Gene) <- with(SampleInfo_human[colnames(PFs_mito_human_Gene), ], paste(sub("^Human", "", BioGroup), Composition, sub("scCLTdegenNuc", "#", SampleID), sep = "_"))
PFs_mito_human_Gene_df <- data.frame(EnsFeaturesGene_human_mito_df[EnsSymbolsChrM$human, c("start", "end", "width", "strand")], PFs_mito_human_Gene[EnsSymbolsChrM$human, ], check.names = FALSE, stringsAsFactors = FALSE)
dirname <- sprintf("Report/Mito/%s_%s", qualOutInPair, mapqTh)
filename <- sprintf("%s/PFs_mito_human_Gene.csv", dirname)
write.csv(PFs_mito_human_Gene_df, file = filename)

PFs_mito_mouse_Gene <- PFs$mouse$Gene[EnsSymbolsChrM$mouse, subset(SampleInfoFull, BioGroup %in% bioGroups_mouse & Species == "mouse" & IsOut == "N" & IsNegCtrl == "N" & CompType == "Biol")[["SampleID"]]]
PFs_mito_mouse_Gene <- PFs_mito_mouse_Gene[, colSums(PFs_mito_mouse_Gene > 0) > 0]
sampleIDs_mouse <- colnames(PFs_mito_mouse_Gene)
SampleInfo_mouse <- SampleInfoFull[sampleIDs_mouse, ]
SampleInfo_mouse$BioGroup <- factor(SampleInfo_mouse$BioGroup, levels = bioGroups_mouse)
SampleInfo_mouse$Composition <- factor(SampleInfo_mouse$Composition, levels = c("Single", "Pooled"))

colnames(PFs_mito_mouse_Gene) <- with(SampleInfo_mouse[colnames(PFs_mito_mouse_Gene), ], paste(sub("^Mouse", "", BioGroup), Composition, sub("scCLTdegenNuc", "#", SampleID), sep = "_"))
PFs_mito_mouse_Gene_df <- data.frame(EnsFeaturesGene_mouse_mito_df[EnsSymbolsChrM$mouse, c("start", "end", "width", "strand")], PFs_mito_mouse_Gene[EnsSymbolsChrM$mouse, ], check.names = FALSE, stringsAsFactors = FALSE)
dirname <- sprintf("Report/Mito/%s_%s", qualOutInPair, mapqTh)
filename <- sprintf("%s/PFs_mito_mouse_Gene.csv", dirname)
write.csv(PFs_mito_mouse_Gene_df, file = filename)

###########################################################################
## Heavy- vs. Light-strand priming rate comparison
###########################################################################
qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
filename <- sprintf("Data/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)

GRs <- readRDS(filename)
GRsVirtual <- sapply(sampleIDsVirtual, function(x) {
    message(x)
    sourceIDs <- SampleInfoFull[x, "SourceIDsNoOut"]
    SIDs <- strsplit(sourceIDs, ",")[[1]]
    grs <- GRs[SIDs]
    Reduce(c, grs)
}, simplify = FALSE)
GRsFull <- c(GRs, GRsVirtual)[sampleIDsFull]

PFs_mt_perbase_human <- sapply(bioGroups_human, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        message(SID)
        grs <- GRsFull[[SID]]
        h <- grs[strand(grs) == "+"]
        l <- grs[strand(grs) == "-"]
        h <- coverage(h)[["chrM"]]
        l <- coverage(l)[["chrM"]]
        data.frame(Heavy = h, Light = l)
    }, simplify = FALSE)
}, simplify = FALSE)

totcnts_mt_perbase_human <- sapply(bioGroups_human, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        X <- PFs_mt_perbase_human[[bioGroup]][[SID]]
        colSums(X)
    })
}, simplify = FALSE)

## We need to exclude bases without too few counts.
## The cutoff set to 10
heavytolight_mt_perbase_human <- sapply(bioGroups_human, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        x <- totcnts_mt_perbase_human[[bioGroup]][, SID]
        unname(ifelse(x[1] + x[2] >= 10, x["Heavy"] / x["Light"], NA))
    })
}, simplify = FALSE)

heavytolight_mt_perbase_human_df <- do.call(rbind, lapply(bioGroups_human, function(x) data.frame(SampleID = names(heavytolight_mt_perbase_human[[x]]), BioGroup = x, Heavy_to_Light = heavytolight_mt_perbase_human[[x]], row.names = NULL)))

PFs_mt_perbase_mouse <- sapply(bioGroups_mouse, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        message(SID)
        grs <- GRsFull[[SID]]
        h <- grs[strand(grs) == "+"]
        l <- grs[strand(grs) == "-"]
        h <- coverage(h)[["chrM"]]
        l <- coverage(l)[["chrM"]]
        data.frame(Heavy = h, Light = l)
    }, simplify = FALSE)
}, simplify = FALSE)

totcnts_mt_perbase_mouse <- sapply(bioGroups_mouse, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        X <- PFs_mt_perbase_mouse[[bioGroup]][[SID]]
        colSums(X)
    })
}, simplify = FALSE)

heavytolight_mt_perbase_mouse <- sapply(bioGroups_mouse, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        x <- totcnts_mt_perbase_mouse[[bioGroup]][, SID]
        unname(ifelse(x[1] + x[2] >= 10, x["Heavy"] / x["Light"], NA))
    })
}, simplify = FALSE)

heavytolight_mt_perbase_mouse_df <- do.call(rbind, lapply(bioGroups_mouse, function(x) data.frame(SampleID = names(heavytolight_mt_perbase_mouse[[x]]), BioGroup = x, Heavy_to_Light = heavytolight_mt_perbase_mouse[[x]], row.names = NULL)))

heavytolight_mt_perbase_df <- rbind(heavytolight_mt_perbase_human_df, heavytolight_mt_perbase_mouse_df)
filename <- sprintf("%s/totcnts_heavytolight_mt_perbase_boxplot.pdf", dirname)
pdf(filename, height = 7, width = 3.5)
ggplot(heavytolight_mt_perbase_df, aes(x = BioGroup, y = log2(Heavy_to_Light))) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.4) + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("") + ylab("Total counts heavy-to-light ratio (log2)") + geom_hline(yintercept = 0, lwd = 0.5, color = "blue", linetype = 2)
dev.off()

###########################################################################
## Inside vs. outside D-loop priming rate comparison
###########################################################################
## human D-loop is 16024-576, according to Sharma et al. 2005
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1352382/
cnts_dloop_perbase_human <- sapply(bioGroups_human, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        X <- PFs_mt_perbase_human[[bioGroup]][[SID]]
        x <- rowSums(X)
        c(inside = mean(x[c(1:576, 16024:16569)]), outside = mean(x[-c(1:576, 16024:16569)]))
    })
}, simplify = FALSE)

cnts_dloop_perbase_human_df <- melt(cnts_dloop_perbase_human)
colnames(cnts_dloop_perbase_human_df) <- c("Dloop", "SampleID", "Perbase", "BioGroup")
cnts_dloop_perbase_human_df$BioGroup <- factor(cnts_dloop_perbase_human_df$BioGroup, levels = bioGroups_human)

## Mouse D-loop is 15417-16295, according to Bibb et al. 1981
## https://pubmed.ncbi.nlm.nih.gov/7332926/
cnts_dloop_perbase_mouse <- sapply(bioGroups_mouse, function(bioGroup) {
    SIDs <- sampleIDsFullByBioGroup[[bioGroup]]
    sapply(SIDs, function(SID) {
        X <- PFs_mt_perbase_mouse[[bioGroup]][[SID]]
        x <- rowSums(X)
        c(inside = mean(x[15417:16295]), outside = mean(x[-(15417:16295)]))
    })
}, simplify = FALSE)

cnts_dloop_perbase_mouse_df <- melt(cnts_dloop_perbase_mouse)
colnames(cnts_dloop_perbase_mouse_df) <- c("Dloop", "SampleID", "Perbase", "BioGroup")
cnts_dloop_perbase_mouse_df$BioGroup <- factor(cnts_dloop_perbase_mouse_df$BioGroup, levels = bioGroups_mouse)

cnts_dloop_perbase_df <- rbind(cnts_dloop_perbase_human_df, cnts_dloop_perbase_mouse_df)
cnts_dloop_perbase_df$BioGroup <- factor(cnts_dloop_perbase_df$BioGroup, levels = bioGroups)
dirname <- sprintf("Report/Mito/%s_%s", qualOutInPair, mapqTh)
filename <- sprintf("%s/totcnts_dloop_perbase.pdf", dirname)
pdf(filename, height = 7, width = 5, useDingbats = FALSE)
ggplot(cnts_dloop_perbase_df, aes(x = BioGroup, fill = Dloop, y = Perbase)) + geom_boxplot(size = 0.4, outlier.color = NA) + geom_jitter(size = 0.1) + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylab("Per-base density") + xlab("") + scale_y_continuous(trans = "log10")
dev.off()

###########################################################################
## Per-strand per-base coverages
###########################################################################
dirname <- sprintf("Report/Mito/%s_%s", qualOutInPair, mapqTh)
filename <- sprintf("%s/cnts_mito_coverage_stranded_human.pdf", dirname)
pdf(filename, height = 3.5, width = 7)
for (bioGroup in bioGroups_human) {
    Graphics$plot_mito_coverage_stranded(GRsFull[[sampleIDsVirtmaxByBioGroup[[bioGroup]]]], as.data.frame(EnsFeatures$Gene$human), main = sprintf("%s", bioGroup), ylab = "", mar = c(2.1, 2.1, 6.0, 0.1), texthjust1 = 60, texthjust2 = -60, textcex = 0.8)
}
dev.off()
filename <- sprintf("%s/cnts_mito_coverage_stranded_mouse.pdf", dirname)
pdf(filename, height = 3.5, width = 7)
for (bioGroup in bioGroups_mouse) {
    Graphics$plot_mito_coverage_stranded(GRsFull[[sampleIDsVirtmaxByBioGroup[[bioGroup]]]], as.data.frame(EnsFeatures$Gene$mouse), main = sprintf("%s", bioGroup), ylab = "", mar = c(2.1, 2.1, 6.0, 0.1), texthjust1 = 100, texthjust2 = -100, textcex = 0.8)
}
dev.off()

###########################################################################
## Strand merged per-bin coverages
###########################################################################
Species <- c("human", "mouse")
binsizes <- c(50, 100, 500, 1000)
binsizeLabs <- c("50", "100", "500", "1k")
names(binsizes) <- binsizeLabs
mitoBinsList <- lapply(binsizes, function(binsize) { 
    list(human = Genome$getTilesFromSeqLengths(Genome$MainSeqLengths$human["chrM"], width = binsize, seqInfo = Genome$MainSeqInfo$human),
         mouse = Genome$getTilesFromSeqLengths(Genome$MainSeqLengths$mouse["chrM"], width = binsize, seqInfo = Genome$MainSeqInfo$mouse)
    )
})

PFs_mito_perbin_human <- sapply(binsizeLabs, function(binsizeLab) {
    bins <- mitoBinsList[[binsizeLab]][["human"]]
    sapply(c(bioGroups_human), function(bioGroup) {
        sampleID <- sampleIDsVirtmaxByBioGroup[bioGroup]
        grs <- GRsFull[[sampleID]]
        grs <- grs[seqnames(grs) == "chrM"]
        cvgs <- coverage(grs)
        cvg <- cvgs[["chrM"]]
        PFs <- unlist(viewSums(Views(cvg, ranges(bins[["chrM"]]))))
        binnames <- unlist(lapply(bins, as.character))
        names(PFs) <- binnames
        PFs
    })
}, simplify = FALSE)

PFs_mito_perbin_mouse <- sapply(binsizeLabs, function(binsizeLab) {
    bins <- mitoBinsList[[binsizeLab]][["mouse"]]
    sapply(c(bioGroups_mouse), function(bioGroup) {
        sampleID <- sampleIDsVirtmaxByBioGroup[bioGroup]
        grs <- GRsFull[[sampleID]]
        grs <- grs[seqnames(grs) == "chrM"]
        cvgs <- coverage(grs)
        cvg <- cvgs[["chrM"]]
        PFs <- unlist(viewSums(Views(cvg, ranges(bins[["chrM"]]))))
        binnames <- unlist(lapply(bins, as.character))
        names(PFs) <- binnames
        PFs
    })
}, simplify = FALSE)


for (binsizeLab in binsizeLabs) {
    filename <- sprintf("%s/PFs_mito_zscore-perbin_human_%s.pdf", dirname, binsizeLab)
    pdf(filename, height = 3.5, width = 7)
    for (bioGroup in bioGroups_human) {
        Graphics$plot_mito_cnts_bin(scale(PFs_mito_perbin_human[[binsizeLab]][, bioGroup], center = TRUE), ann = as.data.frame(EnsFeatures$Gene$human), binsize = binsizes[binsizeLab], main = bioGroup, xlab = "", ylab = "", mar = c(2.1, 2.1, 6.0, 0.1), hline = c(1, 2), hlinelty = c(2, 2), textcex = 0.8)
    }
    dev.off()
}

for (binsizeLab in binsizeLabs) {
    filename <- sprintf("%s/PFs_mito_zscore-perbin_mouse_%s.pdf", dirname, binsizeLab)
    pdf(filename, height = 3.5, width = 7)
    for (bioGroup in bioGroups_mouse) {
        Graphics$plot_mito_cnts_bin(scale(PFs_mito_perbin_mouse[[binsizeLab]][, bioGroup], center = TRUE), ann = as.data.frame(EnsFeatures$Gene$mouse), binsize = binsizes[binsizeLab], main = bioGroup, xlab = "", ylab = "", mar = c(2.1, 2.1, 6.0, 0.1), hline = c(1, 2), hlinelty = c(2, 2), textcex = 0.8)
    }
    dev.off()
}
