## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/release/functions.R")

EnsFeatures <- readRDS("Data/release/GenomicFeatures/EnsFeatures.RDS")
EnsFlanksByGene <- readRDS("Data/release/GenomicFeatures/EnsFlanksByGene.RDS")
EnsFeatureIDsProt <- readRDS("Data/release/GenomicFeatures/EnsFeatureIDsProt.RDS")
Species <- c("human", "mouse")
featureTypesByGene <- c("FlankUp5kByGene", "FlankDn5kByGene", "ExonByGene", "IntronByGene")
EnsFeaturesByGene <- sapply(Species, function(species) sapply(featureTypesByGene, function(featureTypeByGene) { message(species, " ", featureTypeByGene); EnsFeatures[[featureTypeByGene]][[species]] }, simplify = FALSE), simplify = FALSE)
EnsGeneIDsProt <- lapply(EnsFeatureIDsProt, function(X) X[["Gene"]])
EnsIDToSymbolMap <- sapply(Species, function(species) structure(mcols(EnsFeatures[["Gene"]][[species]])[["symbol"]], names = names(EnsFeatures[["Gene"]][[species]])), simplify = FALSE)

SampleInfoFull <- read.csv("Data/release/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[["SampleID"]]

## limit to the samples we finally decided to include
Species <- c("human", "mouse")
bioGroups <- c(
   "K562", 
   "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture", 
   "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice"
)

trtGroups <- c("PositiveSingle", "PositivePooled")

SampleInfoPositive <- subset(SampleInfoFull, CompType == "Biol" & IsOut == "N" & IsNegCtrl == "N" & (Species %in% .GlobalEnv$Species) & BioGroup %in% bioGroups)
sampleIDsPositive <- SampleInfoPositive[["SampleID"]]

bioGroupCols <- structure(if (length(bioGroups) <= 9) { brewer.pal(n = length(bioGroups), "Set1") } else { colorRampPalette(brewer.pal(n = 9, "Set1"))(length(bioGroups)) }, names = bioGroups)
bioGroupsBySpecies <- sapply(Species, function(species) bioGroups[bioGroups %in% unique(subset(SampleInfoFull[sampleIDsPositive, ], Species == species)[["BioGroup"]])], simplify = FALSE)
bioGroupColsBySpecies <- sapply(Species, function(species) bioGroupCols[bioGroupsBySpecies[[species]]], simplify = FALSE)

trtGroupCols <- structure(brewer.pal(n = length(trtGroups), "Greys")[seq(trtGroups)], names = trtGroups)
trtGroupsBySpecies <- sapply(Species, function(species) trtGroups[trtGroups %in% unique(subset(SampleInfoFull[sampleIDsPositive, ], Species == species)[["TrtGroup"]])], simplify = FALSE)
trtGroupColsBySpecies <- sapply(Species, function(species) trtGroupCols[trtGroupsBySpecies[[species]]], simplify = FALSE)

Strandedness <- c("sense", "antisense")

qualOutInPairs <- c( "ABreadCmate5End", "ABread5End", "Aread5End", "Bread5End", "Cmate5End", "Deither5End" )[1]
## Remember that Deither has no barcode strand information.
mapqThs <- c("ge30_le0.1", "ge20_le0.1_strict", "ge20_le0.1", "ge10_le0.1")[2]

GRsFiltered <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, th)
        message(filename)
        readRDS(filename)
    }, simplify = FALSE)
}, simplify = FALSE)

###########################################################################
## Get 5End read counts per featureTypeByGene
## We cannot test intergenic, because it is not stranded.
###########################################################################
dirname <- "Data/release/Strandedness"
dir.create(dirname, FALSE, TRUE)
for (mapqTh in mapqThs) {
    for (qualOutInPair in qualOutInPairs) {
        filename <- sprintf("Data/release/Strandedness/PFs%sFilteredEnsIDStranded_%s.RDS", qualOutInPair, mapqTh)
        X <- sapply(Species, function(species) {
            SIDs <- subset(SampleInfoFull[sampleIDsPositive, ], Species == species)[["SampleID"]]
            sapply(featureTypesByGene, function(featureTypeByGene) {
                features <- EnsFeaturesByGene[[species]][[featureTypeByGene]]
                sapply(Strandedness, function(strandedness) {
                    sapply(SIDs, function(SID) {
                    message(mapqTh, " ", qualOutInPair, " ", species, " ", featureTypeByGene, " ", strandedness, " ", SID)
                    GRs <- GRsFiltered[[mapqTh]][[qualOutInPair]][[SID]]
                    Genome$getCountsByStrandedness(GRs, features = features, strandedness = strandedness)
                    })
                }, simplify = FALSE)
            }, simplify = FALSE)
		}, simplify = FALSE)
        saveRDS(X, file = filename)
	}
}

PFsFilteredEnsIDStranded <- sapply(mapqThs, function(mapqTh) 
    sapply(qualOutInPairs, function(qualOutInPair) {
        filename <- sprintf("Data/release/Strandedness/PFs%sFilteredEnsIDStranded_%s.RDS", qualOutInPair, mapqTh)
        message(filename)
        readRDS(filename)
    }, simplify = FALSE)
, simplify = FALSE)

###########################################################################
## We first do it to the overall expressions (sum all features up).
## This can avoid the issue of different feature numbers between human and mouse.
###########################################################################
cutoff <- 50
PFsFilteredEnsIDStranded_TotCnts <- sapply(mapqThs, function(mapqTh) 
    sapply(qualOutInPairs, function(qualOutInPair) {
        sapply(featureTypesByGene, function(featureTypeByGene) {
            sapply(sampleIDsPositive, function(SID) { 
                message(featureTypeByGene, " ", SID)
                species <- SampleInfoFull[SID, "Species"]
                PFs <- PFsFilteredEnsIDStranded[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]]
                sense <- PFs[["sense"]][, SID, drop = FALSE]
                antisense <- PFs[["antisense"]][, SID, drop = FALSE]
                if ((sum(sense) >= cutoff) && (sum(antisense) >= cutoff)) {
                    return(c(sense = sum(sense), antisense = sum(antisense)))
                } else {
                    return(c(sense = 0, antisense = 0))
                }
            })
        }, simplify = FALSE)
    }, simplify = FALSE)
, simplify = FALSE)

pesudocnt <- 0
PFsFilteredEnsIDStranded_TotCntRatios <- sapply(mapqThs, function(mapqTh) 
    sapply(qualOutInPairs, function(qualOutInPair) {
        sapply(featureTypesByGene, function(featureTypeByGene) {
            message(featureTypeByGene)
                totCnts <- PFsFilteredEnsIDStranded_TotCnts[[mapqTh]][[qualOutInPair]][[featureTypeByGene]]
                A2S = (pesudocnt + totCnts["antisense", ]) / (pesudocnt + totCnts["sense", ])
        })
    }, simplify = FALSE)
, simplify = FALSE)

dirname <- "Report/release/Strandedness"
dir.create(dirname, FALSE, TRUE)
for (qualOutInPair in qualOutInPairs)
    for (mapqTh in mapqThs) {
        filename <- sprintf("%s/PFs%sFilteredEnsIDStranded_TotCntRatios_%s.csv", dirname, qualOutInPair, mapqTh)
        write.csv(PFsFilteredEnsIDStranded_TotCntRatios[[mapqTh]][[qualOutInPair]], file = filename)
}

for (qualOutInPair in qualOutInPairs) {
    for (mapqTh in mapqThs) {
        mat <- PFsFilteredEnsIDStranded_TotCntRatios[[mapqTh]][[qualOutInPair]]
        df <- melt(mat)
        colnames(df) <- c("SampleID", "FeatureTypeByGene", "AntisenseToSense")
        df$SampleID <- as.character(df$SampleID)
        df$BioGroup <- factor(SampleInfoFull[df[["SampleID"]], "BioGroup"], levels = bioGroups)
        df$TrtGroup <- factor(c("Single" = "Single", "Pooled" = "Population")[sub("^Positive", "", SampleInfoFull[df[["SampleID"]], "TrtGroup"])], levels = c("Single", "Population"))
        filename <- sprintf("Report/release/Strandedness/PFs%sFilteredEnsIDStranded_TotCntRatios_%s_violin.pdf",qualOutInPair, mapqTh)
        fig <- ggplot(df, aes(x = FeatureTypeByGene, y = log2(AntisenseToSense))) + geom_violin(lwd = 0.2) + geom_signif(comparisons = list(c("FlankUp5kByGene", "FlankDn5kByGene")), test = "wilcox.test", test.args = list(alternative = "two.sided"), textsize = 2, margin_top = +0.05, map_signif_level = FALSE)  + geom_signif(comparisons = list(c("ExonByGene", "IntronByGene")), test = "wilcox.test", test.args = list(alternative = "two.sided"), textsize = 2, margin_top = +0.05, map_signif_level = FALSE) + theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.key.size = unit(0.1, "in"), legend.text= element_text(size = 6)) + geom_jitter(data = df, aes(color = BioGroup, shape = TrtGroup), size = 0.2) + xlab("") + ylab("Log2 antisense-to-sense ratio (total counts)") + scale_x_discrete(limits = featureTypesByGene, labels = sub("ByGene$", "", featureTypesByGene))  + geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) + guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) 
        ggsave(fig, filename = filename, width = 3, height = 3, useDingbats = FALSE)
    }
}
