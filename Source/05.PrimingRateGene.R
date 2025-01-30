## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/functions.R")

.ISONLINE <- Tools$is_online()
.PWD <- ifelse(.ISONLINE, "~/Workspace.cache/Chex-seq", ".")
EnsFeatures <- readRDS("Data/GenomicFeatures/EnsFeatures.RDS")
EnsFlanksByGene <- readRDS("Data/GenomicFeatures/EnsFlanksByGene.RDS")

## Load grouping information
SampleInfo <- read.csv("Data/SampleInfo20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfo) <- sampleIDs <- SampleInfo[, "SampleID"]

Species <- c("human", "mouse", "rat")
EnsIDToSymbolMap <- sapply(Species, function(species) structure(mcols(EnsFeatures[["Gene"]][[species]])[["symbol"]], names = names(EnsFeatures[["Gene"]][[species]])), simplify = FALSE)
EnsIDToSymbolMap$none <- EnsIDToSymbolMap$mouse
## Calculate priming frequency in each gene
## 0. in GeneExt
## 1. in Gene
## 2. in TSS +/-5kb flanking region
## 3. in promoters
## 4. in 5'UTR
## 5. in exons
## 6. in CDS
## 7. in introns
## 8. in 3'UTR
## 9. in Downstream
featureTypesByGene <- c("GeneExt", "Gene", "Flank5kByGene", "PromoterByGene", "FiveUTRByGene", "ExonByGene", "IntronByGene", "ThreeUTRByGene", "DownstreamByGene", "CDSByGene")
EnsFeaturesByGene <- sapply(Species, function(species) sapply(featureTypesByGene, function(featureTypeByGene) { message(species, " ", featureTypeByGene); EnsFeatures[[featureTypeByGene]][[species]] }, simplify = FALSE), simplify = FALSE)
## We aligned the blank controls to mouse, so we should treat them as the mm10 genome.
EnsFeaturesByGene$none <- EnsFeaturesByGene$mouse
EnsFlanksByGene$none <- EnsFlanksByGene$mouse

## Add Species == "none"
Species <- c(Species, "none")

mapqThs <- c("ge20_le0.1_strict", "ge30_le0.1", "ge20_le0.1", "ge10_le0.1")
qualOutInPairs <- c( "ABreadCmate5End", "ABread5End", "Aread5End", "Bread5End", "Cmate5End", "Deither5End" )
###########################################################################
## Contaminants filtered priming frequency, before hotspot filtering
###########################################################################
## Genic regions
for (i in seq_along(mapqThs)) {
    mapqTh <- mapqThs[i]
    for (qualOutInPair in qualOutInPairs) {
        filename <- sprintf("Data/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
        message(filename)
        GRs <- readRDS(filename)
        PFs <- sapply(Species, function(species) { 
            message(mapqTh, " ", qualOutInPair, " ", species)
            featuresByGene <- EnsFeaturesByGene[[species]]            
            SIDs <- subset(SampleInfo, Species == species)[, "SampleID"]
            Genome$getPrimingRatesPerFeature(featureRangesList = featuresByGene, featureTypesByGene = featureTypesByGene, GRsList = GRs, sampleIDs = SIDs)
        }, simplify = FALSE) 
		dirname <- sprintf("%s/Data/PrimingRateGene", .PWD)
        dir.create(dirname, FALSE, TRUE)
        filename <- sprintf("%s/PFs%sFilteredEnsID_%s.RDS", dirname, qualOutInPair, mapqTh)
        message(filename)
        saveRDS(PFs, file = filename)
    }
}

## TSS flanking regions
flankDistLabs <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
for (i in seq_along(mapqThs)) {
    mapqTh <- mapqThs[i]
    for (qualOutInPair in qualOutInPairs) {
        infile <- sprintf("Data/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
        GRs <- readRDS(infile)
        outfile <- sprintf("Data/PrimingRateGene/FlankPFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
        FlankPFsEnsID <- sapply(Species, function(species) {
            message(infile, " ", species)
            FlanksByGene <- EnsFlanksByGene[[species]] 
            SIDs <- subset(SampleInfo, Species == species)[, "SampleID"]
            Genome$getPrimingRatesPerFeature(featureRangesList = FlanksByGene, featureTypesByGene = flankDistLabs, GRsList = GRs, sampleIDs = SIDs)
        }, simplify = FALSE)
        saveRDS(FlankPFsEnsID, file = outfile)
    }
}

ncores <- 10
for (i in seq_along(mapqThs)) {
    mapqTh <- mapqThs[i]
    for (qualOutInPair in qualOutInPairs) {
        filename <- sprintf("Data/PrimingRateGene/PFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
        message(filename)
        PFsEnsID <- readRDS(filename)
        PFs <- sapply(Species, function(species) {
            sapply(featureTypesByGene, function(featureTypeByGene) { 
                message(mapqTh, " ", qualOutInPair, " ", species, " ", featureTypeByGene)
                mat <- PFsEnsID[[species]][[featureTypeByGene]]
                mat <- Stats$summarizeMatrix(mat, by = EnsIDToSymbolMap[[species]][rownames(mat)], ncores = ncores)
                idx <- rownames(mat) == "" | is.na(rownames(mat))
                mat[!idx, ]
            }, simplify = FALSE)
        }, simplify = FALSE)
		dirname <- sprintf("%s/Data/PrimingRateGene", .PWD)
        dir.create(dirname, FALSE, TRUE)
        filename <- sprintf("%s/PFs%sFiltered_%s.RDS", dirname, qualOutInPair, mapqTh)
        message(filename)
        saveRDS(PFs, file = filename)
    }
}

flankDistLabs <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
for (i in seq_along(mapqThs)) {
    mapqTh <- mapqThs[i]
    for (qualOutInPair in qualOutInPairs) {
        infile <- sprintf("Data/PrimingRateGene/FlankPFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
        PFsEnsID <- readRDS(infile)
        PFs <- sapply(Species, function(species) { 
            IDMap <- EnsIDToSymbolMap[[species]]
            sapply(flankDistLabs, function(flankDistLab) {
                message(infile, " ", species, " ", flankDistLab)
                X <- PFsEnsID[[species]][[flankDistLab]]
                IDs <- rownames(X)
                mat <- Stats$summarizeMatrix(X, by = IDMap[IDs], ncores = ncores)
                idx <- rownames(mat) == "" | is.na(rownames(mat))
                mat[!idx, ]
            }, simplify = FALSE)
        }, simplify = FALSE)
		dirname <- sprintf("%s/Data/PrimingRateGene", .PWD)
        dir.create(dirname, FALSE, TRUE)
        outfile <- sprintf("%s/FlankPFs%sFiltered_%s.RDS", dirname, qualOutInPair, mapqTh)
        saveRDS(PFs, file = outfile)
    }
}
 
