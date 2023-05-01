## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/release/functions.R")

SampleInfoFull <- read.csv("Data/release/SampleInfoFull20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfoFull) <- sampleIDsFull <- SampleInfoFull[, "SampleID"]
sampleIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual")[, "SampleID"]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
sampleIDs <- SampleInfo[, "SampleID"]
###########################################################################
## Load biological priming loci from E.chex
## We finally don't generate BED for combined data (AB, etc.) any more.
## Instead, we combine them here after loading each component class.
###########################################################################
qualOutInPairs <- c( 
    "Aread5End", "AreadAndFrag", "AbothFrag", 
    "Bread5End", "BreadAndFrag", "BbothFrag", 
    "Cmate5End", "CmateReadAndFrag", "CbothFrag", 
    "Deither5End", "DeitherReadAndFrag", "DbothFrag"
)
qualOutInPairNoAdjBases <- c(
    "Aeither.NoAdj", "Aeither.NoAdj", "Aeither.NoAdj", 
    "Beither.NoAdj", "Beither.NoAdj", "Beither.NoAdj", 
    "Ceither.NoAdj", "Ceither.NoAdj", "Ceither.NoAdj", 
    "Deither", "Deither", "Deither"
)
outTypes <- c( 
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag"
)
QualOutInPairConfig <- data.frame(QualOutInPair = qualOutInPairs, QualOutInPairNoAdjBase = qualOutInPairNoAdjBases, OutType = outTypes, stringsAsFactors = FALSE)

GRsQualOutInPairList <- lapply(1:nrow(QualOutInPairConfig), function(i) {
    qualOutInPair <- QualOutInPairConfig[i, "QualOutInPair"]
    qualOutInPairNoAdjBase <- QualOutInPairConfig[i, "QualOutInPairNoAdjBase"]
    outType <- QualOutInPairConfig[i, "OutType"]
    GRs <- sapply(sampleIDs, function(sampleID) {
        species <- SampleInfo[sampleID, "Species"]
        ## We treat water samples as mouse
        if (species == "none") { species <- "mouse" }
        SeqInfo <- Genome$MainSeqInfo[[species]]
        seqLevels <- Genome$MainSeqLevels[[species]]
        seqLengths <- Genome$MainSeqLengths[[species]]
        sampleName <- paste0("Sample_", sampleID)
        filename <- sprintf("Data/E.chex/analyzed/%s/star/%s.star.primaryNoDup.%s.%s.NoBlacklisted.bed", sampleName, sampleName, qualOutInPairNoAdjBase, outType)
        message(filename)
        Genome$getGRsFromBed(filename, SeqInfo = SeqInfo, seqLevels = seqLevels, seqLengths = seqLengths)
    }, simplify = FALSE)
})
names(GRsQualOutInPairList) <- qualOutInPairs

dirname <- "Data/release/PrimingRate"
dir.create(dirname, FALSE, TRUE)
for (qualOutInPair in qualOutInPairs) {
    GRs <- GRsQualOutInPairList[[qualOutInPair]]
    filename <- sprintf("%s/GRs%s.RDS", dirname, qualOutInPair)
    message(filename)
    saveRDS(GRs, file = filename)
}

## Load the data
GRsQualOutInPairList <- sapply(qualOutInPairs, function(qualOutInPair) {
    message(qualOutInPair)
    filename <- sprintf("Data/release/PrimingRate/GRs%s.RDS", qualOutInPair)
    readRDS(filename)
}, simplify = FALSE)

## merge Aread5End Bread5End to form ABread5End
GRsABread5End <- sapply(sampleIDs, function(x) {
    message(x)
    c(GRsQualOutInPairList$Aread5End[[x]], 
      GRsQualOutInPairList$Bread5End[[x]])
}, simplify = FALSE)
filename <- sprintf("Data/release/PrimingRate/GRs%s.RDS", "ABread5End")
saveRDS(GRsABread5End, file = filename)

## merge Aread5End Bread5End and Cmate5End to form ABreadCmate5End
GRsABreadCmate5End <- sapply(sampleIDs, function(x) {
    message(x)
    c(GRsQualOutInPairList$Aread5End[[x]], 
      GRsQualOutInPairList$Bread5End[[x]],
      GRsQualOutInPairList$Cmate5End[[x]])
}, simplify = FALSE)
filename <- sprintf("Data/release/PrimingRate/GRs%s.RDS", "ABreadCmate5End")
saveRDS(GRsABreadCmate5End, file = filename)

qualOutInPairs_merge_map <- list(
    "ABread5End" = c("Aread5End", "Bread5End"), 
    "ABreadAndFrag" = c("AreadAndFrag", "BreadAndFrag"), 
    "ABbothFrag" = c("AbothFrag", "BbothFrag"), 
    "ABreadCmate5End" = c("Aread5End", "Bread5End", "Cmate5End"), 
    "ABCreadAndFrag" = c("AreadAndFrag", "BreadAndFrag", "CmateReadAndFrag"), 
    "ABCbothFrag" = c("AbothFrag", "BbothFrag", "CbothFrag")
)
for (x in names(qualOutInPairs_merge_map)) {
    GRs_each <- sapply(qualOutInPairs_merge_map[[x]], function(y) {
        filename <- sprintf("Data/release/PrimingRate/GRs%s.RDS", y)
            message(x, " < ", y, " ", filename)
            readRDS(filename)
    }, simplify = FALSE)
    GRs_merged <- sapply(sampleIDs, function(x) {
        unname(unlist(as(lapply(GRs_each, function(GRs) GRs[[x]]), "GRangesList")))
    }, simplify = FALSE)
    filename <- sprintf("Data/release/PrimingRate/GRs%s.RDS", x)
    message(x, " > ", filename)
    saveRDS(GRs_merged, file = filename)
}
###########################################################################
## Filter mouse samples (biological and virtual) to get rid of human contaminants.
## Note, these alignments are unique reads aligned to either species
###########################################################################
contamReadIDs <- readRDS("Data/release/AlignmentArtifacts/Contam/primaryNoDupContamReadIDs.RDS")
sharedReadIDs <- readRDS("Data/release/AlignmentArtifacts/SharedReads/primaryNoDupSharedReadIDs.RDS")
rescuedReadIDs_ge20_lt20 <- readRDS("Data/release/AlignmentArtifacts/MapQual/rescuedPrimaryNoDupMapQualReadIDs_ge20_lt30_le0.1.RDS")
rescuedReadIDs_ge30 <- readRDS("Data/release/AlignmentArtifacts/MapQual/rescuedPrimaryNoDupMapQualReadIDs_ge30_le0.1.RDS")

mapqThs <- c("ge30_le0.1", "ge20_le0.1", "ge20_le0.1_strict", "ge10_le0.1")
for (mapqTh in c(mapqThs)) {
    if (mapqTh != "ge20_le0.1_strict") {
        goodMapQualReadIDs <- readRDS(file = sprintf("Data/release/AlignmentArtifacts/MapQual/goodPrimaryNoDupMapQualReadIDs_%s.RDS", mapqTh))
        GRsQualOutInPairListFiltered <- sapply(qualOutInPairs, function(qualOutInPair) {
            GRsList <- GRsQualOutInPairList[[qualOutInPair]]
            sapply(sampleIDs, function(sampleID) {
                message(mapqTh, " ", qualOutInPair, " ", sampleID)
                GRs <- GRsList[[sampleID]]
                if (length(GRs) > 0) {
                    readNames <- mcols(GRs)[, "name"]
                    readIDs <- sapply(strsplit(readNames, "-"), "[", 1)
                    contam <- contamReadIDs[[sampleID]]
                    good <- goodMapQualReadIDs[[sampleID]]
                    shared <- sharedReadIDs[[sampleID]]
                    idx <- (!readIDs %in% union(contam, shared)) & (readIDs %in% good)
                    return(GRs[idx])
                } else {
                    return(GRs)
                }
            }, simplify = FALSE)
        }, simplify = FALSE)
        for (qualOutInPair in qualOutInPairs) {
            GRs <- GRsQualOutInPairListFiltered[[qualOutInPair]]
            filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
            saveRDS(GRs, file = filename)
        }
    } else {
        goodMapQualReadIDs_30 <- readRDS(file = sprintf("Data/release/AlignmentArtifacts/MapQual/goodPrimaryNoDupMapQualReadIDs_%s.RDS", "ge30_le0.1"))
        goodMapQualReadIDs_20 <- readRDS(file = sprintf("Data/release/AlignmentArtifacts/MapQual/goodPrimaryNoDupMapQualReadIDs_%s.RDS", "ge20_le0.1"))
        GRsQualOutInPairListFiltered <- sapply(qualOutInPairs, function(qualOutInPair) {
            GRsList <- GRsQualOutInPairList[[qualOutInPair]]
            sapply(sampleIDs, function(sampleID) {
                message(mapqTh, " ", qualOutInPair, " ", sampleID)
                GRs <- GRsList[[sampleID]]
                if (length(GRs) > 0) {
                    readNames <- mcols(GRs)[, "name"]
                    readIDs <- sapply(strsplit(readNames, "-"), "[", 1)
                    contam <- contamReadIDs[[sampleID]]
                    good_30 <- goodMapQualReadIDs_30[[sampleID]]
                    good_20 <- goodMapQualReadIDs_20[[sampleID]]
                    good_20 <- good_20[!good_20 %in% good_30]
                    shared <- sharedReadIDs[[sampleID]]
                    rescued_30 <- rescuedReadIDs_ge30[[sampleID]]
                    rescued_20 <- rescuedReadIDs_ge20_lt20[[sampleID]]
                    good_20_rescued <- good_20[good_20 %in% rescued_20]
                    good_30_rescued <- good_30[good_30 %in% rescued_30]
                    good_rescued <- union(good_30_rescued, good_20_rescued)
                    idx <- (!readIDs %in% union(contam, shared)) & (readIDs %in% good_rescued)
                    return(GRs[idx])
                } else {
                    return(GRs)
                }
            }, simplify = FALSE)
        }, simplify = FALSE)
        for (qualOutInPair in qualOutInPairs) {
            GRs <- GRsQualOutInPairListFiltered[[qualOutInPair]]
            filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
            saveRDS(GRs, file = filename)
        }
    }
}
###########################################################################
## merge Aread5End Bread5End to form ABread5End
## merge Aread5End Bread5End and Cmate5End to form ABreadCmate5End
###########################################################################
for (mapqTh in mapqThs) {
    for (x in names(qualOutInPairs_merge_map)) {
        GRs_each <- sapply(qualOutInPairs_merge_map[[x]], function(y) {
            filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", y, mapqTh)
            message(mapqTh, " ", x, " < ", y, " ", filename)
            readRDS(filename)
        }, simplify = FALSE)
        GRs_merged <- sapply(sampleIDs, function(x) {
            unname(unlist(as(lapply(GRs_each, function(GRs) GRs[[x]]), "GRangesList")))
        }, simplify = FALSE)
        filename <- sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", x, mapqTh)
        message(mapqTh, " ", x, " > ", filename)
        saveRDS(GRs_merged, file = filename)
    }
}
