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
FeatureIDsMainNoMY <- readRDS("Data/GenomicFeatures/FeatureIDsMainNoMY.RDS")
UCSCCpGIslands <- readRDS("Data/GenomicFeatures/UCSCCpGIslands.RDS")

SampleInfoFull <- read.csv("Data/SampleInfoFull20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[["SampleID"]]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
sampleIDs <- SampleInfo[, "SampleID"]

## Load raw priming rate data after mouse filtered
GRsABreadCmate5End <- readRDS("Data/PrimingRate/GRsABreadCmate5EndFiltered_ge20_le0.1_strict.RDS")
## Extend 1kb toward both sides (so that the point becomes 2kb region)
GRsABreadCmate5EndExt <- lapply(GRsABreadCmate5End, function(GRs) flank(GRs, width = 1000, start = TRUE, both = TRUE))
CvgsABreadCmate5EndExt <- lapply(GRsABreadCmate5EndExt, function(GRs) coverage(GRs))

###########################################################################
## For each biological and computationally merged sample, including negative samples.
###########################################################################
ncores <- 60
CvgsPerFeatureConfig <- tibble(
    DB = c(rep("EnsFeatures", 8), "UCSCCpGIslands", rep("EnsFeatures", 5), "UCSCCpGIslands"),
    FeatureTypeLab = c("TSS", "Gene", "Intergenic", rep(c("Exon", "Intron", "FiveUTR", "ThreeUTR", "CDS", "CpG"), 2)),
    FeatureType = c("Gene", "Gene", "Intergenic", rep(c("ExonByTranscript", "IntronByTranscript", "FiveUTRByTranscript", "ThreeUTRByTranscript", "CDSByTranscript", ""), 2)),
    FeatureIDsMainNoMY = c(rep("Gene", 2), "Intergenic", rep("Transcript", 5), "CpG", rep("Transcript", 5), "CpG"),
    ExtType = c("point", rep("range", 14)),
    nbinsUpstream = c(250, 250, 250, rep(25, 6), rep(150, 6)),
    nbinsDownstream = c(250, 250, 250, rep(25, 6), rep(150, 6)),
    nbinsBody = c("", 750, 750, rep(75, 6), rep(75, 6)),
    upstream = c(5000, 5000, 5000, rep(500, 6), rep(3000, 6)),
    downstream = c(5000, 5000, 5000, rep(500, 6), rep(3000, 6)),
    upstreamLab = c("5k", "5k", "5k", rep("500", 6), rep("3k", 6)),
    downstreamLab = c("5k", "5k", "5k", rep("500", 6), rep("3k", 6)),
    FeatureIDType = c("gene_id", "gene_id", "name", rep("", 5), "", rep("", 5), ""), 
    ncores = c(ncores, ncores, ncores, rep(ncores, 6), rep(ncores, 6)),
)

for (i in (1:nrow(CvgsPerFeatureConfig))[10:15]) {
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    db <- config[["DB"]]
    if (nbinsBody == "") { nbinsBody <- NULL } else { nbinsBody <- as.integer(nbinsBody) }
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    featureIDType <- config[["FeatureIDType"]]
    if (featureIDType == "") { featureIDType <- NULL }
    n <- config[["ncores"]]
    mclapply(sampleIDs, function(sampleID) {
        Cvgs <- CvgsABreadCmate5EndExt[[sampleID]]
        species <- SampleInfoFull[sampleID, "Species"]
        if (species == "none") { species <- "mouse" }
        if (featureType != "") {
            features <- get(db)[[featureType]][[species]]
        } else {
            features <- get(db)[[species]]
        }
        if (inherits(features, "GRangesList")) {
            features <- unique(unlist(features))
        }
        if (extType == "point") {
            outDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
        } else {
            outDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
        }
        dir.create(outDir, FALSE, TRUE)
        outFilename <- sprintf("%s/Sample_%s.csv.gz", outDir, sampleID)
        message(sampleID, " ", species, " ", outFilename)
        X <- Genome$getCvgsPerFeature(Cvgs, features = features, nbinsUpstream = nbinsUpstream, nbinsDownstream = nbinsDownstream, nbinsBody = nbinsBody, upstream = upstream, downstream = downstream, featureIDType = featureIDType, extType = extType, allFeatures = TRUE, ncores = 1)
        write.csv(as.matrix(X), file = gzfile(outFilename))
        return(NULL)
    }, mc.cores = n)
}

###########################################################################
## Remove chrM, chrY and non-major chromosomes and compute 
## (1) average, (2) SEM for each genomic feature type for each sample
###########################################################################
CvgsPerFeatureConfig[["ncores"]] <- c(rep(ncores, 3), rep(ceiling(ncores/2), 6), rep(ceiling(ncores/2), 6))
for (i in (1:nrow(CvgsPerFeatureConfig))[10:15]) {
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    db <- config[["DB"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    featureIDsMainNoMY <- config[["FeatureIDsMainNoMY"]]
    n <- config[["ncores"]]
    mclapply(sampleIDs, function(sampleID) {
        if (extType == "point") {
            inDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
            AvgDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
            SEMDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
        } else {
            inDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
            AvgDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
            SEMDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
        }
        dir.create(AvgDir, FALSE, TRUE)
        dir.create(SEMDir, FALSE, TRUE)
        inFilename <- sprintf("%s/Sample_%s.csv.gz", inDir, sampleID)
        AvgFilename <- sprintf("%s/Sample_%s.csv.gz", AvgDir, sampleID)
        SEMFilename <- sprintf("%s/Sample_%s.csv.gz", SEMDir, sampleID)
        message(featureTypeLab, " ", extType, " ", upstreamLab, " ", downstreamLab, " ", nbinsUpstream, " ", nbinsDownstream, " ", nbinsBody, " ", sampleID, " ")
        species <- SampleInfoFull[sampleID, "Species"]
        if (species == "none") { species <- "mouse" }
        featureIDs <- FeatureIDsMainNoMY[[species]][[featureIDsMainNoMY]]
        X <- read.csv(inFilename, as.is = TRUE, check.names = FALSE)
        IDs <- X[, 1]
        idx <- IDs %in% featureIDs
        Mat <- Matrix(as.matrix(X[idx, -1]), sparse = TRUE) # Remove rownames
        Avg <- Matrix::colMeans(Mat, na.rm = TRUE)
        SEM <- apply(Mat, 2, function(x) { sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))) })
        write.csv(Avg, file = gzfile(AvgFilename), row.names = TRUE)
        write.csv(SEM, file = gzfile(SEMFilename), row.names = TRUE)
    }, mc.cores = n)
}

###########################################################################
## Load individual Avg and SEM file and assemble them into a matrix. Rows are genomic locations in bin index, column sample IDs. 
###########################################################################
for (i in (1:nrow(CvgsPerFeatureConfig))[10:15]) {
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    db <- config[["DB"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    featureIDsMainNoMY <- config[["FeatureIDsMainNoMY"]]
    ncores <- config[["ncores"]]
    outDir <- sprintf("%s/Report/Profiles", .PWD)
    if (extType == "point") {
        AvgDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
        SEMDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
    } else {
        AvgDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
        SEMDir <- sprintf("%s/Report/Profiles/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", .PWD, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
    }
    AvgSEM <- lapply(sampleIDs, function(sampleID) {
        AvgFilename <- sprintf("%s/Sample_%s.csv.gz", AvgDir, sampleID)
        SEMFilename <- sprintf("%s/Sample_%s.csv.gz", SEMDir, sampleID)
        message(featureTypeLab, " ", extType, " ", upstreamLab, " ", downstreamLab, " ", nbinsUpstream, " ", nbinsDownstream, " ", nbinsBody, " ", sampleID, " ")
        Avg <- read.csv(AvgFilename, row.names = 1)
        SEM <- read.csv(SEMFilename, row.names = 1)
        list(Avg, SEM)
    })
    Avg <- do.call(cbind, lapply(AvgSEM, function(x) x[[1]])); colnames(Avg) <- sampleIDs
    SEM <- do.call(cbind, lapply(AvgSEM, function(x) x[[2]])); colnames(SEM) <- sampleIDs
    if (extType == "point") {
        AvgFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s.csv", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
        SEMFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s.csv", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
    } else {
        AvgFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s.csv", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
        SEMFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s.csv", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
    }
    write.csv(Avg, file = AvgFilename)
    write.csv(SEM, file = SEMFilename)
}
