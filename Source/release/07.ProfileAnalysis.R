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
SampleInfofull <- subset(SampleInfoFull, IsOut == "N")
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[, "SampleID"]
sampleIDsVirtual <- subset(SampleInfoFull, CompType == "Virtual")[, "SampleID"]
sourceIDsVirtual <- strsplit(SampleInfoFull[sampleIDsVirtual, "SourceIDs"], ',')
names(sourceIDsVirtual) <- sampleIDsVirtual
sampleIDs <- subset(SampleInfoFull, CompType == "Biol")[, "SampleID"]
SampleInfo <- SampleInfoFull[sampleIDs, ]

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
trtGroups <- c(
    "NegCtrlNone", "NegCtrlSingle", "NegCtrlPooled", "NegCtrlNil", "NegCtrlMerged", "NegCtrlAll", 
    "PositiveNone", "PositiveSingle", "PositivePooled", "PositiveNil", "PositiveMerged", "PositiveAll"
)
bioGroupCols <- structure(if (length(bioGroups) > 9) { colorRampPalette(brewer.pal(n = 9, "Set1"))(length(bioGroups)) } else { brewer.pal(n = length(bioGroups), "Set1") }, names = bioGroups)
trtGroupCols <- structure(c(brewer.pal(n = 6, "Reds"), brewer.pal(n = 6, "Greys")), names = trtGroups)
trtGroupPchs <- structure(c(2, 1, 0, 3, 4, 5, 17, 19, 15, 12, 13, 11), names = trtGroups)
annotationColors <- list(BioGroup = bioGroupCols, TrtGroup = trtGroupCols)

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
    FeatureIDType = c("gene_id", "gene_id", "name", rep("", 5), "name", rep("", 5), "name"), 
    ncores = c(59, 59, 59, rep(59, 6), rep(59, 6)), 
    xat = c(list(c(0, 250, 500)), replicate(2, list(c(0, 250, 1000, 1250))), replicate(6, list(c(0, 25, 100, 125))), replicate(6, list(c(0, 150, 225, 375)))),
    xTickLabel = c(list(c("-5kb", "TSS", "5kb")), list(c("-5kb", "Gene start", "Gene end", "5kb")), list(c("-5kb", "Intergenic start", "Intergenic end", "5kb")), lapply(c("Exon", "Intron", "5'UTR", "3'UTR", "CDS", "CGI"), function(l) c("-500bp", paste(l, "start"), paste(l, "end"), "500bp")), lapply(c("Exon", "Intron", "5'UTR", "3'UTR", "CDS", "CGI"), function(l) c("-3kb", paste(l, "start"), paste(l, "end"), "3kb"))), 
    ylab = c(replicate(2, list("Average coverage per gene")), "Average coverage per intergenic region", lapply(c("Exon", "Intron", "5'UTR", "3'UTR", "CDS", "CGI"), function(l) { paste("Average coverage per", l) }), lapply(c("Exon", "Intron", "5'UTR", "3'UTR", "CDS", "CGI"), function(l) { paste("Average coverage per", l) })), 
    rankCols = c(list(c(150, 350)), list(c(200, 1050)), list(c(200, 1050)), replicate(6, list(c(1, 125))), replicate(6, list(c(100, 275)))), 
    xlab = c(
        list(c("-5kb", rep("", 248), "TSS", rep("", 249), "5kb")), 
        list(c("-5kb", rep("", 249), "Gene start", rep("", 748), "Gene end", rep("", 249), "5kb")), 
        list(c("-5kb", rep("", 249), "Intergenic start", rep("", 748), "Intergenic end", rep("", 249), "5kb")),
        list(c("-500bp", rep("", 24), "Exon start", rep("", 73), "Exon end", rep("", 24), "500bp")), 
        list(c("-500bp", rep("", 24), "Intron start", rep("", 73), "Intron end", rep("", 24), "500bp")),
        list(c("-500bp", rep("", 24), "5'UTR start", rep("", 73), "5'UTR end", rep("", 24), "500bp")), 
        list(c("-500bp", rep("", 24), "3'UTR start", rep("", 73), "3'UTR end", rep("", 24), "500bp")), 
        list(c("-500bp", rep("", 24), "CDS start", rep("", 73), "CDS end", rep("", 24), "500bp")), 
        list(c("-500bp", rep("", 24), "CGI start", rep("", 73), "CGI end", rep("", 24), "500bp")), 
        list(c("-3kb", rep("", 149), "Exon start", rep("", 73), "Exon end", rep("", 149), "3kb")),
        list(c("-3kb", rep("", 149), "Intron start", rep("", 73), "Intron end", rep("", 149), "3kb")),
        list(c("-3kb", rep("", 149), "5'UTR start", rep("", 73), "5'UTR end", rep("", 149), "3kb")), 
        list(c("-3kb", rep("", 149), "3'UTR start", rep("", 73), "3'UTR end", rep("", 149), "3kb")),
        list(c("-3kb", rep("", 149), "CDS start", rep("", 73), "CDS end", rep("", 149), "3kb")), 
        list(c("-3kb", rep("", 149), "CGI start", rep("", 73), "CGI end", rep("", 149), "3kb"))
        ),
)

###########################################################################
## Load Avg and SEM matrices
###########################################################################
inDir <- "Report/release/Profiles"
CvgsABreadCmate5EndExtFullList <- list(Avg = vector("list", nrow(CvgsPerFeatureConfig)), SEM = vector("list", nrow(CvgsPerFeatureConfig)))
for (i in 1:nrow(CvgsPerFeatureConfig)) {
    with(CvgsPerFeatureConfig[i, ], message(DB, " ", FeatureType, " ", nbinsUpstream, " ", nbinsDownstream, " ", nbinsBody))
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    if (nbinsBody == "") {
        AvgFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s.csv", inDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
        SEMFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s.csv", inDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
    } else {
        AvgFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s.csv", inDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
        SEMFilename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredSEM%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s.csv", inDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
    }
    CvgsABreadCmate5EndExtFullList[["Avg"]][[i]] <- as.matrix(read.csv(AvgFilename, as.is = TRUE, check.names = FALSE, row.names = 1))
    CvgsABreadCmate5EndExtFullList[["Avg"]][[i]] <- cbind(CvgsABreadCmate5EndExtFullList[["Avg"]][[i]], sapply(sourceIDsVirtual, function(SIDs) rowSums(CvgsABreadCmate5EndExtFullList[["Avg"]][[i]][, SIDs, drop = FALSE])))
    CvgsABreadCmate5EndExtFullList[["SEM"]][[i]] <- as.matrix(read.csv(SEMFilename, as.is = TRUE, check.names = FALSE, row.names = 1))
    CvgsABreadCmate5EndExtFullList[["SEM"]][[i]] <- cbind(CvgsABreadCmate5EndExtFullList[["SEM"]][[i]], sapply(sourceIDsVirtual, function(SIDs) rowMeans(CvgsABreadCmate5EndExtFullList[["SEM"]][[i]][, SIDs, drop = FALSE])))
}
 
###########################################################################
## Coverage heatmap averaging features
###########################################################################
## Scale profiles to unit SD and zero mean.
CvgsABreadCmate5EndExtFullListAvgZscore <- lapply(CvgsABreadCmate5EndExtFullList[["Avg"]], scale)
outDir <- "Report/release/ProfileAnalysis/zscoreHeatmap"
dir.create(outDir, FALSE, TRUE)

range(CvgsABreadCmate5EndExtFullListAvgZscore[[14]], na.rm = TRUE)
## [1] -3.786245 19.313277
range(CvgsABreadCmate5EndExtFullListAvgZscore[[1]], na.rm = TRUE)
## [1] -3.047219  4.628826
breaks <- c(
    min(CvgsABreadCmate5EndExtFullListAvgZscore[[14]], na.rm = TRUE), 
    seq(min(CvgsABreadCmate5EndExtFullListAvgZscore[[1]], na.rm = TRUE),
        max(CvgsABreadCmate5EndExtFullListAvgZscore[[1]], na.rm = TRUE),
        length.out = 256 - 2
    ),
    max(CvgsABreadCmate5EndExtFullListAvgZscore[[14]], na.rm = TRUE)
)
pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(255)

for (i in 1:nrow(CvgsPerFeatureConfig)) {
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    xat <- config[["xat"]][[1]]
    xTickLabel <- config[["xTickLabel"]][[1]]
    ylab <- config[["ylab"]]
    xlab <- config[["xlab"]][[1]]
    M <- CvgsABreadCmate5EndExtFullListAvgZscore[[i]] 
    M <- apply(M, 2, function(x) if (all(is.na(x))) return(rep(0, length(x))) else return(x))
    rownames(M) <- xlab
    for (bioGroup in bioGroups[1]) {
        SIDs <- subset(SampleInfoFull, BioGroup == bioGroup & IsNegCtrl == "N" & IsOut == "N")[, "SampleID"]
        n <- length(SIDs)
        m <- t(M)[SIDs, , drop = FALSE]
        message(featureTypeLab, " ", extType, " ", upstreamLab, " ", downstreamLab, " ", nbinsUpstream, " ", nbinsDownstream, " ", nbinsBody, " ", bioGroup)
        if (nbinsBody == "") {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, bioGroup)
        } else {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody, bioGroup)
        }
        pdf(filename, width = 10, height = max(1/6 * n, 3))
        if (!all(m == 0)) {
            pheatmap(m, cluster_row = length(SIDs) > 1, cluster_col = FALSE, show_colnames = TRUE, annotation_row = SampleInfoFull[SIDs, c("TrtGroup"), drop = FALSE], annotation_colors = annotationColors, scale = "none", main = bioGroup, labels_row = paste0("E.", SampleInfoFull[SIDs, "ExptID"], " #", sub("scCLTdegenNuc", "", SIDs)), treeheight_row = 0, treeheight_col = 0)
        } else {
            pheatmap(m, cluster_row = length(SIDs) > 1, cluster_col = FALSE, show_colnames = TRUE, annotation_row = SampleInfoFull[SIDs, c("TrtGroup"), drop = FALSE], annotation_colors = annotationColors, scale = "none", main = bioGroup, labels_row = paste0("E.", SampleInfoFull[SIDs, "ExptID"], " #", sub("scCLTdegenNuc", "", SIDs)), breaks = c(-0.01, 0.01), treeheight_row = 0, treeheight_col = 0)
        }
        dev.off()
    }
}

###########################################################################
## Coverage profiles per sample
###########################################################################
SampleInfoFull$TrtGroup <- factor(SampleInfoFull$TrtGroup, levels = trtGroups)
outDir <- "Report/release/ProfileAnalysis/Curve"
dir.create(outDir, FALSE, TRUE)
for (i in 1:nrow(CvgsPerFeatureConfig)) {
    message(i)
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    xat <- config[["xat"]][[1]]
    xTickLabel <- config[["xTickLabel"]][[1]]
    ylab <- config[["ylab"]]
    xlab <- config[["xlab"]][[1]]
    M <- CvgsABreadCmate5EndExtFullList[["Avg"]][[i]]
    M <- apply(M, 2, function(x) if (all(is.na(x))) return(rep(0, length(x))) else return(x))
    for (bioGroup in bioGroups) {
        SIDs <- subset(SampleInfoFull, BioGroup == bioGroup)[["SampleID"]]
        idx <- with(SampleInfoFull[SIDs, ], order(ProbeType, TrtGroup, ExptID, SampleID))
        n <- length(SIDs)
        if (nbinsBody == "") {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, bioGroup)
        } else {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredAvg%sMainNoMY_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody, bioGroup)
        }
        pdf(filename, width = 8, height = 3.5)
        for (SID in SIDs[idx]) {
            bioGroup <- SampleInfoFull[SID, "BioGroup"]
            trtGroup <- SampleInfoFull[SID, "TrtGroup"]
            exptID <- SampleInfoFull[SID, "ExptID"]
            sampleName <- SampleInfoFull[SID, "SampleName"]
            x <- M[, SID]
                Graphics$plotProfile(y = x, xat = xat, xTickLabels = xTickLabel, main = sprintf("%s %s E.%s %s", bioGroup, trtGroup, exptID, SID), ylab = ylab)
        }
        dev.off()
    }
}

###########################################################################
## Coverage curves and matrix
###########################################################################
inDirbase <- "Report/release/Profiles"
outDir <- "Report/release/ProfileAnalysis/ProfileAndMatrix"
FeatureIDsMainNoMY <- readRDS("Data/release/GenomicFeatures/FeatureIDsMainNoMY.RDS")
dir.create(outDir, FALSE, TRUE)
for (i in 1:nrow(CvgsPerFeatureConfig)) {
    config <- CvgsPerFeatureConfig[i, ]
    featureType <- config[["FeatureType"]]
    featureTypeLab <- config[["FeatureTypeLab"]]
    extType <- config[["ExtType"]]
    upstream <- config[["upstream"]]
    downstream <- config[["downstream"]]
    upstreamLab <- config[["upstreamLab"]]
    downstreamLab <- config[["downstreamLab"]]
    nbinsUpstream <- config[["nbinsUpstream"]]
    nbinsDownstream <- config[["nbinsDownstream"]]
    nbinsBody <- config[["nbinsBody"]]
    xat <- config[["xat"]][[1]]
    xTickLabel <- config[["xTickLabel"]][[1]]
    ylab <- config[["ylab"]]
    featureIDsMainNoMY <- config[["FeatureIDsMainNoMY"]]
    rankCols = config[["rankCols"]][[1]]
    if (nbinsBody == "") {
        inDir <- sprintf("%s/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%s", inDirbase, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream)
    } else {
        inDir <- sprintf("%s/CvgsABreadCmate5EndExtFilteredPer%s_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s", inDirbase, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody)
    }
    Avg <- CvgsABreadCmate5EndExtFullList[["Avg"]][[i]] 
    SEM <- CvgsABreadCmate5EndExtFullList[["SEM"]][[i]]
    for (bioGroup in bioGroups[1]) {
        message(featureTypeLab, " ", extType, " ", upstream, " ", downstream, " ", nbinsUpstream, " ", nbinsDownstream, " ", nbinsBody, " ", bioGroup)
        if (nbinsBody == "") {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredPer%sMainNoMYtop5k_Up%sDn%snbinsUp%snbinsDn%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, bioGroup)
        } else {
            filename <- sprintf("%s/CvgsABreadCmate5EndExtFilteredPer%sMainNoMYtop5k_Up%sDn%snbinsUp%snbinsDn%snbinsBody%s_%s.pdf", outDir, featureTypeLab, upstreamLab, downstreamLab, nbinsUpstream, nbinsDownstream, nbinsBody, bioGroup)
        }
        pdf(filename, width = 2, height = 4)
        sampleIDVirtmax <- tail(subset(SampleInfoFull, BioGroup == bioGroup & IsNegCtrl == "N" & IsOut == "N")[, "SampleID"], 1)
        species <- SampleInfoFull[sampleIDVirtmax, "Species"]
        featureIDs <- FeatureIDsMainNoMY[[species]][[featureIDsMainNoMY]]
        if (SampleInfoFull[sampleIDVirtmax, "CompType"] == "Biol") {
            cvgFilename <- sprintf("%s/Sample_%s.csv.gz", inDir, sampleIDVirtmax)
            M <- read.csv(cvgFilename, as.is = TRUE, check.names = FALSE)
            rids <- as.matrix(M[, 1])
            M <- as.matrix(M[, -1])
            M <- M[rids %in% featureIDs, ]
            avg <- as.numeric(Avg[, sampleIDVirtmax])
            sem <- as.numeric(SEM[, sampleIDVirtmax])
        } else {
            SIDs <- strsplit(SampleInfoFull[sampleIDVirtmax, "SourceIDsNoOut"], ",")[[1]]
            cvgs <- sapply(SIDs, function(SID) {
                cvgFilename <- sprintf("%s/Sample_%s.csv.gz", inDir, SID)
                message(cvgFilename)
                M <- read.csv(cvgFilename, as.is = TRUE, check.names = FALSE)
                rids <- M[, 1]
                M <- as.matrix(M[, -1])
                M <- M[rids %in% featureIDs, ]
                Matrix::Matrix(M, sparse = TRUE)
            }, simplify = FALSE)
            M <- Reduce("+", cvgs)
            avg <- Matrix::colMeans(M, na.rm = TRUE)
            sem <- apply(M, 2, function(x) { sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))) })
        }
        Graphics$plotProfileAndMatrix(M, y = avg, ci = sem,  xat = xat, xTickLabels = xTickLabel, ylab = ylab, breaks = seq(quantile(M, 0.001, na.rm = TRUE), quantile(M, 0.999, na.rm = TRUE), length.out = 256), topRow = 500, cex.axis = 1, panelRatio = c(1, 4), rankCols = rankCols, main = "", xlas = ifelse(i > 1, 2, 1), cex.main = 1)
		M1 <- M[order(rowSums(M[, 150:350]), decreasing = TRUE), ]
		m <- M1[1:581, ]
		avg <- colMeans(m, na.rm = TRUE)
		sem <- apply(m, 2, function(x) { sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))) })
		Graphics$plotProfileAndMatrix(log10(1+m), y = avg, ci = sem,  xat = xat, xTickLabels = xTickLabel, ylab = ylab, breaks = seq(quantile(log10(1+m), 0.0, na.rm = TRUE), quantile(log10(1+m), 1, na.rm = TRUE), length.out = 256), topRow = 581, cex.axis = 1, panelRatio = c(1, 4), rankCols = rankCols, main = "", xlas = ifelse(i > 1, 2, 1), cex.main = 1)
        dev.off()
    }
}

###########################################################################
## Negative controls' TSS coverage profiles (z-scores)
###########################################################################
sampleIDs_K562Positive <- subset(SampleInfoFull, BioGroup == "K562" & IsNegCtrl == "N" & CompType == "Biol")[, "SampleID"]
sampleIDs_K562NegCtrl <- subset(SampleInfoFull, BioGroup == "K562" & IsNegCtrl == "Y" & CompType == "Biol")[, "SampleID"]
sampleIDs_K562MungBean <- subset(SampleInfoFull, BioGroup == "K562MungBean" & IsNegCtrl == "Y" & CompType == "Biol")[, "SampleID"]
sampleID_virtmax_K562Positive <- "K562PositiveAll"
sampleID_virtmax_K562NegCtrl <- "K562NegCtrlAll"
sampleID_virtmax_K562MungBean <- "K562MungBeanPositiveMerged"
M <- scale(CvgsABreadCmate5EndExtFullList$Avg[[1]])
pdf("Report/release/ProfileAnalysis/zscoreCurve/CvgsABreadCmate5EndExtFiltered_K562_Positive-NegCtrl-MungBean.pdf", width = 3, height = 6, useDingbats = FALSE)
par(mfcol = c(3, 1), mar = c(2, 2, 1, 1))
matplot(x = 1:500, y = M[, sampleIDs_K562Positive], lty = 1, lwd = 0.1, col = "black", axes = FALSE, type = "l", ylim = c(-3, 3), main = "Positive")
lines(M[, sampleID_virtmax_K562Positive], lty = 1, lwd = 1, col = "black")
axis(at = c(0, 250, 500), side = 1, labels = c("-5kb", "TSS", "5kb"), las = 1)
axis(side = 2)
matplot(x = 1:500, y = M[, sampleIDs_K562NegCtrl], lty = 1, lwd = 0.1, col = "red", axes = FALSE, type = "l", ylim = c(-3, 3), main = "NegCtrl")
lines(M[, sampleID_virtmax_K562NegCtrl], lty = 1, lwd = 1, col = "red")
axis(at = c(0, 250, 500), side = 1, labels = c("-5kb", "TSS", "5kb"), las = 1)
axis(side = 2)
matplot(x = 1:500, y = M[, sampleIDs_K562MungBean], lty = 1, lwd = 0.1, col = "purple", axes = FALSE, type = "l", ylim = c(-3, 3), main = "MungBean")
lines(M[, sampleID_virtmax_K562MungBean], lty = 1, lwd = 1, col = "purple")
axis(at = c(0, 250, 500), side = 1, labels = c("-5kb", "TSS", "5kb"), las = 1)
axis(side = 2)
dev.off()
