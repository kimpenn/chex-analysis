## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
#######################################################################
if (!exists("NGS") || is.environment(NGS)) { 
    NGS <- new.env(parent = emptyenv())
}

local({
    .VERSION = "0.1"
    library("ggplot2")

###########################################################################
## DNA sequences
###########################################################################
    getRC <- function(strings) {
        stringsC <- chartr("atgcrymkswATGCRYMKSW", "tacgyrkmswTACGYRKMSW", strings)
        charsC <- strsplit(stringsC, "")
        charsRC <- lapply(charsC, function(x) rev(x))
        stringsRC <- sapply(charsRC, function(x) paste0(x, collapse = ""))
    }

###########################################################################
## Primers and barcodes
###########################################################################
    filterReadPrimerLocs <- function(ReadPrimerLocs, bcPrimer = "2p", bcLength = 20, endType = "PE") {
        nbPrimer <- ifelse(bcPrimer == "2p", "pC", "2p")
        bcCols <- paste(c("R1", "R2"), bcPrimer, "fw", sep = "_")
        nbCols <- paste(c("R1", "R2"), nbPrimer, "fw", sep = "_")
        bcIdxR1 <- ReadPrimerLocs[, bcCols[1]] >= 6 & ReadPrimerLocs[, bcCols[1]] <= bcLength
        bcIdxR2 <- ReadPrimerLocs[, bcCols[2]] >= 6 & ReadPrimerLocs[, bcCols[2]] <= bcLength
        if (endType == "PE") {
            nbLocsR2 <- ReadPrimerLocs[bcIdxR1, nbCols[2]]
            nbLocsR1 <- ReadPrimerLocs[bcIdxR2, nbCols[1]]
            c(nbLocsR2, nbLocsR1)
        } else if (endType == "SE") {
            nbLocsR1 <- ReadPrimerLocs[, nbCols[1]]
        }
    }
    
###########################################################################
## STAR alignMatesGapMax  option is varied and we compare resultant insert sizes
###########################################################################
    getMatesGapMaxInclusionRates <- function(sampleID, MatesGapMaxLimits, baseDir = "Data/E.chex/6bpPrefixBC1TrimmedIsize800", alignedDir = "cutadaptSTARIsize") {
        n <- length(MatesGapMaxLimits)
        Mat <- matrix(NA, ncol = n, nrow = n, dimnames = list(MatesGapMaxLimits, MatesGapMaxLimits))
        sampleName <- paste0("Sample_", sampleID)
        for (i in 1:n) {
            maxLimit <- MatesGapMaxLimits[i]
            starDir <- file.path(baseDir, sampleName, paste0(alignedDir, maxLimit))
            for (j in 1:i) {
                limit <- MatesGapMaxLimits[j]
                fname <- paste0(sampleName, ".star.primaryR1MatesGapMax", limit, "IncludedLineNum.txt")
                filename <- file.path(starDir, fname)
                lineNum <- as.numeric(strsplit(readLines(filename, n = 1), " ")[[1]][1])
                Mat[j, i] <- lineNum
            }
        }
        Mat
    }

    plotMatesGapMaxInclusionRates <- function(Mat, sampleID, MatesGapMaxLimits) {
        colnames(Mat) <- rownames(Mat) <- NULL
        Df <- melt(Mat)
        colnames(Df) <- c("UpperLimit", "LowerLimit", "InclusionRate")
        ggplot(Df, aes(x = UpperLimit, y = LowerLimit, fill = InclusionRate)) + geom_tile() + scale_x_continuous(breaks = seq(MatesGapMaxLimits), labels = MatesGapMaxLimits) + scale_y_continuous(breaks = seq(MatesGapMaxLimits), labels = MatesGapMaxLimits) + ggtitle(sampleID)
    }

    plotIsizeMatesGapMaxLeJaccard <- function(Mat, sampleID, MatesGapMaxLimits) {
        colnames(Mat) <- rownames(Mat) <- NULL
        Df <- melt(Mat)
        colnames(Df) <- c("MatesGapMaxLimit", "InsertSizeLimit", "Jaccard")
        ggplot(Df, aes(x = MatesGapMaxLimit, y = InsertSizeLimit, fill = Jaccard)) + geom_tile() + scale_x_continuous(breaks = seq(MatesGapMaxLimits), labels = MatesGapMaxLimits) + scale_y_continuous(breaks = seq(MatesGapMaxLimits), labels = MatesGapMaxLimits) + ggtitle(sampleID)
    }

    getIsizeMatesGapMaxLeOverlap <- function(sampleID, MatesGapMaxLimits, baseDir = "Data/E.chex/6bpPrefixBC1TrimmedIsize800", alignedDir = "cutadaptSTARIsize", overlapType = c("Intersection", "Union")) {
        n <- length(MatesGapMaxLimits)
        overlapType <- match.arg(overlapType)
        Mat <- matrix(NA, ncol = n, nrow = n, dimnames = list(MatesGapMaxLimits, MatesGapMaxLimits))
        sampleName <- paste0("Sample_", sampleID)
        for (i in 1:n) {
            maxLimit <- MatesGapMaxLimits[i]
            starDir <- file.path(baseDir, sampleName, paste0(alignedDir, maxLimit))
            for (j in 1:i) {
                limit <- MatesGapMaxLimits[j]
                fname <- paste0(sampleName, ".star.primaryR1IsizeMatesGapMaxLe", limit, overlapType, "LineNum.txt")
                filename <- file.path(starDir, fname)
                lineNum <- as.numeric(strsplit(readLines(filename, n = 1), " ")[[1]][1])
                Mat[j, i] <- lineNum
            }
        }
        Mat
    }

    plotIsizesByMatesGapMax <- function(IsizesDf, sampleID = "", alpha = 1, size = 1, breaks = c(1, seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), xLab = "Insert size inferred by STAR", yLab = "# of primarily aligned read pairs") {
        Df <- with(IsizesDf, table(MatesGapMax, abs(Isize)))
        Df <- data.frame(Df)
        colnames(Df) <- c("MatesGapMax", "Isize", "Freq")
        Df[, "MatesGapMax"] <- ordered(Df[, "MatesGapMax"], levels = sort(as.numeric(levels(Df[, "MatesGapMax"])), decreasing = TRUE))
        Df[, "Isize"] <- as.numeric(Df[, "Isize"])
        ggplot2::ggplot(Df, aes(x = Isize, y = Freq, group = MatesGapMax, col = MatesGapMax)) + geom_line(alpha = alpha, size = size) + scale_x_continuous(trans = "log10", breaks = breaks) + scale_y_continuous(trans = "log10") + ggtitle(sampleID) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab(xLab) + ylab(yLab)
    }

###########################################################################
## RNA-seq normalization methods
###########################################################################
    RLE <- function(counts, locfunc = median) {
        loggeomeans <- rowMeans(log(counts))
        if (all(is.infinite(loggeomeans))) {
            warning("No genes without one zero")
            scaleFactors <- rep(NA, ncol(counts))
        } else {
            scaleFactors <- apply(counts, 2, function(count)
              exp(locfunc((log(count) - loggeomeans)[is.finite(loggeomeans) & count > 0])))
        }
        num_nonzero <- apply(counts, 2, function(count) sum(count != 0))
        normalized <- t(t(counts) / scaleFactors)
        attr(normalized, "NumNonzeros") <- num_nonzero
        attr(normalized, "ScaleFactors") <- scaleFactors
        normalized
    }

    ## Modified RLE
    ## @threshold: the minimum allowed percentage of samples with nonzero expression
    ##             when equal to 1.0, the results fall back to original RLE
    RLE2 <- function(counts, locfunc = median, geoMeans, controlGenes, threshold = 1.0) {
        nSamples <- ncol(counts)
        nGenes <- nrow(counts)
        nSamplesTh <- threshold * ncol(counts)
        if (missing(geoMeans)) {
            loggeomeans <- apply(counts, 1, function(crow) {
                if (sum(crow != 0) < nSamplesTh) { -Inf } 
                else { mean(log(crow[crow != 0])) }
            })
        } else {
            if (length(geoMeans) != nGenes) {
                stop("geoMeans should be as long as the number of rows of counts")
            }
            loggeomeans <- log(geoMeans)
        }
        if (all(is.infinite(loggeomeans))) {
            stop("every gene contains at least one zero, cannot compute log geometric means")
        }
        isFinite <- is.finite(loggeomeans)
        if (missing(controlGenes)) { controlGenes <- seq(nGenes) }
        else {
            if (!(is.numeric(controlGenes) | is.logical(controlGenes) | is.character(controlGenes))) {
                stop("controlGenes should be either a numeric or logical vector")
            }
        }
        isFiniteSub <- isFinite[controlGenes]
        loggeomeansSub <- loggeomeans[controlGenes]
        countsSub <- counts[controlGenes, ]
        nFiniteGenesSub <- sum(isFiniteSub)
        scaleFactors <- apply(countsSub, 2, function(cnts) {
            exp(locfunc((log(cnts) - loggeomeansSub)[isFiniteSub & cnts > 0]))
        })
        normalized <- t(t(counts)/scaleFactors)
        attr(normalized, "ScaleFactors") <- scaleFactors
        attr(normalized, "NumNonzeros") <- nFiniteGenesSub
        normalized
    }

    for (.obj in ls()) assign(.obj, get(.obj), envir = NGS)
})
