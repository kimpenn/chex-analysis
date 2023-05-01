## Author: Youtao Lu <luyoutao@sas.upenn.edu>

## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
#######################################################################
if (!exists("Graphics") || is.environment(Graphics)) { 
    Graphics <- new.env(parent = emptyenv()) 
}

local({
    .VERSION = "0.4"
    library("ape")
    library("GenomicRanges")
    library("RColorBrewer")
    library("viridis")
    library("pheatmap")
    library("gridExtra")
    library("ggplot2")
    library("ggsignif")
    library("ggridges")

    ## Mitochondrial genes coverage
    plot_mito_coverage_stranded <- function(GRs, ann, ylab = "", main = "", rectcolor = gray.colors(10)[c(1, 6)], textcolor = c("red", "blue"), textcolor1 = c("red4", "red"), textcolor2 = c("midnightblue", "blue"), textfield = "symbol", textcex = 0.5, texthjust1 = 50, texthjust2 = -50, srt = 60, mar = c(5, 4, 4, 2) + 0.1) {
        ##  ann must have at least these columns
        ##                  seqnames start   end strand   symbol
        ##  ENSG00000210049     chrM   577   647      +    MT-TF
        ##  ENSG00000211459     chrM   648  1601      +  MT-RNR1
        ##  ENSG00000210077     chrM  1602  1670      +    MT-TV
        ##  ENSG00000210082     chrM  1671  3229      +  MT-RNR2
        ##  ENSG00000209082     chrM  3230  3304      +   MT-TL1
        ##  ENSG00000198888     chrM  3307  4262      +   MT-XXX
        ##  ENSG00000210100     chrM  4263  4331      +   MT-XXX
        par(mar = mar)
        GRs <- GRs[seqnames(GRs) == "chrM"]
        gr1 <- GRs[strand(GRs) == "+"]
        gr2 <- GRs[strand(GRs) == "-"]
        cvg1 <- as.vector(coverage(gr1)[["chrM"]])
        cvg2 <- as.vector(coverage(gr2)[["chrM"]])
        cvgs <- cbind(cvg1, -cvg2)

        range1 <- range(cvg1)
        range2 <- range(cvg2)
        range <- c(-max(range2), max(range1))

        ann <- subset(ann, seqnames == "chrM")
        ann1 <- subset(ann, strand == "+")
        ann2 <- subset(ann, strand == "-")
        ann1$textcolor <- ifelse(seq(nrow(ann1)) %% 2 == 1, textcolor1[1], textcolor1[2])
        ann2$textcolor <- ifelse(seq(nrow(ann2)) %% 2 == 1, textcolor2[1], textcolor2[2])
        ann <- rbind(ann1, ann2)[rownames(ann), ]

        ybase <- range[2] * 2.8
        yheight <- 0.04 * (range1[2] + range2[2]) 
        
        matplot(cvgs, type = "h", axes = FALSE, col = "black", lty = 1, ylab = "", ylim = c(range[1], range[2] * 3))
        axis(1)
        axis(2)
        rect(xleft = ann1[, "start"], xright = ann1[, "end"], ybottom = ybase, ytop = ybase + yheight, col = rectcolor, xpd = TRUE, border = NA)
        rect(xleft = ann2[, "start"], xright = ann2[, "end"], ybottom = ybase, ytop = ybase - yheight, col = rectcolor, xpd = TRUE, border = NA)
        title(main = main, ylab = ylab)

        idx1 <- seq(1, nrow(ann), by = 2)
        idx2 <- seq(2, nrow(ann), by = 2)
        text(x = -60, y = ybase + yheight, label = "H", col = textcolor[1], cex = textcex * 0.8, xpd = TRUE)
        text(x = -60, y = ybase - yheight, label = "L", col = textcolor[2], cex = textcex * 0.8, xpd = TRUE)

        text(x = 0.5 * (ann[idx1, "start"] + ann[idx1, "end"]) + texthjust1, y = ybase + yheight * 7, label = ann[idx1, textfield], xpd = TRUE, srt = 60, col = ann[idx1, "textcolor"], cex = textcex)
        text(x = 0.5 * (ann[idx2, "start"] + ann[idx2, "end"]) + texthjust2, y = ybase - yheight * 7, label = ann[idx2, textfield], xpd = TRUE, srt = 60, col = ann[idx2, "textcolor"], cex = textcex)

    }

    plot_mito_cnts_bin <- function(PF, ann, binsize = 50, xlab = "Coordinate (bp)", ylab = "Priming counts", main = "", sub = "", rectcolor = gray.colors(10)[c(1, 6)], textfield = "symbol", textcolor = c("red", "blue"), textcolor1 = c("red4", "red"), textcolor2 = c("midnightblue", "blue"), textcex = 0.5, vline = NULL, hline = NULL, vlinelty = NULL, hlinelty = NULL, vlinecol = "red", hlinecol = "green", srt = 60, mar = c(5, 4, 4, 2) + 0.1, texthjust1 = 5, texthjust2 = -5) {
        ## ann must have at least these columns
        ##                 seqnames start   end strand   symbol
        ## ENSG00000210049     chrM   577   647      +    MT-TF
        ## ENSG00000211459     chrM   648  1601      +  MT-RNR1
        ## ENSG00000210077     chrM  1602  1670      +    MT-TV
        ## ENSG00000210082     chrM  1671  3229      +  MT-RNR2
        ## ENSG00000209082     chrM  3230  3304      +   MT-TL1
        ## ENSG00000198888     chrM  3307  4262      +   MT-XXX
        ## ENSG00000210100     chrM  4263  4331      +   MT-XXX
        par(mar = mar)
        range <- range(PF)
        ann <- subset(ann, seqnames == "chrM")
        ann1 <- subset(ann, strand == "+")
        ann2 <- subset(ann, strand == "-")
        ann1$textcolor <- ifelse(seq(nrow(ann1)) %% 2 == 1, textcolor1[1], textcolor1[2])
        ann2$textcolor <- ifelse(seq(nrow(ann2)) %% 2 == 1, textcolor2[1], textcolor2[2])
        ann <- rbind(ann1, ann2)[rownames(ann), ]
        ybase <- range[2] * 1.00
        yheight <- 0.04 * (range[2] - range[1]) 
        
        plot(PF, type = "h", axes = FALSE, col = "black", lty = 1, xlim = c(0, length(PF)), ylab = ylab, xlab = xlab)
        axis(1, at = seq(0, as.integer(strsplit(tail(names(PFs_mito_perbin_mouse[[1]][, 1]), 1), "-")[[1]][2]), by = 5000) / binsize, labels = seq(0, as.integer(strsplit(tail(names(PFs_mito_perbin_mouse[[1]][, 1]), 1), "-")[[1]][2]), by = 5000))
        axis(2)
        rect(xleft = ann1[, "start"]/binsize, xright = ann1[, "end"]/binsize, ybottom = ybase, ytop = ybase + yheight, col = rectcolor, xpd = TRUE, border = NA)
        rect(xleft = ann2[, "start"]/binsize, xright = ann2[, "end"]/binsize, ybottom = ybase - yheight, ytop = ybase, col = rectcolor, xpd = TRUE, border = NA)
        title(main = main, sub = sub)

        idx1 <- seq(1, nrow(ann), by = 2)
        idx2 <- seq(2, nrow(ann), by = 2)
        text(x = -4, y = ybase + yheight/2, label = "H", col = textcolor[1], cex = textcex * 0.5, xpd = TRUE)
        text(x = -4, y = ybase - yheight/2, label = "L", col = textcolor[2], cex = textcex * 0.5, xpd = TRUE)

        text(x = 0.5 * (ann[idx1, "start"] + ann[idx1, "end"])/binsize + texthjust1, y = ybase + yheight * 5, label = ann[idx1, textfield], xpd = TRUE, srt = 60, col = ann[idx1, "textcolor"], cex = textcex)
        text(x = 0.5 * (ann[idx2, "start"] + ann[idx2, "end"])/binsize + texthjust2, y = ybase - yheight * 5, label = ann[idx2, textfield], xpd = TRUE, srt = 60, col = ann[idx2, "textcolor"], cex = textcex)
        if (!is.null(vline)) { abline(v = vline, col = vlinecol, lty = vlinelty) }
        if (!is.null(hline)) { abline(h = hline, col = hlinecol, lty = hlinelty) }
    }

    plotProfile <- function(y, x = seq(length(y)), ci = NULL, xaxis = TRUE, yaxis = TRUE, xat = c(0, 250, 500), xlab = "", xTickLabels = c("-5kb", "TSS", "5kb"), type = "l", yat = NULL, ylab = "", xlas = 2, ylas = 2, axes = FALSE, col = "black", lty = 1, lwd = 1, pch = 1, ciShade = "#EBEBEB", legend = FALSE, legendPos = "topright", legendLabel = NULL, legendColor = col, legendLty = lty, legendLwd = lwd, cex.axis = 1, main = "", cex.main = 1, ...) {
        n <- length(y)
        if (!is.null(ci)) {
            if (is.list(ci) && length(ci) == 2) {
                ciU <- ci[[1]]; ciL <- ci[[2]]
            } else {
                if (n != length(ci)) {
                    stop("ci and x have different length!")
                }
                ciU <- ciL <- ci
            }
            plot(x = x, y = y, col = col, lty = lty, lwd = lwd, xlab = xlab, ylab = ylab, axes = FALSE, type = type, panel.first = polygon(c(1:n, n:1), c(y + ciU, rev(y - ciL)), border = NA, col = ciShade), main = main, cex.main = cex.main, ...)
        } else {
            plot(x = x, y = y, col = col, lty = lty, lwd = lwd, xlab = xlab, ylab = ylab, axes = FALSE, type = type, main = main, cex.main = cex.main, ...)
        }
        if (xaxis) { axis(at = xat, side = 1, labels = xTickLabels, las = xlas, cex.axis = cex.axis) }
        if (yaxis) { axis(side = 2, cex.axis = cex.axis) }

        if (legend) {
            legend(legendPos, legend = legendLabel, col = legendColor, lty = legendLty, lwd = legendLwd)
        }
    }

    plotProfileAndMatrix <- function(X, y = NULL, rows = NULL, ordered = TRUE, ci = NULL, ciShade = "#EBEBEB", plotci = TRUE, xaxis = TRUE, yaxis = TRUE, xat = c(0, 250, 500), xlab = "", xTickLabels = c("-5kb", "TSS", "5kb"), yat = NULL, ylab = "", xlas = 1, ylas = 1, ycol = "black", type = "l", lty = 1, lwd = 1, nbreaks = 256, Xcol = colorRampPalette(colors = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"))(nbreaks-1), breaks = seq(quantile(X, 0.001, na.rm = TRUE), quantile(X, 0.999, na.rm = TRUE), length.out = nbreaks), useRaster = TRUE, panelRatio = c(1, 4), topRow = 1000, rankCols = c(200, 300), cex.axis = 1, main = "", cex.main = 1, ...) {
        xcol <- ncol(X); yrow <- nrow(X)
        if (is.null(y)) {
            #' we have to compute average profile ourselves
            y <- colMeans(X, na.rm = TRUE)
        }
        if (!is.null(ci)) {
            if (is.list(ci) && length(ci) == 2) {
                ciU <- ci[[1]]; ciL <- ci[[2]]
            } else {
                if (xcol != length(ci)) {
                    stop("ci and x have different number of bins!")
                }
                ciU <- ciL <- ci
            }
        } else {
            ciU <- ciL <- ci <- apply(X, 2, function(x) sd(x, na.rm = TRUE) / sqrt(yrow))
        }
        layoutMat <- matrix(rep(c(1, 2), panelRatio))
        layout(layoutMat)
        par(mar = c(2,2,0.5,0.5), oma = c(0.5, 0.5, 1, 0.5))
        par(mar = c(2,2,0.5,0.5), oma = c(0.5, 0.5, 1, 0.5))
        if (plotci) {
            plotProfile(y = y, ci = list(ciU, ciL), xaxis = FALSE, yaxis = TRUE, type = type, col = ycol, lwd = lwd, lty = lty, ylas = ylas, cex.axis = cex.axis, main = main, cex.main = cex.main)
        } else {
            plotProfile(y = y, ci = NULL, xaxis = FALSE, yaxis = TRUE, type = type, col = ycol, lwd = lwd, lty = lty, ylas = ylasi, cex.axis = cex.axis, main = main, cex.main = cex.main)
        }
        if (ordered) {
            X <- X[order(rowSums(X[, rankCols[1]:rankCols[2]], na.rm = TRUE), decreasing = TRUE), ]
        } else if (!is.null(rows)) {
            X <- X[rows, ]
        }
        X[is.na(X)] <- 0
        if (!is.null(topRow)) {
            X <- X[1:topRow, ]
        }
        plotMatrix(X = X, xaxis = TRUE, yaxis = TRUE, yat = 0, yTickLabels = "", breaks = breaks, xat = xat, xTickLabels = xTickLabels, col = Xcol, xlas = xlas, cex.axis = cex.axis, ...)
    }

    for (obj in c( ls() )) assign(obj, get(obj), envir = Graphics)
})
