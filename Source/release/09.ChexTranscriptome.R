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
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[["SampleID"]]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
sampleIDs <- SampleInfoFull[, "SampleID"]
bioGroups <- c(
    "K562", "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture",
    "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice"
)

EnsFeatures <- readRDS("Data/release/GenomicFeatures/EnsFeatures.RDS")
EnsFlanksByGene <- readRDS("Data/release/GenomicFeatures/EnsFlanksByGene.RDS")
EnsGeneIDsByGencodeLevel <- readRDS("Data/release/GenomicFeatures/EnsGeneIDsByGencodeLevel.RDS")
EnsFeatureLengthsByGene <- readRDS(file = "Data/release/GenomicFeatures/EnsFeatureLengthsByGene.RDS")
EnsFeatureLengthsByGene$none <- EnsFeatureLengthsByGene$mouse

K562RNAExonExprs <- readRDS("Data/release/Transcriptome/K562RNAExonExprs.RDS")
K562RNAIntronExprs <- readRDS("Data/release/Transcriptome/K562RNAIntronExprs.RDS")
K562RNAGeneExprs <- readRDS("Data/release/Transcriptome/K562RNAGeneExprs.RDS")

K562GROExonExprs <- readRDS("Data/release/Transcriptome/K562GROExonExprs.RDS")
K562GROIntronExprs <- readRDS("Data/release/Transcriptome/K562GROIntronExprs.RDS")
K562GROGeneExprs <- readRDS("Data/release/Transcriptome/K562GROGeneExprs.RDS")

K562scRNAExonExprs <- readRDS("Data/release/Transcriptome/K562scRNAExonExprs.RDS")
K562scRNAGeneExprs <- readRDS("Data/release/Transcriptome/K562scRNAGeneExprs.RDS")

HumanAstroRNAExonExprs <- readRDS("Data/release/Transcriptome/HumanAstroRNAExonExprs.RDS")
HumanAstroRNAIntronExprs <- readRDS("Data/release/Transcriptome/HumanAstroRNAIntronExprs.RDS")
HumanAstroRNAGeneExprs <- readRDS("Data/release/Transcriptome/HumanAstroRNAGeneExprs.RDS")

HumanNeuronRNAExonExprs <- readRDS("Data/release/Transcriptome/HumanNeuronRNAExonExprs.RDS")
HumanNeuronRNAIntronExprs <- readRDS("Data/release/Transcriptome/HumanNeuronRNAIntronExprs.RDS")
HumanNeuronRNAGeneExprs <- readRDS("Data/release/Transcriptome/HumanNeuronRNAGeneExprs.RDS")

MouseAstroRNAExonExprs <- readRDS("Data/release/Transcriptome/MouseAstroRNAExonExprs.RDS")
MouseAstroRNAIntronExprs <- readRDS("Data/release/Transcriptome/MouseAstroRNAIntronExprs.RDS")
MouseAstroRNAGeneExprs <- readRDS("Data/release/Transcriptome/MouseAstroRNAGeneExprs.RDS")

MouseNeuronRNAExonExprs <- readRDS("Data/release/Transcriptome/MouseNeuronRNAExonExprs.RDS")
MouseNeuronRNAIntronExprs <- readRDS("Data/release/Transcriptome/MouseNeuronRNAIntronExprs.RDS")
MouseNeuronRNAGeneExprs <- readRDS("Data/release/Transcriptome/MouseNeuronRNAGeneExprs.RDS")

bioGroupCols <- structure(if (length(bioGroups) <= 9) { brewer.pal(n = length(bioGroups), "Set1") } else { colorRampPalette(brewer.pal(n = 9, "Set1"))(length(bioGroups)) }, names = bioGroups)
trtGroups <- c("PositiveSingle", "PositivePooled", "PositiveMerged", "PositiveAll")
trtGroupCols <- structure(brewer.pal(n = length(trtGroups), "Greys"), names = trtGroups)
annotationColors <- list(BioGroup = bioGroupCols, TrtGroup = trtGroupCols)

Species <- c("human", "mouse")
bioGroupsBySpecies <- sapply(Species, function(species) bioGroups[bioGroups %in% unique(subset(SampleInfoFull, (IsNegCtrl == "N" & BioGroup %in% bioGroups) & ProbeType == "Positive" & Species == species)[["BioGroup"]])], simplify = FALSE)
bioGroupColsBySpecies <- sapply(Species, function(species) bioGroupCols[bioGroupsBySpecies[[species]]], simplify = FALSE)
trtGroupsBySpecies <- sapply(Species, function(species) trtGroups[trtGroups %in% unique(subset(SampleInfoFull, (IsNegCtrl == "N" & BioGroup %in% bioGroups) & ProbeType == "Positive" & Species == species)[["TrtGroup"]])], simplify = FALSE)
trtGroupColsBySpecies <- sapply(Species, function(species) trtGroupCols[trtGroupsBySpecies[[species]]], simplify = FALSE)

annotationColorsList <- sapply(Species, function(species) list(BioGroup = bioGroupColsBySpecies[[species]], TrtGroup = trtGroupColsBySpecies[[species]]), simplify = FALSE)

bioGroup2CellType <- sapply(bioGroups, function(bioGroup) subset(SampleInfoFull, (IsNegCtrl == "N" & BioGroup %in% bioGroups) & ProbeType == "Positive" & BioGroup == bioGroup)[1, "CellType"])
bioGroup2Species <- sapply(bioGroups, function(bioGroup) subset(SampleInfoFull, (IsNegCtrl == "N" & BioGroup %in% bioGroups) & ProbeType == "Positive" & BioGroup == bioGroup)[1, "Species"])
sampleIDsVirtmaxByBioGroup <- sapply(bioGroups, function(bioGroup) tail(subset(SampleInfoFull, IsNegCtrl == "N" & ProbeType == "Positive" & BioGroup == bioGroup), 1)[["SampleID"]])

EnsIDToSymbolMap <- sapply(Species, function(species) structure(mcols(EnsFeatures[["Gene"]][[species]])[["symbol"]], names = names(EnsFeatures[["Gene"]][[species]])), simplify = FALSE)

featureTypesByGene <- c("GeneExt", "Gene", "Flank5kByGene", "PromoterByGene", "FiveUTRByGene", "ExonByGene", "IntronByGene", "ThreeUTRByGene", "DownstreamByGene", "CDSByGene")
flankDistLabs <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")

## We finanlly decide on the quality class: ABreadCmate5End and ge20_le0.1_strict
qualOutInPairs <- c("ABreadCmate5End", "ABread5End", "Aread5End", "Bread5End", "Cmate5End", "Deither5End")[1]
mapqThs <- c("ge20_le0.1_strict", "ge30_le0.1", "ge20_le0.1", "ge10_le0.1")[1]

PFsFiltered <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        filename <- sprintf("Data/release/PrimingRateGene/PFs%sFiltered_%s.RDS", qualOutInPair, th)
        message(filename)
        readRDS(filename)
    }, simplify = FALSE)
}, simplify = FALSE)

FlankPFsFiltered <- sapply(mapqThs, function(th) {
    sapply(qualOutInPairs, function(qualOutInPair) {
        filename <- sprintf("Data/release/PrimingRateGene/FlankPFs%sFiltered_%s.RDS", qualOutInPair, th)
        message(filename)
        readRDS(filename)
    }, simplify = FALSE)
}, simplify = FALSE)

ChexExprsMap <- list(
    K562 = c("K562RNAGeneExprs", "K562RNAIntronExprs", "K562GROGeneExprs", "K562GROIntronExprs", "K562scRNAGeneExprs"),
    HumanAstroCulture = c("HumanAstroRNAGeneExprs", "HumanAstroRNAIntronExprs"), 
    HumanNeuronCulture = c("HumanNeuronRNAGeneExprs", "HumanNeuronRNAIntronExprs"), 
    HumanInterneuronCulture = c("HumanNeuronRNAGeneExprs", "HumanNeuronRNAIntronExprs"), 
    MouseAstroCulture = c("MouseAstroRNAGeneExprs", "MouseAstroRNAIntronExprs"), 
    MouseNeuronCulture = c("MouseNeuronRNAGeneExprs", "MouseNeuronRNAIntronExprs"), 
    MouseNeuronSlice = c("MouseNeuronRNAGeneExprs", "MouseNeuronRNAIntronExprs"), 
    MouseInterneuronSlice = c("MouseInterneuronRNAGeneExprs", "MouseInterneuronRNAIntronExprs")
)
TxMetaGroups <- unique(unname(unlist(ChexExprsMap)))
TxMetaGroupsBySpecies <- list(
    human = TxMetaGroups[grepl("^K562|^Human", TxMetaGroups)], 
    mouse = TxMetaGroups[grepl("^Mouse", TxMetaGroups)]
)

ChexExprsMapComplete <- sapply(bioGroups, function(bioGroup) {
    if (grepl("^K562|^Human", bioGroup)) { return(TxMetaGroupsBySpecies[["human"]]) } 
    else if (grepl("^Mouse", bioGroup)) { return(TxMetaGroupsBySpecies[["mouse"]]) }
}, simplify = FALSE)

BioGroupRNAFeatureTypeMap <- as.data.frame(tibble(
    BioGroup = rep(names(ChexExprsMap), lengths(ChexExprsMap)),
    Transcriptome = str_extract(unlist(ChexExprsMap), "RNA|GRO|scRNA"), 
    FeatureType = str_extract(unlist(ChexExprsMap), "Gene|Intron"), 
    DataType = str_extract(unlist(ChexExprsMap), "Exprs"), 
    MetaGroup = str_replace(unlist(ChexExprsMap), "(RNA|GRO|scRNA).*$", ""), 
))

BioGroupRNAFeatureTypeMapComplete <- as.data.frame(tibble(
    BioGroup = rep(names(ChexExprsMapComplete), lengths(ChexExprsMapComplete)),
    Transcriptome = str_extract(unlist(ChexExprsMapComplete), "RNA|GRO|scRNA"), 
    FeatureType = str_extract(unlist(ChexExprsMapComplete), "Gene|Intron"), 
    DataType = str_extract(unlist(ChexExprsMapComplete), "Exprs"), 
    MetaGroup = str_replace(unlist(ChexExprsMapComplete), "(RNA|GRO|scRNA).*$", ""), 
))

###########################################################################
## Venn diagram to show overlap between primed genes and expressed genes
## CHEX use 0 as cutoff; Trascriptome use zero or median
###########################################################################
for (qualOutInPair in qualOutInPairs) {
    for (mapqTh in mapqThs) {
        Ycutoffs <- c("median", "zero")
        gList <- list()
        dir.create(sprintf("Report/release/ChexTranscriptome/ChexExprsBinaryVenn/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
        for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
            bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
            species <- bioGroup2Species[bioGroup]
            Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
            TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
            DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
            MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
            TxLab <- paste0(MetaGroup, Tx, TxFeatureType, DataType)
            sampleID <- sampleIDsVirtmaxByBioGroup[bioGroup]
            SIDs <- strsplit(SampleInfoFull[sampleID, "SourceIDsNoOut"], ",")[[1]]
            gList[[bioGroup]][[TxLab]] <- NULL
            for (Ycutoff in Ycutoffs) {
                gList[[bioGroup]][[TxLab]][[Ycutoff]] <- NULL
                for (featureTypeByGene in featureTypesByGene) {
                    message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                    gList[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- NULL
                    Chex <- as.matrix(rowSums(PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SIDs, drop = FALSE]))
                    colnames(Chex) <- sampleID
                    g <- CHEX$vennChexTranscriptome(X = Chex, Y = get(TxLab), assayX = "CHEX", assayY = Tx, featureX = featureTypeByGene, featureY = TxFeatureType, Ycutoff = Ycutoff, filenameFormat = NULL, sampleID = sampleID, height = 2 * 300, width = 2 * 300, res = 300, hyper.test = TRUE, lower.tail = FALSE, cex = c(1, 1, 1) * 0.6, sub.cex = 1, main.cex = 1, cat.cex = 0, fill = c("#984EA3", "#FFFFE0"), alpha = 0.5, lwd = 0.4, fontfamily = "Helvetica")
                    gList[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- g
                }
            }
        }
        for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
            bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
            species <- bioGroup2Species[bioGroup]
            Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
            TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
            DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
            MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
            TxLab <- paste0(MetaGroup, Tx, TxFeatureType, DataType)
            grobTrees <- vector("list", length(Ycutoffs) * length(featureTypesByGene))
            i <- 0
            filename <- sprintf("Report/release/ChexTranscriptome/ChexExprsBinaryVenn/%s_%s/%s_%s.pdf", qualOutInPair, mapqTh, bioGroup, TxLab)
            pdf(filename, width = 15, height = 3)
            for (Ycutoff in Ycutoffs) {
                for (featureTypeByGene in featureTypesByGene) {
                    i <- i + 1
                    message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                    grobTrees[[i]] <- grobTree(gList[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]])
                    if (Ycutoff == "median") {
                        if (featureTypeByGene == "GeneExt") {
                            grobTrees[[i]] <- arrangeGrob(grobTrees[[i]], left = Ycutoff, top = featureTypeByGene)
                        } else {
                            grobTrees[[i]] <- arrangeGrob(grobTrees[[i]], top = featureTypeByGene)
                        }
                    } else {
                        if  (featureTypeByGene == "GeneExt") {
                            grobTrees[[i]] <- arrangeGrob(grobTrees[[i]], left = Ycutoff)
                        }
                    }
                }
            }
            do.call(grid.arrange, c(grobTrees, list(ncol = length(featureTypesByGene), widths = rep(1, length(featureTypesByGene)), heights = rep(1, 2))))
            dev.off()
        }
    }
}
###########################################################################
## Correlate CHEX priming density and gene expression levels
###########################################################################
dir.create("Report/release/ChexTranscriptome/ExprsByChexCellNum", FALSE, TRUE)
for (qualOutInPair in qualOutInPairs[1]) {
    for (mapqTh in mapqThs[1]) {
        for (bioGroup in bioGroups[11]) {
            Exprs <- ChexExprsMap[[bioGroup]]
            species <- bioGroup2Species[bioGroup]
            sampleIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & Composition == "Single" & IsOut == "N")[, "SampleID"]
            if (length(sampleIDs) > 1) {
                for (Expr in Exprs[2]) {
                    dir.create(sprintf("Report/release/ChexTranscriptome/ExprsByChexCellNum/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
                    filename <- sprintf("Report/release/ChexTranscriptome/ExprsByChexCellNum/%s_%s/%s_%s.pdf", qualOutInPair, mapqTh, bioGroup, Expr)
                    exprs <- get(Expr)
                    exprsInt <- rowMeans(exprs)
                    pdf(filename, height = 3, width = 3)
                    for (featureTypeByGene in featureTypesByGene) {
                        Chex <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sampleIDs, drop = FALSE]
                        ChexInt <- rowSums(Chex > 0)
                        genes <- intersect(names(ChexInt), names(exprsInt))
                        ## overlapgenes <- readLines("Report/release/ChexTranscriptome/FilterOverlapGenes/EnsGeneMouse_OverlappingGenes.txt")
                        ## genes <- setdiff(genes, overlapgenes)
                        message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", Expr, " ", featureTypeByGene, "...")
                        X <- data.frame(Chex = as.integer(ChexInt[genes]), Exprs = exprsInt[genes])
                        print(ggplot(subset(X, Chex > 0 & Exprs > 0), aes(x = Chex, group = Chex, y = Exprs)) + geom_boxplot(outlier.shape = NA) + scale_x_continuous() + scale_y_continuous(trans = "log10") + xlab("#. of cells primed") + ylab("Expression") + ggtitle(featureTypeByGene) + theme_classic(base_size = 6) + theme(panel.grid = element_blank(), panel.background = element_blank(), title = element_text(size = 6)) + geom_jitter(alpha = 1/20, size = 0.1))
                    }
                    dev.off()
                }
            }
        }
    }
}

###########################################################################
## Chex-seq number of primed cells vs Expression variability
###########################################################################
dir.create("Report/release/ChexTranscriptome/ExprsVarByChexCellNum", FALSE, TRUE)
for (qualOutInPair in qualOutInPairs) {
    for (mapqTh in mapqThs) {
        dir.create(sprintf("Report/release/ChexTranscriptome/ExprsVarByChexCellNum/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
        for (bioGroup in bioGroups[12]) {
            Exprs <- ChexExprsMap[[bioGroup]]
            species <- bioGroup2Species[bioGroup]
            sampleIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & Composition == "Single" & IsOut == "N")[, "SampleID"]
            for (Expr in Exprs) {
                exprs <- get(Expr)
                if (length(sampleIDs) > 1 & ncol(exprs) > 1) {
                    filename <- sprintf("Report/release/ChexTranscriptome/ExprsVarByChexCellNum/%s_%s/%s_%s.pdf", qualOutInPair, mapqTh, bioGroup, Expr)
                    exprsCV <- apply(exprs, 1, Stats$cv)
                    pdf(filename, height = 3, width = 3)
                    for (featureTypeByGene in featureTypesByGene) {
                        Chex <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sampleIDs, drop = FALSE]
                        Chex <- rowSums(Chex > 0)
                        genes <- intersect(names(Chex), names(exprsCV))
                        ## overlapgenes <- readLines("Report/release/ChexTranscriptome/FilterOverlapGenes/EnsGeneMouse_OverlappingGenes.txt")
                        ## genes <- setdiff(genes, overlapgenes)
                        message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", featureTypeByGene, " ", Expr, "...")
                        X <- data.frame(Chex = Chex[genes], Exprs = exprsCV[genes])
                        print(ggplot(subset(X, Chex > 0 & Exprs > 0), aes(x = Chex, group = Chex, y = Exprs)) + geom_boxplot(outlier.shape = NA) + scale_x_continuous() + scale_y_continuous(trans = "log1p") + xlab("#. of primed cells") + ylab("Expression variability (CV)") + ggtitle(featureTypeByGene) + theme_classic(base_size = 6) + theme(panel.background = element_blank(), panel.grid = element_blank(), title = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank()) + geom_jitter(alpha = 1/20, size = 0.1))
                    }
                    dev.off()
                }
            }
        }
    }
}

###########################################################################
## We hypothesize that genes with >0 CHEX counts close to TSS have higher experssion than those with CHEX counts at distance
###########################################################################
dir.create("Report/release/ChexTranscriptome/ExprsByChexTSSDist", FALSE, TRUE)
for (qualOutInPair in qualOutInPairs) {
    for (mapqTh in mapqThs) {
        dir.create(sprintf("Report/release/ChexTranscriptome/ExprsByChexTSSDist/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
        for (bioGroup in bioGroups) {
            Exprs <- ChexExprsMap[[bioGroup]]
            species <- bioGroup2Species[bioGroup]
            sampleID <- sampleIDsVirtmaxByBioGroup[bioGroup]
            SIDs <- strsplit(SampleInfoFull[sampleID, "SourceIDsNoOut"], ",")[[1]]
            for (Expr in Exprs) {
                filename <- sprintf("Report/release/ChexTranscriptome/ExprsByChexTSSDist/%s_%s/%s_%s.pdf", qualOutInPair, mapqTh, bioGroup, Expr)
                exprs <- get(Expr)
                avgExprs <- rowMeans(exprs)
                pdf(filename, height = 3, width = 3)
                exprsByDist <- sapply(flankDistLabs, function(flankDistLab) {
                    message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", sampleID, " ", flankDistLab, " ", Expr, "...")
                    Chex <- FlankPFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[flankDistLab]][, SIDs, drop = FALSE]
                    Chex <- as.matrix(rowSums(Chex))
                    g <- intersect(rownames(Chex), rownames(exprs))
                    g <- g[Chex[g, 1] > 0]
                    avgExprs[g]
                }, simplify = FALSE)
                X <- do.call(rbind, lapply(flankDistLabs, function(d) if (length(exprsByDist[[d]]) == 0) { data.frame(Dist = character(0), Expr = numeric(0)) } else { data.frame(Dist = d, Expr = exprsByDist[[d]]) }))
                X$Dist <- factor(sub("Flank", "", X$Dist), levels = sub("Flank", "", flankDistLabs))
                if (bioGroup == "K562" && Expr == "K562scRNAGeneExprs") {
                    print(ggplot(X, aes(x = Dist, y = Expr)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.02, alpha = 1/20) + xlab("Distance to TSS with CHEX priming") + ylab("Expression") + theme_classic(base_size = 12) + theme(panel.grid = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()) + ggtitle(Expr) + scale_y_continuous(trans = "log1p", limits = c(0, quantile(X$Expr, 0.9)))) ## 90%-tile capping to remove the outliers.
                } else {
                    print(ggplot(X, aes(x = Dist, y = Expr)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.02, alpha = 1/20) + xlab("Distance to TSS with CHEX priming") + ylab("Expression") + theme_classic(base_size = 12) + theme(panel.grid = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()) + ggtitle(Expr) + scale_y_continuous(trans = "log1p"))
                }
                dev.off()
            }
        }
    }
}

###########################################################################
## Test for significance in the overlap between genes with
## CHEX > 0 and genes with expression > median
###########################################################################
## sampleIDVirtmax in each CHEX bioGroup x each TranscriptomeComplete
## #. of bioGroups x #. of ChexExprsMapComplete
dir.create("Data/release/ChexTranscriptome", FALSE, TRUE)
for (qualOutInPair in qualOutInPairs[1]) {
    for (mapqTh in mapqThs[1]) {
        ChexExprsVirmaxCompleteBinTables <- sapply(bioGroups, function(bioGroup) {
            Exprs <- ChexExprsMapComplete[[bioGroup]]
            species <- bioGroup2Species[bioGroup]
            sampleID <- sampleIDsVirtmaxByBioGroup[bioGroup]
            SIDs <- strsplit(SampleInfoFull[sampleID, "SourceIDsNoOut"], ",")[[1]]
            sapply(Exprs, function(Expr) {
                exprs <- get(Expr)
                exprs <- rowMeans(exprs)
                exprs <- t(t(exprs))
                BinTables <- lapply(featureTypesByGene, function(featureTypeByGene) {
                    Chex <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SIDs, drop = FALSE]
                    Chex <- as.matrix(rowSums(Chex))
                    message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", sampleID, " ", Expr, " ", featureTypeByGene, "...")
                    CHEX$tabulateExprs(Chex, exprs, cutoff2 = "median")
                })
                names(BinTables) <- featureTypesByGene
                BinTables
            }, simplify = FALSE)
        }, simplify = FALSE)
        saveRDS(ChexExprsVirmaxCompleteBinTables, file = sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsVirmaxCompleteBinTables.RDS", qualOutInPair, mapqTh))
    }
}

ChexExprsVirmaxCompleteBinTables <- readRDS(sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsVirmaxCompleteBinTables.RDS", qualOutInPair, mapqTh))

delta <- 1
ChexExprsVirtmaxCompleteBinFETpvals <- rapply(ChexExprsVirmaxCompleteBinTables, function(x) { fisher.test(x + delta, alternative = "greater")$p.value}, how = "replace")
dir.create(sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
ChexExprsVirtmaxCompleteBinFETpvalsTab <- sapply(bioGroups, function(bioGroup) {
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    sapply(Exprs, function(Expr) {
        sapply(featureTypesByGene, function(featureTypeByGene) {
            ChexExprsVirtmaxCompleteBinFETpvals[[bioGroup]][[Expr]][[featureTypeByGene]]
        })
    }, simplify = FALSE)
}, simplify = FALSE)
ChexExprsVirtmaxCompleteBinFETpvalsByBioGroup <- sapply(bioGroups, function(bioGroup) {
    pvalsList <- ChexExprsVirtmaxCompleteBinFETpvalsTab[[bioGroup]]
    TxLabs <- names(pvalsList)
    pvalsDf <- do.call(cbind, pvalsList)
    colnames(pvalsDf) <- TxLabs
    pvalsDf
}, simplify = FALSE)
log10pBase <- 1e-200
ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup <- lapply(bioGroups, function(bioGroup) {
    X <- t(ChexExprsVirtmaxCompleteBinFETpvalsByBioGroup[[bioGroup]])
    rownames(X) <- sub("Exprs", "", rownames(X))
    colnames(X) <- sub("ByGene", "", colnames(X))
    pheatmap(-log10(X + log10pBase), cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 8, cellheight = 8, fontsize = 6, fontsize_col = 6, fontsize_row = 6, angle_col = 45, main = bioGroup, border_color = NA)[[4]]
})
ggsave(filename = sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup.pdf", qualOutInPair, mapqTh), plot = grid.arrange(grobs = ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup, ncol = 1), height = dim(BioGroupRNAFeatureTypeMapComplete) / 5, width = 5, unit = "in", limitsize = FALSE)
ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup_SameScale <- lapply(bioGroups, function(bioGroup) {
    X <- t(ChexExprsVirtmaxCompleteBinFETpvalsByBioGroup[[bioGroup]])
    rownames(X) <- sub("Exprs", "", rownames(X))
    colnames(X) <- sub("ByGene", "", colnames(X))
    pheatmap(-log10(X + log10pBase), cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 8, cellheight = 8, fontsize = 6, fontsize_col = 6, fontsize_row = 6, angle_col = 45, main = bioGroup, border_color = NA, col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(255), breaks = seq(0, 200, length.out = 256))[[4]]
})
ggsave(filename = sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup_SameScale.pdf", qualOutInPair, mapqTh), plot = grid.arrange(grobs = ChexExprsVirtmaxCompleteBinFETpvalHmsByBioGroup_SameScale, ncol = 1), height = dim(BioGroupRNAFeatureTypeMapComplete) / 5, width = 5, unit = "in", limitsize = FALSE)
 
ChexExprsVirtmaxCompleteBinFEToddRatios <- rapply(ChexExprsVirmaxCompleteBinTables, function(x) { fisher.test(x + delta, alternative = "greater")$estimate}, how = "replace")
ChexExprsVirtmaxCompleteBinFEToddRatiosTab <- sapply(bioGroups, function(bioGroup) {
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    sapply(Exprs, function(Expr) {
        oddRatioTab <- sapply(featureTypesByGene, function(featureTypeByGene) {
            unname(ChexExprsVirtmaxCompleteBinFEToddRatios[[bioGroup]][[Expr]][[featureTypeByGene]])
        })
    }, simplify = FALSE)
}, simplify = FALSE)
ChexExprsVirtmaxCompleteBinFEToddRatiosByBioGroup <- sapply(bioGroups, function(bioGroup) {
    oddRatiosList <- ChexExprsVirtmaxCompleteBinFEToddRatiosTab[[bioGroup]]
    TxLabs <- names(oddRatiosList)
    oddRatiosDf <- do.call(cbind, oddRatiosList)
    colnames(oddRatiosDf) <- TxLabs
    oddRatiosDf
}, simplify = FALSE)
ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup <- lapply(bioGroups, function(bioGroup) {
    X <- t(ChexExprsVirtmaxCompleteBinFEToddRatiosByBioGroup[[bioGroup]])
    rownames(X) <- sub("Exprs", "", rownames(X))
    colnames(X) <- sub("ByGene", "", colnames(X))
    pheatmap(log2(X), cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 8, cellheight = 8, fontsize = 6, fontsize_col = 6, fontsize_row = 6, angle_col = 45, border_color = NA, main = bioGroup)[[4]]
})
ggsave(filename = sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup.pdf", qualOutInPair, mapqTh), plot = grid.arrange(grobs = ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup, ncol = 1), height = dim(BioGroupRNAFeatureTypeMapComplete) / 5, width = 5, unit = "in", limitsize = FALSE)
ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup_SameScale <- lapply(bioGroups, function(bioGroup) {
    X <- t(ChexExprsVirtmaxCompleteBinFEToddRatiosByBioGroup[[bioGroup]])
    rownames(X) <- sub("Exprs", "", rownames(X))
    colnames(X) <- sub("ByGene", "", colnames(X))
    pheatmap(log2(X), cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 8, cellheight = 8, fontsize = 6, fontsize_col = 6, fontsize_row = 6, angle_col = 45, border_color = NA, main = bioGroup, col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(255), breaks = seq(-2, 4, length.out = 256))[[4]]
})
ggsave(filename = sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup_SameScale.pdf", qualOutInPair, mapqTh), plot = grid.arrange(grobs = ChexExprsVirtmaxCompleteBinFEToddRatioHmsByBioGroup_SameScale, ncol = 1), height = dim(BioGroupRNAFeatureTypeMapComplete) / 5, width = 5, unit = "in", limitsize = FALSE)

###########################################################################
## Each individual CHEX SampleFull x each TranscriptomeComplete
###########################################################################
dir.create("Data/release/ChexTranscriptome", FALSE, TRUE)
ncores <- 10
for (qualOutInPair in qualOutInPairs[1]) {
    for (mapqTh in mapqThs[1]) {
        for (bioGroup in bioGroups) {
            Exprs <- ChexExprsMapComplete[[bioGroup]]
            species <- bioGroup2Species[bioGroup]
            SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
            BinTables <- mclapply(SIDs, function(SID) {
                sapply(Exprs, function(Expr) {
                    exprs <- get(Expr)
                    exprs <- rowMeans(exprs)
                    exprs <- t(t(exprs))
                    sapply(featureTypesByGene, function(featureTypeByGene) {
                        if (SampleInfoFull[SID, "CompType"] == "Virtual") {
                            sids <- strsplit(SampleInfoFull[SID, "SourceIDsNoOut"], ",")[[1]]
                            Chex <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sids, drop = FALSE]
                            Chex <- as.matrix(rowSums(Chex))
                        } else {
                            Chex <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SID, drop = FALSE]
                        }
                        message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", SID, " ", Expr, " ", featureTypeByGene, "...")
                        CHEX$tabulateExprs(Chex, exprs, cutoff2 = "median")
                    }, simplify = FALSE)
                }, simplify = FALSE)
            }, mc.cores = ncores)
            names(BinTables) <- SIDs
            filename <- sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsCompleteBinTables_%s.RDS", qualOutInPair, mapqTh, bioGroup)
            saveRDS(BinTables, file = filename)
        }
    }
}
ChexExprsCompleteBinTables <- sapply(bioGroups, function(bioGroup) {
    filename <- sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsCompleteBinTables_%s.RDS", qualOutInPair, mapqTh, bioGroup)
    readRDS(file = filename)
}, simplify = FALSE)

###########################################################################
## Instead of heatmaps by BioGroup, we plot by each sample.
## setup palette for pheatmap annotation
## all transcriptomes, per cell association test
###########################################################################
delta <- 1
ChexExprsCompleteBinFETpvals <- rapply(ChexExprsCompleteBinTables, function(x) { fisher.test(x + delta, alternative = "greater")$p.value}, how = "replace")
ChexExprsCompleteBinFETpvalsTab <- sapply(bioGroups, function(bioGroup) {
    SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    sapply(Exprs, function(Expr) {
        sapply(SIDs, function(SID) {
            sapply(featureTypesByGene, function(featureTypeByGene) {
                message(bioGroup, " ", SID, " ", Expr, " ", featureTypeByGene, "...")
                ChexExprsCompleteBinFETpvals[[bioGroup]][[SID]][[Expr]][[featureTypeByGene]]
            })
        })
    }, simplify = FALSE)
}, simplify = FALSE)

ChexExprsCompleteBinFEToddsRatio <- rapply(ChexExprsCompleteBinTables, function(x) unname(fisher.test(x + delta, alternative = "greater")$estimate), how = "replace")
ChexExprsCompleteBinFEToddsRatioTab <- sapply(bioGroups, function(bioGroup) {
    SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    sapply(Exprs, function(Expr) {
        sapply(SIDs, function(SID) {
            sapply(featureTypesByGene, function(featureTypeByGene) {
                message(bioGroup, " ", SID, " ", Expr, " ", featureTypeByGene, "...")
                ChexExprsCompleteBinFEToddsRatio[[bioGroup]][[SID]][[Expr]][[featureTypeByGene]][[1]]
            })
        })
    }, simplify = FALSE)
}, simplify = FALSE)

dirname <- sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsBinFEToddsRatioHms", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)
for (bioGroup in bioGroups) {
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    species <- bioGroup2Species[bioGroup]
    SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
    hmList <- list()
    for (Expr in Exprs) {
        X <- ChexExprsCompleteBinFEToddsRatioTab[[bioGroup]][[Expr]]
        if (ncol(X) > 1) {
            filename <- sprintf("%s/%s_%s.pdf", dirname, bioGroup, Expr)
            pdf(filename, height = 6, width = 30)
            pheatmap(log2(X), cluster_rows = FALSE, cluster_cols = TRUE, annotation_col = SampleInfoFull[SIDs, c("BioGroup", "TrtGroup")], scale = "none", fontsize_row = 10, fontsize_col = 10, show_colnames = TRUE, annotation_colors = annotationColorsList[[species]], labels_col = with(SampleInfoFull[SIDs, ], paste0("E.", ExptID, " #", sub("scCLTdegenNuc", "", SampleID))), border_color = NA, angle_col = 45, cellwidth = 11, cellheight = 11, main = Expr)
            dev.off()
        }
    }
}

dirname <- sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsBinFEToddsRatioHms_NoClust", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)
for (bioGroup in bioGroups) {
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    species <- bioGroup2Species[bioGroup]
    SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
    hmList <- list()
    for (Expr in Exprs) {
        X <- ChexExprsCompleteBinFEToddsRatioTab[[bioGroup]][[Expr]]
        if (ncol(X) > 1) {
            filename <- sprintf("%s/%s_%s.pdf", dirname, bioGroup, Expr)
            pdf(filename, height = 6, width = 30)
            Y <- log2(X[c(1, 4, 6, 7), ])
            Y1 <- Y[, SampleInfoFull[colnames(Y), "CompType"] == "Biol"]
            Y2 <- Y[, SampleInfoFull[colnames(Y), "CompType"] == "Virtual"]
            Y1 <- Y1[, hclust(dist(t(Y1)))[["order"]]]
            Z <- cbind(Y1, Y2)
            pheatmap(Z, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = SampleInfoFull[SIDs, c("TrtGroup"), drop = FALSE], scale = "none", fontsize_row = 10, fontsize_col = 10, show_colnames = TRUE, annotation_colors = annotationColorsList[[species]], labels_col = with(SampleInfoFull[SIDs, ], paste0("E.", ExptID, " #", sub("scCLTdegenNuc", "", SampleID))), border_color = NA, angle_col = 45, cellwidth = 11, cellheight = 11, main = Expr, treeheight_col = 0, treeheight_row = 0)
            dev.off()
        }
    }
}
 
dirname <- sprintf("Report/release/ChexTranscriptome/ChexExprsCompleteBinTable/%s_%s/ChexExprsBinFETpvalHms_NoClust", qualOutInPair, mapqTh)
dir.create(dirname, FALSE, TRUE)
for (bioGroup in bioGroups) {
    Exprs <- ChexExprsMapComplete[[bioGroup]]
    species <- bioGroup2Species[bioGroup]
    SIDs <- subset(SampleInfoFull, IsNegCtrl == "N" & BioGroup == bioGroup & ProbeType == "Positive" & IsOut == "N")[, "SampleID"]
    hmList <- list()
    for (Expr in Exprs) {
        X <- ChexExprsCompleteBinFETpvalsTab[[bioGroup]][[Expr]]
        if (ncol(X) > 1) {
            filename <- sprintf("%s/%s_%s.pdf", dirname, bioGroup, Expr)
            pdf(filename, height = 6, width = 30)
            pheatmap(-log10(X + log10pBase), cluster_rows = FALSE, cluster_cols = TRUE, annotation_col = SampleInfoFull[SIDs, c("BioGroup", "TrtGroup")], scale = "none", fontsize_row = 10, fontsize_col = 10, show_colnames = TRUE, annotation_colors = annotationColorsList[[species]], labels_col = with(SampleInfoFull[SIDs, ], paste0("E.", ExptID, " #", sub("scCLTdegenNuc", "", SampleID))), border_color = NA, angle_col = 45, cellwidth = 11, cellheight = 11, main = Expr)
            dev.off()
        }
    }
}

###########################################################################
## Given the CHEX-RNA overlap, we want to know the functions enriched in 
## 1. genes shared by CHEX and RNA
## 2. genes unique to CHEX
## 3. genes unique to RNA
###########################################################################
## GO enrichment analysis of overlapping genes
OrgDbGOSymbolList <- readRDS(file = "Data/release/GenomicFeatures/Annotations/OrgDbGOSymbolList.RDS")
OrgDbGOSymbolList$none <- OrgDbGOSymbolList$mouse
OrgDbGONameList <- readRDS(file = "Data/release/GenomicFeatures/Annotations/OrgDbGONameList.RDS")
OrgDbGONameList$none <- OrgDbGONameList$mouse
GOonts <- c("MF", "BP", "CC")
Ycutoffs <- c("median")

for (qualOutInPair in qualOutInPairs) {
    for (mapqTh in mapqThs) {
        CommDiffGenes <- NULL
        for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
            bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
            species <- bioGroup2Species[bioGroup]
            Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
            TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
            DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
            TxLab <- paste0(Tx, TxFeatureType, DataType)
            MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
            SID <- sampleIDsVirtmaxByBioGroup[bioGroup]
            Y <- get(paste0(MetaGroup, Tx, TxFeatureType, DataType))
            for (Ycutoff in Ycutoffs) {
                for (featureTypeByGene in featureTypesByGene) {
                    message(qualOutInPair, " ", mapqTh, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                    CommDiffGenes[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- NULL
                    if (SampleInfoFull[SID, "CompType"] == "Biol") {
                        X <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SID, drop = FALSE]
                    } else {
                        sids <- strsplit(SampleInfoFull[SID, "SourceIDsNoOut"], ",")[[1]]
                        X <- as.matrix(rowSums(PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sids, drop = FALSE]))
                        colnames(X) <- SID
                    }
                    genesList <- CHEX$ChexTxCommDiff(X = X, Y = Y, Ycutoff = Ycutoff)
                    names(genesList) <- c("Comm", "ChexUniq", "TxUniq")
                    CommDiffGenes[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- genesList
                }
            }
        }
        saveRDS(CommDiffGenes, file = sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsCommDiffGenes.RDS", qualOutInPair, mapqTh))
    }
}
CommDiffGenes <- readRDS(file = sprintf("Data/release/ChexTranscriptome/%s_%s_ChexExprsCommDiffGenes.RDS", qualOutInPair, mapqTh))
 
for (GOont in GOonts) {
    filename <- sprintf("Data/release/ChexTranscriptome/%s_%s_GOenrichChexTxOverlap_%s.RDS", qualOutInPair, mapqTh, GOont)
    GOenrichChexTxOverlap <- list()
    for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
        bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
        species <- bioGroup2Species[bioGroup]
        GOSymbol <- OrgDbGOSymbolList[[species]][[GOont]]
        GOName <- OrgDbGONameList[[species]][[GOont]]
        Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
        TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
        DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
        TxLab <- paste0(Tx, TxFeatureType, DataType)
        MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
        SID <- sampleIDsVirtmaxByBioGroup[bioGroup]
        Y <- get(paste0(MetaGroup, Tx, TxFeatureType, DataType))
        for (Ycutoff in Ycutoffs) {
            for (featureTypeByGene in featureTypesByGene[1]) {
                message(GOont, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                if (SampleInfoFull[SID, "CompType"] == "Biol") {
                    X <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SID, drop = FALSE]
                } else {
                    sids <- strsplit(SampleInfoFull[SID, "SourceIDsNoOut"], ",")[[1]]
                    X <- as.matrix(rowSums(PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sids, drop = FALSE]))
                    colnames(X) <- SID
                }
                GOenrichChexTxOverlap[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- NULL
                GOenrichChexTxOverlap[[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- CHEX$enrichChexTxOverlap(X = X, Y = Y, Ycutoff = Ycutoff, TERM2GENE = GOSymbol, TERM2NAME = GOName)
            }
        }
    }
   saveRDS(GOenrichChexTxOverlap, file = filename)
}

for (GOont in GOonts) {
    filename <- sprintf("Data/release/ChexTranscriptome/%s_%s_GOenrichChexUniq_%s.RDS", qualOutInPair, mapqTh, GOont)
    GOenrichChexUniq <- list()
    for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
        bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
        species <- bioGroup2Species[bioGroup]
        GOSymbol <- OrgDbGOSymbolList[[species]][[GOont]]
        GOName <- OrgDbGONameList[[species]][[GOont]]
        Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
        TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
        DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
        TxLab <- paste0(Tx, TxFeatureType, DataType)
        MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
        SID <- sampleIDsVirtmaxByBioGroup[bioGroup]
        Y <- get(paste0(MetaGroup, Tx, TxFeatureType, DataType))
        for (Ycutoff in Ycutoffs) {
            for (featureTypeByGene in featureTypesByGene) {
                message(GOont, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                if (SampleInfoFull[SID, "CompType"] == "Biol") {
                    X <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SID, drop = FALSE]
                } else {
                    sids <- strsplit(SampleInfoFull[SID, "SourceIDs"], ",")[[1]]
                    X <- as.matrix(rowSums(PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sids, drop = FALSE]))
                    colnames(X) <- SID
                }
                GOenrichChexUniq[[GOont]][[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- NULL
                GOenrichChexUniq[[GOont]][[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- CHEX$enrichChexUniq(X = X, Y = Y, Ycutoff = Ycutoff, TERM2GENE = GOSymbol, TERM2NAME = GOName)
            }
        }
    }
    saveRDS(GOenrichChexUniq, file = filename)
}

for (GOont in GOonts) {
    filename <- sprintf("Data/release/ChexTranscriptome/%s_%s_GOenrichTxUniq_%s.RDS", qualOutInPair, mapqTh, GOont)
    GOenrichTxUniq <- list()
    for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
        bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
        species <- bioGroup2Species[bioGroup]
        GOSymbol <- OrgDbGOSymbolList[[species]][[GOont]]
        GOName <- OrgDbGONameList[[species]][[GOont]]
        Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
        TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
        DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
        TxLab <- paste0(Tx, TxFeatureType, DataType)
        MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
        SID <- sampleIDsVirtmaxByBioGroup[bioGroup]
        Y <- get(paste0(MetaGroup, Tx, TxFeatureType, DataType))
        for (Ycutoff in Ycutoffs) {
            for (featureTypeByGene in featureTypesByGene) {
                message(GOont, " ", bioGroup, " ", TxLab, " ", Ycutoff, " ", featureTypeByGene)
                if (SampleInfoFull[SID, "CompType"] == "Biol") {
                    X <- PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, SID, drop = FALSE]
                } else {
                    sids <- strsplit(SampleInfoFull[SID, "SourceIDs"], ",")[[1]]
                    X <- as.matrix(rowSums(PFsFiltered[[mapqTh]][[qualOutInPair]][[species]][[featureTypeByGene]][, sids, drop = FALSE]))
                    colnames(X) <- SID
                }
                GOenrichTxUniq[[GOont]][[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- NULL
                GOenrichTxUniq[[GOont]][[bioGroup]][[TxLab]][[Ycutoff]][[featureTypeByGene]] <- CHEX$enrichTxUniq(X = X, Y = Y, Ycutoff = Ycutoff, TERM2GENE = GOSymbol, TERM2NAME = GOName)
            }
        }
    }
    saveRDS(GOenrichTxUniq, file = filename)
}

###########################################################################
## Figures for GO enrichment of CHEX-RNA overlapping genes
###########################################################################
dir.create(sprintf("Report/release/ChexTranscriptome/EnrichChexExprsOverlap/%s_%s", qualOutInPair, mapqTh), FALSE, TRUE)
GOontLabs <- c(MF = "GO: Molecular Function", BP = "GO: Biological Process", CC = "GO: Cellular Components")
for (GOont in GOonts) {
    GOenrichChexTxOverlap <- readRDS(sprintf("Data/release/ChexTranscriptome/%s_%s_GOenrichChexTxOverlap_%s.RDS", qualOutInPair, mapqTh, GOont))
    for (i in seq(nrow(BioGroupRNAFeatureTypeMap))) {
        bioGroup <- BioGroupRNAFeatureTypeMap[i, "BioGroup"]
        Tx <- BioGroupRNAFeatureTypeMap[i, "Transcriptome"]
        TxFeatureType <- BioGroupRNAFeatureTypeMap[i, "FeatureType"]
        DataType <- BioGroupRNAFeatureTypeMap[i, "DataType"]
        TxLab <- paste0(Tx, TxFeatureType, DataType)
        MetaGroup <- BioGroupRNAFeatureTypeMap[i, "MetaGroup"]
        MetaTxLab <- paste0(MetaGroup, TxLab)
        pdf(sprintf("Report/release/ChexTranscriptome/EnrichChexExprsOverlap/%s_%s/GO_%s_%s_%s_median_GeneExt.pdf", qualOutInPair, mapqTh, GOont, bioGroup, MetaTxLab), width = 4, height = 2)
        par(mar = c(2, 10, 1, 2))
        top <- 20
        cols = colorRampPalette(c("lightblue", "steelblue"))(256)
        with(as.data.frame(GOenrichChexTxOverlap[[bioGroup]][[TxLab]]$median$GeneExt)[1:top, ], barplot(rev(-log10(p.adjust)), horiz = TRUE, names.arg = rev(Description), las = 1, cex.main = 0.8, cex.names = 0.4, cex.axis = 0.5, cex.sub = 0.6, col = cols[(rev(-log10(p.adjust)) - min(rev(-log10(p.adjust)))) / max(rev(-log10(p.adjust)) - min(rev(-log10(p.adjust)))) * 255 + 1], border = NA, xlab = "", main = GOontLabs[GOont]))
        text(x = max(-log10(GOenrichChexTxOverlap[[bioGroup]][[TxLab]]$median$GeneExt[1:top, "p.adjust"])), y = 1, label = "-log10(BH p-value)", xpd = TRUE, cex = 0.5)
        dev.off()
    }
}
