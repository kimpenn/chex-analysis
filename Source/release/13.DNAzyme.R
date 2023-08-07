## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/release/functions.R")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("data.table")
library("openxlsx")
library("rtracklayer")

webpage <- read.csv("Data/DNAmoreDB/webpage.csv", as.is = TRUE, check.names = FALSE)
dnazymes <- openxlsx::read.xlsx("Data/DNAmoreDB/dnazymes.xlsx")
dim(dnazymes)
## [1] 1782   26
dnazymes <- dnazymes[!duplicated(dnazymes[, 2]), ]
dim(dnazymes)
## [1] 1649   26

dnazymes[, 2] <- gsub("U", "T", dnazymes[, 2])

len <- nchar(dnazymes[, 2])
par(mar = c(5, 25, 1, 1))
barplot(sort(table(dnazymes[, 10])), horiz = TRUE, las = 1)

dnazymes <- dnazymes[len > 0, ]
len <- len[len > 0]
dim(dnazymes)
## [1] 1648   26


short <- dnazymes[len <= 10, ]
long <- dnazymes[len > 10, ]
short$id <- paste(rownames(short), short$name, sep = "_")
long$id <- paste(rownames(long), long$name, sep = "_")

Tools$write_xlsx(list(long = long, short = short), file = "Report/release/DNAzyme/dnazymes.xlsx", row.names = FALSE)

dnazyme_short <- short[, c("id", "e")]
dnazyme_long <- long[, c("id", "e")]

dir.create("Data/release/DNAzyme", FALSE, TRUE)
filename <- "Data/release/DNAzyme/dnazyme_short.fa"
cat("", sep = "", file = filename)
for (i in 1:nrow(dnazyme_short)) {
    cat(">", dnazyme_short[i, "id"], "\n", sep = "", file = filename, append = TRUE)
    cat(dnazyme_short[i, "e"], "\n", sep = "", file = filename, append = TRUE)
}

dir.create("Data/release/DNAzyme", FALSE, TRUE)
filename <- "Data/release/DNAzyme/dnazyme_long.fa"
cat("", sep = "", file = filename)
for (i in 1:nrow(dnazyme_long)) {
    cat(">", dnazyme_long[i, "id"], "\n", sep = "", file = filename, append = TRUE)
    cat(dnazyme_long[i, "e"], "\n", sep = "", file = filename, append = TRUE)
}

EnsFeatures <- readRDS("Data/release/GenomicFeatures/EnsFeatures.RDS")
EnsFeatures_Transcript_human <- EnsFeatures$Transcript$human
EnsFeatures_Transcript_mouse <- EnsFeatures$Transcript$mouse

EnsFeatures_Transcript_human <- EnsFeatures_Transcript_human[seqnames(EnsFeatures_Transcript_human) %in% Genome$MainSeqLevels$human]
EnsFeatures_Transcript_mouse <- EnsFeatures_Transcript_mouse[seqnames(EnsFeatures_Transcript_mouse) %in% Genome$MainSeqLevels$mouse]
seqlevels(EnsFeatures_Transcript_human) <- Genome$MainSeqLevels$human
seqlevels(EnsFeatures_Transcript_mouse) <- Genome$MainSeqLevels$mouse
seqinfo(EnsFeatures_Transcript_human) <- Genome$MainSeqInfo$human
seqinfo(EnsFeatures_Transcript_mouse) <- Genome$MainSeqInfo$mouse

blast_long_hg38 <- read.csv("Data/release/DNAzyme/blast_long_hg38.sam", head = FALSE, sep = "\t", as.is = TRUE)
colnames(blast_long_hg38) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", "score", "evalue", "mismatches", "matches", "gaps")
blast_long_hg38 <- data.frame(blast_long_hg38, long[sapply(strsplit(blast_long_hg38$qname, "_"), '[', 1), c("length", "reaction", "main_article_title", "main_article_pub_date")], row.names = NULL)
blast_long_hg38 <- blast_long_hg38[, -match(c("qual", "rnext", "pnext", "tlen"), colnames(blast_long_hg38))]
blast_long_hg38[, c("score", "evalue", "mismatches", "matches", "gaps")] <- sapply(c("score", "evalue", "mismatches", "matches", "gaps"), function(x) as.numeric(sapply(strsplit(blast_long_hg38[, x], ":"), "[", 3)))
blast_long_hg38 <- subset(blast_long_hg38, !grepl("_", blast_long_hg38$rname))
dim(blast_long_hg38)
## [1] 5638   16
grs_blast_long_hg38 <- GRanges(blast_long_hg38$rname, IRanges(blast_long_hg38$pos, blast_long_hg38$pos), strand = "*", seqinfo = Genome$MainSeqInfo$human)
ann_blast_long_hg38 <- annotatePeak(peak = grs_blast_long_hg38, TxDb = EnsFeatures_Transcript_human, tssRegion=c(-2000, 2000))
ann_blast_long_hg38_df <- as.data.frame(ann_blast_long_hg38)
all(as.character(grs_blast_long_hg38) == with(ann_blast_long_hg38_df, paste(seqnames, start, sep = ":")))
## [1] TRUE
blast_long_hg38 <- cbind(blast_long_hg38, ann_blast_long_hg38_df[, -match(c("seqnames", "start", "end", "width", "strand", "geneChr", "tx_cds_seq_start", "tx_cds_seq_end", "tx_name"), colnames(ann_blast_long_hg38_df))])
write.csv(blast_long_hg38, file = "Report/release/DNAzyme/blast_long_hg38.csv", row.names = FALSE)
blast_long_hg38 <- data.table::data.table(openxlsx::read.xlsx("Report/release/DNAzyme/blast_long.xlsx", sheet = "hg38"))

blast_long_mm10 <- read.csv("Data/release/DNAzyme/blast_long_mm10.sam", head = FALSE, sep = "\t", as.is = TRUE)
colnames(blast_long_mm10) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", "score", "evalue", "mismatches", "matches", "gaps")
blast_long_mm10 <- data.frame(blast_long_mm10, long[sapply(strsplit(blast_long_mm10$qname, "_"), '[', 1), c("length", "reaction", "main_article_title", "main_article_pub_date")], row.names = NULL)
blast_long_mm10 <- blast_long_mm10[, -match(c("qual", "rnext", "pnext", "tlen"), colnames(blast_long_mm10))]
blast_long_mm10[, c("score", "evalue", "mismatches", "matches", "gaps")] <- sapply(c("score", "evalue", "mismatches", "matches", "gaps"), function(x) as.numeric(sapply(strsplit(blast_long_mm10[, x], ":"), "[", 3)))
blast_long_mm10 <- subset(blast_long_mm10, !grepl("_", blast_long_mm10$rname))
dim(blast_long_mm10)
## [1] 6197   16
grs_blast_long_mm10 <- GRanges(blast_long_mm10$rname, IRanges(blast_long_mm10$pos, blast_long_mm10$pos), strand = "*", seqinfo = Genome$MainSeqInfo$mouse)
ann_blast_long_mm10 <- annotatePeak(peak = grs_blast_long_mm10, TxDb = EnsFeatures_Transcript_mouse, tssRegion=c(-2000, 2000))
ann_blast_long_mm10_df <- as.data.frame(ann_blast_long_mm10)
all(as.character(grs_blast_long_mm10) == with(ann_blast_long_mm10_df, paste(seqnames, start, sep = ":")))
## [1] TRUE
blast_long_mm10 <- cbind(blast_long_mm10, ann_blast_long_mm10_df[, -match(c("seqnames", "start", "end", "width", "strand", "geneChr", "tx_cds_seq_start", "tx_cds_seq_end", "tx_name"), colnames(ann_blast_long_mm10_df))])
write.csv(blast_long_mm10, file = "Report/release/DNAzyme/blast_long_mm10.csv", row.names = FALSE)
blast_long_mm10 <- data.table::data.table(openxlsx::read.xlsx("Report/release/DNAzyme/blast_long.xlsx", sheet = "mm10"))

pdf("Report/release/DNAzyme/blast_long_map_quality.pdf")
par(ps = 16)
with(blast_long_hg38, smoothScatter((matches + mismatches)/length*100, mismatches/(matches + mismatches)*100, xlab = "% of the catalytic core being mapped to genome", ylab = "% of mapped based being mismatches", main = "hg38"))
with(blast_long_mm10, smoothScatter((matches + mismatches)/length*100, mismatches/(matches + mismatches)*100, xlab = "% of the catalytic core being mapped to genome", ylab = "% of mapped based being mismatches", main = "mm10"))
dev.off()

###########################################################################
#' process the short group
###########################################################################
homer_short_hg38 <- data.table::fread("Data/release/DNAzyme/homer_short_hg38.bed.gz", header = FALSE, sep = "\t")
homer_short_hg38 <- homer_short_hg38[V1 %in% Genome$MainSeqLevels$human]
homer_short_hg38[V4 %like% '_\\(', V4 := gsub("C22_\\(", "C22 (", V4)]
grs_homer_short_hg38 <- homer_short_hg38[, list(.(GRanges(V1, IRanges(V2, V3), strand = V6, seqinfo = Genome$MainSeqInfo$human))), by = "V4"]
x <- grs_homer_short_hg38[, V4]
y <- grs_homer_short_hg38[, V1]
names(y) <- x
y <- as(y, "GRangesList")
grs_homer_short_hg38 <- y

## reshape(data = homer_short_hg38[, .N, .(V4, V1)], timevar = "V1", idvar = "V4", v.names = "N", times = Genome$MainSeqLevels$human, direction = "wide", sep = "_")
homer_short_hg38_cnts_bychr <- as.data.frame.matrix(xtabs(N ~ V4 + V1, homer_short_hg38[, .N, .(V4, V1)])[, Genome$MainSeqLevels$human])
homer_short_hg38_cnts_bychr <- data.frame(id = rownames(homer_short_hg38_cnts_bychr), homer_short_hg38_cnts_bychr, stringsAsFactors = FALSE)
homer_short_hg38_cnts_bychr <- merge(short, homer_short_hg38_cnts_bychr, by.x = "id", by.y = "id", all.x = TRUE)

grs_homer_short_hg38_byid <- by(homer_short_hg38, homer_short_hg38[, V4], function(X) with(X, GRanges(V1, IRanges(V2, V3), strand = V6, seqinfo = Genome$MainSeqInfo$human)))
grs_homer_short_hg38_byid <- as(lapply(grs_homer_short_hg38_byid, function(X) X), "GRangesList")

ann_homer_short_hg38_byid <- lapply(grs_homer_short_hg38_byid, function(X) { annotatePeak(peak = X, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, tssRegion = c(-2000, 2000)); })
ann_homer_short_hg38_byid <- lapply(ann_homer_short_hg38_byid, function(X) { X@annoStat })
ann_homer_short_hg38_byid <- t(do.call(cbind, lapply(ann_homer_short_hg38_byid, function(x) structure(x[["Frequency"]], names = as.character(x[["Feature"]])))))
ann_homer_short_hg38_byid <- data.frame(id = rownames(ann_homer_short_hg38_byid), ann_homer_short_hg38_byid, stringsAsFactors = FALSE, check.names = FALSE)

homer_short_hg38_cnts_bychr <- merge(homer_short_hg38_cnts_bychr, ann_homer_short_hg38_byid, by.x = "id", by.y = "id", all.x = TRUE)
write.csv(homer_short_hg38_cnts_bychr, file = "Report/release/DNAzyme/homer_short_hg38.csv", row.names = FALSE)

homer_short_mm10 <- data.table::fread("Data/release/DNAzyme/homer_short_mm10.bed.gz", header = FALSE, sep = "\t")
homer_short_mm10 <- homer_short_mm10[V1 %in% Genome$MainSeqLevels$mouse]
homer_short_mm10[V4 %like% '_\\(', V4 := gsub("C22_\\(", "C22 (", V4)]
grs_homer_short_mm10 <- homer_short_mm10[, list(.(GRanges(V1, IRanges(V2, V3), strand = V6, seqinfo = Genome$MainSeqInfo$mouse))), by = "V4"]
x <- grs_homer_short_mm10[, V4]
y <- grs_homer_short_mm10[, V1]
names(y) <- x
y <- as(y, "GRangesList")
grs_homer_short_mm10 <- y


homer_short_mm10_cnts_bychr <- as.data.frame.matrix(xtabs(N ~ V4 + V1, homer_short_mm10[, .N, .(V4, V1)])[, Genome$MainSeqLevels$mouse])
homer_short_mm10_cnts_bychr <- data.frame(id = rownames(homer_short_mm10_cnts_bychr), homer_short_mm10_cnts_bychr, stringsAsFactors = FALSE)
homer_short_mm10_cnts_bychr <- merge(short, homer_short_mm10_cnts_bychr, by.x = "id", by.y = "id", all.x = TRUE)

grs_homer_short_mm10_byid <- by(homer_short_mm10, homer_short_mm10[, V4], function(X) with(X, GRanges(V1, IRanges(V2, V3), strand = V6, seqinfo = Genome$MainSeqInfo$mouse)))
grs_homer_short_mm10_byid <- as(lapply(grs_homer_short_mm10_byid, function(X) X), "GRangesList")

ann_homer_short_mm10_byid <- lapply(grs_homer_short_mm10_byid, function(X) { annotatePeak(peak = X, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion = c(-2000, 2000)); })
ann_homer_short_mm10_byid <- lapply(ann_homer_short_mm10_byid, function(X) { X@annoStat })
ann_homer_short_mm10_byid <- t(do.call(cbind, lapply(ann_homer_short_mm10_byid, function(x) structure(x[["Frequency"]], names = as.character(x[["Feature"]])))))
ann_homer_short_mm10_byid <- data.frame(id = rownames(ann_homer_short_mm10_byid), ann_homer_short_mm10_byid, stringsAsFactors = FALSE, check.names = FALSE)

homer_short_mm10_cnts_bychr <- merge(homer_short_mm10_cnts_bychr, ann_homer_short_mm10_byid, by.x = "id", by.y = "id", all.x = TRUE)
write.csv(homer_short_mm10_cnts_bychr, file = "Report/release/DNAzyme/homer_short_mm10.csv", row.names = FALSE)
###########################################################################
#' test if long catalytic sequences are enriched
###########################################################################
pdf("Report/release/DNAzyme/blast_long_nhits_vs_length.pdf")
par(ps = 16)
x <- table(blast_long_hg38$qname)
y <- structure(blast_long_hg38[duplicated(blast_long_hg38$qname), "length"], names = blast_long_hg38[duplicated(blast_long_hg38$qname), "qname"])
# i <- as.vector(x[names(y)]) > 500 | y > 100
plot(y = as.vector(x[names(y)]), x = jitter(y), log = "y", ylab = "# hits (>80% homology) in the genome", xlab = "core catalytic seq. length", main = "name")
# text(y = as.vector(x[names(y)])[i], x = jitter(y[i]), labels = long[match(names(y)[i], long$id), "e"])
x <- table(blast_long_mm10$qname)
y <- structure(blast_long_mm10[duplicated(blast_long_mm10$qname), "length"], names = blast_long_mm10[duplicated(blast_long_mm10$qname), "qname"])
plot(y = as.vector(x[names(y)]), x = jitter(y), log = "y", ylab = "# hits (>80% homology) in the genome", xlab = "core catalytic seq. length", main = "mm10")
dev.off()

## blast_long_hg38 <- subset(blast_long_hg38, !grepl("_", blast_long_hg38$rname))
## grs_blast_long_hg38 <- GRanges(blast_long_hg38$rname, IRanges(blast_long_hg38$pos, blast_long_hg38$pos), strand = "*", seqinfo = Genome$MainSeqInfo$human)

grs_blast_long_hg38 <- GRanges(blast_long_hg38$rname, IRanges(blast_long_hg38$pos, blast_long_hg38$pos), strand = "*", seqinfo = Genome$MainSeqInfo$human)
grs_blast_long_mm10 <- GRanges(blast_long_mm10$rname, IRanges(blast_long_mm10$pos, blast_long_mm10$pos), strand = "*", seqinfo = Genome$MainSeqInfo$human)
grs_dnazyme_hg38 <- c(as(list(longcores = grs_blast_long_hg38), "GRangesList"), grs_homer_short_hg38)
grs_dnazyme_mm10 <- c(as(list(longcores = grs_blast_long_mm10), "GRangesList"), grs_homer_short_mm10)
###########################################################################
SampleInfoFull <- read.csv("Data/release/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[, "SampleID"]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
SampleInfoVirtual <- subset(SampleInfoFull, CompType == "Virtual")
sampleIDs <- SampleInfo[, "SampleID"]
sampleIDsVirtual <- SampleInfoVirtual[, "SampleID"]
Species <- c("human", "mouse", "rat", "none")

PFs_blast_long_hg38 <- sapply(binsizeLabs, function(binsizeLab) {
    X <- unlist(viewSums(Views(cvgs_blast_long_hg38, BinsList[[binsizeLab]][["human"]])))
    names(X) <- unlist(lapply(BinsList[[binsizeLab]][["human"]], as.character))
    X
}, simplify = FALSE)

qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
GRs <- readRDS(sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh))
GRsVirtual <- sapply(sampleIDsVirtual, function(sampleID) {
    message(sampleID)
    sourceIDs <- SampleInfoFull[sampleID, "SourceIDsNoOut"]
    SIDs <- strsplit(sourceIDs, ",")[[1]]
    grs <- GRs[SIDs]
    Reduce(c, grs)
}, simplify = FALSE)
GRsFull <- c(GRs, GRsVirtual)[sampleIDsFull]


GRs_blacklisted_human <- rtracklayer::import("Data/ENCODE/Blacklist/Kundaje/hg38.blacklist.bed")
GRs_blacklisted_human <- Genome$standardizeSeqInfo(GRs_blacklisted_human, seqInfo = Genome$MainSeqInfo$human, prune = TRUE, seqLevels = Genome$MainSeqLevels$human)
GRs_blacklisted_mouse <- rtracklayer::import("Data/ENCODE/Blacklist/Kundaje/mm10.blacklist.bed")
GRs_blacklisted_mouse <- Genome$standardizeSeqInfo(GRs_blacklisted_mouse, seqInfo = Genome$MainSeqInfo$mouse, prune = TRUE, seqLevels = Genome$MainSeqLevels$mouse)

GRs_human <- GRanges(seqnames = seqnames(Genome$MainSeqInfo$human), IRanges(start = 1, end = seqlengths(Genome$MainSeqInfo$human)), seqinfo = Genome$MainSeqInfo$human)
GRs_mouse <- GRanges(seqnames = seqnames(Genome$MainSeqInfo$mouse), IRanges(start = 1, end = seqlengths(Genome$MainSeqInfo$mouse)), seqinfo = Genome$MainSeqInfo$mouse)
GRsNoBlacklisted_human <- setdiff(GRs_human, GRs_blacklisted_human)
GRsNoBlacklisted_mouse <- setdiff(GRs_mouse, GRs_blacklisted_mouse)

SIDsVirtmax_human <- c("K562PositiveAll", "HumanAstroCulturePositiveAll", "HumanNeuronCulturePositiveMerged", "HumanInterneuronCulturePositiveMerged")
bioGroups_human <- c("K562", "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture")
names(SIDsVirtmax_human)  <- bioGroups_human

ChexDnazymeAssoc <- sapply(bioGroups_human, function(x) {
    s <- SIDsVirtmax_human[x]
    k <- GRsFull[[s]]
    sapply(names(grs_dnazyme_hg38), function(y) {
        message(x, " ", y)
        g <- grs_dnazyme_hg38[[y]]
        Genome$testGRsOverlap(query = k, subject = g, universe = GRsNoBlacklisted_human, ignore.strand = TRUE)
    }, simplify = FALSE)
}, simplify = FALSE)
ChexDnazymeAssoc_Log2OddsRatios <- sapply(bioGroups_human, function(x) {
    sapply(names(grs_dnazyme_hg38), function(y) {
        ChexDnazymeAssoc[[x]][[y]][["log2OR"]]
    })
})
ChexDnazymeAssoc_pvals <- sapply(bioGroups_human, function(x) {
    sapply(names(grs_dnazyme_hg38), function(y) {
        ChexDnazymeAssoc[[x]][[y]][["pval"]]
    })
})

write.csv(ChexDnazymeAssoc_Log2OddsRatios, file = "Report/release/DNAzyme/ChexDnazymeAssc_Log2OddsRatios_human.csv") 
write.csv(ChexDnazymeAssoc_pvals, file = "Report/release/DNAzyme/ChexDnazymeAssc_pvals_human.csv") 

 
SIDsVirtmax_mouse <- c("MouseAstroCulturePositiveAll", "MouseNeuronCulturePositiveAll", "MouseNeuronSlicePositiveAll", "MouseInterneuronSlicePositiveMerged")
bioGroups_mouse <- c("MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice")
names(SIDsVirtmax_mouse)  <- bioGroups_mouse

ChexDnazymeAssoc <- sapply(bioGroups_mouse, function(x) {
    s <- SIDsVirtmax_mouse[x]
    k <- GRsFull[[s]]
    sapply(names(grs_dnazyme_mm10), function(y) {
        message(x, " ", y)
        g <- grs_dnazyme_mm10[[y]]
        Genome$testGRsOverlap(query = k, subject = g, universe = GRsNoBlacklisted_mouse, ignore.strand = TRUE)
    }, simplify = FALSE)
}, simplify = FALSE)
ChexDnazymeAssoc_Log2OddsRatios <- sapply(bioGroups_mouse, function(x) {
    sapply(names(grs_dnazyme_mm10), function(y) {
        ChexDnazymeAssoc[[x]][[y]][["log2OR"]]
    })
})
ChexDnazymeAssoc_pvals <- sapply(bioGroups_mouse, function(x) {
    sapply(names(grs_dnazyme_mm10), function(y) {
        ChexDnazymeAssoc[[x]][[y]][["pval"]]
    })
})

write.csv(ChexDnazymeAssoc_Log2OddsRatios, file = "Report/release/DNAzyme/ChexDnazymeAssc_Log2OddsRatios_mouse.csv") 
write.csv(ChexDnazymeAssoc_pvals, file = "Report/release/DNAzyme/ChexDnazymeAssc_pvals_mouse.csv") 

###########################################################################
#' take the top porphyrin metalation DNAzymes for further experiment
###########################################################################
long <- data.table::setDT(long)
## fasta got a problem that 469_PS5 (PS5.M) was converted to 469_PS5 because the space and whatever follows is dropped
blast_long_hg38[qname == "469_PS5", qname := "469_PS5 (PS5.M)"]
blast_long_hg38[, c("perc_mapped", "perc_mismatch") := list(100 * (mismatches + matches) / length, 100 * mismatches / (mismatches + matches))]
blast_long_hg38_porphyrin <- blast_long_hg38[perc_mapped >= 80 & perc_mismatch <= 1 & reaction == "Porphyrin metalation"][order(-perc_mapped, perc_mismatch)]
blast_long_hg38_porphyrin <- long[blast_long_hg38_porphyrin, on = c(id = "qname")]
write.csv(blast_long_hg38_porphyrin, file = "Report/release/DNAzyme/blast_long_hg38_filtered_porphyrin.csv", row.names = FALSE)

blast_long_mm10[, c("perc_mapped", "perc_mismatch") := list(100 * (mismatches + matches) / length, 100 * mismatches / (mismatches + matches))]
blast_long_mm10_porphyrin <- blast_long_mm10[perc_mapped >= 80 & perc_mismatch <= 1 & reaction == "Porphyrin metalation"][order(-perc_mapped, perc_mismatch)]
blast_long_mm10_porphyrin <- long[blast_long_mm10_porphyrin, on = c(id = "qname")]
write.csv(blast_long_mm10_porphyrin, file = "Report/release/DNAzyme/blast_long_mm10_filtered_porphyrin.csv", row.names = FALSE)

###########################################################################
#' take the top RNA-cleavage DNAzymes for further experiment
###########################################################################
blast_long_hg38_rnacleavage <- blast_long_hg38[perc_mapped >= 80 & perc_mismatch <= 1 & reaction == "RNA cleavage"][order(-perc_mapped, perc_mismatch)]
blast_long_hg38_rnacleavage <- long[blast_long_hg38_rnacleavage, on = c(id = "qname")]
write.csv(blast_long_hg38_rnacleavage, file = "Report/release/DNAzyme/blast_long_hg38_filtered_rnacleavage.csv", row.names = FALSE)

blast_long_mm10_rnacleavage <- blast_long_mm10[perc_mapped >= 80 & perc_mismatch <= 1 & reaction == "RNA cleavage"][order(-perc_mapped, perc_mismatch)]
blast_long_mm10_rnacleavage <- long[blast_long_mm10_rnacleavage, on = c(id = "qname")]
write.csv(blast_long_mm10_rnacleavage, file = "Report/release/DNAzyme/blast_long_mm10_filtered_rnacleavage.csv", row.names = FALSE)

###########################################################################
## how many porphyrin DNAzyme sites overlap with CHEX
###########################################################################
SampleInfoFull <- read.csv("Data/release/SampleInfoFullOutAnnotated20201221CV2b.csv", as.is = TRUE, check.names = FALSE)
sampleIDsFull <- rownames(SampleInfoFull) <- SampleInfoFull[, "SampleID"]
SampleInfo <- subset(SampleInfoFull, CompType == "Biol")
SampleInfoVirtual <- subset(SampleInfoFull, CompType == "Virtual")
sampleIDs <- SampleInfo[, "SampleID"]
sampleIDsVirtual <- SampleInfoVirtual[, "SampleID"]
Species <- c("human", "mouse", "rat", "none")
qualOutInPair <- "ABreadCmate5End"
mapqTh <- "ge20_le0.1_strict"
GRs <- readRDS(sprintf("Data/release/PrimingRate/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh))
GRsVirtual <- sapply(sampleIDsVirtual, function(sampleID) {
    message(sampleID)
    sourceIDs <- SampleInfoFull[sampleID, "SourceIDsNoOut"]
    SIDs <- strsplit(sourceIDs, ",")[[1]]
    grs <- GRs[SIDs]
    Reduce(c, grs)
}, simplify = FALSE)
GRsFull <- c(GRs, GRsVirtual)[sampleIDsFull]

grs_chex_human <- c(GRsFull$K562PositiveAll, GRsFull$HumanAstroCulturePositiveAll, GRsFull$HumanNeuronCulturePositiveMerged, GRsFull$HumanInterneuronCulturePositiveMerged)
grs_chex_mouse <- c(GRsFull$MouseAstroCulturePositiveAll, GRsFull$MouseNeuronCulturePositiveAll, GRsFull$MouseNeuronSlicePositiveAll, GRsFull$MouseInterneuronSlicePositiveMerged)

grs_chex_human_ext <- flank(grs_chex_human, width = 1000, start = TRUE, both = TRUE)
grs_chex_mouse_ext <- flank(grs_chex_mouse, width = 1000, start = TRUE, both = TRUE)

blast_long_porphyrin_hg38 <- data.table(openxlsx::read.xlsx("Report/release/DNAzyme/blast_long_filtered_porphyrin.xlsx", sheet = "hg38"))
blast_long_porphyrin_mm10 <- data.table(openxlsx::read.xlsx("Report/release/DNAzyme/blast_long_filtered_porphyrin.xlsx", sheet = "mm10"))

blast_long_porphyrin_hg38[, ins := sapply(cigar, function(x) Genome$parse_cigar(x, "I"))]
blast_long_porphyrin_hg38[, del := sapply(cigar, function(x) Genome$parse_cigar(x, "D"))]

blast_long_porphyrin_mm10[, ins := sapply(cigar, function(x) Genome$parse_cigar(x, "I"))]
blast_long_porphyrin_mm10[, del := sapply(cigar, function(x) Genome$parse_cigar(x, "D"))]

blast_long_porphyrin_hg38[, strand := ifelse(flag == 0, "+", "-")]
blast_long_porphyrin_mm10[, strand := ifelse(flag == 0, "+", "-")]

blast_long_porphyrin_hg38[, geneStrand := ifelse(geneStrand == 1, "+", "-")]
blast_long_porphyrin_mm10[, geneStrand := ifelse(geneStrand == 1, "+", "-")]

grs_long_porphyrin_hg38 <- blast_long_porphyrin_hg38[, GRanges(rname, IRanges(start = pos, end = pos + matches + del - 1), strand = strand, seqinfo = Genome$MainSeqInfo$human, name = id)]
grs_long_porphyrin_mm10 <- blast_long_porphyrin_mm10[, GRanges(rname, IRanges(start = pos, end = pos + matches + del - 1), strand = strand, seqinfo = Genome$MainSeqInfo$mouse, name = id)]

hits_dnazyme_porphyrin_hg38_chex <- findOverlaps(query = grs_long_porphyrin_hg38, subject = grs_chex_human_ext, ignore.strand = FALSE)
hits_dnazyme_porphyrin_mm10_chex <- findOverlaps(query = grs_long_porphyrin_mm10, subject = grs_chex_mouse_ext, ignore.strand = FALSE)

hits_dnazyme_porphyrin_hg38_chexcnts <- tapply(subjectHits(hits_dnazyme_porphyrin_hg38_chex), queryHits(hits_dnazyme_porphyrin_hg38_chex), length)
hits_dnazyme_porphyrin_mm10_chexcnts <- tapply(subjectHits(hits_dnazyme_porphyrin_mm10_chex), queryHits(hits_dnazyme_porphyrin_mm10_chex), length)

hits_dnazyme_porphyrin_hg38_chexcnts <- hits_dnazyme_porphyrin_hg38_chex_cnts[as.character(1:length(grs_long_porphyrin_hg38))]
hits_dnazyme_porphyrin_hg38_chexcnts[is.na(hits_dnazyme_porphyrin_hg38_chexcnts)] <- 0
hits_dnazyme_porphyrin_hg38_chexcnts <- unname(hits_dnazyme_porphyrin_hg38_chexcnts)
mcols(grs_long_porphyrin_hg38)$chexcnts <- hits_dnazyme_porphyrin_hg38_chexcnts

hits_dnazyme_porphyrin_mm10_chexcnts <- hits_dnazyme_porphyrin_mm10_chex_cnts[as.character(1:length(grs_long_porphyrin_mm10))]
hits_dnazyme_porphyrin_mm10_chexcnts[is.na(hits_dnazyme_porphyrin_mm10_chexcnts)] <- 0
hits_dnazyme_porphyrin_mm10_chexcnts <- unname(hits_dnazyme_porphyrin_mm10_chexcnts)
mcols(grs_long_porphyrin_mm10)$chexcnts <- hits_dnazyme_porphyrin_mm10_chexcnts

blast_long_porphyrin_hg38[, start := start(grs_long_porphyrin_hg38)]
blast_long_porphyrin_hg38[, end := end(grs_long_porphyrin_hg38)]

blast_long_porphyrin_mm10[, start := start(grs_long_porphyrin_mm10)]
blast_long_porphyrin_mm10[, end := end(grs_long_porphyrin_mm10)]

blast_long_porphyrin_hg38[, chexcnts := hits_dnazyme_porphyrin_hg38_chexcnts]
blast_long_porphyrin_mm10[, chexcnts := hits_dnazyme_porphyrin_mm10_chexcnts]

rtracklayer::export(grs_long_porphyrin_hg38, con = "Data/release/DNAzyme/grs_long_porphyrin_hg38_core.bed")
rtracklayer::export(grs_long_porphyrin_mm10, con = "Data/release/DNAzyme/grs_long_porphyrin_mm10_core.bed")

grs_long_porphyrin_hg38_upstream <- flank(grs_long_porphyrin_hg38, start = TRUE, width = 100, both = FALSE)
grs_long_porphyrin_hg38_dnstream <- flank(grs_long_porphyrin_hg38, start = FALSE, width = 100, both = FALSE)
rtracklayer::export(grs_long_porphyrin_hg38_upstream, con = "Data/release/DNAzyme/grs_long_porphyrin_hg38_upstream.bed")
rtracklayer::export(grs_long_porphyrin_hg38_dnstream, con = "Data/release/DNAzyme/grs_long_porphyrin_hg38_dnstream.bed")

grs_long_porphyrin_mm10_upstream <- flank(grs_long_porphyrin_mm10, start = TRUE, width = 100, both = FALSE)
grs_long_porphyrin_mm10_dnstream <- flank(grs_long_porphyrin_mm10, start = FALSE, width = 100, both = FALSE)
rtracklayer::export(grs_long_porphyrin_mm10_upstream, con = "Data/release/DNAzyme/grs_long_porphyrin_mm10_upstream.bed")
rtracklayer::export(grs_long_porphyrin_mm10_dnstream, con = "Data/release/DNAzyme/grs_long_porphyrin_mm10_dnstream.bed")

seq_long_porphyrin_hg38_core <- readLines("Data/release/DNAzyme/grs_long_porphyrin_hg38_core.fa")
seq_long_porphyrin_hg38_core <- seq_long_porphyrin_hg38_core[seq(2, length(seq_long_porphyrin_hg38_core), by = 2)]
seq_long_porphyrin_hg38_ext <- readLines("Data/release/DNAzyme/grs_long_porphyrin_hg38_200nt.fa")
seq_long_porphyrin_hg38_ext <- seq_long_porphyrin_hg38_ext[seq(2, length(seq_long_porphyrin_hg38_ext), by = 2)]

blast_long_porphyrin_hg38[, dna_core_seq := seq_long_porphyrin_hg38_core]
blast_long_porphyrin_hg38[, dna_ext_seq := seq_long_porphyrin_hg38_ext]

seq_long_porphyrin_mm10_core <- readLines("Data/release/DNAzyme/grs_long_porphyrin_mm10_core.fa")
seq_long_porphyrin_mm10_core <- seq_long_porphyrin_mm10_core[seq(2, length(seq_long_porphyrin_mm10_core), by = 2)]
seq_long_porphyrin_mm10_ext <- readLines("Data/release/DNAzyme/grs_long_porphyrin_mm10_200nt.fa")
seq_long_porphyrin_mm10_ext <- seq_long_porphyrin_mm10_ext[seq(2, length(seq_long_porphyrin_mm10_ext), by = 2)]

blast_long_porphyrin_mm10[, dna_core_seq := seq_long_porphyrin_mm10_core]
blast_long_porphyrin_mm10[, dna_ext_seq := seq_long_porphyrin_mm10_ext]

fwrite(blast_long_porphyrin_hg38, file = "Report/release/DNAzyme/blast_long_filtered_porphyrin_hg38.csv")
fwrite(blast_long_porphyrin_mm10, file = "Report/release/DNAzyme/blast_long_filtered_porphyrin_mm10.csv")
