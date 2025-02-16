## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
source("Source/functions.R")

## Genome assembly defaults to hg38 and mm10
.RERUN <- FALSE
if (.RERUN) { 
    warning("Rerunning everything...")
}

## I. Organism level annotations to get gene ID, Gene Ontology mapping table
OrgDbs <- list(human = org.Hs.eg.db, mouse = org.Mm.eg.db, rat = org.Rn.eg.db)

## II. Ensembl gene model; has symbols in Genes
## EnsDb.Rnorvegicus.v79 is rn5, so we cannot use what bioconductor provides. 
## Let's build our own from the GTF file downloaded from ENSEMBL.
if (.RERUN) {
    EnsGTFRn6 <- ensDbFromGtf(gtf = "Data/Annotation/Ensembl/Rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.97.gtf.gz", outfile = "Data/Annotation/Ensembl/Rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.97.sqlite", organism = "Rattus_norvegicus", genomeVersion = "rn6", version = 97)
} else {
    EnsGTFRn6 <- "Data/Annotation/Ensembl/Rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.97.sqlite"
}
EnsHg38 <- EnsDb.Hsapiens.v86
EnsMm10 <- EnsDb.Mmusculus.v79
EnsRn6 <- EnsDb(EnsGTFRn6)

seqlevelsStyle(EnsHg38) <- "UCSC"
EnsHg38SeqLevels <- seqlevels(EnsHg38)
EnsHg38SeqLengths <- seqlengths(EnsHg38)
EnsHg38SeqInfo <- Seqinfo(seqnames = EnsHg38SeqLevels, seqlengths = EnsHg38SeqLengths[EnsHg38SeqLevels], isCircular = ifelse(EnsHg38SeqLevels == "chrM", TRUE, FALSE), genome = "hg38")

seqlevelsStyle(EnsMm10) <- "UCSC"
EnsMm10SeqLevels <- seqlevels(EnsMm10)
EnsMm10SeqLengths <- seqlengths(EnsMm10)
EnsMm10SeqInfo <- Seqinfo(seqnames = EnsMm10SeqLevels, seqlengths = EnsMm10SeqLengths[EnsMm10SeqLevels], isCircular = ifelse(EnsMm10SeqLevels == "chrM", TRUE, FALSE), genome = "mm10")

seqlevelsStyle(EnsRn6) <- "UCSC"
EnsRn6SeqLevels <- seqlevels(EnsRn6)
EnsRn6SeqLengths <- seqlengths(EnsRn6)
EnsRn6SeqInfo <- Seqinfo(seqnames = EnsRn6SeqLevels, seqlengths = EnsRn6SeqLengths[EnsRn6SeqLevels], isCircular = ifelse(EnsRn6SeqLevels == "chrM", TRUE, FALSE), genome = "rn6")

if (.RERUN) {
###########################################################################
## Human
###########################################################################
## Make chrM circular, otherwise GenomicRange will complain later
    EnsHg38Genes <- genes(EnsHg38)
    EnsHg38Genes <- Genome$standardizeSeqInfo(EnsHg38Genes, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)

## Transcripts
    EnsHg38Transcripts <- transcripts(EnsHg38)
    EnsHg38Transcripts$symbol <- mcols(EnsHg38Genes[mcols(EnsHg38Transcripts)[["gene_id"]]])[["symbol"]]
    EnsHg38Transcripts <- Genome$standardizeSeqInfo(EnsHg38Transcripts, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)

## Exons
    EnsHg38ExonsByTranscript <- exonsBy(EnsHg38, by = "tx")
    EnsHg38ExonsByTranscript <- Genome$standardizeSeqInfo(EnsHg38ExonsByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38ExonsByGene <- exonsBy(EnsHg38, "gene")
    EnsHg38ExonsByGene <- Genome$standardizeSeqInfo(EnsHg38ExonsByGene, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
# EnsHg38ExonsByGene[!duplicated(EnsHg38ExonsByGene)] is same as unique(EnsHg38ExonsByGene)
    EnsHg38ExonsByGene <- unique(EnsHg38ExonsByGene)

## Introns
    EnsHg38IntronsByTranscript <- Genome$getIntronsByTranscript(EnsHg38ExonsByTranscript, transcripts = EnsHg38Transcripts)
    EnsHg38IntronsByGene <- split(unlist(EnsHg38IntronsByTranscript, use.names = FALSE), f = mcols(EnsHg38Transcripts[names(unlist(EnsHg38IntronsByTranscript))])[["gene_id"]])
    EnsHg38IntronsByGene <- unique(EnsHg38IntronsByGene)
length(EnsHg38IntronsByGene)
#" [1] 39639
## We lost ~16000 genes because they have no intron at all.

## Promoters; Note `promoters` returns a GRanges object, not a GRangesList
## When `split` works on a GRanges, it returns a GRangesList object, hence no need to further convert it to GRangesList
    EnsHg38PromotersByTranscript <- trim(promoters(EnsHg38, upstream = 2000, downstream = 200))
    EnsHg38PromotersByTranscript <- Genome$standardizeSeqInfo(EnsHg38PromotersByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38PromotersByGene <- split(EnsHg38PromotersByTranscript, f = mcols(EnsHg38PromotersByTranscript)[["gene_id"]])
    EnsHg38PromotersByGene <- unique(EnsHg38PromotersByGene)
    EnsHg38PromotersByTranscript <- split(EnsHg38PromotersByTranscript, f = mcols(EnsHg38PromotersByTranscript)[["tx_id"]])

## 5'UTRs
    EnsHg38FiveUTRsByTranscript <- fiveUTRsByTranscript(EnsHg38)
    EnsHg38FiveUTRsByTranscript <- Genome$standardizeSeqInfo(EnsHg38FiveUTRsByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38FiveUTRsByGene <- split(unlist(EnsHg38FiveUTRsByTranscript), f = mcols(EnsHg38Transcripts[names(unlist(EnsHg38FiveUTRsByTranscript))])[["gene_id"]])
    EnsHg38FiveUTRsByGene <- unique(EnsHg38FiveUTRsByGene)

## 3'UTRs
    EnsHg38ThreeUTRsByTranscript <- threeUTRsByTranscript(EnsHg38)
    EnsHg38ThreeUTRsByTranscript <- Genome$standardizeSeqInfo(EnsHg38ThreeUTRsByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38ThreeUTRsByGene <- split(unlist(EnsHg38ThreeUTRsByTranscript), f = mcols(EnsHg38Transcripts[names(unlist(EnsHg38ThreeUTRsByTranscript))])[["gene_id"]])
    EnsHg38ThreeUTRsByGene <- unique(EnsHg38ThreeUTRsByGene)

## CDS
    EnsHg38CDSsByTranscript <- cdsBy(EnsHg38, "tx")
    EnsHg38CDSsByTranscript <- Genome$standardizeSeqInfo(EnsHg38CDSsByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38CDSsByGene <- cdsBy(EnsHg38, "gene")
    EnsHg38CDSsByGene <- Genome$standardizeSeqInfo(EnsHg38CDSsByGene, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38CDSsByGene <- unique(EnsHg38CDSsByGene)

## Downstream
    EnsHg38DownstreamByTranscript <- Genome$getDownstream(EnsHg38Transcripts)
    EnsHg38DownstreamByTranscript <- Genome$standardizeSeqInfo(EnsHg38DownstreamByTranscript, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38DownstreamByGene <- split(EnsHg38DownstreamByTranscript, f = mcols(EnsHg38Transcripts[names(EnsHg38DownstreamByTranscript)])[["gene_id"]])
    EnsHg38DownstreamByGene <- unique(EnsHg38DownstreamByGene)
    EnsHg38DownstreamByTranscript <- split(EnsHg38DownstreamByTranscript, f = mcols(EnsHg38DownstreamByTranscript)[["tx_id"]])

## Intergenic
    EnsHg38Intergenic <- Genome$getGenomicFeaturesIntergenic(EnsHg38Genes)
    EnsHg38Intergenic <- Genome$standardizeSeqInfo(EnsHg38Intergenic, seqInfo = EnsHg38SeqInfo, seqLevels = EnsHg38SeqLevels)
    EnsHg38Intergenic$name <- paste0("EnsHg38Int", seq(length(EnsHg38Intergenic)))

## TSS
    EnsHg38TSSsByTranscript <- resize(EnsHg38Transcripts, width = 1)
## EnsHg38TSSsByGene <- resize(EnsHg38Genes, width = 1) 
## Note, this method only gets one TSS per gene, which is the most upstream TSS, even it is from an less confident transcript.
## Here we just have TSS from all transcript and decide later which TSS as the representative for a gene. 
    EnsHg38TSSsByGene <- split(EnsHg38TSSsByTranscript, f = mcols(EnsHg38TSSsByTranscript)[["gene_id"]])
    EnsHg38TSSsByGene <- unique(EnsHg38TSSsByGene)
    EnsHg38TSSsByTranscript <- split(EnsHg38TSSsByTranscript, f = mcols(EnsHg38TSSsByTranscript)[["tx_id"]])

## Flank 5kb both sides
    EnsHg38Flank5kByTranscript <- trim(flank(EnsHg38Transcripts, both = TRUE, width = 2500))
    EnsHg38Flank5kByGene <- split(EnsHg38Flank5kByTranscript, f = mcols(EnsHg38Flank5kByTranscript)[["gene_id"]])
    EnsHg38Flank5kByGene <- unique(EnsHg38Flank5kByGene)
    EnsHg38Flank5kByTranscript <- split(EnsHg38Flank5kByTranscript, f = mcols(EnsHg38Flank5kByTranscript)[["tx_id"]])

## Flank 5kb upstream
    EnsHg38FlankUp5kByTranscript <- trim(flank(EnsHg38Transcripts, both = FALSE, width = 5000))
    EnsHg38FlankUp5kByGene <- split(EnsHg38FlankUp5kByTranscript, f = mcols(EnsHg38FlankUp5kByTranscript)[["gene_id"]])
    EnsHg38FlankUp5kByGene <- unique(EnsHg38FlankUp5kByGene)
    EnsHg38FlankUp5kByTranscript <- split(EnsHg38FlankUp5kByTranscript, f = mcols(EnsHg38FlankUp5kByTranscript)[["tx_id"]])

## Flank 5kb downstream
    EnsHg38FlankDn5kByTranscript <- trim(flank(EnsHg38Transcripts, both = FALSE, width = -5000))
    EnsHg38FlankDn5kByGene <- split(EnsHg38FlankDn5kByTranscript, f = mcols(EnsHg38FlankDn5kByTranscript)[["gene_id"]])
    EnsHg38FlankDn5kByGene <- unique(EnsHg38FlankDn5kByGene)
    EnsHg38FlankDn5kByTranscript <- split(EnsHg38FlankDn5kByTranscript, f = mcols(EnsHg38FlankDn5kByTranscript)[["tx_id"]])

###########################################################################
## Mouse
###########################################################################
    EnsMm10Genes <- genes(EnsMm10)
    EnsMm10Genes <- Genome$standardizeSeqInfo(EnsMm10Genes, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)

    EnsMm10Transcripts <- transcripts(EnsMm10)
    EnsMm10Transcripts$symbol <- mcols(EnsMm10Genes[mcols(EnsMm10Transcripts)[["gene_id"]]])[["symbol"]]
    EnsMm10Transcripts <- Genome$standardizeSeqInfo(EnsMm10Transcripts, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)

## Exons
    EnsMm10ExonsByTranscript <- exonsBy(EnsMm10, by = "tx")
    EnsMm10ExonsByTranscript <- Genome$standardizeSeqInfo(EnsMm10ExonsByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10ExonsByGene <- exonsBy(EnsMm10, "gene")
    EnsMm10ExonsByGene <- Genome$standardizeSeqInfo(EnsMm10ExonsByGene, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
# EnsMm10ExonsByGene[!duplicated(EnsMm10ExonsByGene)] is same as unique(EnsMm10ExonsByGene)
    EnsMm10ExonsByGene <- unique(EnsMm10ExonsByGene)

## Introns
    EnsMm10IntronsByTranscript <- Genome$getIntronsByTranscript(EnsMm10ExonsByTranscript, transcripts = EnsMm10Transcripts)
    EnsMm10IntronsByGene <- split(unlist(EnsMm10IntronsByTranscript, use.names = FALSE), f = mcols(EnsMm10Transcripts[names(unlist(EnsMm10IntronsByTranscript))])[["gene_id"]])
    EnsMm10IntronsByGene <- unique(EnsMm10IntronsByGene)

## Promoters
    EnsMm10PromotersByTranscript <- trim(promoters(EnsMm10, upstream = 2000, downstream = 200))
    EnsMm10PromotersByTranscript <- Genome$standardizeSeqInfo(EnsMm10PromotersByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10PromotersByGene <- split(EnsMm10PromotersByTranscript, f = mcols(EnsMm10PromotersByTranscript)[["gene_id"]])
    EnsMm10PromotersByGene <- unique(EnsMm10PromotersByGene)
    EnsMm10PromotersByTranscript <- split(EnsMm10PromotersByTranscript, f = mcols(EnsMm10PromotersByTranscript)[["tx_id"]])

## 5'UTRs
    EnsMm10FiveUTRsByTranscript <- fiveUTRsByTranscript(EnsMm10)
    EnsMm10FiveUTRsByTranscript <- Genome$standardizeSeqInfo(EnsMm10FiveUTRsByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10FiveUTRsByGene <- split(unlist(EnsMm10FiveUTRsByTranscript), f = mcols(EnsMm10Transcripts[names(unlist(EnsMm10FiveUTRsByTranscript))])[["gene_id"]])
    EnsMm10FiveUTRsByGene <- unique(EnsMm10FiveUTRsByGene)

## 3'UTRs
    EnsMm10ThreeUTRsByTranscript <- threeUTRsByTranscript(EnsMm10)
    EnsMm10ThreeUTRsByTranscript <- Genome$standardizeSeqInfo(EnsMm10ThreeUTRsByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10ThreeUTRsByGene <- split(unlist(EnsMm10ThreeUTRsByTranscript), f = mcols(EnsMm10Transcripts[names(unlist(EnsMm10ThreeUTRsByTranscript))])[["gene_id"]])
    EnsMm10ThreeUTRsByGene <- unique(EnsMm10ThreeUTRsByGene)

## CDS
    EnsMm10CDSsByTranscript <- cdsBy(EnsMm10, "tx")
    EnsMm10CDSsByTranscript <- Genome$standardizeSeqInfo(EnsMm10CDSsByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10CDSsByGene <- cdsBy(EnsMm10, "gene")
    EnsMm10CDSsByGene <- Genome$standardizeSeqInfo(EnsMm10CDSsByGene, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10CDSsByGene <- unique(EnsMm10CDSsByGene)

## Downstream
    EnsMm10DownstreamByTranscript <- Genome$getDownstream(EnsMm10Transcripts)
    EnsMm10DownstreamByTranscript <- Genome$standardizeSeqInfo(EnsMm10DownstreamByTranscript, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10DownstreamByGene <- split(EnsMm10DownstreamByTranscript, f = mcols(EnsMm10Transcripts[names(EnsMm10DownstreamByTranscript)])[["gene_id"]])
    EnsMm10DownstreamByGene <- unique(EnsMm10DownstreamByGene)
    EnsMm10DownstreamByTranscript <- split(EnsMm10DownstreamByTranscript, f = mcols(EnsMm10DownstreamByTranscript)[["tx_id"]])

## Intergenic
    EnsMm10Intergenic <- Genome$getGenomicFeaturesIntergenic(EnsMm10Genes)
    EnsMm10Intergenic <- Genome$standardizeSeqInfo(EnsMm10Intergenic, seqInfo = EnsMm10SeqInfo, seqLevels = EnsMm10SeqLevels)
    EnsMm10Intergenic$name <- paste0("EnsMm10Int", seq(length(EnsMm10Intergenic)))

## TSS
    EnsMm10TSSsByTranscript <- resize(EnsMm10Transcripts, width = 1)
    EnsMm10TSSsByGene <- split(EnsMm10TSSsByTranscript, f = mcols(EnsMm10TSSsByTranscript)[["gene_id"]])
    EnsMm10TSSsByGene <- unique(EnsMm10TSSsByGene)
    EnsMm10TSSsByTranscript <- split(EnsMm10TSSsByTranscript, f = mcols(EnsMm10TSSsByTranscript)[["tx_id"]])

## Flank 5kb
    EnsMm10Flank5kByTranscript <- trim(flank(EnsMm10Transcripts, both = TRUE, width = 2500))
    EnsMm10Flank5kByGene <- split(EnsMm10Flank5kByTranscript, f = mcols(EnsMm10Flank5kByTranscript)[["gene_id"]])
    EnsMm10Flank5kByGene <- unique(EnsMm10Flank5kByGene)
    EnsMm10Flank5kByTranscript <- split(EnsMm10Flank5kByTranscript, f = mcols(EnsMm10Flank5kByTranscript)[["tx_id"]])

## Flank 5kb upstream
    EnsMm10FlankUp5kByTranscript <- trim(flank(EnsMm10Transcripts, both = FALSE, width = 5000))
    EnsMm10FlankUp5kByGene <- split(EnsMm10FlankUp5kByTranscript, f = mcols(EnsMm10FlankUp5kByTranscript)[["gene_id"]])
    EnsMm10FlankUp5kByGene <- unique(EnsMm10FlankUp5kByGene)
    EnsMm10FlankUp5kByTranscript <- split(EnsMm10FlankUp5kByTranscript, f = mcols(EnsMm10FlankUp5kByTranscript)[["tx_id"]])

## Flank 5kb downstream
    EnsMm10FlankDn5kByTranscript <- trim(flank(EnsMm10Transcripts, both = FALSE, width = -5000))
    EnsMm10FlankDn5kByGene <- split(EnsMm10FlankDn5kByTranscript, f = mcols(EnsMm10FlankDn5kByTranscript)[["gene_id"]])
    EnsMm10FlankDn5kByGene <- unique(EnsMm10FlankDn5kByGene)
    EnsMm10FlankDn5kByTranscript <- split(EnsMm10FlankDn5kByTranscript, f = mcols(EnsMm10FlankDn5kByTranscript)[["tx_id"]])
###########################################################################
## Rat
###########################################################################
    EnsRn6Genes <- genes(EnsRn6)
    EnsRn6Genes <- Genome$standardizeSeqInfo(EnsRn6Genes, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)

    EnsRn6Transcripts <- transcripts(EnsRn6)
    EnsRn6Transcripts$symbol <- mcols(EnsRn6Genes[mcols(EnsRn6Transcripts)[["gene_id"]]])[["symbol"]]
    EnsRn6Transcripts <- Genome$standardizeSeqInfo(EnsRn6Transcripts, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)

## Exons
    EnsRn6ExonsByTranscript <- exonsBy(EnsRn6, by = "tx")
    EnsRn6ExonsByTranscript <- Genome$standardizeSeqInfo(EnsRn6ExonsByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6ExonsByGene <- exonsBy(EnsRn6, "gene")
    EnsRn6ExonsByGene <- Genome$standardizeSeqInfo(EnsRn6ExonsByGene, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
# EnsRn6ExonsByGene[!duplicated(EnsRn6ExonsByGene)] is same as unique(EnsRn6ExonsByGene)
    EnsRn6ExonsByGene <- unique(EnsRn6ExonsByGene)

## Introns
    EnsRn6IntronsByTranscript <- Genome$getIntronsByTranscript(EnsRn6ExonsByTranscript, transcripts = EnsRn6Transcripts)
    EnsRn6IntronsByGene <- split(unlist(EnsRn6IntronsByTranscript, use.names = FALSE), f = mcols(EnsRn6Transcripts[names(unlist(EnsRn6IntronsByTranscript))])[["gene_id"]])
    EnsRn6IntronsByGene <- unique(EnsRn6IntronsByGene)

## Promoters
    EnsRn6PromotersByTranscript <- trim(promoters(EnsRn6, upstream = 2000, downstream = 200))
    EnsRn6PromotersByTranscript <- Genome$standardizeSeqInfo(EnsRn6PromotersByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6PromotersByGene <- split(EnsRn6PromotersByTranscript, f = mcols(EnsRn6PromotersByTranscript)[["gene_id"]])
    EnsRn6PromotersByGene <- unique(EnsRn6PromotersByGene)
    EnsRn6PromotersByTranscript <- split(EnsRn6PromotersByTranscript, f = mcols(EnsRn6PromotersByTranscript)[["tx_id"]])

## 5'UTRs
    EnsRn6FiveUTRsByTranscript <- fiveUTRsByTranscript(EnsRn6)
    EnsRn6FiveUTRsByTranscript <- Genome$standardizeSeqInfo(EnsRn6FiveUTRsByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6FiveUTRsByGene <- split(unlist(EnsRn6FiveUTRsByTranscript), f = mcols(EnsRn6Transcripts[names(unlist(EnsRn6FiveUTRsByTranscript))])[["gene_id"]])
    EnsRn6FiveUTRsByGene <- unique(EnsRn6FiveUTRsByGene)

## 3'UTRs
    EnsRn6ThreeUTRsByTranscript <- threeUTRsByTranscript(EnsRn6)
    EnsRn6ThreeUTRsByTranscript <- Genome$standardizeSeqInfo(EnsRn6ThreeUTRsByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6ThreeUTRsByGene <- split(unlist(EnsRn6ThreeUTRsByTranscript), f = mcols(EnsRn6Transcripts[names(unlist(EnsRn6ThreeUTRsByTranscript))])[["gene_id"]])
    EnsRn6ThreeUTRsByGene <- unique(EnsRn6ThreeUTRsByGene)

## CDS
    EnsRn6CDSsByTranscript <- cdsBy(EnsRn6, "tx")
    EnsRn6CDSsByTranscript <- Genome$standardizeSeqInfo(EnsRn6CDSsByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6CDSsByGene <- cdsBy(EnsRn6, "gene")
    EnsRn6CDSsByGene <- Genome$standardizeSeqInfo(EnsRn6CDSsByGene, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6CDSsByGene <- unique(EnsRn6CDSsByGene)

## Downstream
    EnsRn6DownstreamByTranscript <- Genome$getDownstream(EnsRn6Transcripts)
    EnsRn6DownstreamByTranscript <- Genome$standardizeSeqInfo(EnsRn6DownstreamByTranscript, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6DownstreamByGene <- split(EnsRn6DownstreamByTranscript, f = mcols(EnsRn6Transcripts[names(EnsRn6DownstreamByTranscript)])[["gene_id"]])
    EnsRn6DownstreamByGene <- unique(EnsRn6DownstreamByGene)
    EnsRn6DownstreamByTranscript <- split(EnsRn6DownstreamByTranscript, f = mcols(EnsRn6DownstreamByTranscript)[["tx_id"]])

## Intergenic
    EnsRn6Intergenic <- Genome$getGenomicFeaturesIntergenic(EnsRn6Genes)
    EnsRn6Intergenic <- Genome$standardizeSeqInfo(EnsRn6Intergenic, seqInfo = EnsRn6SeqInfo, seqLevels = EnsRn6SeqLevels)
    EnsRn6Intergenic$name <- paste0("EnsRn6Int", seq(length(EnsRn6Intergenic)))

## TSS
    EnsRn6TSSsByTranscript <- resize(EnsRn6Transcripts, width = 1)
    EnsRn6TSSsByGene <- split(EnsRn6TSSsByTranscript, f = mcols(EnsRn6TSSsByTranscript)[["gene_id"]])
    EnsRn6TSSsByGene <- unique(EnsRn6TSSsByGene)
    EnsRn6TSSsByTranscript <- split(EnsRn6TSSsByTranscript, f = mcols(EnsRn6TSSsByTranscript)[["tx_id"]])

## Flank 5kb
    EnsRn6Flank5kByTranscript <- trim(flank(EnsRn6Transcripts, both = TRUE, width = 2500))
    EnsRn6Flank5kByGene <- split(EnsRn6Flank5kByTranscript, f = mcols(EnsRn6Flank5kByTranscript)[["gene_id"]])
    EnsRn6Flank5kByGene <- unique(EnsRn6Flank5kByGene)
    EnsRn6Flank5kByTranscript <- split(EnsRn6Flank5kByTranscript, f = mcols(EnsRn6Flank5kByTranscript)[["tx_id"]])

## Flank 5kb upstream
    EnsRn6FlankUp5kByTranscript <- trim(flank(EnsRn6Transcripts, both = FALSE, width = 5000))
    EnsRn6FlankUp5kByGene <- split(EnsRn6FlankUp5kByTranscript, f = mcols(EnsRn6FlankUp5kByTranscript)[["gene_id"]])
    EnsRn6FlankUp5kByGene <- unique(EnsRn6FlankUp5kByGene)
    EnsRn6FlankUp5kByTranscript <- split(EnsRn6FlankUp5kByTranscript, f = mcols(EnsRn6FlankUp5kByTranscript)[["tx_id"]])

## Flank 5kb downstream
    EnsRn6FlankDn5kByTranscript <- trim(flank(EnsRn6Transcripts, both = FALSE, width = -5000))
    EnsRn6FlankDn5kByGene <- split(EnsRn6FlankDn5kByTranscript, f = mcols(EnsRn6FlankDn5kByTranscript)[["gene_id"]])
    EnsRn6FlankDn5kByGene <- unique(EnsRn6FlankDn5kByGene)
    EnsRn6FlankDn5kByTranscript <- split(EnsRn6FlankDn5kByTranscript, f = mcols(EnsRn6FlankDn5kByTranscript)[["tx_id"]])

###########################################################################
## Merge human and mouse to one object
###########################################################################
    EnsGenes <- list(human = EnsHg38Genes, mouse = EnsMm10Genes, rat = EnsRn6Genes)
    EnsTranscripts <- list(human = EnsHg38Transcripts, mouse = EnsMm10Transcripts, rat = EnsRn6Transcripts)
    EnsExonsByTranscript <- list(human = EnsHg38ExonsByTranscript, mouse = EnsMm10ExonsByTranscript, rat = EnsRn6ExonsByTranscript)
    EnsExonsByGene <- list(human = EnsHg38ExonsByGene, mouse = EnsMm10ExonsByGene, rat = EnsRn6ExonsByGene)
    EnsCDSsByTranscript <- list(human = EnsHg38CDSsByTranscript, mouse = EnsMm10CDSsByTranscript, rat = EnsRn6CDSsByTranscript)
    EnsCDSsByGene <- list(human = EnsHg38CDSsByGene, mouse = EnsMm10CDSsByGene, rat = EnsRn6CDSsByGene)
    EnsIntronsByTranscript <- list(human = EnsHg38IntronsByTranscript, mouse = EnsMm10IntronsByTranscript, rat = EnsRn6IntronsByTranscript)
    EnsIntronsByGene <- list(human = EnsHg38IntronsByGene, mouse = EnsMm10IntronsByGene, rat = EnsRn6IntronsByGene)
    EnsPromotersByTranscript <- list(human = EnsHg38PromotersByTranscript, mouse = EnsMm10PromotersByTranscript, rat = EnsRn6PromotersByTranscript)
    EnsPromotersByGene <- list(human = EnsHg38PromotersByGene, mouse = EnsMm10PromotersByGene, rat = EnsRn6PromotersByGene)
    EnsFiveUTRsByTranscript <- list(human = EnsHg38FiveUTRsByTranscript, mouse = EnsMm10FiveUTRsByTranscript, rat = EnsRn6FiveUTRsByTranscript)
    EnsFiveUTRsByGene <- list(human = EnsHg38FiveUTRsByGene, mouse = EnsMm10FiveUTRsByGene, rat = EnsRn6FiveUTRsByGene)
    EnsThreeUTRsByTranscript <- list(human = EnsHg38ThreeUTRsByTranscript, mouse = EnsMm10ThreeUTRsByTranscript, rat = EnsRn6ThreeUTRsByTranscript)
    EnsThreeUTRsByGene <- list(human = EnsHg38ThreeUTRsByGene, mouse = EnsMm10ThreeUTRsByGene, rat = EnsRn6ThreeUTRsByGene)
    EnsDownstreamByTranscript <- list(human = EnsHg38DownstreamByTranscript, mouse = EnsMm10DownstreamByTranscript, rat = EnsRn6DownstreamByTranscript)
    EnsDownstreamByGene <- list(human = EnsHg38DownstreamByGene, mouse = EnsMm10DownstreamByGene, rat = EnsRn6DownstreamByGene)
    EnsIntergenic <- list(human = EnsHg38Intergenic, mouse = EnsMm10Intergenic, rat = EnsRn6Intergenic)
    EnsTSSsByTranscript <- list(human = EnsHg38TSSsByTranscript, mouse = EnsMm10TSSsByTranscript, rat = EnsRn6TSSsByTranscript)
    EnsTSSsByGene <- list(human = EnsHg38TSSsByGene, mouse = EnsMm10TSSsByGene, rat = EnsRn6TSSsByGene)
    EnsFlank5kByTranscript <- list(human = EnsHg38Flank5kByTranscript, mouse = EnsMm10Flank5kByTranscript, rat = EnsRn6Flank5kByTranscript)
    EnsFlank5kByGene <- list(human = EnsHg38Flank5kByGene, mouse = EnsMm10Flank5kByGene, rat = EnsRn6Flank5kByGene)
    EnsFlankUp5kByTranscript <- list(human = EnsHg38FlankUp5kByTranscript, mouse = EnsMm10FlankUp5kByTranscript, rat = EnsRn6FlankUp5kByTranscript)
    EnsFlankUp5kByGene <- list(human = EnsHg38FlankUp5kByGene, mouse = EnsMm10FlankUp5kByGene, rat = EnsRn6FlankUp5kByGene)
    EnsFlankDn5kByTranscript <- list(human = EnsHg38FlankDn5kByTranscript, mouse = EnsMm10FlankDn5kByTranscript, rat = EnsRn6FlankDn5kByTranscript)
    EnsFlankDn5kByGene <- list(human = EnsHg38FlankDn5kByGene, mouse = EnsMm10FlankDn5kByGene, rat = EnsRn6FlankDn5kByGene)

## Add another genic region to EnsFeatures, "GeneExt", which is defined as 2.5kb upstream a gene's TSS and 2.5kb downstream a genes TES (i.e. end of last exon)
    EnsGenesExt <- list(human = EnsHg38Genes + 2500, mouse = EnsMm10Genes + 2500, rat = EnsRn6Genes + 2500)

## Combine all Ens features into an object
    EnsFeatures <- list(
        Gene = EnsGenes, Transcript = EnsTranscripts, 
        CDSByTranscript = EnsCDSsByTranscript, 
        CDSByGene = EnsCDSsByGene, 
        ExonByTranscript = EnsExonsByTranscript,
        ExonByGene = EnsExonsByGene,
        IntronByTranscript = EnsIntronsByTranscript, 
        IntronByGene = EnsIntronsByGene, 
        PromoterByTranscript = EnsPromotersByTranscript, 
        PromoterByGene = EnsPromotersByGene, 
        FiveUTRByTranscript = EnsFiveUTRsByTranscript, 
        FiveUTRByGene = EnsFiveUTRsByGene, 
        ThreeUTRByTranscript = EnsThreeUTRsByTranscript, 
        ThreeUTRByGene = EnsThreeUTRsByGene, 
        DownstreamByTranscript = EnsDownstreamByTranscript,
        DownstreamByGene = EnsDownstreamByGene,
        Intergenic = EnsIntergenic, 
        TSSByTranscript = EnsTSSsByTranscript, 
        TSSByGene = EnsTSSsByGene,
        Flank5kByTranscript = EnsFlank5kByTranscript,
        Flank5kByGene = EnsFlank5kByGene,
        FlankUp5kByTranscript = EnsFlankUp5kByTranscript,
        FlankUp5kByGene = EnsFlankUp5kByGene,
        FlankDn5kByTranscript = EnsFlankDn5kByTranscript,
        FlankDn5kByGene = EnsFlankDn5kByGene,
        GeneExt = EnsGenesExt
    )
    saveRDS(EnsFeatures, file = "Data/GenomicFeatures/EnsFeatures.RDS")
} else {
    EnsFeatures <- readRDS("Data/GenomicFeatures/EnsFeatures.RDS")
}

###########################################################################
## Combind human, mouse and rat into one object for EnsDb
###########################################################################
EnsSeqInfo <- list(human = EnsHg38SeqInfo, mouse = EnsMm10SeqInfo, rat = EnsRn6SeqInfo)
EnsSeqLengths <- list(human = EnsHg38SeqLengths, mouse = EnsMm10SeqLengths, rat = EnsRn6SeqLengths)
EnsSeqLevels <- list(human = EnsHg38SeqLevels, mouse = EnsMm10SeqLevels, rat = EnsRn6SeqLevels)

###########################################################################
#" Load cytoband information
###########################################################################
if (.RERUN) {
    UCSCHg38Cytobands <- import("Data/Genome/Human/UCSC/hg38.cytoBandIdeo.fromTable.gtf")
    UCSCMm10Cytobands <- import("Data/Genome/Mouse/UCSC/mm10.cytoBandIdeo.fromTable.gtf")
    UCSCRn6Cytobands <- import("Data/Genome/Rat/UCSC/rn6.cytoBandIdeo.fromTable.gtf")
    UCSCCytobands <- list(human = UCSCHg38Cytobands, mouse = UCSCMm10Cytobands, rat = UCSCRn6Cytobands)
    UCSCCytobands <- sapply(c("human", "mouse", "rat"), function(species) Genome$standardizeSeqInfo(UCSCCytobands[[species]], seqInfo = Genome$MainSeqInfo[[species]], seqLevels = Genome$MainSeqLevels[[species]], prune = TRUE), simplify = FALSE)
    saveRDS(UCSCCytobands, file = "Data/GenomicFeatures/UCSCCytobands.RDS")
} else {
    UCSCCytobands <- readRDS("Data/GenomicFeatures/UCSCCytobands.RDS")
}

###########################################################################
## Load ENCODE blacklisted regions
###########################################################################
if (.RERUN) {
    BlacklistedGRs <- list(
        human = sort(import(con = "Data/ENCODE/Blacklist/Kundaje/hg38.blacklist.bed", format = "bed")),
        mouse = sort(import(con = "Data/ENCODE/Blacklist/Kundaje/mm10.blacklist.bed", format = "bed")),
        # Rat has an empty bed because ENCODE does not have a blacklist.
        rat = sort(import(con = "Data/ENCODE/Blacklist/Kundaje/rn6.blacklist.bed", format = "bed"))
    )
    BlacklistedGRs <- sapply(c("human", "mouse", "rat"), function(species) Genome$standardizeSeqInfo(BlacklistedGRs[[species]], seqInfo = Genome$MainSeqInfo[[species]], seqLevels = Genome$MainSeqLevels[[species]], prune = TRUE), simplify = FALSE)
    saveRDS(BlacklistedGRs, file = "Data/GenomicFeatures/BlacklistedGRs.RDS")
} else {
    BlacklistedGRs <- readRDS("Data/GenomicFeatures/BlacklistedGRs.RDS")
}

###########################################################################
## Load UCSC repeat masked regions
###########################################################################
if (.RERUN) {
    UCSCHg38Repeats <- import("Data/Genome/Human/UCSC/hg38.rmsk.fromTable.gtf")
    UCSCMm10Repeats <- import("Data/Genome/Mouse/UCSC/mm10.rmsk.fromTable.gtf")
    UCSCRn6Repeats <- import("Data/Genome/Rat/UCSC/rn6.rmsk.fromTable.gtf")
    UCSCRepeats <- list(human = UCSCHg38Repeats, mouse = UCSCMm10Repeats, rat = UCSCRn6Repeats)
    UCSCRepeats <- sapply(c("human", "mouse", "rat"), function(species) Genome$standardizeSeqInfo(UCSCRepeats[[species]], seqInfo = Genome$MainSeqInfo[[species]], seqLevels = Genome$MainSeqLevels[[species]], prune = TRUE), simplify = FALSE)
    saveRDS(UCSCRepeats, file = "Data/GenomicFeatures/UCSCRepeats.RDS")
} else {
    UCSCRepeats <- readRDS("Data/GenomicFeatures/UCSCRepeats.RDS")
}

###########################################################################
## Load CpG islands
###########################################################################
if (.RERUN) {
    UCSCCpGIslands <- list(
        human = import("Data/Genome/Human/UCSC/hg38.cpgIslandExt.fromTable.gtf"), 
        mouse = import("Data/Genome/Mouse/UCSC/mm10.cpgIslandExt.fromTable.gtf"),
        rat = import("Data/Genome/Rat/UCSC/rn6.cpgIslandExt.fromTable.gtf")
    )
    UCSCCpGIslands <- sapply(c("human", "mouse", "rat"), function(species) Genome$standardizeSeqInfo(UCSCCpGIslands[[species]], seqInfo = Genome$MainSeqInfo[[species]], seqLevels = Genome$MainSeqLevels[[species]], prune = TRUE), simplify = FALSE)
    saveRDS(UCSCCpGIslands, file = "Data/GenomicFeatures/UCSCCpGIslands.RDS")
} else {
    UCSCCpGIslands <- readRDS("Data/GenomicFeatures/UCSCCpGIslands.RDS")
}

###########################################################################
## Get a list of EnsIDs (geneID, transcriptID, and intergenicID) with 
## unplaced contigs, chrM, chrY removed
###########################################################################
if (.RERUN) {
    EnsHg38GeneIDsMainNoMY <- names(EnsFeatures$Gene$human)[as.character(seqnames(EnsFeatures$Gene$human)) %in% Genome$MainNoMYSeqLevels$human]
    EnsHg38TranscriptIDsMainNoMY <- names(EnsFeatures$Transcript$human)[as.character(seqnames(EnsFeatures$Transcript$human)) %in% Genome$MainNoMYSeqLevels$human]
    EnsHg38IntergenicIDsMainNoMY <- mcols(EnsFeatures$Intergenic$human)[["name"]][as.character(seqnames(EnsFeatures$Intergenic$human)) %in% Genome$MainNoMYSeqLevels$human]
    EnsMm10GeneIDsMainNoMY <- names(EnsFeatures$Gene$mouse)[as.character(seqnames(EnsFeatures$Gene$mouse)) %in% Genome$MainNoMYSeqLevels$mouse]
    EnsMm10TranscriptIDsMainNoMY <- names(EnsFeatures$Transcript$mouse)[as.character(seqnames(EnsFeatures$Transcript$mouse)) %in% Genome$MainNoMYSeqLevels$mouse]
    EnsMm10IntergenicIDsMainNoMY <- mcols(EnsFeatures$Intergenic$mouse)[["name"]][as.character(seqnames(EnsFeatures$Intergenic$mouse)) %in% Genome$MainNoMYSeqLevels$mouse]
    EnsRn6GeneIDsMainNoMY <- names(EnsFeatures$Gene$rat)[as.character(seqnames(EnsFeatures$Gene$rat)) %in% Genome$MainNoMYSeqLevels$rat]
    EnsRn6TranscriptIDsMainNoMY <- names(EnsFeatures$Transcript$rat)[as.character(seqnames(EnsFeatures$Transcript$rat)) %in% Genome$MainNoMYSeqLevels$rat]
    EnsRn6IntergenicIDsMainNoMY <- mcols(EnsFeatures$Intergenic$rat)[["name"]][as.character(seqnames(EnsFeatures$Intergenic$rat)) %in% Genome$MainNoMYSeqLevels$rat]

    UCSCHg38CpGIDsMainNoMY <- as.character(UCSCCpGIslands$human)[as.character(seqnames(UCSCCpGIslands$human)) %in% Genome$MainNoMYSeqLevels$human]
    UCSCMm10CpGIDsMainNoMY <- as.character(UCSCCpGIslands$mouse)[as.character(seqnames(UCSCCpGIslands$mouse)) %in% Genome$MainNoMYSeqLevels$mouse]
    UCSCRn6CpGIDsMainNoMY <- as.character(UCSCCpGIslands$rat)[as.character(seqnames(UCSCCpGIslands$rat)) %in% Genome$MainNoMYSeqLevels$rat]

    FeatureIDsMainNoMY <- list(
        human = list(Gene = EnsHg38GeneIDsMainNoMY, Transcript = EnsHg38TranscriptIDsMainNoMY, Intergenic = EnsHg38IntergenicIDsMainNoMY, CpG = UCSCHg38CpGIDsMainNoMY),
        mouse = list(Gene = EnsMm10GeneIDsMainNoMY, Transcript = EnsMm10TranscriptIDsMainNoMY, Intergenic = EnsMm10IntergenicIDsMainNoMY, CpG = UCSCMm10CpGIDsMainNoMY),
        rat = list(Gene = EnsRn6GeneIDsMainNoMY, Transcript = EnsRn6TranscriptIDsMainNoMY, Intergenic = EnsRn6IntergenicIDsMainNoMY, CpG = UCSCRn6CpGIDsMainNoMY)
    )
    saveRDS(FeatureIDsMainNoMY, file = "Data/GenomicFeatures/FeatureIDsMainNoMY.RDS")
} else {
    FeatureIDsMainNoMY <- readRDS("Data/GenomicFeatures/FeatureIDsMainNoMY.RDS")
}

###########################################################################
## Get a list of EnsIDs (geneID, transcriptID) with protein coding genes only
###########################################################################
if (.RERUN) {
    EnsHg38GeneIDsProt <- names(EnsFeatures$Gene$human)[grepl("protein_coding", mcols(EnsFeatures$Gene$human)[["gene_biotype"]])]
    EnsHg38TranscriptIDsProt <- names(EnsFeatures$Transcript$human)[grepl("protein_coding", mcols(EnsFeatures$Transcript$human)[["tx_biotype"]])]
    EnsMm10GeneIDsProt <- names(EnsFeatures$Gene$mouse)[grepl("protein_coding", mcols(EnsFeatures$Gene$mouse)[["gene_biotype"]])]
    EnsMm10TranscriptIDsProt <- names(EnsFeatures$Transcript$mouse)[grepl("protein_coding", mcols(EnsFeatures$Transcript$mouse)[["tx_biotype"]])]
    EnsRn6GeneIDsProt <- names(EnsFeatures$Gene$rat)[grepl("protein_coding", mcols(EnsFeatures$Gene$rat)[["gene_biotype"]])]
    EnsRn6TranscriptIDsProt <- names(EnsFeatures$Transcript$rat)[grepl("protein_coding", mcols(EnsFeatures$Transcript$rat)[["tx_biotype"]])]
    EnsFeatureIDsProt <- list(
        human = list(Gene = EnsHg38GeneIDsProt, Transcript = EnsHg38TranscriptIDsProt), 
        mouse = list(Gene = EnsMm10GeneIDsProt, Transcript = EnsMm10TranscriptIDsProt), 
        rat = list(Gene = EnsRn6GeneIDsProt, Transcript = EnsRn6TranscriptIDsProt)
    )
    saveRDS(EnsFeatureIDsProt, file = "Data/GenomicFeatures/EnsFeatureIDsProt.RDS")
} else {
    EnsFeatureIDsProt <- readRDS("Data/GenomicFeatures/EnsFeatureIDsProt.RDS")
}

###########################################################################
## Get a list of EnsIDs (geneID, transcriptID) with mito genes only
###########################################################################
if (.RERUN) {
    EnsHg38GeneIDsChrM <- names(EnsFeatures$Gene$human)[as.character(seqnames(EnsFeatures$Gene$human)) == "chrM" ]
    EnsHg38TranscriptIDsChrM <- names(EnsFeatures$Transcript$human)[as.character(seqnames(EnsFeatures$Transcript$human)) == "chrM"]
    EnsMm10GeneIDsChrM <- names(EnsFeatures$Gene$mouse)[as.character(seqnames(EnsFeatures$Gene$mouse)) == "chrM" ]
    EnsMm10TranscriptIDsChrM <- names(EnsFeatures$Transcript$mouse)[as.character(seqnames(EnsFeatures$Transcript$mouse)) == "chrM"]
    EnsRn6GeneIDsChrM <- names(EnsFeatures$Gene$rat)[as.character(seqnames(EnsFeatures$Gene$rat)) == "chrM" ]
    EnsRn6TranscriptIDsChrM <- names(EnsFeatures$Transcript$rat)[as.character(seqnames(EnsFeatures$Transcript$rat)) == "chrM"]
    EnsFeatureIDsChrM <- list(
        human = list(Gene = EnsHg38GeneIDsChrM, Transcript = EnsHg38TranscriptIDsChrM), 
        mouse = list(Gene = EnsMm10GeneIDsChrM, Transcript = EnsMm10TranscriptIDsChrM), 
        rat = list(Gene = EnsRn6GeneIDsChrM, Transcript = EnsRn6TranscriptIDsChrM)
    )
    saveRDS(EnsFeatureIDsChrM, file = "Data/GenomicFeatures/EnsFeatureIDsChrM.RDS")
} else {
    EnsFeatureIDsChrM <- readRDS("Data/GenomicFeatures/EnsFeatureIDsChrM.RDS")
}

###########################################################################
## Get a list of EnsIDs (geneID, transcriptID) with chrY genes only
###########################################################################
if (.RERUN) {
    EnsHg38GeneIDsChrY <- names(EnsFeatures$Gene$human)[as.character(seqnames(EnsFeatures$Gene$human)) == "chrY" ]
    EnsHg38TranscriptIDsChrY <- names(EnsFeatures$Transcript$human)[as.character(seqnames(EnsFeatures$Transcript$human)) == "chrY"]
    EnsMm10GeneIDsChrY <- names(EnsFeatures$Gene$mouse)[as.character(seqnames(EnsFeatures$Gene$mouse)) == "chrY" ]
    EnsMm10TranscriptIDsChrY <- names(EnsFeatures$Transcript$mouse)[as.character(seqnames(EnsFeatures$Transcript$mouse)) == "chrY"]
    EnsRn6GeneIDsChrY <- names(EnsFeatures$Gene$rat)[as.character(seqnames(EnsFeatures$Gene$rat)) == "chrY" ]
    EnsRn6TranscriptIDsChrY <- names(EnsFeatures$Transcript$rat)[as.character(seqnames(EnsFeatures$Transcript$rat)) == "chrY"]
    EnsFeatureIDsChrY <- list(
        human = list(Gene = EnsHg38GeneIDsChrY, Transcript = EnsHg38TranscriptIDsChrY), 
        mouse = list(Gene = EnsMm10GeneIDsChrY, Transcript = EnsMm10TranscriptIDsChrY), 
        rat = list(Gene = EnsRn6GeneIDsChrY, Transcript = EnsRn6TranscriptIDsChrY)
    )
    saveRDS(EnsFeatureIDsChrY, file = "Data/GenomicFeatures/EnsFeatureIDsChrY.RDS")
} else {
    EnsFeatureIDsChrY <- readRDS("Data/GenomicFeatures/EnsFeatureIDsChrY.RDS")
}

###########################################################################
## More TSS flanking regions
###########################################################################
if (.RERUN) {
    EnsFlanksHg38ByGene <- lapply(c(50, 100, 250, 500, 1000, 1500, 2000, 2500), function(d) Genome$getTSSFlanksFromTxGRs(EnsFeatures[["Transcript"]][["human"]], width = d, both = TRUE))
    names(EnsFlanksHg38ByGene) <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
    EnsFlanksMm10ByGene <- lapply(c(50, 100, 250, 500, 1000, 1500, 2000, 2500), function(d) Genome$getTSSFlanksFromTxGRs(EnsFeatures[["Transcript"]][["mouse"]], width = d, both = TRUE))
    names(EnsFlanksMm10ByGene) <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
    EnsFlanksRn6ByGene <- lapply(c(50, 100, 250, 500, 1000, 1500, 2000, 2500), function(d) Genome$getTSSFlanksFromTxGRs(EnsFeatures[["Transcript"]][["rat"]], width = d, both = TRUE))
    names(EnsFlanksRn6ByGene) <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
    EnsFlanksByGene <- list(human = EnsFlanksHg38ByGene, mouse = EnsFlanksMm10ByGene, rat = EnsFlanksRn6ByGene)
    saveRDS(EnsFlanksByGene, file = "Data/GenomicFeatures/EnsFlanksByGene.RDS")
} else {
    EnsFlanksByGene <- readRDS("Data/GenomicFeatures/EnsFlanksByGene.RDS")
}
