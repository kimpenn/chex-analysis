## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
#######################################################################
sed 1d Data/DNAmoreDB/dnazymes.csv 
seqtk seq -F "I" Data/DNAzyme/dnazyme_long.fa > Data/DNAzyme/dnazyme_long.fq

blastn -query Data/DNAzyme/dnazyme_long.fa -db Data/blast/UCSC/Homo_sapiens/hg38 -num_alignments 1000 -word_size 6 -perc_identity 80 -gapopen 3 -gapextend 1 -evalue 1 -num_threads 12 -out Data/DNAzyme/blast_long_hg38.txt
blastn -query Data/DNAzyme/dnazyme_long.fa -db Data/blast/UCSC/Mus_musculus/mm10 -num_alignments 1000 -word_size 6 -perc_identity 80 -gapopen 3 -gapextend 1 -evalue 1 -num_threads 12 -out Data/DNAzyme/blast_long_mm10.txt
perl Source/functions/blast2sam.pl -s -d Data/DNAzyme/blast_long_hg38.txt > Data/DNAzyme/blast_long_hg38.sam
perl Source/functions/blast2sam.pl -s -d Data/DNAzyme/blast_long_mm10.txt > Data/DNAzyme/blast_long_mm10.sam

perl -e 'while(<>) { chomp; if (/^>(.+)/) { $n=$1 } else { $s=$_; print "$n\t$s\n" } }' Data/DNAzyme/dnazyme_short.fa | parallel --colsep "\t" "seq2profile.pl {2} 0 {1}" > Data/DNAzyme/dnazyme_short.motif

zcat < Data/Genome/UCSC/Mus_musculus/mm10.fa.gz > /tmp/mm10.fa
zcat < Data/Genome/UCSC/Homo_sapiens/hg38.fa.gz > /tmp/hg38.fa

scanMotifGenomeWide.pl Data/DNAzyme/dnazyme_short.motif /tmp/mm10.fa -p 10 -bed | gzip -c > Data/DNAzyme/homer_short_mm10.bed.gz
scanMotifGenomeWide.pl Data/DNAzyme/dnazyme_short.motif /tmp/hg38.fa -p 10 -bed | gzip -c > Data/DNAzyme/homer_short_hg38.bed.gz