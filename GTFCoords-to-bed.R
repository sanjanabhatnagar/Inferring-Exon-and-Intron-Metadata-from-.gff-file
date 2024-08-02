library("bedr")
library("rtracklayer")
suppressPackageStartupMessages(library('GenomicFeatures'))
args <- commandArgs(trailingOnly = TRUE)
intron_coord <- args[1]

Intron_processed <- read.table(intron_coord, sep="\t")


Intron_processed <- Intron_processed[-1, ]
Intron_processed$"V4" <- as.integer(Intron_processed$"V4")
Intron_processed$"V3" <- as.integer(Intron_processed$"V3")

Intron_processed_pv <- Intron_processed[Intron_processed[, "V5"] == "+", ]

print(head(Intron_processed_pv, 5))
Intron_processed_pv <- bedr.sort.region(Intron_processed_pv[c("V2","V3","V4")], check.chr = FALSE)
Intron_processed_pv.bed <- convert2bed(Intron_processed_pv[c("V2","V3","V4")], check.zero.based = FALSE,
                                       check.chr = FALSE, verbose = TRUE )

write.table(Intron_processed_pv.bed, file='Introns_positiveStrn.bed', col.names = TRUE, row.names = FALSE)

Intron_processed_nv <- Intron_processed[Intron_processed[, "V5"] == "-", ]

print(head(Intron_processed_nv, 5))
Intron_processed_nv <- bedr.sort.region(Intron_processed_nv[c("V2","V3","V4")], check.chr = FALSE)
Intron_processed_nv.bed <- convert2bed(Intron_processed_nv[c("V2","V3","V4")], check.zero.based = FALSE,
                                       check.chr = FALSE, verbose = TRUE )

write.table(Intron_processed_nv.bed, file='Introns_negativeStrn.bed', col.names = TRUE, row.names = FALSE)