#!/bin/Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("qqman"))

option_list <- list(
        make_option(c('-r', '--regions'), action='store', type='character', default='regions.bed', help='bed file of regions of interest'),
        make_option(c('-d', '--dnv'), action='store', type='character', default='dnvs.bed', help='de novo variants in a bed file'),
        make_option(c('-p', '--prior'), action='store', type='character', default='file.txt', help='Prior file generated from genomic_mutation_rate.py'),
        make_option(c('-n', '--nIndividuals'), action='store', type='numeric', default='1', help='number of individuals'),
        make_option(c('-w', '--weightCADD'), action='store', type='character', default='no', help='use CADD scores as weights (by default this is no and only the baseline mutation rate is used in the p-value calculation)'),
        make_option(c('-m', '--manhattan'), action='store', type='character', default='manhattan.pdf', help='manhattan plot pdf name'),
        make_option(c('-o', '--output'), action='store', type='character', default='output.txt', help='Output text file name')
)
opt <- parse_args(OptionParser(option_list = option_list))

regions = opt$regions
dnv = opt$dnv
priors = opt$prior
individuals = opt$nIndividuals
use_weights = opt$weightCADD
manhattan_plot = opt$manhattan
output = opt$output

regionFile <- import(regions, format = "BED")  
dnvFile <- import(dnv, format = "BED")

overlaps <- findOverlaps(regionFile, dnvFile)
counts <- countOverlaps(regionFile, dnvFile)

result <- data.frame(
  chr = seqnames(regionFile),
  start = start(regionFile),
  end = end(regionFile),
  count = counts,
  region = mcols(regionFile)$name
)

result <- result[,c("region", "count")]

if (grepl("\\.gz$", priors)) {
  probs <- read.delim(gzfile(priors))
} else {
  probs <- read.delim(priors)
}
colnames(probs)[4] <- "region"

m <- merge(probs, result, by="region")
m <- m[which(m$count > 0),]
m <- m[!is.na(m$constrained_mutation_rate),]

for(i in 1:nrow(m)){
  if (use_weights == "yes") {
    m$pValue[i] <- binom.test(x = m$count[i], n = individuals * 70, p = m$constrained_mutation_rate[i])$p.value
  } else {
    m$pValue[i] <- binom.test(x = m$count[i], n = individuals * 70, p = m$mutation_rate[i])$p.value
  }
}

m$chr <- gsub("chr", "", m$chr)
m$chr <- gsub("X", "23", m$chr)
m$chr <- gsub("Y", "24", m$chr)
m$chr <- as.numeric(m$chr)

m <- m[which(m$pValue != 0),]

pdf(manhattan_plot, width=7, height=6)
manhattan(m, chr="chr", bp="start", p="pValue", snp="region", suggestiveline = F, genomewideline = F)
abline(h=-log10(5e-8), col="blue")
abline(h=-log10(8.333333e-12), col="red")
qq(m$pValue)
dev.off()

write.table(m, file=output, sep="\t", quote=F, row.names=F)