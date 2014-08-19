#
# orenogb.R
#
library("ggbio")
library("GenomicRanges")
library("Mus.musculus")
library("BSgenome.Mmusculus.UCSC.mm10")

# range
args <- commandArgs(TRUE)
# chr <- "chr2"
# start.bp <- 1e8
# end.bp   <- 1e8+3000000
chr         <- args[1]
start.bp    <- as.numeric(args[2])
end.bp      <- as.numeric(args[3])
output.file <- args[4]
cat(paste(chr, start.bp, end.bp), "\n")

range <- GRanges(chr, IRanges(start.bp, end.bp))

# An ideogram track
p.ideo <- Ideogram(genome = "mm10", subchr = chr) + xlim(range)

# genes
p.txdb <- autoplot(Mus.musculus, which = range)

# background
bg <- BSgenome.Mmusculus.UCSC.mm10
p.bg <- autoplot(bg, which = range)

# draw my track
tks <- tracks(p.ideo, ref = p.bg, gene = p.txdb, heights = c(2, 1, 4)) + xlim(range)
# tks <- tks + theme_tracks_sunset(bg = "#DFDFDF")

# output
pdf(output.file)
print(tks)
dev.off
