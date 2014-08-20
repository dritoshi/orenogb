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
# output.file <- "demo.pdf"
chr         <- args[1]
start.bp    <- eval(parse(text = args[2]))
end.bp      <- eval(parse(text = args[3]))
zoom.power  <- eval(parse(text = args[4]))
# output.file <- args[4]
output.file <- args[5]

cat(paste0(chr, ":", start.bp, "-", end.bp), "\n")
cat("Zoom:",   zoom.power,  "\n")
cat("Output:", output.file, "\n")

# make a GenomicRanges object
range <- GRanges(chr, IRanges(start.bp, end.bp))

# ideogram track
p.ideo <- Ideogram(genome = "mm10", subchr = chr) + xlim(range)

# genes
p.txdb <- autoplot(Mus.musculus, which = range)

# background
bg   <- BSgenome.Mmusculus.UCSC.mm10
p.bg <- autoplot(bg, which = range)

# draw my track
tks <- tracks(p.ideo, ref = p.bg, gene = p.txdb, heights = c(2, 1, 4)) + xlim(range)
tks <- tks + zoom(zoom.power)
# tks <- tks + theme_tracks_sunset(bg = "#DFDFDF")

# output
pdf(output.file)
print(tks)
dev.off()
