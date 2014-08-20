#
# orenogb.R
#
library("ggbio")
library("GenomicRanges")
library("GenomicAlignments")
# library("ShortRead")
# library("parallel")

library("Mus.musculus")
library("BSgenome.Mmusculus.UCSC.mm10")

#
# Functions
#
getGRangesfromBam <- function(bamFile, wh, ...) {
  param <- ScanBamParam(
  	what  = c("pos", "qwidth"),
    which = wh,
    flag  = scanBamFlag(isUnmappedQuery = FALSE)
  )
  scanBam(bamFile, ..., param = param)[[1]]
} 

#
# main
#

# range
args <- commandArgs(TRUE)

if ( is.na(args[1]) ) {
	rm(list = ls())
  chr         <- "chr17"
  start.bp    <- 35502880
  end.bp      <- 35516079
  zoom.power  <- 1
  bam.files    <- "~/Sources/Quartz_01.th.rmrRNA.bam,~/Sources/Quartz_02.th.rmrRNA.bam"
  output.file <- "demo.pdf"
} else {  
  chr         <- args[1]
  start.bp    <- eval(parse(text = args[2]))
  end.bp      <- eval(parse(text = args[3]))
  zoom.power  <- eval(parse(text = args[4]))
  bam.files   <- args[5]
  output.file <- args[6]
}  

cat(paste0(chr, ":", start.bp, "-", end.bp), "\n")
cat("Zoom:",   zoom.power,  "\n")
cat("bam files:", bam.files, "\n")
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

# bam
bam.files  <- unlist(strsplit(bam.files, ","))
cat("bam files:", bam.files, "\n")
bam.views  <- BamViews(bam.files, bamRanges = range)
bam.galign <- readGAlignmentsFromBam(bam.views)

# if (0) {
# draw my track
tks <- tracks(
  p.ideo,
  ref      = p.bg,  
  gene     = p.txdb,  
  heights  = c(1, 1, 4)
)
for (i in seq_along(bam.files)) {
	p.mis <- autoplot(bam.files[i], bsgenome = bg, which = range, stat = "mismatch")
  tks <- tks + tracks(p.mis, heights = 2)
}  
tks <- tks + xlim(range)
tks <- tks + ggbio:::zoom(zoom.power)
# tks <- tks + theme_tracks_sunset(bg = "#DFDFDF")

# output
pdf(output.file)
print(tks)
dev.off()
# }