#
# orenogb.R
#
library("ggbio")
library("GenomicRanges")
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
  chr         <- "chr17"
  start.bp    <- 35502880
  end.bp      <- 35516079
  zoom.power  <- 1
  bam.file    <- "~/Sources/Quartz_01.th.rmrRNA.bam" 
  output.file <- "demo.pdf"
} else {  
  chr         <- args[1]
  start.bp    <- eval(parse(text = args[2]))
  end.bp      <- eval(parse(text = args[3]))
  zoom.power  <- eval(parse(text = args[4]))
  bam.file    <- args[5]
  output.file <- args[6]
}  

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

# bam
p.mis <- autoplot(bam.file, bsgenome = bg, which = range, stat = "mismatch")

# draw my track
tks <- tracks(
  p.ideo,
  coverage = p.mis,
  ref      = p.bg,
  gene     = p.txdb,
  heights  = c(2, 3, 1, 4)
)
tks <- tks + xlim(range)
tks <- tks + zoom(zoom.power)
# tks <- tks + theme_tracks_sunset(bg = "#DFDFDF")

# output
pdf(output.file)
print(tks)
dev.off()
