##
## orenogb.R
##
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
# Input
#
args <- commandArgs(TRUE)

if ( is.na(args[1]) ) {

	rm(list = ls())
  search.mode <- "gene"
  # search.mode <- "coordination"
  chr         <- "chr17"
  start.bp    <- 35502880
  end.bp      <- 35516079
  zoom.power  <- 1
  output.file <- "demo.pdf"
  bam.files   <- "~/Dropbox_RIKEN/Public_ACCCBiT/Data/Quartz-Seq/bam/Quartz_01.th.rmrRNA.bam,~/Dropbox_RIKEN/Public_ACCCBiT/Data/Quartz-Seq/bam/Quartz_02.th.rmrRNA.bam"
  gene.name   <- "Pou5f1"

} else {  
	search.mode = args[1]
	if (search.mode == "coordination") {
    chr         <- args[2]
    start.bp    <- eval(parse(text = args[3]))
    end.bp      <- eval(parse(text = args[4]))
    zoom.power  <- eval(parse(text = args[5]))
    bam.files   <- args[6]
    output.file <- args[7]     
  } else if (search.mode == "gene") {
    gene.name = args[2]
    zoom.power  <- eval(parse(text = args[3]))
    bam.files   <- args[4]
    output.file <- args[5]    
  } else {
  	cat("Usage: ...", "\n")
  }     
}

bam.files  <- unlist(strsplit(bam.files, ","))

if (search.mode == 'gene') {
	library(biomaRt)

  ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  gene <- getBM(
       attributes = c('chromosome_name', 'start_position', 'end_position'),
       filters = 'mgi_symbol',
       values  = gene.name,
       mart    = ensembl
  )
  chr      <- paste0('chr', as.character(gene[,1]))
  start.bp <- gene[,2]
  end.bp   <- gene[,3]
}

cat(paste0(chr, ":", start.bp, "-", end.bp), "\n")
cat("Zoom:",   zoom.power,  "\n")
cat( paste(c("bam files:", bam.files), "\n") )
cat("Output:", output.file, "\n")

#
# Plot
#

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