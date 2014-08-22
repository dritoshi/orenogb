##
## orenogb.R
##
library("ggbio")
library("GenomicAlignments")

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
getPositionBySymbol <- function(db, gene.name) {
  
  if (db == 'Mus.musculus') {
    # SYMBOL, CHR, TXSTART, TXEND
    idx <- c(18, 11, 49, 50)
    symbol.idx <- 18
  } else if (db == 'Homo.sapiens') {
    idx <- c(19, 11, 51, 52)
    symbol.idx <- 19
  } 

  cls <- columns(eval(parse(text = db)))[idx]
  kt  <- keytypes(eval(parse(text = db)))[symbol.idx]  # SYMBOL

  res <- select(eval(parse(text = db)), keys = gene.name, columns = cls, keytype = kt)
  res <- data.frame(
    SYMBOL  = res[1,]$SYMBOL[1],
    CHR     = paste0("chr", res[1,]$CHR[1]),
    TXSTART = min(res[1,]$TXSTART),
    TXEND   = max(res[1,]$TXEND)
  )
  return(res)
}
chooseGenome <- function(genome.version) {
  if (genome.version == 'mm10') {
      return(c('Mus.musculus', 'BSgenome.Mmusculus.UCSC.mm10'))
  } else if (genome.version == 'hg19') {
      return(c('Homo.sapiens', 'BSgenome.Hsapiens.UCSC.hg19'))
  }
}

#
# Input
#
args <- commandArgs(TRUE)

if ( is.na(args[1]) ) {
  search.mode <- "gene"
  # search.mode <- "coordination"  

  genome.ver  <- 'mm10'    
  chr         <- "chr17"
  start.bp    <- 35502880
  end.bp      <- 35516079
  zoom.power  <- 1
  output.file <- "demo.pdf"
  bam.files   <- "~/Dropbox_RIKEN/Public_ACCCBiT/Data/Quartz-Seq/bam/Quartz_01.bam,~/Dropbox_RIKEN/Public_ACCCBiT/Data/Quartz-Seq/bam/Quartz_02.bam"
  gene.name   <- "Pou5f1"

} else {  
	search.mode = args[1]

	if (search.mode == "coordination") {   
    genome.ver  <- args[2]
    chr         <- args[3]
    start.bp    <- eval(parse(text = args[4]))
    end.bp      <- eval(parse(text = args[5]))
    zoom.power  <- eval(parse(text = args[6]))
    bam.files   <- args[7]
    output.file <- args[8]
  } else if (search.mode == "gene") {
    genome.ver  <- args[2]
    gene.name   <- args[3]
    zoom.power  <- eval(parse(text = args[4]))
    bam.files   <- args[5]
    output.file <- args[6]
  } else {
  	cat("Usage: ...", "\n")
  }     
}

bam.files  <- unlist(strsplit(bam.files, ","))

# choose genome and load libraries
genome <- chooseGenome(genome.ver)
orgdb    <- genome[1]
bsgenome <- genome[2]
library(orgdb,    character.only = TRUE)
library(bsgenome, character.only = TRUE)

# get genomic coordination of gene
if (search.mode == 'gene') {
  gene <- getPositionBySymbol(orgdb, gene.name)
  chr      <- as.character(gene$CHR)
  start.bp <- gene$TXSTART
  end.bp   <- gene$TXEND
}

# message
cat("### Mode:",      search.mode, "\n")
cat("### Genome:",    genome.ver,  "\n")
cat("### Position:",  paste0(chr, ":", start.bp, "-", end.bp), "\n")
cat("### Zoom:",      zoom.power,  "\n")
cat("### bam files:", bam.files,   "\n")
cat("### Output:",    output.file, "\n")

#
# Plot
#
cat("### Start plot...\n")
# make a GenomicRanges object
range <- GRanges(chr, IRanges(start.bp, end.bp))

# ideogram track
cat("### build ideogram...\n")
p.ideo <- Ideogram(genome = genome.ver, subchr = chr) + xlim(range)

# genes
cat("### build genes...\n")
p.txdb <- autoplot(eval(parse(text = orgdb)), which = range)

# background
# bg   <- BSgenome.Mmusculus.UCSC.mm10
cat("### build background...\n")
bg <- eval(parse(text = bsgenome))
p.bg <- autoplot(bg, which = range)

# bam
cat("### build reads from bam data...\n")
bam.views  <- BamViews(bam.files, bamRanges = range)
bam.galign <- readGAlignmentsFromBam(bam.views)

# draw my track
cat("### draw all tracks...\n")
tks <- tracks(
  p.ideo,
  ref      = p.bg,  
  gene     = p.txdb,  
  heights  = c(1, 1, 4)
)
for (i in seq_along(bam.files)) {
  sample.name = sub('.bam', '', basename(bam.files[i]))
	p.mis <- autoplot(bam.files[i], bsgenome = bg, which = range, stat = "mismatch")

  sample.name <- sub('-', '', sample.name)
  tracks.str <- paste0('tracks(', sample.name, ' = p.mis, heights = 2)')
  cat('### ', tracks.str, "\n")
  tks <- tks + eval(parse(text = tracks.str))
  #tks <- tks + tracks(hoge = p.mis, heights = 2)
 
}
tks <- tks + xlim(range)
tks <- tks + ggbio:::zoom(zoom.power)

# output
pdf(output.file)
print(tks)
dev.off()
