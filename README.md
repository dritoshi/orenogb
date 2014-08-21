Orenogb
====

Visualization command for genomic data

## Description
- just type simple command on your terminal
- beautiful visualization by ggbio and ggplot2
- semantic zoom
- eaily loading of transcriptome annotation
- search by gene Symbol
- choose species

## Demo

    $ R --slave --vanilla -f orenogb.R --args coordination mm10 chr17 35400000 35600000 1 Quartz_01.th.rmrRNA.bam,Quartz_02.th.rmrRNA.bam demo/demo.pdf

![demo](demo/demo.png)

### Exponential notation

    $ R --slave --vanilla -f orenogb.R --args coordination mm10 chr17 3.55e7+2880 3.55e7+16079 1 Quartz_01.th.rmrRNA.bam,Quartz_02.th.rmrRNA.bam demo/demo2.pdf

![demo](demo/demo2.png)

### Semantic Zoom

    $ R --slave --vanilla -f orenogb.R --args coordination mm10 chr17 35502880 35516079 1/200 Quartz_01.th.rmrRNA.bam,Quartz_02.th.rmrRNA.bam demo/demo3.pdf

![demo](demo/demo3.png)    

### Search by Gene Symbol

    $ R --slave --vanilla -f orenogb.R --args gene mm10 Pou5f1 1 Quartz_01.th.rmrRNA.bam,Quartz_02.th.rmrRNA.bam demo/demo4.pdf

![demo](demo/demo4.png)

### Human Genome

    $ R --slave --vanilla -f orenogb.R --args gene hg19 POU5F1 1 Quartz_01.th.rmrRNA.bam,Quartz_02.th.rmrRNA.bam demo/demo4.pdf

![demo](demo/demo5.png)

## Requirement
- R
- Bioconductor Software Packages
    - ggbio
    - GenomicRanges
    - GenomicAlignments
- Bioconductor Annotation Packages
    - Mus.musculus, Homo.sapiens, ...
    - BSgenome.Mmusculus.UCSC.mm10, BSgenome.Hsapiens.UCSC.hg19, ...

## Usage

    $ R --slave --vanilla -f orenogb.R --args [mode] [genome] [chr] [start bp] [end bp] [zoom] [bam1,bam2,...] [output file]

## Install

    $ git clone git@github.com:dritoshi/orenogb.git
    $ cd orenogb
    $ sudo R
    R> source("http://bioconductor.org/biocLite.R")
    R> biocLite(c("ggbio", "GenomicRanges", "GenomicAlignments")
    R> biocLite(c("Mus.musculus", "BSgenome.Mmusculus.UCSC.mm10"))
    R> biocLite(c("Homo.sapiens", "BSgenome.Hsapiens.UCSC.hg19"))    

## ToDo
- choose species
- wapper by shell script
- unit test
- deamonaize R script
- draw sample name

## Contribution

## Licence

[MIT](https://github.com/dritoshi/orenogb/blob/master/LICENCE)

## Author

[dritoshi](https://github.com/dritoshi)
