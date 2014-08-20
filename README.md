Orenogb
====

Visualization command for genomic data

## Description


## Demo

    $ R --slave --vanilla -f orenogb.R --args chr17 35400000 35600000 1 ~/Sources/Quartz_01.th.rmrRNA.bam demo.pdf

![demo](demo.png)

### Exponential notation

    $ R --slave --vanilla -f orenogb.R --args chr17 3.55e7+2880 3.55e7+16079 1 ~/Sources/Quartz_01.th.rmrRNA.bam demo2.pdf

![demo](demo2.png)

### Semantic Zoom

    $ R --slave --vanilla -f orenogb.R --args chr17 35502880 35516079 1/200 ~/Sources/Quartz_01.th.rmrRNA.bam demo3.pdf

![demo](demo3.png)

## Requirement
- R
- Bioconductor Software Packages
    - ggbio
    - GenomicRanges
- Bioconductor Annotation Packages
    - Mus.musculus
    - BSgenome.Mmusculus.UCSC.mm10

## Usage

    $ R --slave --vanilla -f orenogb.R --args [chr] [start bp] [end bp] [zoom] [bam] [output file]

## Install

    $ git clone git@github.com:dritoshi/orenogb.git
    $ cd orenogb
    $ sudo R
    R> source("http://bioconductor.org/biocLite.R")
    R> biocLite(c("ggbio", "GenomicRanges")
    R> biocLite(c("Mus.musculus", "BSgenome.Mmusculus.UCSC.mm10"))

## ToDo
- load multiple data files
- choose species
- wapper by shell script
- unit test

## Contribution

## Licence

[MIT](https://github.com/dritoshi/orenogb/blob/master/LICENCE)

## Author

[dritoshi](https://github.com/dritoshi)
