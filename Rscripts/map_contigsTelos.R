#!/usr/bin/env Rscript

################################################################################
################################################################################
# -- ad hoc functions
find_kmerPos <- function(dnass, kmers, jumpSize, minWidth, kmerProp){

  # -- for each chromosome, find and aggregate all kmers
  acrossChrs <- rbindlist(lapply(names(dnass), function(i){
    mlist <- lapply(kmers, function(x) vmatchPattern(x, subject = dnass[i])[[1]])
    mergem <- reduce(do.call(c, mlist), min.gapwidth = 0, drop.empty.ranges = T)
    if(length(mergem) == 0)
      return(NULL)
    mergem <- reduce(
      mergem,
      min.gapwidth = jumpSize,
      drop.empty.ranges = TRUE,
      with.revmap = TRUE)
    if(length(mergem) == 0)
      return(NULL)
    out <- as.data.table(mergem)
    out[,nbp := sapply(revmap, length) * max(nchar(kmers))]
    out[,`:=`(propSeq = nbp / width, chr = i)]
    out <- subset(out, width >= minWidth & propSeq >= kmerProp)
    out <- out[,c("chr", "start", "end")]
    return(out)
  }))

  if(!is.null(acrossChrs)){
    if(nrow(acrossChrs) < 1){
      acrossChrs <- NULL
    }
  }

  return(acrossChrs)
}

pull_chrLens <- function(dnass){
  data.table(
    chr = names(dnass),
    chrStart = 1,
    chrEnd = width(dnass))
}

add_nogapBlocks <- function(chrLens, kmerBed){
  # -- combine kmer-blocks with chromosome lengths
  tp <- merge(chrLens, kmerBed, by = "chr")
  setkey(tp, chr, start, end)
  tp[,index := 1:.N, by = "chr"]

  # -- get lagged values for the left block
  tp[,`:=`(startLeft = c(1, end[-.N]+1),
           endLeft = c(start-1)),
     by = "chr"]
  tpo <- with(tp, data.table(
    chr = c(chr, chr),
    type = rep(c(tp$type, rep("noGap", nrow(tp)))),
    start = c(start, startLeft),
    end = c(end, endLeft)))

  tpo <- subset(tpo, complete.cases(tpo))
  tpe <- tp[,list(start = ifelse(all(is.na(end)), 1L, as.integer(max(end, na.rm = T)+1)),
                  end = chrEnd[1]), by = "chr"]
  tpe[,type := "noGap"]

  tp <- rbind(tpo, tpe)
  setkey(tp, chr, start, end)
  tp[,width := end - start]
  setcolorder(tp, colnames(kmerBed))

  tmp <- subset(chrLens, !chr %in% tp$chr)
  if(nrow(tmp > 0)){
    chrMiss <- with(tmp, data.table(
      chr = chr, start = chrStart, end = chrEnd, width = chrEnd, type = "noGap"))
    tp <- rbind(tp, chrMiss)
  }

  return(tp)
}

get_teloPos <- function(bed, chrLens, teloBuff){
  bed <- merge(chrLens, bed, by = "chr", all.x = T)
  bed[,pos :=
        ifelse(is.na(start), "none",
               ifelse(start < teloBuff, "left",
                      ifelse(end > (chrEnd - teloBuff), "right", "mid")))]
  bedMd <- bed[,list(teloPosition = paste(unique(pos), collapse = ",")),
               by = c("chr", "chrEnd")]
  setnames(bedMd, "chrEnd", "width")
  return(bedMd)
}

rename_chrs <- function(ss){
  ns <- as.data.table(tstrsplit(trimws(names(ss)), " "))
  hasna <- sapply(ns, function(x) any(is.na(x)))
  if(all(hasna)){
    return(ss)
  }

  ns <- ns[,!hasna, with = F]
  if(ncol(ns) == 1){
    # -- if one column (no whitespace) return names
    chrns <- ns[[1]]
  }else{
    # -- if the first column is all unique and the first entries are "chrx" or numeric, use it
    u1 <- ns[[1]]
    if(!any(duplicated(u1)) && gsub("chr|chromosome|scaffold|chr0","",tolower(u1)[[1]]) == "1"){
      chrns <- ns[[1]]
    }else{
      # -- if there is a column with "chr" in it and it is all unique use that
      wh <- grepl("chr", tolower(unlist(ns[1,])))
      us <- apply(ns, 2, function(x) all(!duplicated(x)))
      if(sum(wh & us) == 1){
        chrns <- as.character(ns[[wh & us]])
      }else{
        # -- if there is a column with "chr", followed by one totally unique, use that
        if(sum(wh) == 1){
          index <- which(wh) + 1
          if(all(!duplicated(ns[[index]]))){
            chrns <- ns[[index]]
          }else{
            # -- otherwise just return the names
            chrns <- names(ss)
          }
        }else{
          chrns <- names(ss)
        }
      }
    }
  }
  names(ss) <- chrns
  return(ss)
}

classify_engine <- function(faFile,
                            teloKmers,
                            gapSize,
                            teloGapSize,
                            minTeloSize,
                            teloRepDens,
                            verbose,
                            minChrSize){

  if(!requireNamespace("Biostrings", quietly = TRUE))
    stop("you need to install Biostrings from Bioconductor before running\n")

  if(!requireNamespace("data.table", quietly = TRUE))
    stop("you need to install data.table before running\n")

  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))

  # -- read in the fasta as a DNAStringSet
  if(verbose)
    cat(sprintf("Reading: %s\n", faFile))
  ss <- readDNAStringSet(faFile)
  ss <- ss[width(ss) > minChrSize]
  ss <- rename_chrs(ss)

  chrLens <- pull_chrLens(ss)
  chrLens <- subset(chrLens, complete.cases(chrLens))
  chrLens[,`:=`(v1 = frank(chr), v2 = frank(as.numeric(gsub("[^0-9.-]", "", chr))))]
  setkey(chrLens, v2, v1, chr)
  chrLens[,`:=`(v1 = NULL, v2 = NULL)]
  chrv <- chrLens$chrEnd; names(chrv) <- chrLens$chr

  if(verbose)
    cat(strwrap(sprintf(
      "\tN. scaffolds = %s (N50 = %sMb); total length = %sMb\n",
      nrow(chrLens),
      round(N50(chrLens$chrEnd / 1e6), 2),
      round(sum(chrLens$chrEnd / 1e6, na.rm = T), 2)),
      indent = 8, exdent = 16), sep = "\n")

  # -- make a bed with the start and end positions of the gap kner (N)
  # then cluster into gap coordinates
  gapNs <- find_kmerPos(
    dnass = ss,
    kmers = "N",
    jumpSize = 1,
    minWidth = gapSize,
    kmerProp = .99)
  if(is.null(gapNs)){
    contigCoords <- with(chrLens, data.table(
      chr = chr, start = chrStart, end = chrEnd, type = "noGap"))
    noGap <- data.table(contigCoords)
  }else{
    gapNs[,type := "gap"]

    # -- pull ungapped (contig) coordinates
    seqi <- Seqinfo(chrLens$chr, seqlengths = chrLens$chrEnd)
    gapr <- with(gapNs, GRanges(
      chr, IRanges(start = start, end = end), seqinfo = seqi))

    noGap <- gaps(gapr)
    noGap <- subset(noGap, strand == "*")
    noGap <- with(as.data.table(noGap), data.table(
      chr = seqnames, start = start, end = end, type = "noGap"))

    contigCoords <- rbind(noGap, gapNs)
  }

  setkey(contigCoords, chr, start, end)
  cat(strwrap(sprintf(
    "\tN. contigs = %s (N50 = %sMb); total length = %sMb\n",
    nrow(noGap),
    round(N50((noGap$end - noGap$start) / 1e6), 2),
    round(sum((noGap$end - noGap$start) / 1e6, na.rm = T), 2)),
    indent = 8, exdent = 16), sep = "\n")

  # -- make a bed with the start and end positions of the telomere kmers
  # then cluster into gap coordinates
  teloBuff <- 1e5
  teloKmers <- unique(
    c(teloKmers, as.character(reverseComplement(DNAStringSet(teloKmers)))))
  telos <- find_kmerPos(
    dnass = ss,
    kmers = teloKmers,
    jumpSize = teloGapSize,
    kmerProp = teloRepDens,
    minWidth = minTeloSize)
  if(!is.null(telos)){
    telos[,type := "telo"]
    telos[,chrLen := chrv[chr]]
    telos[,pos := ifelse(start <= teloBuff, "L", ifelse((chrLen - end) < teloBuff, "R", "M"))]

    telopos <- telos[,list(pos = paste(unique(pos), collapse = "")), by = "chr"]
    cat(strwrap(sprintf(
      "\tN. telomeres = %s (%skb) on chrs: %s\n",
      nrow(telos),
      round(N50((telos$end - telos$start) / 1e3), 2),
      paste(sprintf("%s%s", telopos$chr, telopos$pos), collapse = ", ")),
      indent = 8, exdent = 16), sep = "\n")
    telos[,`:=`(chrLen = NULL, pos = NULL)]
    out <- rbind(telos, contigCoords)
  }else{
    out <- data.table(contigCoords)
    cat("\tCould not find any telomeres\n")
  }

  return(out)
}

plot_classes <- function(blockCalls, palette){

  mxn <- with(blockCalls, max(tapply(
    type, chr, FUN = function(x) sum(x == "noGap"))))
  nColors <- floor(ifelse(mxn/3 < 5, 5, mxn / 3))
  if(nColors > 50)
    nColors <- 50
  cols <- palette(nColors)

  tp <- data.table(blockCalls)
  tp[,genome := factor(genome, levels = unique(genome))]
  if(!"genome" %in% colnames(tp))
    tp[,genome := "genome"]
  tpb <- subset(tp, type == "noGap")
  tel <- subset(tp, type == "telo")
  tel[,x := ((end + start)/2)/1e6]

  tpb[,index := rep(1:nColors, nrow(tp))[1:.N], by = c("genome", "chr")]
  tpb[,y1 := as.numeric(factor(paste(genome, chr), levels = unique(paste(genome, chr))))]
  tpb[,y2 := y1 - .6]
  tpb[,y := (y1 + y2)/2]
  tpb[,genome := factor(genome, levels = unique(genome))]
  ng <- tpb[,list(n = .N, x = max(end)/1e6), by = c("genome", "chr", "y")]
  tmp <- subset(tpb, !duplicated(paste(chr, genome)))
  yv <- tmp$y2; names(yv) <- with(tmp, paste(chr, genome))
  tel[,y2 := yv[paste(chr, genome)]]
  tpl <- subset(tpb, !duplicated(paste(genome, chr)))
  p <- ggplot()+
    geom_rect(data = tpb,
              aes(xmin = start / 1e6, xmax = end / 1e6,
                  ymin = y2, ymax = y1,
                  fill = factor(index)))+
    scale_y_reverse(
      expand = c(0.01, 0.01),
      breaks = tpl$y,
      labels = tpl$chr,
      name = "Chromosome")+
    scale_x_continuous(
      name = "Physical position (Mb); * indicate telomere sequence")+
    geom_text(data = ng,
              aes(x = x, y = y, label = n),
              hjust = -.25, size = 2)+
    geom_point(data = tel, aes(x = x, y = y2), shape = "*", col = "red", size = 4)+

    scale_fill_manual(values = cols, guide = "none")+
    theme(panel.background = element_rect(fill = "darkgrey"),
          panel.grid = element_blank())+
    facet_wrap(~genome, scale = "free")+
    ggtitle(sprintf(
      "Contig positions: each cycle through colors is %s contigs (%s gaps)",
      nColors, nColors - 1))
  print(p)
}

classify_genomes <- function(faFiles,
                             toStrip = ".fa$|.fa.gz$",
                             teloKmers = c("CCCGAAA", "CCCTAAA"),
                             gapSize = 1e4,
                             teloGapSize = 1e4,
                             minTeloSize = 1e3,
                             teloRepDens = 0.1,
                             binScaffolds = TRUE,
                             minChrSize = 1e6,
                             palette = viridisLite::viridis,
                             verbose = TRUE){

  blockCoords <- rbindlist(lapply(faFiles, function(i){
    bc <- classify_engine(
      faFile = i,
      teloKmers = teloKmers,
      gapSize = gapSize,
      teloGapSize = teloGapSize,
      minTeloSize = minTeloSize,
      teloRepDens = teloRepDens,
      minChrSize = minChrSize,
      verbose = verbose)
    bc[,genome := gsub(toStrip, "", basename(i))]
    return(bc)
  }))

  # -- make the plot
  blockCoords[,genome := factor(genome, levels = gsub(toStrip, "", basename(faFiles)))]
  blockCoords[,chr := gsub("chr|chromosome|scaffold|contig", "", tolower(chr))]
  blockCoords[,`:=`(chri = as.numeric(gsub('\\D+','', chr)),
                    chrn = frank(chr, ties.method = "dense")), by = "genome"]
  setorder(blockCoords, genome, chri, chrn, na.last = T)

  plot_classes(
    blockCalls = data.table(blockCoords),
    palette = palette)

  return(blockCoords)
}
################################################################################

if (!require('argparse', quietly = T))
  install.packages('argparse')
suppressPackageStartupMessages(library("argparse"))

################################################################################
################################################################################
# --create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument(
  "-f", "--faFiles", type = "character", default = "",
  help = "comma-separated string of paths to fasta files")
parser$add_argument(
  "-t", "--telokmers", type = "character", default = "CCCGAAA,CCCTAAA",
  help = "comma-separated string of paths to fasta files")
parser$add_argument(
  "-z", "--teloclust", type = "character", default = "50,500,0.75",
  help = "comma-separated cluster params (max gap, min width, min %telo sequence) [default = %(default)s]")
parser$add_argument(
  "-p", "--pdfFile", type = "character", default = "contigPlot.pdf",
  help = "output pdf file name [default= %(default)s]")
parser$add_argument(
  "-o", "--outFile", type = "character", default = "contigCoords.txt",
  help = "output text file with contig coordinates [default= %(default)s]")
parser$add_argument(
  "-m", "--minChrLen", type = "integer", default = 1e6,
  help = "the smallest chromosome to show [default= %(default)s]")
parser$add_argument(
  "-g", "--gapSize", type = "integer", default = 1e3,
  help = "window size to average over to merge gaps [default= %(default)s]")
parser$add_argument(
  "-r", "--regex", type = "character", default = ".fa.gz$|.fa$",
  help = "quoted character, a regex to strip from the file name [default= %(default)s]")
parser$add_argument(
  "-q", "--quietly", action = "store_false",
  dest = "verbose", help = "Print little output")
parser$add_argument(
  "-c", "--print2stdout", action = "store_true",
  dest = "print2stdout", help = "print the summary file to stdout instead of outfile")


################################################################################
################################################################################
# -- get dependencies in order
if(!requireNamespace("Biostrings", quietly = TRUE))
  stop("you need to install Biostrings from Bioconductor before running\n")
if(!requireNamespace("data.table", quietly = TRUE))
  stop("you need to install data.table before running\n")
if(!requireNamespace("ggplot2", quietly = TRUE))
  stop("you need to install ggplot2 before running\n")
if(!requireNamespace("viridisLite", quietly = TRUE))
  stop("you need to install viridisLite before running\n")
if(!requireNamespace("GenomicRanges", quietly = TRUE))
  stop("you need to install GenomicRanges before running\n")
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicRanges))
################################################################################
################################################################################
# -- get the fafiles organized
args <- parser$parse_args()
faf <- as.character(args$faFiles)
fs <- NULL

if(faf == "")
  faf <- getwd()

if(dir.exists(faf)){
  fs <- list.files(path = faf, pattern = ".fa$|.fa.gz$", full.names = T)
}
if(length(fs) == 0){
  fs <- strsplit(faf, ",")[[1]]
  fs <- fs[file.exists(fs)]
}

if(length(fs) == 0)
  stop("could not find faFiles specified\n")

faFilesIn <- fs

# -- get arguments in order
toStripIn <- as.character(args$regex)
teloKmersIn <- as.character(strsplit(args$telokmers, ",")[[1]])
gapSizeIn <- as.numeric(args$gapSize)
teloParamIn <- as.numeric(strsplit(args$teloclust, ",")[[1]])
minChrSizeIn <- args$minChrLen
paletteIn <- viridisLite::viridis
verboseIn <- args$verbose
pdfFileIn <- args$pdfFile

# -- write the plot to file
plotx <- ceiling(sqrt(length(faFilesIn))) * 4 + 1.5
ploty <- floor(sqrt(length(faFilesIn))) * 3 + 1
pdf(pdfFileIn, height = ploty, width = plotx)

test <- classify_genomes(
  faFiles = faFilesIn,
  toStrip = toStripIn,
  teloKmers = teloKmersIn,
  gapSize = gapSizeIn,
  teloGapSize = teloParamIn[1],
  minTeloSize = teloParamIn[2],
  teloRepDens = teloParamIn[3],
  binScaffolds = binScaffoldsIn,
  minChrSize = minChrSizeIn,
  palette = paletteIn,
  verbose = verboseIn)

nope <- dev.off()
if(args$print2stdout){
  fwrite(test, sep = "\t", file = "", quote = FALSE)
}else{
  fwrite(test, file = args$outFile, sep = "\t")
}

# -- write the coordinates to the coordinates file


if(args$verbose)
  cat("Done!\n")
