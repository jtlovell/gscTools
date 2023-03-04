#!/usr/bin/env Rscript

################################################################################
################################################################################
# -- ad hoc functions
find_kmerPos <- function(dnass, kmer){
  m <- vmatchPattern(pattern = kmer, subject = dnass)
  o <- rbindlist(lapply(names(m), function(i)
    if(nrow(as.data.table(m[[i]])) > 0)
      data.table(chr = i, start = start(m[[i]])-1, end = end(m[[i]]))))
  if(nrow(o) > 0){
    return(o)
  }
}

cluster_kmerPos <- function(bed, minGap, minDens, minWidth){
  bed <- data.table(bed)
  setkey(bed, chr, start, end)
  bed[,tmp := c(0, diff(end) - 1)]
  bed$tmp[bed$tmp < minGap] <- 0
  bed[,jump := cumsum(tmp), by = "chr"]
  bed[,gapID := add_rle(jump, which = "id"), by = "chr"]
  bed[,nCovInGap := sum(end - start), by = c("chr", "gapID")]
  bed <- bed[,list(
    start = min(start),  end = max(end),
    nCovInGap = sum(end - start),
    gapWidth = max(end) - min(start)), 
    by = c("chr", "gapID")]
  bed[,propInGap := nCovInGap / gapWidth]
  bed <- subset(bed, propInGap >= minDens)
  bed[,width := end - start]
  bed <- subset(bed, width >= minWidth)
  return(bed[,c("chr", "start", "end", "width")])
}

add_rle <- function(x, which = "n"){
  if(which == "n"){
    rep(rle(x)$lengths, rle(x)$lengths)
  }else{
    rep(1:length(rle(x)$lengths), rle(x)$lengths)
  }
}

round_toInteger <- function(x, to){
  round(x/to, 0) * to
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

classify_engine <- function(faFile, teloKmers, gapSize, teloGapSize, minTeloSize, teloRepDens, binScaffolds, verbose){
  
  if(!requireNamespace("Biostrings", quietly = TRUE))
    stop("you need to install Biostrings from Bioconductor before running\n")
  
  if(!requireNamespace("data.table", quietly = TRUE))
    stop("you need to install data.table before running\n")
  
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))
  
  # -- read in the fasta as a DNAStringSet
  ss <- readDNAStringSet(faFile)
  chrLens <- pull_chrLens(ss)
  contigN50 <- N50(chrLens$chrEnd)
  if(verbose)
    cat(sprintf("\tN. scaffolds = %s (N50 = %sMb); total length = %sMb\n",
                nrow(chrLens), 
                round(N50(chrLens$chrEnd) / 1e6, 2), 
                round(N50(sum(chrLens$chrEnd)) / 1e6, 2)))
  
  # -- make a bed with the start and end positions of the gap kner (N)
  # then cluster into gap coordinates
  gaps <- find_kmerPos(dnass = ss, kmer = "N")
  gapsClus <- cluster_kmerPos(
    bed = gaps, minGap = gapSize, minDens = 1, minWidth = gapSize)
  
  # -- make a bed with the start and end positions of the telomere kmers
  # then cluster into gap coordinates
  teloKmers <- c(teloKmers, 
                 as.character(reverseComplement(DNAStringSet(teloKmers))))
  telos <- rbindlist(lapply(teloKmers, function(i)
    find_kmerPos(dnass = ss, kmer = i)))
  telos <- subset(telos, !duplicated(telos))
  telosClus <- cluster_kmerPos(
    bed = telos, 
    minGap = teloGapSize,
    minDens = teloRepDens, 
    minWidth = minTeloSize)
  
  # -- compile data
  telosClus[,type := "telomere"]
  gapsClus[,type := "gap"]
  kmerClus <- rbind(telosClus, gapsClus)
  
  # -- add the no-gap, no-telomere blocks
  outBed <- add_nogapBlocks(
    chrLens = chrLens, 
    kmerBed = kmerClus)
  
  # -- compile contig metadata
  contigBed <- add_nogapBlocks(
    chrLens = chrLens, 
    kmerBed = gapsClus)
  if(binScaffolds)
    contigBed$chr[grepl("scaff", contigBed$chr)] <- "scaffolds"
  
  contigs <- subset(contigBed, type != "gap")
  if(binScaffolds){
    contigs <- subset(contigs, chr != "scaffolds")
    if(verbose)
      cat(sprintf("\t%s contigs exl. scaffs (N50 = %s, total = %sMb); %s gaps (%sMb)\n",
                  nrow(contigs), 
                  round(N50(contigs$width) / 1e6, 2), 
                  round(N50(sum(contigs$width)) / 1e6, 2),
                  sum(outBed$type == "gap"),
                  round(sum(outBed$width[outBed$type == "gap"])/1e6, 2)))
  }else{
    if(verbose)
      cat(sprintf("\t%s contigs incl. scaffs (N50 = %s, total = %sMb); %s gaps (%sMb)\n",
                  nrow(contigs), 
                  round(N50(contigs$width)/1e6, 2), 
                  round(N50(sum(contigs$width))/1e6, 2),
                  sum(outBed$type == "gap"),
                  round(sum(outBed$width[outBed$type == "gap"])/1e6, 2)))
  }
  
  # -- compile telomere metadata
  teloBuff <- 5e4
  teloMd <- get_teloPos(
    bed = subset(outBed, type == "telomere"), 
    chrLens = chrLens, 
    teloBuff = teloBuff)
  nTelos <- sum(c(grepl("right", teloMd$teloPosition), 
                  grepl("left", teloMd$teloPosition)))
  nTeloChrs <- sum(grepl("right|left", teloMd$teloPosition))
  
  # -- compile full metadata
  
  contigMd <- contigBed[,list(
    nContigs = sum(type == "noGap"), 
    contigN50 = N50(width[type == "noGap"]),
    contigTotal = sum(width[type == "noGap"]),
    nGaps = sum(type == "gap"),
    gapTotal = sum(width[type == "gap"])), by = "chr"]
  
  outMd <- merge(contigMd, teloMd, by = "chr", all.x = T)
  
  if(binScaffolds){
    outTot <- data.table(
      chr = "tot. (no scafs)", 
      nContigs = nrow(contigs), 
      contigN50 = N50(contigs$width), 
      contigTotal = N50(sum(contigs$width)), 
      nGaps = sum(outBed$type == "gap"), 
      gapTotal = sum(outBed$width[outBed$type == "gap"]), 
      width = sum(chrLens$chrEnd[chrLens$chr %in% contigs$chr]), 
      teloPosition = sprintf(
        "%s (%sMb on %s Chrs)", 
        nTelos, 
        round(sum(subset(outBed, type == "telomere")$width) / 1e6, 2),
        sum(outBed$type == "telomere")))
  }else{
    outTot <- data.table(
      chr = "total", 
      nContigs = nrow(contigs), 
      contigN50 = N50(contigs$width), 
      contigTotal = N50(sum(contigs$width)), 
      nGaps = sum(outBed$type == "gap"), 
      gapTotal = sum(outBed$width[outBed$type == "gap"]), 
      width = sum(chrLens$chrEnd), 
      teloPosition = sprintf(
        "n = %s / %sMb / %s Chrs)", 
        nTelos, 
        round(sum(subset(outBed, type == "telomere")$width) / 1e6, 2),
        sum(outBed$type == "telomere")))
  }
  
  outMd <- rbind(outMd, outTot)
  
  return(list(metadata = outMd, blockCalls = outBed)) 
}

plot_classes <- function(blockCalls, minChrSize2plot = 1e6, palette = viridisLite::viridis, nColors = NULL, barThickness = .6,  backgroundColor = "darkgrey", ...){
  
  if(is.null(nColors)){
    mxn <- with(blockCalls, max(tapply(
      type, chr, FUN = function(x) sum(x == "noGap"))))
    nColors <- ifelse(mxn/3 < 5, 5, mxn / 3)
  }
    
  cols <- palette(nColors, ...)
  
  tp <- subset(blockCalls, chr != "scaffolds")
  if(!"genome" %in% colnames(tp))
    tp[,genome := "genome"]
  cl <- with(tp, tapply(end, chr, max))
  tp <- subset(tp, chr %in% names(cl)[cl > minChrSize2plot])
  
  tp[,index := rep(1:nColors, 100)[1:.N], by = c("genome", "chr")]
  tp[,y1 := as.numeric(as.factor(chr))]
  tp[,y2 := y1 - barThickness]
  tp[,y := (y1 + y2)/2]
  tp[,genome := factor(genome, levels = unique(genome))]
  ng <- tp[,list(n = .N, x = max(end)/1e6), by = c("genome", "chr", "y")]
  
  tel <- subset(tp, type == "telomere")
  tel[,x := ((end + start)/2)/1e6]
  p <- ggplot()+
    geom_rect(data = tp, 
              aes(xmin = start / 1e6, xmax = end / 1e6, 
                  ymin = y2, ymax = y1, 
                  fill = factor(index)))+
    scale_y_reverse(
      expand = c(0.01, 0.01), 
      breaks = unique(tp$y), 
      labels = unique(tp$chr),
      name = "Chromosome")+
    scale_x_continuous(
      name = "Physical position (Mb); * indicate telomere sequence")+
    geom_text(data = ng, 
              aes(x = x, y = y, label = n), 
              hjust = -.25, size = 2)+
    geom_point(data = tel, aes(x = x, y = y2), shape = "*", col = "red", size = 6)+
    
    scale_fill_manual(values = cols, guide = "none")+
    theme(panel.background = element_rect(fill = backgroundColor),
          panel.grid = element_blank())+
    facet_wrap(~genome, scale = "free")+
    ggtitle(sprintf(
      "Contig positions: each cycle through colors is %s contigs (%s gaps)", 
      nColors, nColors - 1))
  print(p)
}

classify_genomes <- function(faFiles, toStrip = ".fa$|.fa.gz$", teloKmers = c("CCCGAAA", "CCCTAAA"), gapSize = 1e4, teloGapSize = 1e4, minTeloSize = 1e3, teloRepDens = 0.1, binScaffolds = TRUE, minChrSize2plot = 1e6, palette = viridisLite::viridis, nColors = 25, barThickness = .6, backgroundColor = "darkgrey", verbose = TRUE, ...){
  
  classes <- lapply(faFiles, function(i){
    cat(i,"\n")
    test <- classify_engine(
      faFile = i, 
      teloKmers = teloKmers,
      gapSize = gapSize,
      teloGapSize = teloGapSize,
      minTeloSize = minTeloSize,
      teloRepDens = teloRepDens, 
      binScaffolds = binScaffolds, 
      verbose = verbose)
    test$metadata[,genome := gsub(toStrip, "", basename(i))]
    test$blockCalls[,genome := gsub(toStrip, "", basename(i))]
    return(test)
  })
  
  # -- parse the results
  plotThis <- rbindlist(lapply(classes, function(x) x$blockCalls))
  metaThis <- rbindlist(lapply(classes, function(x) x$metadata))
  
  # -- make the plot
  plot_classes(
    blockCalls = plotThis,
    minChrSize2plot = minChrSize2plot,
    palette = palette, 
    nColors = nColors,
    barThickness = barThickness,
    backgroundColor = backgroundColor,
    ...)
  
  return(list(metadata = metaThis, plotdata = plotThis))
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
  "-z", "--teloclust", type = "character", default = "1e4,1e3,0.1", 
  help = "comma-separated cluster params (max gap, min width, min %telo sequence) [default = %(default)s]")
parser$add_argument(
  "-b", "--binScaffolds", action="store_false", dest = "binScaffolds", 
  help = "Should any chrs with scaff in their name be binned [default = %(default)s]")

parser$add_argument(
  "-p", "--pdfFile", type = "character", default = "contigPlot.pdf",
  help = "output pdf file name [default= %(default)s]")
parser$add_argument(
  "-o", "--outFile", type = "character", default = "contigCoords.txt",
  help = "output text file with contig coordinates [default= %(default)s]")

parser$add_argument(
  "-m", "--minChrLen2plot", type = "integer", default = 1e6,
  help = "the smallest chromosome to show [default= %(default)s]")
parser$add_argument(
  "-g", "--gapSize", type = "integer", default = 1e4,
  help = "window size to average over to merge gaps [default= %(default)s]")
parser$add_argument(
  "-d", "--height", type = "double", default = 3,
  help = "numeric, the height of the plot [default= %(default)s]")
parser$add_argument(
  "-w", "--width", type = "double", default = 3,
  help = "numeric, the width of the plot / genome [default= %(default)s]")
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
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

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
binScaffoldsIn <- args$binScaffolds
minChrSize2plotIn <- args$minChrLen2plot
paletteIn <- viridisLite::viridis
nColorsIn <- NULL
barThicknessIn <- .6
backgroundColorIn <- "darkgrey"
verboseIn <- args$verbose
pdfFileIn <- args$pdfFile

# -- write the plot to file
plotx <- floor(sqrt(length(faFilesIn))) * 4 + .5
ploty <- ceiling(sqrt(length(faFilesIn))) * 3
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
  minChrSize2plot = minChrSize2plotIn, 
  palette = paletteIn, 
  nColors = nColorsIn, 
  barThickness = barThicknessIn, 
  backgroundColor = backgroundColorIn, 
  verbose = verboseIn)

nope <- dev.off()

fwrite(test$plotdata, file = args$outFile, sep = "\t")

# -- write the coordinates to the coordinates file
if(args$print2stdout)
  fwrite(test$metadata, sep = "\t", file = "", quote = FALSE)

if(args$verbose)
  cat("Done!\n")
