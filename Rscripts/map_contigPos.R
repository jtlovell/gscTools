#!/usr/bin/env Rscript
################################################################################
################################################################################
# ad hoc functions to do the work
################################################################################
find_contigGaps <- function(faFiles, 
                            gapKmer = "N",
                            gapSize = 1e4,
                            verbose = T){
  
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
  
  if(!requireNamespace("Biostrings", quietly = TRUE))
    stop("you need to install Biostrings from Bioconductor before running\n")
  
  if(!requireNamespace("data.table", quietly = TRUE))
    stop("you need to install data.table before running\n")
  
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))
  
  if(is.null(names(faFiles))){
    warning("no names given to the genome faFiles. Assigning genome_1:N\n")
    names(faFiles) <- sprintf("genome_%s", 1:length(faFiles))
  }
  
  if(any(duplicated(names(faFiles))))
    stop("some names of the genomes (faFiles vector names) are duplicated\n")
  
  gps <- rbindlist(lapply(names(faFiles), function(i){
    dnass <- readDNAStringSet(faFiles[i])
    m <- vmatchPattern(pattern = gapKmer, subject = dnass)
    o <- rbindlist(lapply(names(m), function(i)
      if(nrow(as.data.table(m[[i]])) > 0)
        data.table(chr = i, start = start(m[[i]])-1, end = end(m[[i]]))))
    o[,tmp := c(0, diff(end) - 1)]
    o$tmp[o$tmp < gapSize] <- 0
    o[,jump := cumsum(tmp), by = "chr"]
    o[,gapID := add_rle(jump, which = "id"), by = "chr"]
    o <- o[,list(start = min(start),  end = max(end)), by = c("chr", "gapID")]
    
    if(verbose)
      cat(sprintf("%s: found %s gaps containing %s Mbp on %s sequences\n",
                  i, nrow(o), round(with(o, sum(end - start))/1e6, 3), uniqueN(o$chr)))
    if(nrow(o) > 0){
      return(data.table(o, genome = i))
    }
  }))
  return(gps)
}

################################################################################
pull_chrLens <- function(faFiles){
  gps <- rbindlist(lapply(names(faFiles), function(i){
    dnass <- readDNAStringSet(faFiles[i])
    cl <- data.table(
      genome = i, chr = names(dnass), chrStart = 1, chrEnd = width(dnass))
    return(cl)
  }))
  return(gps)
}

################################################################################
map_contigPos <- function(faFiles, 
                                palette = rainbow,
                                nColors = 10, 
                                barThickness = .6,
                                gapSize = 1e4,
                                verbose = TRUE,
                                gapKmer = "N",
                                minChrLen2plot = 1e6, 
                                backgroundColor = "darkgrey", 
                                ...){
  
  
  
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("you need to install ggplot before running\n")
  
  if(!requireNamespace("data.table", quietly = TRUE))
    stop("you need to install data.table before running\n")
  
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ggplot2))
  
  gaps <- find_contigGaps(
    faFiles = faFiles, 
    gapSize = gapSize,
    gapKmer = "N",
    verbose = verbose)
  seqLens <- pull_chrLens(faFiles)
  
  tp <- merge(seqLens, gaps, all.x = T, by = c("genome", "chr"))
  setkey(tp, genome, chr, start, end)
  tp[,index := 1:.N, by = c("genome", "chr")]
  tp[,`:=`(startLeft = c(1, end[-.N]+1),
           endLeft = c(start-1)), 
     by = c("genome", "chr")]
  tpo <- with(tp, data.table(
    genome = c(genome, genome),
    chr = c(chr, chr),
    type = rep(c("gap", "noGap"), each = nrow(tp)),
    start = c(start, startLeft),
    end = c(end, endLeft)))
  tpo <- subset(tpo, complete.cases(tpo))
  tpe <- tp[,list(start = ifelse(all(is.na(end)), 1L, as.integer(max(end, na.rm = T)+1)), 
                  end = chrEnd[1]), by = c("genome", "chr")]
  tpe[,type := "noGap"]
  tp <- rbind(tpo, tpe)
  setkey(tp, genome, chr, start, end)
  
  nogap <- subset(tp, type = "noGap")
  nogap[,index := 1:.N]
  
  cols <- palette(nColors, ...)
  nogap[,index := rep(1:nColors, 100)[1:.N], by = c("genome", "chr")]
  nogap[,y1 := as.numeric(as.factor(chr))]
  nogap[,y2 := y1 - barThickness]
  nogap[,y := (y1 + y2)/2]
  nogap[,genome := factor(genome, levels = names(faFiles))]
  goodchrs <- subset(seqLens, chrEnd > minChrLen2plot)
  goodchrs[,n := uniqueN(genome), by = "chr"]
  goodchrs <- subset(goodchrs, n > (uniqueN(genome) / 2))
  nogapsub <- subset(nogap, chr %in% unique(goodchrs$chr))
  ng <- nogapsub[,list(n = .N, x = max(end)/1e6), by = c("genome", "chr", "y")]
  p <- ggplot()+
    geom_rect(data = nogapsub, 
              aes(xmin = start / 1e6, xmax = end / 1e6, 
                  ymin = y2, ymax = y1, 
                  fill = factor(index)))+
    scale_y_continuous(
      expand = c(0.01, 0.01), 
      breaks = unique(nogap$y), 
      labels = unique(nogap$chr),
      name = "Chromosome")+
    scale_x_continuous(
      name = "Physical position (Mb)")+
    geom_text(data = ng, 
              aes(x = x, y = y, label = n), 
              hjust = -.1, size = 2)+
    scale_fill_manual(values = cols, guide = "none")+
    theme(panel.background = element_rect(fill = backgroundColor),
          panel.grid = element_blank())+
    facet_grid(.~genome)+
    ggtitle(sprintf(
      "Contig positions: each cycle through colors is %s contigs (%s gaps)", 
      nColors, nColors - 1))
  print(p)
  
  return(tp)
}

################################################################################
################################################################################
# -- parse the arguments
################################################################################

if (!require('optparse', quietly = T))
  install.packages('optparse')

suppressPackageStartupMessages(library("optparse"))


library("optparse")

option_list = list(
  make_option(
    c("-f", "--faFiles"), 
    type = "character", 
    default = "", 
    help = "comma-separated quoted string of paths to fasta files",
    metavar = "character"),
  
  make_option(
    c("-p", "--pdfFile"), 
    type = "character", 
    default = "contigPlot.pdf",
    help = "output pdf file name [default= %default]",
    metavar = "character"),
  
  make_option(
    c("-c", "--outFile"), 
    type = "character", 
    default = "contigCoords.txt",
    help = "output text file with contig coordinates [default= %default]",
    metavar = "character"),
  
  make_option(
    c("-m", "--minChrLen2plot"), 
    type = "integer", 
    default = "1e6",
    help = "the smallest chromosome to show [default= %default]",
    metavar = "integer"),
  
  make_option(
    c("-g", "--gapSize"), 
    type = "integer", 
    default = "1e4",
    help = "window size to average over to merge gaps [default= %default]",
    metavar = "integer"),
  
  make_option(
    c("-h", "--height"), 
    type = "double", 
    default = 3,
    help = "numeric, giving the height of the plot [default= %default]",
    metavar = "double"),
  
  make_option(
    c("-w", "--width"), 
    type = "double", 
    default = 3,
    help = "numeric, giving the width of the plot per genome [default= %default]",
    metavar = "double"),
  
  make_option(
    c("-r", "--regex"), 
    type = "character", 
    default = ".fa.gz",
    help = "quoted character, giving [default= %default]")
  )

opt_parser <- OptionParser(
  usage = "%prog [options] -f <list,of,fasta,files> ...",
  option_list = option_list,
  add_help_option=FALSE)
opt <- parse_args(opt_parser)

if ("help" %in% names(opt)){
  print_help(opt_parser)
  stop("Please provide an input file")
}

if(!"faFiles" %in% names(opt))
  stop("must provide path to atleast one fasta file\n")
faFiles <- strsplit(gsub(" ", "", opt$faFiles), ",")[[1]]
if(!all(file.exists(faFiles)))
  stop("the following fasta files could not be found: ",
       paste(faFiles[!file.exists(faFiles)], collapse = ", "))
names(faFiles) <- gsub(opt$regex, "", basename(faFiles))

if(opt$faFiles == ""){
  print_help(opt_parser)
  stop("Please provide an input file")
}

pal <- colorRampPalette(colors = c(
  "#440154FF", "#471164FF", "#481F70FF", "#472D7BFF", "#443A83FF", "#404688FF",
  "#3B528BFF", "#365D8DFF", "#31688EFF", "#2C728EFF", "#287C8EFF", "#24868EFF", 
  "#21908CFF", "#1F9A8AFF", "#20A486FF", "#27AD81FF", "#35B779FF", "#47C16EFF",
  "#5DC863FF", "#75D054FF", "#8FD744FF", "#AADC32FF", "#C7E020FF", "#E3E418FF", 
  "#FDE725FF"))

pdf(opt$pdfFile, height = opt$height, width = opt$width * length(faFiles))
dat <- map_contigPos(
  faFiles = faFiles, 
  palette = pal, 
  gapSize = as.numeric(opt$gapSize),
  minChrLen2plot = as.numeric(opt$minChrLen2plot),
  backgroundColor = "darkgrey",
  nColors = 25)
nope <- dev.off()


fwrite(dat, sep = "\t", file = opt$outFile, quote = FALSE)
cat("Done!\n")
