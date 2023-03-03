#!/usr/bin/env Rscript

if (!require('argparse'))
  install.packages('argparse')

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument(
  "-w", "--wd", default = "NULL",
  help = "file.path, where are your cov.bed files stored?")
parser$add_argument(
  "-d", "--mxDepth", type = "double", default = 160,
  help = "maximum per-site depth")
parser$add_argument(
  "-c", "--maxCorDepth", type = "double", default = 4,
  help = "maximum genome-standardized depth")
parser$add_argument(
  "-p", "--pdfFile", default = "NULL",
  help = "file.path to where the pdf file should be written to")
parser$add_argument(
  "-t", "--pdfHt", type = "integer", default = 16,
  help = "plot window height")
parser$add_argument(
  "-f", "--pdfWd", type = "integer", default = 16,
  help = "plot window width")
parser$add_argument(
  "-x", "--excludeChrString", default = "NULL",
  help = "character string - if the chr ID contains this, exclude it")
parser$add_argument(
  "-v", "--verbose", default = "TRUE",
  help = "TRUE or FALSE, should updates be printed to the console?")

args <- parser$parse_args()

if(args$wd == "NULL"){
  wd <- getwd()
}else{
  if(!dir.exists(wd))
    stop(sprintf("specified wd %s does not exist\n", wd))
}

covs <- list.files(pattern = "cov.bed", path = wd, full.names = T)
if(length(covs) == 0)
  stop(sprintf("coult not find any file names containing cov.bed in wd", wd))


mxDepth <- as.numeric(args$mxDepth)
if(is.na(mxDepth))
  stop(sprintf("could not coerce max depth = %s to integer", args$mxDepth))

maxCorDepth <- as.numeric(args$maxCorDepth)
if(is.na(mxDepth))
  stop(sprintf("could not coerce max cor. depth = %s to integer", args$maxCorDepth))

pdfFile <- args$pdfFile
if(pdfFile == "NULL")
  pdfFile <- file.path(wd, "ccsCoverage.pdf")
if(!dir.exists(dirname(pdfFile)))
  stop(sprintf("directory for pdf file: %s does not exist", args$pdfFile))

pdfHt <- as.numeric(args$pdfHt)
if(is.na(pdfHt))
  stop(sprintf("could not coerce pdfHt = %s to integer", args$pdfHt))

pdfWd <- as.numeric(args$pdfWd)
if(is.na(pdfWd))
  stop(sprintf("could not coerce pdfWd = %s to integer", args$pdfWd))

excludeChrString <- as.character(args$excludeChrString)

verbose <- args$verbose
verbose <- as.logical(verbose)
if(!is.na(verbose))
  verbose <- TRUE

if (!require('data.table'))
  install.packages('data.table')
if (!require('ggplot2'))
  install.packages('ggplot2')
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))

if(verbose)
  cat(sprintf("Reading %s coverage files ... ", length(covs)))

d <- rbindlist(lapply(covs, function(i){
  x <- fread(i, col.names = c("tmp", "pos", "depth"))
  x[,genome := gsub("_cov.bed", "", basename(i))]
  setnames(x, 1, "chr")
  return(x)
}))

d$depth[d$depth > mxDepth] <- mxDepth
d[,corDepth := depth/mean(depth), by = "genome"]
d$corDepth[d$corDepth > maxCorDepth] <- maxCorDepth

d <- subset(d, !grepl(excludeChrString, chr))

if(verbose)
  cat(sprintf("Done!\nWriting plots to %s ... ", pdfFile))
pdf(pdfFile, height = pdfHt, width = pdfWd)
for(i in unique(d$genome)){
  x <- subset(d, genome == i)
  p <- ggplot(d, aes(x = pos/1e6, col =  corDepth, y =  corDepth))+
    geom_point(pch = ".")+
    facet_wrap(~ chr, scale = "free")+
    scale_color_gradient(low = "dodgerblue", high = "gold")+
    stat_smooth(span = .2)+
    # scale_y_log10()+
    labs(x = "Physical position (Mbp)", y = "Coverage (X)",
         main = sprintf("%s CCS coverage", i))+
    theme(panel.background = element_rect(fill = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "lightgrey", size = .1, linetype = 2),
          panel.grid.minor.y = element_blank())
  suppressMessages(print(p))
}
tmp <- dev.off()

if(verbose)
  cat("Done!\n")
