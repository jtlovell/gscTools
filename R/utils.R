#' @title Generic internal functions used by genespace
#' @description
#' \code{utils} Convenience functions for gscTools, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param bed data.table or data.frame with at least chr, start, end columns
#' @param dnass DNAStringSet
#' @param seqInfo seqInfo ranges extracted from dnass
#' @param gr granges object
#' @param existingGr existing granges to mask with
#' @param newGr new granges to mask
#' @param kmers DNAStringSet with kmers to find
#' @param nCores integer, then number of cores to use
#' @param x single-value parameter, string, integer, numeric, list, vector
#' @param ... additional arguments passed on
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'
#' @title startup messages
#' @description
#' \code{.onAttach} startup messages
#' @rdname utils
#' @export
.onAttach <- function(...) {
  packageStartupMessage(paste(strwrap(
    "gscTools v0.0.1: tools for genome QC and viz\n",
    indent = 0, exdent = 8), collapse = "\n"))
}

#' @title join_ranges
#' @description
#' \code{join_ranges} join_ranges
#' @rdname utils
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
join_ranges <- function(bed, seqInfo){
  gr <- with(bed, GRanges(
    chr,
    IRanges(start = start, end = end),
    seqinfo = seqInfo))
  return(reduce(gr))
}

#' @title pull_seqInfo
#' @description
#' \code{pull_seqInfo} pull_seqInfo
#' @rdname utils
#' @importFrom GenomeInfoDb Seqinfo
#' @export
pull_seqInfo <- function(dnass){
  ss <- dnass
  fai <- data.table(
    chr = names(ss),
    start = 1,
    end = width(ss))
  seqi <- Seqinfo(fai$chr, seqlengths = fai$end)
  return(seqi)
}

#' @title find_gaps
#' @description
#' \code{find_gaps} find_gaps
#' @rdname utils
#' @importFrom GenomicRanges gaps reduce
#' @export
find_gaps <- function(gr){
  gapr <- gaps(reduce(gr))
  return(reduce(subset(gapr, strand == "*")))
}

#' @title subset_bedInGap
#' @description
#' \code{subset_bedInGap} subset_bedInGap
#' @rdname utils
#' @importFrom GenomicRanges intersect reduce
#' @export
subset_bedInGap <- function(gr, gapGr){
  return(intersect(
    x = reduce(gr),
    y = reduce(gapGr)))
}

#' @title get_gappedRanges
#' @description
#' \code{get_gappedRanges} get_gappedRanges
#' @rdname utils
#' @export
get_gappedRanges <- function(existingGr, newGr){
  ingap <- subset_bedInGap(
    gr = newGr,
    gapGr = find_gaps(existingGr))
  return(ingap)
}

#' @title convert_si2gr
#' @description
#' \code{convert_si2gr} convert_si2gr
#' @rdname utils
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths
#' @export
convert_si2gr <- function(seqInfo){
  gr <- GRanges(
    names(seqInfo),
    IRanges(start = 1, end = seqlengths(seqInfo)),
    seqinfo = seqInfo)
  return(reduce(gr))
}

#' @title rename chromosomes
#' @description
#' \code{rename_chrs} Parse common formats for chromosomes, specifically
#' splitting whitespace separated fields and finding the field like to give
#' the chromosome ID
#' @rdname utils
#' @import data.table
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths
#' @export
rename_chrs <- function(dnass){
  ss <- dnass
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
    hasChr <- gsub("chr|chromosome|scaffold|chr0", "", tolower(u1)[[1]]) == "1"
    if(!any(duplicated(u1)) && hasChr){
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

