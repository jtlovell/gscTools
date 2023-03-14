#' @title Generic internal functions used by genespace
#' @description
#' \code{find_kmers} Convenience functions for gscTools, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name find_kmers
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
#' If called, \code{find_kmers} returns its own arguments.
#'
#'
#' @title Fast method to map many kmers
#' @description
#' \code{find_manyKmers} Finds exact DNA-DNA matches and aggregates across all
#' kmers in the kmer DNAstring
#' @rdname find_kmers
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings PDict matchPDict
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom BiocGenerics start end do.call
#' @export
find_manyKmers <- function(dnass,
                           kmers,
                           nCores = 1){
  pdictKmer <- PDict(kmers)
  grlist <- mclapply(names(dnass), mc.cores = nCores, function(i){
    tmp <- GRanges(matchPDict(
      pdict = pdictKmer,
      subject = dnass[[i]]))
    tmp <- reduce(
      tmp,
      drop.empty.ranges = T)
    if(length(tmp) > 1){
      out <- GRanges(
        seqnames = i,
        IRanges(start = start(tmp),
                end = end(tmp)))
      return(out)
    }
  })
  if(length(grlist) == 0){
    stop("could not find any matches to the kmers\n")
  }else{
    suppressWarnings(outgr <- BiocGenerics::do.call(
      c,
      grlist[sapply(grlist, length) > 0]))
    return(reduce(outgr))
  }
}

#' @title Fast method to map many kmers
#' @description
#' \code{find_manyKmers} Finds exact DNA-DNA matches and aggregates across all
#' kmers in the kmer DNAstring
#' @rdname find_kmers
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings PDict matchPDict
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom BiocGenerics start end do.call
#' @export
find_fewKmers <- function(dnass,
                          kmers,
                          nCores = 1){
  si <- pull_seqInfo(dnass)

  kmerposList <- mclapply(kmers, mc.cores = nCores, function(i)
    vmatchPattern(pattern = as.character(i), dnass))

  # -- remove any chromosomes without any Ns
  kmerposList <- kmerposList[sapply(kmerposList, length) > 0]

  # -- combine grs with Ns
  kmerpos <- BiocGenerics::do.call(c, lapply(kmerposList, function(kmeri){
    kmeri <- kmeri[sapply(kmeri, length) > 0]
    kmeri <- BiocGenerics::do.call(c, lapply(names(kmeri), function(i)
      GRanges(i, kmeri[[i]], seqinfo = si)))
    return(kmeri)
  }))
  return(kmerpos)
}

#' @title Find contig gap positions as a run of N's
#' @description
#' \code{find_runsOfNs} Exact string matching of Ns to a DNAStringSet to find
#' the position of gaps between contigs
#' @rdname find_kmers
#' @import data.table
#' @importFrom Biostrings vmatchPattern
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom BiocGenerics start end do.call
#' @export
find_runsOfNs <- function(dnass, minRunLength){
  # -- get sequence lengths
  si <- pull_seqInfo(dnass)

  # -- find all occurances of N
  nposList <- vmatchPattern(pattern = "N", dnass)

  # -- remove any chromosomes without any Ns
  nposList <- nposList[sapply(nposList, length) > 0]

  # -- combine grs with Ns
  npos <- BiocGenerics::do.call(c, lapply(names(nposList), function(i)
    GRanges(i, nposList[[i]], seqinfo = si)))

  npos <- reduce(
    npos,
    min.gapwidth = 2,
    ignore.strand = TRUE,
    drop.empty.ranges = TRUE)
  return(npos)
}

#' @title find_telomeres
#' @description
#' \code{find_telomeres} find_telomeres
#' @rdname find_kmers
#' @import data.table
#' @importFrom Biostrings vmatchPattern
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom BiocGenerics start end do.call
#' @export
find_telomeres <- function(dnass,
                           kmers,
                           maxDistBtwTelo = 20,
                           minTeloSize = 200,
                           minTeloDens = 0.75,
                           maxDist2end = 1e4,
                           nCores = 1){
  kmerpos <- find_fewKmers(
    dnass = dnass,
    kmers = kmers,
    nCores = nCores)

  if(is.null(kmerpos))
    return(NULL)

  if(uniqueN(nchar(kmers)) != 1)
    warning("variable kmer lengths ... using longest for density calculations\n")
  kmerout <- reduce(
    kmerpos,
    min.gapwidth = maxDistBtwTelo,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE,
    with.revmap = TRUE)

  dens <- sapply(kmerout@elementMetadata$revmap, length) * max(nchar(kmers))
  wid <- width(kmerout)
  dens <- dens / wid
  kmersub <- subset(kmerout, dens >= minTeloDens & wid >= minTeloSize)
  isLeft <- start(kmersub) <= maxDist2end
  rightBound <- (seqlengths(si) - maxDist2end)[as.character(seqnames(kmersub))]
  isRight <- end(kmersub) >= rightBound
  teloclass <- ifelse(isLeft, "left", ifelse(isRight, "right", "inter"))
  kmerout <- GRanges(
    seqnames = seqnames(kmersub),
    IRanges(start = start(kmersub), end = end(kmersub)),
    position = teloclass,
    seqinfo = si)
  return(kmerout)
}




