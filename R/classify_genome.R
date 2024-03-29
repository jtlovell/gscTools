#' @title classify the genome based on bed coordinates
#'
#' @description
#' \code{classify_genome} Hierarchically classify the genome by coordinates in
#' the list of bed files. Sequences in the first bed file are masked in the
#' second and so on.
#'
#' @param seqInfo seqinfo extracted from the genome assembly DNAstringset
#' @param listOfBeds list of bed-like data.frames/data.tables with at least
#' the columns chr, start and end.
#' @param verbose logical, should updates be printed to the console?
#'
#' @details intersects, gaps and reduces overlapping intervals
#'
#' @return a list of granges matchin the listOfBeds
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges reduce union
#' @export
classify_genome <- function(seqInfo, listOfBeds, verbose){

  if(length(listOfBeds) < 1)
    stop("listOfBeds must contain at least one element\n")

  if(!is.list(listOfBeds))
    stop("listOfBeds must be a list of data.frames\n")

  if(!all(sapply(listOfBeds, is.data.frame)))
    stop("all elements in listOfBeds must be a data.frame\n")

  if(!all(sapply(listOfBeds, function(x)
    all(c("start", "end", "chr") %in% colnames(x)))))
    stop("all elements in listOfBeds must have column names: chr, start, end\n")

  if(is.null(names(listOfBeds)))
    stop("listOfBeds must be a named list\n")

  if(sum(duplicated(names(listOfBeds))) > 0)
    stop("listOfBeds had duplicated names\n")

  ##############################################################################
  # 1. Read the sequence lengths
  si <- seqInfo
  if(verbose)
    cat(sprintf(
      " (%s Mb in %s chrs)\n",
      round(sum(seqlengths(si))/1e6, 2),
      length(si)))

  ##############################################################################
  # 2. Loop through the bedfiles
  if(verbose)
    cat(sprintf("Classifying the genome by %s bed files ... \n", length(listOfBeds)))
  maskThis <- NULL
  outList <- list()
  for(i in names(listOfBeds)){
    inGr <- join_ranges(listOfBeds[[i]], seqInfo = si)
    if(is.null(maskThis)){
      maskThis <- inGr
      outGr <- inGr
    }else{
      maskThis <- reduce(maskThis)
      outGr <- get_gappedRanges(
        existingGr = maskThis, newGr = inGr)
      tmp <- GenomicRanges::union(outGr, maskThis)
      maskThis <- reduce(tmp)
    }
    outList[[i]] <- outGr
    if(verbose)
      cat(sprintf("\t%sMb %s sequence (%sMb masked)\n",
                  sum(round(sum(width(outGr))/1e6, 2)),
                  i,
                  sum(round(sum(width(maskThis))/1e6, 2))))
  }

  ##############################################################################
  # 3. Classify un-annotated sequence
  unk <- find_gaps(maskThis)
  if(verbose)
    cat(sprintf("\t%sMb un-annotated sequence\n",
                sum(round(sum(width(unk))/1e6, 2))))
  outList[["unknown"]] <- unk
  return(outList)
}
