#' @title Parse and plot bam files
#'
#' @description
#' \code{strandBam} Parses bam file(s) by the position of genes of interest, then
#' categorize the reads by their flags and plot them against the gene model
#'
#' @param gff.file Character string indicating the name and file location of the
#' gff3 formatted annotation file
#' @param bam.files Character vector of the name(s) of the bam files to parse and concatenate. If the filese
#' are not in the working directory, the file location must also be included in
#' each bam.files name
#' @param genes2test Character string of the genes to examine. These must match the
#' ids in the gff file. See below
#' @param geneID.start The position of the first character of the gene id specified in
#' the gff3 file. This is typically 4, because the leading 3 characters are "ID=".
#' @param geneID.end The position of the last character of the gene id. Depends on
#' the format of the gene IDs. Typically a good idea to include both the gene identifier
#' and the splice variant so that all gene models are unique.
#' @param min.dist The minimum distance (bp) between any two genes. If all genes (regardless
#' whether they overlap) are desired, set to NA
#' @param buffer The bp distance around a gene model to plot and pull reads
#' @param plotit Logical, should plots be made and printed to a pdf?
#' @param plot.pdf The name and file location of the output pdf (if plotit)
#' @param verbose Logical, should updates be printed?
#' @details Here, we do some stuff. More here soon.
#' @return A dataframe with the positions and flags of reads relating to each gene
#'
#' @examples
#' # more here soon
#'
#' @import IRanges Rsamtools plyr
#' @importFrom grDevices dev.off pdf rgb
#' @importFrom graphics arrows axis par plot rect segments title
#' @importFrom utils read.delim
#' @export
strandBam<-function(gff.file, bam.files, genes2test,
                    min.dist = 120, buffer = 100,
                    plotit = T, verbose = T,
                    geneID.start = 4, geneID.end = 17,
                    plot.pdf = "geneModel_reads.pdf"){

  add.min.dist = function(gff){
    h1g<-gff[gff$V3 == "gene",]
    h1gs = gsub(".v3.1","",gsub(".v2.1","",
                                gsub("ID=","",
                                     sapply(h1g$V9,function(x)
                                       strsplit(x,";", fixed = T)[[1]][1])),
                                fixed = T),
                fixed = T)

    h1g$geneID = h1gs
    h1g$min = apply(h1g[,4:5],1,min)
    h1g$max = apply(h1g[,4:5],1,max)

    h1g$dist1<-c(sapply(1:(nrow(h1g)-1), function(x){
      (min(h1g[x+1,c("min","max")])-max(h1g[x,c("min","max")]))
    }),1000)

    h1g$dist2<-c(1000,sapply(2:(nrow(h1g)), function(x){
      (min(h1g[x,c("min","max")])-max(h1g[x-1,c("min","max")]))
    }))

    h1g$min.dist = apply(h1g[,c("dist1","dist2")],1,min)
    return(h1g)
  }


  readIndexBam<-function(bamfiles, index.it = F,
                         chr, start, end, geneID,
                         what = Rsamtools::scanBamWhat(), verbose = T){

    which <- IRanges::RangedData(space=chr,
                                 IRanges::IRanges(as.numeric(start),
                                                  as.numeric(end)))
    param <- Rsamtools::ScanBamParam(which = which,
                                     what= what)

    if(index.it){
      if(verbose) cat("indexing",length(bamfiles),"bamfiles\n")
      Rsamtools::indexBam(bamfiles)
    }

    if(verbose) cat("scanning",length(bamfiles),"bamfiles ... \n")
    bamlist = lapply(bamfiles, function(x){
      if(verbose) cat(" ...",x)

      dat <- Rsamtools::scanBam(file=x, param = param)
      names(dat)<-geneID
      if(verbose) cat("\tn.regions = ",length(dat),"and total reads =")

      ldat = lapply(names(dat), function(y){
        if(length(dat[[y]]$flag)>0){
          with(dat[[y]],
               data.frame(library = x,
                          geneID = y,
                          flag = flag,
                          start = as.numeric(pos),
                          end = as.numeric(pos)+as.numeric(qwidth),
                          stringsAsFactors = F))
        }
      })

      tp = do.call(rbind, ldat)
      if(verbose) cat(nrow(tp),"\n")

      return(tp)
    })
    cat("done")
    return(do.call(rbind, bamlist))
  }

  plotGeneModel<-function(gff3, geneID, downstreamBuffer = 100, upstreamBuffer = 100){

    # --- 1. Format gff
    if(!ncol(gff3)==9) stop("gff3 file must be in standard format")

    gff3 = gff3[,c(1,3,4,5,7,9)]
    colnames(gff3)<-c("chr","feature","start","end","strand","id")
    gff3$start=as.numeric(as.character(gff3$start))
    gff3$end=as.numeric(as.character(gff3$end))

    gff3<-gff3[grep(geneID,gff3$id),]

    if(nrow(gff3)==0) {
      stop("geneID not found in the gff file\n")
    }
    gff3<-gff3[order(gff3$end),]

    # --- 2. Get orientation
    orientation <- ifelse(gff3$strand[1]=="+","forward","reverse")
    dir <- ifelse(orientation == "forward",1,-1)

    # --- 3. calculate standardized positions relative to the transcription start site
    if(orientation == "forward"){
      xrange<-c(min(c(gff3$start,gff3$end))-upstreamBuffer,max(c(gff3$start,gff3$end))+downstreamBuffer)
    }else{
      xrange<-c(min(c(gff3$start,gff3$end))-downstreamBuffer,max(c(gff3$start,gff3$end))+upstreamBuffer)
    }

    # --- 4. Make the basic plot
    ### make the plot window
    chr = gff3$chr[1]
    plot(1,type = "n", xlim=xrange, ylim=c(-1,1),
         yaxt="n", bty="n",
         xlab = paste(chr,"position (bp)"), ylab="")

    ### Draw the strand as an arrow
    arrows(x0=xrange[1], x1=xrange[2],
           y0=0, y1=0,
           length = .1, code = ifelse(orientation == "forward",2,1),
           col = rgb(0,0,0,.4))

    ### Add the gene as a dark line segment
    gene = gff3[gff3$feature == "gene", ]
    with(gene, segments(x0 = start, x1 = end, y0 = 0, y1 = 0,
                        col = "black", lwd = 2))

    ### Add exons as dark boxes
    exons = gff3[gff3$feature == "CDS", ]
    apply(exons[,c("start","end")],1,function(x){
      rect(xleft = x[1],xright = x[2],ytop = .02,ybottom = -.02,
           border = "dodgerblue4", col = "lightsteelblue1")
    })

    # --- 5. Form the introns

    if(nrow(exons)>1){
      introns = exons[-1,]
      for(i in 1:nrow(introns)){
        introns$start[i]=exons$end[i]
        introns$end[i]=exons$start[i+1]
      }
    }
    ### Add intron line segments
    introns$mids<-with(introns, (start+end)/2)
    with(introns, segments(x0=end,x1=mids,y0=0.02,y1=0.04, col="dodgerblue4", lwd=1))
    with(introns, segments(x0=mids,x1=start,y0=0.04,y1=0.02, col="dodgerblue4", lwd=1))
  }



  if(verbose) cat("reading and parsing the gff...")
  gff.file = "/global/projectb/scratch/jtlovell/hallii_genomes/HAL2/PhalliiHALv2.1.gene.gff3"
  gff = read.delim(gff.file,
                   comment.char = "#",
                   header = F,
                   stringsAsFactors = F)
  gff$geneID = substr(gff$V9,geneID.start,geneID.end)
  if(verbose) cat(" done.\n")

  if(!is.na(min.dist)){
    if(verbose) cat("calculating minimum distance between genes... ")
    gff.mindist = add.min.dist(gff = gff)

    gff.mindist <- gff.mindist[gff.mindist$min.dist>min.dist,]

    genes2test = genes2test[genes2test %in% gff.mindist$geneID]
    if(verbose) cat("retaining", length(genes2test), "genes <", min.dist, "bp from nearest neighbor\n")
  }

  bed = gff.mindist[gff.mindist$geneID %in% genes2test, ]
  bed<-bed[!duplicated(bed[,c(1,4,5)]),]
  gff.mindist<-gff.mindist[!duplicated(gff.mindist[,c(1,4,5)]),]
  good.flags = data.frame(flag = c(83,99,147,163),
                          pair = c(1,1,2,2),
                          which.reverse = c("read","mate","read","mate"),
                          stringsAsFactors = F)

  wh <- RangedData(space=bed[,1],
                   IRanges(as.numeric(bed[,4])-buffer,
                           as.numeric(bed[,5])+buffer))
  param.hal <- ScanBamParam(which = wh,
                            what=scanBamWhat())

  reads<-readIndexBam(bamfiles = bam.files,
                      index.it = F,
                      chr = bed[,1],
                      start = as.numeric(bed[,4])-buffer,
                      end = as.numeric(bed[,5])+buffer,
                      geneID = bed$geneID)

  hbam = reads[reads$flag %in% good.flags$flag,]
  hbam$start = as.numeric(hbam$start)
  hbam$end = as.numeric(hbam$end)

  if(plotit){
    pdf(plot.pdf,
        height = 4, width = 6)

    if(verbose) cat("\nplotting gene models and reads for", length(genes2test),"genes\n\t... this might take a while ...")
    for(i in genes2test){
      par(mfrow = c(1,2))
      strand = ifelse(gff$V7[grep(i,gff$V9)][1] == "-",-1,1)

      hb = hbam[hbam$geneID == i,]
      hb$dir = ifelse(hb$flag %in% c(99,147),1,-1)*strand

      hb$read = ifelse(hb$flag %in% c(99,83),"r1","r2")
      hb<-ddply(hb, .(geneID, flag), mutate,
                y = (0.05+(rank(start)/length(start)))*dir,
                lwd = ifelse(length(start)<10,3,
                             ifelse(length(start)<50,2,
                                    ifelse(length(start)<200,1,.5))))

      mod = gff[grep(i, gff[,9]),]
      mod = mod[mod$V3 == "CDS",]
      hb$inexon = apply(hb[,c("start","end")],1, function(x){
        has<-apply(mod[,c("V4","V5")],1,function(y){
          (x[1]>y[1] & x[1]<y[2]) |
            (x[2]<y[2] & x[2]>y[1]) |
            x[1]<y[1] & x[2]>y[2]
        })
        return(any(has))
      })
      hb$col = rgb(0,0,0,.5)
      hb$col[hb$inexon]<-rgb(1,0,0,1)

      hb = split(hb,as.factor(hb$read))


      plotGeneModel(gff3 = gff[,c(1:9)], geneID = i)
      with(hb[[1]], segments(x0 = start, x1 = end,
                             y0 = y, y1 = y,
                             lwd = lwd, col = col))
      axis(2, at = c(.5,-.5),
           labels = c("anti-sense reads","sense reads"), tick = F)
      title(paste(i,"read #1"))

      plotGeneModel(gff3 = gff[,c(1:9)], geneID = i)
      with(hb[[2]], segments(x0 = start, x1 = end,
                             y0 = y, y1 = y,
                             lwd = lwd, col = col))
      title(paste(i,"read #2"))
      axis(2, at = c(.5,-.5),
           labels = c("anti-sense reads","sense reads"), tick = F)

    }
    if(verbose) cat("done.\n")
    dev.off()
  }
  return(hbam)
}
