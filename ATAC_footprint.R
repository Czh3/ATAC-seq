## atac-seq footprinting



## ATACseqqc
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg19)

## generate fragement size distribution
a = list.files("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/", "P*.sort.rmdup.bam$", full.names = T)
b = list.files("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/", "P*.sort.rmdup.bam$", full.names = T)



  bamfile = a[3]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize <- fragSizeDist(bamfile, file_name)
  plot(as.numeric(names(fragSize$P4)), as.numeric(fragSize$P4)*1000/sum(as.numeric(fragSize$P4)), xlim=c(0,600), type = "l")
  
  bamfile = a[2]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize1 <- fragSizeDist(bamfile, file_name)
  lines(as.numeric(names(fragSize1$P26)), as.numeric(fragSize1$P26)*1000/sum(as.numeric(fragSize1$P26)), xlim=c(0,600), type = "l", col = "red")
  
  bamfile = a[1]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize2 <- fragSizeDist(bamfile, file_name)
  lines(as.numeric(names(fragSize2$P14)), as.numeric(fragSize2$P14)*1000/sum(as.numeric(fragSize2$P14)), xlim=c(0,600), type = "l", col = "blue")
  
  bamfile = b[3]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize3 <- fragSizeDist(bamfile, file_name)
  plot(as.numeric(names(fragSize3$P5)), as.numeric(fragSize3$P5)*1000/sum(as.numeric(fragSize3$P5)), xlim=c(0,600), type = "l", col = "black")
  
  bamfile = b[2]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize4 <- fragSizeDist(bamfile, file_name)
  lines(as.numeric(names(fragSize4$P36)), as.numeric(fragSize4$P36)*1000/sum(as.numeric(fragSize4$P36)), xlim=c(0,600), type = "l", col = "red")
  
  bamfile = b[1]
  file_name = gsub(".sort.rmdup.bam", "", basename(bamfile))
  fragSize5 <- fragSizeDist(bamfile, file_name)
  lines(as.numeric(names(fragSize5$P26)), as.numeric(fragSize5$P26)*1000/sum(as.numeric(fragSize5$P26)), xlim=c(0,600), type = "l", col = "blue")

  
  
  tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
  ## files will be output into outPath
  outPath <- paste0(dirname(bamfile), "/", file_name, "_split")
  dir.create(outPath)
  ## shift the coordinates of 5'ends of alignments in the bam file

  #seqlev <- "chr1" ## subsample data for quick run
  #which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
  gal <- readBamFile(bamfile, tag=tags, asMates=TRUE)
  gal1 <- shiftGAlignmentsList(gal)
  shiftedBamfile <- file.path(outPath, "shifted.bam")
  export(gal1, shiftedBamfile)


if (0)
{
## foot prints
library(MotifDb)

d <- read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/AP-1/FOSL1.sites.IARs.bed" )
gr <- GRanges(seqnames=Rle(d$V1),
              ranges = IRanges(d$V2, end=d$V3),
              strand = Rle(strand(d$V6)),
              score = d$V5,
              seqlengths = seqlengths(Hsapiens)[unique(d$V1)])
seqlevels(gr) = levels(seqnames(gr))

### IARs
FOSL1 <- query(MotifDb, c("FOSL1"))
FOSL1 <- as.list(FOSL1)
print(FOSL1[[1]], digits=2)
bams = c("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P4_split/shifted.bam",
         "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P14_split/shifted.bam",
         "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P26_split/shifted.bam")

par(cex=1.3)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P4.footprint_FOSL1.pdf",5,5)
my_factorFootprints(bams[1], pfm=FOSL1[[1]], 
                         genome=Hsapiens,
                         bindingSites = gr,
                         min.score="90%", seqlev=seqlevels(gr),
                         upstream=100, downstream=100)
dev.off()
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P14.footprint_FOSL1.pdf",5,5)
my_factorFootprints(bams[2], pfm=FOSL1[[1]], 
                 genome=Hsapiens,
                 bindingSites = gr,
                 min.score="90%", seqlev=seqlevels(gr),
                 upstream=100, downstream=100)
dev.off()
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P26.footprint_FOSL1.pdf",5,5)
my_factorFootprints(bams[3], pfm=FOSL1[[1]], 
                 genome=Hsapiens,
                 bindingSites = gr,
                 min.score="90%", seqlev=seqlevels(gr),
                 upstream=100, downstream=100)
dev.off()


## DARs
d <- read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/ETV1/ERG.sites.DARs.bed")
gr <- GRanges(seqnames=Rle(d$V1),
              ranges = IRanges(d$V2, end=d$V3),
              strand = Rle(strand(d$V6)),
              score = d$V5,
              seqlengths = seqlengths(Hsapiens)[unique(d$V1)])
seqlevels(gr) = levels(seqnames(gr))


ERG <- query(MotifDb, c("ERG"))
ERG <- as.list(ERG)
print(ERG[[1]], digits=2)
bams = c("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P4_split/shifted.bam",
         "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P14_split/shifted.bam",
         "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/align/P26_split/shifted.bam")

par(cex=1.3)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P4.footprint_ERG.pdf",5,5)
my_factorFootprints(bams[1], pfm=ERG[[3]], 
                    genome=Hsapiens,
                    bindingSites = gr,
                    min.score="90%", seqlev=seqlevels(gr),
                    upstream=100, downstream=100)
dev.off()
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P14.footprint_ERG.pdf",5,5)
my_factorFootprints(bams[2], pfm=ERG[[3]], 
                    genome=Hsapiens,
                    bindingSites = gr,
                    min.score="90%", seqlev=seqlevels(gr),
                    upstream=100, downstream=100)
dev.off()
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/P26.footprint_ERG.pdf",5,5)
my_factorFootprints(bams[3], pfm=ERG[[3]], 
                    genome=Hsapiens,
                    bindingSites = gr,
                    min.score="90%", seqlev=seqlevels(gr),
                    upstream=100, downstream=100)
dev.off()




}
  
  
library("Rsamtools")
library("GenomicAlignments")
library("ChIPpeakAnno")
library("ATACseqQC")
  library("motifStack")
my_factorFootprints = function (bamfiles, index = bamfiles, pfm, genome, min.score = "95%", 
            bindingSites, seqlev = paste0("chr", c(1:22, "X", "Y")), 
            upstream = 100, downstream = 100) 
  {
    stopifnot(is(genome, "BSgenome"))
    stopifnot(all(round(colSums(pfm), digits = 4) == 1))
    stopifnot(upstream > 10 && downstream > 10)
    if (missing(bindingSites)) {
      pwm <- motifStack::pfm2pwm(pfm)
      maxS <- maxScore(pwm)
      if (!is.numeric(min.score)) {
        if (!is.character(min.score)) {
          stop("'min.score' must be a single number or string")
        }
        else {
          nc <- nchar(min.score)
          if (substr(min.score, nc, nc) == "%") {
            min.score <- substr(min.score, 1L, nc - 1L)
          }
          else {
            stop("'min.score' can be given as a character string containing a percentage", 
                 "(e.g. '85%') of the highest possible score")
          }
          min.score <- maxS * as.double(min.score)/100
        }
      }
      else {
        min.score <- min.score[1]
      }
      predefined.score <- maxS * as.double(0.85)
      suppressWarnings({
        mt <- matchPWM(pwm, genome, min.score = min(predefined.score, 
                                                    min.score), with.score = TRUE, exclude = names(genome)[!names(genome) %in% 
                                                                                                             seqlev])
      })
      if (min.score <= predefined.score) {
        mt$userdefined <- TRUE
      }
      else {
        mt$userdefined <- FALSE
        mt$userdefined[mt$score >= min.score] <- TRUE
      }
      if (length(mt) > 10000 && sum(mt$userdefined) < 10000) {
        set.seed(seed = 1)
        mt.keep <- seq_along(mt)[!mt$userdefined]
        n <- 10000 - sum(mt$userdefined)
        if (length(mt.keep) > n) {
          mt.keep <- mt.keep[order(mt[mt.keep]$score, decreasing = TRUE)]
          mt.keep <- mt.keep[seq.int(n)]
          mt.keep <- sort(mt.keep)
          mt.keep <- seq_along(mt) %in% mt.keep
          mt <- mt[mt$userdefined | mt.keep]
        }
      }
    }
    else {
      stopifnot(is(bindingSites, "GRanges"))
      stopifnot(length(bindingSites) > 1)
      stopifnot(length(bindingSites$score) == length(bindingSites))
      mt <- bindingSites
      mt$userdefined <- TRUE
    }
    wid <- ncol(pfm)
    seqlevels(mt) <- seqlev
    seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
    which <- as(seqinfo(mt), "GRanges")
    param <- ScanBamParam(which = which)
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, 
                                                     param = param), bamfiles, index, SIMPLIFY = FALSE)
    bamIn <- lapply(bamIn, as, Class = "GRanges")
    if (class(bamIn) != "GRangesList") 
      bamIn <- GRangesList(bamIn)
    bamIn <- unlist(bamIn)
    bamIn <- promoters(bamIn, upstream = 0, downstream = 1)
    libSize <- length(bamIn)
    coverageSize <- sum(as.numeric(width(reduce(bamIn, ignore.strand = TRUE))))
    libFactor <- libSize/coverageSize
    bamIn <- split(bamIn, strand(bamIn))
    cvglist <- sapply(bamIn, coverage)
    cvglist <- cvglist[c("+", "-")]
    cvglist <- lapply(cvglist, function(.ele) .ele[names(.ele) %in% 
                                                     seqlev])
    cvgSum <- cvglist[["+"]] + cvglist[["-"]]
    mt.s <- split(mt, seqnames(mt))
    seqlev <- intersect(names(cvgSum), names(mt.s))
    cvgSum <- cvgSum[seqlev]
    mt.s <- mt.s[seqlev]
    mt.s.ext <- promoters(mt.s, upstream = wid, downstream = wid + 
                            wid)
    stopifnot(all(lengths(mt.s.ext) == lengths(mt.s)))
    mt.v <- Views(cvgSum, mt.s.ext)
    mt.s <- mt.s[viewSums(mt.v) > 0]
    mt <- unlist(mt.s)
    sigs <- featureAlignedSignal(cvglists = cvglist, feature.gr = reCenterPeaks(mt, 
                                                                                width = 1), upstream = upstream + floor(wid/2), downstream = downstream + 
                                   ceiling(wid/2), n.tile = upstream + downstream + wid)
    cor <- lapply(sigs, function(sig) {
      sig.colMeans <- colMeans(sig)
      windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, 
                                step = 1)[[1]]
      windows <- windows[end(windows) <= upstream | start(windows) >= 
                           upstream + wid]
      sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
      windows.sel <- windows[which.max(sig.windowMeans)][1]
      highest.sig.windows <- rowMeans(sig[, start(windows.sel):end(windows.sel)])
      predictedBindingSiteScore <- mt$score
      suppressWarnings({
        cor <- cor.test(x = predictedBindingSiteScore, y = highest.sig.windows, 
                        method = "spearman")
      })
      cor
    })
    sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, 
                                             ])
    mt <- mt[mt$userdefined]
    mt$userdefined <- NULL
    my_plotFootprints(colMeans(do.call(cbind, sigs) * 2/libFactor, 
                            na.rm = TRUE), Mlen = wid, motif = ATACseqQC:::pwm2pfm(pfm))
    return(invisible(list(signal = sigs, spearman.correlation = cor, 
                          bindingSites = mt)))
}

my_plotFootprints = function (Profile, Mlen = 0, xlab = "Dist. to motif (bp)", ylab = "Cut-site probability", 
                              legTitle, newpage = TRUE, motif) 
{
  stopifnot(is(motif, "pfm"))
  if (newpage) 
    grid.newpage()
  S <- length(Profile)
  W <- ((S/2) - Mlen)/2
  vp <- plotViewport(name = "plotRegion")
  pushViewport(vp)
  vp1 <- viewport(y = 0.4, height = 0.8, xscale = c(0, S/2 + 
                                                      1), yscale = c(0, 0.18 * 1.12), name = "footprints")
  pushViewport(vp1)
  grid.lines(x = 1:(S/2), y = Profile[1:(S/2)], default.units = "native", 
             gp = gpar(lwd = 2, col = "darkblue"))
  grid.lines(x = 1:(S/2), y = Profile[(S/2 + 1):S], default.units = "native", 
             gp = gpar(lwd = 2, col = "darkred"))
  grid.xaxis(at = c(seq(1, W, length.out = 3), W + seq(1, Mlen), 
                    W + Mlen + seq(1, W, length.out = 3)), label = c(-(W + 
                                                                         1 - seq(1, W + 1, length.out = 3)), rep("", Mlen), seq(0, 
                                                                                                                                W, len = 3)))
  grid.yaxis()
  grid.lines(x = c(W, W, 0), y = c(0, 0.18, 0.18 * 
                                     1.12), default.units = "native", gp = gpar(lty = 2))
  grid.lines(x = c(W + Mlen + 1, W + Mlen + 1, S/2), y = c(0, 
                                                           0.18, 0.18 * 1.12), default.units = "native", 
             gp = gpar(lty = 2))
  upViewport()
  vp2 <- viewport(y = 0.9, height = 0.2, xscale = c(0, S/2 + 
                                                      1), name = "motif")
  pushViewport(vp2)
  motifStack::plotMotifLogoA(motif)
  upViewport()
  upViewport()
  grid.text(xlab, y = unit(1, "lines"))
  grid.text(ylab, x = unit(1, "line"), rot = 90)
  if (missing(legTitle)) {
    legvp <- viewport(x = unit(1, "npc") - convertX(unit(1, 
                                                         "lines"), unitTo = "npc"), y = unit(1, "npc") - convertY(unit(1, 
                                                                                                                       "lines"), unitTo = "npc"), width = convertX(unit(14, 
                                                                                                                                                                        "lines"), unitTo = "npc"), height = convertY(unit(3, 
                                                                                                                                                                                                                          "lines"), unitTo = "npc"), just = c("right", "top"), 
                      name = "legendWraper")
    pushViewport(legvp)
    grid.legend(labels = c("For. strand", "Rev. strand"), 
                gp = gpar(lwd = 2, lty = 1, col = c("darkblue", "darkred")))
    upViewport()
  }
  else {
    grid.text(legTitle, y = unit(1, "npc") - convertY(unit(1, 
                                                           "lines"), unitTo = "npc"), gp = gpar(cex = 1.2, fontface = "bold"))
  }
  return(invisible())
}
