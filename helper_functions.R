PeakLinkageCor <- function(deg, object, nn) {
  
  peak.ranges <- granges(object)
  
  pro.ranges <- GetPromoterRange(object, deg)
  
  if (is.null(pro.ranges)) {
    
    return(NULL)
    
  }
  
  pro.peaks <- FindProximalPeaks(peak.ranges, pro.ranges)
  
  DefaultAssay(object) <- "ATAC"
  
  pro.peaks.name <- rownames(object)[pro.peaks$idx]
  
  if (length(pro.peaks.name) == 0){
    
    return(NULL)
    
  }
  
  
  
  gene.expr <- object@assays$RNA@counts[deg, ]
  
  gene.smooth <- SmoothData(gene.expr, nn)
  
  
  
  peak.mat <- object@assays$ATAC@counts[pro.peaks.name, , drop=FALSE]
  
  peak.smooth.mat <- apply(peak.mat, 1, SmoothData, nn=nn)
  
  return(cor(peak.smooth.mat, gene.smooth))
  
}
GetPromoterRange <- function(object, gene, upstream = 5e4) {
  
  exons <- object@assays$ATAC@annotation %>% as_data_frame %>% 
    
    filter(type == "exon", gene_name == gene)
  
  exon1 <- head(exons, n = 1)
  
  exonn <- tail(exons, n = 1)
  
  if (nrow(exon1) == 0) {
    
    return (NULL)
    
  }
  
  return(data.frame(seqnames = exon1$seqnames, 
                    
                    strand = exonn$strand, start = exon1$start - upstream, end = exonn$end + upstream))
  
}
FindProximalPeaks <- function(peak.ranges, promoter.range) {
  
  idx <- which(seqnames(peak.ranges) == as.character(promoter.range$seqnames) & 
                 
                 start(peak.ranges) > promoter.range$start & 
                 
                 start(peak.ranges) < promoter.range$end)
  
  return(list(idx = idx, ranges = peak.ranges[idx, ]))
  
}
SmoothData <- function(data, nn) {
  
  smooth.data <- sapply(1:length(data), function(i) {
    
    nn.idx <- nn$nn.idx[i, ]
    
    nn.dist <- nn$nn.dist[i, ]
    
    weights <- compute_weight(nn.dist)
    
    return(sum(weights * data[nn.idx]))
    
  })
  
  return (as.numeric(smooth.data))
  
}
compute_weight <- function(distances, sigma.idx = 10) {
  
  sigma <- distances[sigma.idx]
  
  if (sigma == 0)
    
    sigma <- 1
  
  w <- exp(-0.5* (distances^2) / (sigma^2))
  
  w.sum <- sum(w)
  
  if (w.sum == 0) 
    
    w.sum <- 1
  
  return(w/w.sum)
  
}
TFGeneCorrelation <- function(tf, peaks, object, nn) {
  
  peak_values <- object@assays$ATAC@data[peaks, , drop=FALSE]
  
  peaks.smooth <- apply(peak_values, 1, SmoothData, nn)
  
  tf.values <- as.numeric(object@assays$SCT@data[tf, ])
  
  tf.smooth <- SmoothData(tf.values, nn)
  
  tf.cor <- cor(peaks.smooth, tf_values)
  
  return(mean(tf.cor))
  
}