getCorr <- function(SE_RNA,SE_ATAC,max.dist= 250000){
  o <- data.frame(findOverlaps(resize(SE_RNA, 2 * max.dist + 
                                        1, "center"), resize(rowRanges(SE_ATAC_filt), 1, "center"), 
                               ignore.strand = TRUE))
  
  o$distance <- IRanges::distance(rowRanges(SE_RNA)[o[, 1]], 
                                  rowRanges(SE_ATAC_filt)[o[, 2]])
  
  colnames(o) <- c("gene_idx", "peak_idx", "distance")
  nullCor <- .getNullCorrelations(SE_ATAC_filt, SE_RNA, o, 1000)
  
  df <- rowRanges(SE_ATAC_filt)[o$peak_idx, ]
  o$gene <- rowData(SE_RNA)[o$gene_idx, ]$SYMBOL
  o$peak <- paste0(df@seqnames, "-", as.data.frame(df@ranges)$start, 
                   "-", as.data.frame(df@ranges)$end)
  o$Correlation <- rowCorr:::rowCorCpp(as.integer(o$peak_idx), as.integer(o$gene_idx), 
                                       assay(SE_ATAC_filt), assay(SE_RNA))
  o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
                                       1e-17, na.rm = TRUE))/(ncol(SE_ATAC_filt) - 2)))
  o$Pval <- 2 * pt(-abs(o$TStat), ncol(SE_ATAC_filt) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$EmpPval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  o$EmpFDR <- p.adjust(o$EmpPval, method = "fdr")
  o <- o[!is.na(o$FDR), ]
  return(o)
}
.getNullCorrelations <- function(seA, seB, o, n){

  o$seq <- seqnames(seA)[o$peak_idx]

  nullCor <- lapply(seq_along(unique(as.character(o$seq))), function(i){

    #Get chr from olist
    chri <- unique(as.character(o$seq))[i]
    #message(chri, " ", appendLF = FALSE)

    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[as.character(o$seq)==chri,]$gene_idx))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
    out <- rowCorr::rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)

    return(out)

  }) %>% SimpleList
  #message("")

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

}

