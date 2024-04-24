.getNullCorrelations <- function(seA, seB, o, n){
  
  o$seq <- seqnames(seA)[o$A]
  
  nullCor <- lapply(seq_along(unique(o$seq)), function(i){
    
    #Get chr from olist
    chri <- unique(o$seq)[i]
    #message(chri, " ", appendLF = FALSE)
    
    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }
    
    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))
    
    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])
    
    seSubA <- seA[idxA]
    seSubB <- seB[idxB]
    
    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)
    
    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)
    
    return(out)
    
  }) %>% SimpleList
  #message("")
  
  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)
  
  return(list(summaryDF, unlist(nullCor)))
  
}

