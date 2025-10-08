filter_integer_alleles <- function(haplotypes_unpacked){
  lapply(haplotypes_unpacked, function(h){
    h_int <- as.integer(h[as.character(suppressWarnings(as.integer(h)))==h])
    h_int <- h_int[h_int>0]
  })
}

pack_haplotypes <- function(haplotypes_unpacked){
  sapply(haplotypes_unpacked, function(x) if((length(x)>0) && (!all(is.na(x)))) paste0(x, collapse = ",") else NA_character_)
}

unpack_haplotypes <- function(haplotypes){
  strsplit(haplotypes,split = ",")
}

clean_and_reorder <- function(x, x_summarised, include_2_loci, verbose = 0){
  
  one_loci <- x_summarised$one_loci
  two_loci <- if (include_2_loci) x_summarised$two_loci else character()
  
  loci <- c(one_loci, two_loci)
  
  # prepare output matrix
  x_cleaned <- matrix(data = character(), nrow = nrow(x), ncol = length(loci))
  colnames(x_cleaned) <- loci
  
  if (verbose>=1) verbose_print("Setting any haplotype that is not 1 or 2 integers at the approriate locus to NA")
  removed_non_standard_dfs <- list()
  
  for(locus in one_loci){
    ind_1 <- which(x_summarised$x_ind_12_other[,locus] == "1")
    not_ind_1 <- which(sapply(x_summarised$x_ind_12_other[,locus] != "1", isTRUE))
    
    x_cleaned[ind_1, locus] <- x[[locus]][ind_1]
    
    number_removed <- length(not_ind_1)
    removed_non_standard_dfs[[1 + length(removed_non_standard_dfs)]] <- 
      data.frame(row = not_ind_1, 
                 locus = rep(locus, number_removed), 
                 haplotype = x[[locus]][not_ind_1])
  }
  
  for(locus in two_loci){
    ind_2 <- which(x_summarised$x_ind_12_other[,locus] == "2")
    not_ind_2 <- which(sapply(x_summarised$x_ind_12_other[, locus] != "2", isTRUE))
    
    x_cleaned[ind_2, locus] <- x[[locus]][ind_2]
    
    number_removed <- length(not_ind_2)
    
    removed_non_standard_dfs[[1+length(removed_non_standard_dfs)]] <- 
      data.frame(row = not_ind_2, 
                 locus = rep(locus, number_removed), 
                 haplotype = x[[locus]][not_ind_2])
  }

  # list out each removed non-standard haplotype
  removed_non_standard_df <- do.call(rbind, removed_non_standard_dfs)
  
  # list out the non-standard haplotypes
  non_standard_df <- unique(removed_non_standard_df[, c("locus", "haplotype")])
  non_standard_df$locus <- factor(non_standard_df$locus, levels = loci)
  rownames(non_standard_df) <- NULL
  
  # assign index
  non_standard_df$index <- -seq_len(nrow(non_standard_df))
  
  if (nrow(non_standard_df) > 0){
    removed_non_standard_df$index <- non_standard_df$index[match(
      paste0(removed_non_standard_df$locus, "___", removed_non_standard_df$haplotype),
      paste0(non_standard_df$locus, "___", non_standard_df$haplotype))]
  }
  
  list(
    x_before_cleaning = x,
    x_cleaned = x_cleaned,
    removed_non_standard_df = removed_non_standard_df,
    non_standard_df = non_standard_df
  )
  
}