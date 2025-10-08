summarise_db <- function(x){
  loci <- names(x)
  tab_12other <- matrix(data=numeric(), nrow=5, ncol=length(loci))
  x_is_12other <- matrix(data=character(), nrow=nrow(x), ncol=length(loci))
  
  rownames(tab_12other) <- c("fraction 1 integer allele", "fraction 2 integer alleles", 
                             "fraction other", "fraction NA", "# not NA")
  colnames(x_is_12other) <- colnames(tab_12other) <- loci
  
  x_unpacked <- list()

  for(locus in loci){
    x_locus <- x[[locus]]
    
    # split by comma
    x_locus_unpacked <- unpack_haplotypes(x_locus)

    # grab only integer alleles
    x_locus_not_NA <- which(!is.na(x_locus))
    x_locus_int <- rep(NA_integer_, length(x_locus))
    x_locus_int[x_locus_not_NA] <- filter_integer_alleles(x_locus_unpacked[x_locus_not_NA])
    
    ind_na <- is.na(x_locus)
    ind_one_int <- !ind_na & 
      sapply(x_locus_int,length)==1 & sapply(x_locus_unpacked,length)==1 &
      (!sapply(x_locus_int, anyNA))
    ind_two_int <- !ind_na & 
      sapply(x_locus_int,length)==2 & sapply(x_locus_unpacked,length)==2 & 
      (!sapply(x_locus_int, anyNA))
    ind_other <- (!ind_one_int) & (!ind_two_int) * (!ind_na)
    
    tab_12other[1,locus] <- mean(ind_one_int)
    tab_12other[2,locus] <- mean(ind_two_int)
    tab_12other[3,locus] <- mean(ind_other)
    tab_12other[4,locus] <- mean(ind_na)
    
    x_is_12other[ind_one_int,locus] <- "1"
    x_is_12other[ind_two_int,locus] <- "2"
    x_is_12other[ind_na,locus] <- NA
    x_is_12other[ind_other,locus] <- "other"
  }
  
  tab_12other[5,] <- sapply(x, function(x0) sum(!is.na(x0)))
  
  # determine suitable 1-loci and 2-loci
  one_loci <- loci[sapply(tab_12other["fraction 1 integer allele", ]/
                            (1 - tab_12other["fraction NA",]) >= 0.5, isTRUE)]
  two_loci <- loci[sapply(tab_12other["fraction 2 integer alleles", ]/
                            (1 - tab_12other["fraction NA",]) >= 0.5, isTRUE)]
  
  list(x_ind_12_other = x_is_12other, one_loci = one_loci, two_loci = two_loci, tab = tab_12other)
}
