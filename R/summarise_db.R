#' @export
summarise_db <- function(x){
  loci <- names(x)
  tab_12other <- matrix(data=numeric(), nrow=3, ncol=length(loci))
  x_is_12other <- matrix(data=character(), nrow=nrow(x), ncol=length(loci))
  
  rownames(tab_12other) <- c("fraction 1 integer allele", "fraction 2 integer alleles", "fraction other")
  colnames(x_is_12other) <- colnames(tab_12other) <- loci
  
  x_unpacked <- list()
  for(locus in loci){
    x_locus <- x[[locus]]
    
    # split by comma
    x_locus_unpacked <- unpack_haplotypes(x_locus)
    x_unpacked[[locus]] <- x_locus_unpacked
    
    # grab only integer alleles
    x_locus_int <- filter_integer_alleles(x_locus_unpacked)
    
    ind_one_int <- sapply(x_locus_int,length)==1 & sapply(x_locus_unpacked,length)==1
    ind_two_int <- sapply(x_locus_int,length)==2 & sapply(x_locus_unpacked,length)==2
    ind_other <- (!ind_one_int) & (!ind_two_int)
    
    tab_12other[1,locus] <- mean(ind_one_int)
    tab_12other[2,locus] <- mean(ind_two_int)
    tab_12other[3,locus] <- mean(ind_other)
    
    x_is_12other[ind_one_int,locus] <- "1"
    x_is_12other[ind_two_int,locus] <- "2"
    x_is_12other[ind_other,locus] <- "other"
  }
 
  # determine suitable 1-loci and 2-loci
  one_loci <- loci[tab_12other[1,] > 0.8]
  two_loci <- loci[tab_12other[2,] > 0.8]
  
  list(x_ind_12_other = x_is_12other, one_loci = one_loci, two_loci = two_loci, tab = tab_12other)
}
