filter_integer_alleles <- function(haplotypes_unpacked){
  lapply(haplotypes_unpacked, function(h){
    h_int <- as.integer(h[as.character(as.integer(h))==h])
    h_int[h_int>0] # do not include the NULL allele
  })
}

# h <- haplotypes_unpacked[[46]]
# h
filter_non_integer_alleles <- function(haplotypes_unpacked, as_numeric = FALSE){
  lapply(haplotypes_unpacked, function(h){
    h_filtered <- h[as.character(as.integer(h))!=h] # also does not include the NULL allele
    if (as_numeric) as.numeric(h_filtered) else h_filtered
  })
}

pack_haplotypes <- function(haplotypes_unpacked){
  sapply(haplotypes_unpacked, function(x) if((length(x)>0) && (!all(is.na(x)))) paste0(x, collapse = ",") else NA_character_)
}

unpack_haplotypes <- function(haplotypes){
  strsplit(haplotypes,split = ",")
}

haplotypes_unpacked_to_df <- function(haplotypes_unpacked, min_number_of_columns = 1L, colname_prefix = "allele", factor_levels){
  number_of_columns <- max(c(min_number_of_columns,sapply(haplotypes_unpacked, length)))
  
  by_column <- lapply(seq_len(number_of_columns), function(i_column) sapply(haplotypes_unpacked, get("["),i_column))
  names(by_column) <- paste0(colname_prefix, seq_len(number_of_columns))
  
  if (!missing(factor_levels)){
    by_column <- lapply(by_column, factor, levels=factor_levels)
  }
  
  data.frame(by_column, stringsAsFactors = FALSE, check.names = FALSE)
}

sort_alleles <- function(alleles){
  if (!is.character(alleles)) stop("expected a character vector")
  if (anyNA(alleles)) stop("NAs encountered among alleles")
  
  alleles_numeric <- suppressWarnings(as.numeric(alleles))
  
  c(sort(alleles_numeric[!is.na(alleles_numeric)]), sort(alleles[is.na(alleles_numeric)]))
}

unpack_profiles <- function(profiles){
  if (!inherits(profiles, "data.frame")) stop("expected a data.frame")
  
  r <- setNames(lapply(profiles, disclapmixExtended:::unpack_haplotypes), nm = names(profiles))
  class(r) <- "data.frame"
  rownames(r) <- rownames(profiles)
  r
}

pack_profiles <- function(profiles){
  if (!inherits(profiles, "data.frame")) stop("expected a data.frame")
  
  packed <- lapply(profiles, function(profiles_locus){
    sapply(profiles_locus, pack_haplotypes)
  })
  
  data.frame(setNames(packed, nm = names(profiles)), stringsAsFactors = FALSE, check.names = FALSE)
}

validate_packed_profiles <- function(profiles){
  if (!inherits(profiles, "data.frame")) stop("expected a data.frame")
  if (!all(sapply(profiles, is.character))) stop("columns of profiles need to be character vectors")
  
  for(locus in names(profiles)){
    idx <- grep(pattern = "[^0-9\\.\\,]", profiles[[locus]])
    if(length(idx)>0) stop("Unexpected input at locus ", locus, ": ", profiles[[locus]][idx])
  }
}

validate_unpacked_profiles <- function(profiles_unpacked, require_integer_alleles = FALSE){
  if (!is.list(profiles_unpacked)) stop("unpacked profiles need to be a list")
  if (!all(sapply(profiles_unpacked, is.list))) stop("unpacked profiles need to be a list of lists")
  
  if (!all(sapply(profiles_unpacked, function(profiles_locus){
    all(sapply(profiles_locus, is.character))
  }))) stop("unpacked profiles need to be a list (by locus) of lists (by profile) of character vectors (alleles)")
  
  if (require_integer_alleles){
    for(locus in names(profiles_unpacked)){
      idx <- grep(pattern = "[^0-9]", profiles_unpacked[[locus]])
      if(length(idx)>0) stop("Non-integer input at locus ", locus, ": ", profiles_unpacked[[locus]][idx])
    }
  }
}

profiles_filter_1_integer_allele <- function(profiles, known_haplotype_summary){
  loci <- names(profiles)
  
  filtered <- list()
  for(locus in loci){
    haplotypes <- profiles[[locus]]
    # convert factor to integers
    haplotype_idx <- as.integer(haplotypes)
    
    haplotype_has_1_integer_allele <- (known_haplotype_summary[[locus]]$number_of_integer_alleles==1)[haplotype_idx]
    haplotypes[!haplotype_has_1_integer_allele] <- NA
    filtered[[locus]] <- haplotypes
  } 
  data.frame(filtered, check.names = FALSE)
}
