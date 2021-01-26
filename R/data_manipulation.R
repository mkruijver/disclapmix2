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
