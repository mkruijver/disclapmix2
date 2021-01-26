#' Count the number of times each haplotype occurs
#'
#' @param x DataFrame (by locus) of character vectors containing haplotypes (rows) where alleles are separated by comma's, e.g. "13,14.2" is a haplotype
#' @return Integer vector 
#' @export
haplotype_counts <- function(x){
  check_input_db(x)

  profiles_char <- apply(as.matrix(x),1,paste0,collapse="|||")
  
  as.integer(table(profiles_char)[profiles_char])
}

x_raw <- list()
for(i in 1:3){
  x_raw_i <- readxl::read_excel(system.file("extdata", "South_Australia.xlsx", package="disclapmix2"), col_types = "text", sheet = i)
  x_raw[[i]] <- x_raw_i 
}

x_raw_bound <- do.call(rbind, x_raw)
x <- x_raw_bound[-c(1,2)] # remove 'Sample Name' and 'Population' columns


#' List unique haplotypes with their counts
#'
#' @param x DataFrame (by locus) of character vectors containing haplotypes (rows) where alleles are separated by comma's, e.g. "13,14.2" is a haplotype
#' @return DataFrame
#' @export
haplotype_table <- function(x){
  check_input_db(x)
  
  profiles_char <- apply(as.matrix(x),1,paste0,collapse="|||")
  unique_profiles_char <- unique(profiles_char)
  idx <- match(unique_profiles_char, profiles_char)
  
  x_out <- x[idx,]
  x_out$Count <- haplotype_counts(x)[idx]
  
  x_out <- x_out[order(x_out$Count, decreasing = TRUE),]
  
  x_out
}
