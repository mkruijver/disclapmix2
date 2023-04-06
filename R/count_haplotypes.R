#' Count the number of times each haplotype occurs
#'
#' @param x DataFrame (by locus) of character vectors containing haplotypes (rows) where alleles are separated by comma's, e.g. "13,14.2" is a haplotype
#' @return Integer vector with count for each row in DataFrame
#' @export
#' @examples 
#' # read haplotypes
#' h <- readxl::read_excel(system.file("extdata","South_Australia.xlsx",
#' package = "disclapmix2"), 
#' col_types = "text")[-c(1,2)]
#' 
#' # obtain counts
#' counts <- disclapmix2::haplotype_counts(h)
#' 
#' # all haplotypes in the dataset are unique
#' stopifnot(all(counts == 1))
haplotype_counts <- function(x){
  check_input_db(x)

  if (any(sapply(x, anyNA))){
    stop("Haplotype counts cannot be determined because an NA value is present")
  }
  
  profiles_char <- apply(as.matrix(x),1,paste0,collapse="|||")
  
  as.integer(table(profiles_char)[profiles_char])
}
NULL
#' List unique haplotypes with their counts
#'
#' @param x DataFrame (by locus) of character vectors containing haplotypes (rows) where alleles are separated by comma's, e.g. "13,14.2" is a haplotype
#' @return DataFrame with unique rows and a Count column added at the end
#' @export
#' @examples 
#' # read haplotypes
#' h <- readxl::read_excel(system.file("extdata","South_Australia.xlsx",
#' package = "disclapmix2"), 
#' col_types = "text")[-c(1,2)]
#' 
#' # obtain counts
#' unique_counts <- disclapmix2::unique_haplotype_counts(h)
#' 
#' # all haplotypes in the dataset are unique
#' stopifnot(all(unique_counts$Count == 1))
unique_haplotype_counts <- function(x){
  check_input_db(x)
  
  if (any(sapply(x, anyNA))){
    stop("Haplotype counts cannot be determined because an NA value is present")
  }
  
  profiles_char <- apply(as.matrix(x),1,paste0,collapse="|||")
  unique_profiles_char <- unique(profiles_char)
  idx <- match(unique_profiles_char, profiles_char)
  
  x_out <- x[idx,]
  x_out$Count <- haplotype_counts(x)[idx]
  
  x_out <- x_out[order(x_out$Count, decreasing = TRUE),]
  
  x_out
}
