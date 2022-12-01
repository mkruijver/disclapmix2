#' @title Compute Profile Probability from fit
#'
#' @param x DataFrame. Columns should be one character vector for each locus
#' @param fit	Output from disclapmix2
#' @description Compute the profile probability for a new profile that was not used in the original fit.
#' @returns Numeric.
#' @examples
#' require(disclapmix)
#'
#' data(danes) 
#'
#' x <- as.data.frame(sapply(danes[rep(seq_len(nrow(danes)), danes$n), -ncol(danes)], as.character))
#'
#' dlm2_fit <- disclapmix2(x, number_of_clusters = 3)
#' 
#' 
#' new_profile <- structure(list(DYS19 = "14", DYS389I = "13", DYS389II = "29", 
#'                               DYS390 = "22", DYS391 = "9", DYS392 = "15", DYS393 = "13", 
#'                               DYS437 = "14", DYS438 = "11", DYS439 = "12"), row.names = 1L, class = #' "data.frame")
#' 
#' profile_pr_by_locus_and_cluster(x = new_profile, dlm2_fit)
#' @export
profile_pr_by_locus_and_cluster <- function(x, fit){
  
  ## check x
  if (!is.data.frame(x)){
    stop("x needs to be a data frame")
  }
  
  if (!all(sapply(x, is.character))){
    stop("columns of x need to be of class character")
  }
  
  ## check loci
  fit_loci <- c(fit$one_loci, fit$two_loci)
  
  if (length(x) != length(fit_loci)){
    stop("x has ", length(x), " loci while fit has ", length(fit_loci), " loci")
  }
  
  if (!all(fit_loci %in% names(x))){
    stop("not all loci in fit are present in x")
  }
  
  # reorder x loci to be in fit order
  x_ordered <- x[fit_loci]
  
  # determine which haplotypes are 1, 2 or non-standard
  x_summary <- summarise_db(x_ordered)
  
  # get array indices of non-standard haplotypes
  one_two <- matrix(c(rep("1", length(fit$one_loci)), 
                      rep("2", length(fit$two_loci))), 
                    nrow = nrow(x), ncol = length(fit_loci), byrow = TRUE)
  
  ns_idx <- which(x_summary$x_ind_12_other != one_two, arr.ind = TRUE)
  
  x_ordered_character_filtered <- x_ordered
  x_ordered_character_filtered[ns_idx] <- NA
  
  # convert to integer matrix
  x_int <- to_int_db(x_ordered_character_filtered, 
                                   fit$one_loci, fit$two_loci)
  
  
  number_of_1_loci <- length(fit$one_loci)
  number_of_2_loci <- length(fit$two_loci)
  
  pr <- compute_profiles_pr_by_cluster_and_locus(fit$p, x_int, fit$y, number_of_1_loci, number_of_2_loci)
  
  ## now deal with the non-standard haplotypes
  ## which so far have pr. 1
  
  # first we set any non-standard haplotypes to pr. 0
  for (i in seq_len(nrow(ns_idx))){
    pr[[ns_idx[i, 1]]][, ns_idx[i, 2]] <- 0.
  }
  
  # put a factor of 1 - sum(q) in the pr.
  if (!is.null(fit$pi)){
    for (i_profile in seq_along(pr)){
      pr[[i_profile]] <- pr[[i_profile]] * (1 - fit$pi)
    }
  }
  
  # set pr to q for all non-standard haplotypes
  for (i in seq_len(nrow(ns_idx))){
    
    i_profile <- ns_idx[i, 1]
    locus <- fit_loci[ns_idx[i, 2]]
    g <- x_ordered[i_profile, ns_idx[i, 2]]
    
    freq_row <- match(paste0(locus, "_", g), rownames(fit$q))
    
    if (!is.na(freq_row)){
      pr[[i_profile]][, ns_idx[i,2]] <- fit$q[freq_row,]
    }
  }
  
  if (nrow(x) == 1){
    return(pr[[1]])
  }
  else{
    return(pr)
  }
}

