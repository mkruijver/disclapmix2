#' @title Discrete Laplace mixture inference using Numerical Optimisation
#'
#' @param x DataFrame. Columns should be one character vector for each locus
#' @param number_of_clusters	The number of clusters to fit the model for.
#' @param include_2_loci Should duplicated loci be included or excluded from the analysis? 
#' @param remove_non_standard_haplotypes Should observations that are not single integer alleles be removed?
#' @param use_stripped_data_for_initial_clustering Should non_standard data be removed for the initial clustering?
#' @param initial_y_method Which cluster method to use for finding initial central haplotypes, y: pam (recommended) or clara.
#' @param verbose Set to 1 (or higher) to print optimisation details. Default is 0.
#' @description An extension to the *disclapmix* method in the *disclapmix* package that supports duplicated loci and other non-standard haplotypes.
#' @returns List.
#' @examples
#' require(disclapmix)
#'
#' data(danes) 
#'
#' x <- as.matrix(danes[rep(seq_len(nrow(danes)), danes$n), -ncol(danes)])
#' x2 <- as.data.frame(sapply(danes[rep(seq_len(nrow(danes)), danes$n), -ncol(danes)], as.character))
#'
#'
#' dlm_fit <- disclapmix(x, clusters = 3L)
#' dlm2_fit <- disclapmix2(x2, number_of_clusters = 3)
#'
#' stopifnot(all.equal(dlm_fit$logL_marginal, dlm2_fit$log_lik))
#' @export
disclapmix2 <- function(x, number_of_clusters, include_2_loci = FALSE, remove_non_standard_haplotypes = TRUE, 
                        use_stripped_data_for_initial_clustering = FALSE, initial_y_method = "pam",
                        verbose=0L){
  
  check_input_db(x)
  if (!(isTRUE(initial_y_method %in% c("pam", "clara")))){
    stop("initial_y_method needs to be \"pam\" or \"clara\"")
  }
  
  if (verbose>=1) verbose_print("Determining loci with suitable data")
  x_summarised <- summarise_db(x)
  
  # clean the data such that 1-loci only have 1 integer or NA
  # and 2-loci only have 2 integers or NA
  number_of_1_loci <- length(x_summarised$one_loci)
  number_of_2_loci <- if (include_2_loci) length(x_summarised$two_loci) else 0

  one_loci <- x_summarised$one_loci
  two_loci <- if (include_2_loci) x_summarised$two_loci else character()
  loci <- c(one_loci, two_loci)
  number_of_loci <- length(loci)
  
  locuslabel <- function(i) if (i==1) "locus" else "loci"
  if (verbose>=1) verbose_print("Found", number_of_1_loci, paste0("suitable 1-", locuslabel(number_of_1_loci)),
                                              "and", number_of_2_loci,paste0("suitable 2-", locuslabel(number_of_2_loci)),
                                              "among", ncol(x), "columns")
  
  x_cleaned <- matrix(data = character(), nrow = nrow(x), ncol = number_of_loci)
  colnames(x_cleaned) <- c(one_loci, if(number_of_2_loci>0) two_loci else NULL)
  
  if (verbose>=1) verbose_print("Setting any haplotype that is not 1 or 2 integers at the approriate locus to NA")
  removed_non_standard_dfs <- list()
  
  # locus = "DYS458"
  for(locus in one_loci){
    ind_1 <- which(x_summarised$x_ind_12_other[,locus]=="1")
    not_ind_1 <- which(sapply(x_summarised$x_ind_12_other[,locus]!="1", isTRUE))
    
    x_cleaned[ind_1,locus] <- x[[locus]][ind_1]
    
    number_removed <- length(not_ind_1)
    removed_non_standard_dfs[[1+length(removed_non_standard_dfs)]] <- data.frame(row = not_ind_1, locus=rep(locus, number_removed), haplotype=x[[locus]][not_ind_1])
  }
  
  if(number_of_2_loci > 0){
    for(locus in two_loci){
      ind_2 <- which(x_summarised$x_ind_12_other[,locus]=="2")
      not_ind_2 <- which(sapply(x_summarised$x_ind_12_other[,locus]!="2", isTRUE))

      x_cleaned[ind_2,locus] <- x[[locus]][ind_2]
      
      number_removed <- length(not_ind_2)
      
      removed_non_standard_dfs[[1+length(removed_non_standard_dfs)]] <- data.frame(row = not_ind_2, locus=rep(locus, number_removed), haplotype=x[[locus]][not_ind_2])
    }
  }
  
  # list out each removed non-standard haplotype
  removed_non_standard_df <- do.call(rbind, removed_non_standard_dfs)
  
  # list out the non-standard haplotypes
  non_standard_df <- unique(removed_non_standard_df[,c("locus","haplotype")])
  non_standard_df$locus <- factor(non_standard_df$locus, levels=loci)
  rownames(non_standard_df) <- NULL
  
  # assign index
  non_standard_df$index <- -seq_len(nrow(non_standard_df))
  if (nrow(non_standard_df) > 0){
    removed_non_standard_df$index <- non_standard_df$index[match(
      paste0(removed_non_standard_df$locus, "___", removed_non_standard_df$haplotype),
      paste0(non_standard_df$locus, "___", non_standard_df$haplotype))]
  }
  
  if (verbose>=1){
    if (nrow(removed_non_standard_df>1)){
      verbose_print("Set", nrow(removed_non_standard_df), "non-standard haplotypes to NA:")
      print(removed_non_standard_df)
    }
    else{
      verbose_print("All haplotypes are integer valued")
    }
  }
  
  # create a stripped down version of the db with only the easy profiles
  x_stripped <- x_cleaned[rowSums(is.na(x_cleaned))==0,]
  x_stripped_int <- to_int_db(x_stripped, one_loci, two_loci)
  
  if (remove_non_standard_haplotypes){
    if (verbose>=1) verbose_print("Removing profiles with any haplotype that is not 1 or 2 integers at the approriate locus")

    x_int <- x_stripped_int
  }else{
    if (verbose>=1) verbose_print("Not removing profiles with any haplotype that is not 1 or 2 integers at the approriate locus")
    x_int <- to_int_db(x_cleaned, one_loci, two_loci)
    
    observations_columns <- c(one_loci, if(number_of_2_loci > 0) paste0(two_loci,"(1)") else character())
    
    if (verbose>=1) verbose_print(sum(is.na(x_int[,observations_columns])), "NA observations in x_int")
  }
  
  if (verbose>=1) verbose_print(nrow(x_int) , paste0("profiles used in analysis (started with ", nrow(x),")"))
  if (verbose>=1) verbose_print("Determining initial clustering using", initial_y_method)
  
  if (!use_stripped_data_for_initial_clustering){
    
    if (initial_y_method=="pam"){
      initial_clustering <- cluster::pam(x_int, k = number_of_clusters, metric = "manhattan",diss = FALSE, keep.diss=FALSE, keep.data=TRUE)
    }
    else{
      initial_clustering <- cluster::clara(x_int, k = number_of_clusters, metric = "manhattan",pamLike = TRUE,samples=100, correct.d = TRUE)
    }
    
    if (anyNA(initial_clustering$medoids)){
      use_stripped_data_for_initial_clustering <- TRUE
      if (verbose>=1) verbose_print("NAs in initial clustering, retrying with stripped data")    }
  }else{
    if (verbose>=1) verbose_print("Using stripped data for initial clustering")
  }
  
  if (use_stripped_data_for_initial_clustering){
    
    if (initial_y_method == "pam"){
      initial_clustering <- cluster::pam(x_stripped_int, k = number_of_clusters, 
                                         metric = "manhattan",diss = FALSE, keep.diss=FALSE, keep.data=TRUE)
    }
    else{
      initial_clustering <- cluster::clara(x_stripped_int, k = number_of_clusters,
                                           metric = "manhattan", pamLike = TRUE, samples=100, correct.d = TRUE)
      
    }
    
  }
  
  initial_pam_tau <- (tabulate(initial_clustering$clustering, nbins = number_of_clusters)/nrow(x_int))
  initial_pam_theta_tau <- if(number_of_clusters==1L) numeric() else log( initial_pam_tau[1:(number_of_clusters-1)])
  
  # now find the variances
  y <- y0 <- initial_clustering$medoids
  theta <- theta0 <- if(number_of_clusters>1) c(initial_pam_theta_tau, rep(-1,number_of_clusters), rep(0,number_of_loci-1)) else rep(-1,number_of_loci)
  theta <- theta0 <- if(number_of_clusters>1) c(initial_pam_theta_tau, rep(-0.5,number_of_clusters), rep(-0.5,number_of_loci-1)) else rep(-1,number_of_loci)
  
  if (verbose>=1) verbose_print("Starting optimisation")
  
  cluster_labels <- paste0("cluster", seq_len(number_of_clusters))
  
  y_iterations <- list()
  theta_iterations <- list()
  
  optim_trace <- if (verbose>0) verbose - 1 else 0
  
  repeat{
    rownames(y) <- cluster_labels
    y_iterations[[1 + length(y_iterations)]] <- y
    
    f <- function(theta) {
      negll <- disclapmix2::neg_loglik_theta(theta, x_int, y, number_of_1_loci, number_of_2_loci)
      
      negll
    }
    
    opt <- stats::optim(par = theta, fn = f, method = "BFGS", control = list(maxit=500, reltol = 1e-9, trace=optim_trace))
    if (opt$convergence != 0) stop("BFGS failed to converge")

    theta_opt <- opt$par
    theta_iterations[[1+length(theta_iterations)]] <- theta_opt
    
    tau_opt <- get_tau(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
    p_opt <- get_P(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
    
    rownames(p_opt) <- cluster_labels
    colnames(p_opt) <- loci
    
    x_profile_pr_by_cluster <- compute_profile_prs(p_by_cluster_and_locus = p_opt, db = x_int, y = y, number_of_1_loci, number_of_2_loci)
    v_matrix <- compute_posterior_cluster_prs(profile_pr = x_profile_pr_by_cluster, tau = tau_opt)
    
    # see if we need to move centers
    y_new <- move_centers(x_int, y, v_matrix)
    
    theta <- theta_opt
    
    dist_new_y <- sum(abs(y - y_new))
    if (dist_new_y == 0) {
      if (verbose >= 1L) {
        verbose_print("Current central haplotypes are optimal")
      }
      break
    } else if (any(duplicated(y_new))) { # new case introduced in disclapmix version 1.6.3
      if (verbose >= 1L) {
        verbose_print("New central haplotypes had at least one duplicated haplotype, change rejected")
      }
      break
    } else {
      if (verbose >= 2L) {      
        verbose_print("Current central haplotypes are not otimal")
        print(y)
        verbose_print("New central haplotypes:")
        print(y_new)
        verbose_print("Differences:")
        print(y_new - y)
        verbose_print("Number of stepwise mutations between center configurations = ", dist_new_y)
        
      } else if (verbose >= 1) {
        verbose_print("Current central haplotypes are not optimal, moving and restarting optimisation")
      }
      
      y <- y_new
    }
  }
  
  if((!remove_non_standard_haplotypes) & (nrow(removed_non_standard_df)>0)){
    
    repeat{
      verbose_print(nrow(removed_non_standard_df), "observations of",
                                  nrow(non_standard_df),
                                  "non-standard haplotypes were ignored, second stage estimation begins")
      
      # relabel the matrix such that non-standard alleles get negative indices
      x_int_ns <- x_int
      for(i in seq_len(nrow(removed_non_standard_df))){
        locus <- removed_non_standard_df$locus[i]
        column <- if (locus %in% one_loci) locus else paste0(locus, "(1)")
        x_int_ns[removed_non_standard_df$row[i], column] <- removed_non_standard_df$index[i]
      }
      if (verbose>=2){
        verbose_print("Relabeled data has", sum(x_int_ns<0, na.rm = TRUE), "negative indices")
      }
      ns_locus_haplotype <- paste0(non_standard_df$locus,"_", non_standard_df$haplotype)
      dimnames_q <- list(ns_locus_haplotype, cluster_labels)
      
      # keep trying to estimate pi,q and p until convergence
      pi_opt <- estimate_pr_ns(x = x_int_ns, v = v_matrix, number_of_1_loci = number_of_1_loci, number_of_2_loci = number_of_2_loci, locus_names = loci)
      q_opt <- estimate_q(x = x_int_ns, v = v_matrix, non_standard_haplotypes = non_standard_df, number_of_1_loci = number_of_1_loci, number_of_2_loci = number_of_2_loci)
      dimnames(q_opt) <- dimnames_q
      
      conv <- FALSE
      while(!conv){
        
        # optimise the variances again including the non-standard haplotypes
        f_ns <- function(theta) {
          negll <- disclapmix2::neg_loglik_theta_ns(theta, x_int_ns, y, pi_opt, q_opt, number_of_1_loci, number_of_2_loci)
          
          # if (is.nan(negll) || is.infinite(negll)){
          #   theta_fail <<- theta
          #   browser()
          # }
          negll
        }
        
        opt_ns <- stats::optim(par = theta, fn = f_ns, method = "BFGS", control = list(reltol = 1e-16, trace = optim_trace))
        if (opt$convergence != 0) stop("BFGS failed to converge")
        theta <- opt_ns$par
        
        # make v-matrix again
        tau_opt <- get_tau(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
        p_opt <- get_P(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
        rownames(p_opt) <- cluster_labels
        colnames(p_opt) <- loci
        
        # recreate the v matrix
        x_profile_pr_by_cluster <- compute_profile_prs_ns(p_by_cluster_and_locus = p_opt, db = x_int, y = y, pi_opt, q_opt, number_of_1_loci, number_of_2_loci)
        v_matrix <- compute_posterior_cluster_prs(profile_pr = x_profile_pr_by_cluster, tau = tau_opt)
        
        # estimate pi
        pi_opt_next <- estimate_pr_ns(x = x_int_ns, v = v_matrix, number_of_1_loci = number_of_1_loci, number_of_2_loci = number_of_2_loci, locus_names = loci)
        
        # estimate q
        q_opt_next <- estimate_q(x = x_int_ns, v = v_matrix, non_standard_haplotypes = non_standard_df, number_of_1_loci = number_of_1_loci, number_of_2_loci = number_of_2_loci)
        dimnames(q_opt_next) <- dimnames_q
        
        abs_diff_pi <- sum(abs(pi_opt_next - pi_opt))
        abs_diff_q <- sum(abs(q_opt_next - q_opt))
        
        if (verbose >=2){
          verbose_print("Abs difference in pi:", abs_diff_pi)
          verbose_print("Abs difference in q:", abs_diff_q)
        }
        
        conv <- abs_diff_pi<1e-12 & abs_diff_q <1e-12
        pi_opt <- pi_opt_next
        q_opt <- q_opt_next
      }
      
      # check if the centers are optimal
      y_new <- move_centers(x_int, y, v_matrix) # note that we need NAs for non-standard haplotypes here
      
      dist_new_y <- sum(abs(y - y_new))
      
      if (dist_new_y > 0){
        verbose_print("Centers need to be moved after second stage")
        verbose_print("Number of stepwise mutations between center configurations = ", dist_new_y)
        
        y <- y_new
      } 
      else{
        verbose_print("Second stage finished: centers need not be moved")
        break
      }
    }
  }
  
  number_of_iterations <- length(y_iterations)
  
  if (verbose >= 1) {
    iterations_label <- function(i) if (i==1) "iteration" else "iterations"
    verbose_print("Finished after", number_of_iterations, iterations_label(number_of_iterations))
  }
  
  x_profile_pr = rowSums(rep(tau_opt, each= nrow(x_int)) * x_profile_pr_by_cluster)
  
  log_lik = sum(log(x_profile_pr))
  
  ret <-  list(
    y_iterations = y_iterations, theta_iterations = theta_iterations,
    x_int = x_int, profile_pr_by_cluster = x_profile_pr_by_cluster, 
    profile_pr = x_profile_pr, posterior_cluster_pr = v_matrix, 
    log_lik = log_lik,
    y = y, p = p_opt
  )
  
  if(exists("x_int_ns")) ret$x_int_ns <- x_int_ns
  if(exists("x_stripped_int")) ret$x_stripped_int <- x_stripped_int
  if(exists("pi_opt")) ret$pi <- pi_opt
  if(exists("q_opt")) ret$q <- q_opt
  
  if(!exists("q_opt")){
    ret$number_of_parameters <- length(theta) + length(y)
    ret$BIC <- -2 * log_lik + ret$number_of_parameters * log(length(x_int))
  }else{
    ret$number_of_parameters <- length(theta) + length(y) + length(q_opt)
    ret$BIC <- -2 * log_lik + ret$number_of_parameters * log(length(x_int_ns))
  }
  
  ret$tau <- tau_opt
  
  ret$one_loci <- one_loci
  ret$two_loci <- two_loci
  
  ret
}

check_input_db <- function(x){
  if (!is.data.frame(x)) stop("x should be a data frame")
  if (!all(sapply(x,is.character))) stop("columns of x should be character vectors")
  
  x_unlist <- unlist(x)
  illegal_haplotypes <- x_unlist[grep(" ", x_unlist)]
  if (length(illegal_haplotypes) > 1){
    stop("Found haplotypes containing space: ", illegal_haplotypes)
  }
}

to_int_db <- function(x, one_loci, two_loci){
  
  # convert the stripped down db to an integer matrix
  x_int_one_loci <- matrix(as.integer(unlist(x[, one_loci])), 
                           nrow = nrow(x))
  colnames(x_int_one_loci) <- one_loci
  
  x_int_two_loci <- matrix(integer(), nrow=nrow(x), ncol=2 * length(two_loci))
  if(length(two_loci) > 0) colnames(x_int_two_loci) <- paste0(rep(two_loci, each=2),c("(1)","(2)"))
  
  for( locus in two_loci){
    x_locus <- x[,locus]
    locus_split <- ifelse(is.na(x_locus), yes = list(c(NA_character_,NA_character_)), no =strsplit(x_locus, split = ","))
    
    x_int_two_loci[,paste0(locus,c("(1)","(2)"))] <- matrix(as.integer(unlist(locus_split)), nrow = nrow(x), ncol=2,byrow = TRUE)
  }
  
  x_int <- cbind(x_int_one_loci, x_int_two_loci)
  
  x_int
}
