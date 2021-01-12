# rm(list=ls())
check_input_db <- function(x){
  if (!is.data.frame(x)) stop("x should be a data frame")
  if (!all(sapply(x,is.character))) stop("columns of x should be character vectors")
}

# x_raw <- readxl::read_excel(system.file("extdata", "South Australia.xlsx", package = "disclapmix2"), col_types = "text", sheet = 2)
# x <- x_raw[-c(1,2,21,26)]
# x_raw <- readxl::read_excel(system.file("extdata", "Chinese Han.xlsx", package = "disclapmix2"), col_types = "text")
# # x_raw
# x <- x_raw[rep(seq_len(nrow(x_raw)), times=as.integer(x_raw$Count)),-c(1:3)]
# include_2_loci=FALSE
# verbose=1
# number_of_clusters = 2L

#' @export
disclapmix2 <- function(x, number_of_clusters, include_2_loci = FALSE, verbose=0L){
  
  check_input_db(x)
  
  if (verbose>=1) disclapmix2:::verbose_print("Determining loci with suitable data")
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
  if (verbose>=1) disclapmix2:::verbose_print("Found", number_of_1_loci, paste0("suitable 1-", locuslabel(number_of_1_loci)),
                                              "and", number_of_2_loci,paste0("suitable 2-", locuslabel(number_of_2_loci)),
                                              "among", ncol(x), "columns")
  
  x_cleaned <- matrix(data = character(), nrow = nrow(x), ncol = number_of_loci)
  colnames(x_cleaned) <- c(one_loci, if(number_of_2_loci>0) two_loci else NULL)
  
  for(locus in one_loci){
    ind_1 <- x_summarised$x_ind_12_other[,locus]=="1"
    x_cleaned[ind_1,locus] <- x[[locus]][ind_1]
  }
  
  if(number_of_2_loci > 0){
    for(locus in two_loci){
      ind_2 <- x_summarised$x_ind_12_other[,locus]=="2"
      x_cleaned[ind_2,locus] <- x[[locus]][ind_2]
    }
  }
  
  if (verbose>=1) disclapmix2:::verbose_print("Removing profiles with any haplotype that is not 1 or 2 integers at the approriate locus")
  # for the first stage estimation, we remove any profiles with NAs
  idx_included <- rowSums(is.na(x_cleaned))==0
  
  # create a stripped down version of the db with only the easy profiles
  x_stripped <- x_cleaned[idx_included,]
  
  # convert the stripped down db to an integer matrix
  x_int_one_loci <- matrix(as.integer(x_stripped[, one_loci]), nrow=nrow(x_stripped))
  colnames(x_int_one_loci) <- one_loci
  
  x_int_two_loci <- matrix(integer(), nrow=nrow(x_stripped), ncol=2 * number_of_2_loci)
  if(number_of_2_loci > 0) colnames(x_int_two_loci) <- paste0(rep(two_loci, each=2),c("(1)","(2)"))
  
  for( locus in two_loci){
    x_int_two_loci[,paste0(locus,c("(1)","(2)"))] <- matrix(as.integer(unlist(strsplit(x_stripped[,locus], split = ","))), nrow = nrow(x_stripped),ncol=2,byrow = TRUE)
  }
  
  x_int <- cbind(x_int_one_loci, x_int_two_loci)
  
  if (verbose>=1) disclapmix2:::verbose_print(nrow(x_int) , paste0("profiles left (started with ", nrow(x),")"))
  if (verbose>=1) disclapmix2:::verbose_print("Determining initial clustering using PAM")
  
  initial_pam <- cluster::pam(x_int, k = number_of_clusters, metric = "manhattan",diss = FALSE, keep.diss=FALSE, keep.data=TRUE)
  
  initial_pam_tau <- (tabulate(initial_pam$clustering, nbins = number_of_clusters)/nrow(x_int))
  initial_pam_theta_tau <- if(number_of_clusters==1L) numeric() else log( initial_pam_tau[1:(number_of_clusters-1)])
  
  # now find the variances
  y <- y0 <- initial_pam$medoids
  theta <- theta0 <- if(number_of_clusters>1) c(initial_pam_theta_tau, rep(-1,number_of_clusters), rep(0,number_of_loci-1)) else rep(-1,number_of_loci)
  
  if (verbose>=1) disclapmix2:::verbose_print("Starting optimisation")
  
  cluster_labels <- paste0("cluster", seq_len(number_of_clusters))
  
  y_iterations <- list()
  theta_iterations <- list()
  
  repeat{
    rownames(y) <- cluster_labels
    y_iterations[[1 + length(y_iterations)]] <- y
    
    # optim
    f <- function(theta) {
      negll <- disclapmix2:::neg_loglik_theta(theta, x_int, y, number_of_1_loci, number_of_2_loci)
      
      # if (is.nan(negll) || is.infinite(negll)){
      #   theta_fail <<- theta
      #   browser()
      # }
      negll
    }
    
    opt <- stats::optim(par = theta, fn = f, method = "BFGS", control = list(reltol = 1e-16))
    if (opt$convergence != 0) stop("BFGS failed to converge")
    
    # dput(opt$value)
    
    theta_opt <- opt$par
    
    theta_iterations[[1+length(theta_iterations)]] <- theta_opt
    
    tau_opt <- disclapmix2:::get_tau(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
    p_opt <- disclapmix2:::get_P(theta = theta_opt, number_of_loci = number_of_loci, number_of_clusters = number_of_clusters)
    
    rownames(p_opt) <- cluster_labels
    colnames(p_opt) <- loci
    
    x_int_profile_pr_by_cluster <- disclapmix2:::compute_profile_prs(p_by_cluster_and_locus = p_opt, db = x_int, y = y, number_of_1_loci, number_of_2_loci)
    v_matrix <- disclapmix2:::compute_posterior_cluster_prs(profile_pr = x_int_profile_pr_by_cluster, tau = tau_opt)
    
    # see if we need to move centers
    y_new <- disclapmix:::move_centers(x_int, y, v_matrix)
    
    theta <- theta_opt
    
    dist_new_y <- sum(abs(y - y_new))
    if (dist_new_y == 0) {
      if (verbose >= 1L) {
        disclapmix2:::verbose_print("Current central haplotypes are optimal")
      }
      break
    } else if (any(duplicated(y_new))) { # new case introduced in disclapmix version 1.6.3
      if (verbose >= 1L) {
        disclapmix2:::verbose_print("New central haplotypes had at least one duplicated haplotype, change rejected")
      }
      break
    } else {
      if (verbose >= 2L) {      
        disclapmix2:::verbose_print("Current central haplotypes are not otimal")
        print(y)
        disclapmix2:::verbose_print("New central haplotypes:")
        print(y_new)
        disclapmix2:::verbose_print("Differences:")
        print(y_new - y)
        disclapmix2:::verbose_print("Number of stepwise mutations between center configurations = ", dist_new_y)
        
      } else if (verbose >= 1) {
        disclapmix2:::verbose_print("Current central haplotypes are not optimal, moving and restarting optimisation")
      }
      
      y <- y_new
    }
  }
  
  number_of_iterations <- length(y_iterations)
  
  if (verbose >= 1) {
    iterations_label <- function(i) if (i==1) "iteration" else "iterations"
    disclapmix2:::verbose_print("Finished after", number_of_iterations, iterations_label(number_of_iterations))
  }
  
  x_int_profile_pr = rowSums(rep(tau_opt, each= nrow(x_int)) * x_int_profile_pr_by_cluster)
  
  log_lik = sum(log(x_int_profile_pr))
  
  list(
    y_iterations = y_iterations, theta_iterations = theta_iterations,
    x_int = x_int, profile_pr_by_cluster = x_int_profile_pr_by_cluster, 
    profile_pr = x_int_profile_pr, posterior_cluster_pr = v_matrix, 
    log_lik = log_lik,
    y = y, p = p_opt, tau = tau_opt
  )
}
