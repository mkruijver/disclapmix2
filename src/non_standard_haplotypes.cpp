#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix estimate_pr_ns(IntegerMatrix x, NumericMatrix v, 
                             int number_of_1_loci, int number_of_2_loci,
                             CharacterVector locus_names) {
  // estimates the probability of a non-standard allele
  // in cluster (row) at locus (column)
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = v.ncol();
  
  int n = x.nrow();
  
  if (v.nrow() != n) Rcpp::stop("v needs as many rows as x");
  if (x.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("the number of columns of x needs to be number_of_1_loci + 2*number_of_2_loci");
  } 
  if (locus_names.length() != number_of_loci){
    Rcpp::stop("length of locus_names needs to be equal to number_of_1_loci + number_of_2_loci");
  }
  
  NumericMatrix pi(number_of_clusters, number_of_loci);
  
  for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
    for(int i_locus = 0; i_locus < number_of_loci; i_locus++){
      int i_column = i_locus < number_of_1_loci ? 
        i_locus : number_of_1_loci + (i_locus - number_of_1_loci) * 2;
      
      double mass_total = 0.0;
      double mass_ns = 0.0;
      
      for(int i = 0; i < n; i++){
        // ignore actual NAs (not negative integers)
        if (x(i, i_column) == NA_INTEGER) continue;
        
        mass_total += v(i, i_cluster);
        if (x(i,i_column) < 0){
          mass_ns += v(i, i_cluster);
        }
      }
      
      if (mass_total > 0){
        pi(i_cluster, i_locus) = mass_ns / mass_total;
      }
      else{
        pi(i_cluster, i_locus) = 1.0;
      }
    }
  }
  
  Function R_paste0("paste0");   
  pi.attr("dimnames") = List::create(R_paste0("cluster", Rcpp::seq_len(number_of_clusters)), 
          locus_names );
  
  return pi;
}

// [[Rcpp::export]]
NumericMatrix estimate_q(IntegerMatrix x, NumericMatrix v, 
                         DataFrame non_standard_haplotypes,
                         int number_of_1_loci, int number_of_2_loci) {
  // estimates the probability of each non-standard allele
  // in cluster (row) at locus (column)
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = v.ncol();
  
  int n = x.nrow();
  
  if (v.nrow() != n) Rcpp::stop("v needs as many rows as x");
  if (x.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("the number of columns of x needs to be number_of_1_loci + 2*number_of_2_loci");
  } 

  int n_ns = non_standard_haplotypes.nrow();
  IntegerVector locus_idx_by_h = non_standard_haplotypes["locus"];
  
  NumericMatrix mass_cluster_locus(number_of_clusters, number_of_loci);
  NumericMatrix q(n_ns, number_of_clusters);
  
  for(int i_locus = 0; i_locus < number_of_loci; i_locus++){
    int i_column = i_locus < number_of_1_loci ? 
    i_locus : number_of_1_loci + (i_locus - number_of_1_loci) * 2;
    
    for(int i = 0; i < n; i++){
      // ignore actual NAs (not negative integers)
      if (x(i, i_column) == NA_INTEGER) continue;
      
      for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
        mass_cluster_locus(i_cluster, i_locus) += v(i, i_cluster);
        
        if (x(i, i_column) < 0){
          // the first haplotype has index -1, then -2, etc.
          // so the corresponding 0-based row in the df is -1 - x
          int i_haplotype = -1 - x(i, i_column);
          q(i_haplotype, i_cluster) += v(i, i_cluster);
        }
      }
    }
  }
  
  // normalise to obtain probabilities
  for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
    for(int i_haplotype = 0; i_haplotype < n_ns; i_haplotype++){
      
      int i_locus = locus_idx_by_h[i_haplotype] - 1;
      
      if (mass_cluster_locus(i_cluster, i_locus) > 0){
        q(i_haplotype, i_cluster) /= mass_cluster_locus(i_cluster, i_locus);
      }
    }
  }
  
  return q;  
}
