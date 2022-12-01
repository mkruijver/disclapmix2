#include <Rcpp.h>
using namespace Rcpp;

const int number_of_precomputed_powers = 32;

// [[Rcpp::export]]
NumericMatrix get_P(NumericVector theta, int number_of_loci, int number_of_clusters) {
  // obtains the the discrete Laplace parameters (p) by cluster (row) and locus (column)
  // from the parametrisation
  
  int r = number_of_loci;
  int c = number_of_clusters;
  
  if (theta.size()!=(c-1 + r+c-1)) Rcpp::stop("theta needs length number_of_clusters-1 + number_of_loci+number_of_clusters-1");
  
  // theta is the full parameter vector:
  // (theta_tau, theta_pjk)
  
  // we use theta_pjk to parametrise the locus and cluster effects:
  // theta_pjk = (omega_1, ..., omega_c, lambda_1, ..., lambda_{r-1})
  
  // log p_{jk} = omega_j + lambda_k
  //              (age)     (mutation rate)
  auto p = NumericMatrix(c, r);
  
  if (number_of_clusters==1){
    for(int i_locus=0; i_locus < number_of_loci; i_locus++){
      p[i_locus] = std::exp(theta[i_locus]);
      
      if (p[i_locus] > 0.99) p[i_locus] = 0.99;
    }
    
    return p;
  }
  
  for(auto i_cluster = 0; i_cluster < c; i_cluster++){   // cluster
    for(auto i_locus = 0; i_locus < r; i_locus++){ // locus
      
      if (i_locus==0){
        p(i_cluster, i_locus) = std::exp(theta[c - 1 + i_cluster]);
      }else{
        p(i_cluster, i_locus) = std::exp(theta[c - 1 + i_cluster] + theta[c - 1 + c + i_locus - 1]);
      }
      
      if (p(i_cluster, i_locus)>0.99) p(i_cluster, i_locus) = 0.99;
    }
    
  }
  
  return p;
}

// [[Rcpp::export]]
NumericVector get_tau(NumericVector theta, int number_of_loci, int number_of_clusters){
  int r = number_of_loci;
  int c = number_of_clusters;
  
  if (theta.size()!=(c-1 + r+c-1)) Rcpp::stop("theta needs length number_of_clusters-1 + number_of_loci+number_of_clusters-1");

  NumericVector tau(number_of_clusters);
  
  // extract tau from parameter vector
  double tau_cum = 0.0;
  
  for(int i_cluster = 0; i_cluster < number_of_clusters - 1; i_cluster++){
    
    double tau_i = exp(theta[i_cluster]);
    
    tau[i_cluster] = tau_i;
    
    tau_cum += tau_i;
  }
  
  tau[number_of_clusters - 1] = std::max(0.0, 1.0 - tau_cum);
  
  return tau;
}

std::vector<NumericMatrix> precompute_dlm_powers(NumericMatrix p_by_cluster_and_locus){
  std::vector<NumericMatrix> prs_by_cluster;
  
  int number_of_clusters = p_by_cluster_and_locus.nrow();
  int number_of_loci = p_by_cluster_and_locus.ncol();
  
  for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
    NumericMatrix prs(number_of_precomputed_powers, number_of_loci);
    
    for(int i_locus = 0; i_locus < number_of_loci; i_locus++){
      
      double p = p_by_cluster_and_locus(i_cluster, i_locus);
      
      prs(0, i_locus) = (1.0-p)/(1.0+p);
      
      for(int i = 1; i < number_of_precomputed_powers; i++){
        prs(i, i_locus) = (1.0-p)/(1.0+p) * std::pow(p,i);//prs(i-1, i_locus) * p;
        // prs(i, i_locus) = prs(i-1, i_locus) * p;
      }
    }
    
    prs_by_cluster.push_back(prs);
  }
  
  return prs_by_cluster;
}

void range_error(int x, int y, int i_profile){
      std::string error_message = "range outside of pre-computations: x=" + 
        std::to_string(x) + " y=" + std::to_string(y) + 
		" row = " + std::to_string(i_profile + 1);
      
      Rcpp::stop(error_message);
}

void range_error(int xa, int xb, int ya, int yb, int i_profile){
      std::string error_message = "range outside of pre-computations: x=" + 
        std::to_string(xa) + "," +  std::to_string(xb) +
		" y=" + std::to_string(ya) + "," +  std::to_string(yb)+ 
		" row = " + std::to_string(i_profile + 1);

      Rcpp::stop(error_message);
}

double compute_profile_pr_locus(int i_profile, int i_cluster, int i_locus,
                        std::vector<NumericMatrix> &prs_by_cluster, 
                        IntegerMatrix &db, IntegerMatrix &y, 
                        int number_of_1_loci,
                        int number_of_2_loci){
  
  if (i_locus < number_of_1_loci){
    // pr. for the 1 loci
    
    int x = db(i_profile, i_locus);
    if (x == NA_INTEGER) return(1.0);
    
    // Rcpp::Rcout << "i: " << i << " i_locus: " << i_locus << "\n";
    int delta = std::abs(x - y(i_cluster, i_locus));
    if (delta >= number_of_precomputed_powers){
      range_error(x, y(i_cluster, i_locus), i_profile);
    } 
    
    return(prs_by_cluster[i_cluster](delta, i_locus));
  }
  else{
    
    // pr. for the 2-loci
    int i_2_locus = i_locus - number_of_1_loci;
    
    int col_a = number_of_1_loci + 2 * i_2_locus;
    int col_b = number_of_1_loci + 2 * i_2_locus + 1;
    
    int x_a = db(i_profile, col_a);
    int x_b = db(i_profile, col_b);
    
    if (x_a == NA_INTEGER || x_b == NA_INTEGER) {
      return(1.0);
    };
    
    int y_a = y(i_cluster, col_a);
    int y_b = y(i_cluster, col_b);
    
    int delta_a_a = std::abs(x_a - y_a);
    int delta_b_b = std::abs(x_b - y_b);
    
    int delta_b_a = std::abs(x_b - y_a);
    int delta_a_b = std::abs(x_a - y_b);
    
    if (delta_a_a >= number_of_precomputed_powers || delta_b_b >= number_of_precomputed_powers 
          || delta_b_a >= number_of_precomputed_powers || delta_a_b >= number_of_precomputed_powers){
      range_error(x_a, x_b, y_a, y_b, i_profile);
    } 
    
    // Rcpp::Rcout << "i: " << i << " i_locus: " << i_locus << "\n";
    double pr = 0.5 * 
      (prs_by_cluster[i_cluster](delta_a_a, i_locus) * prs_by_cluster[i_cluster](delta_b_b, i_locus) +
      prs_by_cluster[i_cluster](delta_b_a, i_locus) * prs_by_cluster[i_cluster](delta_a_b, i_locus));
    
    return pr;
  }
}

double compute_profile_pr(int i, int i_cluster, std::vector<NumericMatrix> &prs_by_cluster, 
                          IntegerMatrix &db, IntegerMatrix &y, int number_of_1_loci,
                          int number_of_2_loci){
  double pr_cluster = 1.0;
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  
  for (int i_locus = 0; i_locus < number_of_loci; i_locus++){
    double pr_cluster_locus = compute_profile_pr_locus(i, i_cluster, 
                                i_locus, prs_by_cluster, db, y, 
                                number_of_1_loci, number_of_2_loci);
    
    pr_cluster *= pr_cluster_locus;
  }
  
  return pr_cluster;
} 
  
double compute_profile_pr_ns(int i, int i_cluster, std::vector<NumericMatrix> &prs_by_cluster, 
                            IntegerMatrix &db, IntegerMatrix &y, 
                            NumericMatrix &pi, NumericMatrix &q,
                            int number_of_1_loci,
                            int number_of_2_loci){
    double pr_cluster = 1.0;
    
    // pr. for the 1-loci
    for(int i_locus = 0; i_locus < number_of_1_loci; i_locus++){
      
      int x = db(i,i_locus);
      
      if (x == NA_INTEGER) continue;
      
      if (x > 0){
        // standard haplotype
        int delta = std::abs(x - y(i_cluster, i_locus));
        if (delta >= number_of_precomputed_powers) range_error(x, y(i_cluster, i_locus), i);
        
        pr_cluster *= (1.0 - pi(i_cluster, i_locus)) *
          prs_by_cluster[i_cluster](delta, i_locus);
      }
      else{
        // non-standard haplotype
        int i_haplotype = -1 - x;
        pr_cluster *= q(i_haplotype, i_cluster);
        // Rcpp::Rcout << "Row " << i+1 << " haplotype " << x <<
        //             " (idx is " << i_haplotype << ")" <<
        //           " pr in cluster "  << i_cluster +1 << " is " <<
        //             q(i_haplotype, i_cluster) << "\n";
        
      }
    }
    
    // pr. for the 2-loci
    for(int i_2_locus = 0; i_2_locus < number_of_2_loci; i_2_locus++){
      int i_locus = number_of_1_loci + i_2_locus;
      
      int col_a = number_of_1_loci + 2 * i_2_locus;
      int col_b = number_of_1_loci + 2 * i_2_locus + 1;
      
      int x_a = db(i, col_a);
      int x_b = db(i, col_b);
      
      if (x_a==NA_INTEGER || x_b==NA_INTEGER) continue;
      
      if (x_a > 0){
        // standard haplotype
        
        int y_a = y(i_cluster, col_a);
        int y_b = y(i_cluster, col_b);
        
        int delta_a_a = std::abs(x_a - y_a);
        int delta_b_b = std::abs(x_b - y_b);
        
        int delta_b_a = std::abs(x_b - y_a);
        int delta_a_b = std::abs(x_a - y_b);
        
        if (delta_a_a >= number_of_precomputed_powers || delta_b_b >= number_of_precomputed_powers 
              || delta_b_a >= number_of_precomputed_powers || delta_a_b >= number_of_precomputed_powers){	        
			range_error(x_a, x_b, y_a, y_b, i);
        } 
        
        // Rcpp::Rcout << "i: " << i << " i_locus: " << i_locus << "\n";
        pr_cluster *= (1.0 - pi(i_cluster, i_locus)) * 0.5 * 
          (prs_by_cluster[i_cluster](delta_a_a, i_locus) * prs_by_cluster[i_cluster](delta_b_b, i_locus) +
          prs_by_cluster[i_cluster](delta_b_a, i_locus) * prs_by_cluster[i_cluster](delta_a_b, i_locus));
        
      }
      else{
        // non-standard haplotype
        int i_haplotype = -1 - x_a;
        pr_cluster *= q(i_haplotype, i_cluster);
      }
      
    }
    
    return pr_cluster;
  } 
  
// [[Rcpp::export]]
double loglik_tau_p(NumericVector tau, NumericMatrix p_by_cluster_and_locus, IntegerMatrix db, IntegerMatrix y,
                int number_of_1_loci, int number_of_2_loci) {
  
  int n = db.nrow();

  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = tau.length();
  
  if (p_by_cluster_and_locus.nrow() != number_of_clusters){
    Rcpp::stop("p should have as many rows as length of tau");
  }
  if (p_by_cluster_and_locus.ncol() != number_of_loci){
    Rcpp::stop("p should have as many columns as number of loci");
  }
  if (db.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("db should have as many columns as number_of_1_loci + 2*number_of_2_loci");
  }
  if (y.nrow() != number_of_clusters){
    Rcpp::stop("y should have as many rows as length of tau");
  }
  if (y.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("y should have as many columns as number_of_1_loci + 2 * number_of_2_loci");
  }
  
  // make sure tau is valid
  double tau_sum = 0.0;
  double tau_penalty = 0.0;
  for(int i = 0; i < tau.size(); i++){
    tau_sum += tau[i];
    if (tau[i] < 0) return R_NegInf;
  }
  if (tau_sum > 1){
    tau_penalty -= (tau_sum-1) * 1e7;  
  }
  
  // pre-compute part of the discrete Laplace pmf for each cluster and locus
  std::vector<NumericMatrix> prs_by_cluster = precompute_dlm_powers(p_by_cluster_and_locus);

  double loglik = 0.0;
  
  // for each profile
  for(int i = 0; i < n; i++){
    double pr = 0.0;
    
    // compute the pr. in each cluster
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      
      double pr_cluster = compute_profile_pr(i, i_cluster, prs_by_cluster, 
                                             db, y, number_of_1_loci, number_of_2_loci);
      
      pr += tau[i_cluster] * pr_cluster;
    }
    
    loglik += std::log(pr);
  }
  
  return loglik + tau_penalty;
}

// [[Rcpp::export]]
double loglik_tau_p_ns(NumericVector tau, NumericMatrix p_by_cluster_and_locus, 
                       IntegerMatrix db, IntegerMatrix y,
                       NumericMatrix pi, NumericMatrix q,
                       int number_of_1_loci, int number_of_2_loci) {
  
  int n = db.nrow();
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = tau.length();
  
  if (p_by_cluster_and_locus.nrow() != number_of_clusters){
    Rcpp::stop("p should have as many rows as length of tau");
  }
  if (p_by_cluster_and_locus.ncol() != number_of_loci){
    Rcpp::stop("p should have as many columns as number of loci");
  }
  if (db.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("db should have as many columns as number_of_1_loci + 2*number_of_2_loci");
  }
  if (y.nrow() != number_of_clusters){
    Rcpp::stop("y should have as many rows as length of tau");
  }
  if (y.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("y should have as many columns as number_of_1_loci + 2 * number_of_2_loci");
  }
  if (pi.nrow() != number_of_clusters){
    Rcpp::stop("pi should have as many rows as length of tau");
  }
  if (pi.ncol() != number_of_loci){
    Rcpp::stop("pi should have as many columns as number of loci");
  }
  
  // make sure tau is valid
  for(int i = 0; i < tau.size(); i++){
    if (tau[i] < 0 || tau[i] >1) return R_NegInf;
  }
  
  // pre-compute part of the discrete Laplace pmf for each cluster and locus
  std::vector<NumericMatrix> prs_by_cluster = precompute_dlm_powers(p_by_cluster_and_locus);
  
  double loglik = 0.0;
  
  // for each profile
  for(int i = 0; i < n; i++){
    double pr = 0.0;
    
    // compute the pr. in each cluster
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      
      double pr_cluster = compute_profile_pr_ns(i, i_cluster, prs_by_cluster, 
                                             db, y, pi, q,
                                             number_of_1_loci, number_of_2_loci);
      
      pr += tau[i_cluster] * pr_cluster;
    }
    
    loglik += std::log(pr);
  }
  
  return loglik;
}

// [[Rcpp::export]]
double neg_loglik_theta(NumericVector theta, IntegerMatrix db, IntegerMatrix y,
                    int number_of_1_loci, int number_of_2_loci) {
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = y.nrow();
  
  // obtain model parameters from theta
  NumericVector tau = get_tau(theta, number_of_loci, number_of_clusters);
  NumericMatrix p = get_P(theta, number_of_loci, number_of_clusters);
  
  double ll = loglik_tau_p(tau, p, db, y, number_of_1_loci, number_of_2_loci);
  
  return -ll;
}


// [[Rcpp::export]]
NumericMatrix compute_profile_prs(NumericMatrix p_by_cluster_and_locus, IntegerMatrix db, IntegerMatrix y,
                    int number_of_1_loci, int number_of_2_loci) {
  
  int n = db.nrow();
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = p_by_cluster_and_locus.nrow();
  
  NumericMatrix pr(n, number_of_clusters);
    
  if (p_by_cluster_and_locus.ncol() != number_of_loci){
    Rcpp::stop("p should have as many columns as number of loci");
  }
  if (db.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("db should have as many columns as number_of_1_loci + 2*number_of_2_loci");
  }
  if (y.nrow() != number_of_clusters){
    Rcpp::stop("y should have as many rows as length of tau");
  }
  if (y.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("y should have as many columns as number_of_1_loci + 2 * number_of_2_loci");
  }
  
  // pre-compute part of the discrete Laplace pmf for each cluster and locus
  std::vector<NumericMatrix> prs_by_cluster = precompute_dlm_powers(p_by_cluster_and_locus);
  
  // for each profile
  for(int i = 0; i < n; i++){
    
    // compute the pr. in each cluster
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      
      double pr_cluster = compute_profile_pr(i, i_cluster, prs_by_cluster, 
                                             db, y, number_of_1_loci, number_of_2_loci);
      
      
      pr(i, i_cluster) = pr_cluster;
    }
  }
  
  return pr;
}

// [[Rcpp::export]]
List compute_profiles_pr_by_cluster_and_locus(
    NumericMatrix p_by_cluster_and_locus, IntegerMatrix x, IntegerMatrix y,
    int number_of_1_loci, int number_of_2_loci) {
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = p_by_cluster_and_locus.nrow();
  
  int number_of_profiles = x.nrow();
  
  List prs = List(number_of_profiles);
  
  if (p_by_cluster_and_locus.ncol() != number_of_loci){
    Rcpp::stop("p should have as many columns as number of loci");
  }
  if (x.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("db should have as many columns as number_of_1_loci + 2*number_of_2_loci");
  }
  if (y.nrow() != number_of_clusters){
    Rcpp::stop("y should have as many rows as length of tau");
  }
  if (y.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("y should have as many columns as number_of_1_loci + 2 * number_of_2_loci");
  }
  
  // pre-compute part of the discrete Laplace pmf for each cluster and locus
  std::vector<NumericMatrix> prs_by_cluster = precompute_dlm_powers(p_by_cluster_and_locus);
  
  for (int i_profile = 0; i_profile < number_of_profiles; i_profile++){
    NumericMatrix pr(number_of_clusters, number_of_loci);
    
    // compute the pr. in each cluster
    for (int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      for (int i_locus = 0; i_locus < number_of_loci; i_locus++){
        
        double pr_cluster_locus = compute_profile_pr_locus(i_profile, i_cluster, 
                                                           i_locus, prs_by_cluster, x, y, 
                                                           number_of_1_loci, number_of_2_loci);
        
        pr(i_cluster, i_locus) = pr_cluster_locus;
      }
      
    }
    
    pr.attr("dimnames") = p_by_cluster_and_locus.attr("dimnames");
   
    prs[i_profile] = pr; 
  }
  
  return prs;
}

// [[Rcpp::export]]
NumericMatrix compute_profile_prs_ns(NumericMatrix p_by_cluster_and_locus, IntegerMatrix db, IntegerMatrix y,
                                     NumericMatrix pi, NumericMatrix q,
                                  int number_of_1_loci, int number_of_2_loci) {
  
  int n = db.nrow();
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = p_by_cluster_and_locus.nrow();
  
  NumericMatrix pr(n, number_of_clusters);
  
  if (p_by_cluster_and_locus.ncol() != number_of_loci){
    Rcpp::stop("p should have as many columns as number of loci");
  }
  if (db.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("db should have as many columns as number_of_1_loci + 2*number_of_2_loci");
  }
  if (y.nrow() != number_of_clusters){
    Rcpp::stop("y should have as many rows as number of clusters");
  }
  if (y.ncol() != number_of_1_loci + 2*number_of_2_loci){
    Rcpp::stop("y should have as many columns as number_of_1_loci + 2 * number_of_2_loci");
  }
  if (pi.nrow() != number_of_clusters){
    Rcpp::stop("pi should have as many rows as number of clusters");
  }
  if (pi.ncol() != number_of_loci){
    Rcpp::stop("pi should have as many columns as number of loci");
  }  
  
  // pre-compute part of the discrete Laplace pmf for each cluster and locus
  std::vector<NumericMatrix> prs_by_cluster = precompute_dlm_powers(p_by_cluster_and_locus);
  
  // for each profile
  for(int i = 0; i < n; i++){
    
    // compute the pr. in each cluster
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      
      double pr_cluster = compute_profile_pr_ns(i, i_cluster, prs_by_cluster, 
                                             db, y, pi, q,
                                             number_of_1_loci, number_of_2_loci);
      
      pr(i, i_cluster) = pr_cluster;
    }
  }
  
  return pr;
}

// [[Rcpp::export]]
NumericMatrix compute_posterior_cluster_prs(NumericMatrix profile_pr, NumericVector tau){
  int n = profile_pr.nrow();
  int number_of_clusters = profile_pr.ncol();
  
  if (tau.length() != number_of_clusters) Rcpp::stop("tau should have length equal to number of columns in profile_pr");
  
  NumericMatrix posterior(n, number_of_clusters);
  
  for(int i = 0; i < n; i++){
    
    double total_pr = 0.0;
    
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      total_pr += tau[i_cluster] * profile_pr(i, i_cluster);
    }
    
    for(int i_cluster = 0; i_cluster < number_of_clusters; i_cluster++){
      posterior(i, i_cluster) = (1.0/total_pr) * tau[i_cluster] * profile_pr(i, i_cluster);
    }
    
  }
  
  return posterior;
}

// [[Rcpp::export]]
double neg_loglik_theta_ns(NumericVector theta, IntegerMatrix db, IntegerMatrix y,
                           NumericMatrix pi, NumericMatrix q,
                        int number_of_1_loci, int number_of_2_loci) {
  // theta is the parameter vector (parametrises the variances)
  // pi contains the pr. of a non-standard haplotype by cluster (row) and locus (column)
  // q  contains the pr's of each non-standard haplotype (row) by cluster (column)
  
  int number_of_loci = number_of_1_loci + number_of_2_loci;
  int number_of_clusters = y.nrow();
  
  // obtain model parameters from theta
  NumericVector tau = get_tau(theta, number_of_loci, number_of_clusters);
  NumericMatrix p = get_P(theta, number_of_loci, number_of_clusters);
  
  double ll = loglik_tau_p_ns(tau, p, db, y, pi, q, number_of_1_loci, number_of_2_loci);
  
  return -ll;
}
