test_that("danes", {

  require(disclapmix)
  
  data(danes) 
  
  x <- as.matrix(danes[rep(seq_len(nrow(danes)), danes$n), -ncol(danes)])
  x2 <- as.data.frame(sapply(danes[rep(seq_len(nrow(danes)), danes$n), -ncol(danes)], as.character))
  
  for(K in 1:5){
    dlm_fit <- disclapmix(x, K)
    
    dlm2_fit <- disclapmix2(x2, number_of_clusters = K)
    
    testthat::expect_equal(dlm_fit$logL_marginal, dlm2_fit$log_lik)
  }
  
})