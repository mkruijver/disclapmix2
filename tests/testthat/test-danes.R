test_that("danes", {

  require(disclapmix)
  
  data(danes) 
  x <- as.matrix(danes)
  x2 <- as.data.frame(sapply(danes, as.character))

  for(K in 1:5){
    dlm_fit <- disclapmix(x, K)
    
    dlm2_fit <- disclapmix2(x2, number_of_clusters = K)
    
    testthat::expect_equal(dlm_fit$logL_marginal, dlm2_fit$log_lik)
  }
  
})