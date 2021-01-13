
move_centers <- function(x, y, v_matrix) {
  
  y_candidate_range_by_locus <- apply(x, 2L, function(x) range(x, na.rm = TRUE))
  
  new_y <- matrix(data=integer(), nrow = nrow(y), ncol=ncol(y), dimnames = dimnames(y))
  
  for(i_locus in seq_len(ncol(y))){
    for(i_cluster in seq_len(nrow(y))){
      
      y_candidates <- seq(from=y_candidate_range_by_locus[1,i_locus], 
                            to=y_candidate_range_by_locus[2,i_locus])
      
      loss <- sapply(y_candidates, function(y_candidate){
        sum(v_matrix[, i_cluster] * abs(x[, i_locus] - y_candidate), na.rm = TRUE) })
      
      new_y[i_cluster, i_locus] <- y_candidates[which.min(loss)]
    }
  }
  
  colnames(new_y) <- colnames(y)
  return(new_y)
}