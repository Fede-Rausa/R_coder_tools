
contingency_test = function(A, B){
  #A and B are 2 vectors of the same length
  #this function return the statistics of the contigency table
  tab_freq_oss = table(A,B)
  N = length(A)
  r = length(unique(A))
  c = length(unique(B))
  freq_mar_r = rowSums(tab_freq_oss)
  freq_mar_c = colSums(tab_freq_oss)
  tab_freq_ind = (freq_mar_r %*% t(freq_mar_c))/N
  rownames(tab_freq_ind) = names(freq_mar_r)
  cont_table = tab_freq_oss - tab_freq_ind
  test_chi = sum(c(cont_table)^2/c(tab_freq_ind))
  df = (r-1)*(c-1)
  cont_coeff = sqrt(test_chi/(test_chi+N))
  pvalue = 1 - pchisq(test_chi, df=df)
  return(list(contingency_coefficent = cont_coeff,
              test_chi2 = test_chi,
              pvalue = pvalue,
              contingency_table = cont_table,
              observed_joint_frequencies = tab_freq_oss,
              joint_frequencies_under_independence = tab_freq_ind))
}

