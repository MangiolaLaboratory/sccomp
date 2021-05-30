library(dplyr)
library(tibble)
data_generation_intercept_only = function(N, Tm, M, p_n) {
  # Sample N proportions for categories
  mu = runif(N, min = 0, max = 10)
  mu = mu / sum(mu)
  
  # Generate proportions
  I = matrix(c(rep(1, N)), nrow = N)
  beta = (p_n * (I - mu))
  alpha = (p_n * mu)
  proportion_begin = 0
  for (i in 1:length(I)) {
    proportion_i = rbeta(M, alpha[i], beta[i])
    proportion_begin = cbind(proportion_begin, proportion_i)
  }
  proportion = as.matrix(proportion_begin)[,1:N + 1]

  # Generate beta binomial cell count
  Cmn_md = matrix(c(rep(rep(0, N), M)), nrow = M)
  for (i in 1:M) {
    for (j in 1:N) {
      Cmn_md[i,j] =rbinom(1, Tm, proportion[i,j])
    }
  }
  count = as.data.frame(Cmn_md)
  colname=character(N)
  for (j in 1:N){
    col_id=j
    colname[j]=paste("categpry",col_id)
  }
  colnames(count)=colname
  rowname=character(M)
  for (i in 1:M){
    row_id=i
    rowname[i]=paste("subject",row_id)
  }
  rownames(count)=rowname
  return(count)
}
#example
data_generation_intercept_only(N = 10,
               Tm = 300,
               M = 5,
               p_n = exp(4))

