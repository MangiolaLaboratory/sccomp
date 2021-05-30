#generate sum-to-one random numbers
#input N (# of groups)
N = 10
#input of Tm
Tm = 200
#input of M (nmubers of subjects in each group)
M = 5
# input precision coefficients
p_n = exp(3)
#Tm input total counts
#function
data_generation = function(X, coeff, Tm,  p_n){
  M = nrow(X)
  N = ncol(coeff)
  
  # X * coeff (-Inf, Inf)
  # mu (0,q)
  mu[m]  = softmax( X[m] * coeff )
}


library(dplyr)
library(tibble)
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}
#input coefficients must be column vector
#input design matrix must be matrix that can left-multiple coefficient
data_generation = function(X, Tm, coeff, p_n) {
  # Sample N proportions for categories
  M=nrow(coeff)
  N=nrow(X)
  mu=softmax(X%*%coeff)
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
trial=matrix(runif(50,min = 0,max=10),nrow=5,ncol=10)
data_generation(X=t(trial),
               Tm = 300,
               coeff=matrix(rep(0.3,5),nrow=5),
               p_n = exp(4))



