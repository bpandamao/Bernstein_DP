# A matrix of Beta parameters and use Stick-breaking process to come up with the weights
# L: Truncation parameter
# M: Shape parameter for breaks
# K: the degree of Bernstein polynomial
# G0: Base distribution (Default: Uniform) 
# ...: Additional arguments in G0. 
# e.g., min and max for runif, mean and sd for rnorm, etc

Bernstein_DP=function(L,M=1,K,G0=runif,...){
  # the Stick-breaking process to have Zl
  # Sample V
  V = rbeta(L, 1, M)
  
  # Construct stick-breaking probabilities (p from V)
  p = rep(NA, L)
  p[1] = V[1]
  for (i in 2:L) {
    p[i] = prod(1 - V[1:(i - 1)]) * V[i]
  }
  p = c(1 - sum(p), p)  # Attach p0 as p[1]
  
  # Sample U from base measure G0
  U = G0(L + 1, ...)

  # matrix of beta density
  Beta_alpha=seq(1,K)
  Beta_beta=K-seq(1,K)+1

  #band partition
  Zbound = matrix(0,K,2)
  for (i in 1:K) {
    Zbound[i,] = c((i-1)/K,i/K)
  }

  # weights of each beta density
  weight = rep(0, K) 
  for (i in 1:nrow(Zbound)) {
    for (j in 1:L+1){
      if ((U[j]>Zbound[i,1]) & (U[j]<=Zbound[i,2])){
        weight[i]=weight[i]+p[j]
        }
    }
  }

  # Output
  return(list(Beta_alpha = Beta_alpha,
              Beta_beta = Beta_beta,
              weight = weight))
}

# example
L=10000
M=1
K=10
Bernstein_DP(L, M, K, runif, min = 0, max = 1)      # Base measure: U[0, 1]
Bernstein_DP(L, M, K, rnorm, mean = 0, sd = 10)     # Base measure: N(0, 100)
Bernstein_DP(L, M, K, rgamma, shape = 1, rate = 1)  # Base measure: Gamma(1, 1)
