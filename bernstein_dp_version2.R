# version2
# get the weights
# add the weights to the prior curve
# plot the prior curve
set.seed(42)
x = seq(0, 1, length=10000)
K=12 # 12 beta mixture
L=10000
M=1
G0=runif

# the Stick-breaking process to have Zl
# Sample V
V1 = rbeta(L, 1, M)

# Construct stick-breaking probabilities (p from V)
p = rep(NA, L)
p[1] = V1[1]
for (i in 2:L) {
  p[i] = prod(1 - V1[1:(i - 1)]) * V1[i]
}
p = c(1 - sum(p), p)  # Attach p0 as p[1]

# Sample U from base measure G0
Z = G0(L + 1)

# get the beta density function
density_beta=function(x,j,K){
  alpha=j 
  beta=K-j+1
  out=dbeta(x,alpha,beta)
}

# come out the weights
weights=function(j,K,p,Z){
  # get Xbound
  boundd=(j-1)/K
  boundu=j/K
  # get according q_sum
  out=0
  for (i in 1:length(p)){
    if (Z[i]<=boundu&Z[i]>boundd){out=out+p[i]}
  }
  return(out)
}

# get the mixture beta density
betam=function(x,K){
  out=matrix(0,length(x),K)
  for(i in 1:K){
    out[,i]=density_beta(x,i,K)*weights(i,K,p,Z)
  }
  return(out)
}

# Plot functions
betamatrix=betam(x,K)
betamatrixcs=apply(betamatrix,1,sum)
plot(x, betamatrix[,1], type = "l", ylim = c(min(betamatrix), max(betamatrix)), main = "beta mixture prior")
for (i in 2:ncol(betamatrix)){lines(x, betamatrix[,i], col = i)}
lines(x,betamatrixcs,col='black')
