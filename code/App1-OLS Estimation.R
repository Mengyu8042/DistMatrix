library(pracma)

N = 10000
p_arr = c(20, 200, 500, 1000)
L = length(p_arr)

for (l in 1:L) {
  p = p_arr[l]
  gamma = p/N
  
  #create a test covariance matrix
  set.seed(12345)
  v = runif(p) + 1
  Sigma = diag(v)
  iSigma = solve(Sigma)
  
  #generate a random matrix with prescribed covariance matrix
  data = rnorm(N*p, 0, 1)
  X = matrix(data, N, p) %*% pracma::sqrtm(Sigma)$B
  
  #find the full-sample covariance matrix of beta_hat
  Sigma_h = t(X) %*% X / N
  iSigma_h = solve(Sigma_h)
  Gamma = iSigma_h / N
  
  #find the distributed covariance matrix of beta_hat
  dist_iSigma_h = list()
  dist_Gamma = list()
  K = floor(0.75*N/(p+1))
  for (k in 1:K) {
    n = floor(N/k)
    Sigma_h_k = list()
    iSigma_h_k = list()
    dist_iSigma_h[[k]] = matrix(0, p, p)
    
    for(j in 1:k) {
      Sigma_h_k[[j]] = t(X[((j-1)*n+1):(j*n), ]) %*% X[((j-1)*n+1):(j*n), ] / n
      iSigma_h_k[[j]] = solve(Sigma_h_k[[j]])
      dist_iSigma_h[[k]] = dist_iSigma_h[[k]] + iSigma_h_k[[j]] / k
      dist_Gamma[[k]] = dist_iSigma_h[[k]] / N
    }
  }
  
  #calculate and print relative efficiency
  re = rep(0, length = K)
  re_th = c()
  for (k in 1:K) {
    for (i in 1:p) {
      re[k] = re[k] + (Gamma[i,i] / dist_Gamma[[k]][i,i]) / p
    }
    re_th[k] = (N - k * p) / (N - p)
    print(c(p, k, re[k], re_th[k]))
  }
  
  #plot: relative efficiency of OLS estimation ~ number of machines
  plot(1:K, re, type = "l", lwd = 2, col = "blue", xlab = "Number of Machines", 
       ylab = "Relative Efficiency of OLS") 
  lines(1:K, re_th, lty = 2, lwd = 2, col = "darkorange")
  legend("topright", legend=c("Numerical","Theoretical"), lty=c(1,2), lwd = c(2,2),
         col=c("blue","darkorange"))
  
}
