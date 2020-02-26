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
  
  #find the full-sample inverse covariance matrix
  Sigma_h = t(X) %*% X / N
  iSigma_h = solve(Sigma_h)
  
  #find the distributed inverse covariance matrix
  dist_iSigma_h = list()
  K = floor(N/(p+1))
  for (k in 1:K) {
    n = floor(N/k)
    Sigma_h_k = list()
    iSigma_h_k = list()
    dist_iSigma_h[[k]] = matrix(0, p, p)
    
    for(j in 1:k) {
      Sigma_h_k[[j]] = t(X[((j-1)*n+1):(j*n), ]) %*% X[((j-1)*n+1):(j*n), ] / n
      iSigma_h_k[[j]] = solve(Sigma_h_k[[j]])
      dist_iSigma_h[[k]] = dist_iSigma_h[[k]] + iSigma_h_k[[j]] / k
    }
  }
  
  #calculate and print relative efficiency
  phi = 0
  for (i in 1:p) {
    phi = phi + (iSigma_h[i,i] / iSigma[i,i]) / p
  }
  
  dist_phi = rep(0, length = K)
  re = c()
  re_th = c()
  for (k in 1:K) {
    for (i in 1:p) {
      dist_phi[k] = dist_phi[k] + (dist_iSigma_h[[k]][i,i] / iSigma[i,i]) / p
    }
    re[k] = phi/dist_phi[k]
    re_th[k] = (N - k * p) / (N - p)
    print(c(p, k, re[k], re_th[k]))
  }

  #plot: relative efficiency ~ number of machines
  plot(1:K, re, type = "l", lwd = 2, col = "blue", xlab = "Number of Machines", 
       ylab = "Relative Efficiency") 
  lines(1:K, re_th, lty = 2, lwd = 2, col = "darkorange")
  legend("topright", legend=c("Numerical","Theoretical"), lty=c(1,2), lwd = c(2,2), 
         col=c("blue","darkorange"))
  
}
