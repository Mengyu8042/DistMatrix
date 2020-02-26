library(pracma)

N = 10000
p_arr = c(20, 200, 500, 1000)
#p_arr = c(10, 20, 50, 100)
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
  X_bar = colMeans(X)
  
  #find the full-sample covariance matrix of beta_hat
  Sigma_h = t(X) %*% X / N
  iSigma_h = solve(Sigma_h)
  T2 = (N-1) * t(X_bar) %*% iSigma_h %*% X_bar

  #find the distributed covariance matrix of beta_hat
  K = floor(N/(p+1))
  dist_iSigma_h = list()
  dist_T2 = rep(0, length = K)
  
  for (k in 1:K) {
    n = floor(N/k)
    Sigma_h_k = list()
    iSigma_h_k = list()
    dist_iSigma_h[[k]] = matrix(0, p, p)
    X_bar_k = list()
 
    for (j in 1:k) {
      Sigma_h_k[[j]] = t(X[((j-1)*n+1):(j*n), ]) %*% X[((j-1)*n+1):(j*n), ] / n
      iSigma_h_k[[j]] = solve(Sigma_h_k[[j]])
      dist_iSigma_h[[k]] = dist_iSigma_h[[k]] + iSigma_h_k[[j]] / k
    }
    for (i in 1:k) {
      X_bar_k[[i]] = colMeans(X[((i-1)*n+1):(i*n), ])
      #dist_T2[k] = dist_T2[k] + (n-1)/k * t(X_bar_k[[i]]) %*% iSigma_h_k[[i]] %*% X_bar_k[[i]]
      dist_T2[k] = dist_T2[k] + (n-1)/k * t(X_bar) %*% iSigma_h_k[[i]] %*% X_bar
    }

  }
  
  #calculate and print relative efficieny
  re = c()
  re_th = c()
  for (k in 1:K) {
    re[k] =  dist_T2[k] / T2
    re_th[k] = (N-k)*(N-p) / (k*(N-1)*(N-k*p))
    print(c(p, k, re[k], re_th[k]))
  }
  
  #plot: relative efficieny of hypothesis test ~ number of machines
  plot(1:K, re, type = "l", lwd = 2, col = "blue", xlab = "Number of Machines", 
       ylab = "Relative Efficiency of Test") 
  lines(1:K, re_th, lty = 2, lwd = 2, col = "darkorange")
  legend("topright", legend=c("Numerical","Theoretical"), lty=c(1,2), lwd = c(2,2), 
         col=c("blue","darkorange"))
  
}

#------------------------------------------------
N = 10000
p_arr = c(20, 200, 500, 1000)
L = length(p_arr)

for (l in 1:L) {
  p = p_arr[l]
  K = floor(N/(p+1))
  re_loss_th = c()
  for (k in 1:K) {
    re_loss_th[k] = (N-k)*(N-p) / (k*(N-1)*(N-k*p))
    print(re_loss_th[k])
  }
  plot(1:K, re_loss_th, type = "l") 
}
