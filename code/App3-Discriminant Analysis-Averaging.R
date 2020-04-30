library(pracma)

n = 10000
p_arr = c(20, 200, 500, 1000)
m = length(p_arr)

for (j in 1:m) {
  p = p_arr[j]
  #create a test covariance matrix
  set.seed(12345)
  v = runif(p) + 1
  Sigma = diag(v)
  iSigma = solve(Sigma)
  
  #generate a random matrix with prescribed covariance matrix
  data1 = rnorm(n*p/2, 0, 1)
  X1 = matrix(data1, n/2, p) %*% pracma::sqrtm(Sigma)$B
  data2 = rnorm(n*p/2, 1, 1)
  X2 = matrix(data2, n/2, p) %*% pracma::sqrtm(Sigma)$B
  X = rbind(X1, X2)
  
  #whole sample discriminant analysis
  X1_bar = colMeans(X1)
  X2_bar = colMeans(X2)
  Sigma_h = t(X) %*% X / n
  iSigma_h = solve(Sigma_h)
  delta2 = t(X1_bar-X2_bar) %*% iSigma_h %*% (X1_bar-X2_bar)
  
  #distributed sample discriminant analysis
  K = floor(0.8*n/(p+1))
  delta2_dist = c()
  re = c()
  for (k in 1:K) {
    m = floor(n/(2*k))  #sample size of each machine = 
    delta2_l = c()
    for (l in 1:k) {
      X_l = rbind(X1[((l-1)*m+1):(l*m),],X2[((l-1)*m+1):(l*m),])
      X_bar = colMeans(X_l)
      X1_l_bar = colMeans(X1[((l-1)*m+1):(l*m),])
      X2_l_bar = colMeans(X2[((l-1)*m+1):(l*m),])
      Sigma_l_h = t(X_l) %*% X_l / (2*m)
      iSigma_l_h = solve(Sigma_l_h)
      delta2_l[l] = t(X1_bar-X2_bar) %*% iSigma_l_h %*% (X1_bar-X2_bar)
      #delta2_l[l] = t(X1_l_bar-X2_l_bar) %*% iSigma_l_h %*% (X1_l_bar-X2_l_bar)
    }
    delta2_dist[k] = mean(delta2_l)
    re[k] = delta2_dist[k] / delta2
  }
  
  #calculate the ARE_DA
  re_th = c()
  for (k in 1:K) {
    re_th[k] = (n-p) / (n-k*p)
  }

  plot(1:K, re, type = "l", col = "blue", xlab = "Number of Machines", 
       ylab = "RE of Discriminant Analysis", main = "p/n=0.05") 
  lines(1:K, re_th, lty = 2, col = "darkorange")
  legend("topleft", legend=c("Numerical","Theoretical"), lty=c(1,2), col=c("blue","darkorange"))
}
