library(pracma)

n = 10000   #训练样本量
N = 10000  #测试样本量
p_arr = c(20, 200, 500, 1000)
m = length(p_arr)
G = rep(1,N)  #实际分类
for(i in (N/2+1):N) G[i] = 0

for (j in 1:m) {
  p = p_arr[j]
  #create a covariance matrix
  #set.seed(12345)
  v = runif(p) + 1
  Sigma = diag(v)
  iSigma = solve(Sigma)
  
  #generate a random matrix with prescribed covariance matrix
  data1 = rnorm(n*p/2, 0, 1)
  X1 = matrix(data1, n/2, p) %*% pracma::sqrtm(Sigma)$B
  miu1 = rep(0, length = p)
  data2 = rnorm(n*p/2, 1, 1)
  X2 = matrix(data2, n/2, p) %*% pracma::sqrtm(Sigma)$B
  miu2 = rep(1, length = p)
  X = rbind(X1, X2)
  
  data3 = rnorm(N*p/2, 0, 1)
  test1 = matrix(data3, N/2, p) %*% pracma::sqrtm(Sigma)$B
  data4 = rnorm(N*p/2, 1, 1)
  test2 = matrix(data4, N/2, p) %*% pracma::sqrtm(Sigma)$B
  test = rbind(test1, test2)
  
  #whole sample discriminant analysis
  R = c()  #分类结果
  W = c()  #判别函数
  X_bar = colMeans(X) 
  X1_bar = colMeans(X1)
  X2_bar = colMeans(X2)
  Sigma_h = t(X) %*% X / n
  iSigma_h = solve(Sigma_h)
  sum = 0
  for (i in 1:N) {
    W[i] = t(test[i,]-X_bar) %*% iSigma_h %*% (X1_bar-X2_bar)
    if (W[i] > 0) R[i] = 1
    else R[i] = 0
    if (R[i] == G[i]) sum = sum + 1
  }
  P = sum / N
  
  #distributed sample discriminant analysis
  K = floor(0.8*n/(p+1))
  sum = rep(0, length = K)
  P_dist = c()
  re = c()
  for (k in 1:K) {
    W2 = matrix(NA, N, k) #判别函数，N行k列，一列是一台机器的结果
    R_l = matrix(NA, N, k)  #本地分类结果
    m = floor(n/(2*k))  #每半块的样本量
    for (l in 1:k) {
      X_l = rbind(X1[((l-1)*m+1):(l*m),],X2[((l-1)*m+1):(l*m),])
      X_l_bar = colMeans(X_l)
      X1_l_bar = colMeans(X1[((l-1)*m+1):(l*m),])
      X2_l_bar = colMeans(X2[((l-1)*m+1):(l*m),])
      Sigma_l_h = t(X_l) %*% X_l / (2*m)
      iSigma_l_h = solve(Sigma_l_h)
      for (i in 1:N) {
        W2[i,l] = t(test[i,]-X_l_bar) %*% iSigma_l_h %*% (X1_l_bar-X2_l_bar)
        if (W2[i,l]>0) R_l[i,l] = 1
        else R_l[i,l] = 0
      }
    }
    R_dist = c()  #汇总结果
    R = c()  #最终结果
    for (i in 1:N) {
      R_dist[i] = sum(R_l[i,])
      if (R_dist[i] > k/2) R[i] = 1
      else R[i] = 0
      if (R[i] == G[i]) sum[k] = sum[k] + 1
    }
    P_dist[k] = sum[k] / N
    re[k] = P_dist[k] / P
    print(c(k, P_dist[k]))
  }
  
  #calculate the ARE^DA
  d = sum(diag(iSigma))*sqrt(t(miu1-miu2)%*%(miu1-miu2))
  re_th = c()
  for (k in 1:K) {
    re_th[k] = 0
    if (k%%2 == 0) {
      for (i in (k/2+1):k) {
        re_th[k] = re_th[k] + choose(k,i)*((1-pnorm(-sqrt(n/(n-k*p))*d/2))^i)*(pnorm(-sqrt(n/(n-k*p))*d/2)^(k-i))/(1-pnorm(-sqrt(n/(n-p))*d/2))
      }
    }
    else {
      for (i in ((k+1)/2):k) {
        re_th[k] = re_th[k] + choose(k,i)*((1-pnorm(-sqrt(n/(n-k*p))*d/2))^i)*(pnorm(-sqrt(n/(n-k*p))*d/2)^(k-i))/(1-pnorm(-sqrt(n/(n-p))*d/2))
      }
    }
  }
  
  plot(1:K, re, type = "l", col = "blue", xlab = "Number of Machines", ylab = "RE of Discriminant Analysis", ylim = c(0, 1.2), main = "p/n=0.002") 
  lines(1:K, re_th, lty = 2, col = "darkorange")
  legend("bottomright", legend=c("Numerical","Theoretical"), lty=c(1, 2), col=c("blue","darkorange"))

}

