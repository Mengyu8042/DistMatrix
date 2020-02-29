#X = as.matrix(rbind(iris[71:80,1:4], iris[111:120,1:4]))
#test = as.matrix(rbind(iris[51:70,1:4],iris[81:100,1:4], iris[101:110,1:4],iris[121:150,1:4]))

new = sample(100, 100, replace = FALSE)
data = as.matrix(iris[51:150,1:4])
newdata = matrix(NA, 100, 4)
for (i in 1:100) {
  newdata[i,] = data[new[i],]
}
X = as.matrix(rbind(newdata[1:15,1:4], newdata[51:65,1:4]))
test = as.matrix(rbind(newdata[16:50,1:4], newdata[66:100,1:4]))
n = length(X[,1])  #训练样本量n=30
N = length(test[,1])  #测试样本量N=70
G = rep(1,N)  #实际分类
for(i in (N/2+1):N) G[i] = 0

#whole sample discriminant analysis
R = c()  #分类结果
W = c()  #判别函数
X_bar = colMeans(X) 
X1_bar = colMeans(X[1:n/2,])
X2_bar = colMeans(X[(n/2+1):n,])
Sigma_h = t(X) %*% X / n
iSigma_h = solve(Sigma_h)
sum = 0
for (i in 1:N) {
  W[i] = t(test[i,]-X_bar)%*%iSigma_h%*%(X1_bar-X2_bar)
  if (W[i]>0) R[i] = 1
  else R[i] = 0
  if (R[i]==G[i]) sum =sum + 1
}
(P = sum / N)

#distributed sample discriminant analysis
K = floor(n/(4+1))
sum = rep(0, length = K)
P_dist = c()
re = c()
for (k in 1:K) {
  W2 = matrix(NA, N, k) #判别函数，N行k列，一列是一台机器的结果
  R_l = matrix(NA, N, k)  #分类结果
  m = floor(n/(2*k))  #每半块的样本量
  for (l in 1:k) {
    X_l = rbind(X[((l-1)*m+1):(l*m),],X[((n/2+1)+(l-1)*m):(n/2+l*m),])
    X_bar = colMeans(X_l)
    X1_bar = colMeans(X[((l-1)*m+1):(l*m),])
    X2_bar = colMeans(X[((n/2+1)+(l-1)*m):(n/2+l*m),])
    Sigma_h = t(X_l) %*% X_l / (2*m)
    iSigma_h = solve(Sigma_h)
    for (i in 1:N) {
      W2[i,l] = t(test[i,]-X_bar)%*%iSigma_h%*%(X1_bar-X2_bar)
      if (W2[i,l]>0) R_l[i,l] = 1
      else R_l[i,l] = 0
    }
  }
  R_dist = c()
  R = c()  #最终结果
  for (i in 1:N) {
    R_dist[i] = sum(R_l[i,])
    if (R_dist[i]>k/2) R[i] = 1
    else R[i] = 0
    if (R[i]==G[i]) sum[k] =sum[k] + 1
  }
  P_dist[k] = sum[k] / N
  re[k] = P_dist[k] / P
  print(c(k, P_dist[k]))
}

plot(1:K, re, type = "l", lwd = 2, col = "blue", xlab = "Number of Machines", 
     ylab = "RE of Discriminant Analysis") 
legend("topright", legend=c("Numerical"), lty=1, lwd = 2, col="blue")

