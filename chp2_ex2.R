# Introduction #### 
library(zoo)
# This code replicates code from Blake and Mumtaz in R.

# The code estimates a VAR model using a Gibbs sampled. Beta coefficients are then used to forecast 2 years ahead.

# The Gibbs sampler algorithm does the following:
# Step1: 
library(readxl)
data_usgdp_2 <- read_excel("Code2017/CHAPTER2/DATA/dataUS.XLS", col_names = FALSE)
N = ncol(data_usgdp_2) # number of series?
L = 2 # number of lags
Y = data_usgdp_2
t = nrow(data_usgdp_2)
lag0 = data_usgdp_2[3:t,]
lag1 = data_usgdp_2[2:(t-1),]
lag2 = data_usgdp_2[1:(t-2),]
x_0 = rep(1, length(lag0))
colnames(lag1)[2] = "lag1"
colnames(lag2)[2] = "lag2"


x = as.matrix(cbind(x_0, lag1, lag2))
t = nrow(x)
y= as.matrix(Y[(3:nrow(Y)),])
check_data = cbind(y,x)



y_1 = y[,1]
x_1 = x[,1:2]
b0 = solve(t(x_1) %*% x_1)  %*% (t(x_1) %*%y_1)
s1 = sqrt((t( y_1 - x_1 %*% b0) %*% ( y_1 - x_1 %*% b0))/ (length(y_1)-2))


y_2 = y[,2]
x_2 = x[,c(1,3)]
b0 = solve(t(x_2) %*% x_2)  %*% (t(x_2) %*% y_2)
s2 = sqrt((t( y_2 - x_2 %*% b0) %*% ( y_2 - x_2 %*% b0))/ (length(y_2)-2))

y_3 = y[,1]
x_3 = x[,c(1,4)]
b0 = solve(t(x_3) %*% x_3)  %*% (t(x_3) %*%y_3)
s3 = sqrt((t( y_3 - x_3 %*% b0) %*% ( y_3 - x_3 %*% b0))/ (length(y_3)-2))


y_4 = y[,2]
x_4 = x[,c(1,5)]
b0 = solve(t(x_4) %*% x_4)  %*% (t(x_4) %*% y_4)
s4 = sqrt((t( y_4 - x_4 %*% b0) %*% ( y_4 - x_4 %*% b0))/ (length(y_4)-2))


#### Param for the minnesota priors

lambda_1 = 0.1
lambda_2 = 1
lambda_3 = 0.03
lambda_4 = 1

B0 = matrix(0, N*L+1,N)
for(i in 1:N){
  B0[i+1,i]=0.95
}
B0 = as.vector(B0)

# Prior variance of vec(B)
H = diag(1, N*(N*L+1),N*(N*L+1))

H[3,3] = 1e-9
H[4,4] = 1e-9
H[5,5] = 1e-9
H[7,7] = 1e-9
H[8,8] = 1e-9
H[9,9] = 1e-9

H[1,1] = (s1 * lambda_4)^2
H[2,2] = lambda_1^2
H[6,6] = (lambda_1 / (2^lambda_3))^2

H[10,10] = (s2 * lambda_4)^2
H[11,11] = ((s2*lambda_1)/s1)^2
H[12,12] = lambda_1^2
H[13,13] = ((s2 * lambda_1)/s3)^2
H[14,14] = ((s2 * lambda_1)/s4)^2
H[15,15] = ((s2 * lambda_1)/ (s1 * (2^lambda_3)))^2
H[16,16] = (lambda_1/(2^lambda_3))^2
H[17,17] = ((s2*lambda_1) / (s3 * (2^lambda_3)))^2
H[18,18] = ((s2*lambda_1) / (s4 * (2^lambda_3)))^2


H[19,19] = (s3 * lambda_4)^2
H[20,20] = ((s3*lambda_1)/s1)^2
H[21,21] = ((s3 * lambda_1)/s2)^2
H[22,22] = lambda_1^2
H[23,23] = ((s3 * lambda_1)/s4)^2
H[24,24] = ((s3 * lambda_1)/ (s1 * (2^lambda_3)))^2
H[25,25] = ((s3*lambda_1) / (s2 * (2^lambda_3)))^2
H[26,26] = (lambda_1/(2^lambda_3))^2
H[27,27] = ((s3*lambda_1) / (s4 * (2^lambda_3)))^2



H[28,28] = (s4 * lambda_4)^2
H[29,29] = ((s4*lambda_1)/s1)^2
H[30,30] = ((s4 * lambda_1)/s2)^2
H[31,31] = ((s4 * lambda_1)/s3)^2
H[32,32] = lambda_1^2
H[33,33] = ((s4 * lambda_1)/ (s1 * (2^lambda_3)))^2
H[34,34] = ((s4*lambda_1) / (s2 * (2^lambda_3)))^2
H[35,35] = ((s4*lambda_1) / (s3 * (2^lambda_3)))^2
H[36,36] = (lambda_1/(2^lambda_3))^2

s=(diag(1,N,N))
alpha = N+1

Sigma = (diag(1,N,N))
betaols= as.vector(solve(t(x) %*% x) %*%(t(x) %*%y))



Reps = 40000
burn = 5000

beta=as.numeric(Reps)
out1=matrix(0,Reps,58)
out2=matrix(0,Reps,58)
out3=matrix(0,Reps,58)
out4=matrix(0,Reps,58)


i =1
for(j in 1:Reps){
  
  M = solve( solve(H) + kronecker(solve(Sigma),t(x)%*%x)  ) %*% (solve(H) %*% B0 + kronecker( solve(Sigma), t(x) %*%x) %*%betaols)
  V = solve(solve(H) + kronecker(solve(Sigma),t(x)%*%x) )
  beta = M + t(rnorm(N * (N * L +1)) %*% chol(V))
  e = y - x %*% matrix( beta, nrow= N*L +1, ncol= N)
  scale = t(e) %*% e + s
  k_s = nrow(solve(scale))
  z_s = matrix(0, t +alpha, k_s)
  mu_s = matrix(0,k_s,1)
  for( i in 1:(t+alpha)){
    z_s[i,] = (t(chol(solve(scale))) %*% rnorm(k_s))
  }
  Sigma = solve(t(z_s) %*% z_s)
 A0 = chol(Sigma)
 v = matrix(0,60,N)
 v[L+1,2]=-1
  yhat = matrix(0,60,N)
  
  for( i in 3:60){
    yhat[i,] = as.matrix(cbind(0, t(yhat[i-1,]), t(yhat[i-2,]))) %*% matrix( beta, N*L+1,N) + v[i,] %*%A0
  }
  out1[j,] = (c(yhat[3:nrow(yhat),1]))
  out2[j,] = (c(yhat[3:nrow(yhat),2]))   
  out3[j,] = (c(yhat[3:nrow(yhat),3]))
  out4[j,] = (c(yhat[3:nrow(yhat),4]))   
  
}

quants <- c(0.50,0.16,0.84)
result <- t(apply( out2[burn:Reps,] , 2 , quantile , probs = quants , na.rm = TRUE ))
period <- 1:58
result <- cbind(period, result)

plot(0,0,xlim = c(0,60),ylim = c(-0.0002,0.0002),type = "n")
cl <- rainbow(5)
invisible(lapply(2:4, function(i) lines(result[,i],col = cl[i-1],type = 'l')))


# check the convergence of beta
for(i in 1:10){
  plot(cumsum(beta[i,])/1:Reps,type="l")
}



