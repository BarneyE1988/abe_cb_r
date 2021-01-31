# Introduction #### 
library(zoo)
# This code replicates code from Blake and Mumtaz in R.

# The code estimates a VAR model using a Gibbs sampled. Beta coefficients are then used to forecast 2 years ahead.

# The Gibbs sampler algorithm does the following:
# Step1: 
library(readxl)
data_usgdp <- read_excel("Code2017/CHAPTER2/DATA/DATAIN.XLS", col_names = FALSE)
N = ncol(data_usgdp) # number of series?
L = 2 # number of lags
Y = data_usgdp
t = nrow(inflation)
lag0 = data_usgdp[3:t,]
lag1 = data_usgdp[2:(t-1),]
lag2 = data_usgdp[1:(t-2),]
x_0 = rep(1, length(lag0))
colnames(lag1)[2] = "lag1"
colnames(lag2)[2] = "lag2"

x = as.matrix(cbind(x_0, lag1, lag2))
t = nrow(x)
y= as.matrix(Y[(3:nrow(Y)),])
check_data = cbind(y,x)

#### Use OLS for priors ####

lm(y ~x)

y_1 = y[,1]
x_1 = x[,1:2]
b0 = solve(t(x_1) %*% x_1)  %*% (t(x_1) %*%y_1)
s1 = sqrt((t( y_1 - x_1 %*% b0) %*% ( y_1 - x_1 %*% b0))/ (length(y_1)-2))


y_2 = y[,2]
x_2 = x[,c(1,3)]
b0 = solve(t(x_2) %*% x_2)  %*% (t(x_2) %*% y_2)
s2 = sqrt((t( y_2 - x_2 %*% b0) %*% ( y_2 - x_2 %*% b0))/ (length(y_2)-2))


#### Param for the minnesota priors

lambda_1 = 1
lambda_2 = 1
lambda_3 = 1
lambda_4 = 1

# Prior means of the coefficients of the two equations of the VAR

B0_1 = as.matrix(c(0,1,0,0,0),nrow=1,ncol =5)
B0_2 = as.matrix(c(0,0,1,0,0),nrow=1,ncol =5)
B0 = as.matrix(c(B0_1,B0_2))

#prior variance of vec(B)

H = matrix(0, nrow = 10, ncol = 10)
H[1,1] = (s1 *lambda_4)^2
H[2,2] = (lambda_1)^2
H[3,3] = ((s1 * lambda_1 * lambda_2) / s2)^2
H[4,4] = (lambda_1/(2^lambda_3))^2
H[5,5] = ((s1 * lambda_1 *lambda_2) / (s2 * (2^lambda_3))) ^2

H[6,6] = (s2 *lambda_4)^2
H[7,7] = ((s2 * lambda_1 * lambda_2) / s1)^2
H[8,8] = (lambda_1)^2
H[9,9] = ((s2 * lambda_1 *lambda_2) / (s1 * (2^lambda_3))) ^2
H[10,10] = (lambda_1/(2^lambda_3))^2

# prior scale matrix for sigma the VAR covariance

s= diag(N)
alpha = N+1
Sigma= diag(N)
betaols= as.matrix(c( solve(t(x) %*% x) %*%(t(x) %*%y)), nrow =6, ncol =1 )

Reps = 10000
burn = 5000
#preallocate 
beta=matrix(0,10, Reps)
out1=matrix(0,Reps,nrow(y)+12)
out2=matrix(0,Reps,nrow(y)+12)
i =1
for(j in 1:Reps){
  
  M = solve( solve(H) + kronecker(solve(Sigma),t(x)%*%x)  ) %*% (solve(H) %*% B0 + kronecker( solve(Sigma), t(x) %*%x) %*%betaols)
  V = solve(solve(H) + kronecker(solve(Sigma),t(x)%*%x) )
  beta[,j] = M + t(rnorm(N * (N * L +1)) %*% chol(V))
  e = y - x %*% matrix( beta[,j], nrow= N*L +1, ncol= N)
  scale = t(e) %*% e + s
  k_s = nrow(solve(scale))
  z_s = matrix(0, t +alpha, k_s)
  mu_s = matrix(0,k_s,1)
  for( i in 1:(t+alpha)){
    z_s[i,] = (t(chol(solve(scale))) %*% rnorm(k_s))
  }
  Sigma = solve(t(z_s) %*% z_s)
  
 
    yhat = matrix(0,14,2)
    yhat[1:2,] = y[(nrow(y)-1):nrow(y),]
    for( i in 3:14){
      yhat[i,] = as.matrix(cbind(1, t(yhat[i-1,]), t(yhat[i-2,]))) %*% matrix( beta[,j], N*L+1,N) + t(rnorm(N)) %*% chol(Sigma)
    }
    out1[j,] = (c(y[,1],yhat[3:nrow(yhat),1]))
    out2[j,] = (c(y[,2],yhat[3:nrow(yhat),2]))   
  
  
}

quants <- c(0,0.05,0.25,0.50,0.75,0.90,0.95,0.99,1)
result <- t(apply( out1 , 2 , quantile , probs = quants , na.rm = TRUE ))
period <- as.yearqtr(1948 + seq(3, 4*(2014-1948))/4)
result <- cbind(period, result)
plot(result[,1], result[,9], type="l")

#check the convergence of beta
for(i in 1:10){
plot(cumsum(beta[i,])/1:Reps,type="l")
}