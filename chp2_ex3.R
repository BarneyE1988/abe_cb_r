# Introduction #### 
library(zoo)
# This code replicates code from Blake and Mumtaz in R.

#The code estimate a bi-variate VAR with Minnesota prior and Gibbs sampling and steady state priors

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

x = as.matrix(cbind(lag1, lag2,x_0))
t = nrow(x)
y= as.matrix(Y[(3:nrow(Y)),])
check_data = cbind(y,x)

# Compute the standard deviation of each series residual via an OLS regression to be used in setting the prior

# first variable

y_1 = y[,1]
x_1 = x[,c(1,3)]
b0 = solve(t(x_1) %*% x_1)  %*% (t(x_1) %*%y_1)
s1 = sqrt((t( y_1 - x_1 %*% b0) %*% ( y_1 - x_1 %*% b0))/ (length(y_1)-2))


y_2 = y[,2]
x_2 = x[,c(2,3)]
b0 = solve(t(x_2) %*% x_2)  %*% (t(x_2) %*% y_2)
s2 = sqrt((t( y_2 - x_2 %*% b0) %*% ( y_2 - x_2 %*% b0))/ (length(y_2)-2))


#### Param for the minnesota priors

lambda_1 = 1
lambda_2 = 1
lambda_3 = 1
lambda_4 = 1


# Prior means of the coefficients of the two equations of the VAR

B0_1 = as.matrix(c(1,0,0,0),nrow=1,ncol =4)
B0_2 = as.matrix(c(0,1,0,0),nrow=1,ncol =4)
B0 = as.matrix(c(B0_1,B0_2))


H = matrix(0, nrow = 8, ncol = 8)

H[1,1] = (lambda_1)^2
H[2,2] = ((s1 * lambda_1 * lambda_2) / s2)^2
H[3,3] = (lambda_1/(2^lambda_3))^2
H[4,4] = (( s1 * lambda_1 * lambda_2) / (s2 * (2^lambda_3)))^2

H[5,5] = ((s2 * lambda_1 * lambda_2) / s1)^2
H[6,6] = (lambda_1)^2
H[7,7] = (( s2 * lambda_1 * lambda_2) / (s1 * (2^lambda_3)))^2
H[8,8] = (lambda_1/(2^lambda_3))^2


s= diag(N)
alpha = N+1

#Prior mean
M0 = c(1,1)
V0 = diag(N)*0.001 #prior variance

#Starting value via OLS
betaols=  solve(t(x) %*% x) %*%(t(x) %*%y)
f = rbind(t(betaols[1:(N*L),]), diag(1,N*(L-1), N*L))

C = matrix(0,nrow = nrow(f),1)
C[1:N,] <- t(betaols[N*L+1,])
MU= solve(diag(1, nrow = nrow(f))-f)%*% C

e = as.matrix(Y - x%*%betaols)

Sigma = (t(e) %*% e)/t

Reps= 10000
burn = 5000
out1=matrix(0,Reps,294)
out2=matrix(0,Reps,294)

# start loop 

i =1
for(j in 1:Reps){
  
Y0 = data_usgdp - kronecker(matrix(1,nrow = nrow(data_usgdp),1), t(MU[1:N,]))
t=nrow(Y0)
  X0_lag0 = Y0[3:t,]
  X0_lag1 = Y0[2:(t-1),]
  X0_lag2 = Y0[1:(t-2),]
  X0 = as.matrix(cbind(X0_lag1, X0_lag2))
  Y0= as.matrix(Y0[(3:nrow(Y0)),])
  check_data = cbind(Y0,X0)
  
bols =    as.vector(solve(t(X0) %*% X0) %*%(t(X0) %*%Y0))
M = solve( solve(H) + kronecker(solve(Sigma),t(X0) %*%X0)) %*% ( solve(H)%*%B0 + kronecker(solve(Sigma),t(X0)%*%X0)%*%bols)
V = solve( solve(H) + kronecker(solve(Sigma),t(X0) %*%X0))
beta = M + t(t(rnorm(N*(N*L))) %*%chol(V))
beta_1 = matrix(beta,N*L, N )

e = Y0 - X0 %*% beta_1

scale = t(e) %*% e + s

k_s = nrow(solve(scale))
z_s = matrix(0, t +alpha, k_s)
mu_s = matrix(0,k_s,1)
for( i in 1:(t+alpha)){
  z_s[i,] = (t(chol(solve(scale))) %*% rnorm(k_s))
}
Sigma = solve(t(z_s) %*% z_s)

Y1 = y - x[,(1:(ncol(x)-1))]%*%beta_1
U= diag(1, N,N)
U=rbind(U,beta_1)
t=nrow(Y1)
D=cbind(matrix(1,t,1),- matrix(1,t,L))
vstar1 = solve(t(U)%*% kronecker(t(D)%*%D,solve(Sigma))%*%U + solve(V0))
mstar1 = vstar1 %*%(t(U)%*%as.vector(solve(Sigma)%*%t(Y1)%*%D) + solve(V0)%*%matrix(M0))

MU = mstar1 + t(rnorm(N)%*%chol(vstar1))

#forecast

f= rbind( t(beta_1[(1:(N*L)),]),diag(1,N*(L-1),N*L))   
mu=numeric()

for(i in 1:L){
  mu=rbind(mu,MU)
}

C= (diag(1,nrow(f))-f)%*%mu
Y=as.matrix(Y)
yhat = matrix(0,44,2)
yhat[1:2,1:2] = as.matrix(Y[(nrow(Y)-1):nrow(Y),])

for (i in 3:44) {
yhat[i,]=t(C[1:N]) +  c(yhat[i-1,], yhat[i-2,]) %*% matrix( beta, N*L,N) + rnorm(N)%*%chol(Sigma)  
}

out1[j,]=t((c(Y[,1],yhat[(3:nrow(yhat)),1])))
out2[j,]=((c(Y[,2],yhat[(3:nrow(yhat)),2])))
}




quants <- c(0.50,0.1,0.2,0.3,0.7,0.8,0.9)
result <- t(apply( out2[burn:Reps,] , 2 , quantile , probs = quants , na.rm = TRUE ))
period <- 1:294


out2[244,]
result <- cbind(period, result)
plot(result)