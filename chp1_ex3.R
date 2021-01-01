set.seed(123)
# Data set up ####
library(readxl)
inflation <- read_excel("Code2017/CHAPTER1/DATA/inflation.xls")
colnames(inflation)[1] = "period"

x_ts = as.ts(inflation)
x_ts_1 = lag(x_ts,1)
x_ts_2 = lag(x_ts,2)
x_ts_3 = lag(x_ts,3)

#bring together the new lags removing the time period 


new_x = as.ts(cbind(x_ts, x_ts_1[,2],x_ts_2[,2],x_ts_3[,2]))

# remove all rows with NAs.
new_x = new_x[complete.cases(new_x),]
x_0 = rep(1, nrow(new_x))
#remove ts
x = cbind(x_0, new_x[,3:4])
y = new_x[,2]
t = nrow(x)

#Set priors and starting values ####

## Priors for B 

B0 = (c(0,0,0))
sigma0 = diag(1,3,3)

## Priors for Sigma2 

T0 = 1 
D0 = 0.1

## Priors for rho

rho0 = 0
sigma0r = 1
## Starting values

B = B0
rho = rho0
sigma2 =1
reps = 15000
burn = 12000
out1 = matrix(0,reps,3)
out2 = c()
out3 = c()

# Loop ####
for(i in 1:reps){
  
  ystar = y[2:t] - y[1:(t-1)]*rho
  xstar = x[2:t,] - x[1:(t-1),]*drop(rho)
  
  
  M = solve( solve(sigma0) + drop((1 / sigma2)) * (t(xstar) %*% xstar)) %*% ( solve(sigma0) %*% B0 + drop((1 / sigma2)) *( t(xstar) %*% ystar) )
  V = solve(solve(sigma0) + drop((1 / sigma2))*(t(xstar) %*% xstar))
  chck= -1
  
  B = M + t((rnorm(3)) %*% chol(V))
  
  b= rbind(c(B[2],B[3]), c(1,0) )
  
  ee = max(abs(eigen(b)$values))
  
  # Step 3 compute rho
  
  y_1 =  y[2:t] -x[2:t,] %*% B
  y_2 = c(0, y_1[1:(length(y_1)-1)]) #not sure this is correct , I think it should lag in the opposite direction?
  y_1 = y_1[1:(length(y_1)-1)]
  y_2 = y_2[1:(length(y_2)-1)]
  
  MM = solve( solve(sigma0r) + drop((1 / sigma2)) * (t(y_2) %*% y_2)) %*% ( solve(sigma0r) %*% rho0 + drop((1 / sigma2)) *( t(y_2) %*% y_1) )
  VV = solve(solve(sigma0r) + drop((1 / sigma2))*(t(y_2) %*% y_2))
  
  rho = MM +t((rnorm(1) * chol(VV)))
  
  # compute residuals
  
  resid = ystar - xstar %*% B
  
  T1 = t + T0
  D1 = D0 + t(resid)%*% resid
  
  # draw from IG
  
  z0 = rnorm(T1)
  z0_z0 = t(z0) %*% z0
  
  sigma2 = D1 /z0_z0
  out1[i,] = t(B)
  out2[i] = sigma2
  out3[i] = rho 
  
} 

plot(out1[12000:15000,1], type = "l")
plot(out1[12000:15000,2], type = "l")
plot(out1[12000:15000,3], type = "l")
plot(out2[12000:15000], type = "l")
plot(out3[12000:15000], type = "l")

mean(out1[12000:15000,1])
mean(out1[12000:15000,2])
mean(out1[12000:15000,3])
