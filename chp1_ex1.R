set.seed(123)
# Data set up ####
library(readxl)
inflation <- read_excel("Code2017/CHAPTER1/DATA/inflation.xls")
colnames(inflation)[1] = "period"

t = nrow(inflation)
lag0 = inflation[3:t,]
lag1 = inflation[2:(t-1),]
lag2 = inflation[1:(t-2),]
x_0 = rep(1, length(lag0))
colnames(lag1)[2] = "lag1"
colnames(lag2)[2] = "lag2"
y= as.matrix(lag0[,2])
x = as.matrix(cbind(x_0, lag1[,2],lag2[,2]))
t = nrow(x)

#Set priors and starting values ####

## Priors for B 

B0 = (c(0,0,0))
sigma0 = diag(1,3,3)

## Priors for Sigma2 

T0 = 1 
D0 = 0.1

## Starting values

B = B0
sigma2 =1
reps = 5000
burn = 4000
out1 = matrix(0,reps,3)
out2 = c()

# Loop ####
for(i in 1:reps){
  
  M = solve( solve(sigma0) + drop((1 / sigma2)) * (t(x) %*% x)) %*% ( solve(sigma0) %*% B0 + drop((1 / sigma2)) *( t(x) %*% y) )
  V = solve(solve(sigma0) + drop((1 / sigma2))*(t(x) %*% x))
  chck= -1
  
  B = M + t((rnorm(3)) %*% chol(V))
  
  b= rbind(c(B[2],B[3]), c(1,0) )
  
  ee = max(abs(eigen(b)$values))
  
  resid = y - x %*% B
  
  T1 = t + T0
  D1 = D0 + t(resid)%*% resid
  
  # draw from IG
  
  z0 = rnorm(T1)
  z0_z0 = t(z0) %*% z0
  
  sigma2 = D1 /z0_z0
  out1[i,]=t(B)
  out2[i]=sigma2
  
  
} 