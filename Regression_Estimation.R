#Hubers weight function for location
W = function(x){
  k = 1.03
  return(min(1,k/abs(x)))
}
#Bisquare weight function for location
W = function(x){
  k = 4.68
  return((1-(x/k)^2)^2*(abs(x)<=k))
}
#Bisquare weight function for scale
W = function(x){
  #k is set to 1
  return(min(3-3*x^2+x^4,1/x^2))
}

bisquare = function(x){
  k = 1
  if(abs(x)>k){
    return(1)
  }else{
    return(1-(1-(x/k)^2)^3)
  }
}
MAD = function(x){
  #x is the vector of data
  #return median absolute deviation
  return (median(abs(x-median(x))))
}
MADN = function(x){
  #x is a vector of data
  #return MADN dispersion estimator
  return (MAD(x)/.675)
}


#Location robust estimator
# n = 100
# x = 16 + rnorm(n)
# for (u in 1:10){
# sum_w = 0
# sum_wx = 0
# mu = median(x)
# for(i in 1:n){
#   w = W(x[i]-mu)
#   sum_w = sum_w+w
#   sum_wx = sum_wx+w*x[i]
# }
# mu = sum_wx/sum_w
# }
# (mu)


#Scale robust estimator
# n = 100
# final = 1:100
# for(j in 1:100){
#   x = 4*rnorm(n)
#   sigma = MADN(x)
#   for(u in 1:10){
#     sum_wx2 = 0
#     sum_rho = 0
#     for(i in 1:n){
#       w = W(x[i]/sigma)
#       sum_wx2 = sum_wx2 + w*x[i]^2
#       sum_rho = sum_rho+bisquare(x[i]/sigma)
#     }
#     delta = sum_rho/n
#     sigma=sqrt(1/(n*delta)*sum_wx2)
#   }
#   final[j] = sigma
# }
# hist(final)


#Joint Location and Scale Robust Estimation
#Bisquare weight function for location
W1 = function(x){
  k = 4.68
  return((1-(x/k)^2)^2*(abs(x)<=k))
}
#Bisquare weight function for scale
W2 = function(x){
  return(min(3-3*x^2+x^4,1/x^2))
}

final_u = 1:100
final_sigma = 1:100
base_u = 1:100
base_sigma = 1:100
for(q in 1:100){
  n = 100
  p = .97
  not_outlier = rbinom(n,1,p)
  x = 59 + 25*(rnorm(n)*not_outlier+(1-not_outlier)*rnorm(n,mean=4,sd=15))
  
  u=median(x)
  sigma=MADN(x)
  for(u in 1:50){
    sum_rho = 0
    sum_w1 = 0
    sum_w1x = 0
    sum_w2r2 = 0
    for(i in 1:n){
      ri = (x[i]-u)/sigma
      w1 = W1(ri)
      w2 = W2(ri)
      sum_w1 = sum_w1 + w1
      sum_w1x = sum_w1x + w1*x[i]
      sum_w2r2 = sum_w2r2 + w2*ri^2
      sum_rho = sum_rho+bisquare(ri)
    }
    u = sum_w1x/sum_w1
    delta = 1/n*sum_rho
    sigma = sqrt(sigma^2/(n*delta)*sum_w2r2)
  }
  final_u[q] = u
  final_sigma[q] = sigma
  base_u[q] = mean(x)
  base_sigma[q] = sd(x)
}
hist(x,nclass = 30,main='Mixed Normal Data')
hist(final_sigma,xlab='Scale',main='Robust M-Estimator of Scale')
hist(base_sigma,xlab='Scale',main='Standard Deviation Estimate of Scale')
hist(final_u,xlab = 'Location',main='Robust M-Estimator of Location')
hist(base_u,xlab = 'Location',main='Sample Mean Estimator of Location')




#M-Estimation of Regression
#Bisquare weight function for location
W = function(x){
  k = 4.68
  return((1-(x/k)^2)^2*(abs(x)<=k))
}
base_b0 = 1:500
base_b1 = 1:500
robust_b0 = 1:500
robust_b1 = 1:500
for(q in 1:500){
#Generate data
n = 100
x = cbind(1,rnorm(n))
y = x%*%c(5,-4) + rt(n,1)

B = solve(t(x)%*%x)%*%t(x)%*%y
#plot(x[,2],y,xlab = 'x',main="OLS v M-Estimate")
#abline(a=B[1],b=B[2],col='red')
base_b0[q] = B[1]
base_b1[q] = B[2]
sigma = median(abs(y-x%*%B))/.675
old_B = 1+B
while(max(old_B-B)>.01){
  #sigma = median(abs(y-x%*%B))/.675
  w = c(W((y-x%*%B)/sigma))
  w = diag(x=w,nrow = n)
  old_B = B
  B = solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y
}
#abline(a=B[1],b=B[2],col='blue')
#legend(1, -18, legend=c("OLS", "Robust M"),
#       col=c("red", "blue"),lty=1)
robust_b0[q]=B[1]
robust_b1[q]=B[2]
}
hist(base_b0,nclass = 30,xlab = 'Intercept',main='OLS Estimator of Intercept')
hist(robust_b0,nclass = 30, xlab='Intercept',main="Robust M-Estimator of Intercept")
hist(base_b1,nclass = 30,xlab = 'Slope',main = 'OLS Estimator of Slope')
hist(robust_b1,nclass = 30, xlab = 'Slope',main="Robust M-Estimator of Slope")
