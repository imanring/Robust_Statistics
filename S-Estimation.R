bisquare = function(x){
  k = 1
  return((x^2/2-x^4/(2*k^2)+x^6/(6*k^4))*(abs(x)<=k)+(abs(x)>k)*k^2/6)
}
W = function(x){
  k = 1
  return((1-(x/k)^2)^2*(abs(x)<=k))
}

base_b0 = 1:500
base_b1 = 1:500
robust_b0 = 1:500
robust_b1 = 1:500
for(q in 1:500){

n = 100
x = cbind(1,rnorm(n))
y = x%*%c(5,-40) + rt(n,1)

B = coef(lm(y~x[,2]))
base_b0[q] = B[1]
base_b1[q] = B[2]
old_B = B + 1
first = TRUE
while(max(old_B-B)>.05){
  resid = y-x%*%B
  if (first){
    sigma = median(abs(resid))/.6745
  }else {
    sigma = c(wi%*%resid^2)/(n*.5)
  }
  u = resid/sigma
  if (first){
    wi = c(W(u))
    first = FALSE
  }else{
    wi = c(bisquare(u)/u^2)
  }
  w = diag(x=wi,nrow = n)
  old_B = B
  B = solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y
}
robust_b0[q]=B[1]
robust_b1[q]=B[2]
}
hist(base_b0,nclass = 30,xlab = 'Intercept',main='OLS Estimator of Intercept')
hist(robust_b0,nclass = 30, xlab='Intercept',main="Robust S-Estimator of Intercept")
hist(base_b1,nclass = 30,xlab = 'Slope',main = 'OLS Estimator of Slope')
hist(robust_b1,nclass = 30, xlab = 'Slope',main="Robust S-Estimator of Slope")

