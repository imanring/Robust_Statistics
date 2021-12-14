#Actual values
n = 100#data points
z = 3 #number of classes
#Each column contains the parameters for a model associated with a single class
true_B = replicate(z,rnorm(2))
x = cbind(1,rnorm(n)*5)
z_sample = sample(1:z,n,replace = T)
y = replicate(n,1)
for (i in 1:n) {
  y[i]=x[i,]%*%true_B[,z_sample[i]]+rnorm(1,sd=.5)
}

#distance from a point to a line. This is a 2d case, but it can be generalized to n-dimensions
dist = function(a,b,x0,y0){
  #a and b are parameters for a line. x0 and y0 are coordinates of a point
  #minimize (x-x0)^2+(a*x+b-y0)^2
  #derivative = 2*(x-x0)+2*a*(a*x+b-y0) = 0
  #2x+2a^2*x = 2*x0-2ab+2a*y0
  x = (x0-a*b+a*y0)/(1+a^2)
  y = a*x+b
  return(sqrt((x-x0)^2+(y-y0)^2))
}
initial_B = replicate(z,rnorm(2))
dist_all = 1:z
z_guess = 1:n
plot(x[,2],y,xlab = 'x',main='Data with random estimates')
for(i in 1:z){
  abline(a=initial_B[2,i],b=initial_B[1,i])
}
old_B = initial_B+1
while(max(old_B-initial_B)>.5){
  for(j in 1:n){
    for(i in 1:z){
      #find distance from point to each line
      dist_all[i] = dist(initial_B[2,i],initial_B[1,i],x[j,2],y[j])
    }
    #Assign point to a line
    z_guess[j]=which.min(dist_all)
  }
  old_B=initial_B
  #Find new line based on new data
  for(i in 1:z){
    initial_B[,i]=coef(lm(y[z_guess==i]~x[,2][z_guess==i]))
  }
}

initial_B
true_B
plot(x[,2],y,xlab='x',main='Lines After Clustering Algorithm')
for(i in 1:z){
  abline(a=initial_B[1,i],b=initial_B[2,i],col='red')
}
plot(x[,2][z_guess==3],y[z_guess==3])
plot(x[,2][z_guess==2],y[z_guess==2])
plot(x[,2][z_guess==1],y[z_guess==1])