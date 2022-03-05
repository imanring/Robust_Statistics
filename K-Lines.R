#distance from a point to a line. This is a 2d case, but it can be generalized to n-dimensions
dist_2d = function(a,b,x0,y0){
  #a and b are parameters for a line. x0 and y0 are coordinates of a point
  #minimize (x-x0)^2+(a*x+b-y0)^2
  #derivative = 2*(x-x0)+2*a*(a*x+b-y0) = 0
  #2x+2a^2*x = 2*x0-2ab+2a*y0
  x = (x0-a*b+a*y0)/(1+a^2)
  y = a*x+b
  return(sqrt((x-x0)^2+(y-y0)^2))
}
#Generalized Distance function
dist = function(B,r0){
  #B is the parameter vector for the linear regression
  #r0 is the point formatted (y,x1,x2,x3,x4,...)
  #vector normal to the regression plane
  n = c(1,-B[2:length(B)])
  return(abs(B[1]-n%*%r0)/sqrt(n%*%n))
}


#Generate data
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


#find initial guesses
initial_B = replicate(z,rnorm(2))
dist_all = 1:z
z_guess = 1:n

#Plot random guesses
#plot(x[,2],y,xlab = 'x',main='Data with random estimates')
#for(i in 1:z){
#  abline(a=initial_B[2,i],b=initial_B[1,i])
#}

old_B = initial_B+1
count = 0
while(max(old_B-initial_B)>.05){
  for(j in 1:n){
    for(i in 1:z){
      #find distance from point to each line
      dist_all[i] = dist(initial_B[,i],c(y[j],x[j,2]))
    }
    #Assign point to a line
    z_guess[j]=which.min(dist_all)
  }
  count = count + 1
  old_B=initial_B
  #Find new line based on new data
  for(i in 1:z){
    mod = lm(y[z_guess==i]~x[,2][z_guess==i])
    #Add noise to the estimates based on how well the algorithm is doing 
    #in order to search the space parameter space
    decay_factor = 1.5
    initial_B[,i]=coef(mod)+rnorm(2,sd=(1-summary(mod)$r.squared)/(decay_factor^count))
  }
  if (any(is.na(initial_B))){
    print(q)
  }
}
pick = 1:z
tot = 0
for(i in 1:z){
  difs = 1:z
  for(j in 1:z){
    difs[j] = sum((initial_B[,i]-true_B[,j])^2)
  }
  if (i!=1){
    pick[i] = which(difs==min(difs[-pick[1:(i-1)]]))
  }else{
    pick[i] = which.min(difs)
  }
  tot = tot + difs[pick[i]]
}

#check clustering algorithm
plot(x[,2],y,xlab='x',main='Lines After Clustering Algorithm')
for(i in 1:z){
 abline(a=initial_B[1,i],b=initial_B[2,i],col='red')
}


for(i in 1:z){
 abline(a=true_B[1,i],b=true_B[2,i],col='blue')
}

#plot lines points assigned to a line
plot(x[,2][z_guess==3],y[z_guess==3])
plot(x[,2][z_guess==2],y[z_guess==2])
plot(x[,2][z_guess==1],y[z_guess==1])
