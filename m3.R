rm(list=ls(all=TRUE))
##################################################################
## read data & hyperparameters
##################################################################
library(MASS)
library(mvtnorm)
cal.pi <- function(alpha,beta,x)
{
pi_val <- 1/(1+exp(alpha+beta*x))
return(pi_val)
}
cal.poster <- function(alpha,beta,x,n,y)
{
likeli <- 1
nx <- length(x)
for (i in 1:nx)
likeli = likeli*dbinom(y[i], n[i], cal.pi(alpha,beta,x[i]), log = FALSE)
likeli = likeli*dmvnorm(c(alpha,beta),c(0,0),sqrt(100)*diag(2))
return(likeli)
}
alpha_list <- seq(from=-4,to=2,by = 0.1)
beta_list <- seq(from=-8,to=1,by = 0.1)
n1 <- length(alpha_list) ## length for alpha and beta list
n2 <- length(beta_list)
## read data
dta <- scan("/Users/zk794/Box Sync/Monte Carlo Methods in Stats/M3/sage.dta")
data = matrix(dta,4,4)
x = data[2,]
n = data[3,]
y = data[4,]
poster = matrix(1,n1,n2)
for (j in 1:n1)
for (k in 1:n2)
poster[j,k] = cal.poster(alpha_list[j],beta_list[k],x,n,y)
poster = poster/sum(poster)
#lled.contour(alpha_list, beta_list, poster,color = function(x)rev(rainbow(x)), xlab = "alpha", ylab =
"beta")
contour(alpha_list, beta_list, poster, xlab = "alpha", ylab = "beta")
######1b
###potential energy
potU = function(q){
poster = cal.poster(q[1],q[2],x,n,y)
return(-log(poster))
}
##gradient of potential energy
gradU = function(q){
tmp = exp(q[1] + q[2]*x)
grad_alpha = q[1]/100 + sum(y-n) + sum(n*tmp/(1+tmp))
grad_beta = q[2]/100 + sum(x*(y-n)) + sum(n*tmp*x/(1+tmp))
return(c(grad_alpha,grad_beta))
}
###Hamiltonian function
HMC = function (epsilon=0.01, L=100, q0)
{
q = q0
p0 = rnorm(length(q),0,1) # independent standard normal variates
p = p0
# Full steps for position and momentum
for (i in 1:L)
{
grad_U = gradU(q)
p = p - epsilon * grad_U / 2
q = q + epsilon * p
grad_U = gradU(q)
p = p - epsilon * grad_U/2
}
# Evaluate potential and kinetic energies at start and end of trajectory
U0 = potU(q0)
K0 = sum(p0^2) / 2
U1 = potU(q)
K1 = sum(p0^2) / 2
# Accept or reject the state at end of trajectory
rate = exp(U0-U1+K0-K1)
if (runif(1) < rate)
{
return (q) # accept
}
else
{
return (q0) # reject
}
}
#####initialization###
theta_b = matrix(c(rep(2*1000)),1000,2) ####storing alpha and beta
q0 = c(-5,-10)
niter = 1000
for(i in 1:niter){
####iterations
if (runif(1) < 0.5)
{
theta_b[i,] = HMC(0.01,100,q0) ###forward
}
else
{
theta_b[i,] = HMC(-0.01,100,q0) ###backward
}
q0 = theta_b[i,]
}
points(theta_b[,1],theta_b[,2],col = "red",pch=20, cex=.5)
###########################
#########1c ###############
mh.transition <- function(q0)
{
s = 10
q1 = mvrnorm(1,q0,s*diag(2))
likeli0 = cal.poster(q0[1],q0[2],x,n,y)
likeli1 = cal.poster(q1[1],q1[2],x,n,y)
acc_prob <- likeli1/likeli0
if (runif(1) < acc_prob)
{ # accept
qstar <- q1 # replace
}
else
{
qstar = q0
}
return(qstar)
}
mh <- function(n.iter=100)
{
## initialize the parameters
q0 = c(-5,-10)
qm = matrix(c(rep(NA,n.iter*2)),n.iter,2) ###store all the alpha and beta
for(iter in 1:n.iter){
q1 = mh.transition(q0)
qm[iter,] = q1
q0 = q1
}
return(qm)
}
theta_c = mh(1000)
points(theta_c[,1],theta_c[,2],col = "blue",pch=19,cex=.5)
legend("topright", c("HMC", "MH"), pch = c(20,19), col = c("red","blue"))