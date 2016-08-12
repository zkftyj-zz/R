
require(gtools)  # I found a Dir r.v. generator in this package
## if needed install it with
## install.packages("gtools")  

## DATA
rm(list=ls(all=TRUE))
##################################################################
## read data & hyperparameters
##################################################################
library(MASS)

## hyper parameters
Kmx <- 15                # max number of frequencies
b0 <- rep(0,2*Kmx+1)    # p(beta) = N(b0, A0^-1) (A0 precision matrix)
A0 <- diag(1,2*Kmx+1)
lambda <- 7
asig <- 1               # 1/sig2 ~ Ga(asig, bsig)
bsig <- 1

## read data
dta <- read.table("C:/Users/zk794/Box Sync/DISC/Austin Courses/2016 Spring/Monte Carlo Methods in Stats/M2/cepheid.dta",header=T,skip=15)
o <- order(dta$phase)   # sort the data points, just for easier plotting
dta <- dta[o,]
n <- nrow(dta)          # sample size
x <- dta$phase          # phase
y <- dta$velocity       # velocity

make.x <- function(xx,k=1)
  { ## make columns in the design matrix for k-th harmonic
    tm <- xx*(k*2*pi)
    return(cbind(sin(tm),cos(tm)))
  }

## make design matrix for up to Kmx trig polys
## we will use it later by selecting the first (1+2K) columns 
X <- rep(1,n)
for(k in 1:Kmx)
  X <- cbind(X,make.x(x,k))  # see make.x(.) below

## design matrix for plotting fitted curve on a grid
n0 <- 100
X0 <- rep(1,n0)
x0 <- seq(from=0,to=1.1,length=n0)
for(k in 1:Kmx)
  X0 <- cbind(X0,make.x(x0,k))
##################################################################
 plots
##################################################################
plt.dta <- function(plt.spline=F)
  { # plots data and adds smoothing spline (if plt.spine=T)
    plot(x,y, pch=19,bty="l",xlab="PHASE",ylab="VELOCITY",xlim=c(0,1.1),
         ylim=c(-25,50))
    if (plt.spline){ # add smoothing spline
      fit <- smooth.spline(x,y)
      lines(fit,col=3,type="l",lwd=2,lty=3)
    }
  }

########################Gibbs transition probs

sample.b <- function(K,sig2)
  { # generate b ~ p(b | K, sig2, y)
    idx <- 1:(2*K+1)  # select columns (elements) for
    ## K harmonics
    Xk <- X[,idx]
    
    V = solve(0.1 * diag(2*K+1) + (1/sig2) * t(Xk)%*%Xk)
    M = (1/sig2) * V %*% (t(Xk)%*%y)
    b = mvrnorm(1, M, V)
    return(b)
  }

sample.sig2 <- function(K,b)
  { # generate 1/sig2 ~ p(1/sig2 | K,b,y)
  	idx  <- 1:(2*K +1)
    Xk <- X[,idx]
    tmp = sum((y-Xk%*%b)^2)/2

    sig2 <- 1/rgamma(1, shape = 1+n/2,rate = 1+tmp)
    return(sig2)
  }

##################################################################
## RJ

rj <- function(K,b,sig2)
  { # RJ move.
    ## returns th=list(b=b,K=K)
    q <- qbirth(K)       # prob of move up ("birth") - see below
    u <- runif(1)        # flip coin 
    if (u < q)           #      birth
      th <- rj.birth(K,b,sig2)
    else                 #      death
      th <- rj.death(K,b,sig2)
    return(th) 
  }

rj.birth <- function(K,b,sig2)
  { ## birth move -- add one harmonic
    
    ## 1. generate auxiliaries (u1,u2) for the new regression coefficients
    ##    we use a normal linear regression of the residual on the
    ##    (K+1)-st harmonics
    u <- rnorm(2)
    idx <- 1:(2*K+1)
    Xk <- X[,idx]
    xi = X[,c(2*K+2,2*K+3)]
    epsi = y - Xk%*%b

    ## 2. b1 = T(b,u); save Jacobian = |L| in J.
    Hinv = solve(t(xi)%*%xi)
    bhat = Hinv%*%t(xi)%*%epsi
    ##sighat2 = sum((epsi - xi%*%bhat)^2)/(n-2)
    Shat = Hinv*sig2
    L = t(chol(Shat))
    J = abs(det(L))
    
    baux = bhat + L%*%u 
    b1 = c(b,baux)

    ## 2. acc prob
    r <- rho(K,b,u,b1,J,sig2)
    coin <- runif(1)
    ## 3. accept (or not)
    if (coin < r){ # accept with pr min(1,rho)
      b <- b1
      K <- K+1
    }
    ## else reject (do nothing :-)
    return(list(b=b,K=K))
  }

rj.death <- function(K1,b1,sig2)
  { ## death move -- drop last harmonic
    ## NOTE: it is convenient for notation to label now
    ##    current pars     K1,b1
    ##    proposed pars    K, b
    ## this keeps the notation compatible with the birth move
    ## with b1 being the larger model vector

    ## 1. T inv mapping (b,u) = Tinv(b1)
    ## and save J=Jacobian of mapping T(b,u)
    K = K1-1
    idx <- 1:(2*K+1)
    b = b1[idx]
    
    Xk = X[,idx]
    xi = X[,c(2*K+2,2*K+3)]
    epsi = y - Xk%*%b
    Hinv = solve(t(xi)%*%xi)
    bhat = Hinv%*%t(xi)%*%epsi
    Shat = Hinv*sig2
    L = t(chol(Shat))
    J = abs(det(L))
    u = solve(L)%*%(b1[c(2*K+2,2*K+3)]-bhat)

    ## 2. acc ratio (for opposite birth move)
    r <- rho(K,b,u,b1,J,sig2)
    ## 3. accept (or not)
    coin = runif(1)
    if (coin < 1/r){  # accept with pr min(1, 1/rho)
      b1 <- b
      K1 <- K
    }
    ## else reject (do nothing :-)
    return(list(b=b1,K=K1))
}    

qbirth <- function(K)
  { # prob of proposing a birth move = 0.5, except when K=1
    return(ifelse(K==1,1,0.5))
  }

rho <- function(K,b,u,b1,J,sig2)
  { ##  acceptance ratio for birth move,
    ##  moving from b -> (b,u)

    if(length(b) != 2*K+1)        # check
      cat("\n *** Error 1: ",
          "rho(.) should be called with rho(K,b1,b,u,J).\n")
    if(length(b1) != length(b)+2) # check
      cat("\n *** Error 2: ",
          "rho(.) should be called with rho(K,b1,b,u,J).\n")
    K1 <- K+1

    r = exp(-(1/(2*sig2)) *(sum((y-X[,1:(2*K+3)]%*%b1)^2)-sum((y-X[,1:(2*K+1)]%*%b)^2))) * (lambda/K1) * dnorm(b1[2*K+2],mean=0,sd=sqrt(10)) *
dnorm(b1[2*K+3],mean=0,sd=sqrt(10)) * (1-qbirth(K+1)) * J / (qbirth(K) * dnorm(u[1],mean=0,sd=1) * dnorm(u[2],mean=0,sd=1))
    return(r)
  }

##################################################################
## MCMC
##################################################################

mcmc <- function(niter=100)
  { # main loop for mcmc
    K <- 4     # initial values
    sig2 <- 1
    ## initize lists to save imputed par values and mean function
    sig2list <- NULL
    Klist <- NULL
    flist <- NULL
    for(iter in 1:niter){
      b    <- sample.b(K,sig2)   # Gibbs transition for beta
      sig2 <- sample.sig2(K,b)   # Gibbs for sig2
      th <-   rj(K,b,sig2)       # RJ for changing K
      K <- th$K
      b <- th$b
      idx <- 1:(2*K+1)
      ## update lists
      sig2list <- c(sig2list,sig2)
      Klist <-    c(Klist,K)
      f <- X0[,idx]%*%b
      flist <- cbind(flist,f)
    }
    return(list(K=Klist,f=flist,sig2=sig2list))
  }
result = mcmc(10000)
Klist = result$K
f = result$f
sig2list = result$sig2

Ef = apply(f,1,mean)  
plt.dta(plt.spline=T)
points(seq(from=0,to=1.1,length=n0),Ef,type = "l",col=4,pch=16)

plot(Klist, type="l")
hist(Klist, freq = F,breaks=3:12, xlab="K",ylab="p(K/y)",main = paste("Histogram of K"))

plot(sig2list, type="l")
hist(sig2list,freq = F,xlab="sig2",ylab="p(sig2/y)",main = paste("Histogram of sig2"))

### problem 6
Sk.det = rep(0,10)
betak.bar = list()
p.ybetabar = rep(0,10)
p.betabar = rep(0,10)
p.yk = rep(0,10)
for(i in 1:10){
  Sk.det[i] = 1/det(0.1 * diag(2*i+1) + 2 * t(X[,1:(2*i+1)]) %*% X[,1:(2*i+1)])
  betakbar = 2 * solve(0.1 * diag(2*i+1) + 2 * t(X[,1:(2*i+1)]) %*% X[,1:(2*i+1)]) %*% t(X[,1:(2*i+1)]) %*% y
  betak.bar[[i]] = list(beta = betakbar)
  p.ybetabar[i] = prod(dnorm(y, mean = X[,1:(2*i+1)] %*% betak.bar[[i]]$beta, sd = sqrt(0.5)))
  p.betabar[i] = prod(dnorm(betak.bar[[i]]$beta, mean = 0, sd = sqrt(10)))
  p.yk[i] = sqrt(Sk.det[i]) * p.ybetabar[i] * p.betabar[i]
}

barplot(p.yk/sum(p.yk),names.arg = 1:10,xlab="K",ylab="p(y/K)")


