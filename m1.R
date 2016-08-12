
require(gtools)  # I found a Dir r.v. generator in this package
## if needed install it with
## install.packages("gtools")  

## DATA
x <- scan("/Users/zk794/Box Sync/DISC/Austin Courses/2016 Spring/Monte Carlo Methods in Stats/mid-term/sage.dta")    # raw data
y <- table(x)                        # counts
N <- length(y)       
names(y) <- 1:N
n <- length(x)

## HYPERPARS:
rho <- 0.1
as <- 0.9
bs <- 0.1
a1 <- .1
a0 <- 10

## initialize
## this function creates a list with
##   z=(z1,.. zN); pis=pi*, r=(r~[1],..r~[M0])
##   q=(q~[1],..q~[M1])
## Youc an use it to initialize the state of the MC
init <- function()
   { # initialize parameters
     ## z
     z <- ifelse(y<10,0,1)
     ## pi*: empirical frquency
     A0 <- which(z==0);     A1 <- which(z==1)
     M0 <- length(A0);      M1 <- length(A1)
     Y0 <- sum(y[A0]);      Y1 <- sum(y[A1])
     pis <- sum(Y1)/n
    ## r and q: empirical fequencies
    q <- y[A1]/Y1  # this is the q~ of the text
    r <- y[A0]/Y0  # this is the r~ of the text
    return(th=list(z=z,pis=pis,r=r,q=q))
  }

## function for sampling z 
# 1. z ~ p(z | pis, y)
sample.z <- function(pis,z)
  {
    znew <- z # initialize - will save z here
    
    for(j in 1:N){ # loop over all zs
    	A0 <- which(znew==0);     A1 <- which(znew==1)
    	M0 <- length(A0);      M1 <- length(A1)
      M0_j <- length(A0[-j]);      M1_j <- length(A1[-j])
     	Y0 <- sum(y[A0]);      Y1 <- sum(y[A1])
      Y0_j <- sum(y[A0[-j]]);      Y1_j <- sum(y[A1[-j]])

      logpz1 = log(pis^y[j]*rho)+lgamma((1+M1_j)*a1)+lgamma(M0_j*a0)+lgamma(a1+y[j])-lgamma(a1)-lgamma((1+M1_j)*a1+Y1_j+y[j])-lgamma(M0_j*a0+Y0_j)
      logpz0 = log((1-pis)^y[j]*(1-rho))+lgamma(M1_j*a1)+lgamma((1+M0_j)*a0)+lgamma(a0+y[j])-lgamma(a0)-lgamma(M1_j*a1+Y1_j)-lgamma((1+M0_j)*a0+Y0_j+y[j])
    	pz1=1/(1+exp(logpz0-logpz1))
      znew[j] <- rbinom(1,1,pz1) 	# update
    }
    return(array(znew))
  }
  
## function for sampling q 
# 2. q ~ p(q | pis, z, y)
sample.q <- function(pis,z)
  {
    q = rep(0,N)
    A0 <- which(z==0);     A1 <- which(z==1)
    q[A0] = 0
    q[A1] = rdirichlet(1, a1+y[A1])
    return(array(q))
  }

## function for sampling r 
# 3. r ~ p(r | pis,z,y)
sample.r <- function(pis,z)
  {
    r = rep(0,N)
    A0 <- which(z==0);     A1 <- which(z==1)
    r[A1] = 0
    r[A0] = rdirichlet(1, a0+y[A0])
    return(array(r))
  }

## function for sampling pi 
# 4. pi
sample.pis <- function(z)
  {
    A0 <- which(z==0);     A1 <- which(z==1)
  	Y0 <- sum(y[A0]);      Y1 <- sum(y[A1])
    pis = rbeta(1, as+Y1, bs+Y0)
    return(pis)
  }
  
## main function for MCMC
gibbs <- function(n.iter=100, verbose=0)
  {
    TH <- NULL # initialize - will save pi*,z here
               ##             for each iteration
    PI <- NULL # similar - will save (pi1,.., piN) here
    th <-  init()  
    pis <- th$pis        # initialize pis = pi* 
    z <- th$z            # initialize z
    for(it in 1:n.iter){ # loop over iterations
      z <- sample.z(pis, z)   # 1. z ~ p(z | pis, y)
      q <- sample.q(pis,z)    # 2. q ~ p(q | pis,z,y)
      r <- sample.r(pis,z)    # 3. r ~ p(r | pis,z,y)
      pis <- sample.pis(z)    # 4. pi
      if (verbose > 0){
        if (it %% 10 ==0)       # print short summary
          prt.summary(z,q,r,pis)
      }
      ## save iteration
      TH <- rbind(TH, c(pis,z))

      pi0 <- pis*q+(1-pis)*r
      PI <- rbind(PI, pi0)
    }
    return(list(TH=TH, PI=PI))
  }

## run the MCMC :-)
ex <- function()
  { ## RUN these lines to get the plots
    n.iter <- 500
    gbs <- gibbs(n.iter)
    ## assume gbs returns a list with elements
    ## TH = (niter x p) matrix with each row being the
    ##      state (pi, z)
    ## PI = (niter x 1) vector with pi
    TH <- gbs$TH
    PI <- gbs$PI
    its <- 1:n.iter

    ## trajectory plot
    plot(its, TH[,1],xlab="ITER",ylab="PI*",bty="l",type="l")

    ## boxplot
    boxplot(log(PI))

    ## plotting posterior means vs. mle's
    pibar <- apply(PI,2,mean) # posterior means
    pihat <- as.numeric(y)/n
    plot(pihat, pibar, type="p", 
         pch=19, bty="l",xlab="MLE pihat", ylab="E(pi | y)")
    abline(0,1)

    ## same thing, zoom in to left lower corner
    plot(pihat, pibar, type="p", xlim=c(0,0.03), ylim=c(0,0.03),
         pch=19, bty="l",xlab="MLE pihat", ylab="E(pi | y)")
    abline(0,1)
  }

