##Group member:
##1. Aditya Prabaswara Mardjikoen(s2264710)
##2. Xiangtian Duan(s2248742)
##3. Huiying Chen(s2264943)

## Github repositories link: 
## https://github.com/aprabaswara/BFGS-Optimizer.git

##overview:
##---------------------------------------------------------------------------------------------------------------------------------------------
##Finite differencing approximation reduces reliability, accuracy, and efficiency a little,
##but for many cases it is good enough.

##Sys.setenv("_R_USE_PIPEBIND_" = "true")
##test function that professor provided to test the BFGS method
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by 'bfgs'
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

bfgs <- function(theta, f, ..., tol, fscale, maxit){
  ##function to implementing the BFGS Quasi-Newton minimization method.
  
  ##input:
  ##theta = vector of initial values for the optimization parameter; f = objective function to minimize;
  ##... = any argument of f after the parameter vector and gradient logical indicator; tol = the convergence tolerance;
  ##fscale = a rough estimate of the magnitude of f at the optimum (used in covergence testing); maxit = the maximum number of BFGS iterations.
  
  ##output:
  ##list containing this following element:
  ##f = scalar value of the objective function at the minimum; theta = vector of values of the parameters at the minimum;
  ##iter = number of iterations taken to reach the minimum; g = gradient vector at the minimum;
  ##H = approximate Hessian matrix (obtained by finite difference method) at the minimum.
  
  g <- function(theta,f,...){
    ##function to calculate gradient given the value of theta and function f
    
    ##input: 
    ##theta = vector of initial values for the optimization parameter; f = objective function to minimize;
    ##... = any argument of f after the parameter vector and gradient logical indicator;
    
    ##output: grad = gradient of theta
    
    #compute gradient by finite difference approximation if the user didn't supply it
    
    if (is.null(attr(f(theta,...),"gradient"))){
      #if (f[2] == TRUE){
      f0 <- f(theta,...) ## f at theta
      eps <- 1e-7 ## finite difference interval
      grad <- c() ##initialize vector to store gradient values
      
      for (i in 1:length(theta)) { ## loop over parameters
        
        th <- theta; th[i] <- th[i] + eps ## increase theta by eps
        f1 <- f(th,...) ## compute f at theta+eps
        grad[i] <- (f1 - f0) / eps ## approximate first order derivative
      }
    }
    
    #get the gradient from f if the user supply it
    else{
      grad <- attr(f(theta,...),'gradient')
    }
    
    return(grad)
  }
  #
  param_num <- length(theta) ##length of vector theta
  B <- I <- diag(param_num) ##initialize identity matrix as the initial value of the inverse of Hessian matrix
  iter<-0
  f0<-f(theta,...)
  f1<-f0+abs(f0)/2
  grad <- g(theta,f,...)
  delta <- - B %*% g(theta,f,...)
  for (i in 1:maxit){
    iter <- iter +1
    #initial value
    grad_index <- which(is.finite(g(theta,f,...))==FALSE)
    ##check if objective function is finite or not
    
    if (is.finite(f(theta,...))==FALSE){
      warning("Objective function is not finite at given initial theta value!")
      break
    }
    
    ##check if derivative is finite or not
    
    #else if (g(theta,f,...) %in% c(-Inf,Inf) | is.na(g(theta,f,...)))
    else if (length(grad_index)!=0){
      warning("Derivative is not finite at given initial theta value!")
      break
    }
    
    if (max(abs(grad)) < (abs(f(theta,...)) + fscale) * tol){
      
      ##condition for convergence
      break
    }
    n <- 0
    while (!is.finite(f1)|f1>f0){
      
      n <- n +1
      #if(n > maxit){
      # warning("Steps failed to reduce objective before convergence occured!")
      #break
      
      #}
      delta <- delta / 2#make f1 get smaller
      
      f0 <- f(theta,...)
      f1 <- f(theta + delta,...)
    }
    #now reduce step we get decrese on objective
    #condition 2
    #delta <- - B %*% g(theta,f,...)
    grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
    grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
    c2 <- 0.9
    m <- 0.1
    
    n <- 0
    while(grad_delta2 < c2 * grad_delta1){
      m <- m / 2
      n<- n+1
      if (n > maxit){
        warning("Steps failed to reduce objective before convergence occured!")
        break
      }else if(!is.finite(f1) | f1 > f0){
        delta <-  delta / 2
        
      }else{
        delta<-(1 + m) * delta
        
      }
      
      grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
      grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
      
      f0 <- f(theta,...) 
      f1 <- f(theta + delta,...)
      
    }
    
    s <- delta
    theta_new <- theta + s
    grad <- g(theta_new,f,...)
    # get delta satisfy condition 2
    #iteration of B
    y <- g(theta_new,f,...) - g(theta,f,...)
    inv_rho <- drop(crossprod(delta,y)) ##rho^(-1)=s^T*y
    rho <- 1/inv_rho ##1/rho^(-1)
    
    rho_syt <- rho * tcrossprod(s,y) ##rho*s*y^T
    rho_yst <- rho * tcrossprod(y,s) ##rho*y*s^T
    matrix1 <- I - rho_syt ##I-rho*s*y^T
    matrix2 <- I-rho_yst ##I-rho*y*s^T
    B <- (matrix1 %*% B) %*% matrix2 + rho * tcrossprod(s) ##update B
    
    theta <- theta_new
    grad <- g(theta,f,...)
    delta <- drop(- B %*% g(theta,f,...))
    
  }
  
  #if (iter==maxit & max(abs(grad)) >= (abs(f(theta,...)) + fscale) * tol){
  #warning('Convergence not reached until maximum iteration!')
  #}
  
  
  H <- matrix(0, nrow = param_num, ncol = param_num) ##initialize Hessian matrix
  
  g0 <- g(theta,f,...) ##gradient value at the minimum theta value
  eps <- 1e-7 ## finite difference interval
  
  for (i in 1:param_num) { ## loop over parameters
    th <- theta; th[i] <- th[i] + eps ##increase theta by eps
    g1 <- g(th,f,...) ##gradient value at theta+eps
    
    H[i,] <- (g1 - g0)/eps ## approximate second derivatives in Hessian matrix
    
  }
  
  ##convert asymmetric Hessian matrix into symmetric Hessian matrix
  if (isSymmetric(H)==FALSE){
    H <- 0.5 * (t(H) + H)
  }
  
  f_optimum <- f(theta,...) ##scalar value of objective function at the minimum
  
  list(f = f_optimum, theta = theta, iter = iter, g = g0, H = H)
}
bfgs(theta=c(-1,2),f=rb,getg=FALSE,k=10,tol=1e-5,fscale=1,maxit=100)