##Group member:
##1. Aditya Prabaswara Mardjikoen(s2264710)
##2. Xiangtian Duan(s2248742)
##3. Huiying Chen(s2264943)

## Github repositories link: 
## https://github.com/aprabaswara/BFGS-Optimizer.git

##overview:
##---------------------------------------------------------------------------------------------------------------------------------------------
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
    
    #compute gradient by forward difference if the user didn't supply it
    if (is.null(attr(f(theta,...),"gradient"))){
      f0 <- f(theta,...) ## f at theta
      eps <- 1e-7 ## finite difference interval
      grad <- c() ##initialize vector to store gradient values
      
      for (i in 1:length(theta)) { ## loop over parameters
        
        th <- theta; th[i] <- th[i] + eps ## increase theta by eps
        f1 <- f(th,...) ## compute f at theta+eps
        grad[i] <- (f1 - f0)/eps ## approximate first order derivative
      }
    }
    
    #get the gradient from f if the user supply it
    else{
      grad <- attr(f(theta,...),"gradient")
    }
    
    return(grad)
  }
  
  #check if derivative is finite or not
  if (g(theta,f,...) %in% c(-Inf,Inf) | is.na(g(theta,f,...))){
    warning("Derivative is not finite at given initial theta value!")
  }
  
  c1 = 1e-4 #constant for checking the first Wolfe condition
  c2 = 0.9 #constant for checking the second Wolfe condition
  n <- length(theta) #length of vector theta
  B <- diag(n) #initial value of the inverse of Hessian matrix
  
  iter = 1 #initialize first iteration
  
  while(iter <= maxit){
    
    #exit the loop if the maximum of iteration is one
    if (maxit == 1){
      break
    }
    
    #take the absolute value of each element in f and its gradient
    abs_grad <- abs(g(theta,f,...)) 
    abs_f <- abs(f(theta,...)) 
    
    delta <- -B %*% g(theta,f,...) ##compute step length
    delta <- drop(delta) ##return delta as a vector
    iter = iter+1 ##increase iteration by one
  }
  
  #check if maximum iteration is reached without convergence
  if (iter == maxit & max(abs_grad) => (abs_f+fscale)*tol){
    warning("Maximum iteration is reached without convergence!")
  } 
  
  H <- matrix(0, nrow = n, ncol = n) ##initialize Hessian matrix
  g0 <- g(theta,f,...) ##gradient value at the minimum theta value
  eps <- 1e-7 ## finite difference interval
  
  for (i in 1:length(theta)) { ## loop over parameters
    th <- theta; th[i] <- th[i] + eps ##increase theta by eps
    g1 <- g(th,f,...) ##gradient value at theta+eps
    H[i,] <- (g1 - g0)/eps ## approximate second derivatives in Hessian matrix
  }
  
  #convert asymmetric Hessian matrix into symmetric Hessian matrix
  if (max(abs(H-t(H))) != 0){
    H <- 0.5 * (t(H) + H)
  }
  
  f_optimum <- f(theta,...) ##scalar value of objective function at the minimum
  
  list(f = f_optimum, theta = theta, iter = iter, g = g0, H = H)
}