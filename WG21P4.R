##Group member:
##1. Aditya Prabaswara Mardjikoen(s2264710)
##2. Xiangtian Duan(s2248742)
##3. Huiying Chen(s2264943)

## Github repositories link: 
## https://github.com/aprabaswara/BFGS-Optimizer.git

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
  #function to 
  #input:
  #theta = vector of initial values for the optimization parameter; f = objective function to minimize;
  #... = any argument of f after the parameter vector and gradient logical vector; tol = the convergence tolerance
  #fscale = a rough estimate of the magnitude of f at the optimum (used in covergence testing); maxit = the maximum number of BFGS iterations to try before giving up
  n <- length(theta) #length of vector theta
  s <- c()
  s[k] <- theta[k+1]-theta[k] 
  if (is.null(attr(f,"gradient"))){
    
  }
  list(f, theta, iter, g, H)
}