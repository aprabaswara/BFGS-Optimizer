##Group member:
##1. Aditya Prabaswara Mardjikoen(s2264710)
##2. Xiangtian Duan(s2248742)
##3. Huiying Chen(s2264943)

## Github repositories link: 
## https://github.com/aprabaswara/BFGS-Optimizer.git

##overview:
##---------------------------------------------------------------------------------------------------------------------------------------------
##Finite differencing approximation reduces reliability, accuracy, and efficiency a little,
##but for many cases it is good enough. In this code, we will implement the BFGS method by
##minimizing the initial value using this algorithm. One advantage of using this method is 
##we mostly dealing with first order derivative instead of second order derivative. Here 
##is a short overview how the algorithm works:
##1. Initialize initial value, gradient, and inverse of Hessian matrix, which usually set as 
##   identity matrix or any matrix that approximate the inverse of Hessian.
##
##2. Check whether the gradient approximate zero, if yes terminate the process and if not continue
##   the process until the gradient approximate zero.
##
##3. compute delta by multiplying the negative value of inverse Hessian with its gradient at the 
##   current initial value.
##
##4. Check whether the delta reduces the objective. If not, decrease the delta value until it reduces the objective. 
##
##5. Check whether the delta satisfies the second wolfe condition. If not, increase delta slightly until it satisfies this 
##   condition. 
##
##6. Set the upper limit of iteration times of delta, otherwise we will stuck in this process for too long. 
##   If we still can not obtain appropriate delta when iteration times exceed the limit, we move to the next step.
##
##7. Update the theta using the new delta by adding it to the previous theta.
##
##8. Compute the new matrix B.
##
##9. Repeat the whole process until the gradient approaches zero.
##
##10. For the whole process, we issue error if objective or derivative is not finite at some 
##    theta, objective function values not decreased before convergence, and maximum iteration
##    is reached before convergence.
##
##11. If no error in point 10, compute the Hessian matrix for optimal values using finite difference on the
##    gradient.
##
##12. Compute the optimal values for objective functions, gradient, and the number of iteration 
##    to reach convergence if no error in point 10.
##----------------------------------------------------------------------------------------------------------

g <- function(theta,f,...){
  
  ## function: calculate gradient given the value of theta and function f
  ##
  ## input:     theta: vector of initial values for the optimization parameter
  ##           f    : objective function
  ##           ...  : any argument of f after the parameter vector and gradient logical indicator;
  ## output:    grad:  gradient of theta
  
  if (is.null(attr(f(theta,...),"gradient"))){## compute gradient by finite difference approximation if the user didn't supply it
    f0 <- f(theta,...) ## f at theta
    eps <- 1e-7 ## finite difference interval
    grad <- c() ##initialize vector to store gradient values
    
    for (i in 1:length(theta)) { ## loop over parameters
      
      th <- theta; th[i] <- th[i] + eps ## increase theta by eps
      f1 <- f(th,...) ## compute f at theta+eps
      grad[i] <- (f1 - f0) / eps ## approximate first order derivative
    }
  }
  else{ ## if the user supplies gradient
    grad <- attr(f(theta,...),'gradient') ## get the gradient from f 
  }
  return(grad)
}

bfgs <- function(theta, f, ..., tol = 1e-5, fscale = 1, maxit = 100){
  
  ## function: implementing the BFGS Quasi-Newton minimization method.
  
  ## input:
  ## theta  = vector of initial values for the optimization parameter
  ## f      = objective function to minimize
  ##          Its first argument is the vector of optimization parameters
  ##          Its second argument is a logical indicating whether or not gradients of the objective 
  ## ...    = any argument of f after the parameter vector and gradient logical indicator
  ## tol    = the convergence tolerance;
  ## fscale = a rough estimate of the magnitude of f at the optimum (used in convergence testing)
  ## maxit  = the maximum number of BFGS iterations to try before giving up
  ##         
  ## output:(named list containing following element)
  ##
  ## f     = scalar value of the objective function at the minimum
  ## theta = vector of values of the parameters at the minimum
  ## iter  = number of iterations taken to reach the minimum
  ## g     = gradient vector at the minimum
  ## H     = approximate Hessian matrix (obtained by finite difference method) at the minimum.
  
  
  c2 <- 0.9                  ## constant for checking second Wolfe condition
  param_num <- length(theta) ## length of vector theta
  B <- I <- diag(param_num)  ## initialize identity matrix as the initial value of the inverse of Hessian matrix
  iter <- 0                  ## initialize iteration
  f0 <- f(theta,...)         ## initialize objective function values
  grad <- g(theta,f,...)     ## initialize gradient
  
  
  ## Before starting the whole process we have to check weather the initial objective function
  ## values or gradient is finite. If not, we could not continue the BFGS algorithm.
  
  grad_index <- which(is.finite(grad)==FALSE)## find gradient element that isn't finite
  if (is.finite(f0)==FALSE){                 ## issue error if objective function values is not finite
    stop("Objective function is not finite at given initial theta value!")
  }
  else if (length(grad_index) != 0){           ## issue error if derivative is not finite
    stop("Derivative is not finite at given initial theta value!")
  }
  
  while(iter <= maxit & max(abs(grad)) >= (abs(f(theta,...)) + fscale) * tol){## when the iteration is less than 
    ## the limit and it is not convergent
    ## we implement BFGS algorithm
    iter <- iter + 1          ## record the number of iterations
    f0 <- f(theta,...)          ## objective function value with initial value of parameters
    grad <- g(theta,f,...)    ## gradient of objective function with initial value of parameters
    
    
    ## In this part, we will check is there any objective function values or 
    ## gradient that is not finite during the BFGS algorithm iteration. If there
    ## exist objective function values or gradient that is not finite, we could 
    ## not continue the BFGS algorithm during this iteration
    
    grad_index <- which(is.finite(grad) == FALSE) ## find gradient element that isn't finite
    if (is.finite(f0) == FALSE){                  ## issue error if objective function values is not finite
      stop("Objective function is not finite at given initial theta value!")
      break
    }
    else if (length(grad_index) != 0){            ## issue error if derivative is not finite
      stop("Derivative is not finite at given initial theta value!")
      break
    }
    
    ## initial value for checking second wolfe condition 
    delta <- drop(- B %*% g(theta,f,...))           
    f1 <- f(theta + delta,...)
    grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
    grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
    
    count <- 0                                                          ## record the number of loops in order to limit the number of cycles
    
    while(!is.finite(f1) | f1 > f0 | grad_delta2 < c2 * grad_delta1) {  ## cases where delta adjustment is required
      count<- count + 1                                                 ## update trial
      
      if(count > maxit){                                               ## terminate loop if its times exceed the upper limit
        break
      }
      
      else if(!is.finite(f1) | f1 > f0){                                ## if the objective is infinite and if delta
        ## can not decrease the objective function value
        delta <- delta / 2                                              ## reduce delta 
      }
      
      
      else if(grad_delta2 < c2 * grad_delta1){                          ## if delta makes objective function decrease and 
        ## it do not satisfy the second Wolfe condition
        delta <- (1 + 0.1) * delta                                    ## increase delta slightly
      }
      
      
      else{                                                             ## terminate process if delta make objective function 
        ## decrease and it satisfies the second Wolfe condition
        break
      }
      
      f1 <- f(theta + delta,...)                                        ## update f1 using new delta to prepare for the loop
      grad_delta1 <- drop(crossprod(g(theta,f,...),delta))              ## the right side of the Wolfe condition 2
      grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))      ## the left side of the Wolfe condition 2
    }
    
    s <- delta                                                          ## s equals to the final delta in the iteration                                
    theta_new <- theta + s                                              ## update the new theta with s
    
    
    ## issue error if steps failed to reduce objective function values before convergence
    if (f(theta_new,...) > f(theta,...)){
      stop('Steps failed to reduced objective values but convergence not occured!')
      break
    }
    
    ## compute the B matrix (approximate inverse Hessian)
    y <- g(theta_new,f,...) - g(theta,f,...)                                                         ## difference in gradients at old and new parameters
    inv_rho <- drop(crossprod(delta,y))                                                              ## rho^(-1) = s^T*y
    rho <- 1 / inv_rho                                                                               ## 1 / rho^(-1)
    s_yt <- tcrossprod(s,y)                                                                          ## s*y^T
    y_st <- tcrossprod(y,s)                                                                          ## y*s^T
    
    B <- B - rho * B %*% y_st - rho * s_yt %*% B + rho^2 * s_yt %*% B %*% y_st + rho * tcrossprod(s) ## update B, unpack the formula to calculate
    ## in order to have the cost O(p^2) rather than
    ## O(p^3)
    theta <- theta_new                                                                               ## update parameter value
    grad <- g(theta,f,...)                                                                           ## the current gradient
    delta <- drop(- B %*% g(theta,f,...))                                                            ## delta computed with new B matrix
  }
  
  ##issue error if convergence not reached until maximum iteration
  if (max(abs(grad)) >= (abs(f(theta,...)) + fscale) * tol){
    stop('Convergence not reached until maximum iteration!')
  }
  
  H <- matrix(0, nrow = param_num, ncol = param_num)               ## initialize Hessian matrix
  g0 <- g(theta,f,...)                                             ## gradient value at the minimum theta value
  eps <- 1e-7                                                      ## finite difference interval
  
  for (i in 1:param_num) {                                         ## loop over parameters
    th <- theta; th[i] <- th[i] + eps                              ## increase theta by eps
    g1 <- g(th,f,...)                                              ## gradient value at theta + eps
    H[i,] <- (g1 - g0) / eps                                       ## approximate second derivatives in Hessian matrix
  }
  
  ##convert asymmetric Hessian matrix into symmetric Hessian matrix
  if (isSymmetric(H) == FALSE){
    H <- 0.5 * (t(H) + H)
  }
  
  f_optimum <- f(theta,...)                                                      ## scalar value of objective function at the minimum
  bfgs_summary <- list(f = f_optimum, theta = theta, iter = iter, g = g0, H = H) ## create the named list
  return(bfgs_summary)                                                           ## return the list
}