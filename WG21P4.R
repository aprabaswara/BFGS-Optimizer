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
      grad <- attr(f(theta,...),"gradient")
       
    }
    
    return(grad)
  }
 
  
  #step_reduce <- function(theta, delta, g, f, ...){
    ##function to reduce the step until satisfying the second Wolfe condition
    
    ##input: 
    ##theta = vector of initial values for the optimization parameter; f = objective function to minimize;
    ##... = any argument of f after the parameter vector and gradient logical indicator; g = gradient function
    ##delta = the perturbation vector that we want to reduce
    
    ##output: delta = the perturbation vector that satisfying the second Wolfe condition and reduce the values of objective function
    
    #c2 = 0.9 ##constant for checking the second Wolfe condition
    #grad_delta1 <- crossprod(g(theta,f,...),delta) ##grad(theta)^T delta
    #grad_delta1 <- drop(grad_delta1)
    #grad_delta2 <- crossprod(g(theta+delta,f,...),delta) ##grad(theta+delta)^T delta
    #grad_delta2 <- drop(grad_delta2)
    
    ##find delta that satisfying second Wolfe condition and decreases the objective function
    #while (grad_delta2 < c2 * grad_delta1){
      #delta <- delta/2 ##halve the delta
      #grad_delta1 <- crossprod(g(theta,f,...),delta) ##update grad(theta)^T delta
      #grad_delta1 <- drop(grad_delta1)
      #grad_delta2 <- crossprod(g(theta+delta,f,...),delta) ##update grad(theta+delta)^T delta
      #grad_delta2 <- drop(grad_delta2)
      #return(delta)
    #}
    

    step_length <- function(theta, f, delta, B, g,...){
      
      step_direc = - B %*% g(theta,f,...)#grad?
      grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
      grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
      c2 <- 0.9
      #ieration according to the second condition of Wolfe conditions                    
      n <- 1
      while (grad_delta2 < c2 * grad_delta1){
        delta <- (1 + 1 / (n + 1)) * delta 
        grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
        grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
        n <- n + 1
        
      }
      return(delta / step_direc)
    }
      #iteration for delta
    
      #method 1(new material)
      #initial value
      #step_direc = - B %*% g(theta,f,...)
      
      #delta <- step_direc * step_len
      #c <- 3 / 2
      #while(f(theta + delta_initial,...) > (f(theta,...) + c1 * crossprod(g(theta), delta))){
        #delta <- c * delta_initial  
      #}
      
      # method 2
      B <- I 
      step_direc = - B %*% g(theta,f,...)
      delta <- step_direc * 2
      iter <- 0
      step_len <- step_length(theta, f, delta, B, g,...)
      s <- step_len * step_direc
      theta_new <- theta + s
      grad <- g(theta_new,f,...)

      while (iter <= maxit | max(abs(grad)) <= (abs(f(theta_new)) + fscale) * tol){
        iter <- iter + 1
        f0 <- f(theta,...)
        step_len <- step_length(theta, f, delta, B, g,...)
        s <- step_len * step_direc
        theta_new <- theta + s
        f1 <- f(theta_new,...)
        
        if (f1 < f0 | f1 %in% c(-Inf,Inf) | is.na(f1)){#f0 na?
          delta <- delta / 2
          step_len <- step_length(theta, f, delta, B, g,...)
          
        }else{#havn't finish,out put error
          
        }
        y <- g(theta_new,f,...) - g(theta,f,...)
        inv_rho <- drop(crossprod(delta,y)) ##rho^(-1)=s^T*y
        rho <- 1/inv_rho ##1/rho^(-1)
        rho_syt <- rho * tcrossprod(s,y) ##rho*s*y^T
        rho_yst <- rho * tcrossprod(y,s) ##rho*y*s^T
        matrix1 <- I - rho_syt ##I-rho*s*y^T
        matrix2 <- I-rho_yst ##I-rho*y*s^T
        B <- (matrix1 %*% B) %*% matrix2 + rho * tcrossprod(s) ##update B
        theta_new <- theta
        step_direc = - B %*% g(theta,f,...)
        step_len <- step_length(theta, f, delta, B, g,...)
        s <- step_len * step_direc
        theta_new <- theta + s
        grad <- g(theta_new,f,...)
        
      }
       
      
      
      
      
  ##update start from here (find update that is numbered)
  n <- length(theta) ##length of vector theta
  I <- diag(n) ##initialize identity matrix
  B <- I ##initial value of the inverse of Hessian matrix
  
  iter <- 0 ##initialize first iteration
  
  while(iter <= maxit){
    iter <- iter+1 ##increase iteration by one
    
    ##update one
    grad_index <- which(g(theta,f,...) %in% c(-Inf,Inf) | is.na(g(theta,f,...)))
    ##check if objective function is finite or not
    
    if (f(theta,...) %in% c(-Inf,Inf) | is.na(f(theta,...))){
      warning("Objective function is not finite at given initial theta value!")
      break
    }
    
    ##check if derivative is finite or not
    
    #else if (g(theta,f,...) %in% c(-Inf,Inf) | is.na(g(theta,f,...)))
    else if (length(grad_index)!=0){
      warning("Derivative is not finite at given initial theta value!")
      break
    }
    ##-----------------------------------------------------------------
    
    ##take the absolute value of each element in f and its gradient
    abs_grad <- abs(g(theta,f,...)) 
    abs_f <- abs(f(theta,...)) 
    
    ##update two
    ##exit the loop if the convergence is reached
    if (max(abs_grad) < (abs_f+fscale)*tol){
      break
    }
    
    ##check if maximum iteration is reached without convergence
    else if (iter == maxit){
      warning("Maximum iteration is reached without convergence!")
      break
    } 
    ##-------------------------------------------------------------
    
    delta <- -B %*% g(theta,f,...) ##compute step length
    delta <- drop(delta) ##return delta as a vector
    delta <- step_reduce(theta, delta, g, f, ...) ##find delta that reduce objective function
    
    ##update 3
    ##check if the step reduces the value of objective function
    if (f(theta+delta,...) >= f(theta,...)){
      warning("Steps failed to reduce the objective before convergence occured!")
      break
    } ##still unsure?
    ##---------------------------------------------------------------------------
    
    theta_new <- theta+delta ##update the theta value
    s <- theta_new-theta
    y <- g(theta_new,f,...)-g(theta,f,...)
    inv_rho <- crossprod(s,y) ##rho^(-1)=s^T*y
    inv_rho <- drop(inv_rho)
    rho <- 1/inv_rho ##1/rho^(-1)
    
    rho_syt <- rho*tcrossprod(s,y) ##rho*s*y^T
    rho_yst <- rho*tcrossprod(y,s) ##rho*y*s^T
    
    ##update 5
    matrix1 <- I-rho_syt ##I-rho*s*y^T
    matrix2 <- I-rho_yst ##I-rho*y*s^T
    matrix3 <- (matrix1 %*% B) %*% matrix2 ##(I-rho*s*y^T)B(I-rho*y*s^T)
    
    B <- matrix3 + rho*tcrossprod(s) ##update B
    ##------------------------------------------------------------------
    theta <- theta_new ##update theta for the next iteration
  }
  
  
  H <- matrix(0, nrow = n, ncol = n) ##initialize Hessian matrix
  g0 <- g(theta,f,...) ##gradient value at the minimum theta value
  eps <- 1e-7 ## finite difference interval
  
  for (i in 1:n) { ## loop over parameters
    th <- theta; th[i] <- th[i] + eps ##increase theta by eps
    g1 <- g(th,f,...) ##gradient value at theta+eps
    H[i,] <- (g1 - g0)/eps ## approximate second derivatives in Hessian matrix
  }
  
  ##convert asymmetric Hessian matrix into symmetric Hessian matrix
  if (isSymmetric(H)==FALSE){
    H <- 0.5 * (t(H) + H)
  }
  
  f_optimum <- f(theta,...) ##scalar value of objective function at the minimum
  
  ##update 4
  optim_summary <- list(f = f_optimum, theta = theta, iter = iter, g = g0, H = H)
  return(optim_summary)
  ##-----------------------------------------------------------------------------
}