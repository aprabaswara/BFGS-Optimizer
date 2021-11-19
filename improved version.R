#测试函数
c2<-0.9
rb <- function(theta,getg=FALSE,k=10) {
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
}

#梯度函数  
g <- function(theta,f,...){
  
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
  
  else{
    grad <- attr(f(theta,...),'gradient')
  }
  return(grad)
}

#计算step
delta_compute <- function(f0,f1,f,delta,theta,...){
  
  count<- 0
  grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
  grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))

  
  while(!is.finite(f1) | f1 > f0 | grad_delta2 < c2 * grad_delta1) {
    
    if(count > maxit){
      warning("spend too much time") 
      break
    }else if(!is.finite(f1) | f1 > f0){
      delta <- delta / 2
    
    }else if(grad_delta2 < c2 * grad_delta1){
      delta<- 1.5  * delta
   
    }else{
      break
    }
    f0 <- f(theta,...)
    f1 <- f(theta + delta,...)
    
    count <- count + 1
  }
 
  print(count)
}


#优化函数
bfgs <- function(theta, f, ..., tol, fscale, maxit){
  
 
  param_num <- length(theta) ##length of vector theta
  B <- I <- diag(param_num) ##initialize identity matrix as the initial value of the inverse of Hessian matrix
  iter<-0
  grad <- g(theta,f,...)

  while(iter <= maxit & max(abs(grad)) >= (abs(f(theta,...)) + fscale) * tol){
    

    iter <- iter +1
    count<-0
    f0<-f(theta,...)
    grad <- g(theta,f,...)
    
    grad_index <- which(is.finite(grad)==FALSE)
    ##check if objective function is finite or not
    
    if (is.finite(f0)==FALSE){
      warning("Objective function is not finite at given initial theta value!")
      break
    }
    
    ##check if derivative is finite or not
    
    #else if (g(theta,f,...) %in% c(-Inf,Inf) | is.na(g(theta,f,...)))
    else if (length(grad_index)!=0){
      warning("Derivative is not finite at given initial theta value!")
      break
    }
    
    delta <- drop(- B %*% g(theta,f,...))
    f1<-f(theta+delta,...)
    grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
    grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
    
    while( f1 > f0 | grad_delta2 < c2 * grad_delta1) {
      count<- count +1
      if(count > maxit){
        
        warning("spend too much time") 
        break
      }else if(!is.finite(f1) | f1 > f0){
        
        delta <- delta / 2}else if(grad_delta2 < c2 * grad_delta1){
          delta<-(1 + 0.2) * delta}else{
            break}
      f1 <- f(theta + delta,...)
      grad_delta1 <- drop(crossprod(g(theta,f,...),delta))
      grad_delta2 <- drop(crossprod(g(theta + delta,f,...),delta))
    }
    print(count)
    s <- delta
    theta_new <- theta + s
    
    if (f(theta_new,...) > f(theta,...)){
      warning('Steps failed to reduced objective values but convergence not occured!')
      break
    }
    
    #iteration of B
    y <- g(theta_new,f,...) - g(theta,f,...)
    inv_rho <- drop(crossprod(delta,y)) ##rho^(-1)=s^T*y
    rho <- 1/inv_rho ##1/rho^(-1)
    s_yt <- tcrossprod(s,y) ##s*y^T
    y_st <- tcrossprod(y,s) ##y*s^T
    B <- B-rho*B%*%y_st -rho*s_yt%*%B+rho^2*s_yt%*%B%*%y_st+ rho * tcrossprod(s) ##update B
    theta <- theta_new
    grad <- g(theta,f,...)
    delta <- drop(- B %*% g(theta,f,...))
    
  }
  
  if (max(abs(grad)) >= (abs(f(theta,...)) + fscale) * tol){
    warning('Convergence not reached until maximum iteration!')
  }
  
  
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