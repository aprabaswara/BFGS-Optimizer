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








step_reduce <- function(theta, f, B, g,...){
  direction_delta = - B %*% g(theta,f,...)#grad?
  print(2)
  alpha=1e-4
  n=1
  c1=1e-4
  c2=0.9
  while(TRUE){
    alpha=alpha*1.1
    print(1)
    if(f(theta+alpha*direction_delta)>f(theta)+c1*crossprod(g(theta,f,...),alpha*direction_delta)){
      return(alpha)
      break
    }
  }
  while(TRUE){
    alpha=(1-(1/2)^n)*alpha
    n=n+1
    print(n)
    if((crossprod(g(theta+alpha*direction_delta,f,...),direction_delta))<c2*crossprod(g(theta,f,...),alpha*direction_delta)){
      delta=alpha*direction_delta
      return(delta)
      break
    }
  }
}