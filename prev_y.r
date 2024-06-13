source("rbeta_function.r")

prev_y<-function(fit,ar=NA,ma=NA,h=6)
{
  
  mlink<-make.link(link)
  linkinv<-mlink$linkinv
  linkfun <- mlink$linkfun
  n <- length(fit$serie)
  ynew <- linkfun(fit$serie)
  ynew_prev <- c(ynew,rep(NA,h))
  errorhat<-rep(0,n) 
  mu_hat<- c(fit$fitted,rep(0,h))
  eta_hat<- linkfun(mu_hat)
  y_prev <- c(fit$serie,rep(0,h))
  eta_hat <- c(eta_hat,rep(0,h))
  
  ##### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F))
  { 
    p1 <- length(ar)
    q1 <- length(ma)
    ar<-1:(length(ar))
    ma<-1:(length(ma))
    
    alpha <- fit$coeff[1]
    phi <- fit$coeff[2:(p1+1)] 
    theta <- fit$coeff[(p1+2):(p1+q1+1)]
    prec <- fit$coeff[p1+q1+2] # precision parameter
    
    
    for(i in 1:h)
    {
      eta_hat[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      mu_hat[n+i]<- linkinv(eta_hat[n+i])
      y_prev[n+i] <- rbeta_function(1, mu_hat[n+i],prec)
      ynew_prev[n+i]<- linkfun(y_prev[n+i])
      errorhat[n+i] <- ynew_prev[n+i]-eta_hat[n+i] # original scale 
    }

    z<-c()
    z$y_prev<-y_prev[(n+1):(n+h)]
    z$mu_hat<-mu_hat[(n+1):(n+h)]
  }
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T))
  {
    p1 <- length(ar)
    ar<-1:(length(ar))
    
    alpha <- fit$coeff[1]
    phi <- fit$coeff[2:(p1+1)] 
    prec <- fit$coeff[p1+2] # precision parameter

    for(i in 1:h)
    {
      eta_hat[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) 
      mu_hat[n+i]<- linkinv(eta_hat[n+i])
      y_prev[n+i] <- rbeta_function(1, mu_hat[n+i],prec)
      ynew_prev[n+i]<- linkfun(y_prev[n+i])
      errorhat[n+i] <- ynew_prev[n+i]-eta_hat[n+i] # original scale
    }

    z<-c()
    z$y_prev<-y_prev[(n+1):(n+h)]
    z$mu_hat<-mu_hat[(n+1):(n+h)]
  
  }

  ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F))
  {
    q1 <- length(ma)
    ma<-1:(length(ma))
    
    alpha <- fit$coeff[1]
    theta <- fit$coeff[2:(q1+1)]
    prec <- fit$coeff[q1+2] # precision parameter
    
    for(i in 1:h)
    {
      eta_hat[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      mu_hat[n+i]<- linkinv(eta_hat[n+i])
      y_prev[n+i] <- rbeta_function(1, mu_hat[n+i],prec)
      ynew_prev[n+i]<- linkfun(y_prev[n+i])
      errorhat[n+i] <- ynew_prev[n+i]-eta_hat[n+i] # original scale  
    }
    
    z<-c()
    z$y_prev<-y_prev[(n+1):(n+h)]
    z$mu_hat<-mu_hat[(n+1):(n+h)]
    
  }

return(z)
}

