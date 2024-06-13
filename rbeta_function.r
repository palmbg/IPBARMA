rbeta_function<-function(n,mu,phi)
{
  
  p<-mu*phi
  
  q<-(1-mu)*phi
  
  y_prev<-rbeta(n,p,q)
  
  return(y_prev)
  
  
}