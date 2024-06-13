source("prev_y.r")
source("rbeta_function.r")

##########################################################
# intervalo de predição BCa
ip_boot<-function(fit,y=y,ym=ym,ar=NA,ma=NA,h=6,B=100,alpha_prec=0.1)
{
  
  mlink<-make.link(link)
  linkinv<-mlink$linkinv
  linkfun <- mlink$linkfun
  mu.eta <-  mlink$mu.eta
  M_resid_bca<-c()
  quantil <-1-alpha_prec/2
  z<-qnorm(quantil)
  
  
  ym <- ts(ym, start=(end(y)+c(0,1)), frequency=frequency(y))
  forecast<- ts(fit$forecast, start=(end(y)+c(0,1)), frequency=frequency(y))
  
  ##### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F))
  { 
    p1 <- length(ar)
    q1 <- length(ma)
    p<- max(ar)
    q<-max(ma)
    m <- max(p,q,na.rm=T)
    n <- length(fit$serie)-m
    n_first<-n+h
    phi1<-ar
    ma1<-ma
    h1 <- h
    
    #########################BCa#########################
    tc_bca<-0
    ystar <- (log(fit$serie/(1-fit$serie)))[(m+1):n]
    mustar <- (digamma(fit$fitted*fit$prec) - digamma((1-fit$fitted)*fit$prec))[(m+1):n]
    nu<- (trigamma(fit$fitted*fit$prec) + trigamma((1-fit$fitted)*fit$prec))[(m+1):n]
    resid1<- (ystar-mustar)/sqrt(nu)
    phi_more<- (fit$prec^3)*(psigamma((fit$forecast*fit$prec),deriv=2) - psigamma(((1-fit$forecast)*fit$prec),deriv=2))
    nu_more<- (fit$prec^2)*(trigamma(fit$forecast*fit$prec)+trigamma((1-fit$forecast)*fit$prec))
    a<- ((1/6) * (phi_more))/(nu_more^(3/2))
    
    #####################################################
    i<-0
    while(i<B)
    {
      
      n <- length(fit$serie)-m
      ### amostra BCa
      resid2<-sample(resid1,(n-m),replace=T)
      ystar_boot_bca<- (exp(mustar+ resid2*sqrt(nu))) / (1 + exp(mustar+ resid2*sqrt(nu)))
      y_boot_bca<-ts(ystar_boot_bca)
      
      if(min(y_boot_bca,na.rm=T)>0.0001 && max(y_boot_bca,na.rm=T)<0.999)
      {
        
        ###modelo BCa
        fit_bca_boot<-try(barma(y_boot_bca,ar=phi1,ma=ma1,diag=0),silent = T)
        #convergenciabca<-fit_bca_boot$conv
        
        phi <- fit_bca_boot$coeff[2:(p1+1)] 
        theta <- (fit_bca_boot$coeff[(p1+2):(p1+q1+1)])*(-1)
        
        names_phi<- rep(0,p)
        for(b in 1:p){
          for(j in 1:(length(ar))){
            
            if(b == ar[j]){
              names_phi[b]<-phi[j]
            }
          }
        }
        
        names_theta<- rep(0,q)
        for(b in 1:q){
          for(j in 1:(length(ma))){
            
            if(b == ma[j]){
              names_theta[b]<-theta[j]
            }
          }
        }
      }
      
      if(m<h){
        
        n<-h-length(names_phi)
        n<-rep(0,n)
        phi<-c(names_phi,n)
        
        n<-h-length(names_theta)
        n<-rep(0,n)
        theta<-c(names_theta,n)
        
        
        
        psi<-rep(NA,h+1)
        psi[1]<-1
        
        for(b in 2:(h+1))
        {
          l<-2+b
          end<-l-(b+1)
          psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end]) - theta[b-1]
        }
        
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          forecast_bca<-prev_y(fit_bca_boot,ar=ar1,ma=ma1,h=h)
          mu_hat_bca1<-forecast_bca$mu_hat
          
          if(min(mu_hat_bca1,na.rm=T) < 0.0001 && max(mu_hat_bca1,na.rm=T) > 0.999 )
          {
            forecast_bca<-prev_y(fit_bca_boot,ma=ma1,h=h)
            mu_hat_bca1<-forecast_bca$mu_hat
          } 
          
          if(min(mu_hat_bca1,na.rm=T) > 0.0001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          {
            mu_hat_bca<-mu_hat_bca1
            
          }
          
          ############ IP BCa
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          psi_new<-psi[2:h]^2
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h
          
        }
      }else{
        
        if(p == q && p>h && q>h){
          theta<-c(names_theta[1:h])
          phi<-c(names_phi[1:h])
        }
        
        if(p != q && p>h){ 
          phi<-c(names_phi[1:h])
        } 
        
        if(p != q && p<h){
          n<-h-length(names_phi)
          n<-rep(0,n)
          phi<-c(names_phi,n)
        }
        
        if(p != q && q>h){
          theta<-c(names_theta[1:h])   
        }
        
        if(p != q && q<h){
          n<-h-length(names_theta)
          n<-rep(0,n)
          theta<-c(names_theta,n)
        }
        
        psi<-rep(NA,h+1)
        psi[1]<-1
        
        for(b in 2:(h+1))
        {
          l<-2+b
          end<-l-(b+1)
          psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end]) - theta[b-1]
        }
        
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          
          forecast_bca<-prev_y(fit_bca_boot,ar=ar1,ma=ma1,h=h)
          mu_hat_bca<-forecast_bca$mu_hat
          
          ############ IP BCa
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          psi_new<-psi[2:h]^2
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h
        }
      }
      
      resid<-sample(resid1,h1,replace=T)
      
      mustar_prev<-digamma(fit$forecast*prec) - digamma((1-fit$forecast)*prec)
      nuboot<-trigamma(fit$forecast*prec)+trigamma((1-fit$forecast)*prec)
      mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      nuboot2<- trigamma(mu_hat_bca*prec) + trigamma((1-mu_hat_bca)*prec)
      y_more<- exp(mustar_prev + resid*sqrt(nuboot))/(1 + exp(mustar_prev+ resid*sqrt(nuboot)))
      ystar_boot<- log(y_more/(1-y_more))
      resid_more<- (ystar_boot-mustar_boot)/sqrt(nuboot2)  
      M_resid_bca <- rbind(M_resid_bca, resid_more)
      tc_bca<- tc_bca + (resid_more<median(resid1)) # se prev boot e menoss que prev conta +1         
      
      
      
    } ### fim do laco boot
    
    
    ############ IP BCa
    phi <- fit$coeff[2:(p1+1)] 
    theta <- (fit$coeff[(p1+2):(p1+q1+1)])*(-1)
    
    names_phi<- rep(0,p)
    for(b in 1:p){
      for(j in 1:(length(ar))){
        
        if(b == ar[j]){
          names_phi[b]<-phi[j]
        }
      }
    }
    
    names_theta<- rep(0,q)
    for(b in 1:q){
      for(j in 1:(length(ma))){
        
        if(b == ma[j]){
          names_theta[b]<-theta[j]
        }
      }
    }
    
    if(m<h){
      
      n<-h-length(names_phi)
      n<-rep(0,n)
      phi<-c(names_phi,n)
      
      n<-h-length(names_theta)
      n<-rep(0,n)
      theta<-c(names_theta,n)
      
      psi<-rep(NA,h+1)
      psi[1]<-1
      
      for(b in 2:(h+1))
      {
        l<-2+b
        end<-l-(b+1)
        psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end]) - theta[b-1]
      }
      
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      psi_new<-psi[2:h]^2
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
    } else{
      
      if(p == q && p>h && q>h){
        theta<-c(names_theta[1:h])
        phi<-c(names_phi[1:h])
      }
      
      if(p != q && p>h){ 
        phi<-c(names_phi[1:h])
      } 
      
      if(p != q && p<h){
        n<-h-length(names_phi)
        n<-rep(0,n)
        phi<-c(names_phi,n)
      }
      
      if(p != q && q>h){
        theta<-c(names_theta[1:h])   
      }
      
      if(p != q && q<h){
        n<-h-length(names_theta)
        n<-rep(0,n)
        theta<-c(names_theta,n)
      }
      
      psi<-rep(NA,h+1)
      psi[1]<-1
      
      for(b in 2:(h+1))
      {
        l<-2+b
        end<-l-(b+1)
        psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end]) - theta[b-1]
      }
      
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      psi_new<-psi[2:h]^2
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
    }
    
    nu_zero<-qnorm(tc_bca/B)
    alpha_i<-pnorm(nu_zero + ((nu_zero-z)/(1-a*(nu_zero-z))))
    alpha_s<-pnorm(nu_zero + ((nu_zero+z)/(1-a*(nu_zero+z))))
    li_bca<-apply(M_resid_bca,2,quantile,alpha_i)
    ls_bca<- apply(M_resid_bca,2,quantile,alpha_s)
    mustar_more<-digamma(fit$forecast*prec_fit) - digamma((1-fit$forecast)*prec_fit)
    numore<-trigamma(fit$forecast*prec_fit)+trigamma((1-fit$forecast)*prec_fit)
    
    LI_bca<- exp(mustar_more+diag(li_bca)*sqrt(numore))/(1+exp(mustar_more+diag(li_bca)*sqrt(numore)))
    LS_bca<- exp(mustar_more+diag(ls_bca)*sqrt(numore))/(1+exp(mustar_more+diag(ls_bca)*sqrt(numore)))
    
    
    z<-c()
    
    z$li_bca<-LI_bca
    z$ls_bca<-LS_bca
    
    prev_presentation <- cbind(round(ym,4),round(fit$forecast,4),round(z$li_bca,4),round(z$ls_bca,4))
    colnames(prev_presentation)<-c("Real","Forecast","Lo","Hi")
    
    z$prev <- prev_presentation
    print(prev_presentation)
    
    fim<-end(y)[1]+end(y)[2]/12
    li<-ts(z$li_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    ls<-ts(z$ls_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    
    par(mfrow=c(1,1))
    plot(y,type="l", ylim=c(min(c(y,li)),max(c(y,ls))))
    abline(v=fim,lty=2)
    lines(li,col="blue")
    lines(ls,col="blue")
    lines(ym,lty=2)
  }
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T))
  {
    p1 <- length(ar)
    p<- max(ar)
    m <- max(p,na.rm=T)
    n <- length(fit$serie)
    n_first<-n+h
    phi1<-ar
    h1 <- h
    
    #########################BCa#########################
    tc_bca<-0
    ystar <- (log(fit$serie/(1-fit$serie)))[(m+1):n]
    mustar <- (digamma(fit$fitted*fit$prec) - digamma((1-fit$fitted)*fit$prec))[(m+1):n]
    nu<- (trigamma(fit$fitted*fit$prec) + trigamma((1-fit$fitted)*fit$prec))[(m+1):n]
    resid1<- (ystar-mustar)/sqrt(nu)
    phi_more<- (fit$prec^3)*(psigamma((fit$forecast*fit$prec),deriv=2) - psigamma(((1-fit$forecast)*fit$prec),deriv=2))
    nu_more<- (fit$prec^2)*(trigamma(fit$forecast*fit$prec)+trigamma((1-fit$forecast)*fit$prec))
    a<- ((1/6) * (phi_more))/(nu_more^(3/2))
    
    #####################################################
    i<-0
    while(i<B)
    {
      
      n <- length(fit$serie)-m
      ### amostra BCa
      resid2<-sample(resid1,(n),replace=T)
      ystar_boot_bca<- (exp(mustar+ resid2*sqrt(nu))) / (1 + exp(mustar+ resid2*sqrt(nu)))
      y_boot_bca<-ts(ystar_boot_bca)
      
      if(min(y_boot_bca,na.rm=T)>0.001 && max(y_boot_bca,na.rm=T)<0.999)
      {
        
        ###modelo BCa
        fit_bca_boot<-try(barma(y_boot_bca,ar=phi1,diag=0),silent = T)
        #convergenciabca<-fit_bca_boot$conv
        
        phi <- fit_bca_boot$coeff[2:(p1+1)] 
        
        names_phi<- rep(0,p)
        
        for(b in 1:p){
          for(j in 1:(length(ar))){
            
            if(b == ar[j]){
              names_phi[b]<-phi[j]
            }
          }
        }
        

      
      if(p<=h){
        
        n1<-h-length(names_phi)
        n1<-rep(0,n1)
        phi<-c(names_phi,n1)
        psi<-rep(NA,h+1)
        psi[1]<-1
        
        for(b in 2:(h+1))
        {
          l<-2+b
          end<-l-(b+1)
          psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end])
          
        }
        
  
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          forecast_bca<-prev_y(fit_bca_boot,ar=phi1,h=h)
          mu_hat_bca1<-forecast_bca$mu_hat
          
          if(min(mu_hat_bca1,na.rm=T) > 0.0001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          {
            forecast_bca<-prev_y(fit_bca_boot,ar=phi1,h=h)
            mu_hat_bca1<-forecast_bca$mu_hat
            
          } 
          
          if(min(mu_hat_bca1,na.rm=T) > 0.001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          {
            mu_hat_bca<-mu_hat_bca1
            
          }
          
          ############ IP BCa
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          psi_new<-psi[2:h]^2
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h
          
        }
      } else{
        phi<-c(names_phi[1:h])
        psi<-rep(NA,h+1)
        psi[1]<-1
        
        for(b in 2:(h+1))
        {
          l<-2+b
          end<-l-(b+1)
          psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end])
          
        }
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          forecast_bca<-prev_y(fit_bca_boot,ar=phi1,h=h)
          mu_hat_bca<-forecast_bca$mu_hat

          
          # if(min(mu_hat_bca1,na.rm=T) < 0.0001 && max(mu_hat_bca1,na.rm=T) > 0.999 )
          # {
          #   forecast_bca<-prev_y(fit_bca_boot,ar=phi1,h=h)
          #   mu_hat_bca1<-forecast_bca$mu_hat
          # } 
          # 
          # if(min(mu_hat_bca1,na.rm=T) > 0.0001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          # {
          #   mu_hat_bca<-mu_hat_bca1
          #   
          # }
          
          ############ IP BCa
          
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          psi_new<-psi[2:h]^2
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h
          
        }
      }
    
      #     mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      #     nuboot<- trigamma(mu_hat_bca*prec) + trigamma((1-mu_hat_bca)*prec)
      #     resid<-sample(resid1,h,replace=T)
      #     
      #     y_more<- exp(mustar_boot+ resid*sqrt(nuboot))/(1 + exp(mustar_boot+ resid*sqrt(nuboot)))
      #     ystar_boot<- log(y_more/(1-y_more))
      #     resid_more<- (ystar_boot-mustar_boot)/sqrt(nuboot) 
      #     M_resid_bca <- rbind(M_resid_bca, resid_more)
      
      #         mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      #         nuboot<- trigamma(mu_hat_bca*prec) + trigamma((1-mu_hat_bca)*prec)
      #         resid<-sample(resid1,h1,replace=T)
      #         
      #         y_more<- exp(mustar_boot+ resid*sqrt(nuboot))/(1 + exp(mustar_boot+ resid*sqrt(nuboot)))
      #         ystar_boot<- log(y_more/(1-y_more))
      #         resid_more<- (ystar_boot-mustar_boot)/sqrt(nuboot) 
      #         M_resid_bca <- rbind(M_resid_bca, resid_more)
      
      #mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      
      resid<-sample(resid1,h,replace=T)
  
      
      mustar_prev<-digamma(fit$forecast*prec) - digamma((1-fit$forecast)*prec)
      nuboot<-trigamma(fit$forecast*prec)+trigamma((1-fit$forecast)*prec)
      
      mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      nuboot2<- trigamma(mu_hat_bca*prec) + trigamma((1-mu_hat_bca)*prec)
      
      #y_more<- exp(mustar_boot+ resid*sqrt(nuboot))/(1 + exp(mustar_boot+ resid*sqrt(nuboot)))
      y_more<- exp(mustar_prev + resid*sqrt(nuboot))/(1 + exp(mustar_prev+ resid*sqrt(nuboot)))
      ystar_boot<- log(y_more/(1-y_more))
      resid_more<- (ystar_boot-mustar_boot)/sqrt(nuboot2) 
      #resid_more<- (ystar_boot-mustar_prev)/sqrt(nuboot) 
      M_resid_bca <- rbind(M_resid_bca, resid_more)
      tc_bca<- tc_bca + (resid_more<median(resid1)) # se prev boot e menoss que prev conta +1        
      
    
      
    } ### fim do laco boot
      
    }    
    
    ############ IP BCa
    p1 <- length(ar)
    p<- max(ar)
    
    phi <- fit$coeff[2:(p1+1)] 
    
    names_phi<- rep(0,p)
    for(b in 1:p){
      for(j in 1:(length(ar))){
        
        if(b == ar[j]){
          names_phi[b]<-phi[j]
        }
      }
    }
    if(p<=h){
      
      n<-h-length(names_phi)
      n<-rep(0,n)
      phi<-c(names_phi,n)
      psi<-rep(NA,h+1)
      psi[1]<-1
      
      for(b in 2:(h+1))
      {
        l<-2+b
        end<-l-(b+1)
        psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end])
        
      }
      
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      psi_new<-psi[2:h]^2
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
      
    }else{
      n<-h
      n<-rep(0,n)
      phi<-names_phi[1:h]
      psi<-rep(NA,h+1)
      psi[1]<-1
      
      for(b in 2:(h+1))
      {
        l<-2+b
        end<-l-(b+1)
        psi[b]<-sum(phi[1:(b-1)]*psi[(b-1):end])
        
      }
      
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      psi_new<-psi[2:h]^2
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
    }
    
    nu_zero<-qnorm(tc_bca/B)
    alpha_i<-pnorm(nu_zero + ((nu_zero-z)/(1-a*(nu_zero-z))))
    alpha_s<-pnorm(nu_zero + ((nu_zero+z)/(1-a*(nu_zero+z))))
    li_bca<-apply(M_resid_bca,2,quantile,alpha_i)
    ls_bca<- apply(M_resid_bca,2,quantile,alpha_s)
    mustar_more<-digamma(fit$forecast*prec_fit) - digamma((1-fit$forecast)*prec_fit)
    numore<-trigamma(fit$forecast*prec_fit)+trigamma((1-fit$forecast)*prec_fit)
    
    if ( h == 1)
    {
      
      LI_bca<- exp(mustar_more+(li_bca)*sqrt(numore))/(1+exp(mustar_more+(li_bca)*sqrt(numore)))
      LS_bca<- exp(mustar_more+(ls_bca)*sqrt(numore))/(1+exp(mustar_more+(ls_bca)*sqrt(numore)))
    }else
    {
      LI_bca<- exp(mustar_more+diag(li_bca)*sqrt(numore))/(1+exp(mustar_more+diag(li_bca)*sqrt(numore)))
      LS_bca<- exp(mustar_more+diag(ls_bca)*sqrt(numore))/(1+exp(mustar_more+diag(ls_bca)*sqrt(numore)))   
    }
    
    
    
    
    z<-c()
    
    z$li_bca<-LI_bca
    z$ls_bca<-LS_bca
    
    prev_presentation <- cbind(round(ym,4),round(fit$forecast,4),round(z$li_bca,4),round(z$ls_bca,4))
    colnames(prev_presentation)<-c("Real","Forecast","Lo","Hi")
    
    z$prev <- prev_presentation
    print(prev_presentation)
    
    fim<-end(y)[1]+end(y)[2]/12
    li<-ts(z$li_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    ls<-ts(z$ls_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    
    par(mfrow=c(1,1))
    plot(y,type="l", ylim=c(min(c(y,li)),max(c(y,ls))))
    abline(v=fim,lty=2)
    lines(li,col="blue")
    lines(ls,col="blue")
    lines(ym,lty=2)
  }
  
  ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F))
  { 
    q1 <- length(ma)
    q  <- max(ma)
    m <- max(q,na.rm=T)
    n <- length(fit$serie)-m
    n_first<-n+h
    ma1<-ma
    h1 <- h
    
    #########################BCa#########################
    tc_bca<-0
    ystar <- (log(fit$serie/(1-fit$serie)))[(m+1):n]
    mustar <- (digamma(fit$fitted*fit$prec) - digamma((1-fit$fitted)*fit$prec))[(m+1):n]
    nu<- (trigamma(fit$fitted*fit$prec) + trigamma((1-fit$fitted)*fit$prec))[(m+1):n]
    resid1<- (ystar-mustar)/sqrt(nu)
    phi_more<- (fit$prec^3)*(psigamma((fit$forecast*fit$prec),deriv=2) - psigamma(((1-fit$forecast)*fit$prec),deriv=2))
    nu_more<- (fit$prec^2)*(trigamma(fit$forecast*fit$prec)+trigamma((1-fit$forecast)*fit$prec))
    a<- ((1/6) * (phi_more))/(nu_more^(3/2))
    
    #####################################################
    i<-0
    while(i<B)
    {
      
      n <- length(fit$serie)-m
      ### amostra BCa
      resid2<-sample(resid1,(n-m),replace=T)
      ystar_boot_bca<- (exp(mustar+ resid2*sqrt(nu))) / (1 + exp(mustar+ resid2*sqrt(nu)))
      y_boot_bca<-ts(ystar_boot_bca)
      
      if(min(y_boot_bca,na.rm=T)>0.0001 && max(y_boot_bca,na.rm=T)<0.999)
      {
        
        ###modelo BCa
        fit_bca_boot<-try(barma(y_boot_bca,ma=ma1,diag=0),silent = T)
        #convergenciabca<-fit_bca_boot$conv
        
        theta <- (fit_bca_boot$coeff[2:(q1+1)])*(-1)
        
        names_theta<- rep(0,q)
        for(b in 1:q){
          for(j in 1:(length(ma))){
            
            if(b == ma[j]){
              names_theta[b]<-theta[j]
            }
          }
        }
      }
      
      if(q<h){
        
        psi<-rep(NA,q)
        
        for(b in 1:q)
        {
          if(names_theta[b] != 0)
          {
            psi[b]<-- names_theta[b]
          }else{
            psi[b]<- 0
          }
          
        }
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          forecast_bca<-prev_y(fit_bca_boot,ma=ma1,h=h)
          mu_hat_bca1<-forecast_bca$mu_hat
          
          if(min(mu_hat_bca1,na.rm=T) < 0.0001 && max(mu_hat_bca1,na.rm=T) > 0.999 )
          {
            forecast_bca<-prev_y(fit_bca_boot,ma=ma1,h=h)
            mu_hat_bca1<-forecast_bca$mu_hat
          } 
          
          if(min(mu_hat_bca1,na.rm=T) > 0.0001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          {
            mu_hat_bca<-mu_hat_bca1
            
          }
          
          ############ IP BCa
          
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          n<-h-length(psi)
          n<-rep(0,n)
          psi_new<-c(psi^2,n)
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h 
        }
      } else{
        psi<-rep(NA,h)
        names_theta<-names_theta[1:h]
        for(b in 1:h)
        {
          if(names_theta[b] != 0)
          {
            psi[b]<-- names_theta[b]
          }else{
            psi[b]<- 0
          }
          
        }
        
        if(fit_bca_boot$conv == 0) {
          i<-i+1
          
          if (i==B*0.10)
            print("10%")
          if (i==B*0.20)
            print("20%")
          if (i==B*0.3)
            print("30%")
          if (i==B*0.4)
            print("40%")
          if (i==B*0.5)
            print("50%")
          if (i==B*0.6)
            print("60%")
          if (i==B*0.7)
            print("70%")
          if (i==B*0.8)
            print("80%")
          if (i==B*0.9)
            print("90%")
          
          ### valores previstos
          forecast_bca<-prev_y(fit_bca_boot,ma=ma1,h=h)
          mu_hat_bca1<-forecast_bca$mu_hat
          
          if(min(mu_hat_bca1,na.rm=T) < 0.0001 && max(mu_hat_bca1,na.rm=T) > 0.999 )
          {
            forecast_bca<-prev_y(fit_bca_boot,ma=ma1,h=h)
            mu_hat_bca1<-forecast_bca$mu_hat
          } 
          
          if(min(mu_hat_bca1,na.rm=T) > 0.0001 && max(mu_hat_bca1,na.rm=T) < 0.999 )
          {
            mu_hat_bca<-mu_hat_bca1
            
          }
          
          ############ IP BCa
          
          v_mu<- mu_hat_bca*(1-mu_hat_bca)
          sigma<- (mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)/(1+fit_bca_boot$prec)) 
          psi_new<-c(psi^2)
          v_h<-rep(0,h)
          for(b in 1:h){
            if(b==1){
              
              v_h[b]<-sigma[b]
            }else{
              v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
            }
          }
          prec<-(((mu.eta(mu_hat_bca)^{-2})*(mu_hat_bca*(1-mu_hat_bca)))-v_h)/v_h 
        }
      }  
      mustar_boot<- digamma(mu_hat_bca*prec) - digamma((1-mu_hat_bca)*prec) 
      nuboot<- trigamma(mu_hat_bca*prec) + trigamma((1-mu_hat_bca)*prec)
      resid<-sample(resid1,h,replace=T)
      
      y_more<- exp(mustar_prev+ resid*sqrt(nuboot))/(1 + exp(mustar_prev+ resid*sqrt(nuboot)))
      ystar_boot<- log(y_more/(1-y_more))
      resid_more<- (ystar_boot-mustar_boot)/sqrt(nuboot) 
      M_resid_bca <- rbind(M_resid_bca, resid_more)
      tc_bca<- tc_bca + (resid_more<median(resid1)) # se prev boot e menoss que prev conta +1         
      
    } ### fim do laco boot
    
    
    ############ IP BCa
    q1 <- length(ma)
    q  <- max(ma)
    
    theta <- (fit$coeff[2:(q1+1)])*(-1)
    
    names_theta<- rep(0,q)
    for(b in 1:q){
      for(j in 1:(length(ma))){
        
        if(b == ma[j]){
          names_theta[b]<-theta[j]
        }
      }
    }
    
    if(q<h){
      
      psi<-rep(NA,q)
      
      for(b in 1:q)
      {
        if(names_theta[b] != 0)
        {
          psi[b]<-- names_theta[b]
        }else{
          psi[b]<- 0
        }
        
      }
      
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      n<-h-length(psi)
      n<-rep(0,n)
      psi_new<-c(psi^2,n)
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
    } else{
      psi<-rep(NA,h)
      names_theta<-names_theta[1:h]
      for(b in 1:h)
      {
        if(names_theta[b] != 0)
        {
          psi[b]<-- names_theta[b]
        }else{
          psi[b]<- 0
        }
      }
      v_mu<- fit$forecast*(1-fit$forecast)
      sigma<- (mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)/(1+fit$prec)) 
      psi_new<-c(psi^2)
      v_h<-rep(0,h)
      for(b in 1:h){
        if(b==1){
          
          v_h[b]<-sigma[b]
        }else{
          v_h[b]<- (1+sum(psi_new[1:(b-1)]))*sigma[b]
        }
      }
      prec_fit<-(((mu.eta(fit$forecast)^{-2})*(fit$forecast*(1-fit$forecast)))-v_h)/v_h
      
    }
    
    nu_zero<-qnorm(tc_bca/B)
    alpha_i<-pnorm(nu_zero + ((nu_zero-z)/(1-a*(nu_zero-z))))
    alpha_s<-pnorm(nu_zero + ((nu_zero+z)/(1-a*(nu_zero+z))))
    li_bca<-apply(M_resid_bca,2,quantile,alpha_i)
    ls_bca<- apply(M_resid_bca,2,quantile,alpha_s)
    mustar_more<-digamma(fit$forecast*prec_fit) - digamma((1-fit$forecast)*prec_fit)
    numore<-trigamma(fit$forecast*prec_fit)+trigamma((1-fit$forecast)*prec_fit)
    
    LI_bca<- exp(mustar_more+diag(li_bca)*sqrt(numore))/(1+exp(mustar_more+diag(li_bca)*sqrt(numore)))
    LS_bca<- exp(mustar_more+diag(ls_bca)*sqrt(numore))/(1+exp(mustar_more+diag(ls_bca)*sqrt(numore)))
    
   
    z<-c()
    
    z$li_bca<-LI_bca
    z$ls_bca<-LS_bca
    
    prev_presentation <- cbind(round(ym,4),round(fit$forecast,4),round(z$li_bca,4),round(z$ls_bca,4))
    colnames(prev_presentation)<-c("Real","Forecast","Lo","Hi")
    
    z$prev <- prev_presentation
    print(prev_presentation)
    
    fim<-end(y)[1]+end(y)[2]/12
    li<-ts(z$li_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    ls<-ts(z$ls_bca, start=(end(y)+c(0,1)), frequency=frequency(y))
    
    par(mfrow=c(1,1))
    plot(y,type="l", ylim=c(min(c(y,li)),max(c(y,ls))))
    abline(v=fim,lty=2)
    lines(li,col="blue")
    lines(ls,col="blue")
    lines(ym,lty=2)
  }
  
  
  
  
  
  return(z)
  
}