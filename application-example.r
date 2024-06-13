###############################################################################################################
#
#  Application presented in:            
#  Palm, B.G. and Bayer, F.M. and Cintra, R. J. (2018) 
#  Prediciton Intervals in the Beta Autoregressive Moving Average models
#
###############################################################################################################

source("barma.fit.r")
source("barma.r")
source("ip-bca.r")

set.seed(5)

y_first<-read.table("canta1.txt", h=T) # data
y_first<-y_first[,2]/100
y_first<-ts(y_first) # time series  
h1<-10 # forecasting
link<-"logit"
n_first<-length(y_first) # sample size
n<-n_first-h1
y<-ts(y_first[1:n],start=c(1991,4),frequency=12) 
ym<- c(y_first[(n+1):n_first]) 

ar1<-c(1,2,3,4)
ma1<-c(NA)
fit<-barma(y,ar=ar1,ma=ma1,link = "logit",diag=0,h=h1) # barma model

ipboot<- ip_boot(fit,y=y,ym=ym,ar=ar1,ma=ma1,h=h1,B=1000) # BCa prediction interval





