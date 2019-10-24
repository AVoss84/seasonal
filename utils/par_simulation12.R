
#------------------------------------##
## Simulate PAR(1) process (S=12):
##-----------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

par1.12 = function(N=150, mu=c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,mu12), beta=0, phi=c(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,phi12), sig=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12), g=TRUE,seed,start=NULL,burnin=500,ur=F,...)
{
  if(!missing(seed)){set.seed(seed)} ;   
  if(is.null(start)){start=begin=1} else{(begin = ifelse(length(start)==2,yes=start[2],no=1))}
  if(missing(sig)){sig=rep(1,12)}                                             #note that estimation procedure does not account for unequal variances
  n = N+burnin+begin ; y=ts(numeric(n),fre=12); y[1]=0
  if(ur){phi[12]=prod(phi[-12])^(-1); message(' Unit root restriction imposed!\n')}     #impose the null restriction (see Boswijk/Franses,p.244,12)
  stopifnot(length(beta)==1) ; print(prod(phi))  
    
  for(tt in 2:n){
    st = 1+trunc((tt-1)%%12); year = 1+trunc((tt-1)/12);
    y[tt] = beta*year + 
      (st==1)*(mu[1]+phi[1]*y[tt-1]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi[2]*y[tt-1]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi[3]*y[tt-1]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi[4]*y[tt-1]+rnorm(1,mean=0,sd=sig[4])) +
      (st==5)*(mu[5]+phi[5]*y[tt-1]+rnorm(1,mean=0,sd=sig[5])) +
      (st==6)*(mu[6]+phi[6]*y[tt-1]+rnorm(1,mean=0,sd=sig[6])) +
      (st==7)*(mu[7]+phi[7]*y[tt-1]+rnorm(1,mean=0,sd=sig[7])) +
      (st==8)*(mu[8]+phi[8]*y[tt-1]+rnorm(1,mean=0,sd=sig[8])) +
      (st==9)*(mu[9]+phi[9]*y[tt-1]+rnorm(1,mean=0,sd=sig[9])) +
      (st==10)*(mu[10]+phi[10]*y[tt-1]+rnorm(1,mean=0,sd=sig[10])) +
      (st==11)*(mu[11]+phi[11]*y[tt-1]+rnorm(1,mean=0,sd=sig[11])) +
      (st==12)*(mu[12]+phi[12]*y[tt-1]+rnorm(1,mean=0,sd=sig[12]))
  }  
  (take.low = which(cycle(y)==begin)[10]) ;       #take the draws for example after 10 cycles
  (take.up = take.low+N-1);
  (yy = ts(y[take.low:take.up],fre=12,start))
 print(abs(prod(phi)))  

  if(g){
    graphics.off()
    par(las=1,cex.main=1.2, font.main=11)
    plot(yy,type="l",ylab=expression(y[t])) ; 
    title(paste(ifelse(abs(prod(phi))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(1) process"),sub=expression(S==12));
  }
  return(list(y=yy,phi=phi,sig=sig,mu=mu,ur.check=abs(prod(phi))))
}
dput(par1.12,"par1.12.R")
#
#set.seed(12)
(out = par1.12(N=150,mu=rpois(n=12,lambda=1.5), beta=0, start=c(1950,4), phi=c(.5,.7,1.5,1,.5,.56,1.55,1.42,1.22,1,1.12,1.33),sig=rep(1,12))$y)  #sig=rgamma(n=12,shape=2)



#----------------------------##
## Simulate PAR(2) processes:
##---------------------------##
## y_{s,\tau} = \alpha_{s,1} * y_{s-1,tau} + \alpha_{s,2} * y_{s-2,tau} + \epsilon_{s,\tau}  

par2.12 = function(N=150, mu=c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,mu12), beta=0, phi=c(phi11,phi12,phi13,phi14,phi15,phi16,phi17,phi18,phi19,phi110,phi111,phi112,
       phi21,phi22,phi23,phi24,phi25,phi26,phi27,phi28,phi29,phi210,phi211,phi212), sig=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12), g=TRUE,seed,start=NULL,burnin=500,...)
{
  if(!missing(seed)){set.seed(seed)} ;   
  if(is.null(start)){start=begin=1} else{(begin = ifelse(length(start)==2, yes=start[2],no=1))}
  if(missing(sig)){sig=rep(1,12)}                                             #note that estimation procedure does not account for unequal variances
  n = N+burnin+begin ; y=ts(numeric(n),fre=12); 
  
  phi11=phi[1];phi12=phi[2];phi13=phi[3];phi14=phi[4]          #lag1
  phi15=phi[5];phi16=phi[6];phi17=phi[7];phi18=phi[8];
  phi19=phi[9];phi110=phi[10];phi111=phi[11];phi112=phi[12]
  phi21=phi[13];phi22=phi[14];phi23=phi[15];phi24=phi[16]          #lag2
  phi25=phi[17];phi26=phi[18];phi27=phi[19];phi28=phi[20]
  phi29=phi[21];phi210=phi[22];phi211=phi[23];phi212=phi[24]

   stopifnot(length(beta)==1) ;y[1:2]=c(0,0)      		#initialize for s=1,2
  for(tt in 3:n){
    st = 1+trunc((tt-1)%%12); year = 1+trunc((tt-1)/12);
    y[tt] = beta*year + 
      (st==1)*(mu[1]+phi11*y[tt-1]+phi21*y[tt-2]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi12*y[tt-1]+phi22*y[tt-2]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi13*y[tt-1]+phi23*y[tt-2]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi14*y[tt-1]+phi24*y[tt-2]+rnorm(1,mean=0,sd=sig[4])) +
      (st==5)*(mu[5]+phi15*y[tt-1]+phi25*y[tt-2]+rnorm(1,mean=0,sd=sig[5])) +
      (st==6)*(mu[6]+phi16*y[tt-1]+phi26*y[tt-2]+rnorm(1,mean=0,sd=sig[6])) +
      (st==7)*(mu[7]+phi17*y[tt-1]+phi27*y[tt-2]+rnorm(1,mean=0,sd=sig[7])) +
      (st==8)*(mu[8]+phi18*y[tt-1]+phi28*y[tt-2]+rnorm(1,mean=0,sd=sig[8])) +
      (st==9)*(mu[9]+phi19*y[tt-1]+phi29*y[tt-2]+rnorm(1,mean=0,sd=sig[9])) +
      (st==10)*(mu[10]+phi110*y[tt-1]+phi210*y[tt-2]+rnorm(1,mean=0,sd=sig[10])) +
      (st==11)*(mu[11]+phi111*y[tt-1]+phi211*y[tt-2]+rnorm(1,mean=0,sd=sig[11])) +
      (st==12)*(mu[12]+phi112*y[tt-1]+phi212*y[tt-2]+rnorm(1,mean=0,sd=sig[12]))
  }  
  (take.low = which(cycle(y)==begin)[10]) ;       #take the draws for example after 10 cycles
  (take.up = take.low+N-1);
  (yy = ts(y[take.low:take.up],fre=12,start))
  
  if(g){
    graphics.off()
    par(las=1,cex.main=1.2, font.main=11)
    plot(yy,type="l",ylab=expression(y[t])) ; 
    title(paste(ifelse(abs(prod(phi[1:12]))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(2) process"),sub=expression(S==12));
  }
 yy;
}
dput(par2.12,"par2.12.R")
#

#(out = par2.12(N=150,mu=rep(0.2,12),start=c(1976,9),beta=0, phi=c(rep(0.2,6),rep(0.35,6),rep(0.4,6),rep(0.35,6)),sig=rep(1,12)))


##----------------------------------------##
## Simulate broken PAR(1) processes (S=12):
##----------------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

broken.par1.12 = function(n=1000,Tb,mu=c(mu1,mu1b,mu2,mu2b,mu3,mu3b,mu4,mu4b,mu5,mu5b,mu6,mu6b,mu7,mu7b,mu8,mu8b,mu9,mu9b,mu10,mu10b,mu11,mu11b,mu12,mu12b), 
                          alpha=c(a1,a1b,a2,a2b,a3,a3b,a4,a4b,a5,a5b,a6,a6b,a7,a7b,a8,a8b,a9,a9b,a10,a10b,a11,a11b,a12,a12b), 
                          phi=c(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,phi12), 
                          sig=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12), g=T,seed,ur=F,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(Tb)){(Tb=trunc(n*.5))} else{stopifnot(Tb<n)}
  if(missing(sig)){sig=rep(1,12)}                      #note that estimation procedure does not account for unequal variances
  y=numeric(n); y[1]=0 ;
  if(ur){phi[12]=prod(phi[-12])^(-1); message(' Unit root restriction imposed!\n')}     #impose the null restriction (see Boswijk/Franses,p.244,12)
  
  for(tt in 2:n){
    st = 1+trunc((tt-1)%%12); year = 1+trunc((tt-1)/12);
    y[tt] =  
      (st==1)*(mu[1]*(tt<=Tb) + mu[2]*(tt>Tb) + alpha[1]*(tt<=Tb)*year + alpha[2]*(tt>Tb)*year +  
      phi[1]*y[tt-1]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[3]*(tt<=Tb) + mu[4]*(tt>Tb) + alpha[3]*(tt<=Tb)*year + alpha[4]*(tt>Tb)*year + 
      phi[2]*y[tt-1]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[5]*(tt<=Tb) + mu[6]*(tt>Tb) + alpha[5]*(tt<=Tb)*year + alpha[6]*(tt>Tb)*year +
      phi[3]*y[tt-1]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[7]*(tt<=Tb) + mu[8]*(tt>Tb) + alpha[7]*(tt<=Tb)*year + alpha[8]*(tt>Tb)*year +
      phi[4]*y[tt-1]+rnorm(1,mean=0,sd=sig[4])) +    
      (st==5)*(mu[9]*(tt<=Tb) + mu[10]*(tt>Tb) + alpha[9]*(tt<=Tb)*year + alpha[10]*(tt>Tb)*year +
      phi[5]*y[tt-1]+rnorm(1,mean=0,sd=sig[5])) +     
      (st==6)*(mu[11]*(tt<=Tb) + mu[12]*(tt>Tb) + alpha[11]*(tt<=Tb)*year + alpha[12]*(tt>Tb)*year +
      phi[6]*y[tt-1]+rnorm(1,mean=0,sd=sig[6])) +     
      (st==7)*(mu[13]*(tt<=Tb) + mu[14]*(tt>Tb) + alpha[13]*(tt<=Tb)*year + alpha[14]*(tt>Tb)*year +
      phi[7]*y[tt-1]+rnorm(1,mean=0,sd=sig[7])) +
      (st==8)*(mu[15]*(tt<=Tb) + mu[16]*(tt>Tb) + alpha[15]*(tt<=Tb)*year + alpha[16]*(tt>Tb)*year +
      phi[8]*y[tt-1]+rnorm(1,mean=0,sd=sig[8])) +
      (st==9)*(mu[17]*(tt<=Tb) + mu[18]*(tt>Tb) + alpha[17]*(tt<=Tb)*year + alpha[18]*(tt>Tb)*year +
      phi[9]*y[tt-1]+rnorm(1,mean=0,sd=sig[9])) +
      (st==10)*(mu[19]*(tt<=Tb) + mu[20]*(tt>Tb) + alpha[19]*(tt<=Tb)*year + alpha[20]*(tt>Tb)*year +
      phi[10]*y[tt-1]+rnorm(1,mean=0,sd=sig[10])) +
      (st==11)*(mu[21]*(tt<=Tb) + mu[22]*(tt>Tb) + alpha[21]*(tt<=Tb)*year + alpha[22]*(tt>Tb)*year +
      phi[11]*y[tt-1]+rnorm(1,mean=0,sd=sig[11])) +
      (st==12)*(mu[23]*(tt<=Tb) + mu[24]*(tt>Tb) + alpha[23]*(tt<=Tb)*year + alpha[24]*(tt>Tb)*year +
      phi[12]*y[tt-1]+rnorm(1,mean=0,sd=sig[12]))
  }  
 print(abs(prod(phi)))  
  if(g){
    graphics.off(); par(las=1,cex.main=1.2, font.main=11)
    plot(y,type="l",ylab=expression(y[t]),col="black"); abline(v=Tb,lty=3,lwd=1,col="blue4") 
    title(paste(ifelse(abs(prod(phi))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(1) process with break"),sub=expression(S==12));
  }
 return(list(y=ts(y,fre=12,...),Tb=Tb,mu=mu,alpha=alpha,phi=phi,sig=sig,ur.check=abs(prod(phi))))
}
dput(broken.par1.12,"broken.par1.12.R")
##
#(out = broken.par1.12(n=70,ur=F,Tb=40,mu=rep(c(1.5,0.2),12), alpha=rep(c(.05,0.02),12), phi=c(.5,.7,1.5,1,.5,.56,1.55,1.42,1.22,1,1.12,0.33), sig=rep(.1,12))$y)



#################################################################################################


##----------------------------------------##
## Simulate broken PAR(1) processes (S=12):
##----------------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)


 broken.par1.12b = function(N=1000, Tb,mu=c(mu1,mu1b,mu2,mu2b,mu3,mu3b,mu4,mu4b,mu5,mu5b,mu6,mu6b,mu7,mu7b,mu8,mu8b,mu9,mu9b,mu10,mu10b,mu11,mu11b,mu12,mu12b), 
                           alpha=c(a1,a1b,a2,a2b,a3,a3b,a4,a4b,a5,a5b,a6,a6b,a7,a7b,a8,a8b,a9,a9b,a10,a10b,a11,a11b,a12,a12b), 
                           phi=c(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,phi12), 
                           sig=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12), g=TRUE,seed,burnin=500,start=1,...)
 {
   if(!missing(seed)){set.seed(seed)}
   if(missing(Tb)){(Tb=trunc(n*.5))} else{stopifnot(Tb<N)}
   if(missing(sig)){sig=rep(1,12)}                      #note that estimation procedure does not account for unequal variances
   n = N+burnin   
   (begin=ifelse(length(start)==1,yes=1,no=start[2]))
   (y = ts(numeric(n),fre=12))
   (seas = c(cycle(y))); y[1]=0 ; 
   Tb = Tb + burnin                        #shift break date by one lag in the future, because below initial observation is cut off (i.e. Tb is shifted to the past by one lag) 
   
     tt=2;
   for(st in seas[-1]){
     year = 1+trunc((tt-1)/12);
     y[tt] =  
       (st==1)*(mu[1]*(tt<=Tb) + mu[2]*(tt>Tb) + alpha[1]*(tt<=Tb)*year + alpha[2]*(tt>Tb)*year +  
       phi[1]*y[tt-1]+rnorm(1,mean=0,sd=sig[1])) + 
       (st==2)*(mu[3]*(tt<=Tb) + mu[4]*(tt>Tb) + alpha[3]*(tt<=Tb)*year + alpha[4]*(tt>Tb)*year + 
       phi[2]*y[tt-1]+rnorm(1,mean=0,sd=sig[2])) + 
       (st==3)*(mu[5]*(tt<=Tb) + mu[6]*(tt>Tb) + alpha[5]*(tt<=Tb)*year + alpha[6]*(tt>Tb)*year +
       phi[3]*y[tt-1]+rnorm(1,mean=0,sd=sig[3])) + 
       (st==4)*(mu[7]*(tt<=Tb) + mu[8]*(tt>Tb) + alpha[7]*(tt<=Tb)*year + alpha[8]*(tt>Tb)*year +
       phi[4]*y[tt-1]+rnorm(1,mean=0,sd=sig[4])) +    
       (st==5)*(mu[9]*(tt<=Tb) + mu[10]*(tt>Tb) + alpha[9]*(tt<=Tb)*year + alpha[10]*(tt>Tb)*year +
       phi[5]*y[tt-1]+rnorm(1,mean=0,sd=sig[5])) +     
       (st==6)*(mu[11]*(tt<=Tb) + mu[12]*(tt>Tb) + alpha[11]*(tt<=Tb)*year + alpha[12]*(tt>Tb)*year +
       phi[6]*y[tt-1]+rnorm(1,mean=0,sd=sig[6])) +     
       (st==7)*(mu[13]*(tt<=Tb) + mu[14]*(tt>Tb) + alpha[13]*(tt<=Tb)*year + alpha[14]*(tt>Tb)*year +
       phi[7]*y[tt-1]+rnorm(1,mean=0,sd=sig[7])) +
       (st==8)*(mu[15]*(tt<=Tb) + mu[16]*(tt>Tb) + alpha[15]*(tt<=Tb)*year + alpha[16]*(tt>Tb)*year +
       phi[8]*y[tt-1]+rnorm(1,mean=0,sd=sig[8])) +
       (st==9)*(mu[17]*(tt<=Tb) + mu[18]*(tt>Tb) + alpha[17]*(tt<=Tb)*year + alpha[18]*(tt>Tb)*year +
       phi[9]*y[tt-1]+rnorm(1,mean=0,sd=sig[9])) +
       (st==10)*(mu[19]*(tt<=Tb) + mu[20]*(tt>Tb) + alpha[19]*(tt<=Tb)*year + alpha[20]*(tt>Tb)*year +
       phi[10]*y[tt-1]+rnorm(1,mean=0,sd=sig[10])) +
       (st==11)*(mu[21]*(tt<=Tb) + mu[22]*(tt>Tb) + alpha[21]*(tt<=Tb)*year + alpha[22]*(tt>Tb)*year +
       phi[11]*y[tt-1]+rnorm(1,mean=0,sd=sig[11])) +
       (st==12)*(mu[23]*(tt<=Tb) + mu[24]*(tt>Tb) + alpha[23]*(tt<=Tb)*year + alpha[24]*(tt>Tb)*year +
       phi[12]*y[tt-1]+rnorm(1,mean=0,sd=sig[12]))
    tt=tt+1;
   }  
   (without.burnin = c(cycle(y))[-c(1:burnin)])
   (start.index = which.min(which(without.burnin==begin)))
   #(take.low = which(cycle(y)==begin)[floor(burnin/12)]) ;       #take the draws for example after some 'burnin cycles'
   #(yy = ts(y[-c(1:(take.low-1))],fre=12,start))
   (take.low = burnin+start.index)
   (take.up = take.low+N-1)
   #(yy = ts(y[-c(1:(burnin+start.index-1))],fre=12,start))
   (yy = ts(na.omit(y[take.low:take.up]),fre=12,start))

   if(g){
     graphics.off()
     par(las=1,cex.main=1.2, font.main=11)
     plot(y[-c(1:(burnin+start.index-1))],type="l",ylab=expression(y[t]),col="black"); abline(v=Tb-burnin,lty=3,lwd=1,col="blue4") 
     title(paste(ifelse(abs(prod(phi))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(1) process with break"),sub=expression(S==12));
   }
   return(list(y=yy,Tb=Tb,mu=mu,alpha=alpha,phi=phi,sig=sig,ur.check=abs(prod(phi))))
 }
#
dput(broken.par1.12b,"broken.par1.12b.R")
##
#(out = broken.par1.12b(N=200,burnin=200,start=c(1950,12),Tb=20,mu=rep(c(1.5,0.2),12), alpha=rep(c(.05,0.05),12), phi=c(.5,.7,1.5,1,.5,.56,1.55,1.42,1.22,1,1.12,0.33), sig=rep(1,12))$y)
# 





















