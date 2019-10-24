
##-----------------------------------##
## Simulate PAR(1) processes (S=4):
##-----------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

par1 = function(N=150, mu=c(mu1,mu2,mu3,mu4), beta=0, phi=c(phi1,phi2,phi3,phi4), sig=c(s1,s2,s3,s4), g=T,seed,start=NULL,burnin=500,ur=F,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(is.null(start)){start=begin=1} else{(begin = ifelse(length(start)==2,yes=start[2],no=1))}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
  n=N+burnin+begin ; y=ts(numeric(n),fre=4); y[1]=0
  if(ur){phi[4]=prod(phi[-4])^(-1); message(' Unit root restriction imposed!\n')}     #impose the null restriction (see Boswijk/Franses,p.244,12)
  stopifnot(length(beta)==1); print(prod(phi))  

  for(tt in 2:n){
    (st = 1+trunc((tt-1)%%4)); (year = 1+trunc((tt-1)/4));
    (y[tt] = beta*year + 
      (st==1)*(mu[1]+phi[1]*y[tt-1]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi[2]*y[tt-1]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi[3]*y[tt-1]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi[4]*y[tt-1]+rnorm(1,mean=0,sd=sig[4])))     
  }  
  (take.low = which(cycle(y)==begin)[10]) ;       #take the draws for example after 10 cycles
  (take.up = take.low+N-1);
  (yy = ts(y[take.low:take.up],fre=4,start))
  
  (adjoint.1 = matrix(c(1,phi[1]*phi[3]*phi[4],phi[1]*phi[4],phi[1],
                        phi[2],1,phi[1]*phi[2]*phi[4],phi[1]*phi[2],
                        phi[2]*phi[3],phi[3],1,phi[1]*phi[2]*phi[3],
                        phi[2]*phi[3]*phi[4],phi[3]*phi[4],phi[4],1)
                      ,byrow=T,nrow=4,ncol=4))                    #see Adjoint matrix at L=1, in Osborn, p.147 (6.19)
  (M = (adjoint.1/(1-prod(phi)))%*%cbind(mu))                       #p.147, (6.20/21)
  (inv.Phi0 = matrix(c(1,0,0,0,phi[2],1,0,0,phi[2]*phi[3],phi[3],1,0,phi[2]*phi[3]*phi[4],phi[3]*phi[4],phi[4],1),ncol=4,byrow=T))  #see own stuff
  (Phi1 = matrix(c(0,0,0,phi[1],0,0,0,0,0,0,0,0,0,0,0,0),ncol=4,byrow=T))     #see Osborn, p.145
  
  if(g){
    graphics.off()    #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(las=1,cex.main=1.4, font.main=11)
    plot(yy,type="l",ylab=expression(y[t])) ; 
    #abline(h=M,lty=3,col="red",lwd=2)
    title(paste(ifelse(abs(prod(phi))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(1) process"),sub=expression(S==4));
  }
  return(list(y=yy,phi=inv.Phi0%*%Phi1,smeans=M))
}
dput(par1,'par1.R')
#
(out = par1(N=100, ur=F,mu=c(1,1,1,1), beta=.1, phi=c(.5,.7,.5,1),sig=c(1,1,1,1)))


##----------------------------------------##
## Simulate broken PAR(1) processes (S=4):
##----------------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

broken.par1 = function(n=1000, Tb,mu=c(mu1,mu1b,mu2,mu2b,mu3,mu3b,mu4,mu4b), ur=F,
			alpha=c(a1,a1b,a2,a2b,a3,a3b,a4,a4b), phi=c(phi1,phi2,phi3,phi4), sig=c(s1,s2,s3,s4), g=TRUE,seed,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(Tb)){(Tb=trunc(n*.5))} else{stopifnot(Tb<n)}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
   y=numeric(n); y[1]=0
  if(ur){phi[4]=prod(phi[-4])^(-1); message(' Unit root restriction imposed!\n')}     #impose the null restriction (see Boswijk/Franses,p.244,12)
  
  for(tt in 2:n){
    st = 1+trunc((tt-1)%%4); year = 1+trunc((tt-1)/4);
    y[tt] =  
      (st==1)*(mu[1]*(tt<=Tb) + mu[2]*(tt>Tb) + alpha[1]*(tt<=Tb)*year + alpha[2]*(tt>Tb)*year +  
       phi[1]*y[tt-1]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[3]*(tt<=Tb) + mu[4]*(tt>Tb) + alpha[3]*(tt<=Tb)*year + alpha[4]*(tt>Tb)*year + 
       phi[2]*y[tt-1]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[5]*(tt<=Tb) + mu[6]*(tt>Tb) + alpha[5]*(tt<=Tb)*year + alpha[6]*(tt>Tb)*year +
       phi[3]*y[tt-1]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[7]*(tt<=Tb) + mu[8]*(tt>Tb) + alpha[7]*(tt<=Tb)*year + alpha[8]*(tt>Tb)*year +
       phi[4]*y[tt-1]+rnorm(1,mean=0,sd=sig[4]))    
  }  
  if(g){
    graphics.off()
    #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(las=1,cex.main=1.4, font.main=11)
    plot(y,type="l",ylab=expression(y[t]),col="black") ; abline(v=Tb,lty=3,lwd=1,col="red4") 
    title(paste(ifelse(abs(prod(phi))<1,yes="Stationary",no="Nonstationary"),"Gaussian PAR(1) process with break"),sub=expression(S==4));
  }
return(y=ts(y,fre=4,...))
}
#
dput(broken.par1,'broken.par1.R')
#
(out = broken.par1(n=100,Tb=50,mu=rep(c(1,0.1),4), alpha=rep(.01,8), phi=c(.1,.1,.1,.1),sig=c(.1,.1,.1,.1)))


## Simulate PAR(2) processes:
##---------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)
## y_{s,\tau} = \alpha_{s,1} * y_{s-1,tau} + \alpha_{s,2} * y_{s-2,tau} + \epsilon_{s,\tau}  

par2 = function(n=1000, mu=c(mu1,mu2,mu3,mu4), beta=0, phi=c(phi11,phi12,phi13,phi14,phi21,phi22,phi23,phi24), sig=c(s1,s2,s3,s4), g=TRUE,seed,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
  phi11=phi[1];phi12=phi[2];phi13=phi[3];phi14=phi[4]          #lag1
  phi21=phi[5];phi22=phi[6];phi23=phi[7];phi24=phi[8]          #lag2
  stopifnot(length(beta)==1) 

   y=numeric(n); y[1:2]=c(0,0)      		#initialize for s=1,2
  for(tt in 3:n){
    st = 1+trunc((tt-1)%%4); year = 1+trunc((tt-1)/4);
    y[tt] = beta*year + 
      (st==1)*(mu[1]+phi11*y[tt-1]+phi21*y[tt-2]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi12*y[tt-1]+phi22*y[tt-2]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi13*y[tt-1]+phi23*y[tt-2]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi14*y[tt-1]+phi24*y[tt-2]+rnorm(1,mean=0,sd=sig[4]))     
  }  
   y=ts(y,fre=4,...)
  (Phi0 = matrix(c(1,0,0,0,-phi12,1,0,0,-phi23,-phi13,1,0,0,-phi24,-phi14,1),byrow=T,ncol=4,nrow=4))
  (Phi1 = matrix(c(0,0,phi21,phi11,0,0,0,phi22,0,0,0,0,0,0,0,0),byrow=T,ncol=4,nrow=4))
  
  if(g){
    graphics.off()
    #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(las=1,cex.main=1.4, font.main=11); plot(y,type="l",ylab=expression(y[t]))
    title(paste("Gaussian quarterly PAR(2) process"),sub=expression(S==4));
  }
 return(list(Phi0=Phi0,Phi1=Phi1,y=y))
}
#dput(par2,'par2.R')

(y = par2(n=1000, mu=c(0,0,0,0), beta=0, phi=c(.1,.1,.1,.1,.5,.5,.5,.5),sig=c(1,1,1,1),seed=533))

#(out = bayes.par(data=y, p=2,deter=NULL, sdeter="const"))

##--------------------------------------##
## Simulate broken PAR(2) processes:
##--------------------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

broken.par2 = function(n=1000, Tb,mu=c(mu1,mu1b,mu2,mu2b,mu3,mu3b,mu4,mu4b), alpha=c(a1,a1b,a2,a2b,a3,a3b,a4,a4b), 
                       phi=c(phi11,phi12,phi21,phi22,phi31,phi32,phi41,phi42), sig=c(s1,s2,s3,s4), g=TRUE,seed,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
  if(missing(Tb)){(Tb=trunc(n*.5))}else{stopifnot(Tb<n)}
   y=numeric(n); y[1:2]=c(0,0) 
   for(tt in 3:n){
    (st = 1+trunc((tt-1)%%4)); (year = 1+trunc((tt-1)/4));    #Trendmodellierung wie in Franses/Paap, S.30 (impliziert Konstanz ¸bers Jahr hinweg!)
    y[tt] =  
      (st==1)*(mu[1]*(tt<=Tb) + mu[2]*(tt>Tb) + alpha[1]*(tt<=Tb)*year + alpha[2]*(tt>Tb)*year + #oder wie in Bayes-Sch‰tzung und partsm: alpha[2]*(tt>Tb)*tt
       phi[1]*y[tt-1]+phi[2]*y[tt-2]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[3]*(tt<=Tb) + mu[4]*(tt>Tb) + alpha[3]*(tt<=Tb)*year + alpha[4]*(tt>Tb)*year + 
       phi[3]*y[tt-1]+phi[4]*y[tt-2]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[5]*(tt<=Tb) + mu[6]*(tt>Tb) + alpha[5]*(tt<=Tb)*year + alpha[6]*(tt>Tb)*year +
       phi[5]*y[tt-1]+phi[6]*y[tt-2]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[7]*(tt<=Tb) + mu[8]*(tt>Tb) + alpha[7]*(tt<=Tb)*year + alpha[8]*(tt>Tb)*year +
       phi[7]*y[tt-1]+phi[8]*y[tt-2]+rnorm(1,mean=0,sd=sig[4]))  
   }  
  if(g){
    graphics.off()
    par(las=1,cex.main=1.4, font.main=11)
    plot(y,type="l",ylab=expression(y[t]))
    abline(v=Tb,lty=3,col="red4") ; title(paste("Gaussian PAR(2) process with break"),sub=expression(S==4));
  }
return(y=ts(y,fre=4,...))
}
#
(y = broken.par2(n=500, Tb=100,mu=c(1,-.1,1,.1,1,.1,1,.1), alpha=rep(.02,8), phi=c(.1,.1,.1,.1,.5,.5,.5,.5),sig=c(.5,.5,.5,.5)))

#require(pear)
#(ans <- pear(y, m=rep(2,times=4), ic="none")$phi)

#(out = bayes.par(data=y, p=2,deter=NULL, sdeter="const"))

##---------------------------##
## Simulate PAR(3) processes:
##---------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

par3 = function(n=1000, mu=c(mu1,mu2,mu3,mu4), beta=0, phi=c(phi11,phi12,phi13,phi14,phi21,phi22,phi23,phi24,phi31,phi32,phi33,phi34), sig=c(s1,s2,s3,s4), g=TRUE,seed,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
 
  phi11=phi[1];phi12=phi[2];phi13=phi[3];phi14=phi[4]          #lag1
  phi21=phi[5];phi22=phi[6];phi23=phi[7];phi24=phi[8]          #lag2
  phi31=phi[9];phi32=phi[10];phi33=phi[11];phi34=phi[12]          #lag3
  stopifnot(length(beta)==1) 

   y=numeric(n); y[1:3]=rep(0,3)      		#initialize for s=1,2,3
  for(tt in 4:n){
    st = 1+trunc((tt-1)%%4); year = 1+trunc((tt-1)/4);
    y[tt] = beta*year + 
      (st==1)*(mu[1]+phi11*y[tt-1]+phi21*y[tt-2]+phi31*y[tt-3]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi12*y[tt-1]+phi22*y[tt-2]+phi32*y[tt-3]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi13*y[tt-1]+phi23*y[tt-2]+phi33*y[tt-3]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi14*y[tt-1]+phi24*y[tt-2]+phi34*y[tt-3]+rnorm(1,mean=0,sd=sig[4]))     
  }  
   y=ts(y,fre=4,...);
  if(g){
    graphics.off() #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(las=1,cex.main=1.4, font.main=11); plot(y,type="l",ylab=expression(y[t]));
    title(paste("Gaussian quarterly PAR(3) process"),sub=expression(S==4));
  }
 return(y=y)
}
dput(par3,"par3.R")
#
(y = par3(n=1000, mu=c(0,0,0,0), beta=0, phi=c(.1,.1,.1,.1,.5,.5,.5,.5,.2,.2,.2,.2),sig=c(1,1,1,1)))


##---------------------------##
## Simulate PAR(4) processes:
##---------------------------##
## Using the univariate representation in Franses, p.78, (4.60)
## (with seasonal dummy variables)

par4 = function(n=1000, mu=c(mu1,mu2,mu3,mu4), beta=0, phi=c(phi11,phi12,phi13,phi14,phi21,phi22,phi23,phi24,phi31,phi32,phi33,phi34,phi41,phi42,phi43,phi44), sig=c(s1,s2,s3,s4), g=TRUE,seed,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
 
  phi11=phi[1];phi12=phi[2];phi13=phi[3];phi14=phi[4]          #lag1
  phi21=phi[5];phi22=phi[6];phi23=phi[7];phi24=phi[8]          #lag2
  phi31=phi[9];phi32=phi[10];phi33=phi[11];phi34=phi[12]          #lag3
  phi41=phi[13];phi42=phi[14];phi43=phi[15];phi44=phi[16]          #lag4
  stopifnot(length(beta)==1) 

   y=numeric(n); y[1:4]=rep(0,4)      		#initialize for s=1,2,3,4
  for(tt in 5:n){
    st = 1+trunc((tt-1)%%4); year = 1+trunc((tt-1)/4);
    y[tt] = beta*year + 
      (st==1)*(mu[1]+phi11*y[tt-1]+phi21*y[tt-2]+phi31*y[tt-3]+phi41*y[tt-4]+rnorm(1,mean=0,sd=sig[1])) + 
      (st==2)*(mu[2]+phi12*y[tt-1]+phi22*y[tt-2]+phi32*y[tt-3]+phi42*y[tt-4]+rnorm(1,mean=0,sd=sig[2])) + 
      (st==3)*(mu[3]+phi13*y[tt-1]+phi23*y[tt-2]+phi33*y[tt-3]+phi43*y[tt-4]+rnorm(1,mean=0,sd=sig[3])) + 
      (st==4)*(mu[4]+phi14*y[tt-1]+phi24*y[tt-2]+phi34*y[tt-3]+phi44*y[tt-4]+rnorm(1,mean=0,sd=sig[4]))     
  }  
   y=ts(y,fre=4,...);
  if(g){
    graphics.off() #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(las=1,cex.main=1.4, font.main=11); plot(y,type="l",ylab=expression(y[t]));
    title(paste("Gaussian quarterly PAR(4) process"),sub=expression(S==4));
  }
 return(y=y)
}
#
dput(par4,"par4.R")
#
(y = par4(n=1000, mu=c(0,0,0,0), beta=0, phi=c(.1,.1,.1,.1,.5,.5,.5,.5,.2,.2,.2,.2,.02,.02,.02,.02),sig=c(1,1,1,1)))


##-------------------------------------------------------------------------##
## Simulate periodic quarterly moving average process of order 1 PMA(1):
##-------------------------------------------------------------------------##

pma1 = function(N=150, mu=c(mu1,mu2,mu3,mu4), beta=0, theta=c(theta1,theta2,theta3,theta4), sig=c(s1,s2,s3,s4), g=T,seed,start=NULL,burnin=500,...)
{
  if(!missing(seed)){set.seed(seed)}
  if(is.null(start)){start=begin=1} else{(begin = ifelse(length(start)==2,yes=start[2],no=1))}
  if(missing(sig)){sig=rep(1,4)}                      #note that estimation procedure does not account for unequal variances
  (n=N+burnin+begin) ; (tau = ceiling(n/4)); S=4;             #rounded up number of years
  (eps1 = rnorm(tau,mean=0,sd=sig[1])); eps2 = rnorm(tau,mean=0,sd=sig[2]);
  eps3 = rnorm(tau,mean=0,sd=sig[3]); eps4 = rnorm(tau,mean=0,sd=sig[4]);    
  y = ts(numeric(n),fre=4); stopifnot(length(beta)==1)  

  for(tt in (S+1):n){
    (st = 1+trunc((tt-1)%%4)); (year = 1+trunc((tt-1)/4));
    (y[tt] = beta*year + 
      (st==1)*(mu[1] + eps1[year] + theta[1]*eps4[year-1]) + 
      (st==2)*(mu[2] + eps2[year] + theta[2]*eps1[year]) + 
      (st==3)*(mu[3] + eps3[year] + theta[3]*eps2[year]) + 
      (st==4)*(mu[4] + eps4[year] + theta[4]*eps3[year]))    
   }  
  (take.low = which(cycle(y)==begin)[5]) ;                          #take the draws for example after 10 cycles
  (take.up = take.low+N-1); (yy = ts(y[take.low:take.up],fre=4,start))
  if(g){graphics.off(); par(las=1,cex.main=1.4, font.main=11)
    plot(yy,type="l",ylab=expression(y[t])) ; title("Gaussian PMA(1) process",sub=expression(S==4));
  }
 yy;
}
#
dput(pma1,"pma1.R")
#
(out = pma1(N=1000,mu=c(1,1,1,1),beta=0, theta=c(.5,.7,.5,1),sig=c(1,1,1,1)))


##--------------------------------------------------------------------------------------##
## Now use the VAR(1) representation of a PAR(1) process (-> VQ(1) representation)
## see for example Franses, p.64, (4.12), also Osborn and Ghysels

## PAR(1) (passt noch nicht ganz):
##-----------##

par1VQ = function(n=1000, mu=c(mu1,mu2,mu3,mu4), beta=0, phi=c(phi1,phi2,phi3,phi4), sig=c(s1,s2,s3,s4), g=TRUE)
{ 
  p11=phi[1];p12=phi[2];p13=phi[3];p14=phi[4] ; (sim=matrix(rnorm(4),ncol=1,nrow=4))  
  (inv.Phi0 = matrix(c(1,0,0,0,p12,1,0,0,p12*p13,p13,1,0,p12*p13*p14,p13*p14,p14,1),ncol=4,byrow=T))  #see own stuff
  (Phi1 = matrix(c(0,0,0,p11,0,0,0,0,0,0,0,0,0,0,0,0),ncol=4,byrow=T))     #Osborn, p.145  
  
  for(ii in 2:n){
    (eps = matrix(rnorm(4),ncol=1))
    sim = cbind(sim, inv.Phi0%*%Phi1%*%cbind(sim[,ii-1]) + mu + eps)
  }
  VQsim = ts(t(sim),fre=4); colnames(VQsim)=paste("S",1:4,sep="")
  
  if(g){
    #layout(matrix(c(1,1,2,3),ncol=2,byrow=T))
    par(mfrow=c(2,2))
    plot(VQsim[,1],type="l");plot(VQsim[,2],type="l");plot(VQsim[,3],type="l");plot(VQsim[,4],type="l")
    #title(paste("Gauﬂscher PAR(1)-Prozess \n"));
    #acf(y,main="SACF");pacf(y,main="SPACF")
  }
  return(list(y=VQsim,phi=inv.Phi0%*%Phi1))
  #return(list(y=y,acf=acf(y,plot=F),pacf=pacf(y,plot=F)))
}
(y = par1VQ(n=550, mu=c(1,1,1,1), beta=.1, phi=c(.1,.5,.3,.6))$y)






