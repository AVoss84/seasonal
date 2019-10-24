
##================================================================##
## Estimate and predict monthly PAR(p) process via Gibbs Sampling:
##================================================================##

rm(list=ls())

##---------------------------------------------------------------------------------------##
## Simulate some monthly PAR(1) data with(out) structural break:
##---------------------------------------------------------------------------------------##
source("par_simulation12.R")

n=200 ;             #The higher n, the more accurate the forecasts!
(y = par1.12(N=n, mu=rep(1,times=12),start=c(1950,5),beta=-.04, phi=c(.5,.7,1.5,1,.5,.56,1.55,1.42,1.22,1,1.12,1.33),sig=rep(1,12))$y)  #sig=rgamma(n=12,shape=2)

(y = par2.12(N=n,mu=rep(1,12), beta=0, phi=c(rep(0.4,6),rep(0.45,6),rep(0.4,6),rep(0.75,6)),sig=rep(1,12)))

(Tb.true = 105)
(y = broken.par1.12(n=n,Tb=Tb.true, mu=rep(c(.1,.9),12), alpha=rep(c(0,0),12), phi=c(.75,.7,1.5,1,.5,.56,1.05,1.42,1.22,1,1.12,0.33), sig=rep(.5,12))$y)

#(y = broken.par1.12b(N=n,Tb=floor((n-1)*.5),burnin=100, mu=rep(c(1.5,0.2),12), alpha=rep(c(.05,0.02),12), phi=c(.5,.7,1.5,1,.5,.56,1.55,1.42,1.22,1,1.12,0.33), sig=rep(.1,12))$y)


##------------------------------------##
## Anpassen vor Beginn der Session!!
##------------------------------------##
folder = "par2_ob"    				                 #Target folder
label = NULL							                     #Label for simulation design
dir.create(path)						                  #create folder in specified directory with given name
#file = folder						                    #Anfang aller Filenamen, ob Plots oder txt-Files
##------------------------------------##



## Plot series, (non)periodic (P)ACF, periodogram and boxplot of seasonal subseries:
##----------------------------------------------------------------------------------##
pdf(paste(path,"nonperiodic_plots_y.pdf",sep=""))

par(mfrow=c(2,2));
plot(y);
acf(as.vector(y),main="");                        #nonperiodic (homogenous) correlograms
pacf(as.vector(y),main="");
spectrum(as.vector(y),main="")

graphics.off()

#------------------------------------------------------------------------------------#

pdf(paste(path,"periodic_plots_y.pdf",sep=""))
require(pear)

par(mfrow=c(2,2),las=1);
plot(y);
peacf(y);                                 #periodic (heterogenous) correlograms
ppacf = pepacf(y,lag.max=7, plot=T);
attributes(ppacf)
ppacf$bic
peboxplot(y,xlab="period")     
ans = pear(y, ic="bic")           #fit possibly individual model orders       
ans$model.orders
ans$phi

graphics.off()

#-------------------------------------------------------------------------------------------#

gibbs12 = function(y, p , n.ahead, M, burnin, restr2stat=F, sbreak=T, sbreak.fix=F, s2.hyp=c(nu,lambda),b.hyp=c(B=1,delta=0,c1=100,c2=100),in.samp.fc=T,le.W=100,low_up=c(1,1),deter=NULL, sdeter="const",isbreak.nsea=c(1,0),isbreak.sea=c(1,0),tau=frequency(y)-1)  
{             
   library(MASS) ;
  bayes.par12=dget("bayes.par12.R"); strbr.compS = dget("strbr.compS.R"); mats01 = dget('mats01.R'); 

  nseason=deter; season=sdeter; coding=c("mean","effect"); S=12 ; 

  ## Save last n-ahead time observations for in-sample forecast control:
  ##---------------------------------------------------------------------##
  if(frequency(y)!=S){stop('Time series has wrong seasonal frequency!\n')}
  
  if(in.samp.fc){
    (y.true.T = ts(y[c((length(y)-n.ahead+1):length(y))],end=end(y),fre=S));
    (y = ts(y[-c((length(y)-n.ahead+1):length(y))],start=start(y),fre=S))       #actually used sample
   } else{y.true.T=NULL}          

  (n = length(y)) ;(n.new = n-p)                       
  (z = rep(mean(y),times=n.ahead))
  (y.aug = ts(c(y,z),fre=S,start=start(y)))
  (k = c(1,floor((n-1)*.5),n.new))                 #Break interval for the single break date
  #(k = c(1,Tb.true,n.new))                      #Attention: If k does not start from 1, you must modify the vector indices in prob[j..], see step for break date draws
  (date = if(sbreak){k[2]}else{NULL})
  if(sbreak.fix && !sbreak){
    (date = floor(runif(1,min=k[1] + tau + p, max=k[3]-1)));
   message('Draw break date from uniform instrumental density instead of posterior.\n')}    #draw some date from a uniform density
  if(sbreak.fix && sbreak){
    message('Option is not supported. Set "sbreak <- FALSE" if uniform density should be used.\nMultinomial posterior density is used instead.\n')}
    
  #M=10000; burnin=200; 
  MCsim = M+burnin;                                    #Monte Carlo runs and burnin draws
  (out1 = bayes.par12(data=y, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau))
  (lhs = out1@lhs); (X.til = out1@W); (X = X.til[,1:(p*S)]);       #only autoregressive terms ; lhs and rhs

   # For future values/only rhs:
  out2 = bayes.par12(data=y.aug, class.est=T, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau) 
  (W.star = out2@W); (W.star = ts(matrix(W.star[(nrow(W.star)-n.ahead+1):nrow(W.star),],nrow=n.ahead),fre=S,end=end(W.star)))
  
  ##-------------------------## 
  ## Initialize some objects:
  ##-------------------------##
  s2.eps = Tb = prod.unit = numeric(MCsim)
  if(sbreak){Tb[1]=date} else{Tb=NULL}
  b = matrix(,ncol=ncol(X.til),nrow=MCsim)
  W.post = matrix(,ncol=n.ahead,nrow=MCsim)         #posterior draws of (unknown/latent) future values

  max.date = trunc(n/2); trigger = -10^10             #initializations for searching the most likely break date for DIC calculation at the end, see line 323 below 
   
  (W.post[1,]=z) ; (s2.eps[1] = out2@sigma2)          #initialize s2 with classical estimates
  b[1,] = out2@est.classic[,1]
  (OLS = out1@est.classic)                                  #frequentist estimates
  (prod.unit[1] = abs(prod(b[1,1:(p*S)])))           #product of AR coefficients to test for a (monthly) periodic unit root
  margloglike=jlogp.values=NULL ; 
  (W_supp = seq(min(y)-low_up[1]*sd(y),max(y)+low_up[2]*sd(y),length.out=le.W))
  print(W_supp) 		#if(is.null(le.W)){le.W = length(W_supp)}            #is truncated when in-sample pred. was chosen

  marg.Ws = matrix(0,nrow=n.ahead,ncol=length(W_supp),dimnames=list(paste("W",1:n.ahead,sep=""),paste(1:length(W_supp),sep="")))  
  if((MCsim)%%le.W){stop(" 'MCsim' must be a multiple of length(W_supp).\n")
   } else{cat("There are",MCsim/le.W,"draws per point in 'W_supp' used.\n")}

  ## Define all hyperparameters for \beta normal prior (as on p.9):
  ##---------------------------------------------------------------##
  c1 = b.hyp[3] ; c2 = b.hyp[4]                                        #hyperparameters prior covariance matrix of \beta
  (b0 = c(rep(b.hyp[1],times=p*S),rep(b.hyp[2],ncol(X.til)-p*S)))       #prior means
  (inv.V = diag(c(rep(1/c1,times=p*S),rep(1/c2,ncol(X.til)-p*S))))    # V^(-1); see p.19 and p.10 for definitions
  a = s2.hyp[1]; lambda = s2.hyp[2]                                    #parameters IG prior for \sigma^2 
                   #existence of the first four moments if \nu>4
                   #Note: large \lambda induces large prior mean and variance

  ## Function to compute the joint loglikelihood function:
  ##------------------------------------------------------##
   loglike = function(X.til,lhs,s2.eps){
     (beta = ginv(t(X.til)%*%X.til)%*%t(X.til)%*%lhs);     
     (err = lhs - X.til%*%beta); nNew=length(lhs);
     LL = -(nNew/2)*log(2*pi)-(nNew/2)*log(s2.eps)-(1/(2*s2.eps))*t(err)%*%err;  
    return(list(LL=LL,err=as.vector(err)))  
    }
  ## Function to compute the joint logposterior distribution (under Jeffreys prior for s2):
  ##----------------------------------------------------------------------------------------##
  jointpost = function(X.til,lhs,s2.eps,y.til, W.til, B, B0, inv.V){
    (err = lhs - X.til%*%B); nNew=length(lhs); (d = length(B)); (k = length(y.til))
    (pred.err = y.til - W.til%*%B) ; (QF.y = t(err)%*%err) ;
    (QF.B = t(B-B0)%*%inv.V%*%(B-B0)); (QF.fut = t(y.til-W.til%*%as.matrix(B))%*%(y.til-W.til%*%as.matrix(B)))
    #QF.fut=k=0
    #nconst = -.5*(nNew+d+k)*log(2*pi) + .5*log(prod(diag(inv.V)) ; 
    nconst=0;
    (jpost = nconst -.5*(nNew+d+k+2)*log(s2.eps)-(1/(2*s2.eps))*(QF.y + QF.B + QF.fut))  
   return(as.numeric(jpost))  
  }

   ii=2; hh=gg=1;  
  #----------------------Gibbs Sampler with Monte Carlo integration step----------------------------------------------------#
  for(ii in 2:MCsim)            #Gibbs cycles
  {
    cat("Draw nr.",ii,"\n");
    
    ## Break date: Draw from posterior m|\beta,\sigma^2, Wk,y (see p.22 own derivations):
    ##-------------------------------------------------------------------------------------##
    if(sbreak){     #should break point be sampled from the multinomial posterior?
  
      (kBeg = k[1] + tau + p) ; (kEnd = k[3]-1); (prob = numeric(kEnd-kBeg+1)); (j = kBeg);
      while(j <= kEnd)
      {
        (X.til = strbr.compS(nseason=nseason,season=season,parX=X,lags=p,y=y, Tb=j,isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau))
        (beta = ginv(t(X.til)%*%X.til)%*%t(X.til)%*%lhs)        #update \beta vector via OLS using the updated design matrix X.til
        (err = lhs[kBeg:kEnd] - X.til[kBeg:kEnd,]%*%beta)
        (prob[j-p-tau] = -length(prob)*log(sqrt(s2.eps[ii-1]))-sum(err^2)/(2*s2.eps[ii-1]))
        j=j+1; print(j)
      }
      prob = exp(prob-max(prob)) ; prob = prob/sum(prob);         #normalize to class probabilities!
      
     #Write output to file:
   # write(prob,file="C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Seasonal_model\\Output\\bdat_eprobs.out", append=T, ncol=length(prob))
      
      if(length(prob)>0){
        mn = rmultinom(n = 1, size = length(prob), prob = prob)
        postki = which.max(mn)   
        Tb[ii] = postki + k[1] + p + tau - 1      #correct for lagging, to get the true date within the sample
      }
      cat("\nDrawn break date:", Tb[ii], "\n\n")
    }   

    ## Update W.star matrix of lagged W-values with respect to deterministic parts -> break date at t: 
    ##-------------------------------------------------------------------------------------------------##
    (y.aug = ts(c(y,W.post[ii-1,]),fre=S,start=start(y)))     
    (date = if(sbreak){Tb[ii]} else{NULL})
    if(sbreak.fix && !sbreak){date = floor(runif(1,min=k[1] + tau + p, max=k[3]-1)); print(date)}
    
    out2 = bayes.par12(data=y.aug, class.est=F, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau)
    (W.star = out2@W); (fc.rows = (nrow(W.star)-n.ahead+1):nrow(W.star))                 #row indices of the future observations  
    (X.til = W.star[-fc.rows,])                 #everything except rows of future observations; autoregressive lags of y and deterministic parts
    (X = X.til[,1:(p*S)]);                      #only autoregressive lags of y
    (W.star = matrix(W.star[fc.rows,],nrow=n.ahead,ncol=ncol(W.star),dimnames=list(paste("k=",1:n.ahead,sep=""),NULL)))
    
    ## Draw from posterior \beta|\sigma^2,W,m,y (see p.19 own derivations):
    ##--------------------------------------------------------------------------##
    (R = t(X.til)%*%X.til + t(as.matrix(W.star))%*%as.matrix(W.star) + inv.V)
      
    Rinv = qr.solve(R,tol = 1e-10)                #use QR- decomposition to compute the inverse
    #Rinv = ginv(R)                               #use Moore Penrose inverse, because with increasing sample size there are too many zeros (column containing only zeros) due to the structural break in the seasonal means (=dummies)
    
    (b.mean = Rinv%*%(t(X.til)%*%lhs + inv.V%*%b0 + t(as.matrix(W.star))%*%as.matrix(W.post[ii-1,])))
    (b.sigma = s2.eps[ii-1]*Rinv) 
    (U = chol(b.sigma,pivot=F))                           #upper triangular matrix from Cholesky decomp.

    ## Restrict posterior draws to stationary parameter region:
    ##-----------------------------------------------------------##
    if(restr2stat)
   { 
   repeat{
    (b.new = b.mean + t(U)%*%as.matrix(rnorm(length(b.mean),0,1)))
     b[ii,] = b.new                                             #draw new beta vector from full cond. posterior
     
     mat = mats01(b[ii,1:(p*S)],p=p,S=S)           #compute \Omega_{0} and \Omega_{1} matrices for the VQ(1) representation, see Franses(2006), p.32
     prod.unit[ii] = abs(-(mat$det01 - 1))          #nonlinear restriction from characteristic equation (for z=1) -> periodic unit root, Franses(2006, p.35) 
      
        #Does determinant not lie in region in the neighborhood of det==0 ? Otherwise, y_t instationary, because then z=1 is solution.
    if(ifelse(restr2stat,yes = !(abs(mat$det01)<.05), no=TRUE)){break} else{
      cat('\nPosterior draws of beta discarded.\n')       #draw again a beta vector from the posterior, until stationarity condition is fullfilled
      }
     }
    } else{(b[ii,] = b.mean + t(U)%*%as.matrix(rnorm(length(b.mean),0,1)))
     prod.unit[ii] = abs(prod(b[ii,1:(p*S)]))        #only equivalent to abs(-(mat$det01 - 1)) if PAR(1) process, use here nevertheless, because this output is not really needed in the subsequent  
   }
    
    ## Draw from posterior \sigma^2|\beta,W,m,y (see p.22 own derivations):
    ##----------------------------------------------------------------------##
    (a.star = n - p + n.ahead + length(b[ii,] + a))
    (b.star = t(as.matrix(lhs)-X.til%*%b[ii,])%*%(as.matrix(lhs)-X.til%*%b[ii,]) + t(as.matrix(W.post[ii-1,])-W.star%*%b[ii,])%*%(as.matrix(W.post[ii-1,])-W.star%*%b[ii,]) + t(b[ii,]-b0)%*%inv.V%*%(b[ii,]-b0) + lambda)
    (s2.eps[ii] = 1/rgamma(1,shape = a.star, rate = b.star))            #Draw from Inverse Gamma Distr. (s2e)
    
    ## Draw from posterior Wk|\beta,\sigma^2, W_k1,m,y (see p.22 own derivations):
    ##---------------------------------------------------------------------------##
     step=kk=lags=1; 
    (W.mean = t(W.star[step,])%*%b[ii,])                       #posterior predictive mean       
    (W.post[ii,step] = rnorm(1,mean=W.mean,sd=sqrt(s2.eps[ii])))   # One step ahead forecast -> kk=1
    (cycle.orig = as.numeric(cycle(y)))                       # phases associated with the original series
    (cycle.aug = as.numeric(cycle(y.aug)))                    # phases associated with the future observation augmented series

     ## 1-step ahead forecast -> Update step:
    if(n.ahead>1){
      repeat{
        if(kk<=(n.ahead-1) && lags<=p){
          (W.star[(kk+1), cycle.aug[length(cycle.orig) + (kk+1)] +(lags-1)*S] = W.post[ii,step]);  #alles um (p-1)*S weiterverschieben
          #(W.star[(kk+1),end(y)[2]+(kk+1)+(lags-1)*S] = W.post[ii,step])    
          kk=kk+1; lags=lags+1;
        } else{break}
      }
      ##Monte Carlo integration -> univariate posterior predictive distributions------------------------------#

     marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + dnorm(W_supp[hh],mean=t(W.star[step,])%*%b[ii,], sd=sqrt(s2.eps[ii])))/gg      #Bauwens, p.86 (3.53)
      #marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + hh)/gg   #step/gg 
      #---------------------------------------#
    }
    #------k>1 step ahead----------------------#
     step=2;
    while(step<=n.ahead){
      kk=step;
      (W.mean = t(W.star[step,])%*%b[ii,])                          
      (W.post[ii,step] = rnorm(1,mean=W.mean,sd=sqrt(s2.eps[ii])))       #draw from the full cond. posterior of W(step)
      
       lags=1;      
      if(n.ahead>1){
        repeat{
          if(kk<=(n.ahead-1) && lags<=p){
            (W.star[(kk+1), cycle.aug[length(cycle.orig) + (kk+1)] +(lags-1)*S] = W.post[ii,step])
            kk=kk+1; lags=lags+1;
          } else{break}
        }
      }
      ##Monte Carlo integration -> univariate posterior predictive distributions-------------------#   
     marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + dnorm(W_supp[hh],mean=t(W.star[step,])%*%b[ii,], sd=sqrt(s2.eps[ii])))/gg      
      #--------------------------------------------------------------------------------------------#
      step=step+1; 
    }   #end while loop
    if(gg < (MCsim/le.W)){gg=gg+1} else{gg=1; if(hh<le.W){hh=hh+1}}    #nächsten punkt in grid wählen, wieder von 1 ab zählen
     
    mll = loglike(X.til=X.til,lhs=lhs,s2.eps=s2.eps[ii])
    (marglog = mll$LL); margloglike = c(margloglike,marglog);        #marginal Loglikelihood       		     
    
    if(marglog > trigger){max.date = date; trigger = marglog}
    
       # Evaluate the joint log-posterior at the current MCMC draws, for DIC, see below:
    jlogp.values = c(jlogp.values, -2*jointpost(X.til=X.til,lhs=lhs,s2.eps=s2.eps[ii],y.til=W.post[ii,], W.til=W.star, B=b[ii,], B0=b0, inv.V=inv.V))

  } #end Gibbs loop
  
      # Compute the DIC:
      #-----------------#
     #(jposts= exp(jlogp.values - max(jlogp.values)))
    (D.overline = mean(jlogp.values))
    (D.postmeans = -2*jointpost(X.til=X.til,lhs=lhs,s2.eps=mean(s2.eps),y.til=colMeans(W.post), W.til=W.star, B=colMeans(b), B0=b0, inv.V=inv.V))
    (p.D = D.overline - D.postmeans) ; (DIC = p.D + D.overline)	
   
     # Compute the IC, see Ando (2012):
     #----------------------------------#
     (IC = DIC + p.D)
   
  #(marglike = exp(margloglike - max(margloglike)))
  #(MLH = 1/mean(1/marglike))                       # Harmonic Mean Estimator of the conditional posterior (given p=p^{\star})
  
  ## Compute marginal predictive densities and corresponding predictive means, variances and HPDs:
  ##------------------------------------------------------------------------------------------------##
   require(Bolstad); erw=vari=hpd95=NULL; 
  for(ii in 1:nrow(marg.Ws)){
    (con = sintegral(x = W_supp, fx=marg.Ws[ii,])$value)
    (erw = c(erw,sintegral(x = W_supp, fx=W_supp*marg.Ws[ii,]/con)$value))
    (vari = c(vari,sintegral(x = W_supp, fx=((W_supp-erw[ii])^2)*marg.Ws[ii,]/con)$value))
    marg.Ws[ii,] = marg.Ws[ii,]/con
    (cdf = sintegral(W_supp,marg.Ws[ii,])$cdf)
    (q0025=cdf$x[min(which(cdf$y>=0.025))])        #compute 0.025-quantile
    (q0975=cdf$x[min(which(cdf$y>=0.975))])        #compute 0.975-quantile
    (hpd95 = rbind(hpd95, c(q0025,q0975)))          #HPD regions
  }
  if(in.samp.fc){(BIAS = erw - y.true.T)} else{BIAS=NULL}

  (s_t = ts(1 + ((seq(length(y) + n.ahead + start(y)[2]-1) - 1)%%S)[-c(1:(start(y)[2]-1))],start=start(y),fre=S))    
   #R2 = 1 - mean(s2.eps)/var(y)   
 return(invisible(list(control=c(M,burnin,n.ahead),OLS=OLS,s2.eps=s2.eps,Tb=Tb,prod.unit=prod.unit,
             b=b,W.post=W.post,ytrunc=y,y.true.T=y.true.T,n.ahead=n.ahead, 
             jlogp.values=jlogp.values, marg.Ws=marg.Ws,W_supp=W_supp,
             hpd95=ts(hpd95,end=end(s_t),fre=S), st=s_t,DIC=DIC,D.overline=D.overline, p.D.=p.D,IC=IC,
             pred.mean=list(Mean=erw,Variance=vari,y.true.T=y.true.T,BIAS=BIAS))))
}  #end function
#----------------------------------------------End-----------------------------------------------------------------------------------------------------#
dput(gibbs12,"gibbs12.R")

(out12 = gibbs12(y=y, p=1, n.ahead=1*12, M=5500, burnin=500, in.samp.fc=F, sbreak=T, sbreak.fix=0, sdeter="const", deter=NULL, isbreak.sea=c(1,0), isbreak.nsea=c(0,0), s2.hyp=c(4.01,4.01),le.W=60,low_up=c(2,2), tau=11))    #existence of the first four moments for \nu>4 large \lambda induces large mean and variance

length(orig)
length(out12$ytrunc)

#(cand1 = data.frame(k=1:12, Wk = out12$pred.mean$Mean, Wk.true = out12$pred.mean$y.true.T, bias=out12$pred.mean$BIAS))
(m1 = median(abs(out12$pred.mean$BIAS)))

(cand2 = data.frame(k=1:12, Wk = out12$pred.mean$Mean, Wk.true = out12$pred.mean$y.true.T, bias=out12$pred.mean$BIAS))
(m2 = median(abs(out12$pred.mean$BIAS)))

sink(paste(path,"comparison.txt",sep=""))

print(cand1);print(cand2);print(mean(cand2[,4]));print(mean(cand1[,4]))

sink()


plot(cumsum(out12$pred$BIAS^2),type="l",ylab="",xlab="Forecast horizon (in months)")



## Plot marginal predictive densities:
##----------------------------------------##

W_supp = out12$W_supp ; pred.dens = t(out12$marg.Ws) ; y.true.T = out12$y.true.T ; hpd95 = out12$hpd95

par(ask=T,las=1,font.main=11,mfrow=c(2,2)) ; 
ii=nr=1;
nr=6

#pdf(paste(path,"margPred_s12_",nr,".pdf",sep=""))  
#par(ask=F,las=1,font.main=11,mfrow=c(2,2)) ;

while(ii<=nr*4)    #nrow(out12$marg.Ws)
{
 plot(W_supp,pred.dens[,ii],type="l",xlim=c(hpd95[ii,1]-10,hpd95[ii,2]+10),ylab="density",xlab=substitute(y[t+k],list(k=ii)),sub=expression(S==12))   #substitute(list(eta[j],eta[j+K]) == group("(",list(x,y),")"),list(x=(ii-5), y=ii))
 abline(v=y.true.T[ii],lty=3,col="red4");
 abline(v=hpd95[ii,],lty=2,col="gray")
 legend("topleft",legend=substitute(y[t+k]^true,list(k=ii)),bty="n",lwd=2,col="red4",lty=3)

ii=ii+1;
}
#dev.off()
#(nr=nr+1);             #number of plots


#=========================================================================================================================================#
# Integrate out all modelparameters via importance sampling
# use the (conditional) likelihood values (under each model) to compute model probabilities
# and 'integrate' out the number of breaks and the number of lags via Bayesian model averaging (Monte Carlo integration over model space)
#=========================================================================================================================================#

#options(show.error.messages = FALSE)

margpostTb12 = dget("margpostTb12.R") ; require(coda)

n.ahead=12*2                         #number of full cycles to forecast
pmax=5                               #maximum number of autoregressive lags
sdeter="const" ;isbreak.sea=c(1,0);
expand=1; M=10000; burnin=400; le.W=40;
pred=mLL=fc.BMA=aic=bic=hq=fpe=Marg=bma.forec=NULL; m=pred.dens=0; 
#pred.dens=matrix(0,nrow=M+burnin,ncol=n.ahead,dimnames=list(NULL,paste(1:n.ahead,"step",sep="")))

while(m<=1){  
  for(pp in 1:pmax)
  {
    print(m); print(pp);     
    (out = gibbs12(y=y, p=pp, n.ahead=n.ahead, M=M, sbreak=F, sbreak.fix=m, s2.hyp=c(2.001,4.001),sdeter=sdeter, isbreak.sea=isbreak.sea, burnin=burnin,le.W=le.W,low_up=c(2,2)))
    (caseTb = margpostTb12(data=y,p=pp, sbreak=m, normc=F, Tb.fix=NULL, deter=NULL,sdeter=sdeter,isbreak.nsea=c(1,0),isbreak.sea=isbreak.sea,coding=c("mean","effect")))
    aic = c(aic,caseTb$info$bic); bic = c(bic,caseTb$info$bic);               #info criteria
    hq = c(hq,caseTb$info$hq); fpe = c(fpe,caseTb$info$fpe);
    Marg = c(Marg,caseTb$marg)                                              #marginal likelihood (analytical) estimates   
    pred.dens = pred.dens + (expand * caseTb$marg * t(out$marg.Ws)[])       #compute weighted sum of marginal predictive density values under each model (finally divide by marginal likelihood value)    
    (bma.forec = rbind(bma.forec, out$pred$Mean))                           #save (marginal) predictive means from Monte Carlo integration for weighting with model posterior probabilities -> BMA point forecasts   
    W.post.new = mcmc(out$W.post[-c(1:burnin),])                                #compare Bayesian in-sample forecast with true value
    (bayes.fc = summary(W.post.new))
    if(n.ahead>1){(fc=bayes.fc$stat[,"Mean"])}else{(fc=bayes.fc$stat["Mean"])}     #mean
    fc.BMA=rbind(fc.BMA,fc)                                                   #läuft auch, Ergebnisse so gut wie identisch                                                            #multiply likelihood values with some large number to avoid numerical problems due to too small numbers which are multiplied by the predictive densities
  }
  m=m+1;
}

#print(.Last.value); 
#options(show.error.messages = TRUE)
#===============================================================================================================================#

pred.dens[] = pred.dens[]/(sum(Marg)*expand); 
W_supp=out$W_supp; y.true.T = out$y.true.T          #used for plots, see above
dim(pred.dens)

## Finally normalize the model averaged marginal predictive densities:
##---------------------------------------------------------------------##
  require(Bolstad); erw.BMA=hpd95=NULL; 
 for(ii in 1:ncol(pred.dens)){
   (con = sintegral(x = W_supp, fx=pred.dens[,ii])$value)
   pred.dens[,ii] = pred.dens[,ii]/con
   (erw.BMA = c(erw.BMA,sintegral(x = W_supp, fx=W_supp*pred.dens[,ii])$value))
   (cdf = sintegral(W_supp,pred.dens[,ii])$cdf)
   (q0025=cdf$x[min(which(cdf$y>=0.025))])        #compute 0.025-quantile
   (q0975=cdf$x[min(which(cdf$y>=0.975))])        #compute 0.975-quantile
   (hpd95 = rbind(hpd95, c(q0025,q0975)))          #HPD regions
 }
erw.BMA;hpd95


## Save results to file:
##-----------------------##
#unlink(paste(path,"forecast_output_s12.out",sep=""))
sink(paste(path,"forecast_output_s12.out",sep=""))

model.postwei = round(Marg/sum(Marg),10);                                                              #with prior model probabilities implictly set to zero
data.frame(AIC=round(aic,6),BIC=round(bic,6),HQ=round(hq,6),FPE=round(fpe,6),ModProb=model.postwei)    #compare with frequentist info criteria

BMA_forecasts = t(model.postwei)%*%bma.forec;
colnames(BMA_forecasts)=colnames(pred.dens);
data.frame(BMA=as.vector(BMA_forecasts),ytrue=as.vector(out$y.true.T),hpd95.BMA=hpd95)       #compare forecasts with true values

#Wp = summary(W.post.new); Wp$stat[,"Mean"];

sink()


## Plot series together with multiple forecasts (and possibly true values):
##==========================================================================##

#pdf(paste(path,"forecasts_vs_true.pdf",sep=""))  

par(font.main=11,cex.main=1.1,mfrow=c(2,1),las=1)
plot(ts(c(out$ytrunc,BMA_forecasts),fre=12,start=start(y)),col="black",type="l",ylab=expression(y[t]),lty=2);
lines(out$ytrunc);
legend("topleft",legend=c(expression(y[t]),expression(tilde(y)[T+k])),col=c("black","black"),bty="n",lwd=3,,lty=c(1,2))
#title("Simulated PAR(p) series\n with k step prediction")
plot(ts(y,fre=12,start=start(y)),col="black",type="l",ylab=expression(y[t]),lty=2,ylim=c(min(y),max(y)));
lines(out$ytrunc,lty=1)
legend("topleft",legend=c(expression(y[t]),expression(y[T+k]^true)),col=c("black","black"),lty=c(1,2),bty="n",lwd=3)

dev.off()

#boa.menu()


##----------------------------------------------------##
## Compare with classical estimates:
##----------------------------------------------------##
require(coda)

b=out12$b ; prod.unit=out12$prod.unit ; s2.eps=out12$s2.eps ;W.post=out12$W.post ; 
(burnin=out12$control[2]);(OLS=out12$OLS); (n.ahead=out12$control[3])

b.new = mcmc(b[-c(1:burnin),]);                              #keep only draws after burnin phase
(bayes = summary(b.new)); attributes(bayes)
hpd=HPDinterval(b.new); #rownames(hpd)=colnames(X.til);hpd

prod.unit.new = mcmc(prod.unit[-c(1:burnin)])
summary(prod.unit.new)

#nvar(b.new); mcpar(b.new)

s2.eps.new = mcmc(s2.eps[-c(1:burnin)])             
summary(s2.eps.new) ; (s2.class = OLS$sig^2)

effectiveSize(s2.eps.new)
#(var.mean = spectrum0(s2.eps.new)$spec/niter(s2.eps.new))       #Estimated variance of mean(s2.eps.new)      

bay = as.matrix(bayes$stat); #rownames(bay)=colnames(X.til);
data.frame(bayes=bay[,1],classic=coef(OLS)[,1])                #compare both

#summary(mcmc(Tb[-c(1:burnin)]))

#W.post.new = mcmc(W.post) 
W.post.new = mcmc(W.post[-c(1:burnin),])        #compare Bayesian in-sample forecast with true value
(bayes.fc = summary(W.post.new))

out12$y.true.T
if(n.ahead>1){(fc=bayes.fc$quant[,"50%"])}else{(fc=bayes.fc$quant["50%"])}
hppd = HPDinterval(W.post.new); rownames(hppd)=paste("y_T+",1:n.ahead,":",sep="");
hppd            #Highest posterior predictive density regions for each point forecast


par(font.main=11,cex.main=1.2)
#plot(ts(c(y,fc),fre=12,start=start(y)),col="red3",type="l",ylab=expression(y[t]),lty=2);
plot(c(y,fc),col="red3",type="l",ylab=expression(y[t]),lty=2, las=1)
lines(c(y))
abline(v=length(y),lty=2)
legend("topleft",legend=c(expression(y[t]),expression(tilde(y)[T+k])),col=c("black","red3"),bty="n",lwd=3,lty=c(1,2))
title("Simulated PAR(p) process with k-step prediction")


##---------------------------------##
## Some convergence diagnostics:
##---------------------------------##

## Regression coefficients (i.e. AR and deterministic):
##------------------------------------------------------##
graphics.off()

par(ask=T);
for(ii in 1:ncol(b.new))
{
  plot(b.new[,ii],main="")
  cumuplot(b.new[,ii],main=paste("Chain of '",colnames(b.new)[ii],"'",sep=""))    
  #hist(b.new[,ii],prob=TRUE,nclass=120,col="wheat2",xlab="",main="")
  #densityplot(b.new[,ii], start = 10)
  acf(b.new[,ii],main=paste("Chain of '",colnames(b.new)[ii],"'",sep="")) 
  print(autocorr.diag(b.new[,ii]))                           
  #gelman.plot(mcmc.list(b.new[,1],b.new[,2]))
}

## Residual autocovariances:
##---------------------------##
graphics.off() ;par(ask=T);

plot(s2.eps.new,main="");cumuplot(s2.eps.new)    
acf(s2.eps.new);print(autocorr.diag(s2.eps.new))     		
#traceplot(s2.eps.new)

## Single Break date:
##--------------------##
dev.off();par(ask=F,las=1,cex.main=1);
hist(Tb,prob=TRUE,col="wheat2",xlab=expression(T[B]),main="",xlim=c(1,n))
title("Posterior probability mass\n function of the break point")

#plot(Tb,main="");
#cumuplot(Tb)    
#acf(Tb);print(autocorr.diag(Tb))

## Product of PAR coefficients for unit root testing:
##---------------------------------------------------##
dev.off()
par(ask=T);
plot(prod.unit.new,main="");
cumuplot(prod.unit.new)    
acf(prod.unit.new);print(autocorr.diag(prod.unit.new))       	

## Model prediction:
##-------------------##
graphics.off()
#par(ask=T,cex.main=1,font.main=12);
for(ii in 1:n.ahead)
{
  if(n.ahead>1){pred=W.post.new[,ii]}else{pred=W.post.new}  
  plot(pred,main="",sub=paste("Posterior prediction k=",ii,sep=""));
  cumuplot(pred,main=paste("Posterior prediction k=",ii,sep=""))    
  acf(pred,main=paste("Posterior prediction k=",ii,sep=""));
  print(autocorr.diag(pred))
}
graphics.off()
#------------------------------------------------------------------------------------------#




















