
##==================================================================##
## Estimate and predict quarterly PAR(p) process via Gibbs Sampling:
##==================================================================##
## In this version, first break dates are drawn and then all the other parameters (was just thought as a test)

rm(list=ls())

## Simulate some quarterly PAR(1)/PMA(1) or PAR(2) data with(out) structural break:
##----------------------------------------------------------------------------------##
source("par_simulation4.R")    #Home

set.seed(123)
TT=500 ;                         #the higher n, the more accurate the forecasts!

(y = par1(N=TT, start=c(1976,2), mu=c(1,1,1,1), beta=0, phi=c(.87,1.2,.95,.85))$y)               #PAR(1) without breaks

(y = par2(n=TT, mu=c(0,0,0,0), beta=0, phi=c(.2,.2,.2,.2,.5,.5,.5,.5))$y)              

(y = par3(n=TT, mu=c(1,1,1,1), beta=0, phi=c(.1,.1,.1,.1,.4,.4,.4,.4,.4,.3,.3,.3), sig=c(1,1,1,1)))

(y = par4(n=TT, mu=c(0,0,0,0), beta=0, phi=c(.1,.1,.1,.4,.1,.1,.1,.4,.2,.2,.2,.2,.4,.4,.4,.4),sig=c(1,1,1,1)))

(y = pma1(N=TT,mu=c(1,1,1,1), beta=0, theta=c(.7,.7,.75,1), sig=c(1,1,1,1)))               #PMA(1)

Tb.true = 55;
(y = broken.par1(n=TT, Tb=Tb.true,mu=c(1,.1,1,.65,1,.65,1,.65), alpha=rep(0,8), phi=c(0.7,1.2,0.75,.55),sig=rep(0.3,4)))

(y = broken.par2(n=TT, Tb=Tb.true,mu=c(1,-.1,1,.1,1,.1,1,.1), alpha=rep(0,8), phi=c(.1,.1,.1,.1,.5,.5,.5,.5),sig=c(.5,.5,.5,.5)))


## Plot series, (non)periodic (P)ACF, periodogram and boxplot of seasonal subseries:
##----------------------------------------------------------------------------------##

par(mfrow=c(2,2));
plot(y);
acf(as.vector(y),main="");          #nonperiodic (homogenous) correlograms
pacf(as.vector(y),main="");
spectrum(as.vector(y),main="")

graphics.off()

par(mfrow=c(2,2),las=1);
plot(y);
pa = peacf(y,plot=T);                                 #periodic (heterogenous) correlograms
ppacf = pepacf(y,lag.max=7, plot=T);
attributes(ppacf)
ppacf$bic
peboxplot(y,xlab="period")     
ans = pear(y, ic="bic")           #fit possibly individual model orders       

graphics.off()



pa$periodicity.test
ans$portmanteau.test
ans$model.orders
ans$phi

sink()
##===============================================================================================##

gibbs4 = function(y,p,n.ahead,M,burnin, restr2stat=F, sbreak=T, sbreak.fix=F, s2.hyp=c(nu,lambda),in.samp.fc=T, le.W=100,low_up=c(4,4),deter=NULL, sdeter="const",isbreak.nsea=c(1,0),isbreak.sea=c(1,0),tau=3)  
{                     
   library(MASS) 
  bayes.par=dget("bayes.par.R"); strbr.comp = dget("strbr.comp.R"); mats01 = dget('mats01.R'); 

  nseason=deter; season=sdeter; coding=c("mean","effect"); S=4 ; 

  ## Save last n-ahead time observations for in-sample forecast control:
  ##---------------------------------------------------------------------##
  if(frequency(y)!=S){stop('Time series has wrong seasonal frequency!\n')}
  if(in.samp.fc){
    (y.true.T = ts(y[c((length(y)-n.ahead+1):length(y))],end=end(y),fre=S));
    (y = ts(y[-c((length(y)-n.ahead+1):length(y))],start=start(y),fre=S))       #actually used sample
  } else{y.true.T=NULL}          
  (n = length(y)); (n.new = n-p)                       
  (z = rep(mean(y),times=n.ahead))
  (y.aug = ts(c(y,z),fre=S,start=start(y)))
  (k = c(1,floor((n-1)*.5),n.new))                 #Break interval for the single break date
  #(k = c(1,Tb.true,n.new))                      #Attention: If k does not start from 1, you must modify the vector indices in prob[j..], see step for break date draws
  (date = if(sbreak){k[2]} else{NULL})
  if(sbreak.fix && !sbreak){
    date = floor(runif(1,min=k[1] + tau + p, max=k[3]-1));
   message('Draw break date from uniform instrumental density instead of posterior.\n')}    #draw some date from a uniform density
  if(sbreak.fix && sbreak){
    message('Option is not supported. Set "sbreak <- FALSE" if uniform density should be used.\nMultinomial posterior density is used instead.\n')}
  
  #M=10000; burnin=200; 
  MCsim = M+burnin;          #Monte Carlo runs and burnin draws
  
  (out1 = bayes.par(data=y, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau))
  (lhs = out1@lhs);(X.til = out1@W)          
  (X = X.til[,1:(p*S)]);                    #only autoregressive terms 
  out2 = bayes.par(data=y.aug, class.est=T, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau)
  (W.star = out2@W); (W.star = ts(matrix(W.star[(nrow(W.star)-n.ahead+1):nrow(W.star),],nrow=n.ahead),fre=S,end=end(W.star)))
  
  ##-------------------------## 
  ## Initialize some objects:
  ##-------------------------##
  s2.eps = Tb = prod.unit = numeric(MCsim)
  if(sbreak){Tb[1]=date}else{Tb=NULL}
  b = matrix(,ncol=ncol(X.til),nrow=MCsim)
  W.post = matrix(,ncol=n.ahead,nrow=MCsim)
  (W.post[1,]=z)  
  (s2.eps[1] = out2@sigma2)          #initialize s2 with classical estimates
  (b[1,] = out2@est.classic)
  (OLS = out1@est.classic)                      #frequentist estimates
  (prod.unit[1] = abs(prod(b[1,1:(p*S)])))           #product of AR coefficients to test for a (monthly) periodic unit root
  margloglike=NULL ; 
  (W_supp = seq(min(y)-low_up[1]*sd(y),max(y)+low_up[2]*sd(y),length.out=le.W))
  print(W_supp)
   #if(is.null(le.W)){le.W = length(W_supp)} 
  marg.Ws = matrix(0,nrow=n.ahead,ncol=length(W_supp),dimnames=list(paste("W",1:n.ahead,sep=""),paste(1:length(W_supp),sep="")))  
  if((MCsim)%%le.W){stop(" 'MCsim' must be a multiple of length(W_supp).\n")
   } else{cat("There are",MCsim/le.W,"draws per point in 'W_supp' used.\n")}

  ## Define all hyperparameters for \beta normal prior (as on p.9):
  ##---------------------------------------------------------------##
  c1 = 100 ; c2 = 100                                        #hyperparameters prior covariance matrix of \beta
  (b0 = c(rep(1,times=p*S),rep(0,ncol(X.til)-p*S)))          #prior mean
  V = diag(c(rep(c1,times=p*S),rep(c2,ncol(X.til)-p*S)))     #see p.19 and p.10 for definitions
  a = s2.hyp[1]; lambda = s2.hyp[2]                          #parameters IG prior for \sigma^2 
                   #existence of the first four moments for \nu>4
                   #large \lambda induces large mean and variance

  ## Function to compute the joint loglikelihood function:
  ##------------------------------------------------------##
  loglike = function(X.til,lhs,s2.eps){
    (beta = ginv(t(X.til)%*%X.til)%*%t(X.til)%*%lhs);     
    (err = lhs - X.til%*%beta); nNew=length(lhs);
    -(nNew/2)*log(2*pi)-(nNew/2)*log(s2.eps)-(1/(2*s2.eps))*t(err)%*%err    
   }
  
   ii=2; hh=gg=1;  
  #--------------Gibbs Sampler with Monte Carlo integration step----------------------------------------------------#

  for(ii in 2:MCsim)
  {
    cat("Draw nr.",ii,"\n")
    
    ## Draw from posterior m|\beta,\sigma^2, Wk,y (see p.22 own derivations):
    ##------------------------------------------------------------------------##
    if(sbreak){
      
      assign("lags",p, envir=environment(strbr.comp))       #assignment is only for S=4 function necessary
      (kBeg = k[1] + tau + p) ; (kEnd = k[3]-1)
      (prob = numeric(kEnd-kBeg+1)); (j = kBeg);
      while(j <= kEnd)
      {
        (X.til = strbr.comp(nseason=nseason,season=season,parX=X,y=y, Tb=j, tau=tau,isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea))  
        (beta = ginv(t(X.til)%*%X.til)%*%t(X.til)%*%lhs)     #update \beta vector via OLS using the updated design matrix X.til
        (err = lhs[kBeg:kEnd] - X.til[kBeg:kEnd,]%*%beta)
        (prob[j-p-tau] = -length(prob)*log(sqrt(s2.eps[ii-1]))-sum(err^2)/(2*s2.eps[ii-1]))
        j=j+1; print(j)
      }
      prob = exp(prob-max(prob)) ; prob = prob/sum(prob);
      
     #Write output to file:
   # write(prob,file="C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Seasonal_model\\Output\\bdat_eprobs.out", append=T, ncol=length(prob))
      
      if(length(prob)>0){
        mn = rmultinom(n = 1, size = length(prob), prob = prob)
        postki = which.max(mn)   
        Tb[ii] = postki + k[1] + p + tau - 1      #correct for lagging, to get the true date within the sample
      }
      cat("\nDrawn break date:", Tb[ii], "\n\n")
    }   
    ## Update W.star matrix of lagged W values with respect to deterministic parts -> break date at t: 
    ##------------------------------------------------------------##
    (y.aug = ts(c(y,W.post[ii-1,]),fre=S,start=start(y)))     
    (date = if(sbreak){Tb[ii]} else{NULL})
    if(sbreak.fix && !sbreak){date = floor(runif(1,min=k[1] + tau + p, max=k[3]-1)); print(date)}
    
    out2 = bayes.par(data=y.aug, class.est=FALSE, p=p, Tb=date, deter=deter, sdeter=sdeter, isbreak.nsea=isbreak.nsea,isbreak.sea=isbreak.sea,tau=tau)
    (W.star = out2@W); (fc.rows = (nrow(W.star)-n.ahead+1):nrow(W.star))   #row indices of the future observations  
    (X.til = W.star[-fc.rows,])                #everything except rows of future observations; autoregressive lags of y and deterministic parts
    (X = X.til[,1:(p*S)]);                      #only autoregressive lags of y
    (W.star = matrix(W.star[fc.rows,],nrow=n.ahead,ncol=ncol(W.star),dimnames=list(paste("k=",1:n.ahead,sep=""),NULL)))
    
    ## Draw from posterior \beta | \sigma^2, W, m, y (see p.19 own derivations):
    ##--------------------------------------------------------------------------##
    (R = t(X.til)%*%X.til + t(as.matrix(W.star))%*%as.matrix(W.star) + solve(V))

      #Use Moore Penrose inverse, because with increasing sample size there are too many zeros (column containing only zeros) due to the structural break in the seasonal means (=dummies)
    
    (b.mean = ginv(R)%*%(t(X.til)%*%lhs + solve(V)%*%b0 + t(as.matrix(W.star))%*%as.matrix(W.post[ii-1,])))   #posterior mean    
    (b.sigma = s2.eps[ii-1]*ginv(R)) 
    (U = chol(b.sigma,pivot=F))                           #upper triangular matrix from Cholesky decomp.
    
     ## Restrict posterior draws to stationary parameter region:
     ##-----------------------------------------------------------##
    if(restr2stat)
   { 
   repeat{
    (b.new = b.mean + t(U)%*%as.matrix(rnorm(length(b.mean),0,1)))
     b[ii,] = b.new                                             #draw new beta vector from full cond. posterior
          #own matrix function:
     mat = mats01( b[ii,1:(p*S)], p=p,S=S)         #compute \Omega_{0} and \Omega_{1} matrices for the VQ(1) representation, see Franses(2006), p.32
     prod.unit[ii] = abs(mat$det01 - 1)          #nonlinear restriction from characteristic equation (for z=1) -> periodic unit root, Franses(2006, p.35) 

        #Does determ. not lie in region in the neighborhood of det==0 ? Otherwise, y_t instationary, because then z=1 is solution.
    if(ifelse(restr2stat,yes = !(abs(mat$det01)<.05), no=TRUE)){break} else{ 
      cat('\nPosterior draws of beta discarded.\n')       #draw again a beta vector from the posterior, until stationarity condition is fullfilled
      } 
     }   #end loop
    } else{ (b[ii,] = b.mean + t(U)%*%as.matrix(rnorm(length(b.mean),0,1)));
     prod.unit[ii] = abs(prod(b[ii,1:(p*S)]))            #only equivalent to abs(-(mat$det01 - 1)) if PAR(1) process, use here nevertheless, because this output is not really needed in the subsequent  
   }
 
    ## Draw from posterior \sigma^2|\beta,W,m,y (see p.22 own derivations):
    ##----------------------------------------------------------------------##
    (a.star = n - p + n.ahead + length(b[ii,] + a))
    (b.star = t(as.matrix(lhs)-X.til%*%b[ii,])%*%(as.matrix(lhs)-X.til%*%b[ii,]) + t(as.matrix(W.post[ii-1,])-W.star%*%b[ii,])%*%(as.matrix(W.post[ii-1,])-W.star%*%b[ii,]) + t(b[ii,]-b0)%*%solve(V)%*%(b[ii,]-b0) + lambda)
    (s2.eps[ii] = 1/rgamma(1,shape = a.star, rate = b.star))                        #Draw from Inverse Gamma Distr. (s2e)
    
    ## Draw from posterior Wk|\beta,\sigma^2, W_k1,m,y (see p.22 own derivations):
    ##---------------------------------------------------------------------------##
    step=kk=lags=1; 
    (W.mean = t(W.star[step,])%*%b[ii,])                          
    (W.post[ii,step] = rnorm(1,mean=W.mean,sd=sqrt(s2.eps[ii])))   # One step ahead forecast -> kk=1
    
    (cycle.orig = as.numeric(cycle(y)))                       # phases associated with the original series
    (cycle.aug = as.numeric(cycle(y.aug)))                    # phases associated with the future observation augmented series

     ## 1-step ahead forcast -> Update step:
    if(n.ahead>1){
      repeat{
        if(kk<=(n.ahead-1) && lags<=p){
          (W.star[(kk+1), cycle.aug[length(cycle.orig) + (kk+1)] +(lags-1)*S] = W.post[ii,step]);  #alles um (p-1)*S weiterverschieben
          #(W.star[(kk+1),end(y)[2]+(kk+1)+(lags-1)*S] = W.post[ii,step])    
          kk=kk+1; lags=lags+1;
        } else{break}
      }
      #------------------------------#
     marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + dnorm(W_supp[hh],mean=t(W.star[step,])%*%b[ii,], sd=sqrt(s2.eps[ii])))/gg
     #marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + hh)/gg   #step/gg 
      #---------------------------------------#
    }
    #------k>1 step ahead----------------------#
     step=2;
    while(step<=n.ahead){
      kk=step;
      (W.mean = t(W.star[step,])%*%b[ii,])                          
      (W.post[ii,step] = rnorm(1,mean=W.mean,sd=sqrt(s2.eps[ii])))       
      
      lags=1;      
      if(n.ahead>1){
        repeat{
          if(kk<=(n.ahead-1) && lags<=p){
            (W.star[(kk+1), cycle.aug[length(cycle.orig) + (kk+1)] +(lags-1)*S] = W.post[ii,step])
            kk=kk+1; lags=lags+1;
          } else{break}
        }
      }
      #---------------------------------#     
     marg.Ws[step,hh] = (marg.Ws[step,hh]*(gg-1) + dnorm(W_supp[hh],mean=t(W.star[step,])%*%b[ii,], sd=sqrt(s2.eps[ii])))/gg  
      #---------------------------------#
      step=step+1; 
      
    }   #end while loop
    if(gg < (MCsim/le.W)){gg=gg+1} else{gg=1; if(hh<le.W){hh=hh+1}}    #nächsten punkt in grid wählen, wieder von 1 ab zählen
    
    (marglog = loglike(X.til=X.til,lhs=lhs,s2.eps=s2.eps[ii]))   #marginal Loglikelihood       		    
    margloglike = c(margloglike,marglog);    					
    
  } #end Monte Carlo runs
  
  (marglike = exp(margloglike - max(margloglike)))
  (MLH = 1/mean(1/marglike))                       # Harmonic Mean Estimator of the conditional likelihood (given p=p^{\star})
  
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

 return(list(control=c(M,burnin,n.ahead),OLS=OLS,s2.eps=s2.eps,Tb=Tb,prod.unit=prod.unit,b=b,
	W.post=W.post,ytrunc=y,y.true.T=y.true.T,n.ahead=n.ahead,MLH=MLH,marg.Ws=marg.Ws,
	W_supp=W_supp, hpd95=hpd95,pred.mean=list(Mean=erw,Variance=vari,y.true.T=y.true.T,BIAS=BIAS)))
}  #end function
#-----------------------------------------------End---------------------------------------------------------------------------------------------------------#
dput(gibbs4,"gibbs4.R")

(out4 = gibbs4(y=y,p=1,n.ahead=4*2,M=1000,restr2stat=T, burnin=100, sbreak=F,sbreak.fix=F, in.samp.fc=T, s2.hyp=c(4.001,5),le.W=50,low_up=c(3,3)))    #existence of the first four moments for \nu>4 large \lambda induces large mean and variance


## Plot marginal predictive densities:
##----------------------------------------##

#W_supp = out4$W_supp; pred.dens = t(out4$marg.Ws); y.true.T = out4$y.true.T; hpd95 = out4$hpd95

#par(ask=T,las=1,font.main=11,mfrow=c(2,2)) ; ii=1; nr=3;

 nr=ii=1;     #first plot

par(ask=F,las=1,font.main=11,mfrow=c(2,2)) ; 

while(ii<=nr*4)    #nrow(out4$marg.Ws)
{
  plot(W_supp,pred.dens[,ii],type="l",ylab="density",xlab=substitute(y[t+k],list(k=ii)),sub=expression(S==4))   #substitute(list(eta[j],eta[j+K]) == group("(",list(x,y),")"),list(x=(ii-5), y=ii))
  abline(v=y.true.T[ii],lty=3,col="red4");
  abline(v=hpd95[ii,],lty=2,col="gray")
  legend("topleft",legend=substitute(y[t+k]^true,list(k=ii)),bty="n",lwd=1,col="red4",lty=3)
  
 ii=ii+1;
}
dev.off()

(nr=nr+1);             #switch to next four units block of plots


#=========================================================================================================================================#
# Integrate out all modelparameters via importance sampling
# use the (conditional) likelihood values (under each model) to compute model probabilities
# and 'integrate' out the number of breaks and the number of lags via Bayesian model averaging (Monte Carlo integration over model space)
#=========================================================================================================================================#

margpostTb = dget("margpostTb.R"); gibbs4 = dget("gibbs4.R"); require(coda)

n.ahead=4*3                         #number of full cycles to forecast
pmax=4                              #maximum number of autoregressive lags
sdeter="const"; isbreak.sea=c(1,0)
expand=1; M=4500; burnin=500; le.W=40;
pred=mLL=fc.BMA=aic=bic=hq=fpe=Marg=bma.forec=NULL; m=pred.dens=0; 
#pred.dens=matrix(0,nrow=M+burnin,ncol=n.ahead,dimnames=list(NULL,paste(1:n.ahead,"step",sep="")))
  m=0;pp=1 
while(m<=1){  
  for(pp in 1:pmax)
  {
    print(m); print(pp);     
    (out = gibbs4(y=y, p=pp, n.ahead=n.ahead, M=M, sbreak=F, sbreak.fix=m, le.W=le.W, s2.hyp=c(2.001,4.001),sdeter=sdeter, isbreak.sea=isbreak.sea, burnin=burnin))
     mLL = c(mLL, out$MLH)           #marginal model likelihood via importance estimate   
    (caseTb = margpostTb(data=y,p=pp, sbreak=m, normc=F, Tb.fix=NULL, deter=NULL,sdeter=sdeter,isbreak.nsea=c(1,0),isbreak.sea=isbreak.sea,coding=c("mean","effect")))
    aic = c(aic,caseTb$info$bic); bic = c(bic,caseTb$info$bic);               #info criteria
    hq = c(hq,caseTb$info$hq); fpe = c(fpe,caseTb$info$fpe);
    Marg = c(Marg,caseTb$marg)                                              #marginal likelihood (analytical) estimates
    pred.dens = pred.dens + (expand * caseTb$marg * t(out$marg.Ws)[])       #compute weighted sum of marginal predictive density values under each model (finally divide by marginal likelihood value)
    (bma.forec = rbind(bma.forec, out$pred$Mean))                            #save (marginal) predictive means from Monte Carlo integration for weighting with model posterior probabilities -> BMA point forecasts
    
    W.post.new = mcmc(out$W.post[-c(1:burnin),])                  #compare Bayesian in-sample forecast with true value
    (bayes.fc = summary(W.post.new))
    if(n.ahead>1){(fc=bayes.fc$stat[,"Mean"])} else{(fc=bayes.fc$stat["Mean"])}     #mean
    fc.BMA=rbind(fc.BMA,fc)                                         #läuft auch, Ergebnisse so gut wie identisch                                                               #multiply likelihood values with some large number to avoid numerical problems due to too small numbers which are multiplied by the predictive densities
  }
  m=m+1;
}
#==============================================================================================#
     
pred.dens[] = pp = pred.dens[]/(sum(Marg)*expand);       #marginal posterior predictive density of the multiple forecasts W
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
#
model.postwei = round(Marg/sum(Marg),10);                                                              #with prior model probabilities implictly set to zero
data.frame(AIC=round(aic,6),BIC=round(bic,6),HQ=round(hq,6),FPE=round(fpe,6),ModProb=model.postwei)    #compare with frequentist info criteria

BMA_forecasts = t(model.postwei)%*%bma.forec;
colnames(BMA_forecasts)=colnames(pred.dens);

data.frame(BMA=as.vector(BMA_forecasts),ytrue=as.vector(out$y.true.T),hpd95.BMA=hpd95)       #compare forecasts with true values

#(model.weights = round(mLL/sum(mLL),5));(BMAs=t(t(fc.BMA)%*%model.weights))
#Wp = summary(W.post.new); Wp$stat[,"Mean"];

sink()



## Plot series together with multiple forecasts (and possibly true values):
##==========================================================================##

par(font.main=11,cex.main=1.1,mfrow=c(2,1))
plot(ts(c(out$ytrunc,BMA_forecasts),fre=4,start=start(y)),col="black",type="l",ylab=expression(y[t]),lty=2);
lines(out$ytrunc);
legend("topright",legend=c(expression(y[t]),expression(tilde(y)[T+k])),col=c("black","black"),bty="n",lwd=1,,lty=c(1,2))
#title("Simulated PAR(p) series\n with k step prediction")
plot(ts(y,fre=4,start=start(y)),col="black",type="l",ylab=expression(y[t]),lty=2);
lines(out$ytrunc,lty=1)
legend("topright",legend=c(expression(y[t]),expression(tilde(y)[T+k]^true)),col=c("black","black"),lty=c(1,2),bty="n",lwd=1)

dev.off()

help(setBreakpoint)


##----------------------------------------------------##
## Compare with classical estimates:
##----------------------------------------------------##
require(coda)

b=out4$b ; prod.unit=out4$prod.unit ; s2.eps=out4$s2.eps ;W.post=out4$W.post ; 
(burnin=out4$control[2]);(OLS=out4$OLS); (n.ahead=out4$control[3])

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

W.post.new = mcmc(W.post[-c(1:burnin),])                  #compare Bayesian in-sample forecast with true value
(bayes.fc = summary(W.post.new))

out4$y.true.T
if(n.ahead>1){(fc=bayes.fc$quant[,"50%"])}else{(fc=bayes.fc$quant["50%"])}
hppd = HPDinterval(W.post.new); rownames(hppd)=paste("y_T+",1:n.ahead,":",sep="");
hppd            #Highest posterior predictive density regions for each point forecast


par(font.main=11,cex.main=1.1)
plot(ts(c(y,fc),fre=4,start=start(y)),col="red3",type="l",ylab=expression(y[t]));lines(y)
legend("topright",legend=c(expression(y[t]),expression(tilde(y)[T+k])),col=c("black","red3"),bty="n",lwd=3)
title("Simulated PAR(p) series\n with k step prediction")


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
par(ask=T,cex.main=1,font.main=12);
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




































