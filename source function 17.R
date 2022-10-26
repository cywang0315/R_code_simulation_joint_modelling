m=1000
q=7
internal_knots=c(5,10,15)/2
boundary_knots=c(0,20)/2

####### B-spline basis #########
## my cubic B-spline basis functions ####
mycubicbs=function(s,internal_knots,boundary_knots){
  df=length(internal_knots)+3+1
  r=numeric(df)
  out_knots=c(seq(from=min(boundary_knots),by=-0.1,length=4),seq(from=max(boundary_knots),by=0.1,length=4))
  knots=sort(c(internal_knots,out_knots))
  temp=function(x){
    my0bs=function(x,k){
      a=ifelse(x<knots[k+1]&x>=knots[k],1,0)
      return(a)
    }
    my1bs=function(x,k){
      (x-knots[k])/(knots[k+1]-knots[k])*my0bs(x,k)+
        (knots[k+2]-x)/(knots[k+2]-knots[k+1])*my0bs(x,k+1)
    }
    
    my2bs=function(x,k){
      (x-knots[k])/(knots[k+2]-knots[k])*my1bs(x,k)+
        (knots[k+3]-x)/(knots[k+3]-knots[k+1])*my1bs(x,k+1)
    }
    my3bs=function(x,k){
      (x-knots[k])/(knots[k+3]-knots[k])*my2bs(x,k)+
        (knots[k+4]-x)/(knots[k+4]-knots[k+1])*my2bs(x,k+1)
    }
    for(k in 1:df){
      r[k]=my3bs(x,k)
    }
    return(r)}
  return=list(mat=t(sapply(s,temp)),knots=knots)# a matrix with ncol=df
}

R.2=function(t,knots){
  R.2=matrix(0,ncol=q-2,nrow=q-2)
  for (l in 1:(q-2)){
    if((t>knots[l+2])&(t<=knots[l+3])){
      R.2[l,l]=(t-knots[l+2])^3/(3*(knots[l+3]-knots[l+2])^2)
    }else{
      if((t>knots[l+3])&(t<knots[l+4])){
        R.2[l,l]=(knots[l+4]-knots[l+2])/3-(knots[l+4]-t)^3/(3*(knots[l+4]-knots[l+3])^2)
        R.2[l,l+1]=(-t^3/3+(knots[l+4]+knots[l+3])/2*t^2-knots[l+4]*knots[l+3]*t-(knots[l+3])^3/6+
                      (knots[l+3])^2*knots[l+4]/2)/(knots[l+4]-knots[l+3])^2
        R.2[l+1,l]= R.2[l,l+1]
      }else{
        if(t>=knots[l+4]){
          R.2[l,l]=(knots[l+4]-knots[l+2])/3
          R.2[l,l+1]=(knots[l+4]-knots[l+3])/6
          R.2[l+1,l]=(knots[l+4]-knots[l+3])/6
          
        }
      }
    }
  }
  R.2[1,1]=ifelse(t>=knots[5],(knots[5]-knots[4])/3,(knots[5]-knots[4])/3-(knots[5]-t)^3/(3*(knots[5]-knots[4])^2))
  return(R.2)
}

# true value of longitudinal biomarker (no measurement error)
longit.true=function(t,fix_eff, ran_eff){
  desmat=mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat# design matrix
  desmat%*%(fix_eff+ran_eff)
}







adv=function(x){
  grad=diag(1,length(x))
  grad2=matrix(0,ncol=length(x),nrow = length(x))
  attr(x,"grad")=grad
  attr(x,"grad2")=grad2
  class(x)="adv"
  x
}


quad.adv=function(A,b){
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=c(t(b)%*%A%*%b) # scalar
  attr(d,"grad")=grad.b%*%(2*A%*%b)
  attr(d,"grad2")=grad.b%*%(2*A)
  class(d)="adv"
  d
}


"^.adv"=function(b,alpha){ # the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=b^alpha
  attr(d,"grad")=alpha*b^(alpha-1)*grad.b
  attr(d,"grad2")=alpha*(alpha-1)*b^(alpha-2)*tcrossprod(grad.b)+alpha*b^(alpha-1)*grad2.b
  class(d)="adv"
  d
}


linear_sca.adv=function(A,b,j){ ### derivative of <Ab>j with respect to b, A is matrix independent of b
  grad.b=attr(b,"grad")
  b=as.numeric(b) ## transform to a vector
  d=c(A%*%b)[j]
  K=numeric(nrow(A))
  K[j]=1
  attr(d,"grad")=t(A)%*%K
  attr(d,"grad2")=matrix(0,ncol=length(b),nrow=length(b))
  class(d)="adv"
  d
}

linear.adv=function(A,b){## derivative of Ab w.r.t b, A is a matrix or a column vector
  grad.b=attr(b,"grad")
  b=as.numeric(b)
  check=length(A)==length(b)
  if(check){
    d=t(A)%*%b
    attr(d,"grad")=A
    attr(d,"grad2")=matrix(0,ncol=length(b),nrow = length(b))
  }else{
    d=A%*%b
    attr(d,"grad")=A
  }
  
  class(d)="adv"
  d
}

exp.adv=function(b){ ## the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=exp(b)
  attr(d,"grad")=exp(b)*grad.b
  attr(d,"grad2")=exp(b)*tcrossprod(grad.b)+exp(b)*grad2.b
  class(d)="adv"
  d
}



"*.adv"=function(a1,a2){#derivative of a1(b)*a2(b) w.r.t. b, where a1 is a scalar function of b and a2 is a scalar or vector function of b.
  
  if(is.null(attr(a1,"grad"))){
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a2=as.numeric(a2)
    d=a1*a2
    attr(d,"grad")=a1*grad.a2
    attr(d,"grad2")=a1*grad2.a2
    
  }else{
    grad.a1=attr(a1,"grad")
    grad2.a1=attr(a1,"grad2")
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a1=as.numeric(a1)
    a2=as.numeric(a2)
    d=a1*a2
    check=length(a2)==1
    if(check){
      attr(d,"grad")=grad.a1*a2+a1*grad.a2
      attr(d,"grad2")=grad.a1%*%t(grad.a2)+grad.a2%*%t(grad.a1)+a2*grad2.a1+a1*grad2.a2
    }else{
      attr(d,"grad")=a2%*%t(grad.a1)+a1*grad.a2
    }
  }
  
  class(d)="adv"
  d
}


"+.adv"=function(a,b){ #derivative of a(b)+b(b) where a and b are vector or scalar functions of b, and a(.) can be not function of b
  if(is.null(attr(a,"grad"))){
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.b
    attr(d,"grad2")=grad2.b
    
  }else{
    grad.a=attr(a,"grad")
    grad2.a=attr(a,"grad2")
    
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    
    a=as.numeric(a)
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.a+grad.b
    attr(d,"grad2")=grad2.a+grad2.b
    
  }
  
  class(d)="adv"
  d
}

##################Generate data 

gendat=function(s,obstim,obsmax,gammatrue,alpha1true,alpha2true,betatrue,D0,sigmatrue,knots,Q.2){
  set.seed(s)
  basetrue=function(t){ifelse(t<10/2,0,exp(-2))} ## 0 when t<10
  b=mvrnorm(m,mu=rep(0,q),Sigma=D0) # random effects in longitudinal submodel
  W=sample(c(0,1),size=m,replace = TRUE)
  desmat=mycubicbs(seq(10/2,obsmax,by=0.005),internal_knots=internal_knots,boundary_knots=boundary_knots)$mat# design matrix  
  set.seed(1)
  TD.fun=function(i,t){# i is a scale, t can be a vector
    tt=function(t){
      c(basetrue(t)*exp(gammatrue*W[i]+alpha1true*desmat[match(t,seq(10/2,obsmax,by=0.005)),]%*%(betatrue+b[i,])+
                          alpha2true*sqrt(t(betatrue+b[i,])%*%Q.2%*%R.2(t,knots)%*%t(Q.2)%*%(betatrue+b[i,]))))
    }
    #hazard=sapply(t,tt)
    #p=1-exp(-hazard*0.01)
    r=rbinom(n=length(t),size=1,prob=1-exp(-sapply(t,tt)*0.005)) ##prob: event occurs
    Time1=obsmax
    delta1=0
    if(max(r)>0){
      Time1=t[min(which(r==1))]
      delta1=1
    }
    return(list(Time1=Time1,delta1=delta1))
  }
  
  TD=sapply(c(1:m),TD.fun,t=seq(10/2,obsmax,by=0.005))
  Time1=as.numeric(TD[1,])
  delta1=as.numeric(TD[2,])
  ####generate censoring 
  #censtim=sample(seq(12,20,by=0.001),size=m,replace = TRUE)
  censtim=sample(seq(12/2,30/2,by=0.001),size=m,replace = TRUE)
  delta=ifelse(Time1<=censtim,delta1,0)
  sum(delta)
  Time=ifelse(Time1<=censtim,Time1,censtim)
  #############generate longitudinal data
  set.seed(1)
  Y=c()
  l=numeric(m)
  time=c()
  for (i in 1:m){
    l[i]=floor(2*Time[i])+1
    for(j in seq(0,Time[i],by=0.5)){
      YY=longit.true(j,fix_eff=betatrue,ran_eff=b[i,])+rnorm(1,sd=sigmatrue)
      Y=rbind(Y,YY)
      time=rbind(time,j)
    }
  }
  
  Y=as.vector(Y)
  time=as.vector(time)
  ### construct data frame. Note here we use "id" not "sub"
  data=data.frame(id=rep(c(1:m),l),Time=rep(Time,l),delta=rep(delta,l),Y=Y,W=rep(W,l),
                  time=time, start=time,stop=time+1/2,event=numeric(length(Y))) ##should adjust "stop" and "event"
  
  data$stop=ifelse(data$stop<=data$Time,data$stop,data$Time)
  data$event=ifelse((data$stop==data$Time)&(data$delta==1),1,data$event)
  return(data)
  
}



############################# Get initial value ##################
inival=function(data,data.id,ctl,knots,Q.2){
  lme.data=lme(fixed=Y~-1+time1+time2+time3+time4+time5+time6+time7,random=list(id=pdDiag(form=~-1+time1+time2+time3+time4+time5+time6+time7)),data = data,control=ctl)
  cox.data=coxph(Surv(Time,delta)~W,data=data.id,x=TRUE)
  beta=lme.data$coefficients[[1]]
  sigma2=lme.data$sigma^2
  D=diag(c(as.numeric(VarCorr(lme.data)[1]),as.numeric(VarCorr(lme.data)[2]),as.numeric(VarCorr(lme.data)[3]),
           as.numeric(VarCorr(lme.data)[4]),as.numeric(VarCorr(lme.data)[5]),as.numeric(VarCorr(lme.data)[6]),as.numeric(VarCorr(lme.data)[7])))
  
  gamma=cox.data$coefficients[1]
  #### two stage ####
  data$stop[data$start==data$stop]=data$stop[data$start==data$stop]+0.01
  data.ts=data
  l=as.data.frame(table(data$id))$Freq
  raneff=apply(random.effects(lme.data),2,rep,times=l)
  data.ts$preY=rowSums(data[,10:16]*t(beta+t(raneff)))
  data.ts$prevar=sapply(c(1:(dim(data)[1])),function(x) sqrt(t(as.numeric(beta+raneff[x,]))%*%(Q.2%*%R.2(data$time[x],knots)%*%t(Q.2))%*%(as.numeric(beta+raneff[x,]))))
  coxts.data=coxph(Surv(start, stop, event)~W+preY+prevar,data=data.ts)
  coxts=coxts.data$coefficients
  cumbase=rbind(c(0,0),basehaz(cox.data, centered=FALSE))
  out=list(beta=beta,sigma2=sigma2,D=D,gamma=gamma,cumbase=cumbase,coxts=coxts)
}


###############Estimation##########

est=function(data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,knots,Q.2){
  Time=data$Time[!duplicated(data$id)]
  delta=data$delta[!duplicated(data$id)]
  W=data$W[!duplicated(data$id)]
  des.Y=as.matrix(data[,c(paste0("time",c(1:q)))])
  des.T=mycubicbs(data.id$Time[data.id$delta==1],internal_knots=c(5,10,15)/2,boundary_knots=c(0,20)/2)$mat
  N=nrow(des.Y)
  l=as.data.frame(table(data$id))$Freq
  diXZ=matrix(0,ncol=q*m,nrow=N) ##block diagnal version of XZ
  for(i in 1:m){
    diXZ[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),c((q*(i-1)+1):(q*i))]=des.Y[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),]
  }
  ### collection of K(t)=Q.2%*%R.2(t)%*%t(Q.2) for t=data.id$Time[data.id$delta==1]
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) Q.2%*%R.2(t,knots)%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  ### EM algorithm  ###
  K=10
  gamma_set=numeric(K)
  alpha1_set=numeric(K)
  alpha2_set=numeric(K)
  beta_set=matrix(0,ncol=q,nrow=K)
  sigma2_set=numeric(K)
  D_set=array(0,dim=c(q,q,K))
  Q.fun_set=numeric(K)
  logLik_set=numeric(K)
  b_set=matrix(0,ncol=q,nrow=m)
  inv.Fish_set=array(0,dim=c(q,q,m)) ## V(bi)~=inv.Fishi
  ## Estimation
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  for (k in 1:K){
    for(i in 1:m){
      Zi=des.Y[data$id==i,]
      Xi=Zi
      Yi=data$Y[data$id==i]
      chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      Q=solve(chol_invSig)
      bi=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta)) #E(bi|Yi,theta),initial value of bi
      a=unique(data.id$Time[(data.id$Time<=Time[i])&(data.id$delta==1)])
      if(length(a)>0){
        if(length(a)>1){
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*t(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y%*%(beta+bi)))
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=tensor(array(apply(ZK,1,tcrossprod),dim=c(q,q,length(a)))+
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y),dim=c(q,q,length(a)))-
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-3/2)*tcrossprod(y%*%(beta+bi))),dim=c(q,q,length(a))),
                         wei,3,1)+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+binew)%*%y%*%(beta+binew)))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
        }else{
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=wei*(tcrossprod(ZK)-alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                         alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])])+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*sqrt(t(beta+binew)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+binew))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
          
        }
        
      }else{
        Fishi=t(Zi)%*%Zi/sigma2+solve(D)
      }
      
      b_set[i,]=bi
      # Fish_set[,,i]=Fishi
      inv.Fish_set[,,i]=solve(Fishi)
    }
    
    tr=0 #sum of trace
    for(i in 1:m){
      tra=sum(diag(t(des.Y[data$id==i,])%*%des.Y[data$id==i,]
                   %*%inv.Fish_set[,,i]))
      tr=tr+tra
    }
    
    ###################### Updata parameters using b_set
    sgamma=0
    salpha1=0
    salpha2=0
    sbeta=numeric(q)
    Igamma=0
    Ialpha1=0
    Ialpha2=0
    Ibeta=matrix(0,ncol=q,nrow=q)
    Igamalp1=0
    Igamalp2=0
    Ialp12=0
    
    
    for(i in 1:m){
      
      if(i %in% unique(data$id[data$delta==1])){
        
        rs=riskset(Time[i])
        Exp.f=numeric(length(rs))
        Exp.fBb=numeric(length(rs))
        Exp.fBb2=numeric(length(rs))
        Exp.fbetaKb=numeric(length(rs))
        sqrbetaKb=numeric(length(rs))
        Exp.sqrfbetaKb=numeric(length(rs))
        Exp.fBbbetaKb=numeric(length(rs))
        Exp.scorebeta=matrix(0,ncol=q,nrow=length(rs))
        
        
        K2=K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]
        desT=des.T[match(Time[i],data.id$Time[data.id$delta==1]),]
        
        betabi=adv(c(beta+b_set[i,]))
        misqrbetaKbKb=c(c(t(beta+b_set[i,])%*%K2%*%(beta+b_set[i,]))^(-1/2)*K2%*%(beta+b_set[i,]))+
          sapply(c(1:q),function(s) sum(diag(inv.Fish_set[,,i]%*%attr((quad.adv(K2,betabi))^(-1/2)*linear_sca.adv(K2,betabi,s),"grad2")))/2)
        
        
        
        for(j in rs){
          f=c(exp(alpha1*desT%*%b_set[j,]+
                    alpha2*sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))))
          inv.Fishj=inv.Fish_set[,,j]
          
          
          fg1=c(f*(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])))
          
          fg2=f*(tcrossprod(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))-
                   alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
                   alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2)
          
          h=c(sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,])))
          
          hg1=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])
          
          hg2=-c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
            c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2
          
          Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
          
          Bb=c(t(desT)%*%b_set[j,])
          
          Exp.fBb[which(rs==j)]=f*Bb+
            sum(diag(inv.Fishj%*%(fg2*Bb+fg1%*%t(desT)+desT%*%t(fg1))))/2
          
          Exp.fBb2[which(rs==j)]=f*Bb^2+
            sum(diag(inv.Fishj%*%(Bb^2*fg2+2*Bb*(fg1%*%t(desT)+desT%*%t(fg1))+2*f*desT%*%t(desT))))/2
          
          betaKb=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))
          
          Exp.fbetaKb[which(rs==j)]=f*betaKb+
            sum(diag(inv.Fishj%*%(betaKb*fg2+2*fg1%*%t(K2%*%(beta+b_set[j,]))+2*(K2%*%(beta+b_set[j,]))%*%t(fg1)+2*f*K2)))/2
          
          sqrbetaKb[which(rs==j)]=h+sum(diag(inv.Fishj%*%hg2))/2
          
          Exp.sqrfbetaKb[which(rs==j)]=f*h+sum(diag(inv.Fishj%*%(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)))/2
          
          Exp.fBbbetaKb[which(rs==j)]=f*h*Bb+sum(diag(inv.Fishj%*%(Bb*(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)+(h*fg1+f*hg1)%*%t(desT)+desT%*%t((h*fg1+f*hg1)))))/2
          
          betabj=adv(c(beta+b_set[j,]))
          Exp.scorebeta[which(rs==j),]=c(f*alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))+
            sapply(c(1:q),function(s) sum(diag(inv.Fishj%*%attr(exp(alpha1*linear.adv(desT,adv(b_set[j,])))*exp(alpha2*(quad.adv(K2,betabj))^(1/2))*(alpha2*(quad.adv(K2,betabj))^(-1/2))*linear_sca.adv(K2,betabj,s),
                                                                "grad2")))/2)
          
        }
        
        ### score of gamma; Fisher of gamma
        d=exp(gamma*W[rs])#*(exp(alpha*xbeta))# a vector, can be cancelled 
        Deno=c(t(d)%*%Exp.f)
        ssgamma=W[i]-t(W[rs]*d)%*%Exp.f/Deno
        sgamma=sgamma+ssgamma
        igamma=t(W[rs]^2*d)%*%Exp.f/Deno-(t(W[rs]*d)%*%Exp.f/Deno)^2
        Igamma=Igamma+igamma
        
        ### score of alpha1; Fisher of alpha1
        xbeta=matrix(rep(c(t(desT)%*%beta),length(rs)),ncol=1)#collection of Xj(Ti)%*%beta
        ssalpha1=desT%*%(beta+b_set[i,])-t(d)%*%(xbeta*Exp.f+Exp.fBb)/Deno
        salpha1=salpha1+ssalpha1
        ialpha1=t(d)%*%(xbeta*Exp.fBb+Exp.fBb2)/Deno-(t(d)%*%(xbeta*Exp.f+Exp.fBb))*(t(d)%*%Exp.fBb)/(Deno^2)
        Ialpha1=Ialpha1+ialpha1
        
        ### score of alpha2; Fisher of alpha2
        ssalpha2=sqrbetaKb[which(rs==i)]-t(d)%*%Exp.sqrfbetaKb/Deno
        salpha2=salpha2+ssalpha2
        ialpha2=t(d)%*%Exp.fbetaKb/Deno-(t(d)%*%Exp.sqrfbetaKb/Deno)^2
        Ialpha2=Ialpha2+ialpha2
        
        ##############
        igamalp1=t(d)%*%(Exp.fBb*W[rs])/Deno-(t(d)%*%Exp.fBb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
        Igamalp1=Igamalp1+igamalp1
        
        igamalp2=t(d)%*%(Exp.sqrfbetaKb*W[rs])/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
        Igamalp2=Igamalp2+igamalp2
        
        ialp12=t(d)%*%(Exp.fBbbetaKb+c(desT%*%beta)*Exp.sqrfbetaKb)/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(d)%*%(Exp.fBb+c(desT%*%beta)*Exp.f))/(Deno^2)
        Ialp12=Ialp12+ialp12
        
        ### part of score of beta; part of Fisher of beta
        
        ssbeta=alpha1*c(desT)+alpha2*misqrbetaKbKb-
          (Deno*alpha1*desT+colSums(d*Exp.scorebeta))/Deno+c(1/sigma2*t(des.Y[data$id==i,])%*%(data$Y[data$id==i]-des.Y[data$id==i,]%*%(beta+b_set[i,])))
        sbeta=sbeta+ssbeta
        ibeta=tcrossprod(ssbeta)
        Ibeta=Ibeta+ibeta
      }else{
        
        ssbeta=c(1/sigma2*t(des.Y[data$id==i,])%*%(data$Y[data$id==i]-des.Y[data$id==i,]%*%(beta+b_set[i,])))
        sbeta=sbeta+ssbeta
        ibeta=tcrossprod(ssbeta)
        Ibeta=Ibeta+ibeta
      }
      
    }
    
    
    Igamalp12=matrix(0,ncol=3,nrow=3)
    Igamalp12[upper.tri(Igamalp12)]=c(Igamalp1,Igamalp2,Ialp12)
    Igamalp12=Igamalp12+t(Igamalp12)
    diag(Igamalp12)=c(Igamma,Ialpha1,Ialpha2)
    
    
    par.sur=c(gamma,alpha1,alpha2)+solve(Igamalp12,c(sgamma,salpha1,salpha2))
    gammanew=par.sur[1]
    alpha1new=par.sur[2]
    alpha2new=par.sur[3]
    betanew=c(beta+solve(Ibeta)%*%sbeta)
    
    
    gamma=gammanew
    alpha1=alpha1new
    alpha2=alpha2new
    beta=betanew
    
    
    sigma2=1/N*(t(data$Y-des.Y%*%beta)%*%(data$Y-des.Y%*%beta-2*diXZ%*%c(t(b_set)))+tr+t(diXZ%*%c(t(b_set)))%*%diXZ%*%c(t(b_set)))
    sigma2=c(sigma2)
    D=1/m*(apply(inv.Fish_set,c(1,2),sum)+t(b_set)%*%b_set)
    D=diag(diag(D))
    
    
    cc=c()
    for(t in sort(unique(data$Time[data$delta==1]))){
      rs=riskset(t)
      Exp.f=numeric(length(rs))
      for(j in rs){
        f=c(exp(alpha1*des.T[match(t,data.id$Time[data.id$delta==1]),]%*%b_set[j,]+
                  alpha2*sqrt(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))))
        
        inv.Fishj=inv.Fish_set[,,j]
        
        fg2=f*(tcrossprod(alpha1*des.T[match(t,data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))-
                 alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))+
                 alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])])
        
        
        Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
      }
      xbeta=matrix(rep(des.T[match(t,data.id$Time[data.id$delta==1]),],length(rs)),ncol=q,byrow=TRUE)%*%beta #collection of Xj(t)%*%beta1
      cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/(c(t(exp(gamma*W[rs]+alpha1*xbeta))%*%Exp.f)))
    }
    
    cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
    cumbase[,1]=cumsum(c(0,cc))
    
    ### Log-likelihood 
    set.seed(1)
    L=2000
    log.p.Y=numeric(m)
    log.s_set=matrix(0,nrow=m,ncol=L)
    ZTimeb=matrix(0,nrow=sum(delta==1),ncol=L) ##Z_Time%*%b
    bKb=matrix(0,nrow=sum(delta==1),ncol=L)
    rmc=rmvnorm(L,mean=rep(0,q))
    for(i in 1:m){
      Zi=des.Y[data$id==i,]
      Xi=Zi
      Yi=data$Y[data$id==i]
      Sigma=solve(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      mu=c(Sigma%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta))
      log.p.Y[i]=dmvnorm(Yi,mean=c(Xi%*%beta),sigma=sigma2*diag(length(Yi))+Zi%*%D%*%t(Zi),log=TRUE)
      rmci=t(mu+t(chol(Sigma))%*%t(rmc))## matrix with ncol=q,nrow=L, MC sample from bi|yi
      if(delta[i]==1){
        ZTimeb[match(i,which(delta==1)),]=des.T[match(i,which(delta==1)),]%*%t(rmci)
        bKb[match(i,which(delta==1)),]=colSums((t(rmci)+beta)*(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
      }
      
      a=unique(data.id$Time[(data.id$Time<=Time[i])&(data.id$delta==1)])
      if(length(a)>0){
        if(length(a)==1){
          btKb=colSums((t(rmci)+beta)*(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
        }else{
          btKb=apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) colSums((t(rmci)+beta)*(y%*%(t(rmci)+beta))))
        }
        log.s_set[i,]=apply(sapply(a,lambda0)*exp(gamma*W[i])*
                              exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(t(rmci)+beta)+
                                    alpha2*t(sqrt(btKb))),2,sum)
        
      }
      
    }
    
    
    log.hazard=matrix(0,ncol=L,nrow=m)
    log.hazard[which(delta==1),]=log(sapply(Time[delta==1],lambda0))+gamma*W[delta==1]+alpha1*(c(des.T%*%beta)+ZTimeb)+
      alpha2*sqrt(bKb)
    log.p.Tb=log.hazard-log.s_set
    p.Tb=exp(log.p.Tb)
    logLikmc=sum(log.p.Y+log(rowMeans(p.Tb)))
    print(logLikmc)
    
    
    if(k>1){
      if((logLikmc<logLik_set[k-1])|((abs(logLikmc-logLik_set[k-1])/abs(logLik_set[k-1]))<10^(-8))){
        #if(((abs(logLikmc-logLik_set[k-1])/abs(logLik_set[k-1]))<10^(-4))){
        break
      }
    }
    
    logLik_set[k]=logLikmc
    gamma_set[k]=gamma
    alpha1_set[k]=alpha1
    alpha2_set[k]=alpha2
    beta_set[k,]=beta
    sigma2_set[k]=sigma2
    D_set[,,k]=D
    #Q.fun_set[k]=Q.fun.up
  }
  
  #out=list(iter.num=k-1,loglike=logLik_set[k-1],gamma=gamma_set[k-1],alpha1=alpha1_set[k-1],alpha2=alpha2_set[k-1],beta=beta_set[k-1,],sigma2=sigma2_set[k-1],D_diag=diag(D_set[,,k-1]))
  return(c(k-1,logLik_set[k-1],gamma_set[k-1],alpha1_set[k-1],alpha2_set[k-1],beta_set[k-1,],sigma2_set[k-1],diag(D_set[,,k-1])))
}
