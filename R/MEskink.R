

###----function------------
logis=function(smooth,yd,thres){
  (1+exp(-smooth*(yd-thres)))^(-1)}


##=====================

MEskink=function(y,x,number,Z,V){
  s=c(Z,5,5,V)

  L=min(x)+0.5*abs(mean(x))
  U=max(x)-0.5*abs(mean(x))
  ###############  entropy
  k=3 # number of parameter
  if (number=="3")
  {
    z=cbind(rep(-s[1],k),rep(0,k),rep(s[1],k))    # support beta
    v1=cbind(rep(L,1),rep(mean(x),1),rep(U,1))       # support kink
    v2=cbind(rep(0,1),rep(1,1),rep(s[3],1))       # support smooth
    v3=cbind(rep(-s[4],n),rep(0,n),rep(s[4],n))   # support error
    c=ncol(z)  # number of support
  }

  if (number=="5")
  {
    z =cbind(rep(-s[1],k),rep((-s[1]/2),k),rep(0,k),rep((s[1]/2),k),rep(s[1],k))    # support beta
    v1=cbind(rep(L,1),rep(L+0.5,1),rep(mean(x),1),rep(U-0.5,1),rep(U,1))  		    	  # support kink
    v2=cbind(rep(0,1),rep((0.5),1),rep(1,1),rep(1.5,1),rep(s[3],1)) 			  # support smooth
    v3=cbind(rep(-s[4],n),rep((-s[4]/2),n),rep(0,n),rep((s[4]/2),n),rep(s[4],n))    # support error
    c=ncol(z)  # number of support
  }

  if (number=="7")
  {
    z =cbind(rep(-s[1],k),rep((-s[1]/2),k),rep((-s[1]/3),k),rep(0,k),rep((s[1]/3),k),rep((s[1]/2),k),rep(s[1],k))   # support beta
    v1=cbind(rep(L,1),rep(L+0.25,1),rep(L+0.5,1),rep(mean(x),1),rep(U-0.5,1),rep(U-0.25,1),rep(U,1))     # support kink
    v2=cbind(rep(0,1),rep(0.5,1),rep(0.75,1),rep(1,1),rep(1.5,1),rep(1.75,1),rep(s[3],1))   # support smooth
    v3=cbind(rep(-s[4],n),rep((-s[4]/2),n),rep((-s[4]/3),n),rep(0,n),rep((s[4]/3),n),rep((s[4]/2),n),rep(s[4],n))    # support error
    c=ncol(z)  # number of support
  }

  obj=function(par)
  {
    par1=matrix(par,(n+k+2),c)
    n1=nrow(par1)
    p=abs(par1[(1:3),])
    w1=abs(par1[(4),])
    w2=abs(par1[(5),])
    w3=abs(par1[(6:n1),])
    H=sum(p*log(p))+sum(w1*log(w1))+sum(w2*log(w2))+sum(w3*log(w3))
    sumH=H
    return(sumH)
  }

  ##--------ME reg==========================
  # Maximize Problem

  n=length(y)

  #initial value for p and w
  p=matrix(1/c,k,c)
  w1=matrix(1/c,1,c)
  w2=matrix(1/c,1,c)
  w3=matrix(1/c,n,c)

  par=rbind(p,w1,w2,w3)
  par=c(par)

  con=function(par){

    par1=matrix(par,(n+k+2),c)
    n1=nrow(par1)
    y=((rowSums(rbind(z[1,]*par1[1,])))+((rowSums(rbind(z[2,]*par1[2,]))))*(x*(1-logis((rowSums(rbind(v2[1,]*par1[(5),]))),x,(rowSums(rbind(v1[1,]*par1[4,]))))))+
         ((rowSums(rbind(z[3,]*par1[3,]))))*(x*(logis((rowSums(rbind(v2[1,]*par1[(5),]))),x,(rowSums(rbind((v1[1,]*par1[4,])))))))+
         (rowSums(rbind(v3*par1[(6:n1),])))
    )          # LHS is y
    p=rowSums(par1[(1:n1),])
    return(c(y,p))
  }


  # 1 is sum(P)=1
  # 4 is z2
  powell=solnp(par, fun = obj, eqfun = con, eqB = c(y,rep(1,n+k+2)),LB = c(rep(0,((n+k+2)*c))), UB =c(rep(1,((n+k+2)*c))))
  param=powell$par
  Maxent=max(powell$value)


  #=================================================
  nn=length(y)
  param1=matrix(param,(nn+k+2),c)
  propp=param1[(1:3),]
  propw1=param1[4,]
  propw2=param1[5,]
  propw3=param1[(6:(nn+k+2)),]


  beta0=sum((propp*z)[1,])
  beta1=sum((propp*z)[2,])
  beta2=sum((propp*z)[3,])
  thres1=sum((propw1*v1))
  gam1=sum((propw2*v2))
  par=c(beta0,beta1,beta2,thres1,abs(gam1))

  result=list(beta=par[1:3], threshold=thres1, smooth=abs(gam1),Maxent=Maxent)
  return(result)
}

############## Example smooth kink regression

#library("Rsolnp")
#set.seed(1)
#thres=3
#gam=1.2
#e=rnorm(n)
#x=rnorm(n,thres,5)
#alpha=c(0.5,1,-1)

#y=alpha[1]+(alpha[2]*(x*(1-logis(gam,x,thres))))+(alpha[3]*(x*(logis(gam,x,thres))))+e
#MEskink(y,x,number="5",Z=10,V=5)
