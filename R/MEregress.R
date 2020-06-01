

###----function------------
MEregress=function(y,x,number,Z,V){
  s=c(Z,5,5,V)                                 #support value
  k=ncol(as.matrix(x))+1

  ############### Belief based entropy

  if (number=="3")
  {
    z=cbind(rep(-s[1],k),rep(0,k),rep(s[1],k))   # support beta
    v3=cbind(rep(-s[4],n),rep(0,n),rep(s[4],n))   # support error
    c=ncol(z)  # number of support
  }

  if (number=="5")
  {
    z =cbind(rep(-s[1],k),rep((-s[1]/2),k),rep(0,k),rep((s[1]/2),k),rep(s[1],k))   # support beta
    v3=cbind(rep(-s[4],n),rep((-s[4]/2),n),rep(0,n),rep((s[4]/2),n),rep(s[4],n))    # support error
    c=ncol(z)  # number of support
  }

  if (number=="7")
  {
    z =cbind(rep(-s[1],k),rep((-s[1]/2),k),rep((-s[1]/3),k),rep(0,k),rep((s[1]/3),k),rep((s[1]/2),k),rep(s[1],k))   # support beta
    v3=cbind(rep(-s[4],n),rep((-s[4]/2),n),rep((-s[4]/3),n),rep(0,n),rep((s[4]/3),n),rep((s[4]/2),n),rep(s[4],n))    # support error
    c=ncol(z)  # number of support
  }

  obj=function(par)
  {
    par1=matrix(par,(n+k),c)
    n1=nrow(par1)
    p=abs(par1[(1:k),])
    w3=abs(par1[((k+1):n1),])
    H=sum(p*log(p))+sum(w3*log(w3))
    sumH=H
    return(sumH)
  }

  ##--------ME reg==========================
  # Maximize Problem
  MElm=function(y,x){
    n=length(y)

    #initial value for p and w
    p=matrix(1/c,k,c)
    w3=matrix(1/c,n,c)

    par=rbind(p,w3)
    par=c(par)

    con=function(par){

      par1=matrix(par,(n+k),c)
      n1=nrow(par1)
      y=cbind(1,x)%*%(rowSums(z*par1[(1:k),]))+(rowSums(v3*par1[((k+1):n1),]))          # LHS is y
      p=rowSums(par1[(1:n1),])
      return(c(y,p))
    }


    # 1 is sum(P)=1
    # 4 is z2
    powell=solnp(par, fun = obj, eqfun = con, eqB = c(y,rep(1,n+k)),LB = c(rep(0,((n+k)*c))), UB =c(rep(1,((n+k)*c))))
    param=powell$par
    value=powell$values

    result=list(param=param,value=value)
    return(result)
  }

  #=================================================
  Entreg=MElm(y=y,x=x)
  nn=length(y)
  prob_ent=Entreg$param
  param=prob_ent
  param1=matrix(param,(nn+k),c)
  propp=param1[(1:k),]
  propw3=param1[((k+1):(nn+k)),]

  B=rep(0,k)

  for ( j in 1:k){
    B[j]=sum((propp*z)[j,])
  }
  ME=length(Entreg$value)
  MaXent=Entreg$value[ME]

  result=list(beta=B, MaXent=MaXent)
  return(result)
}


############## Example regression
#library("Rsolnp")
#set.seed(1)
#n=100
#e=rnorm(n)
#x0=rnorm(n)
#x1=rnorm(n)
#y=1+2*x0+3*x1+e
#x=cbind(x0,x1)
#MEregress(y,x,number="3",Z=10,V=5)
