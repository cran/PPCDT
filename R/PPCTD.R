#' We introduce an optimal subset selection  for distributed hypothesis testing called as PPCDT.
#'
#' @param X is a independent variable
#' @param Y is the response variable
#' @param alpha is significance level
#' @param K is the number of blocks into which variable X is divided
#' @return Xopt,Yopt,Bopt,Eopt,Aopt
#' @examples
#' alpha=0.05
#' t=5;K=10;n=1000;p=5
#' X=matrix(rnorm(n*p,0,1),ncol=p)
#' beta=matrix(runif(p),nrow = p)
#' esp=matrix(rnorm(n),nrow = n)
#' Y=X%*%beta+esp
#' PPCDT(X=X,Y=Y,alpha=alpha,K=K)

PPCDT=function(X,Y,alpha,K){
  n=nrow(X)
  p=ncol(X)
  nk=n/K
  power1=power2=p1=p2=c(1:K)
  Rm=matrix(rep(0, nk*K),ncol=K)
  mr=matrix(rep(0,K*nk), ncol=nk)
  for(i in 1:K )
  {
    mr[i,]=sample(1:n,nk,replace = TRUE)
    r=matrix(c(1:nk,mr[i,]),ncol = nk,byrow = TRUE)
    Rm[,i]=r[2,]
    R=matrix(rep(0,nk*n),ncol=n)
    R[t(r)]=1
    X1=R%*%X
    Y1=R%*%Y
    m=mu0=sdH0=sd=1
    Be=ginv(crossprod(X1))%*%t(X1)%*%Y1
    H=X1%*%ginv(crossprod(X1))%*%t(X1)
    I=diag(rep(1,nk))
    sx=t(Y1)%*%(I-H)%*%Y1/(nk-p)
    A=matrix(rnorm(m*p),nrow=m)
    b=c(1:m)
    BeH0=Be-ginv(crossprod(X1))%*%t(A)%*%ginv(A%*%ginv(crossprod(X1))%*%t(A))%*%(A%*%Be-b)
    sig=sxH0=(t(Y1-X1%*%(BeH0))%*%(Y1-X1%*%(BeH0)))/(n-p)
    power1=pt((mean(X1)-mu0)/(sd(X1)/sqrt(nk)),df=nk-1,lower.tail=TRUE)
    power2=pf(((sig-sx)/(sdH0*m))/(sx/(sd*(nk-p))),df1=m,df2=nk-p)
    T=mean(X1)-mu0/(sd(X1)/sqrt(n))
    p1=pt(T,df=nk-1)
    F=((sig-sx)/(sdH0*m))/(sx/(sd*(nk-p)))
    p2=pf(F,nk-1,nk-1,lower.tail=TRUE)
  }
  t_power=max(power1)
  t_pvalue=min(p1)
  t_power_opt=Rm[,which.max(t_power)]
  t_pvalue_opt=Rm[,which.min(t_pvalue)]

  F_power=max(power2)
  F_pvalue=min(p2)
  F_power_opt=(Rm[,which.max(F_power)])
  F_pvaluer_opt=(Rm[,which.min(F_pvalue)])

  opt=intersect(t_power_opt,t_pvalue_opt)
  Yopt=Y[opt]
  Xopt=X[opt,]
  Bopt=ginv(crossprod(Xopt))%*%t(Xopt)%*%Yopt
  Hopt=Xopt%*%Bopt
  Nopt=length(Yopt)
  Eopt=(t(Yopt-Hopt)%*%(Yopt-Hopt))/Nopt
  Aopt=sum(abs(Yopt-Hopt))/Nopt

  return(list(Xopt=Xopt,Yopt=Yopt,Bopt=Bopt,Eopt=Eopt,Aopt=Aopt))
}
