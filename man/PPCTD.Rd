\name{PPCDT}
\alias{PPCDT}
\title{An Optimal Subset Selection for Distributed Hypothesis Testing}
\description{
We introduce an optimal subset selection  for distributed hypothesis testing called as PPCDT.
}
\usage{
PPCDT(X,Y,alpha,K)
}
\examples{
alpha=0.05
t=5;K=10;n=1000;p=5
X=matrix(rnorm(n*p,0,1),ncol=p)
beta=matrix(runif(p),nrow = p)
esp=matrix(rnorm(n),nrow = n)
Y=X\%*\%beta+esp
PPCDT(X=X,Y=Y,alpha=alpha,K=K)
}
\author{
Guangbao Guo, Jiarui Li
}
\arguments{
\item{X}{A independent variable}

\item{Y}{The response variable}

\item{alpha}{Significance level}

\item{K}{The number of blocks into which variable X is divided}
}

\value{
 \item{Xopt}{optimal subset of selected independent variables}
 \item{Yopt}{optimal subset of selected response variables}
 \item{Bopt}{Regression coefficients}
 \item{Eopt}{The Mean Squared Error of optimal subset}
 \item{Aopt}{The Mean Absolute Error of optimal subset}
}
