\name{plspredtest}
\alias{plspredtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use the built model to predict property values
}
\description{
Use the built model to predict property values
}
\usage{
plspredtest(Wstar,Q,Xtext,xpara1, xpara2, ypara1, ypara2,A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Wstar}{weights for original Sample matrix which use to establish a nipals regression model, n x p matrix(Sample matrix)}
  \item{Q}{weights for original measured property which use to establish a nipals regression model, n x q}
  \item{Xtext}{test spectra, used to test the predictive ability of the model, m x p matrix(Sample matrix)}
  \item{xpara1}{Column mean vector(length=p) of original Sample matrix(n x p)}
  \item{xpara2}{a vector of length p, value of 1}
  \item{ypara1}{Column mean vector(length=q) of original Sample property values}
  \item{ypara2}{a vector of length q, value of 1}
  \item{A}{Number of optimal latent variables}
}

\value{
  \item{ypred}{The value obtained by model prediction}
}

\examples{
data(corn)
X = corn$X.MP5
Y = corn$Y
A = 10
Xtext=corn$X.MP6
result<-pretreat (X)
mx=result$mx
xpara1=result$para1
xpara2=result$para2
result<-pretreat (Y)
my=result$mx
ypara1=result$para1
ypara2=result$para2
res <- pls1_nipals(mx,my,A)
P=res$P
W=res$W
Wstar=W %*% solve(t(P) %*% W)
Q=res$C
xpara1=t(xpara1)
xpara2=t(xpara2)
ypara1=t(ypara1)
ypara2=t(ypara2)
ypred<-plspredtest(Wstar,Q,Xtext,xpara1, xpara2, ypara1, ypara2,A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plspredtest }% use one of  RShowDoc("KEYWORDS")
\keyword{ plspred }% __ONLY ONE__ keyword per line
