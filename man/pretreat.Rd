\name{pretreat}
\alias{pretreat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pretreat dataset before build a model
}
\description{
data pretreatment
}
\usage{
pretreat(x)
}
\arguments{
  \item{x}{Sample matrix(n x p matrix) or measured property(n x q matrix)}
}
\details{
pretreatment method. Contains: autoscaling,pareto,minmax,center or none.Centralized here
}
\value{
Return a list of result
  \item{mx}{matrix after Centralized,(n x p matrix)}
  \item{xpara1 }{Column mean vector of X matrix(length=p)}
  \item{xpara2 }{a vector of length p, value of 1}

}

\examples{
data(corn)
X = corn$X.MP5
result<-pretreat (X)
mx=result$mx
xpara1=result$para1
xpara2=result$para2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{meancenter}% use one of  RShowDoc("KEYWORDS")
\keyword{ pretreat }% __ONLY ONE__ keyword per line
