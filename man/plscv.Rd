\name{plscv}
\alias{plscv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
K-fold Cross-validation for Partial Least Square
}
\description{
Determine the number of optimal latent variables required for partial least squares modeling by interactive test
}
\usage{
plscv(X,Y,A,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{n x p matrix(Sample matrix)}
  \item{Y}{n x m matrix(measured property)}
  \item{A}{The maximal number of latent variables for cross-validation}
  \item{K}{fold. when K=m, it is leave-one-out CV}
}

\value{
Return a list of result
  \item{RMSECV }{Root mean square error obtained under different latent variables}
  \item{RMSEP }{Corresponding minimum root mean square error}
  \item{Q2 }{A parameter for evaluating model prediction performance, the closer Q2 is to 1, the better the model}
  \item{Optlv }{The number of latent variables corresponding to the smallest root mean square error}
%% ...
}

\seealso{\code{\link{plscv}}}
\examples{
data(corn)
X = corn$X.MP5
Y = corn$Y
result<- plscv (X,Y,A=20,K=10)
RMSECV=result$RMSECV
RMSEP=result$RMSEP
Q2=result$Q2
Optlv=result$Optlv
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plscv }% use one of  RShowDoc("KEYWORDS")
\keyword{ K-fold Cross-validation for Partial Least Square }% __ONLY ONE__ keyword per line
