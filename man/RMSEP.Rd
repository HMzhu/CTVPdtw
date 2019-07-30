\name{RMSEP}
\alias{RMSEP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Root-Mean-Square Error of Prediction
}
\description{
Used to estimate the error between the predicted value and the standard value
}
\usage{
RMSEP (ytest, ypred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ytest}{property Value obtained by standard methods}
  \item{ypred}{The value obtained by model prediction}
}
\details{
Used to estimate the error between the predicted value and the standard value, that is, calculate the root mean square error of the two sets of data
}
\value{
Return a list of result
  \item{RMSEP }{Root-Mean-Square Error of Prediction}
  \item{Q2 }{A parameter for evaluating model prediction performance, the closer Q2 is to 1, the better the model}
}

\seealso{\code{\link{RMSEP}}}
\examples{
ytest <- c(1,5,4,3,9,8,5,2,6,5,4)
ypred <- c(1.5,5.2,4.3,3.2,8,9,5.4,2.5,5,5,4.7)
result <- RMSEP (ytest, ypred)
RMSEP=result$RMSEP
Q2=result$Q2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RMSEP }% use one of  RShowDoc("KEYWORDS")
\keyword{ Root-Mean-Square Error of Prediction }% __ONLY ONE__ keyword per line
