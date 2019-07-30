\name{savgol}
\alias{savgol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Savitsky-Golay smoothing and differentiation
}
\description{
Savitsky-Golay smoothing and differentiation~
}
\usage{
savgol(X,width,order,deriv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the matrix of ROW vectors to be smoothed}
  \item{width}{the optional variables specifying the number of points in filter}
  \item{order}{the order of the polynomial}
  \item{deriv}{the derivative}
}

\value{
Return a list of result
  \item{D }{the matrix of coefficients (cm) which can be used to create a new smoothed differentiated matrix, i.e. y_hat = y*cm.}
 \item{y_hat }{the matrix of smoothed and differentiated ROW vectors (y_hat)}
}
\references{A. Savitzky, M.J.E. Golay, Smoothing and Differentiation of Data by Simplified Least Squares Procedures, Analytical Chemistry, 36 (1964) 1627-1639.}

\examples{
data(corn)
x = corn$X.MP5
result <- savgol (x,15,2,2)
D=result$D
y_hat=result$y_hat
}

\keyword{ Savitzky-Golay }% use one of  RShowDoc("KEYWORDS")
\keyword{ savgol }% __ONLY ONE__ keyword per line
