\name{meancenter}
\alias{meancenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean centering
}
\description{
Data preprocessing to reduce the difference caused by different dimensions, self-variation}
\usage{
meancenter(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the matrix of the sample spectra}
}

\value{
\item{result }{result of matrix centered }
}

\author{ZHuiMin}

\examples{
data(corn)
X = corn$X.MP5
h <- meancenter (X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meancenter }% use one of  RShowDoc("KEYWORDS")
