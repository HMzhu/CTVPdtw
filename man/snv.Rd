\name{snv}
\alias{snv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard Normalized Variable
}
\description{
Standard Normalized Variable(snv) is used to correct the spectral difference caused by scattering between samples.
}
\usage{
snv(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data to be corrected
}
}
\details{
The snv algorithm can correct spectral differences caused by sample inhomogeneities.
}
\value{
Return a list of result
  \item{h}{Corrected data}

}

\examples{
data(corn)
X = corn$X.MP5
h<-snv (X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ snv }% use one of  RShowDoc("KEYWORDS")
\keyword{ Standard Normalized Variable }% __ONLY ONE__ keyword per line
