\name{autoscale}
\alias{autoscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{autosacle}
\description{
Autosacle eliminates the impact of dimensions on data
}
\usage{
autoscale(X)
}
\arguments{
  \item{X}{Data to be Processed}
}
\details{
The autosacle eliminates the effect of the dimension on the data, so that the pre-processed spectral data fully reflects the information changes.
}
\value{

  \item{h}{Processed data}
}

\examples{
data(corn)
X = corn$X.MP5
h<-autoscale (X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ auto }% use one of  RShowDoc("KEYWORDS")
\keyword{ autoscale }% __ONLY ONE__ keyword per line
