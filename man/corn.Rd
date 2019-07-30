\name{corn}
\alias{corn}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The classical corn dataset with standards from Cargill.
}
\description{
This dataset is the classical dataset with standards. 80 corn samples were measured on 3 different NIR spectrometers by Cargill Corporation, and they are two FOSS NIR Systems 5000 (abbreviated as MP5 and M5) and one FOSS NIR Systems 6000 (abbreviated as MP6).
}
\usage{
data(corn)
}
\format{
  four matrix with 80 objects and 700 variables:
  \describe{
    \item{\code{X.M5}}{matrix with 80 rows and 700 columns}
    \item{\code{X.MP5}}{matrix with 80 rows and 700 columns}
    \item{\code{X.MP6}}{matrix with 80 rows and 700 columns}
    \item{\code{Y}}{matrix with 80 rows and 1 columns}
  }
}

\details{
The data set can be used for calibration transfer.
}

\references{
http://www.eigenvector.com/data/Corn/
}

\examples{
data(corn)
names(corn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }% use one of  RShowDoc("KEYWORDS")

