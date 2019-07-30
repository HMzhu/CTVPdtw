\name{KennardStone}
\alias{KennardStone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KennardStone sample partitioning
}
\description{
Use the Euclidean distance between samples to do sample partitioning
%  Input
}
\usage{
KennardStone(x,num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the matrix of the sample spectra}
  \item{num}{the number of the sample spectra you want select }
}
\details{
KennardStone selects the samples XSelected which uniformly distributed in the exprimental data x's space
}
\value{
Return a list of result
  \item{vSelectedRowIndex }{the row index of the selected sample in the X matrix }
  \item{xSelected }{the sample spectras was selected from the X}
  \item{xRest }{the sample spectras remain int the X after select}
  \item{vNotSelectedSample }{the row index of remain int the X after select}
}
\references{
R.W. Kennard, L.A. Stone, Computer Aided Design of Experiments, Technometrics, 11 (1969) 137-148.
}
\author{
ZHuiMin
}
\examples{
data(corn)
x = corn$X.MP5
result<-KennardStone (x,num=64)
X1selected=result$xSelected
X1Rest=result$xRest
vSelectedRowIndex=result$vSelectedRowIndex
vNotSelectedSample=result$vNotSelectedSample
}
\keyword{ ks }
\keyword{ KennardStone }
