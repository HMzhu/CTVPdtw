\name{wheat}
\alias{wheat}
\docType{data}
\title{
This dataset is the standard-free dataset. This is the traditional shootout dataset of IDRC 2016.
}
\description{
248 wheat samples were measured on 9 NIR instruments for calibration and test. Another 150 wheat samples were measured on 3 new instruments (one new instrument per manufacturer) for validation. The wheat samples and reference protein values were provided by the U.S. Department of Agriculture’s Grain Inspection, Packers and Stockyards Administration.
}
\usage{
data(wheat)
}
\format{
  fourteen objects are generated:
  \describe{
    \item{A1.CAL}{matrix with 248 rows and 397 columns}
    \item{A2.CAL}{matrix with 248 rows and 397 columns}
    \item{A3.CAL}{matrix with 248 rows and 397 columns}
    \item{A.Test}{matrix with 248 rows and 397 columns}
    \item{B1.CAL}{matrix with 248 rows and 397 columns}
    \item{B2.CAL}{matrix with 248 rows and 397 columns}
    \item{B3.CAL}{matrix with 248 rows and 397 columns}
    \item{B.Test}{matrix with 248 rows and 397 columns}
    \item{C.Test}{matrix with 248 rows and 100 columns}
    \item{vliadA}{matrix with 150 rows and 397 columns}
    \item{vliadB}{matrix with 150 rows and 397 columns}
    \item{vliadC}{matrix with 150 rows and 100 columns}
    \item{Y.CAL}{matrix with 248 rows and 1 columns}
    \item{y.Valid}{matrix with 150 rows and 1 columns}
  }
}
\details{
This dataset is the standard-free dataset. The data set can be used for calibration transfer.
}

\examples{
data(wheat)
names(wheat)
}
\keyword{ datasets }

