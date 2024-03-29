% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{parent_captures}
\alias{parent_captures}
\title{Find the "parents" of a given capture history}
\usage{
parent_captures(k, nlists = 10)
}
\arguments{
\item{k}{An encoded capture history that corresponds to the row number of the capture history data set}

\item{nlists}{The total number of lists}
}
\value{
a vector giving the encoded versions of the parents
}
\description{
Given any encoded capture history and the number of lists, find the encoded capture histories which are obtained by leaving out just one list in turn
}
\examples{
parent_captures(2,10)
parent_captures(1,4)

}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2024).
Bootstrapping Multiple Systems Estimates to Account for Model Selection
\emph{Statistics and Computing}, \strong{34(44)},
Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
}
