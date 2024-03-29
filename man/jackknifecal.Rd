% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{jackknifecal}
\alias{jackknifecal}
\title{Jackknife abundance and Jackknife bic}
\usage{
jackknifecal(z, checkexist = T)
}
\arguments{
\item{z}{Results from \code{assemble_bic} or \code{subsetmat}.}

\item{checkexist}{If \code{checkexist=TRUE}, check for existence in cases where the jackknife introduces an additional zero, else it does
not check for existence.  Note that in the current version it is assume that models for which the fit doesn't exist for the original data
have already been excluded.}
}
\value{
A list with the following components
\describe{
  \item{jackabund}{Jackknife abundance matrix}
  \item{jackbic}{Jackknife BIC matrix}
  \item{countsobserved}{Capture counts in the same order as the columns of \code{jackabund} and \code{jackbic}}
}
}
\description{
This routine takes the output from \code{subsetmat} or from \code{assemble_bic} and returns the jackknife
abundance matrix and jackknife BIC matrix.
}
\examples{
data(Korea)
z=assemble_bic(Korea)
jackknifecal(z,checkexist=T)

}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2024).
Bootstrapping Multiple Systems Estimates to Account for Model Selection
\emph{Statistics and Computing}, \strong{34(44)},
Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
}
