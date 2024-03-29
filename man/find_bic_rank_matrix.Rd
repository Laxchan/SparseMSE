% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ktop_routines.R
\name{find_bic_rank_matrix}
\alias{find_bic_rank_matrix}
\title{Find BIC rank matrix up to a specified number of bic ranks for a given data set}
\usage{
find_bic_rank_matrix(zsortbic, nbicranks = 5)
}
\arguments{
\item{zsortbic}{Output from applying \code{assemble_bic} to the data}

\item{nbicranks}{The number of bic ranks to be propagated}
}
\value{
A bic rank matrix encoding how many steps the models considered need to get to the optimal one.
}
\description{
Find BIC rank matrix up to a specified number of bic ranks for a given data set
}
\examples{
data(Korea)
zsortbic=assemble_bic(Korea, checkexist=T)
find_bic_rank_matrix(zsortbic, nbicranks=5)
}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2024).
Bootstrapping Multiple Systems Estimates to Account for Model Selection
\emph{Statistics and Computing}, \strong{34(44)},
Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
}
