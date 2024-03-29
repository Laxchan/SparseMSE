% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{gethiermodels}
\alias{gethiermodels}
\title{Get a list of all hierarchical models for given number of lists and maximum order}
\usage{
gethiermodels(nlists, maxorder = nlists - 1, modelvec = hiermodels)
}
\arguments{
\item{nlists}{Number of lists}

\item{maxorder}{Maximum order of models to be returned (defaults to nlists-1)}

\item{modelvec}{vector of hierarchical models (defaults to hiermodels)}
}
\value{
A list of models satisfying the given criteria
}
\description{
Extracts from a larger vector of hierarchical models the ones satisfying the given criterion
}
\examples{
data(hiermodels)
# Five lists with maximum order of 4
gethiermodels(nlists=5,maxorder=4)
# Five lists with maximum order of 2
gethiermodels(nlists=5, maxorder=2)

}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2024).
Bootstrapping Multiple Systems Estimates to Account for Model Selection
\emph{Statistics and Computing}, \strong{34(44)},
Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
}
