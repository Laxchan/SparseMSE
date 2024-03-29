% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{modelorder}
\alias{modelorder}
\title{Order of models}
\usage{
modelorder(x)
}
\arguments{
\item{x}{a model character string}
}
\value{
the order of the model character string
}
\description{
This routine returns the order of the model given the model character string
}
\examples{
x="[1,2,3,4,5]"
modelorder(x)
y="[12,24,25,35,45]"
modelorder(y)

}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2024).
Bootstrapping Multiple Systems Estimates to Account for Model Selection
\emph{Statistics and Computing}, \strong{34(44)},
Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
}
