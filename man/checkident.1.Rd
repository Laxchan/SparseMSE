% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{checkident.1}
\alias{checkident.1}
\title{The Fienberg-Rinaldo linear program check for the existence of the estimates}
\usage{
checkident.1(parset, datlist)
}
\arguments{
\item{parset}{Either the hierarchical representation of the model, or the vector of the corresponding capture histories to the model.}

\item{datlist}{The output of \code{\link{ingest_data}} on the data set}
}
\value{
The value of the linear program. The parameter estimates within the extended ML framework exist if and only if this value is nonzero.
}
\description{
This routine performs the Fienberg-Rinaldo linear program check in the framework of extended maximum likelihood estimates,
the parameters estimates exist if and only if the return value of the check is nonzero
}
\examples{
data(Korea)
datlist = ingest_data(Korea)
parset ="[0,0,1]"
checkident.1(parset,datlist)

}
\references{
Silverman, B. W., Chan, L. and  Vincent, K., (2022).
Bootstrapping Multiple Systems Estimates to Account for Model Selection

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
models. \emph{Ann. Statist.} \strong{40}, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
}
