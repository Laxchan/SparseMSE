% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{checkident}
\alias{checkident}
\title{Check a model for the existence and identifiability of the maximum likelihood estimate}
\usage{
checkident(zdat, mX = 0, verbose = FALSE)
}
\arguments{
\item{zdat}{Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
are 0s and 1s defining the capture histories
observed. The last column is the count of cases with that particular capture history.
List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
it is assumed that it has observed count zero.}

\item{mX}{A \eqn{2 \times k} matrix giving the \eqn{k} two-list parameters to be included in the model.
Each column of \code{mX} contains the numbers of the corresponding pair of lists.
If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no two-list parameters are included and
 the main effects model is fitted.}

\item{verbose}{Specifies the output.   If \code{FALSE} then the error code is returned.  If \code{TRUE} then
in addition the routine prints an error message if the model/data fail either of the two tests, and also
returns both the error code and the \code{lp} object.}
}
\value{
If \code{verbose=FALSE}, then return the error code \code{ierr} which is 0 if there are no errors, 1 if the linear program test shows that the maximum likelihood
 estimate does not exist, 2 if it is not identifiable, and 3 if both tests are failed.

If \code{verbose=TRUE}, then return a list with components as below

\code{ierr} As described above.

\code{zlp} Linear programming object, in particular giving the value of the objective function at optimum.
}
\description{
Apply the linear programming test as derived by Fienberg and Rinaldo (2012), and a calculation of the rank of the design
matrix, to check whether a particular model yields an identifiable maximum likelihood estimate
based on the given data.  The linear programming problem is as described on page 3 of Fienberg and Rinaldo (2012), with a typographical error corrected.
Further details are given by Chan, Silverman and Vincent (2021).
}
\examples{
data(Artificial_3)
#Build a matrix that contains all two-list effects
m=dim(Artificial_3)[2]-1
mX = t(expand.grid(1:m, 1:m)); mX = mX[ , mX[1,]<mX[2,]]
# When the model is not identifiable
checkident(Artificial_3,mX=mX, verbose=TRUE)
# When the maximum likelihood estimate does not exist
checkident(Artificial_3, mX=mX[,1],verbose=TRUE)
#Model passes both tests
checkident(Artificial_3, mX=mX[,2:3],verbose=TRUE)

}
\references{
Chan, L., Silverman, B. W., and Vincent, K. (2021).
 Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
\emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
models. \emph{Ann. Statist.} \strong{40}, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
}
