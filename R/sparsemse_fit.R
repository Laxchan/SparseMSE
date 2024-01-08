<<<<<<< HEAD
#' A routine for naive user
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' @param maxorder Maximum order of models to be included.
#' @param nboot Number of bootstrap replications.
#' @param iseed seed for initialisation
#'@param alpha Levels of confidence intervals to be constructed and assessed
#'
#'@return A list with the following components
#'\describe{
#' \item{confvals}{BCa confidence intervals for each considered model}
#'   \item{probests}{Corresponding probabilities of the estimates}
#'  }
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'data(Korea)
#'bootstrap_mse(Korea)
#'
#'@export
bootstrap_mse=function(zdat,maxorder=dim(zdat)[2]-2, nboot=10000, iseed=1234,alpha=c(0.025, 0.1, 0.9, 0.975)){

  if (maxorder == 6){
    z= assemble_bic(zdat, maxorder =2 ,checkexist =T)} else {
      z=assemble_bic(zdat, maxorder=maxorder, checkexist=T)
    }

  if (maxorder == 4){z= subsetmat(z, ntopmodels = 20, maxorder = maxorder)
  } else if(maxorder == 5 | maxorder ==6){
    z =subsetmat(z,ntopmodels = 100, maxorder = maxorder)} else {
      z=subsetmat(z, ntopmodels = Inf, maxorder = maxorder)
    }


  z=bootstrapcal(z, nboot=nboot, iseed=iseed, checkexist=T)
  z=jackknifecal(z,checkexist=T)
  ntopBCa(z, alpha=alpha, maxorder=maxorder)



}
=======
#' A routine for naive user
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' @param maxorder Maximum order of models to be included.
#' @param nboot Number of bootstrap replications.
#' @param iseed seed for initialisation
#'@param alpha Levels of confidence intervals to be constructed and assessed
#'
#'@return A list with the following components
#'\describe{
#' \item{confvals}{BCa confidence intervals for each considered model}
#'   \item{probests}{Corresponding probabilities of the estimates}
#'  }
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'data(Korea)
#'bootstrap_mse(Korea)
#'
#'@export
bootstrap_mse=function(zdat,maxorder=dim(zdat)[2]-2, nboot=10000, iseed=1234,alpha=c(0.025, 0.1, 0.9, 0.975)){

  if (maxorder == 6){
    z= sortmodelsbic(zdat, maxorder =2 ,checkexist =T)} else {
      z=sortmodelsbic(zdat, maxorder=maxorder, checkexist=T)
    }

  if (maxorder == 4){z= subsetmat(z, ntopmodels = 20, maxorder = maxorder)
  } else if(maxorder == 5 | maxorder ==6){
    z =subsetmat(z,ntopmodels = 100, maxorder = maxorder)} else {
      z=subsetmat(z, ntopmodels = Inf, maxorder = maxorder)
    }


  z=bootstrapcal(z, nboot=nboot, iseed=iseed, checkexist=T)
  z=jackknifecal(z,checkexist=T)
  ntopBCa(z, alpha=alpha, maxorder=maxorder)



}
>>>>>>> 60449ffe702bd1b8cd8b68685c88791d5c872b98
