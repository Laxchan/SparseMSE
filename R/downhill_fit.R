<<<<<<< HEAD
#' Conduct downhill search among hierarchical models starting
#'   from the main effects only.
#'
#' Find a local optimum by downhill search among hierarchical models
#'
#' @param counts  Observed counts for the capture histories defined by desmat
#' @param desmat  Incidence matrix defining the capture histories observed with counts given by counts
#' @param maxorder Maximum order of models to be included
#' @param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#' @param niter Number of iterations
#' @param verbose Specifies the output, if F then only returns the best value, if T, returns a more detailed list of objects
#'
#' @return A list with the following components
#'\describe{
#' \item{optimum_hierarchy}{Optimal hierarchical model}
#'   \item{minimum_value}{hierarchical model with the minimum value}
#'   \item{hierarchies_considered}{hierarhical models considered}
#'   \item{function_values}{Values of function}
#'  }
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' xdata=Korea
#' counts = xdata[,dim(xdata)[2]]
#' desmat = xdata[,1:(dim(xdata)[2]-1)]
#' downhill_fit(counts, desmat)
#'
#'@export
downhill_fit = function(counts, desmat, maxorder=dim(desmat)[2]-1, checkid=T, niter=20,verbose=F) {
  # initialise
 nlists = dim(desmat)[2]
 xdata = ingest_data(cbind(desmat, counts))
 fhm = function(hiermod, xdata, checkid) {
      zfit = fit_hier_model(xdatin=xdata, hiermod=hiermod, checkid=checkid)
      zret = matrix(c(zfit$bic, zfit$abundance), nrow=2, dimnames = list( c("BIC", "abundance"), hiermod))
      return(zret)
      }
 # find initial hierarchy
     inithier = paste0("[", paste0(1:nlists, collapse = ","), "]", collapse ="")
  # initialise vector of models considered and values of function
  hiers_considered = inithier
  funvals = fhm(inithier, xdata=xdata, checkid=checkid)
  opthier = inithier
  best_value = funvals
  # now find neighbours excluding those already considered
  for (iter in (1:niter)) {
    newhiers = unlist(find_neighbour_hierarchies(opthier, nlists, keepmaineffects=T,maxorder=maxorder))
    newhiers = setdiff(newhiers, hiers_considered)
    # if no neighbours left, finish search. Otherwise attach new neighbours to those already
    #  considered
    if (length(newhiers) == 0)
      break
    hiers_considered = c(hiers_considered, newhiers)
    #  find minimum among the new neighbours.  If that's more than the current best, finish
    newvals = sapply(newhiers, fhm, xdata=xdata, checkid=checkid)
    funvals = cbind(funvals,newvals)
    znew = min(c(Inf, newvals[1,]), na.rm=T)
    if (znew >= best_value[1])
      break
    #  update minimum to best found so far
    whichmin= which.min(newvals[1,])
    opthier = newhiers[whichmin]
    best_value = newvals[, whichmin]
  }
  if (verbose) { return(
    list(
      optimum_hierarchy = opthier,
      minimum_value = best_value,
      hierarchies_considered = hiers_considered,
      function_values = funvals
    )
  )} else return(best_value[2])
}
#' Bootstrap downhill
#'
#' Construct bootstrap replications and use the downhill fit method to obtain point estimates of total population sizes from each bootstrap sample.
#'
#' @param xdata original data matrix
#' @param nboot number of bootstrap replicates
#' @param iseed random seed
#' @param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#' @param verbose If TRUE, return the list of extra output from \code{downhill_fit}
#' @param maxorder Maximum order of models to be included
#'
#' @return Point estimates of total population sizes from each bootstrap sample.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' downhill_bootstrapcal(Korea)
#' @export
downhill_bootstrapcal <- function(xdata, nboot = 1000, iseed = 1234,
      checkid = T, verbose=F, maxorder=dim(xdata)[2]-2) {
  set.seed(iseed)
  nlists = dim(xdata)[2] - 1
  xdata = tidylists(xdata, includezerocounts = T)
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  nobs = sum(countsobserved)
    # construct all the bootstrap data
  bootreplications = rmultinom(nboot, nobs, countsobserved)
  z = apply(bootreplications, 2, downhill_fit, desmat=desmat,
            maxorder=maxorder, checkid=checkid, niter=20,verbose=verbose)
  return (z)
}
#' Jackknife downhill
#'
#'It uses the downhill approach to calculate the jackknife abundance and returns the estimated
#'acceleration factor
#'
#'@param xdata original data matrix
#'@param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#'@param maxorder Maximum order of models to be included
#'
#'@return the estimated acceleration factor
#'
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' downhill_jackknifecal(Korea)
#'
#' @export
downhill_jackknifecal <- function(xdata,checkid = T, maxorder=dim(xdata)[2]-2) {
  xdata = tidylists(xdata, includezerocounts=F)
  n1= dim(xdata)[1]
  nlists = dim(xdata)[2] - 1
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  jackabund= rep(NA, n1)
  # now the relevant jackknife values for this particular model
  for (j in (1:n1)) {
    yy = countsobserved
    yy[j] = yy[j] - 1
    jackabund[j] = downhill_fit(yy, desmat, maxorder=maxorder, checkid=checkid, niter=20,verbose=F)}
    # now calculate ahat
    jr = sum(countsobserved*jackabund)/sum(countsobserved) - jackabund
    # find estimated acceleration factor by counting each residual the number of times it would occur,
    #   via the count of the corresponding capture history
    ahat = sum(countsobserved*jr^3)/(6 * (sum(countsobserved*jr^2))^{3/2})
  return(ahat)
}
=======
#' Conduct downhill search among hierarchical models starting
#'   from the main effects only.
#'
#' Find a local optimum by downhill search among hierarchical models
#'
#' @param counts  Observed counts for the capture histories defined by desmat
#' @param desmat  Incidence matrix defining the capture histories observed with counts given by counts
#' @param maxorder Maximum order of models to be included
#' @param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#' @param niter Number of iterations
#' @param verbose Specifies the output, if F then only returns the best value, if T, returns a more detailed list of objects
#'
#' @return A list with the following components
#'\describe{
#' \item{optimum_hierarchy}{Optimal hierarchical model}
#'   \item{minimum_value}{hierarchical model with the minimum value}
#'   \item{hierarchies_considered}{hierarhical models considered}
#'   \item{function_values}{Values of function}
#'  }
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' xdata=Korea
#' counts = xdata[,dim(xdata)[2]]
#' desmat = xdata[,1:(dim(xdata)[2]-1)]
#' downhill_fit(counts, desmat)
#'
#'@export
downhill_fit = function(counts, desmat, maxorder=dim(desmat)[2]-1, checkid=T, niter=20,verbose=F) {
  # initialise
 nlists = dim(desmat)[2]
 xdata = ingest_data(cbind(desmat, counts))
 fhm = function(hiermod, xdata, checkid) {
      zfit = fit_hier_model(xdatin=xdata, hiermod=hiermod, checkid=checkid)
      zret = matrix(c(zfit$bic, zfit$abundance), nrow=2, dimnames = list( c("BIC", "abundance"), hiermod))
      return(zret)
      }
 # find initial hierarchy
     inithier = paste0("[", paste0(1:nlists, collapse = ","), "]", collapse ="")
  # initialise vector of models considered and values of function
  hiers_considered = inithier
  funvals = fhm(inithier, xdata=xdata, checkid=checkid)
  opthier = inithier
  best_value = funvals
  # now find neighbours excluding those already considered
  for (iter in (1:niter)) {
    newhiers = unlist(find_neighbour_hierarchies(opthier, nlists, keepmaineffects=T,maxorder=maxorder))
    newhiers = setdiff(newhiers, hiers_considered)
    # if no neighbours left, finish search. Otherwise attach new neighbours to those already
    #  considered
    if (length(newhiers) == 0)
      break
    hiers_considered = c(hiers_considered, newhiers)
    #  find minimum among the new neighbours.  If that's more than the current best, finish
    newvals = sapply(newhiers, fhm, xdata=xdata, checkid=checkid)
    funvals = cbind(funvals,newvals)
    znew = min(c(Inf, newvals[1,]), na.rm=T)
    if (znew >= best_value[1])
      break
    #  update minimum to best found so far
    whichmin= which.min(newvals[1,])
    opthier = newhiers[whichmin]
    best_value = newvals[, whichmin]
  }
  if (verbose) { return(
    list(
      optimum_hierarchy = opthier,
      minimum_value = best_value,
      hierarchies_considered = hiers_considered,
      function_values = funvals
    )
  )} else return(best_value[2])
}
#' Bootstrap downhill
#'
#' Construct bootstrap replications and use the downhill fit method to obtain point estimates of total population sizes from each bootstrap sample.
#'
#' @param xdata original data matrix
#' @param nboot number of bootstrap replicates
#' @param iseed random seed
#' @param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#' @param verbose If TRUE, return the list of extra output from \code{downhill_fit}
#' @param maxorder Maximum order of models to be included
#'
#' @return Point estimates of total population sizes from each bootstrap sample.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' downhill_bootstrapcal(Korea)
#' @export
downhill_bootstrapcal <- function(xdata, nboot = 1000, iseed = 1234,
      checkid = T, verbose=F, maxorder=dim(xdata)[2]-2) {
  set.seed(iseed)
  nlists = dim(xdata)[2] - 1
  xdata = tidylists(xdata, includezerocounts = T)
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  nobs = sum(countsobserved)
    # construct all the bootstrap data
  bootreplications = rmultinom(nboot, nobs, countsobserved)
  z = apply(bootreplications, 2, downhill_fit, desmat=desmat,
            maxorder=maxorder, checkid=checkid, niter=20,verbose=verbose)
  return (z)
}
#' Jackknife downhill
#'
#'It uses the downhill approach to calculate the jackknife abundance and returns the estimated
#'acceleration factor
#'
#'@param xdata original data matrix
#'@param checkid If it is T, then \code{checkident.1} is called and it performs the Fienberg-Rinaldo linear program check for the existence of the estimates
#'@param maxorder Maximum order of models to be included
#'
#'@return the estimated acceleration factor
#'
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' downhill_jackknifecal(Korea)
#'
#' @export
downhill_jackknifecal <- function(xdata,checkid = T, maxorder=dim(xdata)[2]-2) {
  xdata = tidylists(xdata, includezerocounts=F)
  n1= dim(xdata)[1]
  nlists = dim(xdata)[2] - 1
  countsobserved = xdata[, nlists+1]
  desmat = xdata[, 1:nlists]
  jackabund= rep(NA, n1)
  # now the relevant jackknife values for this particular model
  for (j in (1:n1)) {
    yy = countsobserved
    yy[j] = yy[j] - 1
    jackabund[j] = downhill_fit(yy, desmat, maxorder=maxorder, checkid=checkid, niter=20,verbose=F)}
    # now calculate ahat
    jr = sum(countsobserved*jackabund)/sum(countsobserved) - jackabund
    # find estimated acceleration factor by counting each residual the number of times it would occur,
    #   via the count of the corresponding capture history
    ahat = sum(countsobserved*jr^3)/(6 * (sum(countsobserved*jr^2))^{3/2})
  return(ahat)
}
>>>>>>> 60449ffe702bd1b8cd8b68685c88791d5c872b98
