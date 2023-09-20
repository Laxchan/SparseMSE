#' Downhill wrapper
#'
#' The routine carries out pointwise, bootstrap and jackknife BCa confidence interval calculation and it calls
#' the following routines \code{downhill_fit}, \code{downhill_bootstrapcal}, \code{downhill_jackknifecal}
#' and \code{bcaconfvalues}.
#'
#' @param xdata The data matrix
#' @param maxorder Maximum order of models to be included
#' @param checkid if TRUE, do the Fienburg-Renaldo test
#' @param verbose if TRUE, return the list of extra output from \code{downhill_fit}
#' @param nboot The number of bootstrap replications
#' @param alpha Bootstrap quantiles of interests.
#'
#' @return A list with the following components
#'\describe{
#' \item{point_est}{A point estimate using the downhill fit approach}
#'   \item{boot_res}{Point estimates of total population sizes from each bootstrap sample}
#'   \item{jack_res}{The estimated acceleration factor}
#'   \item{BCaconf_int}{BCa confidence intervals}
#'  }
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2022).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#'
#' @examples
#' data(Korea)
#' xdata=Korea
#' downhill_funs(xdata)
#'
#' @export
downhill_funs = function(xdata, maxorder=dim(xdata)[2]-2, checkid=T, verbose=F, nboot=1000, alpha=c(0.025, 0.05, 0.95, 0.975)){

  counts = xdata[,dim(xdata)[2]]
  desmat = xdata[,1:(dim(xdata)[2]-1)]


  fitres_downhill = downhill_fit(counts = counts, desmat =desmat, maxorder = maxorder, verbose = verbose)

  bootres_downhill = downhill_bootstrapcal(xdata, checkid = checkid, verbose=verbose, maxorder=maxorder,
                                                                   nboot=nboot)
  jackres_downhill = downhill_jackknifecal(xdata, checkid = checkid, maxorder = maxorder)

  BCaconf_interval = bcaconfvalues(bootreps = bootres_downhill, popest = fitres_downhill ,ahat= jackres_downhill, alpha=alpha)

  return(list(
    point_est =  fitres_downhill,
    boot_res = bootres_downhill,
    jack_res = jackres_downhill,
    BCaconf_int = BCaconf_interval
    ))
}


