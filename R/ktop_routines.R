#' Return a BIC rank matrix breaking ties with each roles
#'
#' This routine takes the ouput of \code{find_bic_rank_matrix}
#' and produce a BIC rank matrix with a specified order of rank breaking
#' ties with each previous column.
#'
#' @param BICmatrix_prop The output from \code{find_bic_rank_matrix}
#' @param k The order of BIC rank
#'
#' @return A BIC rank matrix with the last column breaking ties with the previous ones.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' z=sortmodelsbic(Korea)
#' z=jackknifecal(z)
#' BICmatrix_prop= find_bic_rank_matrix(z)
#' BICmatrix_break = BICrank_tiebreak(BICmatrix_prop, 2)
#' @export
BICrank_tiebreak <- function(BICmatrix_prop, k){
  BICmatrix_break <- BICmatrix_prop[,1:k]
  BICmatrix_break <- BICmatrix_break[do.call(order, rev(asplit(BICmatrix_break, 2))),]
  return(BICmatrix_break)
}
#' Find BCa confidence intervals using ktop idea.
#'
#'@param z The output of \code{sortmodelsbic}, \code{bootstrapcal} and \code{jackknifecal}
#'@param alpha Levels of confidence intervals to be constructed and assessed
#'@param BICmatrix_break output from \code{BICrank_tiebreak}
#'@param BCaFD If T, return also the probabilities of the estimates
#'@param maxorder Maximum order of models to be included
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
#' data(Korea)
#' z=sortmodelsbic(Korea)
#' z=bootstrapcal(z)
#' z=jackknifecal(z)
#' BICmatrix_prop= find_bic_rank_matrix(z)
#' BICmatrix_break = BICrank_tiebreak(BICmatrix_prop, 2)
#' ktopBCa(z,BICmatrix_break)
#'
#'@export
ktopBCa = function(z,BICmatrix_break,alpha=c(0.025, 0.05, 0.1, 0.16,0.2, 0.5, 0.8, 0.84, 0.9, 0.95, 0.975), maxorder=Inf, BCaFD=FALSE) {
  # set up and calculate order of model scores
  # restrict to models of prescribed order
  choosemods = (z$res[,3] <= maxorder)
  bootabund = z$bootabund[choosemods,]
  bootbic = z$bootbic[choosemods, ]
  #
  nmods=dim(bootabund)[1]
  modelslist = dimnames(bootabund)[[1]]
  nreps = dim(bootabund)[2]
  #modscores = (1:nmods)
  modscores = BICmatrix_break[,1]
  #if (ktop*lweight > 0) modscores = rescoremodels(modelslist, ktop, lweight, verbose=F)
  #modelorder= order(modscores)
  modelorder = as.vector(modscores)
  # now set up to find all bootstrap estimates
  zba=bootabund[modelorder,]
  zbootest=zba
  zbic = bootbic[modelorder,]
  #
  for (j in (1:nreps)) {
    jrec = which.min(zbic[,j])
    zbootest[jrec:nmods,j] = zba[jrec,j]
    while (jrec >1) {
      jrec0=jrec-1
      jrec = which.min(zbic[1:jrec0,j])
      jrec=max(jrec, 1)
      zbootest[jrec:jrec0,j] = zba[jrec,j]
    }
  }
  # now find acceleration and then bias correction (vector, one entry for each value of ntop)
  ahat = bicktopahatcal(z, modelorder)
  popest = z$res[choosemods,][1,1]
  z0 = qnorm(apply(zbootest < popest, 1, mean))
  za = qnorm(alpha)
  za0mat = outer(z0,za, "+")
  zqmat = pnorm( z0 + za0mat/(1-ahat*za0mat))
  # now find confidence values
  bcac = matrix(NA, nrow= nmods, ncol=length(alpha), dimnames=list(dimnames(zbootest)[[1]], alpha))
  for (i in (1:nmods)) bcac[i,]=bcaconfvalues(zbootest[i,],popest,ahat[i],alpha=alpha)
  #
  # now find the estimated BCa confidence levels, based on the full data,
  #   of the confidence values found
  if (BCaFD == TRUE){
    z01 = z0[nmods]
    cvals = qnorm( ((1:nreps)-0.5)/nreps )
    yvals = pnorm( (cvals-z01)/(1+ahat[nmods]*(cvals-z01))  - z01)
    xvals = sort(zbootest[nmods,])
    pfun = splinefun(xvals,yvals, method="monoH.FC")
    pbcac = pfun(as.vector(bcac))
    pbcac = pmin(pmax(pbcac, 0), 1)
    dim(pbcac)=dim(bcac)
    dimnames(pbcac)=dimnames(bcac)
    return(list(confvals=bcac, probests=pbcac))
  }
  return(confvals=bcac)
}
#' Find BIC rank matrix up to a specified number of bic ranks for a given data set
#'
#'
#' @param zsortbic Output from applying \code{sortmodelsbic} to the data
#' @param nbicranks The number of bic ranks to be propagated
#'
#' @return A bic rank matrix encoding how many steps the models considered need to get to the optimal one.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' zsortbic=sortmodelsbic(Korea, checkexist=T)
#' find_bic_rank_matrix(zsortbic, nbicranks=5)
#' @export
find_bic_rank_matrix= function(zsortbic,  nbicranks=5){
  #  set things up; this could be done directly through arguments to the
  # original function call
  models = dimnames(zsortbic$res)[[1]]
  nlists=dim(zsortbic$xdata)[2]-1
  maxorder = zsortbic$maxorder
  nranks = dim(zsortbic$res)[1]
  #
  # now start off the bic_rank_matrix and construct the neighbour list
  bic_rank_matrix = matrix(NA, nrow=nranks, ncol=nbicranks, dimnames=list(models,1:nbicranks))
  bic_rank_matrix[,1] = rank(zsortbic$res[,2], ties.method="f")
  neighbour_list = sapply(models, find_neighbour_hierarchies, nlists=nlists, maxorder=maxorder,
                          simplify=F)
  #  Now do the recursion to find all the bic ranks.  For debugging purposes
  #  you can call it with nbicranks=1 and the output will be different
  #  Perhaps remove this line in the final version or put it differently.
  if (nbicranks==1) return(list(neighbour_list, bic_rank_matrix))
  #
  for (j in (2:nbicranks)) {
    for (jm in models) {
      #brnew = min(bic_rank_matrix[neighbour_list[[jm]], j-1])
      #this is to make sure the neighbour is actually one of the models
      brnew = min(bic_rank_matrix[intersect(neighbour_list[[jm]],models), j-1])
      bic_rank_matrix[jm,j]  = min(bic_rank_matrix[jm, j-1], brnew)
    } }
  return(bic_rank_matrix)

}
