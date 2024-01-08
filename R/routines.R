#' Build model for multiple systems estimation
#'
#' For multiple systems estimation model corresponding to a specified set of two-list effects,
#'    set up the GLM model formula and data matrix.
#'
#'
#'@param zdat  Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero counts.
#'
#'@param mX A \eqn{2 \times k} matrix giving the \eqn{k} two-list effects to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list effects are included. If \code{mX = NULL}, no such effects are included and
#'  the main effects model is fitted.
#'
#' @return A list with components as below.
#'
#'  \code{datamatrix} A matrix with all possible capture histories, other than those equal to or containing non-overlapping pairs
#'      indexed by parameters
#'  that are within the model specified by \code{mX}.  A non-overlapping pair is a pair of lists \eqn{(i,j)} such that
#'  no case is observed in both lists,
#'  regardless of whether it is present on any other lists.   If \eqn{(i,j)} is within the model specified by \code{mX},
#'  all capture histories containing both \eqn{i} and \eqn{j} are
#'  then excluded.
#'
#' \code{modelform} The model formula suitable to be called by the Generalized Linear Model function \code{glm}. Model terms corresponding to
#' non-overlapping pairs are not included, because they are handled by removing appropriate rows from the data matrix supplied to glm. The list
#' of non-overlapping pairs are provided in \code{emptyoverlaps}. See Chan, Silverman and Vincent (2021) for details.
#'
#' \code{emptyoverlaps} A matrix with two rows, whose columns give the indices of non-overlapping pairs of lists where the parameter indexed by
#' the pair is within the specified model. The column names give the names of the lists
#' corresponding to each pair.
#'
#' @references
#' Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' @examples
#' data(NewOrl)
#' buildmodel(NewOrl, mX=NULL)
#' #Build a matrix that contains all two-list effects
#' m=dim(NewOrl)[2]-1
#' mX = t(expand.grid(1:m, 1:m)); mX = mX[ , mX[1,]<mX[2,]]
#' # With one two-list effect
#' buildmodel(NewOrl, mX=mX[,1])
#' #With three two-list effects
#' buildmodel(NewOrl, mX=mX[,1:3])
#'
#' @importFrom stats as.formula
#'
#'@export


buildmodel <- function(zdat, mX) {
  m = dim(zdat)[2] - 1
  #Recover the data matrix with combination of lists which have count zero
  zdfull = tidylists(zdat, includezerocounts = TRUE)
  vn = dimnames(zdfull)[[2]]
  #count
  cn = vn[m+1]
  #names of list excluding count
  vn = vn[ -(m+1)]
  modf = paste(cn, "~ .")
  if (is.null(mX)) return(list(datamatrix = zdfull, modelform = modf, emptyoverlaps = matrix(nrow=2,ncol=0)))
  if (sum(mX) == 0) {
    mX = t(expand.grid(1:m, 1:m))
    mX = mX[ , mX[1,]<mX[2,]]
  }
  mX=as.matrix(mX)
  # ---  find if any of the overlaps specified within the model are empty
  neff = dim(mX)[2]
  effc = rep(TRUE,neff)
  for (j in (1:neff)) {
    i1 = mX[1,j]
    i2 = mX[2,j]
    joverlap = sum(zdfull[,i1]*zdfull[,i2]*zdfull[,m+1])
    if (joverlap == 0) {
      effc[j] = FALSE
      zdfull = zdfull[zdfull[,i1]*zdfull[,i2]==0, ]
    }
  }
  #----construct matrix of two-list parameters where overlap of the corresponding lists is empty
  mempty = as.matrix(mX[, !effc])
  dimnames(mempty)[[1]] = c("","")
  menames =  apply(matrix(vn[mempty],nrow=2),2, paste, collapse=":")
  dimnames(mempty)[[2]] = menames
  #----construct matrix of two-list parameters where overlap of the corresponding lists is not empty
  mX = as.matrix(mX[,effc])
  mXc = vn[mX]
  dim(mXc) = dim(mX)
  #----construct the model formula to include
  #      the two-list effects in the model for which the overlap is not empty
  mXmod = paste(apply(mXc,2, paste, collapse=":"), collapse=" + ")
  if (nchar(mXmod) > 0) modf = paste(modf, mXmod, sep= " + ")
  modf = as.formula(modf)
  return(list(datamatrix = zdfull, modelform = modf, emptyoverlaps = mempty))
}
#' Produce a data matrix with a unique row for each capture history
#'
#' This routine finds rows with the same capture history and consolidates them into a single row whose count is the sum of counts of
#' the relevant rows.  If \code{includezerocounts = TRUE} then it also includes rows for all the capture histories with zero count; otherwise
#' these are all removed.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param includezerocounts  If \code{FALSE} then remove rows corresponding to capture histories with zero count.
#' If \code{TRUE} then include all possible capture histories including those with zero count,
#' excluding the all-zero row corresponding to the dark figure.
#'
#' @return A data matrix in the form specified above, including all capture histories with zero counts if  \code{includezerocounts=TRUE}.
#'
#' @examples
#' data(NewOrl)
#' tidylists(NewOrl,includezerocounts=TRUE)
#'
#'@export
tidylists <- function(zdat, includezerocounts = FALSE) {
  zdat = as.matrix(zdat)
  m = dim(zdat)[2] - 1
  # ----construct full capture history matrix
  zm = NULL
  #-----produce an unordered matrix of all possible capture history including the one corresponds to dark figure
  for (j in (1:m)) zm = rbind(cbind(1, zm), cbind(0, zm))
  # ----calculate the number of 1s in each row of the unordered matrix
  ztot = apply(zm, 1, sum)
  #----order the rows of the matrix according to the number of 1s appearing on each row in ascending order
  zm = zm[order(ztot), ]
  #---remove the row that has 0s for all its entries and add a column of 0s
  zm = cbind(zm[-1, ], 0)
  #---supply column names if necessary
  vn = dimnames(zdat)[[2]]
  if(is.null(vn)) {
    vn = c(LETTERS[1:m], "count")
  }
  dimnames(zm)[[2]] = vn
  # find row of zm corresponding to each row of zmse and update count
  bcode = zdat[, -(m + 1)] %*% (2^(1:m))
  bc = zm[, -(m + 1)] %*% (2^(1:m))
  for (j in (1:dim(zdat)[1])) {
    ij = (1:dim(zm)[1])[bc == bcode[j]]
    zm[ij, m + 1] = zm[ij, m + 1] + zdat[j, m + 1]
  }
  # remove rows with zero counts if includezerocounts is F
  if (!includezerocounts)
    zm = zm[zm[, m + 1] > 0, ]
  zdatr = as.data.frame(zm)
  return(zdatr)
}

#' Build the model matrix based on particular data, as required to check for identifiability and existence of the maximum likelihood estimate
#'
#' This routine builds a model matrix as required by the linear program check \code{\link{checkident}} and checks if the matrix is of full rank.
#'  In addition, for each individual list, and for each pair of lists included in the model,
#'   it returns the total count of individuals appearing on the specific list or lists whether or not in combination with other lists.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list parameters to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list parameters are included. If \code{mX = NULL}, no such parameters are included and
#'  the main effects model is fitted.
#'
#' @return A list with components as below
#'
#' \code{modmat} The matrix that maps the parameters in the model (excluding any corresponding to non-overlapping lists) to the log expected value
#'    of the counts of capture histories that do not contain non-overlapping pairs in the data.
#'
#' \code{tvec} A vector indexed by the parameters in the model, excluding those corresponding to non-overlapping pairs of lists.  For each parameter
#'    the vector contains the total count of individuals in all the capture histories that include both the lists indexed by that parameter.
#'
#' \code{rankdef} The column rank deficiency of the matrix \code{modmat}. If \code{rankdef = 0}, the matrix has full column rank.
#'
#'
#'@examples
#' data(NewOrl)
#' buildmodelmatrix(NewOrl, mX=NULL)
#'
#' @export

buildmodelmatrix <- function(zdat, mX=NULL) {
  m = dim(zdat)[2] - 1
  #Recover the data matrix
  zdfull = tidylists(zdat, includezerocounts =TRUE)
  if (length(mX)==0) { mX= matrix(nrow=2,ncol=0)
  }
  else {
    if (sum(mX) == 0) {
      mX = t(expand.grid(1:m, 1:m))
      mX = mX[ , mX[1,]<mX[2,]]
    }
    mX=as.matrix(mX)
    # ---  find if any of the overlaps specified within the model are empty
    neff = dim(mX)[2]
    effc = rep(TRUE,neff)
    for (j in (1:neff)) {
      i1 = mX[1,j]
      i2 = mX[2,j]
      joverlap = sum(zdfull[,i1]*zdfull[,i2]*zdfull[,m+1])
      if (joverlap == 0) {
        effc[j] = FALSE
        zdfull = zdfull[zdfull[,i1]*zdfull[,i2]==0, ]
      }
    }
    #----construct matrix of two-list parameters where overlap of the corresponding lists is not empty
    mX = as.matrix(mX[,effc])
  }
  counts = zdfull[,m+1]
  nobs = dim(zdfull)[1]
  vn = c(0,dimnames(zdfull)[[2]][-(m+1)])
  modmat = cbind(rep(1,nobs), zdfull[, -(m+1)])
  #----construct the model matrix to include
  #      the two-list parameters in the model for which the overlap of the corresponding lists is not empty
  neffs = dim(mX)[2]
  if (neffs > 0) {for ( j in (1:dim(mX)[2])) {
    i1 = mX[1,j]
    i2 = mX[2,j]
    modmat = cbind(modmat, zdfull[,i1]*zdfull[,i2])
    vn = c(vn, paste(vn[i1+1],vn[i2+1],sep=":"))
  }}
  dimnames(modmat)[[2]] =vn
  tvec = t(modmat)%*%counts
  rankdef = dim(modmat)[2] - qr(modmat)$rank
  return(list(modmat=modmat, tvec=tvec, rankdef=rankdef))
}

#' Check a model for the existence and identifiability of the maximum likelihood estimate
#'
#' Apply the linear programming test as derived by Fienberg and Rinaldo (2012), and a calculation of the rank of the design
#' matrix, to check whether a particular model yields an identifiable maximum likelihood estimate
#' based on the given data.  The linear programming problem is as described on page 3 of Fienberg and Rinaldo (2012), with a typographical error corrected.
#' Further details are given by Chan, Silverman and Vincent (2021).
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list parameters to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no two-list parameters are included and
#'  the main effects model is fitted.
#'
#' @param verbose Specifies the output.   If \code{FALSE} then the error code is returned.  If \code{TRUE} then
#' in addition the routine prints an error message if the model/data fail either of the two tests, and also
#' returns both the error code and the \code{lp} object.
#'
#' @return
#'
#'If \code{verbose=FALSE}, then return the error code \code{ierr} which is 0 if there are no errors, 1 if the linear program test shows that the maximum likelihood
#'  estimate does not exist, 2 if it is not identifiable, and 3 if both tests are failed.
#'
#'If \code{verbose=TRUE}, then return a list with components as below
#'
#'\code{ierr} As described above.
#'
#'\code{zlp} Linear programming object, in particular giving the value of the objective function at optimum.
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
#' models. \emph{Ann. Statist.} \strong{40}, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
#'
#'
#' @examples
#' data(Artificial_3)
#' #Build a matrix that contains all two-list effects
#' m=dim(Artificial_3)[2]-1
#' mX = t(expand.grid(1:m, 1:m)); mX = mX[ , mX[1,]<mX[2,]]
#' # When the model is not identifiable
#' checkident(Artificial_3,mX=mX, verbose=TRUE)
#' # When the maximum likelihood estimate does not exist
#' checkident(Artificial_3, mX=mX[,1],verbose=TRUE)
#' #Model passes both tests
#' checkident(Artificial_3, mX=mX[,2:3],verbose=TRUE)
#'
#' @export
#'

checkident <- function(zdat, mX=0, verbose=FALSE)
{
  #---use LP to check if ML estimate exists
  #--- construct model matrix
  zb = buildmodelmatrix(zdat, mX)
  amat = t(zb$modmat)
  tt = zb$tvec
  npar = dim(amat)[1]
  nobs = dim(amat)[2]
  #----construct objective and constraint matrix
  f.obj = c(rep(0,nobs),1)
  const.rhs = c( tt, rep(0,nobs))
  amat1 = cbind(amat, rep(0,npar))
  amat2 = cbind(diag(nobs), rep(-1, nobs))
  #amat3 = c(rep(0,nobs), 1)
  const.mat = rbind(amat1, amat2)
  const.dir = c( rep("=",npar), rep(">=",nobs))
  #--------------carry out the LP step
  zlp = lpSolve::lp("max",f.obj,const.mat,const.dir,const.rhs)
  ierr = (zlp$objval == 0) + 2*(zb$rankdef > 0 )
  if (!verbose) return(ierr)
  if (zlp$objval == 0) cat ("The estimate is not finite \n")
  if (zb$rankdef > 0 ) cat ("The estimate is not identifiable \n")
  return(list(ierr=ierr, lp= zlp))
}

#' Estimate the total population including the dark figure.
#'  If the user wishes to find bootstrap confidence intervals then the routine \code{\link{estimatepopulation}} should be used instead.
#'
#' This routine estimates the total population size, which includes the dark figure, together with confidence intervals as specified.
#' It also returns the details of the fitted model.  The user can choose whether to fit main effects only,
#' to fit a particular model containing specified two-list parameters, or
#' to choose the model using the stepwise approach described by Chan, Silverman and Vincent (2021).
#'
#' @param zdat  Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param method If \code{method = "stepwise"} the stepwise method implemented in \code{\link{stepwisefit}} is used.  If \code{method = "fixed"} then a specified fixed model is used; the model is then given by \code{mX}.  If \code{method = "main"} then main effects only are fitted.
#'
#' @param quantiles Quantiles of interest for confidence intervals.
#'
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list parameters to be included in the model if \code{method = "fixed"}.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list parameters are included. If \code{mX = NULL}, no two-list parameters are included and
#'  the main effects model is fitted.
#' If only one two-list parameter is to be fitted, it is sufficient to specify it as a vector of length 2, e.g \code{mX=c(1,3)}
#' for the parameter indexed by lists 1 and 3.   If \code{method} is equal to \code{"stepwise"} or \code{"main"}, then \code{mX} is ignored.
#'
#' @param pthresh Threshold p-value used if \code{method = "stepwise"}.
#'
#' @return A list with components as below
#'
#' \code{estimate} Point estimate and confidence interval estimates corresponding to specified quantiles.
#'
#' \code{MSEfit} The model fitted to the data in the format described in \code{\link{modelfit}}.
#'
#'@references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' @examples
#' data(NewOrl)
#' data(NewOrl_5)
#' estimatepopulation.0(NewOrl, method="stepwise", quantiles=c(0.025,0.975))
#' estimatepopulation.0(NewOrl_5, method="main", quantiles=c(0.01, 0.05,0.95, 0.99))
#'
#' @importFrom stats qnorm
#'
#'@export
estimatepopulation.0 <-function(zdat,method="stepwise", quantiles=c(0.025,0.975),mX=NULL, pthresh=0.02){
  if (method=="stepwise"){
    zMSE=stepwisefit(zdat, pthresh)
  }
  else if (method=="fixed"){
    zMSE=modelfit(zdat,mX)
  }
  else if (method=="main") {
    zMSE=modelfit(zdat, mX=NULL)
  }
  else {stop ("Unknown method")}
  zfit=zMSE$fit
  zzdat=zfit$data
  totobs = sum(zzdat[, dim(zzdat)[2]])#summing all observed history
  estdarklogmean = zfit$coefficients[1]#intercept of the fit
  estdarklogsd   = sqrt(summary(zfit)$cov.unscaled[1,1])
  quantiles_with_pt.estimate<-sort(unique(c(quantiles,0.5)))
  estdarkci = exp(qnorm(quantiles_with_pt.estimate, estdarklogmean, estdarklogsd))#estimated dark figure confidence interval
  totci = totobs + estdarkci
  names(totci)<-paste( quantiles_with_pt.estimate*100,"%",sep="")
  for (i in 1:length(totci)){
    if (names(totci)[i]=="50%"){
      names(totci)[i]<-"point est."}
  }
  return(list(estimate=totci, MSEfit = zMSE))
}

#' Fit a specified model to multiple systems estimation data
#'
#' This routine fits a specified model to multiple systems estimation data, taking account of the possibility of empty overlaps between
#'  pairs of observed lists.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.

#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list parameters to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list parameters are included. If \code{mX = NULL}, no two-list parameters are included and
#'  the main effects model is fitted.
#' If only one two-list parameter is to be fitted, it may be specified as a vector of length 2, e.g \code{mX=c(1,3)}
#' for the parameter corresponding to lists 1 and 3.
#'
#' @param check If \code{check = TRUE} check first of all if the maximum likelihood
#' estimate exists and is identifiable, using the routine \code{\link{checkident}}.  If either condition fails, print an appropriate error message
#' and return the error code.
#'
#'
#' @return A list with components as below
#'
#' \code{fit} Details of the fit of the specified model as output by \code{glm}.  The Akaike information criterion is adjusted to take account
#'   of the number of parameters corresponding to non-overlapping pairs.
#'
#' \code{emptyoverlaps} Matrix with two rows, giving the list pairs within the model for which no cases are observed in common.
#'   Each column gives the indices of a pair of lists, with the names of the lists in the column name.
#'
#' \code{poisspempty} the Poisson p-values of the parameters corresponding to non-overlapping pairs.
#'
#' @examples
#' data(NewOrl)
#' modelfit(NewOrl,mX= c(1,3), check=TRUE)
#'
#' @importFrom stats glm poisson
#'
#' @export

modelfit = function (zdat, mX = NULL, check = TRUE){
  if (check) {
    zcheck = checkident(zdat, mX, verbose = TRUE)
    if (zcheck$ierr > 0)
      return(zcheck$ierr)
  }
  zz = buildmodel(zdat, mX)
  eo = zz$emptyoverlaps
  nover = dim(eo)[2]
  fit = glm(zz$modelform, family = poisson, data = zz$datamatrix, x = TRUE)
  fit$aic = fit$aic + 2 * nover
  #
  # find the probabilities of each observed empty overlap in the model being empty, given a model excluding that parameter.
  probzero = NULL
  if (nover > 0) {
    probzero= vector("numeric", nover)
    m = dim(zdat)[2]-1
    if (sum(mX) == 0) {
      mX = t(expand.grid(1:m, 1:m))
      mX = mX[, mX[1, ] < mX[2, ]]
    }
    if (is.vector(mX)) mX= as.matrix(mX, ncol=1)
    names(probzero) = dimnames(eo)[[2]]
    for (iover in (1:nover)) {
      # for each empty overlap, remove that parameter from the model and fit the model without it
      jeocol = (1:dim(mX)[2]) [ mX[1,]==eo[1,iover] & mX[2,]==eo[2,iover]]
      mX1 = mX[,-jeocol]
      if (length(mX1)==0) mX1=NULL
      zz1 = buildmodel(zdat, mX1)
      zfit = glm(zz1$modelform, family=poisson, data=zz1$datamatrix, x=TRUE)
      #  now calculate the estimated probability that the given overlap is empty, summing over all relevant cells
      pred0 = predict(zfit)
      dat0 = zfit$data
      lambdastar = sum(exp(pred0) * dat0[, eo[1,iover]] * dat0[, eo[2,iover]])
      probzero[iover] = exp(-lambdastar)
    }
  }
  return(list(fit = fit, emptyoverlaps = eo, poisspempty = probzero))
}

#'Stepwise fit using Poisson p-values.
#'
#'Starting with a model with main effects only, two-list parameters are added one by one.
#'At each stage the parameter with the lowest p-value is added, provided that p-value is lower than \code{pthresh}, and provided that the resulting
#' model does not fail either of the tests in \code{\link{checkident}}.
#'
#' For each candidate two-list parameter for possible addition to the model, the p-value is calculated as follows.
#' The total of cases occurring on both lists indexed by the parameter (regardless of
#' whether or not they are on any other lists) is calculated.
#' On the null hypothesis that the effect is not included in the model, this statistic has a Poisson
#' distribution whose mean depends on the parameters within the model.    The one-sided Poisson p-value of the observed statistic is calculated.
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#'@param pthresh this is the threshold below which the p-value of the newly added parameter needs to be in order to be included in the model.
#'  If \code{pthresh = 0} then the model with main effects only is returned.
#'
#'@return A list with components as below
#'
#'\code{fit} Details of the fit of the specified model as output by \code{glm}.  The Akaike information criterion is adjusted to take account
#'   of the number of parameters corresponding to non-overlapping pairs.
#'
#' \code{emptyoverlaps} Matrix with two rows, giving the list pairs within the model for which no cases are observed in common.
#'   Each column gives the indices of a pair of lists, with the names of the lists in the column name.
#'
#' \code{poisspempty} the Poisson p-values of the parameters corresponding to non-overlapping pairs.
#'
#'@examples
#'data(NewOrl)
#'stepwisefit(NewOrl, pthresh=0.02)
#'
#'@importFrom stats ppois predict
#'
#'
#'@export

stepwisefit<- function(zdat, pthresh=0.02) {
  m = dim(zdat)[2] - 1
  ierrfullmodel=NA
  if (pthresh == 1){
    ierrfullmodel = checkident(zdat, mX=0)}
  if ((pthresh==1) & (ierrfullmodel==0)) {
    zfit=modelfit(zdat,mX=0, check=FALSE)
  } else if(pthresh == 0){
    zfit=modelfit(zdat,mX=NULL, check=FALSE)}
  else {
    # ------  Set up a list of two-list parameters to be included in the mode
    ints = NULL
    mX = NULL
    zfit = modelfit(zdat, mX=NULL, check=FALSE)
    for (i in (1:(m - 1))) for (j in ((i + 1):m)) {
      ints = cbind(ints,c(i, j))
    }
    nints = dim(ints)[2]
    nstar = rep(0,nints)
    for (j in (1:nints)) {
      nstar[j] = sum( zdat[,ints[1,j]]*zdat[,ints[2,j]]*zdat[,m+1])
    }
    #--------cycle through and find best parameter to add.
    intsincluded = rep(FALSE, nints)
    for (jcycle in (1:nints)) {
      pred0 = predict(zfit$fit)
      dat0 = zfit$fit$data
      pvalnew = rep(1, nints)
      for (j in (1:nints)) {
        if (!intsincluded[j]) {
          #----- find the predicted Poisson intensity of this overlap
          pstar = sum(exp(pred0) * dat0[,ints[1,j]]*dat0[,ints[2,j]] )
          pvalnew[j] = min( ppois(nstar[j],pstar), ppois(nstar[j]-1, pstar, FALSE) )
          #----- test the potential model and set the p-value to 1 if the model
          #---      would yield an infinite result or be non-identifiable
          ierr = checkident(zdat, cbind(mX, ints[,j]))
          if (ierr > 0) pvalnew[j] = 1
        }
      }
      #---------if none of them have p value below threshold, finish, otherwise incorporate and proceed
      pvmin = min(pvalnew)
      if (pvmin >= pthresh) return(zfit)
      jintmax = min((1:nints)[pvalnew == pvmin])
      mX1 = cbind(mX, ints[, jintmax])
      zf1 = modelfit(zdat, mX = mX1, check=FALSE)
      zfit = zf1
      mX = mX1
      intsincluded[jintmax] = TRUE
    }}
  return(zfit)
}
#' Plot of simulation study
#'
#' This routine reproduces Figure 1 of Chan, Silverman and Vincent (2021).
#'
#'   Simulations are carried out for two different three-list models.
#'   In one model, the probabilities of capture are 0.01, 0.04 and 0.2 for the three lists respectively, while in the other the probability
#'   is 0.3 on all three lists.  In both cases,
#'   captures on the lists occur independently of each other,
#'   so there are no two-list effects.
#'   The first model is chosen to be somewhat more typical of the sparse capture
#'   case, of the kind which often occurs in the human trafficking context,
#'   while the second, more reminiscent of a classical mark-recapture study.
#'
#'   The probability of an individual having each possible capture history is first evaluated.
#'   Then these probabilities are multiplied by \code{Nsamp = 1000} and, for each simulation
#'   replicate, Poisson random values with expectations equal to these values are generated
#'   to give a full set of observed capture histories;
#'   together with the null capture history the expected number of counts
#'   (population size) is equal to \code{Nsamp}.
#' Inference was carried out both for the model with main effects only, and for the model with the addition of
#'   a correlation effect between the first two lists.
#' The reduction in deviance between the two models was determined. Ten thousand simulations were carried out.
#'
#' Checking for compliance with the conditions for existence and identifiability of the
#' estimates shows that a very small number of the simulations for the sparse model (two out of
#' ten thousand) fail the checks for existence even within the extended maximum likelihood context.
#' Detailed investigation shows that in neither of these cases is the dark figure itself not estimable;
#' although the parameters themselves cannot all be estimated, there is a maximum likelihood estimate
#' of the expected capture frequencies, and hence the deviance can still be calculated.
#'
#' The routine produces QQ-plots
#'   of the resulting deviance reductions against quantiles of the \eqn{\chi^2_1} distribution,
#'   for \code{nsim} simulation replications.
#'
#' @param nsim The number of simulation replications
#' @param Nsamp The expected value of the total population size within each simulation
#' @param seed  The random number seed
#'
#' @return
#'
#'  An \code{nsim} \eqn{\times} 2 matrix giving the changes in deviance for each replication for each of the two models.
#'
#'
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#'@importFrom stats qchisq rpois ppoints
#'@importFrom graphics matplot abline
#'
#' @export
investigateAIC <- function(nsim=10000, Nsamp= 1000, seed = 1001) {
  #
  # ---- define functions to be called
  #
  simulate_deviances <- function(N, plist=0.1, nsim=100, seed= NULL ) {
    if (length(plist)==1) plist= rep(plist, 3)
    set.seed(seed)
    devdiffs = rep(NA, nsim)
    for (j in (1:nsim)) {
      zdat1 = simulate_capture(N, plist)
      devdiffs[j] = deviance_compare(zdat1)
    }
    return(devdiffs)
  }
  #
  #
  simulate_capture <- function(N=1000, plist=rep(0.1,3)) {
    zdata = rbind(diag(3), 1-diag(3), rep(1,3))
    meancounts = N* exp( zdata %*% log(plist) + (1-zdata)%*% log(1-plist))
    simcounts = rpois(7, meancounts)
    zdata = cbind(zdata,simcounts)
    dimnames(zdata)[[2]] = c("A", "B", "C", "count")
    return(zdata)
  }
  #
  #
  deviance_compare <- function(zdat) {
    if(sum(zdat[4:7,4])==0) return(0)
    devnull =  modelfit(zdat)$fit$deviance
    devalt =   modelfit(zdat, mX=c(1,2),check=FALSE)$fit$deviance
    devdiff = devnull - devalt
    return(devdiff)
  }
  #
  # ----carry out the simulations and plot the results
  #
  zdev0 = simulate_deviances(Nsamp, plist=c(0.01,0.04,0.2),nsim, seed=seed )
  zdev1 = simulate_deviances(Nsamp, plist=0.3,nsim, seed=seed )
  zdev= cbind(sort(zdev0), sort(zdev1))
  zchisq = qchisq(ppoints(nsim), df=1)
  matplot(zchisq, zdev, type="l", xlim=c(0,8), ylim=c(0,8), xaxs="i", yaxs="i",
          xlab="Chi-Squared distribution", ylab="Simulated deviances", main="QQ plot")
  abline( 0,1, lty=3)
  return(zdev)
}

#' Check all possible models for existence and identifiability
#'
#'This routine efficiently checks whether every possible model
#'  satisfies the conditions for the existence and identifiability of an extended maximum likelihood estimate.
#'   For identifiability it
#' is only necessary to check the model containing all two-list effects.   For existence, the approach set out in Chan, Silverman and Vincent (2021)
#' is used. This uses the linear programming approach described in a more general and abstract context by Fienberg and Rinaldo (2012).
#'
#'
#' If the extended maximum likelihood estimator exists for a particular model, then it will still exist if one or more overlapping pairs are removed from the model.   This allows the search to be carried out efficiently, as follows:
#'
#' 1. Search over `top-level' models which contain all overlapping pairs.  The number of such models is \eqn{2^k}, where \eqn{k} is the number of non-overlapping pairs, because there are \eqn{2^k} possible choices of the set of non-overlapping pairs to include in the model.
#'
#' 2. For each top-level model, if the estimate exists there is no need to consider that model further.  If the estimate does not exist, then a branch search is carried out over the tree structure obtained by successively leaving out overlapping pairs.
#'
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#'@param nreport A message is printed for every \code{nreport} models checked.
#'  It gives the number of top level models to be checked altogether, and also the number of models found so far for which the estimate does not exist.
#'
#'@return
#'
#' The routine prints a message if the model with all two-list parameters is not identifiable.   As set out above, it gives regular progress reports in
#'  the case where there are a large number of models to be checked.
#'
#' If all models give estimates which exist, then a message is printed to that effect.
#'
#' If there are models for which the estimate does not exist, the routine reports the number of such models and returns a
#'   matrix whose rows give the parameters included in them.
#'
#'@references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
#' models. \emph{Ann. Statist.} \strong{40}, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
#'
#' @examples
#'data(Artificial_3)
#'data(Western)
#'checkallmodels(Artificial_3)
#'checkallmodels(Western)
#'
#'@export

checkallmodels <-function (zdat, nreport= 1024) {
  m = dim(zdat)[2] - 1
  if (count_triples(zdat) ==0) cat("The model including all two-list parameters is not identifiable. \n")
  varnames = dimnames(zdat)[[2]][1:m]
  zdfull = tidylists(zdat, includezerocounts = TRUE)
  Omegamat = zdfull[, 1:m]
  Countobserved = zdfull[, m+1]
  #--------construct parameter set as capture histories
  Thetamat = rbind(rep(0,m), diag(m), matrix(0, nrow=m*(m-1)/2, ncol=m))
  dimnames(Thetamat)[[1]] = rep("", 1 +m + m*(m-1)/2)
  dimnames(Thetamat)[[1]][2:(m+1)] =  varnames
  n1 = m+2
  for ( i1 in (1:(m-1))) for (i2 in ((i1+1):m)) {
    Thetamat[n1,i1] = 1
    Thetamat[n1,i2] = 1
    dimnames(Thetamat)[[1]][n1] = paste(varnames[c(i1,i2)], collapse=":")
    n1 = n1+1}
  #---- construct full A^T matrix
  Afull = matrix(0, nrow = dim(Thetamat)[1], ncol= dim(Omegamat)[1],
                 dimnames= list(dimnames(Thetamat)[[1]], NULL) )
  for (i in (1:dim(Thetamat)[1])) for (j in (1:dim(Omegamat)[1])) {
    Afull[i,j] = all(Thetamat[i,] <= Omegamat[j,])
  }
  #---construct Nstar vector
  Nstar = Afull%*%Countobserved
  npairsempty = sum(Nstar==0)
  #---find parameters corresponding to non-overlapping pairs and keep the rest
  parkeep= (Nstar > 0)
  Athetadag = Afull[parkeep, ]
  Nstardag = Nstar[parkeep]
  if (npairsempty > 0) {
    Athetanonoverlap = Afull[Nstar==0,]
    if (npairsempty == 1)  Athetanonoverlap = as.matrix(Athetanonoverlap, nrow = 1)
    Namesnonoverlap = dimnames(Thetamat)[[1]][Nstar==0]
  }
  Namesoverlap = dimnames(Thetamat)[[1]][parkeep][-(1:(m+1))]
  npairsnonempty = m*(m-1)/2 - npairsempty
  #----construct each possible subset of rows of Athetanonoverlap and remove corresponding columns
  #----- of Athetadag
  nmodels = 2^npairsempty
  ierr = rep(0, nmodels)
  errlist=list()
  #----formulate the LP for every possible model
  for (j in (1:nmodels)) {
    if (npairsempty==0) jkeep = 0 else jkeep = (intToBits(j)[1:npairsempty] == 1)
    if ( sum(jkeep) == 0 ) {
      amat=Athetadag
    } else {
      if (sum(jkeep) == 1) {
        Omegakeep = (Athetanonoverlap[jkeep, ] == 0)
      } else {
        Omegakeep = (apply( Athetanonoverlap[jkeep,],2,sum) ==0 )
      }
      amat = Athetadag[, Omegakeep]
    }
    #---set up LP
    if ( checkthetasubset(NULL, amat, Nstardag,m) ) {
      if (npairsempty==0) errmodel0=NULL else errmodel0 = Namesnonoverlap[jkeep]
      errmodel = c(errmodel0, Namesoverlap)
      errlist[[1+length(errlist)]] = errmodel
      errlist1 = subsetsearch(npairsnonempty, checkthetasubset,  testnull = FALSE, amat=amat ,tvec=Nstardag , nlists=m)
      if (length(errlist1)>0) for (k in (1:length(errlist1))) { jkeep1 = errlist1[[k]]
      errmodel = c( errmodel0, Namesoverlap[-jkeep1])
      errlist[[1+length(errlist)]] = errmodel}
    }
    #if (ierr[j]>0) errmat = cbind(errmat, jkeep)
    if (j%%nreport == 0) cat(j, " of ", nmodels, " models.  Sum of errors so far ",length(errlist), "\n")
  }
  if (length(errlist)==0) {
    cat("The extended maximum likelihood estimate exists for all models.", "\n")
    invisible()
  } else { cat(length(errlist), " models found for which the extended maximum likelihood estimate does not exist.","\n")
    errmat = matrix(0, nrow=length(errlist), ncol= m*(m-1)/2, dimnames= list(NULL, dimnames(Thetamat)[[1]][-(1:(m+1))]))
    for (j in (1:length(errlist))) {
      errmat[j,errlist[[j]]] = 1
    }
    return(ordercaptures(errmat))
  }
}


#' Search subsets for a property which is inherited in a particular way
#'
#' The following routine returns a list of all subsets of a given set for which a specified property is \code{TRUE}.
#'  It is assumed that if the property is \code{FALSE} for any particular subset, it is also \code{FALSE} for all supersets of that subset.
#'  This enables a branch search strategy to be used to obviate the need to search over supersets of subsets already eliminated from consideration.
#'  It is used within the hierarchical search step of \code{\link{checkallmodels}}.
#'
#'@param n an integer such that the search is over all subsets of \eqn{\{1,...,n\}}
#'
#'@param checkfun a function which takes arguments \code{zset}, a subset of \eqn{\{1,...,n\}} and ...
#'The function returns the value \code{TRUE} or \code{FALSE}. It needs to have the property that if it is \code{FALSE} for any
#'particular \code{zset}, it is also \code{FALSE} for all supersets of \code{zset}.
#'
#'@param testnull If \code{TRUE}, then include the null subset in the search. If \code{FALSE}, do not test the
#'null subset.
#'
#'@param ... other arguments
#'
#'@return A list of all subsets for which the value is \code{TRUE}
#'
#'@references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#'@export
subsetsearch <- function(n, checkfun, testnull = TRUE, ...) {

  #--- test the null set first
  truelist=list()
  if (n==0) return(truelist)
  if (testnull) {
    if (!checkfun(NULL, ...)) {return(truelist)} else {truelist = list(vector("numeric"))}
  }
  #--- set up the initial candidates
  candidates = as.list(1:n)
  ntested = 0
  #--- as long as there are untested candidates, test them and act appropriately
  while (ntested < length(candidates))  {
    ntested = ntested + 1
    zs = candidates[[ntested]]
    if (checkfun(zs, ...))  {
      #---- if true, add zs to truelist and add new candidates obtained by adding one member to zs
      truelist[[1 + length(truelist)]] = zs
      if (max(zs) < n) {
        for (j in ((1+max(zs)):n)) {
          zs1 = c(zs, j)
          candidates[[ 1 + length(candidates)]] = zs1
        }}}}
  return(truelist)
}

#'
#' Check a subset of the parameter set theta
#'
#' This routine leaves out a particular set of parameters corresponding to the two-list effects from the parameter set theta.
#' For the resulting model, it constructs the linear programming
#' problem to check whether the extended maximum likelihood estimates of the parameters exists. It is called internally by \code{checkallmodels}.
#'
#' @param zset set of indices that is not included, numbered among the two-list effects only
#'
#' @param amat a design matrix
#'
#' @param tvec vector of sufficient statistics
#'
#' @param nlists number of lists in the original capture-recapture matrix
#'
#' @return If the return result is \code{TRUE}, the linear program shows that the extended maximum likelihood estimate does not exist.
#' If the return result is \code{FALSE}, the estimate exists.
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#'@export
checkthetasubset <- function(zset, amat, tvec, nlists) {

  #---construct the subsetted a matrix and t vector
  if (length(zset) > 0) {
    zset = 1+nlists+zset
    amat = amat[-zset,]
    tvec = tvec[-zset]
  }
  npar = dim(amat)[1]
  nobs = dim(amat)[2]
  #---construct and solve the LP problem
  f.obj = c(rep(0, nobs), 1)
  const.rhs = c(tvec, rep(0, nobs))
  const.mat = rbind(cbind(amat, rep(0, npar)), cbind(diag(nobs), rep(-1, nobs) ) )
  const.dir = c(rep("=", npar), rep(">=", nobs))
  zlp = lpSolve::lp("max", f.obj, const.mat, const.dir, const.rhs)
  return((zlp$objval==0))
}

#' Order capture histories
#'
#' Given a matrix with capture histories only, the routine orders the capture histories first by the number of
#' 1s in the capture history and then lexicographically by columns.
#'
#' @param zmat Data matrix with \eqn{t} columns. The \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories observed. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @return A data matrix that is ordered first by the number of 1s in the capture history
#' and then lexicographically by columns.
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#'
#' @examples
#' #UK 5 list data
#' data(UKdat_5)
#' # Obtain matrix with capture histories only
#' UKdat_5_cap <- UKdat_5[, 1:dim(UKdat_5)[2]-1]
#' ordercaptures(UKdat_5_cap)
#
#' @export

ordercaptures <- function(zmat){
  #  Given a matrix of capture histories order them first by the number of 1s and then lexicographically by columns
  zz = paste(rep(LETTERS[1:26], each=26), rep(LETTERS[1:26], times=26), sep="")
  zza = zz[1: dim(zmat)[2]]
  zo = character(dim(zmat)[1])
  for (j in (1:dim(zmat)[1])) {
    zo[j] = paste(zz[1+ sum(zmat[j,]) ], zza[(zmat[j,]==1)], collapse=":")
  }
  zmat = zmat[order(zo),]
  return(zmat)
}


#' Count number of triples of overlapping lists
#'
#' The routine counts the number of subsets of size three of lists such that every pair of
#' lists in the triple overlaps. If the number is zero, then the model with all two-list effects
#' is unidentifiable.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @return a count of subsets of size three of lists such that every pair of lists in the triple overlaps.
#'
#' @examples
#' data(Western)
#' data(Artificial_3)
#' count_triples(Western)
#' count_triples(Artificial_3)
#'
#'@export

count_triples <- function(zdat) {

  m = dim(zdat)[2] - 1
  zdat = as.matrix(zdat[zdat[,m+1]>0, -(m+1)])
  zadj = crossprod(zdat)
  zadj[zadj > 0] = 1
  diag(zadj) = 0
  ntriples = sum(zadj*crossprod(zadj))/6
  return(ntriples)
}

#' Bootstrapping to evaluate confidence intervals using BCa methods
#'
#' This routine implements the bootstrapping and jacknife approach as detailed in Section 3.3 of Chan, Silverman and Vincent (2021).
#'  It calls the routine \code{\link{estimatepopulation.0}} and so is the preferred routine to be called if a user wishes to estimate the
#'   population and obtain BCa confidence intervals.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param nboot Number of bootstrap replications.
#'
#' @param pthresh p-value threshold used in \code{\link{estimatepopulation.0}}.
#'
#' @param iseed seed for initialisation.
#'
#' @param alpha Bootstrap quantiles of interests.
#'
#' @param ... other arguments which will be passed to \code{\link{estimatepopulation.0}}
#'
#' @return A list with components as below:
#'
#'\code{popest}  point estimate of the total population of the original data set
#'
#'\code{MSEfit} model fitted to the data, in the format described in \code{\link{modelfit}}
#'
#'\code{bootreps} point estimates of total population sizes from each bootstrap sample
#'
#'\code{ahat} the estimated acceleration factor
#'
#'\code{BCaquantiles} BCa confidence intervals
#'
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}, \strong{40(3)}, 189-228.
#'
#' @examples
#' data(Korea)
#' zdat=Korea
#' estimatepopulation(zdat,nboot=10)
#'
#' @export
#'

estimatepopulation <- function(zdat, nboot=1000, pthresh=0.02, iseed=1234, alpha= c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975), ...) {
  #  find nboot bootstrap estimates of population size
  set.seed(iseed)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  nobs=sum(countsobserved)
  bootreps = rep(NA, nboot)
  #   find point estimate and corresponding model using given pthresh
  populationestimatefromdata = estimatepopulation.0(zdat, quantiles=NULL, pthresh=pthresh, ...)
  popest = populationestimatefromdata$estimate
  MSEfit = populationestimatefromdata$MSEfit
  #   set up bootstrap model
  #   generate nboot multinomial samples.
  #   Then use a multinomial distribution to find the actual realization.  Then set out synthetic data.
  for (j in (1:nboot)) {
    counts = rmultinom(1,nobs,countsobserved)
    zdatboot = cbind(zdat[,-n2], counts)
    #   for each sample, find estimated total population size
    bootreps[j] = estimatepopulation.0(zdatboot, quantiles=NULL, pthresh=pthresh,...)$estimate
  }
  # use jackknife to find acceleration factor
  jackest = rep(0, n1)
  for (j in (1:n1)) {
    nj = zdat[j,n2]
    if (nj > 0) {
      zd1 = zdat
      zd1[j,n2] = nj - 1
      jackest[j] = estimatepopulation.0(zd1,quantiles=NULL, pthresh=pthresh, ...)$estimate
    }
  }
  jr = sum(countsobserved*jackest)/sum(countsobserved) - jackest
  # find estimated acceleration factor by counting each residual the number of times it would occur,
  #   via the count of the corresponding capture history
  ahat = sum(countsobserved*jr^3)/(6 * (sum(countsobserved*jr^2))^{3/2})
  #  find BCa conf limits
  confquantiles = bcaconfvalues(bootreps, popest,ahat, alpha)
  return(list(popest=popest, MSEfit=MSEfit, bootreps=bootreps, ahat=ahat, BCaquantiles=confquantiles))
}

#' BCa confidence intervals
#'
#' The BCa confidence intervals use percentiles of the bootstrap distribution of the population size
#' , but adjust the percentile actually used. The adjusted percentiles depend on an
#' estimated bias parameter, and the quantile function of the estimated bias parameter is the proportion
#' of the bootstrap estimates that fall below the estimate from the original data, and an
#' estimated acceleration factor, which derivation depends on a jackknife approach. This routine is called internally
#' by \code{estimatepopulation}.
#'
#'
#' @param bootreps Point estimates of total population sizes from each bootstrap sample.
#'
#' @param popest A point estimate of the total population of the original data set.
#'
#'@param ahat the estimated acceleration factor
#'
#'@param alpha Bootstrap quantiles of interests
#'
#'@return BCa confidence intervals
#'
#'@references
#' Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}, \strong{40(3)}, 189-228.
#'
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. \emph{Journal of the American Statistical Association}, \strong{82(397)}, 171-185.
#' @export
bcaconfvalues<-function(bootreps, popest, ahat, alpha=c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975) ) {
  # find BCA critical values
  z0 = qnorm(mean(bootreps < popest, na.rm=TRUE))
  za = qnorm(alpha)
  za0 = za+z0
  zq = pnorm( z0 + za0/(1-ahat*za0))
  bcac = quantile(bootreps, probs=zq, type=8, na.rm = TRUE)
  names(bcac) = alpha
  return(bcac)
}

#'Comparison of BIC approach and BCa approach
#'
#' This routine carries out the simulation study as detailed in Section 3.4 of Chan, Silverman and Vincent (2021).
#'  If the original data set has low counts, so that there is a possibility of a simulated data set containing empty lists, then it
#'  may be advisable to use the option \code{noninformativelist=TRUE}.
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has zero count.
#'
#' @param nsims Number of simulations to be carried out.
#'
#' @param nboot Number of bootstrap replications for each simulation
#'
#' @param pthresh p-value threshold used in \code{\link{estimatepopulation.0}}.
#'
#' @param iseed seed for initialization.
#'
#' @param alpha bootstrap quantiles of interests.
#'
#' @param noninformativelist  if \code{noninformativelist=TRUE} then each generated data set in the simulation study (including all bootstrap replications)
#'     will be passed to \code{\link{removenoninformativelists}}.
#'
#' @param verbose If \code{verbose=FALSE}, then the progress of the simulation will not show.
#' If \code{verbose=TRUE}, then the progress of the simulation will be shown.
#'
#'@param ... other arguments.
#'
#'@return A list with components as below
#'
#'\code{popest} Total population point estimate from the original data using
#'\code{estimatepopulation.0} with default threshold.
#'
#'\code{BICmodels} The best model chosen by the BIC at each simulation.
#'
#'\code{BICvals} Point estimates of the total population and standard error of the best model chosen by the BIC at each simulation.
#'
#'\code{simreps} Counts associated to each capture history at each simulation.
#'
#'\code{modelmat} A full capture history matrix excluding the row corresponding to the dark figure.
#'
#'\code{popestsim} Total population estimate given by the BCa method in each simulation.
#'
#'\code{BCaquantiles} bootstrap confidence intervals given by the BCa method.
#'
#'\code{BICconf} confidence interval given by the BIC method.
#'
#'@references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}, \strong{40(3)}, 189-228.
#'
#' Rivest, L-P. and Baillargeon, S. (2014) Rcapture. CRAN package. Available from Available from \url{https://CRAN.R-project.org/package=Rcapture}.
#'
#'
#'@importFrom stats pnorm quantile rbinom rmultinom
#'
#'@examples
#'zdat=UKdat_5
#'BICandbootstrapsim(zdat,nsims=1000, nboot=100, pthresh=0.02, iseed=1234, noninformativelist=T)
#'@export

BICandbootstrapsim <- function(zdat, nsims=100,nboot=100,pthresh=0.02, iseed=1234, alpha=c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975),noninformativelist=F,verbose=F, ...){
  #  simulate data
  set.seed(iseed)
  nobserved = sum(zdat[, dim(zdat)[2]])
  BICmods = vector("character", length=nsims)
  BICvals = matrix(nrow=nsims, ncol=2)
  BICconf = matrix(nrow=nsims, ncol=4, dimnames=list(NULL, c(0.025,0.1,0.9,0.975)))
  #   find point estimate using given pthresh
  fullestimate = estimatepopulation.0(zdat, quantiles=NULL, pthresh=pthresh, ...)
  #   find estimated cell frequencies for fitted model
  popest = fullestimate$estimate
  #   set up bootstrap model
  predfreq = exp(predict(fullestimate$MSEfit$fit))
  modelmat = fullestimate$MSEfit$fit$model[,-1]
  simreps = matrix(NA, nrow = dim(modelmat)[1], ncol=nsims)
  popestr = round(popest)
  #   generate nsims multinomial samples.  First use a binomial distribution to find the number of observations excluding the dark figure.
  #   Then use a multinomial distribution to find the actual realization.  Then set out synthetic data.
  for (j in (1:nsims)) {
    nobs = rbinom(1,popestr,sum(predfreq)/popest)
    simreps[,j] = rmultinom(1,nobs,predfreq)
    # find the BIC estimates at the same time
    zdatsim = cbind(modelmat, simreps[,j])
    if (noninformativelist) zdatsim=removenoninformativelists(zdatsim)
    zallres=Rcapture::closedpMS.t(zdatsim, dfreq=T,maxorder=2)
    zr = zallres$results
    indmin = (1:dim(zr)[1])[zr[,7]==min(zr[,7])]
    BICmods[j] = dimnames(zr)[[1]][indmin]
    BICvals[j,] = zr[indmin, 1:2]
    #  find multinomial 95% and 80% confidence intervals--probably not the most elegant way but it works
    BICconf[j,c(2,3)] = Rcapture::closedpCI.t(zdatsim, dfreq=T,mX=BICmods[j], alpha=0.2)$CI[2:3]
    BICconf[j,c(1,4)] = Rcapture::closedpCI.t(zdatsim, dfreq=T,mX=BICmods[j])$CI[2:3]
    if (verbose) cat(j)
  }
  #  now do stepwise estimates and bootstrap CIs
  popests = rep(0,nsims)
  bootCIs = matrix(nrow=nsims, ncol=length(alpha))
  dimnames(bootCIs)[[2]] = alpha
  for (j in (1:nsims)) {
    count = simreps[,j]
    zdatsim = cbind(modelmat, count)
    if (noninformativelist) zdatsim=removenoninformativelists(zdatsim)
    bootoutput = estimatepopulation(zdatsim, nboot=nboot, alpha=alpha, pthresh=pthresh)
    popests[j]=bootoutput$popest
    bootCIs[j,] = bootoutput$BCaquantiles
    if (verbose)  cat(j)
  }
  return(list(popest=popest, BICmodels = BICmods, BICvals=BICvals, simreps=simreps, modelmat=modelmat, popestsim=popests, BCaquantiles=bootCIs, BICconf=BICconf))
}

#' Remove non-informative list
#'
#' The routine cleans up the data set by removing any lists that contain no data, any lists which contain all the observed data, and
#'  any list whose results duplicate those of another list.
#'   If as a result the original data set has no list left, it returns a matrix with the value of the total count.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#'
#' @return data matrix that contains no duplicate lists, no lists with no data, and no lists that contain all the observed data.
#' If all lists are removed, the total count is returned.
#'
#'@references
#'Chan, L., Silverman, B. W., and Vincent, K. (2021).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists.
#' \emph{Journal of American Statistcal Association}, \strong{116(535)}, 1297-1306,
#' Available from \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748}.
#'
#' @export
removenoninformativelists<-function(zdat){
  zdat=as.matrix(zdat)
  #remove duplicate list
  zdat=unique(zdat,MARGIN=2)
  #remove list that contain no or all data
  m2=dim(zdat)[2]
  countname = dimnames(zdat)[[2]][m2]
  ltot = t( zdat[, -m2])%*% zdat[,m2]
  mtot= sum(zdat[,m2])
  jkeep = (ltot > 0) & (ltot < mtot)
  if (sum(jkeep) > 0 ) zdat = zdat[,c(jkeep, T)] else {
    #  if there are no lists left at all, just return the total count as a matrix
    zdat = matrix(mtot, nrow=1, ncol=1, dimnames= list(NULL, countname))
  }
  return(zdat)
}
#' Decode capture history
#'
#' Given a capture history as a number and the number of lists, decode it into a logical vector giving
#' presence or absence in the capture history.
#'
#' @param k The capture history to be decoded
#' @param nlists The number of lists
#' @return A logical vector of length \code{nlists} giving presence or absence in the capture history
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' decode_capture(2,5)
#' decode_capture(1,4)
#'
#' @export
decode_capture = function(k, nlists) {
  z = as.logical(intToBits(k-1))[1:nlists]
  return(z)
}
#' Find the "descendants" of a given capture history
#'
#' Given any encoded capture history, find all the encoded capture histories that include the original capture history and any other lists
#'
#' @param k An encoded capture history
#' @param nlists The total number of lists
#' @param omitk Determine whether the original capture history is included as a descendant of itself. If \code{omitk=T} it is not.
#' @return a vector giving the encoded versions of the descendants
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' descendants(2,5)
#' descendants(5,10)
#'
#' @export
descendants = function(k,nlists, omitk = FALSE) {
  if (length(k)==0) return(numeric(0))
  z = decode_capture(k, nlists)
  jz = sum(z)
  kdesc = as.vector(k)
  if (omitk) kd1 = numeric(0) else kd1 = kdesc
  for (i in ((jz+1):nlists)) {
    kd2 = NULL
    for (j in (1:length(kdesc))) kd2 = c(kd2, child_captures(kdesc[j], nlists))
    kdesc = unique(kd2)
    kd1 = c(kd1, kdesc) }
  return(sort(kd1))
}
#' Find the "ancestors" of a given capture history
#'
#' Given any encoded capture history and the number of lists, find all the encoded capture histories that are included in the original capture history
#'
#' @param k An encoded capture history
#' @param nlists The total number of lists
#'
#' @return a vector giving the encoded versions of the ancestors
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' ancestors(2,10)
#' ancestors(1,5)
#'
#' @export
ancestors = function(k,nlists=10) {
  z = decode_capture(k, nlists)
  jz = sum(z)
  kanc = as.vector(k)
  kd1 = kanc
  for (i in (1:nlists)) {
    kd2 = NULL
    for (j in (1:length(kanc))) kd2 = c(kd2, parent_captures(kanc[j], nlists))
    kanc = unique(kd2)
    kd1 = c(kd1, kanc) }
  return(sort(kd1))
}
#' Find the "parents" of a given capture history
#'
#' Given any encoded capture history and the number of lists, find the encoded capture histories which are obtained by leaving out just one list in turn
#'
#' @param k An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#' @return a vector giving the encoded versions of the parents
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' parent_captures(2,10)
#' parent_captures(1,4)
#'
#' @export
parent_captures = function(k, nlists=10) {
  z = decode_capture(k, nlists)
  kd = 2^{(0:(nlists-1))}[z]
  return(k - kd)
}
#' Find the "children" of a given capture history
#'
#' Given any encoded capture history that corresponds to the row number of the capture history data set and the number of lists, find the encoded capture histories which are obtained by adding one more list in turn
#'
#' @param k An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#' @return a vector giving the encoded versions of the children
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' child_captures(2,5)
#' child_captures(1,4)
#'
#' @export
child_captures = function(k, nlists) {
  z = decode_capture(k, nlists)
  kd = 2^{(0:(nlists-1))}[!z]
  return(k + kd)
}
#' Set up the inclusion matrix for all possible capture histories
#'
#' This is the master design matrix which maps parameters to observations.
#' Rows correspond to observations and columns to parameters.
#'
#'@param nlists The number of lists
#'
#'@return A matrix whose  \eqn{(i,j)} element is 1 if the expected log of observation \eqn{i} depends on parameter \eqn{j},
#' in other words if \eqn{j} is an ancestor of \eqn{i}.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' #Create master design matrix with 5 lists
#' make_master_design(5)
#' #Create master design matrix with 3 lists
#' make_master_design(3)
#'
#'@export
make_master_design = function(nlists) {
  # make design matrix where rows correspond to observations
  #  and columns to parameters
  ncaps = 2^nlists
  xdes = matrix(0, nrow=ncaps, ncol=ncaps, dimnames= list( 1:ncaps, 1:ncaps))
  for (i in (1:ncaps)) {
    ipars = ancestors(i, nlists)
    xdes[i,ipars] = 1
  }
  xdes=xdes[-1,]
  return(xdes)
}
#' Preliminary processing of a data matrix
#'
#' Perform various preprocessing tasks on the data
#'
#' @param xdat Data matrix of the usual kind
#'
#' @return A list with the following elements
#' \describe{
#' \item{nobs}{Numbers of observations indexed by encoded histories}
#' \item{nstar}{For each capture history, total number of observations for that capture history and all its descendants}
#'  \item{nlists}{Total number of lists}
#'  \item{listnames}{Names of the lists, constructed to be A, B, ... if necessary}
#'  \item{data}{The input data matrix}
#'  \item{notestimable}{A vector indicating which parameters are not estimable, because they are strict descendants of parameters
#'  which would be already estimated to be \eqn{-\infty} if they are included in the model}
#'  \item{masterdesign}{The inclusion matrix as constructed by \code{\link{make_master_design}}}}
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'#Korea data
#'data(Korea)
#'ingest_data(Korea)
#'#Kosovo data
#'data(Kosovo)
#'ingest_data(Kosovo)
#'
#'@export
ingest_data = function(xdat)  {
  nlists = dim(xdat)[2]  - 1
  listnames = dimnames(xdat)[[2]][-(1+nlists)]
  if (is.null(listnames)) listnames = LETTERS[1:nlists]
  ncaps = 2^nlists
  nobs = rep(0, ncaps)
  names(nobs)= (1:ncaps)
  nstar = nobs
  ncount = xdat[, 1+nlists]
  xcap = apply(xdat[, -(1+nlists)], 1, encode_capture)
  for (j in (1:length(ncount))) nobs[xcap[j]] = nobs[xcap[j]]+ ncount[j]
  nobs[1]=0
  for (i in (1:ncaps)) nstar[i] = sum(nobs[descendants(i, nlists)])
  notestimable = rep(F, ncaps)
  notestimable[descendants((1:ncaps)[nstar==0], nlists, omitk=T)] = T
  masterdesign = make_master_design(nlists)
  return(list(nobs=nobs, nstar=nstar, nlists=nlists, listnames=listnames, data=xdat,
              notestimable = notestimable, masterdesign= masterdesign))
}
#' Find the vector of captures corresponding to a given hierarchical model
#'
#' Given a hierarchical model, find the vector of all the corresponding encoded captures
#'
#' @param modelstr A given hierarchical model
#' @param findancestors If T then find all the captures.  If F then just return the encoded defining histories of the hierarchy
#'
#' @return The encoded capture histories that corresponds to the row number of the capture history data set
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'modelstr = "[12,23]"
#'convert_from_hierarchy(modelstr)
#'modelstr = "[12,3]"
#'convert_from_hierarchy(modelstr, findancestors=F)
#'
#'@export
convert_from_hierarchy = function(modelstr, findancestors=T) {
  # first decode to numerical vectors of root capture histories to obtain a list of vectors
  #  each of which gives the captures in the capture history of the particular root
  zz = lapply(strsplit(unlist(strsplit(substring(strsplit(modelstr, split="]"), 2), ",")), split=""),as.numeric)
  nlists = max(unlist(zz))
  encode1 = function(z) 1+ sum(2^(z-1))
  captures = unlist(lapply(zz, encode1))
  # now find all the captures in the hierarchy
  if (findancestors) captures = unique(ancestors(captures, nlists))
  return(captures)
}
#' Encode capture history
#'
#' Given a 0/1 capture history, encode it as number that corresponds to the row number of the capture history data set
#'
#' @param z The capture history to be encoded, as a logical vector or a vector of 0s and 1s
#'
#' @return The capture history encoded as a number that corresponds to the row number of the capture history data set
#'
#' @examples
#' encode_capture(c(1,0,0,0,0))
#' encode_capture(c(1,1,1,1,0))
#' encode_capture(c(T,F,T,F))
#'
#' @export
encode_capture = function(z) {
  nlists = length(z)
  k = 1+sum(z*2^{(0:(nlists-1))})
  return(k)
}
#' Find BCa confidence intervals for all possible ntop
#'
#'@param z The output of \code{assemble_bic}, \code{bootstrapcal} and \code{jackknifecal}
#'@param alpha Levels of confidence intervals to be constructed and assessed
#'@param maxorder Maximum order of models to be included
#'@param BCaFD If T, return also the probabilities of the estimates
#'
#'
#'
#'
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
#'z=assemble_bic(Korea, checkexist=T)
#'z=bootstrapcal(z, checkexist=T)
#'z=jackknifecal(z,checkexist=T)
#'ntopBCa(z)
#'
#'@export
ntopBCa = function(z,alpha=c(0.025, 0.05, 0.1, 0.16,0.2, 0.5, 0.8, 0.84, 0.9, 0.95, 0.975), maxorder=Inf, BCaFD=FALSE) {
  # set up and calculate order of model scores
  # restrict to models of prescribed order
  choosemods = (z$res[,3] <= maxorder)
  bootabund = z$bootabund[choosemods,]
  bootbic = z$bootbic[choosemods, ]
  #
  nmods=dim(bootabund)[1]
  modelslist = dimnames(bootabund)[[1]]
  nreps = dim(bootabund)[2]
  modscores = (1:nmods)
  modelorder= order(modscores)
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
#' Acceleration factor calculation
#'
#' Returns acceleration factors for each value of ntop for given modelorder
#'
#' @param z The output of \code{assemble_bic}, \code{bootstrapcal} and \code{jackknifecal}
#' @param modelorder the order of the considered models
#'
#' @return the estimated acceleration factor
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' z=assemble_bic(Korea, checkexist=T)
#' z=jackknifecal(z, checkexist=T)
#' nmods= dim(z$jackabund)[1]
#' modscores = (1:nmods)
#' modelorder= order(modscores)
#' bicktopahatcal(z,modelorder)
#'
#'
#'@export
bicktopahatcal = function(z, modelorder) {
  # set up
  nmods=length(modelorder)
  zja=z$jackabund[modelorder,]
  zjackest=zja
  zjbic = z$jackbic[modelorder,]
  datobs = z$countsobserved
  mdat = length(datobs)
  # take care because a lot of the entries will be NA but only in columns corresponding to zero data count
  #   in the original data.
  md1 = (1:mdat)[datobs>0]
  ndat = sum(datobs)
  # find all the jackknife estimates
  for (j in md1) {
    jrec = which.min(zjbic[,j])
    zjackest[jrec:nmods,j] = zja[jrec,j]
    while (jrec >1) {
      jrec0=jrec-1
      jrec = which.min(zjbic[1:jrec0,j])
      zjackest[jrec:jrec0,j] = zja[jrec,j]
    }
  }
  # now find the vector of acceleration parameters
  #  replace NAs by 0s to avoid computational errors (which don't affect the result)
  zjackest[is.na(zjackest)] = 0
  jackmeans = as.vector(zjackest%*%datobs/ ndat)
  jackd = jackmeans -zjackest
  ahat = (1/6)* jackd^3 %*% datobs / (jackd^2 %*% datobs)^(3/2)
  ahat=as.vector(ahat)
  # now return ahat
  return(ahat)
}
#' Subset matrix
#'
#' The idea of this routine is to reduce either or both of \code{ntopmodels} and \code{maxorder} without the need to recalculate
#' any actual model fits.
#'
#' The routine subsets the results matrix as part of the output given by \code{assemble_bic}
#' based on specified parameters \code{ntopmodels} and \code{maxorder}. It returns the subsetted matrix,
#' original data matrix with capture histories and counts and the new actual value of \code{maxorder} (reducing it
#' from the input value if necessary or if the default input value of \eqn{\infty} is used).
#'
#'
#' @param z output from \code{assemble_bic} or a list output from applying \code{assemble_bic}, \code{jackknifecal} and \code{bootstrapcal}.
#' @param ntopmodels number of top models.  If (taking into account any change in the maximum order of models) there are fewer than
#' \code{ntopmodels} in the data supplied, then it will be reduced to that value.  If it is not specified then there will be no
#' reduction in the number.
#' @param maxorder the maximum order of the models to be considered.  If not specified, it will be set to the corresponding
#'  value in the input data.
#'
#'@return  A list with the following components
#' \describe{
#'   \item{res}{a matrix containing models being considered, abundance, BIC and their ordered after being subsetted by maxorder and ntopmodels}
#'   \item{xdata}{Original data matrix with counts and capture histories}
#'   \item{maxorder}{The maximum order of models considered after subsetting}
#'   \item{jackabund}{Jackknife abundance matrix, subsetted by maxorder and ntopmodels}
#'   \item{jackbic}{Jackknife BIC matrix, subsetted by maxorder and ntopmodels}
#'   \item{countsobserved}{Capture counts in the same order as the columns of \code{jackabund} and \code{jackbic}}
#'   \item{bootabund}{Bootstrap abundance matrix, subsetted by maxorder and ntopmodels}
#'   \item{bootbic}{Bootstrap BIC matrix, subsetted by maxorder and ntopmodels}
#'   } If the input only has the output from \code{assemble_bic}, the last five items of the list
#'   do not appear.
#'
#' @examples
#' data(UKdat_5)
#' #Example 1
#' z=assemble_bic(UKdat_5)
#' subsetmat(z,ntopmodels=100,maxorder=2)
#' subsetmat(z)
#' #Example 2
#' z=assemble_bic(UKdat_5)
#' z=jackknifecal(z)
#' z=bootstrapcal(z)
#' subsetmat(z,ntopmodels=100,maxorder=2)
#'
#' @export
subsetmat <- function(z,
                      ntopmodels = Inf,
                      maxorder = Inf) {
  nrows = dim(z$res)[1]
  #
  # work out which rows correspond to models of the order to be kept
  #
  originalorder = max(z$res[, "modelsorder"])
  # bug corrected in next two lines...
  z$maxorder = min(maxorder, originalorder)
  keeprows = (1:nrows)[z$res[, "modelsorder"] <= z$maxorder]
  #
  # only keep the ntopmodels number of models
  #
  ntopmodels = min(ntopmodels, length(keeprows))
  keeprows = keeprows[1:ntopmodels]
  #
  # carry out the actual subsetting on those matrices which are present
  #
  z$res = z$res[keeprows, ]
  if (!is.null(z$jackbic)) {
    z$jackbic = z$jackbic[keeprows,]
    z$jackabund = z$jackabund[keeprows,]
  }
  if (!is.null(z$bootbic)) {
    z$bootbic = z$bootbic[keeprows,]
    z$bootabund = z$bootabund[keeprows,]
  }
  return(z)
}
#'Order of models
#'
#'This routine returns the order of the model given the model character string
#'
#'@param x a model character string
#'
#'@return the order of the model character string
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'x="[1,2,3,4,5]"
#'modelorder(x)
#'y="[12,24,25,35,45]"
#'modelorder(y)
#'
#'@export

modelorder = function(x) {
  return(max(nchar(unlist(
    strsplit(gsub("\\[|\\]", "", x), ",")
  ))))
}
#' Get a list of all hierarchical models for given number of lists and maximum order
#'
#' Extracts from a larger vector of hierarchical models the ones satisfying the given criterion
#'
#' @param nlists Number of lists
#' @param maxorder Maximum order of models to be returned (defaults to nlists-1)
#' @param modelvec vector of hierarchical models (defaults to hiermodels)
#'
#' @return A list of models satisfying the given criteria
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(hiermodels)
#' # Five lists with maximum order of 4
#' gethiermodels(nlists=5,maxorder=4)
#' # Five lists with maximum order of 2
#' gethiermodels(nlists=5, maxorder=2)
#'
#'@importFrom stats glm.fit na.omit splinefun
#'@export
gethiermodels=function(nlists, maxorder=nlists-1, modelvec=hiermodels){
  # define extraction functions
  model_order = function(x) max(nchar(unlist(strsplit(gsub("\\[|\\]", "", x), ","))))
  nlistfind  = function(x) as.numeric(max(unlist(strsplit(x, split=""))))
  # find relevant models, firstly extracting by number of lists
  zhier = modelvec[sapply(modelvec, nlistfind) == nlists]
  # now consider order
  zhier1 = zhier[sapply(zhier, model_order) <= maxorder]
  return(zhier1)
}
#' The Fienberg-Rinaldo linear program check for the existence of the estimates
#'
#' This routine performs the Fienberg-Rinaldo linear program check in the framework of extended maximum likelihood estimates,
#' the parameters estimates exist if and only if the return value of the check is nonzero
#'
#'
#' @param parset Either the hierarchical representation of the model, or the vector of the corresponding capture histories to the model.
#' @param datlist The output of \code{\link{ingest_data}} on the data set
#' @return The value of the linear program. The parameter estimates within the extended ML framework exist if and only if this value is nonzero.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2022).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
#' models. \emph{Ann. Statist.} \strong{40}, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
#'
#' @examples
#' data(Korea)
#' datlist = ingest_data(Korea)
#' parset ="[0,0,1]"
#' checkident.1(parset,datlist)
#'
#' @export
checkident.1= function(parset, datlist) {
  if (is.character(parset)) parset = convert_from_hierarchy(parset)
  estneginf = parset[datlist$nstar[parset] == 0]
  if (length(estneginf) > 0) {
    parset = setdiff(parset, estneginf)
    datestzero = descendants(estneginf, datlist$nlists)
    amat = t(datlist$masterdesign[-(datestzero-1),parset])
  } else {
    amat = t(datlist$masterdesign[, parset])
  }
  tt = datlist$nstar[parset]
  npar = dim(amat)[1]
  nobs = dim(amat)[2]
  f.obj = c(rep(0, nobs), 1)
  const.rhs = c(tt, rep(0, nobs))
  amat1 = cbind(amat, rep(0, npar))
  amat2 = cbind(diag(nobs), rep(-1, nobs))
  const.mat = rbind(amat1, amat2)
  const.dir = c(rep("=", npar), rep(">=", nobs))
  zlp = lpSolve::lp("max", f.obj, const.mat, const.dir, const.rhs)
  return(zlp$objval)
}
#'Models BICS, abundance and maxorder.
#'
#'This routine sorts the models in increasing order according to their BICs, returns the sorted models with their
#'corresponding BICs and abundance. The original data as well as the maxorder of the models
#'considered are returned as well.
#'
#'@param xdata The original data matrix with capture histories and counts.
#'@param maxorder Maximum order of models to be included
#'@param checkexist If T then the Fienberg-Rinaldo condition is checked for each model
#'@param removeFRfail If checkexist is T then models which fail the FR condition are removed from the results
#' @param ... Parameters to be fed to \code{gethiermodels}.
#'
#'@return A list with the following components
#'\describe{
#' \item{res}{A matrix with the models considered, their abundance, BIC and their order, sorted
#'  into increasing order of BIC}
#'   \item{xdata}{The original data matrix with capture histories and counts.}
#'   \item{maxorder}{Order parameter that was feed into \code{gethiermodels}}
#'   }
#'
#'
#'@examples
#'data(hiermodels)
#'data(Korea)
#'assemble_bic(Korea, checkexist=T)
#'
#'@export
assemble_bic <-
  function(xdata,maxorder=dim(xdata)[2]-2, checkexist=T, removeFRfail=T, ...){
    # number of lists
    nlists=dim(xdata)[2]-1
    hiermodels_cons=gethiermodels(nlists,maxorder=maxorder, modelvec=hiermodels)
    # Use ingest data for fit_hier_model
    ing_dat= ingest_data(xdata)
    #Obtain fit for each considered hierarchical model
    hiermodfit=lapply(hiermodels_cons,
                      fit_hier_model,xdatin=ing_dat,checkid=checkexist, ...)


    #Create matrix with abundance and BIC for corresponding
    #hierarchical model

    res=matrix(NA, nrow=length(hiermodels_cons),ncol=2)
    rownames(res)=hiermodels_cons
    colnames(res)=c("abundance", "BIC")
    for (i in 1:length(hiermodels_cons)){

      res[i,1]=hiermodfit[[i]]$abundance
      res[i,2]=hiermodfit[[i]]$bic

    }
    if (removeFRfail) res = na.omit(res)
    # arrange rows in order of BIC
    res = res[order(res[, 2]),]
    modelnames = dimnames(res)[[1]]
    modelsorder = sapply(modelnames, modelorder)
    res = cbind(res, modelsorder)
    maxorder = max(modelsorder)
    return(list(
      res = res,
      xdata = xdata,
      maxorder = maxorder
    ))
  }
#' Fit a hierarchical model taking account of possible sparsity
#'
#' @param xdatin data obtained using \code{ingest_data}
#' @param hiermod hierarchical model to fit
#' @param bicRcap  if T then use the Rcapture convention that the BIC sample size is the number of cases observed.  Otherwise use the number of cells in the Poisson log linear model.
#' @param checkid if T then \code{checkident.1} is called inside the routine
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' xdatin = ingest_data(Korea)
#' fit_hier_model(xdatin,"[12,23]")
#' fit_hier_model(xdatin,"[12,3]")
#'
#' @importFrom stats glm.fit na.omit splinefun
#' @export
fit_hier_model= function(xdatin, hiermod, bicRcap=T, checkid=F) {

  # convert hiermod to a vector of encoded parameters
  parvec = convert_from_hierarchy(hiermod)
   #  if appropriate, perform check for identification and existence of model fit
  if (checkid==T) {if (checkident.1(parvec, xdatin)==0){
      zglm = NULL
      zglm$abundance =NA
      zglm$bic =NA
     return(zglm)}}
  #
  npars = length(parvec)
  nlists = xdatin$nlists
  nobs = 2^nlists - 1
  # find sparse parameters and pars to be estimated
  sparsepars = parvec[xdatin$nstar[parvec]==0]
  sparseparest = rep(-Inf, length(sparsepars))
  names(sparseparest) = sparsepars
  estpars = setdiff(parvec,sparsepars)
  # find data points to be removed and to be kept. Note that rows of masterdesign matrix are labelled from 2 upwards,
  # hence the offset in removedat, which is then restored when considering the nobs vector itself

  removedat= sort(unique(unlist(lapply(sparsepars, descendants, nlists=nlists)))) -1
  keepdat = setdiff((1:nobs), removedat)
  # construct design matrix and data vector for fitting
  xdes = xdatin$masterdesign[keepdat, estpars]
  yobs = xdatin$nobs[1+keepdat]
  # carry out glm fit and calculate estimated abundance, BIC and AIC (using the number of cells not the number of cases as the number of data points)
  # note that the "structural zeroes" have fitted value zero and so their log likelihood is zero.
  zglm = glm.fit(xdes, yobs, family=poisson())
  zf = zglm$fitted.values

  zglm$abundance = sum(yobs) + exp(zglm$coefficients[1])
  loglikhat = sum(yobs*log(zf) - zf - lfactorial(yobs))
  if (bicRcap) nobs = sum(yobs)
  zglm$bic = npars*log(nobs) - 2*loglikhat
  zglm$aic = 2*npars - 2*loglikhat
  zglm$neginfpars = sparsepars
  zglm$coefficients = c(zglm$coefficients, sparseparest)
  return(zglm)
}
#' Find all neighbouring hierarchical model to a given one
#'
#'@param modelstr Model string written in hierarchical form
#'@param nlists Number of lists.
#'@param keepmaineffects If T, keep the main effects. If F remove.
#'@param maxorder Maximum order of models to be included
#'
#'@return neighbour hierarchical models
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' modelstr = "[12,23]"
#' find_neighbour_hierarchies(modelstr)
#'@export
find_neighbour_hierarchies = function(modelstr, nlists=NA, keepmaineffects=T, maxorder=nlists-1){
  # First find models obtained by adding a capture history
  if (is.na(nlists)) nlists=as.numeric(max(unlist(strsplit(modelstr, split=""))))
  zhierroots = convert_from_hierarchy(modelstr, F)
  zhier = unique(ancestors(zhierroots, nlists))
  znew = boundary_captures(zhier,nlists)
  newmodels = lapply(znew, function(x) union(zhier, x))
  outerneighbours = lapply(newmodels, convert_to_hierarchy, nlists=nlists)
  #
  order_ok = (lapply(outerneighbours, modelorder) <= maxorder)
  outerneighbours = outerneighbours[order_ok]
  # Now find those by removing a history, necessarily one of the defining histories of the hierarchy
  #   but not any which are a single list
  if (keepmaineffects) zhierroots = setdiff(zhierroots, 1+2^(0:(nlists-1)))
  newmodels1 = lapply(zhierroots, function(x) setdiff(zhier, x))
  innerneighbours = lapply(newmodels1, convert_to_hierarchy, nlists=nlists)
  neighbours = c(innerneighbours,outerneighbours)
  neighbours= unlist(neighbours)
  # nbsorder = sapply(neighbours,modelorder)
  # neighbours = neighbours[nbsorder == nlists-1]
  return(neighbours)
}
#' Find hierarchical representation of a vector of captures
#'
#' Given a vector of encoded captures defining a hierarchical model, re-express
#' it in hierarchical model form
#'
#' @param kcap A vector of captures
#' @param nlists The number of lists
#'
#' @return A hierarchical representation of the vector of encoded captures.
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' kcap=c(1,2,3,5,4)
#' nlists=3
#' convert_to_hierarchy(kcap, nlists)
#'@export
convert_to_hierarchy = function(kcap, nlists) {
  # first find the roots of the hierarchy
  # this will be those in kcap none of whose children are in the set
  nk = length(kcap)
  kroot = vector(length=nk)
  for (jk in (1:nk)) {
    kch = child_captures(kcap[jk], nlists)
    kroot[jk] = (length(intersect(kch, kcap))==0)
  }
  rootcaps = kcap[kroot]
  # now find the captures in each rootcaps
  # also find a weight that will get them into the right order
  nr = length(rootcaps)
  rootdecode = vector("character", length=nr)
  rootweight = vector("numeric", length=nr)
  for (jr in (1:nr)){
    zr = (1:nlists)[decode_capture(rootcaps[jr],nlists)]
    rootweight[jr] = length(zr) + sum(0.5^zr)
    rootdecode[jr] = paste0(zr, collapse="")
  }
  rootdecode = rootdecode[order(rootweight, decreasing=T)]
  # now concatenate to find full hierarchical representation
  zhier = paste0("[", paste(rootdecode, collapse=","), "]", collapse="")
  return(zhier)
}
#' Given a vector of captures, find those which are not in the vector but all of whose parents are
#'
#' Call the resulting set the "boundary".  Supposing that the current set of captures is a hierarchical model, that property
#' will be preserved if a capture in the boundary is added to it. The routine is called internally by \code{find_neighbour_hierarchies}.
#'
#'
#' @param kcap An encoded capture history that corresponds to the row number of the capture history data set
#' @param nlists The total number of lists
#'
#' @return a vector giving the encoded versions of the descendants
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#'modelstr = "[12,23]"
#'nlists=4
#'zhierroots = convert_from_hierarchy(modelstr, F)
#'zhier = unique(ancestors(zhierroots, nlists))
#'boundary_captures(zhier,nlists)
#'
#' @export
boundary_captures = function(kcap, nlists) {
  #  Here kcap is a vector of captures.  Find the captures which are not in kcap but
  #    all of whose parents are.
  #
  # first find all children of captures in kcap
  kchild = NULL
  for (k in kcap) kchild = union(kchild, child_captures(k, nlists))
  # now exclude kcap and then check if all parents are in kcap
  kchild = setdiff(kchild, kcap)
  nkc = length(kchild)
  kinclude = vector(length=nkc)
  for (kc in (1:nkc)) {
    # I think this was the problem....it was 3 and now is nlists
    kpar = parent_captures(kchild[kc], nlists)
    kinclude[kc] = setequal(kpar, intersect(kpar, kcap))}
  kboundary = kchild[kinclude]
  return(kboundary)
}
#' Find unique patterns in matrix columns
#'
#' Given a matrix (for example of bootstrap replications) construct the matrix of unique patterns of non-zeroes,
#' together with a vector of pointers back to that matrix.
#'
#' @param x a matrix
#' @returns The original data \code{x} with the additional components {\describe{
#'   \item{yuniq}{matrix of unique patterns of non-zeroes/zeroes in the columns of \code{x}}
#'   \item{pointers}{vector of length dim(x)[2] giving the column of yuniq corresponding to each column of x}
#' }}
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'data(Korea)
#'z=assemble_bic(Korea)
#'z=bootstrapcal(z)
#'find_unique_patterns(z$bootreplications)
#'
#' @export
find_unique_patterns = function(x) {
  # find non-zeroes and convert resulting matrix to numeric
  xn = (x > 0)
  xn[xn] = 1
  # find unique columns
  yuniq = unique(xn, MARGIN = 2)
  # find pointers
  nuniq = dim(yuniq)[2]
  pointers = numeric(nuniq)
  for (j in (1:(dim(x)[2]))) {
    pointers[j] = which.max(apply(yuniq, 2, identical, y = xn[, j]))
  }
  # construct output and return
  zreturn = list(x=x, yuniq=yuniq, pointers=pointers)
  return(zreturn)
}
#' Carry out the Fienberg-Rinaldo procedure on an array of data vectors and a vector of models
#'
#' Suppose we have a collection of different data outcomes on the same set of capture histories and a vector of models.
#' Typically the data outcomes will be bootstrap replications. This routine finds the unique support patterns among the data
#' and hence economises the task of finding which model/data combinations satisfy the Fienberg-Rinaldo condition
#'
#' @param x a matrix of data observations for a common capture matrix
#' @param xcap the incidence matrix of the capture histories corresponding to the rows of x
#' @param zmods a vector of models
#'
#' @returns a matrix with rows corresponding to the models and columns to the columns of x, with elements
#' taking the value T if the FR linear program for a vector of 0s and 1s with the same zero pattern as the x data yields a strictly
#' positive value
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
#' z=assemble_bic(Korea)
#' z=bootstrapcal(z)
#' n2=dim(z$xdata)[2]
#' checkident.2(z$bootreplications, z$xdata[,-n2], dimnames(z$res)[[1]])
#'
#'@export
checkident.2 = function(x, xcap, zmods) {
  # set up the unique patterns within x and initialise
  zu = find_unique_patterns(x)
  yuniq = zu$yuniq
  nyuniq = dim(yuniq)[2]
  # first convert all models from hierarchy
  nmods = length(zmods)
  zmodsc = lapply(zmods, convert_from_hierarchy)
  # carry out F-R procedure for each column of yuniq for each model
  # set up results matrix
  frmaty = matrix(
    NA,
    nrow = nmods,
    ncol = nyuniq,
    dimnames = list(zmods, 1:nyuniq)
  )
  # consider unique data patterns in turn and solve the FR criterion over the various models

  for (j in (1:nyuniq)) {
    # construct the data matrix and ingest the data
    datlist = ingest_data(cbind(xcap, yuniq[,j]))
    # now consider the various models for this data pattern
    frmaty[,j] = unlist(lapply(zmodsc, checkident.1, datlist=datlist))
  }
  # construct the results matrix for the original data
  frmaty = (frmaty>0)
  frmat = frmaty[, zu$pointers]
  return(frmat)
}
#' Bootstrap abundance and bic
#'
#' This routine takes the output from \code{assemble_bic} or \code{subsetmat} and returns bootstrap
#' abundance matrix and BIC matrix. This version makes use of the \code{checkident.2}.
#'
#' @param z Results from \code{assemble_bic} or \code{subsetmat}
#' @param nboot The number of bootstrap replications.
#' @param iseed Integer seed to allow for replicability.
#' @param checkexist If \code{checkexist=TRUE}, check for existence, else it does
#' not check for existence.
#' @param saveinterval If this is set to a finite value, the output list \code{z}
#' will be saved every time the number of replications is a multiple of it. A message
#' will be printed every time it is saved.
#' @param savefile The file to which the output will be saved if \code{saveinterval}
#' is set to a finite value.
#'
#'
#' @return The original input list \code{z} with the additional components
#' \describe{
#'   \item{bootabund}{Bootstrap abundance matrix}
#'   \item{bootbic}{Bootstrap BIC matrix}
#' }
#' If there are already components with these names they will be overwritten.
#'
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#'@examples
#'z=assemble_bic(Korea)
#'bootstrapcal(z, checkexist=T, saveinterval=50, savefile="Koreabootresults.Rdata")
#'
#' @export
bootstrapcal <- function(z,
                           nboot = 1000,
                           iseed = 1234,
                           checkexist=T,
                           saveinterval = Inf,
                           savefile = "bootout.Rdata") {
  set.seed(iseed)
  zdat = tidylists(z$xdata, includezerocounts = T)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  z$countsobserved = countsobserved
  nobs = sum(countsobserved)
  topmodels = dimnames(z$res)[[1]]
  ntopmodels = length(topmodels)
  z$bootabund = matrix(
    NA,
    nrow = ntopmodels,
    ncol = nboot,
    dimnames = list(topmodels, 1:nboot)
  )
  z$bootbic = z$bootabund
  # construct all the bootstrap data
  z$bootreplications = rmultinom(nboot, nobs, countsobserved)
  # carry out the F-R check
  if (checkexist) frmat = checkident.2(z$bootreplications, zdat[, -n2], topmodels)
  #
  for (j in (1:nboot)) {
    zdat[, n2] = z$bootreplications[, j]
    ing_dat = ingest_data(zdat)
    for (imod in (1:ntopmodels)) {
      if (!checkexist || frmat[imod, j]) {
        zfit = fit_hier_model(xdatin = ing_dat,
                              hiermod = topmodels[imod],
                              checkid = F)
        z$bootabund[imod, j] = zfit$abundance
        z$bootbic[imod, j] = zfit$bic
      }
    }
    # note that z$bootabund and z$bootbic are initialised to NA, so there is no need to do anything if frmat[imod, j] is FALSE
    if (j %% saveinterval == 0) {
      save(z, file = savefile)
      cat("File saved for j = ", j, "\n")
    }
  }
  return (z)
}
#' Jackknife abundance and Jackknife bic
#'
#' This routine takes the output from \code{subsetmat} or from \code{assemble_bic} and returns the jackknife
#' abundance matrix and jackknife BIC matrix.
#'
#' @param z Results from \code{assemble_bic} or \code{subsetmat}.
#' @param checkexist If \code{checkexist=TRUE}, check for existence in cases where the jackknife introduces an additional zero, else it does
#' not check for existence.  Note that in the current version it is assume that models for which the fit doesn't exist for the original data
#' have already been excluded.
#'
#' @return A list with the following components
#' \describe{
#'   \item{jackabund}{Jackknife abundance matrix}
#'   \item{jackbic}{Jackknife BIC matrix}
#'   \item{countsobserved}{Capture counts in the same order as the columns of \code{jackabund} and \code{jackbic}}
#' }
#'
#'@references
#' Silverman, B. W., Chan, L. and  Vincent, K., (2024).
#' Bootstrapping Multiple Systems Estimates to Account for Model Selection
#' \emph{Statistics and Computing}, \strong{34(44)},
#' Available from \url{https://doi.org/10.1007/s11222-023-10346-9}.
#'
#' @examples
#' data(Korea)
#' z=assemble_bic(Korea)
#' jackknifecal(z,checkexist=T)
#'
#' @export
jackknifecal <- function(z, checkexist = T) {

  zdat = tidylists(z$xdata, includezerocounts = T)
  n1 = dim(zdat)[1]
  n2 = dim(zdat)[2]
  countsobserved = zdat[, n2]
  nobs = sum(countsobserved)
  jest = (countsobserved > 0)
  topmodels = dimnames(z$res)[[1]]
  ntopmodels = length(topmodels)
  jackabund = matrix(
    NA,
    nrow = ntopmodels,
    ncol = n1,
    dimnames = list(topmodels, 1:n1)
  )
  jackbic = jackabund

  #hiermod = topmodels[imod]
  # now the relevant jackknife values for this particular model
  for (j in (1:n1)[jest]) {
    yy = countsobserved
    yy[j] = yy[j] - 1
    zdat[, n2] = yy
    ing_dat = ingest_data(zdat)
    # make sure that if a new zero has been introduced the check is carried out
    checkidj = (checkexist & (yy[j]==0))
    for (imod in (1:ntopmodels)) {
      zfit = fit_hier_model(xdatin = ing_dat, hiermod = topmodels[imod],checkid=checkidj)
      jackabund[imod, j] = zfit$abundance
      jackbic[imod, j] = zfit$bic
    }
  }
  z$jackabund = jackabund
  z$jackbic = jackbic
  z$countsobserved = countsobserved
  return(z)
}

