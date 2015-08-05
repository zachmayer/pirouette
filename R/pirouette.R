#' @title pirouetteControl
#' @description Define a random sparse projection
#'
#' @param ntrees number of models to fit
#' @param newdim dimensions to project into
#' @param allowParallel the function to run in parallel (if a foreach backend is loaded)
#' @return a list
#' @export
pirouetteControl <- function(
  ntrees = 100,
  newdim = 2,
  allowParallel = TRUE
  ){
  list(
    ntrees = ntrees,
    newdim = newdim,
    allowParallel = allowParallel
  )
}

#' @title splat
#' @description Splat a matrix into a lower dimension.
#'
#' @param x a sparse matrix
#' @param newdim new dimensionality
#' @param prob the probability each entry of the rotation matrix is non-zero.  1/3 is one heuristic to use (for 3x faster processing), 1/sqrt(ncol(x)) is another, for sqrt(x)-fold faster processing.
#' @param retx If true, return the splatted matrix as well as the rotation.
#'
#' @import methods
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix rsparsematrix
#' @return a splatted Matrix
#'
#' @references \url{http://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf}
#'
#' @export
#' @examples
#' splat(matrix(runif(100), ncol=10), 2)
splat <- function(x, newdim=2, prob=1/sqrt(ncol(x)), retx=TRUE){
  D <- ncol(x)
  nnz <- ceiling(D * newdim * prob)
  r <- rsparsematrix(
    nrow = ncol(x),
    ncol = newdim,
    nnz = nnz,
    rand.x = function(n) sample(c(-1L,1L), nnz, replace=TRUE)
    )
  out <- list(
    r = r,
    x = NULL
  )
  if(retx){
    out[['x']] <- x %*% r
  }
  class(out) <- 'splat'
  return(out)
}

#' @title Make predictions from a splatted object
#' @description Splat a new matrix into the same rotation as an old one.
#'
#' @param object a splat object
#' @param newx A sparse matrix.
#' @param ... ignored
#'
#' @import methods
#' @importMethodsFrom Matrix %*%
#' @method predict splat
#' @return a sparse Matrix
#'
#' @export
#' @examples
#' x <- splat(matrix(runif(100), ncol=10), 2)
#' predict(x, matrix(runif(10), ncol=10))
predict.splat <- function(object, newx, ...){
  newx %*% object[['r']]
}

#' @title Random rotation tree
#' @description This is the workhorse function of the pirouette package.  It
#' splats the dateset and fits an rpart model to it.
#'
#' @param x a sparse matrix of x variables
#' @param y the target variable for classification or regression
#' @param weights option case weights
#' @param newdim dimensions to project into
#' @param gbm_control Control arguments to pass to the gbm call. Use 2 trees.
#' @param prob passed to splat
#' @param ... passed to rpart
#'
#' @importFrom Matrix drop0 rowSums
#' @importFrom gbm gbm.fit
#' @export
#'
#' @return an object of class enpointe
enpointe <- function(
  x, y, weights=NULL,
  newdim=2, prob=1/sqrt(ncol(x)),
  gbm_control = list(
    n.trees=1,
    interaction.depth=2,
    shrinkage=0.1,
    verbose=FALSE
  ), ...){

  #Splat matrix
  x_splat <- splat(x, newdim=newdim, prob=prob)
  x_splat$x <- drop0(x_splat$x)

  #Subset to non-zero rows
  #TODO: If all rows are 0, skip fitting, and always predict median
  #TODO: Warning to splat with more dimensions
  keep <- rowSums(sign(x_splat$x)) != 0
  x_splat$x <- x_splat$x[keep,]
  y <- y[keep]

  #Fit model
  model <- gbm.fit(
    x=data.frame(as.matrix(x_splat$x)),
    y=y,
    w=weights,
    n.trees=gbm_control$n.trees,
    interaction.depth=gbm_control$interaction.depth,
    shrinkage=gbm_control$shrinkage,
    verbose=gbm_control$verbose,
    ...)

  #Todo: trim model more
  model$call <- NULL
  out <- list(
    r = x_splat$r,
    model = model,
    gbm_control = gbm_control
  )
  class(out) <- 'enpointe'
  return(out)
}

#' @title Make predictions from a enpointe object
#' @description Predict from a single rotation tree.
#'
#' @param object an enpointe object
#' @param newx A sparse matrix.
#' @param ... passed to predict.rpart
#'
#' @import methods
#' @importMethodsFrom Matrix %*%
#'
#' @method predict enpointe
#' @return a vector
#' @export
#' @examples
#' a <- matrix(runif(10000), ncol=10)
#' b <- as.vector(a %*% runif(ncol(a)))
#' m <- enpointe(a, b, distribution='gaussian')
#' predict(m, matrix(runif(1000), ncol=10))
predict.enpointe <- function(object, newx, ...){
  newx <- as.matrix(newx %*% object[['r']])
  predict(object[['model']], newx, type = 'link', n.trees=object$gbm_control$n.trees, ...)
}

#' @title pirouette
#' @description Random rotation forests for sparse data
#'
#' @param x a sparse matrix of x variables
#' @param y the target variable for classification or regression
#' @param prob the probability each entry of the rotation matrix is non-zero.  1/3 is one heuristic to use (for 3x faster processing), 1/sqrt(ncol(x)) is another, for sqrt(x)-fold faster processing.  Passed through enpointe to splat.
#' @param ctrl a list of control parameters for the algorithm
#' @param gbm_control a list of control parameters for the gbm
#' @param ... passed through enpointe to gbm
#'
#'@references \itemize{
#'  \item{\url{https://stat.ethz.ch/pipermail/r-sig-hpc/2013-January/001575.html}}
#'  \item{\url{http://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf}}
#'  \item{\url{http://stats.lse.ac.uk/fryzlewicz/rre/rre.pdf}}
#'}
#' @importFrom foreach foreach %do% %dopar%
#' @return an object of class pirouette
#' @export
pirouette <- function(
  x, y,
  prob=1/sqrt(ncol(x)),
  ctrl = pirouetteControl(),
  gbm_control = list(
    n.trees=1,
    interaction.depth=2,
    shrinkage=0.1,
    verbose=FALSE
  ), ...){

  stopifnot(is.vector(y) | is.factor(y))
  stopifnot(nrow(x) == length(y))
  if(is.factor(y)){
    if(length(unique(y)) < 2){
      warning('Less than 2 unique factor levels in classification problem.  Fitting models will probably work, but prediction wont.')
    }
    if(length(unique(y)) > 2){
      stop('More than 2 unique factor levels in classification problem.  Multiclass is not yet supported.')
    }
  }

  `%op%` <- if(ctrl$allowParallel) `%dopar%` else `%do%`

  models <- foreach(i=1:ctrl$ntrees) %op% {
    enpointe(x, y, prob=prob, newdim=ctrl$newdim, ...)
  }

  out <- list(
    models = models,
    ctrl = ctrl,
    prob = prob
  )
  class(out) <- 'pirouette'
  return(out)
}

#' @title Make predictions from a pirouette object
#' @description Predict from a pirouette forest
#'
#' @param object an pirouette object
#' @param newx A sparse matrix.
#' @param allowParallel Allow predictions in parallel
#' @param ... passed to predict.enpointe
#'
#' @method predict pirouette
#' @importFrom foreach foreach %do% %dopar%
#' @return a vector
#' @export
#'
#' @examples
#' nrow <- 10000
#' ncol <- 100
#' m <- pirouette(matrix(runif(nrow * ncol), ncol=ncol), runif(nrow), distribution='gaussian')
#' predict(m, matrix(runif(ncol), ncol=ncol))
predict.pirouette <- function(object, newx, allowParallel = FALSE, ...){

  `%op%` <- if(allowParallel) `%dopar%` else `%do%`

  #Todo: change combine from cbind if multiclass
  p <- foreach(
    i=seq_along(object$models), .combine=cbind, .multicombine=TRUE
  ) %op% {
    predict.enpointe(object$models[[i]], newx, ...)
  }

  #dirty dirty hack
  m <- median(p, na.rm=TRUE)
  p[!is.finite(p)] <- m

  #Return
  #Todo: for multiclass, average matricies, then normalize rows
  return(rowMeans(p))
}
