#' @title pirouetteControl
#' @description Define a random sparse projection
#'
#' @param ntrees number of models to fit
#' @param mtry dimensionality of the random projection
#' @param nobs maximum number of observations used per tree
#' @param allowParallel the function to run in parallel (if a foreach backend is loaded)
#' @return a list
#' @export
pirouetteControl <- function(
  ntrees = 100,
  mtry = 2,
  nobs = 500,
  allowParallel = TRUE
  ){
  list(
    ntrees = ntrees,
    mtry = mtry,
    nobs = nobs,
    allowParallel = allowParallel
  )
}

#' @title splat
#' @description Splat a matrix into a lower dimension.
#'
#' @param x a sparse matrix
#' @param mtry new dimensionality
#' @param retx If true, return the splatted matrix as well as the rotation.
#'
#' @import methods
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix rsparsematrix
#' @return a splatted Matrix
#'
#' @export
#' @examples splat(matrix(runif(100), ncol=10), 2)
splat <- function(x, mtry=2, retx=TRUE){
  D <- ncol(x)
  nnz <- ceiling( D * 1/sqrt(D))
  r <- rsparsematrix(
    nrow = ncol(x),
    ncol = mtry,
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
#' @param maxrows the maximum number of rows to use for the tree fit
#' @param ... passed to rpart
#'
#' @importFrom Matrix drop0 rowSums
#' @importFrom rpart rpart
#' @export
#'
#' @return an object of class enpointe
enpointe <- function(x, y, weights=NULL, maxrows=500, ...){

  #Drop 0s
  x_splat <- splat(x)
  x_splat$x <- drop0(x_splat$x)
  keep <- rowSums(sign(x_splat$x)) >= 0
  x_splat$x <- x_splat$x[keep,]
  y <- y[keep]

  #Subset if needed
  if(nrow(x_splat$x) > maxrows){
    keep <- sample(1:nrow(x_splat$x), maxrows)
    x_splat$x <- x_splat$x[keep,]
    y <- y[keep]
  }

  #Fit model
  tmp <- data.frame('.outcome'=y, as.matrix(x_splat$x))
  if(!is.null(weights)){
    model <- rpart(.outcome ~ ., tmp, y=FALSE, ...)
  } else {
    model <- rpart(.outcome ~ ., tmp, weights=weights, y=FALSE, ...)
  }

  #Todo: trim model more
  model$call <- NULL
  out <- list(
    r = x_splat$r,
    model = model
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
#' m <- enpointe(matrix(runif(10000), ncol=10), runif(100))
#' predict(m, matrix(runif(1000), ncol=10))
predict.enpointe <- function(object, newx, ...){
  newx <- data.frame(as.matrix(newx%*% object[['r']]))
  if(object$model$method == 'anova'){
    return(predict(object[['model']], newx, type = 'vector', ...))
  } else if(object$model$method == 'class'){
    out <- predict(object[['model']], newx, type = 'prob', ...)
    if(ncol(out) == 2){
      return(out[,1])
    } else{
      #Dangerous... not planning to support multiclass yet
      return(out)
    }
  } else{
    stop(paste('rpart method', object$model$method, 'not supported'))
  }
}

#' @title pirouette
#' @description Random rotation forests for sparse data
#'
#' @param x a sparse matrix of x variables
#' @param y the target variable for classification or regression
#' @param ctrl a list of control parameters for the algorithm
#' @param ... passed through enpointe to rpart
#'
#'@references \url{https://stat.ethz.ch/pipermail/r-sig-hpc/2013-January/001575.html}
#'
#' @importFrom foreach foreach %do% %dopar%
#' @return an object of class pirouette
#' @export
pirouette <- function(x, y, ctrl = pirouetteControl(), ...){

  stopifnot(is.vector(y) | is.factor(y))
  stopifnot(nrow(x) == length(y))
  if(is.factor(y)){
    if(length(unique(y)) < 2){
      warning('Less than 2 unique factor levels in classification problem.  Fitting models will probably work, but prediction wont.')
    }
    if(length(unique(y)) > 2){
      warning('More than 2 unique factor levels in classification problem.  Fitting models will probably work, but prediction wont.')
    }
  }

  `%op%` <- if(ctrl$allowParallel) `%dopar%` else `%do%`

  models <- foreach(i=1:ctrl$ntrees) %op% {
    enpointe(x, y)
  }

  out <- list(
    models = models
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
#' m <- pirouette(matrix(runif(nrow * ncol), ncol=ncol), runif(nrow))
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
