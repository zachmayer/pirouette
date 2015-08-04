#' @title pirouetteControl
#' @description Define a random sparse projection
#'
#' @param mtry dimensionality of the random projection
#' @param nobs maximum number of observations used per tree
#' @return a list
#' @export
pirouetteControl <- function(
  mtry = 2,
  nobs = 500
  ){
  list(
    mtry = mtry,
    nobs = nobs
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
    x_splat$x <- x_splat$x[,]
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
#' m <- enpointe(matrix(runif(1000), nrow=100), runif(100))
#' predict(m, matrix(runif(10), ncol=10))
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
#'@references \url{http://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf}
#'
#' @importFrom foreach %do% %dopar%
#' @return an object of class pirouette
#' @export
pirouette <- function(x, y, ctrl = pirouetteControl(), ...){

  stopifnot(is.vector(y))
  stopifnot(nrow(x) == length(y))


}