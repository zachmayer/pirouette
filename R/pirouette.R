#' @title pirouetteControl
#' @description Define a random sparse projection
#'
#' @param mtry dimensionality of the random projection
#' @param nobs number of observations to use to fit each tree
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
#' @importMethodsFrom Matrix rsparsematrix %*%
#' @return a sparse Matrix
splat <- function(x, mtry=2){
  D <- ncol(x)
  nnz <- ceiling( D * 1/sqrt(D))
  r <- rsparsematrix(
    nrow = ncol(x),
    ncol = mtry,
    nnz = nnz,
    rand.x = function(n) sample(c(-1L,1L), nnz, replace=TRUE)
    )
  x %*% r
}

#' @title pirouette
#' @description Random rotation forests for sparse data
#'
#' @param x a sparse matrix of x variables
#' @param y the target variable for classification or regression
#'
#'@references \url{http://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf}
#'
#' @import methods
#' @importMethodsFrom Matrix rsparsematrix %*%
#' @importFrom foreach %do% %dopar%
#'
#' @return an object of class pirouette
#' @export
#'
pirouette <- function(x, y, ctrl = pirouetteControl()){

  stopifnot(is.vector(y))
  stopifnot(nrow(x) == length(y))


}