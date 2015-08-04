#http://stackoverflow.com/questions/31570437/really-fast-word-ngram-vectorization-in-r
set.seed(1)
samplefun <- function(n, x, collapse){
  paste(sample(x, n, replace=TRUE), collapse=collapse)
}
words <- sapply(rpois(10000, 3) + 1, samplefun, letters, '')
sents1 <- sapply(rpois(100000, 5) + 1, samplefun, words, ' ')

find_ngrams <- function(dat, n, verbose=FALSE){
  library(pbapply)
  stopifnot(is.list(dat))
  stopifnot(is.numeric(n))
  stopifnot(n>0)
  if(n == 1) return(dat)
  pblapply(dat, function(y) {
    if(length(y)<=1) return(y)
    c(y, unlist(lapply(2:n, function(n_i) {
      if(n_i > length(y)) return(NULL)
      do.call(paste, unname(as.data.frame(embed(rev(y), n_i), stringsAsFactors=FALSE)), quote=FALSE)
    })))
  })
}

text_to_ngrams <- function(sents, n=2){
  library(stringi)
  library(Matrix)
  tokens <- stri_split_fixed(sents, ' ')
  tokens <- find_ngrams(tokens, n=n, verbose=TRUE)
  token_vector <- unlist(tokens)
  bagofwords <- unique(token_vector)
  n.ids <- sapply(tokens, length)
  i <- rep(seq_along(n.ids), n.ids)
  j <- match(token_vector, bagofwords)
  M <- sparseMatrix(i=i, j=j, x=1L)
  colnames(M) <- bagofwords
  return(M)
}

sparse_text_matrix <- text_to_ngrams(sents1)
sparse_text_matrix <- sparse_text_matrix[rowSums(sparse_text_matrix)>1,]
sparse_text_matrix <- sparse_text_matrix[,colSums(sparse_text_matrix)>1]
devtools::use_data(sparse_text_matrix, overwrite = TRUE, compress='xz')
