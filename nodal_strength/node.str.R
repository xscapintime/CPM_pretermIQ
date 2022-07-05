node.str=function(M) {
  # check that input is square
  if (dim(M)[1]!=dim(M)[1]) stop('the input matrix is not square')
  n=dim(M)[1]
  M[seq(from=1,to=n^2,by=(n+1))]=NaN
  rowMeans(M,na.rm=T)
}


