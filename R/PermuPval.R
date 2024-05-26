`PermuPval` <- function(nperm, Genotype1, Genotype2, Res, Nsets, PermuteIndex = matrix(0,1,1)){

  N = nrow(Res)
  K = ncol(Res)

  QMatrix = matrix(0,nperm,Nsets)
  dims = rep(0,8)
  dims[1] = nperm
  dims[2] = N
  dims[3] = K
  dims[4] = N
  dims[5] = ncol(Genotype1)
  dims[6] = ncol(QMatrix)
  dims[7] = nrow(PermuteIndex)
  dims[8] = ncol(PermuteIndex)


  
  Z = .C("PermuPval", as.double(t(Genotype1)), as.double(t(Genotype2)), as.double(t(Res)), as.integer(dims),
   QMatrix = as.double(t(QMatrix)), as.integer(t(PermuteIndex)), PACKAGE="SubsetRV")

  QMatrix = matrix(Z$QMatrix, nrow = nperm, ncol = Nsets, byrow = TRUE)
  return(QMatrix)
}

