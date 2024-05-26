`testNorm` <- function(nnorm = 100000){


  dims = rep(0,1)
  dims[1] = nnorm
  
 
  Normoutput = rep(0,nnorm)

  
  Z = .C("testNormC", Normoutput=as.double(Normoutput),as.integer(dims), PACKAGE="SubsetRV")

  Normoutput= Z$Normoutput
  return(Normoutput)
}

