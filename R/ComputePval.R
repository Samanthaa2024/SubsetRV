`ComputePval` <- function(nsimchi = 100000, allsets, Cobs, Tsetsvec,Leglambda,tmplambda, IndepIndex){

  #length(tmplambda): number of sets
  LambdaMatrix = matrix(NA,Leglambda,length(tmplambda))
  for(tt in 1:length(tmplambda)){
    LambdaMatrix[,tt] = tmplambda[[tt]]
  }
  nrow_LambdaMatrix = nrow(LambdaMatrix)
  ncol_LambdaMatrix = ncol(LambdaMatrix)# number of columns is the number of sets
  
  
  nrow_allsets = nrow(allsets)


  dims = rep(0,7)
  dims[1] = nsimchi
  dims[2] = nrow_LambdaMatrix
  dims[3] = ncol_LambdaMatrix
  dims[4] = Leglambda
  dims[5] = length(tmplambda)
  dims[6] = nrow_allsets
  dims[7] = IndepIndex
 
  vals = rep(0,2)
  vals[1] = Cobs
  
  Z = .C("ComputePval", as.double(t(LambdaMatrix)),as.integer(dims),
   vals = as.double(vals), PACKAGE="SubsetRV")

  vals= Z$vals
  pval = vals[2]
  return(pval)
}

