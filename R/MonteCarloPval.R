`MonteCarloPval` <- function(nSimu, tL, sqAmatrix_list, Cobs, Seedfixed=1, RandomSimu = 0, singleSim = 10000){
  nrow_tL = nrow(tL) #nrow_tL = ncol_tL
  
  len_sqAmatrix = length(sqAmatrix_list)
  ncol_All_sqAmatrix = nrow_tL*nrow_tL
  Nonzero_sqA = rep(0,len_sqAmatrix)
  
  for(subset in 1:len_sqAmatrix){
    Nonzero_sqA[subset] = sum(sqAmatrix_list[[subset]]!=0)
  }
  MaxCount = max(Nonzero_sqA)
  
  sqAIndex_i = sqAIndex_j = matrix(1,len_sqAmatrix,MaxCount)
  sqAval = matrix(0,len_sqAmatrix,MaxCount)
  for(subset in 1:len_sqAmatrix){
    wh = which(sqAmatrix_list[[subset]]!=0)
    Index = which(sqAmatrix_list[[subset]]!=0,arr.ind=T)
    sqAIndex_i[subset,1:nrow(Index)] = Index[,1]
    sqAIndex_j[subset,1:nrow(Index)] = Index[,2]
    sqAval[subset,1:nrow(Index)] = sqAmatrix_list[[subset]][wh]
    rm(wh,Index)
  }
  
 
  
  nrep = ceiling(nSimu/singleSim)
  if(RandomSimu == 1){
    nsimU = sample(ceiling(singleSim*0.8):(singleSim*1.2),nrep,replace=T)
  }else{
    nsimU = rep(singleSim,nrep)
  }
  dims = rep(0,5)
  dims[1] = nSimu
  dims[2] = nrow_tL
  dims[3] = len_sqAmatrix
  dims[4] = MaxCount
  dims[5] = nrep
  
  
  if(Seedfixed){
    set.seed(12345)
  }
  
  pval = 1
  sqAIndex_i = sqAIndex_i - 1
  sqAIndex_j = sqAIndex_j - 1
  Z = .C("MonteCarloPval", as.double(t(tL)), as.integer(t(sqAIndex_i)), as.integer(t(sqAIndex_j)), as.double(t(sqAval)),as.integer(dims),
   as.double(Cobs), as.integer(nsimU), pval = as.double(pval), as.integer(Nonzero_sqA), PACKAGE="SubsetRV")
  pval = Z$pval
  ll = list()
  ll[[1]] = pval
  return(ll)
}

