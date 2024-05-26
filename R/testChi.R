`testChi` <- function(nchi = 100000){


  dims = rep(0,1)
  dims[1] = nchi
  
 
  Chioutput = rep(0,nchi)

  
  Z = .C("testChiC", Chioutput=as.double(Chioutput),as.integer(dims), PACKAGE="SubsetRV")

  Chioutput= Z$Chioutput
  return(Chioutput)
}

