
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "ComputePval.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void ComputePval(double* RLambdaMatrix,int* dims, double* vals)
{
    int nrow_LambdaMatrix = dims[1], ncol_LambdaMatrix = dims[2];
    double **LambdaMatrix;
    
    /* reorganize vector into matrix */
    reorg(RLambdaMatrix, &LambdaMatrix, nrow_LambdaMatrix, ncol_LambdaMatrix);
        
    ComputePvalC(LambdaMatrix, dims, vals);
} 


void ComputePvalC(double** LambdaMatrix, int* dims, double*vals)
{
  
    int nsimchi = dims[0], nrow_LambdaMatrix = dims[1], ncol_LambdaMatrix = dims[2], Leglambda = dims[3], length_tmplambda = dims[4], nrow_allsets = dims[5];
    int i, j, sim, tt, IndepIndex = dims[6];
    double Cobs = vals[0], Count = 0.0, tmp, *tmpchi;
    
    tmpchi = (double *)calloc(nrow_LambdaMatrix, sizeof(double));


    GetRNGstate();
    Count = 0.0;
    for(sim=0; sim<nsimchi; sim++){
      for(tt=0; tt<length_tmplambda; tt++){
        tmp = 0.0;
        for(j=0; j<nrow_LambdaMatrix; j++){
          if(IndepIndex==1){
            tmpchi[j] = rchisq(1);
          }else{
            if(tt==0){tmpchi[j] = rchisq(1);}
          }
          tmp = tmp + LambdaMatrix[j][tt]*tmpchi[j];
        }
        if(tmp >= Cobs){
          Count = Count + 1;
          break;
        }
      }
    }
    vals[1] = Count/nsimchi;

    PutRNGstate();
 }

 
