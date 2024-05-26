
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
#include "PermuPval.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void PermuPval(double* RG1, double* RG2, double* RRes, int* dims, double* RQMatrix, int* RPermuteIndex)
{
    int nrow_Res = dims[1], ncol_Res = dims[2], nrow_G = dims[3], ncol_G = dims[4];
    int nrow_QMatrix = dims[0], ncol_QMatrix = dims[5], nrow_PermuteIndex = dims[6], ncol_PermuteIndex = dims[7];
    int **PermuteIndex;
    double **Res, **G1, **G2, **QMatrix;
    
    /* reorganize vector into matrix */
    reorg(RRes, &Res, nrow_Res, ncol_Res);
    reorg(RG1, &G1, nrow_G, ncol_G);
    reorg(RG2, &G2, nrow_G, ncol_G);
    reorg(RQMatrix, &QMatrix, nrow_QMatrix, ncol_QMatrix);
    
    reorg_int(RPermuteIndex, &PermuteIndex, nrow_PermuteIndex, ncol_PermuteIndex);
    PermuPvalC(G1, G2, Res, dims, QMatrix, PermuteIndex);
} 


void PermuPvalC(double** G1, double** G2, double** Res, int* dims, double** QMatrix, int** PermuteIndex)
{
  
    int nperm = dims[0], N = dims[1], ncol_Res = dims[2], nrow_G = dims[3], ncol_G = dims[4];
    int nrow_QMatrix = dims[0], ncol_QMatrix = dims[5];;
    int i, j, temp_nperm = 0, nrow_PermuteIndex = dims[6], ncol_PermuteIndex = dims[7];
    double tmp;
    
   
    const int N1 = N;
    int Value[N1];

    GetRNGstate();

    for (i = 0; i < N1; i++) {
      Value[i] = 0;
    }
    //printf("start=%d\n", 0);
    //printf("N1=%d\n", N1);
    if(nrow_PermuteIndex==1){
      permute(Value, N1, 0, temp_nperm, QMatrix, nperm, G1, G2, Res, ncol_G);
    }else{
      double tmp1, tmp2, Q1val, Q2val, Q3val;
      for(temp_nperm = 0; temp_nperm < nperm; temp_nperm++){
      
        Q1val = 0;
        Q2val = 0;
        Q3val = 0;
        for(j = 0; j < ncol_G; j++){
          tmp1 = 0;
          tmp2 = 0;
          for (i = 0; i < N; i++) {
            tmp1 = tmp1 + Res[PermuteIndex[temp_nperm][i]][0]*G1[i][j];
            tmp2 = tmp2 + Res[PermuteIndex[temp_nperm][i]][1]*G2[i][j];
          }
          Q1val = Q1val + tmp1*tmp1;
          Q2val = Q2val + tmp2*tmp2;
          Q3val = Q3val + (tmp1 + tmp2)*(tmp1 + tmp2);
        }
        QMatrix[temp_nperm][0] = Q1val;
        QMatrix[temp_nperm][1] = Q2val;
        QMatrix[temp_nperm][2] = Q3val;
      }
    }
    //visit_test(Value, N1, 0);

    PutRNGstate();
 }

 
