
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
#include "MonteCarloPval.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void MonteCarloPval(double* RtL, int* RsqAIndex_i, int* RsqAIndex_j, double* RsqAval,int* dims, double* Cobs, int* nsimU, double* pval, int* Nonzero_sqA)
{
    int nrow_tL =dims[1],  len_sqAmatrix = dims[2], MaxCount = dims[3];
    
    double **tL, **sqAval;
    int **sqAIndex_i, **sqAIndex_j;
    
    /* reorganize vector into matrix */
    reorg(RtL, &tL, nrow_tL, nrow_tL);
    reorg(RsqAval, &sqAval, len_sqAmatrix, MaxCount);
    
    reorg_int(RsqAIndex_i, &sqAIndex_i, len_sqAmatrix, MaxCount);
    reorg_int(RsqAIndex_j, &sqAIndex_j, len_sqAmatrix, MaxCount);
    
    
    
    //Rprintf("test=%d\n",0);
    MonteCarloPvalC(tL, sqAIndex_i, sqAIndex_j, sqAval, dims, Cobs, nsimU, pval, Nonzero_sqA);
} 


void MonteCarloPvalC(double** tL, int** sqAIndex_i, int** sqAIndex_j, double** sqAval, int* dims, double* Cobs, int* nsimU, double* pval, int* Nonzero_sqA)
{
  
    int nSimu =  dims[0], nrow_tL = dims[1],  len_sqAmatrix = dims[2], MaxCount = dims[3], nrep = dims[4], Maxcol = dims[5];
    int ii, jj, iL, sim, subset, indexij, j;
    double Count = 0.0, tmp, *tmpnorm, *simZtL, *Cvec, **simZtLAq, simTotal;
    
    tmpnorm = (double *)calloc(nrow_tL, sizeof(double));
    simZtL = (double *)calloc(nrow_tL, sizeof(double));
    
    Cvec = (double *)calloc(len_sqAmatrix, sizeof(double));
    
    simZtLAq = (double **) malloc(len_sqAmatrix * sizeof(double*));
    
    simZtLAq[0] = (double *) calloc(len_sqAmatrix*nrow_tL, sizeof(double));
    if(simZtLAq[0] == NULL){ error("fail to allocate memory of simZtLAq"); }
    for(j=0; j<len_sqAmatrix; j++){
      simZtLAq[j] = simZtLAq[0] + j*nrow_tL;
    }

    //Rprintf("test=%d\n",1);
    GetRNGstate();
    Count = 0.0;
    simTotal = 0.0;
    for(ii = 0; ii < nrep; ii++){
      for(sim = 0; sim < nsimU[ii]; sim++){
        simTotal = simTotal + 1;
        //Rprintf("test=%d\n",2);
        for(j = 0; j < nrow_tL; j++){
          tmpnorm[j] = rnorm(0,1);
        }
        for(iL = 0; iL < nrow_tL; iL++){
          tmp = 0;
          //Rprintf("test=%d\n",3);
          for(j = 0; j < nrow_tL; j++){
            tmp = tmp + tmpnorm[j]*tL[j][iL];
          }
          simZtL[iL] = tmp;
        }
        //Rprintf("test=%d\n",4);
        for(subset = 0; subset < len_sqAmatrix; subset++){
          Cvec[subset] = 0;
          for(j = 0; j < nrow_tL; j++){
            simZtLAq[subset][j] = 0;
          }
          //Rprintf("test=%d\n",5);
          for(indexij = 0; indexij < Nonzero_sqA[subset]; indexij++){
            simZtLAq[subset][sqAIndex_j[subset][indexij]] = simZtLAq[subset][sqAIndex_j[subset][indexij]] + simZtL[sqAIndex_i[subset][indexij]]*sqAval[subset][indexij];
          }
          //Rprintf("test=%d\n",6);
          for(j = 0; j < nrow_tL; j++){
            Cvec[subset] = Cvec[subset] + simZtLAq[subset][j]*simZtLAq[subset][j];
          }
          //Rprintf("test=%d\n",7);
          if(Cvec[subset] >= Cobs[0]){
            Count = Count + 1;
            break;
          }
        }
      }
      if(Count > 0){
        break;
      }
    }
   
    pval[0] = Count/simTotal;

    PutRNGstate();
    
    free(simZtLAq[0]);
    free(simZtLAq);
    free(tmpnorm);
    free(simZtL);
    free(Cvec);
 }

 
