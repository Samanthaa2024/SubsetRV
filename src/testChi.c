
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



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void testChiC(double* Chioutput,int* dims)
{
    int nchi = dims[0], j;
    GetRNGstate();
    for(j=0; j<nchi; j++){
      Chioutput[j] = rchisq(1);
    }
    PutRNGstate();
} 

