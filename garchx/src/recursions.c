# include <R.h>
# include <math.h>
# include "recursions.h"

void GARCHXRECURSION( int *iStart, int *iEnd, int *iGARCHorder, double *sigma2, double *parsgarch, double *innov ){

  double garchsum; 
  for(int i=*iStart; i < *iEnd; i++){

    /*GARCH sum:*/
    garchsum = 0;
    for(int j=0; j < *iGARCHorder; j++){
      garchsum = garchsum + parsgarch[j] * sigma2[i-1-j];
    }

    /*recursion:*/
    sigma2[i] =  garchsum + innov[i];

  }

}
