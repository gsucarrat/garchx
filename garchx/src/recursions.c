# include <R.h>
# include <math.h>
# include "recursions.h"

void GARCHXRECURSION( int *iStart, int *iEnd, int *iGARCHorder, double *sigma2, double *parsgarch, double *innov ){

  double garchsum; 
  for(int i=*iStart; i < *iEnd; i++){

    /* GARCH sum */
    garchsum = 0;
    for(int j=0; j < *iGARCHorder; j++){
      garchsum = garchsum + parsgarch[j] * sigma2[i-1-j];
    }

    /* recursion */
    sigma2[i] =  garchsum + innov[i];

  }

}

void GARCHXRECURSIONSIM(int * iStart, int * iEnd, int * iARCHorder, int * iGARCHorder, int * iASYMorder, double * parsarch, double * parsgarch, double * parsasym, double * sigma2, double * z2, double * Ineg, double * xregsum ){
	
  double archsum;
  double garchsum;
  double asymsum;
  
  archsum = 0;
  garchsum = 0;
  asymsum = 0;

  for(int i=*iStart; i < *iEnd; i++){

    /* ARCH sum */
    if(iARCHorder != 0){
      archsum = 0;
      for(int j=0; j < *iARCHorder; j++){
        archsum = archsum + parsarch[j] * z2[i-1-j] * sigma2[i-1-j];
      }
    }

    /* GARCH sum */
    if(iGARCHorder != 0){
      garchsum = 0;
      for(int j=0; j < *iGARCHorder; j++){
        garchsum = garchsum + parsgarch[j] * sigma2[i-1-j];
      }
    }
    
    /* ASYM sum */
    if(iASYMorder != 0){
      asymsum = 0;
      for(int j=0; j < *iASYMorder; j++){
        asymsum = asymsum + parsasym[j] * Ineg[i-1-j] * z2[i-1-j] * sigma2[i-1-j];
      }
    }
    
    /*recursion:*/
    sigma2[i] = archsum + garchsum + asymsum + xregsum[i];
	
	} /* close for loop */

}
