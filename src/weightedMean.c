#include <R.h>
/* Computes the weighted mean of the target.

   Input: 
   - target:   vector of length n with target values
   - w:        vector of length n with nonnegative weights
   - iord:     vector of length n that can be used as index to order the variables,
               missing values are in the last tieblock
   - ties:     vector of lengths nties that elements giving the length of the tie block
   - nties:    number of ties
   
   Output
   - y:        vector of length nties with the weighted mean vector y
   - sumwvec:  vector of length nties with the sum of the weights
   
   Author: Patrick Groenen
   Date:   June 5, 2014
   Version: 1.0
   
   */
void weightedMean(double *y, double *sumwvec, double *target, double *w, 
                   int *iord, int* ties, int *n, int *nties)
{
  int i, k, l, ind, nprevties;
  double sumw, sumwt;

  i = 0;
  nprevties = 0;
	for (k = 0; k<*nties; k++) {
    sumwt = 0;
    sumw  = 0;   
    for (l = 0; l < ties[k]; l++){
      ind = iord[nprevties + l] - 1;
      sumwt += w[ind]*target[ind];
      sumw  += w[ind];
    }
    if (sumw > 1E-10){
      y[k] = sumwt / sumw;
    } else {
      y[k] = 0;
    }
    
    sumwvec[k] = sumw;
    nprevties  += ties[k]; 
	}
}


