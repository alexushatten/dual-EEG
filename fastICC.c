#include "mex.h"
#include "math.h"
/************************************************************/
/*                                                          */
/* This function computes the icc between vectors x and y.  */
/*                                                          */
/************************************************************/
/*                                                          */
/*==========================================================*/
/*                                                          */
/*     Renaud Marquis @FBMlab, June 2018                    */
/*     based on Liber Eleutherios's fastcorr.c              */
/*     and Arash Salarian's ICC.m                           */
/*                                                          */
/*==========================================================*/
/*                                                          */
/************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int i;
double mean_x, mean_y, sweep, delta_x, delta_y, SStotal, MSC, MSR, MSE, ICC, delta_row, sum_sq_row, mean_row, delta_xx, delta_yy;
double *x, *y;
int dimx, dimy;
double *out;

dimx = mxGetNumberOfElements(prhs[0]);
dimy = mxGetNumberOfElements(prhs[1]);

if((dimx == dimy) && (dimx > 1)){

  x = (double *) mxGetPr(prhs[0]);
  y = (double *) mxGetPr(prhs[1]);

  mean_x = *x;
  mean_y = *y;
  mean_row = (*x + *y) / 2;
  sum_sq_row = 0;

  for (i = 1; i < dimx; i++){
    sweep = (double) (i + 0.0) / (i + 1.0);
    delta_x = *(x + i) - mean_x; // pointer arithmetics !
    delta_y = *(y + i) - mean_y;
    mean_x += delta_x / (i + 1.0);
    mean_y += delta_y / (i + 1.0);
    delta_row = ((*(x + i) + *(y + i)) / 2) - mean_row; // OK
    
    sum_sq_row += delta_row * delta_row * sweep; // OK
    mean_row += delta_row / (i + 1.0); // aka global mean, OK
  }

  SStotal = ((*x - mean_row) * (*x - mean_row)) + ((*y - mean_row) * (*y - mean_row));
  
  for (i = 1; i < dimx; i++){
    delta_xx = *(x + i) - mean_row;
    delta_yy = *(y + i) - mean_row;
    SStotal += (delta_xx * delta_xx + delta_yy * delta_yy); // OK
  }

  MSC = ((mean_x - mean_row) * (mean_x - mean_row) + (mean_y - mean_row) * (mean_y - mean_row)) * dimx; // OK
  MSR = sum_sq_row / (dimx - 1) * 2; // OK
  MSE = (SStotal - MSR * (dimx - 1) - MSC * (2 -1) ) / ( (dimx - 1) * (2 - 1) );
  ICC = (MSR - MSE) / (MSR + (2 - 1) * MSE);
}

/***************************************************************************************/
/***************************************************************************************/
plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
out = mxGetPr(plhs[0]);
*(out) = ICC;

// plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
// out = mxGetPr(plhs[1]);
// *(out) = MSE;

}
