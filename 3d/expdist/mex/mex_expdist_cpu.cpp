#include <stdio.h>
#include <cmath>

/*%%=====================================================================
 * %% Project:   Pointset Registration using Gaussian Mixture Model
 * %% Module:    $RCSfile: mex_GaussTransform.c,v $
 * %% Language:  C
 * %% Author:    $Author: bing.jian $
 * %% Date:      $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
 * %% Version:   $Revision: 109 $
 * %%=====================================================================*/

#include "mex.h"

#include "expdist_ref.h"


void mexFunction(int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    int m, n, dim; 
    double *A, *B, *result, *scale_A, *scale_B;

    //prevent clearing from memory
    void mexLock(void);
    
    /* Check for proper number of input and output arguments */    
    if (nlhs > 2){
    mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument */
    if (!(mxIsDouble(prhs[0]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[1]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[2]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[3]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }     
     
    /* Get the number of elements in the input argument */
    /* elements=mxGetNumberOfElements(prhs[0]); */
    /* Get the data */
    A = (double *)mxGetPr(prhs[0]);
    B = (double *)mxGetPr(prhs[1]);
    scale_A = (double *)mxGetPr(prhs[2]);
    scale_B = (double *)mxGetPr(prhs[3]);
      /* Get the dimensions of the matrix input A&B. */
      m = static_cast<int>(mxGetM(prhs[0]));
      n = static_cast<int>(mxGetM(prhs[1]));

      dim = static_cast<int>(mxGetN(prhs[0]));

      if (mxGetN(prhs[1])!=dim)
      {
          mexErrMsgTxt("The two input point sets should have same dimension.");
      }
    
      if (mxGetM(prhs[0])!=mxGetM(prhs[2]))
      {
          mexErrMsgTxt("localizations and uncertainties (A) should have same dimension.");
      }    
    
      if (mxGetM(prhs[1])!=mxGetM(prhs[3]))
      {
          mexErrMsgTxt("localizations and uncertainties (B) should have same dimension.");
      }        
    
    if (dim == 2) {
        if (nrhs != 4) {
            mexErrMsgTxt("Four input arguments required.");
        } 
          if (mxGetN(prhs[2])!=1)
          {
              mexErrMsgTxt("uncertainties should be a m*1 matrix.");
          }
    
          if (mxGetN(prhs[3])!=1)
          {
              mexErrMsgTxt("uncertainties should be a n*1 matrix.");
          }
    }
    else if (dim == 3) {
        if (nrhs != 5) {
            mexErrMsgTxt("Five input arguments required.");
        } 
          if (mxGetN(prhs[2])!=2)
          {
              mexErrMsgTxt("uncertainties should be a m*2 matrix for 3D.");
          }
    
          if (mxGetN(prhs[3])!=2)
          {
              mexErrMsgTxt("uncertainties should be a n*2 matrix for 3D.");
          }

        //rotation matrix
        if (!(mxIsDouble(prhs[4]))) {
          mexErrMsgTxt("Input array must be of type double.");
        }     
          if (mxGetN(prhs[4])!=3 || mxGetM(prhs[4])!=3)
          {
              mexErrMsgTxt("rotation matrix should be 3x3 for 3D.");
          }
    }

    /* Allocate the space for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    result = mxGetPr(plhs[0]);

    if (dim == 3) {
        double *rotated_scales = (double *)malloc(9*n*sizeof(double));
        //double *rotated_B = (double *)malloc(3 * n * sizeof(double));
        double *rotation_matrix = (double *)mxGetPr(prhs[4]);
        rotate_scales(rotated_scales, rotation_matrix, n, scale_B);
        //rotate_B(rotated_B, rotation_matrix, n, B);  //rotation now happens before B is passed to this function
        *result = expdist3D(A, B, m, n, scale_A, rotated_scales);
        //free(rotated_B);
        free(rotated_scales);
    } else {
        //*result = expdist(A, B, m, n, dim, scale_A, scale_B); // not used in alltoall3d
        mexErrMsgTxt("wrong number of dimensions.");
    }
}

