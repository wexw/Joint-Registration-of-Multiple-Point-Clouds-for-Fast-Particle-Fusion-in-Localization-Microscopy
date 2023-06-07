#include <stdio.h>

/*%%=====================================================================
 * %% Project:   Pointset Registration using Gaussian Mixture Model
 * %% Module:    $RCSfile: mex_GaussTransform.c,v $
 * %% Language:  C
 * %% Author:    $Author: bing.jian $
 * %% Date:      $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
 * %% Version:   $Revision: 109 $
 * %%=====================================================================*/

#include "mex.h"
#include "expdist.h"

GPUExpDist *gpu_expdist = (GPUExpDist *)NULL;

/*
 * For 2D this function can be called as:
 *      mex_expdist(Mt.points, S.points, Mt.sigma, S.sigma);
 * For 3D this function can be called using:
 *      mex_expdist(Mt.points, S.points, Mt.sigma, S.sigma, rotation_matrix);
 *   Where:
 *   * The points arrays are Nx2 arrays, which have a data layout of:
 *      x1,x2,x3,...xn followed by y1,y2,y3,...,yn
 *   * The sigma arrays are 2xN arrays for particles of size N with a data layout of:
 *      sxy1,sz1,sxy2,sz2,sxy3,sz3,...,sxyn,szn
 *   * The rotation_matrix is a 3x3 array of doubles, assumed to be stored in row-major
 *      for example,    r11, r12, r13
 *                      r21, r22, r23
 *                      r31, r23, r33 
 *      is assumed to be stored as:
 *           r11,r12,r13,r21,r22,r23,r31,r23,r33
 */
void mexFunction(int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    int m, n, dim; 
    int use_prerotated_scale_B = 0;
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

        if (mxGetN(prhs[3])==2 || mxGetN(prhs[3])==9)
        {
            if (mxGetN(prhs[3])==9) {
                use_prerotated_scale_B = 1;
            }
        } else {
            mexErrMsgTxt("uncertainties should be an n*2 or n*9 matrix for 3D.");
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

    if (gpu_expdist == (GPUExpDist *)NULL) {
        gpu_expdist = new GPUExpDist(100000, dim);
    } else {
    if (gpu_expdist->dim != dim) {
            mexErrMsgTxt("Error GPUExpDist dimension does not match");
        }
        if ((gpu_expdist->max_n < m) || (gpu_expdist->max_n < n)) {
            mexErrMsgTxt("Error GPUExpDist not sufficiently large");
        }
    }

    if (dim == 2) {
        *result = gpu_expdist->compute(A, B, m, n, scale_A, scale_B);
    } else if (dim == 3) {
        *result = gpu_expdist->compute(A, B, m, n, scale_A, scale_B, 
                                       (double *)mxGetPr(prhs[4]), use_prerotated_scale_B);
    }


}

