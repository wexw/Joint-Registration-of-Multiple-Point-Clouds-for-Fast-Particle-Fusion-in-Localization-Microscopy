#include <stdio.h>

/*%%=====================================================================
%% Project:   Pointset Registration using Gaussian Mixture Model
%% Module:    $RCSfile: mex_GaussTransform.c,v $
%% Language:  C
%% Author:    $Author: bing.jian $
%% Date:      $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% Version:   $Revision: 109 $
%%=====================================================================*/

#include "mex.h"
#include "gausstransform.h"

GPUGaussTransform *gpu_gt = (GPUGaussTransform *)NULL;

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    int m, n, dim; 
    double *A, *B, *result, *grad, scale;

    //prevent clearing from memory
    void mexLock(void);
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 3) {
	mexErrMsgTxt("Three input arguments required.");
    } 
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

    /* Get the number of elements in the input argument */
    /* elements=mxGetNumberOfElements(prhs[0]); */
    /* Get the data */
    A = (double *)mxGetPr(prhs[0]);
    B = (double *)mxGetPr(prhs[1]);
    scale = mxGetScalar(prhs[2]);
    /* Get the dimensions of the matrix input A&B. */
    m = static_cast<int>(mxGetN(prhs[0]));
    n = static_cast<int>(mxGetN(prhs[1]));
    dim = static_cast<int>(mxGetM(prhs[0]));
    if (mxGetM(prhs[1])!=dim)
    {
        mexErrMsgTxt("The two input point sets should have same dimension.");
    }
    /* Allocate the space for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(dim,m,mxREAL);
    result = mxGetPr(plhs[0]);
    grad = mxGetPr(plhs[1]);

    if (gpu_gt == (GPUGaussTransform *)NULL) {
        gpu_gt = new GPUGaussTransform(1000000, dim);
    } else {
        if (gpu_gt->dim != dim) {
            mexErrMsgTxt("Error GPUGaussTransform dimension does not match");
        }
        if ((gpu_gt->max_n < m) || (gpu_gt->max_n < n)) {
            mexErrMsgTxt("Error GPUGaussTransform not sufficiently large");
        }
    }

    *result = gpu_gt->compute(A, B, m, n, scale, grad);
}
