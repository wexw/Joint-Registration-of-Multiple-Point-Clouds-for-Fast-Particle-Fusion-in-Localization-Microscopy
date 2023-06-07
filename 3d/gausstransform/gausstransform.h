/*
 * (C) Copyright 2018-2020      
 * Faculty of Applied Sciences
 * Delft University of Technology
 *
 * Ben van Werkhoven, November 2020.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 */
 
#pragma once

#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT 
#endif

#include <stdint.h>
#include <cuda_runtime.h>

class GPUGaussTransform {
    public:
        DLL_EXPORT GPUGaussTransform(int max_n, int argdim);
        ~GPUGaussTransform();
        DLL_EXPORT double compute(const double *A, const double *B, int m, int n, double scale, double *grad);
        int dim;
        int max_n;
    private:
        double *d_A;
        double *d_B;
        double *d_grad;
        double *d_cross_term;

        cudaStream_t stream;
};


