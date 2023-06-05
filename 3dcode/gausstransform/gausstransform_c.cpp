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
 
#include <math.h>

//definition
#include "gausstransform_ref.h"



//use openmp just for measuring time
#include <omp.h>


//extern C interface
extern "C" {

float call_GaussTransform(double* cost, const double* A, const double* B, int m, int n, int dim, double 
        scale, double* grad) {
    *cost = GaussTransform(A, B, m, n, dim, scale, grad);
    return 0.0;
}


float time_GaussTransform(double* cost, const double* A, const double* B, int m, int n, int dim, double 
        scale, double* grad) {
    double start = omp_get_wtime();
    *cost = GaussTransform(A, B, m, n, dim, scale, grad);
    return (float)((omp_get_wtime() - start)*1e3);
}

}

