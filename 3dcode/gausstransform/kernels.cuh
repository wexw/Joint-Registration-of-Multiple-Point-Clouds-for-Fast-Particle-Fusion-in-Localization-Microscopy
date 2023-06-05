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
 
#include <cub/cub.cuh>

#ifndef block_size_x
    #define block_size_x 64
#endif

#ifndef tile_size_x
    #define tile_size_x 2
#endif

#ifndef use_registers_B
    #define use_registers_B 1
#endif

#ifndef cub_algorithm
    #define cub_algorithm BLOCK_REDUCE_WARP_REDUCTIONS
#endif

#define use_cub_algorithm cub::cub_algorithm


template<typename T, int dim>
__device__ __forceinline__ T compute_grad_and_cross(const int i, const int j, const T (&l_A)[dim], const T *B,
        const T scale_sq, T (&grad_i)[dim], const int m, const int n) {

    T cost_ij = 0.0;

    if (i<m) {

            T dist_ij = 0;
            #pragma unroll
            for (int d = 0; d < dim; ++d) {
                dist_ij += (l_A[d] - B[d])*(l_A[d] - B[d]);
            }
            cost_ij = exp(-dist_ij/scale_sq);

            #pragma unroll
            for (int d = 0; d < dim; ++d) {
                grad_i[d] -= cost_ij * 2.0 * (l_A[d] - B[d]);
            }

    }

    return cost_ij;
}



/*
 * This function performs the main body of work for computing the Gauss transform
 * The parallelization is such that one thread block is created
 * for each item in A, which is of size m. This implies that each thread block
 * does n (size of B) work.
 * The gradient computed in this function is reduced to a single value within the
 * thread block. The same is done for the cross term, which then needs to be
 * reduced in a second kernel. 
 */
template<typename T, int dim>
__device__ __forceinline__ void GaussTransform_blocked_i(const T *A, const T *B,
                const int m, const int n, const T scale_sq, T *d_grad, T *d_cross_term) {

    int tx = threadIdx.x;

    // Specialize BlockReduce for a 1D block of block_size_x threads on type T
    typedef cub::BlockReduce<T, block_size_x, use_cub_algorithm> BlockReduce;
    // Allocate shared memory for BlockReduce
    __shared__ typename BlockReduce::TempStorage temp_storage;

    int i = blockIdx.x*tile_size_x; // i-loop is parallelized over thread blocks with tile_size_x stride

    T cross_term = 0.0;
    T grad_i[tile_size_x][dim];
    T l_A[tile_size_x][dim];
    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {
        #pragma unroll
        for (int d = 0; d < dim; d++) {
            l_A[ti][d] = A[(i+ti)*dim+d];
            grad_i[ti][d] = 0.0;
        }
    }

    #if use_registers_B == 1
    T l_B[dim];
    #endif

    //j-loop parallelized over threads within thread block
    for (int j = tx; j<n; j+=block_size_x) {

        #if use_registers_B == 1
        #pragma unroll
        for (int d = 0; d < dim; d++) {
            l_B[d] = B[j*dim+d];
        }
        #endif

        //#pragma unroll
        for (int ti=0; ti<tile_size_x; ti++) {
            #if use_registers_B == 1
            cross_term += compute_grad_and_cross<double, dim>(i+ti, j, l_A[ti], l_B, scale_sq, grad_i[ti], m, n);
            #else
            cross_term += compute_grad_and_cross<double, dim>(i+ti, j, l_A[ti], B+j*dim, scale_sq, grad_i[ti], m, n);
            #endif
        }

    }

    //reduce grad_i for each d, within the block
    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {
        #pragma unroll
        for (int d = 0; d < dim; d++) {
            grad_i[ti][d] = BlockReduce(temp_storage).Sum(grad_i[ti][d]);
            __syncthreads();
        }
    }

    //reduce cross_term within the block, (division by m*n on CPU)
    cross_term = BlockReduce(temp_storage).Sum(cross_term);

    if (tx == 0) {
        #pragma unroll
        for (int ti=0; ti<tile_size_x; ti++) {
            if (i+ti < m) {
                #pragma unroll
                for (int d = 0; d < dim; d++) {
                    d_grad[(i+ti) * dim + d] = grad_i[ti][d] / (scale_sq * m * n);
                }
            }
        }
        d_cross_term[blockIdx.x] = cross_term;
    }
}


extern "C"
__global__ void
GaussTransform(const double* A, const double* B,
    int m, int n, double scale_sq, double *grad, double *cross_term);


extern "C"
__global__ void
GaussTransform3D(const double* A, const double* B,
    int m, int n, double scale_sq, double *grad, double *cross_term);



/*
 * Reduce the per thread block cross terms computed in the GaussTransform kernel to single value
 * and divide by (m*n)
 *
 * This kernel is designed to run as single-thread block, because the number of terms to reduce is
 * of size n or m, which is expected to be around 2000 or so. The number of items to reduce
 * is passed as the last argument 'nblocks', which corresponds to the number of thread blocks used
 * by the first kernel.
 */
extern "C"
__global__ void reduce_cross_term(double *output, double *d_cross_term, const int m, const int n, const int nblocks);


