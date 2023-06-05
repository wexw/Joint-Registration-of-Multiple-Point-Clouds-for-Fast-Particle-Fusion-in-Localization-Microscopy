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

__constant__ double rotation_matrixd[9];
__constant__ double rotation_matrix_transposedd[9];

//defaults tuned for K40
#ifndef block_size_x
    #define block_size_x 32
#endif
#ifndef block_size_y
    #define block_size_y 4
#endif

#ifndef tile_size_x
    #define tile_size_x 2
#endif
#ifndef tile_size_y
    #define tile_size_y 4
#endif

#ifndef use_shared_mem
    #define use_shared_mem 1
#endif

#include "expdist_functions.cuh"


template <int tile_size, int stride, typename T>
__device__ __forceinline__ void fill_shared_mem_tiled_1D(T (&sh_mem)[tile_size*stride], const T *d_mem, int sh_offset, int d_offset, int N) {
    #pragma unroll
    for (int ti=0; ti<tile_size; ti++) {
        if (d_offset+ti*stride < N) {
            sh_mem[sh_offset+ti*stride] = d_mem[d_offset+ti*stride];
        }
    }
}



/*
 * This function performs the main body of work for computing the Bhattacharya
 * cost function for two given point sets.
 * The parallelization is such that a 2D array of 2D thread blocks are created
 * to match the m*n iteration space. The amount of work per thread is controlled
 * through tiling factors tile_size_x and tile_size_y.
 * The cross term is reduced to a single value per thread block, which then needs
 * to be reduced to a single value in a second kernel. 
 */
template<typename T, int dim>
__device__ __forceinline__ void ExpDist_tiled(const T *A, const T *B,
                 int m, int n, const T *scale_A, const T *scale_B, T *d_cross_term) {

    // Specialize BlockReduce for a 1D block of block_size_x threads on type T
    typedef cub::BlockReduce<T, block_size_x, cub::BLOCK_REDUCE_WARP_REDUCTIONS, block_size_y> BlockReduce;

    // Allocate shared memory for BlockReduce
    __shared__ typename BlockReduce::TempStorage temp_storage;

    T cross_term = 0.0;
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int i = tx + blockIdx.x * block_size_x * tile_size_x;
    int j = ty + blockIdx.y * block_size_y * tile_size_y;

    #if use_shared_mem == 1
    __shared__ T sh_A[dim][block_size_x*tile_size_x];
    __shared__ T sh_B[dim][block_size_y*tile_size_y];
    __shared__ T sh_scale_A[block_size_x*tile_size_x];
    __shared__ T sh_scale_B[block_size_y*tile_size_y];

    if (ty == 0) {
        #pragma unroll
        for (int d=0; d<dim; d++) {
            fill_shared_mem_tiled_1D<tile_size_x, block_size_x>(sh_A[d], A+d*m, tx, i, m);
        }
        fill_shared_mem_tiled_1D<tile_size_x, block_size_x>(sh_scale_A, scale_A, tx, i, m);
    }
    if (tx == 0) {
        #pragma unroll
        for (int d=0; d<dim; d++) {
            fill_shared_mem_tiled_1D<tile_size_y, block_size_y>(sh_B[d], B+d*n, ty, j, n);
        }
        fill_shared_mem_tiled_1D<tile_size_y, block_size_y>(sh_scale_B, scale_B, ty, j, n);
    }
    __syncthreads();
    #endif

    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {
        #pragma unroll
        for (int tj=0; tj<tile_size_y; tj++) {

            if ((i+ti*block_size_x < m) && (j+tj*block_size_y < n)) {

                T dist_ij = 0;

                #if use_shared_mem == 0
                    #pragma unroll
                    for (int d=0; d<dim; d++) {
                        int id = (i+ti*block_size_x) + d * m;
                        int jd = (j+tj*block_size_y) + d * n;
                        dist_ij += (A[id]-B[jd])*(A[id]-B[jd]);
                    }
                    cross_term += exp(-dist_ij/(scale_A[i+ti*block_size_x] + scale_B[j+tj*block_size_y]));
                #elif use_shared_mem == 1
                    #pragma unroll
                    for (int d=0; d<dim; d++) {
                        dist_ij += (sh_A[d][tx+ti*block_size_x]-sh_B[d][ty+tj*block_size_y])*(sh_A[d][tx+ti*block_size_x]-sh_B[d][ty+tj*block_size_y]);
                    }
                    cross_term += exp(-dist_ij/(sh_scale_A[tx+ti*block_size_x] + sh_scale_B[ty+tj*block_size_y]));
                #endif

            }
        }
    }

    //reduce cross_term within the block
    cross_term = BlockReduce(temp_storage).Sum(cross_term);

    //write back the per-thread block partial cross term
    if (tx == 0 && ty == 0) {
        d_cross_term[blockIdx.y*gridDim.x+blockIdx.x] = cross_term;
    }
}


/*
 * This function performs the main body of work for computing the Bhattacharya
 * cost function for two given point sets in 3D coordinates.
 * The parallelization is such that a 2D array of 2D thread blocks are created
 * to match the m*n iteration space. The amount of work per thread is controlled
 * through tiling factors tile_size_x and tile_size_y.
 * The cross term is reduced to a single value per thread block, which then needs
 * to be reduced to a single value in a second kernel. 
 */
template<typename T, int dim>
__device__ __forceinline__ void ExpDist3D_tiled(const T *A, const T *B,
                 int m, int n, const T *scale_A, const T *scale_B, T *d_cross_term) {

    // Specialize BlockReduce for a 2D block of block_size_x,block_size_y threads on type T
    typedef cub::BlockReduce<T, block_size_x, cub::BLOCK_REDUCE_WARP_REDUCTIONS, block_size_y> BlockReduce;

    // Allocate shared memory for BlockReduce
    __shared__ typename BlockReduce::TempStorage temp_storage;

    T cross_term = 0.0;
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int i = tx + blockIdx.x * block_size_x * tile_size_x;
    int j = ty + blockIdx.y * block_size_y * tile_size_y;

    #if use_shared_mem == 1
    __shared__ T sh_A[dim][block_size_x*tile_size_x];
    __shared__ T sh_B[dim][block_size_y*tile_size_y];
    __shared__ T sh_scale_A[block_size_x*tile_size_x][2];   // could be a number of 2-way bank conflicts here
    __shared__ T sh_scale_B[block_size_y*tile_size_y][9];   // 9 is a generator of 32, so no bank conflicts

    if (ty == 0) {
        #pragma unroll
        for (int d=0; d<dim; d++) {
            fill_shared_mem_tiled_1D<tile_size_x, block_size_x>(sh_A[d], A+d*m, tx, i, m);
        }
        #pragma unroll
        for (int ti=0; ti<tile_size_x; ti++) {
            #pragma unroll
            for (int d=0; d<2; d++) {
                sh_scale_A[tx+ti*block_size_x][d] = scale_A[(i+ti*block_size_x)*2+d];  //reading from global mem with stride 2 here
            }
        }
    }

    if (tx == 0) {  // a column of threads accessing memory here is not very efficient, should remap to horizontal group of threads of same size
        #pragma unroll
        for (int d=0; d<dim; d++) {
            fill_shared_mem_tiled_1D<tile_size_y, block_size_y>(sh_B[d], B+d*n, ty, j, n);
        }
        #pragma unroll
        for (int tj=0; tj<tile_size_y; tj++) {
            #pragma unroll
            for (int d=0; d<9; d++) {
                sh_scale_B[ty+tj*block_size_y][d] = scale_B[(j+tj*block_size_y)*9+d];  //reading from global mem with stride 9 here
            }
        }
    }

    __syncthreads();
    #endif

    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {

        if (i+ti*block_size_x < m) {

        //register version of shared memory
        T pA[3];
        #pragma unroll
        for (int d=0; d<3; d++) {
            #if use_shared_mem == 1
            pA[d] = sh_A[d][tx+ti*block_size_x];
            #else //use_shared_mem == 0
            pA[d] = A[i+ti*block_size_x+d*m];
            #endif
        }

        T Sigma_i[9];
        zero_matrix(Sigma_i);

        #if use_shared_mem == 0
        Sigma_i[0] = scale_A[i*2+0];   // 0 1 2
        Sigma_i[4] = scale_A[i*2+0];   // 3 4 5
        Sigma_i[8] = scale_A[i*2+1];   // 6 7 8
        #else
        Sigma_i[0] = sh_scale_A[tx+ti*block_size_x][0];   // 0 1 2
        Sigma_i[4] = sh_scale_A[tx+ti*block_size_x][0];   // 3 4 5
        Sigma_i[8] = sh_scale_A[tx+ti*block_size_x][1];   // 6 7 8
        #endif

        #pragma unroll
        for (int tj=0; tj<tile_size_y; tj++) {

            if (j+tj*block_size_y < n) {

                T pB[3];
                #pragma unroll
                for (int d=0; d<3; d++) {
                    #if use_shared_mem == 1
                    pB[d] = sh_B[d][ty+tj*block_size_y];
                    #else //use_shared_mem == 0
                    pB[d] = B[j+tj*block_size_y+d*n];
                    #endif
                }

                T Sigma_j[9];
                #if use_shared_mem == 1
                load_matrix(Sigma_j, sh_scale_B[ty+tj*block_size_y], 0);
                #else
                load_matrix(Sigma_j, scale_B, j);
                #endif

                cross_term += compute_expdist_3D<T, 3>(pA, pB, Sigma_i, Sigma_j);
                
            }
        }

        } //end of bounds check on ti
    }

    //reduce cross_term within the block
    cross_term = BlockReduce(temp_storage).Sum(cross_term);

    //write back the per-thread block partial cross term
    if (tx == 0 && ty == 0) {
        d_cross_term[blockIdx.y*gridDim.x+blockIdx.x] = cross_term;
    }
}


extern "C"
__global__ void
ExpDist(const double *A, const double *B,
    int m, int n, const double *scale_A, const double *scale_B, double *cross_term);

extern "C"
__global__ void
ExpDist3D(const double *A, const double *B,
    int m, int n, const double *scale_A, const double *scale_B, double *cross_term);


template<typename T, int dim>
__device__ __forceinline__ T compute_expdist_block_shared(int tx, int ty, int i, int j,
                                                          T (&sh_A)[dim][block_size_x*tile_size_x], 
                                                          T (&sh_B)[dim][block_size_y*tile_size_y],
                                                          T (&sh_scale_A)[block_size_x*tile_size_x],
                                                          T (&sh_scale_B)[block_size_y*tile_size_y], int m, int n) {
    T cross_term = 0;

    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {
        #pragma unroll
        for (int tj=0; tj<tile_size_y; tj++) {

            if ((i+ti*block_size_x < m) && (j+tj*block_size_y < n)) {

                T dist_ij = 0;

                #pragma unroll
                for (int d=0; d<dim; d++) {
                     dist_ij += (sh_A[d][tx+ti*block_size_x]-sh_B[d][ty+tj*block_size_y])*(sh_A[d][tx+ti*block_size_x]-sh_B[d][ty+tj*block_size_y]);
                }
                cross_term += exp(-dist_ij/(sh_scale_A[tx+ti*block_size_x] + sh_scale_B[ty+tj*block_size_y]));

            }
        }
    }

    return cross_term;
}


template<typename T, int dim>
__device__ __forceinline__ T compute_expdist_block(int i, int j, const T *A, const T *B,
                                                   const T *scale_A, const T* scale_B, int m, int n) {

    T cross_term = 0;

    #pragma unroll
    for (int ti=0; ti<tile_size_x; ti++) {
        #pragma unroll
        for (int tj=0; tj<tile_size_y; tj++) {

            if ((i+ti*block_size_x < m) && (j+tj*block_size_y < n)) {

                T dist_ij = 0;

                #pragma unroll
                for (int d=0; d<dim; d++) {
                    int id = (i+ti*block_size_x) + d * m;
                    int jd = (j+tj*block_size_y) + d * n;
                    dist_ij += (A[id]-B[jd])*(A[id]-B[jd]);
                }
                cross_term += exp(-dist_ij/(scale_A[i+ti*block_size_x] + scale_B[j+tj*block_size_y]));

            }
        }
    }

    return cross_term;
}



/*
 * This function performs the main body of work for computing the Bhattacharya
 * cost function for two given point sets.
 * The parallelization is such that a 1D array of 2D thread blocks is created over
 * the m-dimension. The thread blocks then iterate over n, to process the entire
 * m*n iteration space. The amount of work per thread is controlled
 * through tiling factors tile_size_x and tile_size_y.
 * The cross term is reduced to a single value per thread block, which then needs
 * to be reduced to a single value in a second kernel. 
 */
template<typename T, int dim>
__device__ __forceinline__ void ExpDist_tiled_column(const T *A, const T *B,
                 int m, int n, const T *scale_A, const T *scale_B, T *d_cross_term) {

    // Specialize BlockReduce for a 1D block of block_size_x threads on type T
    typedef cub::BlockReduce<T, block_size_x, cub::BLOCK_REDUCE_WARP_REDUCTIONS, block_size_y> BlockReduce;

    // Allocate shared memory for BlockReduce
    __shared__ typename BlockReduce::TempStorage temp_storage;

    T cross_term = 0.0;
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int i = tx + blockIdx.x * block_size_x * tile_size_x;
    int j = ty + blockIdx.y * block_size_y * tile_size_y;

    #if use_shared_mem == 1
    __shared__ T sh_A[dim][block_size_x*tile_size_x];
    __shared__ T sh_B[dim][block_size_y*tile_size_y];
    __shared__ T sh_scale_A[block_size_x*tile_size_x];
    __shared__ T sh_scale_B[block_size_y*tile_size_y];

    #pragma unroll
    for (int d=0; d<dim; d++) {
        fill_shared_mem_tiled_1D<tile_size_x, block_size_x>(sh_A[d], A+d*m, tx, i, m);
    }
    fill_shared_mem_tiled_1D<tile_size_x, block_size_x>(sh_scale_A, scale_A, tx, i, m);
    __syncthreads();
    #endif

    int step_size = gridDim.y * block_size_y * tile_size_y;
    for (int sj = j; sj < n; sj += step_size) {

        #if use_shared_mem == 1
        for (int d=0; d<dim; d++) {
            fill_shared_mem_tiled_1D<tile_size_y, block_size_y>(sh_B[d], B+d*n, ty, sj, n);
        }
        fill_shared_mem_tiled_1D<tile_size_y, block_size_y>(sh_scale_B, scale_B, ty, sj, n);
        __syncthreads();
        #endif

        #if use_shared_mem == 0
            cross_term += compute_expdist_block<double, dim>(i, sj, A, B, scale_A, scale_B, m, n);
        #elif use_shared_mem == 1
            cross_term += compute_expdist_block_shared<double, dim>(tx, ty, i, sj, sh_A, sh_B, sh_scale_A, sh_scale_B, m, n);
            __syncthreads();
        #endif

    }

    //reduce cross_term within the block
    cross_term = BlockReduce(temp_storage).Sum(cross_term);

    //write back the per-thread block partial cross term
    if (tx == 0 && ty == 0) {
        d_cross_term[blockIdx.y*gridDim.x+blockIdx.x] = cross_term;
    }
}



#ifdef reduce_block_size
 #define block_size reduce_block_size
#else
 #define block_size block_size_x
#endif


extern "C"
__global__ void
ExpDist(const double *A, const double *B,
                 int m, int n, const double *scale_A, const double *scale_B, double *cross_term) {

    //2-dimensional with double precision
    ExpDist_tiled<double, 2>(A, B, m, n, scale_A, scale_B, cross_term);

}

extern "C"
__global__ void
ExpDist3D(const double *A, const double *B,
                 int m, int n, const double *scale_A, const double *scale_B, double *cross_term) {

    //3-dimensional with double precision
    ExpDist3D_tiled<double, 3>(A, B, m, n, scale_A, scale_B, cross_term);

}

extern "C"
__global__ void
ExpDist_column(const double *A, const double *B,
                 int m, int n, const double *scale_A, const double *scale_B, double *cross_term) {

    //2-dimensional with double precision
    //ExpDist_tiled_column<double, 2>(A, B, m, n, scale_A, scale_B, cross_term);

}

extern "C"
__global__ void
ExpDist_column3D(const double *A, const double *B,
                 int m, int n, const double *scale_A, const double *scale_B, double *cross_term) {


    return;
}

/*
 * Reduce the per thread block cross terms computed in the GaussTransform kernel to single value
 *
 * This kernel is designed to run as single-thread block, because the number of terms to reduce is
 * of size n or m, which is expected to be around 2000 or so. The number of items to reduce
 * is passed as the last argument 'nblocks', which corresponds to the number of thread blocks used
 * by the first kernel.
 */
extern "C"
__global__ void reduce_cross_term(double *output, double *d_cross_term, int m, int n, int nblocks) {

    int tx = threadIdx.x;
    // Specialize BlockReduce for a 1D block of block_size threads on type double
    typedef cub::BlockReduce<double, block_size> BlockReduce;
    // Allocate shared memory for BlockReduce
    __shared__ typename BlockReduce::TempStorage temp_storage;

    double cross_term = 0.0;
    for (int i=tx; i<nblocks; i+=block_size) {
        cross_term += d_cross_term[i];
    }

    //reduce to single value within thread block
    cross_term = BlockReduce(temp_storage).Sum(cross_term);

    //thread 0 writes output
    if (tx == 0) {
        output[0] = cross_term;
    }

}


extern "C"
__global__ void rotate_scales_double(double *rotated_scales, const int n, const double *scale_B) {

    int x = blockIdx.x * block_size_x + threadIdx.x;

    if (x < n) {
        rotate_scale(rotated_scales, rotation_matrixd, rotation_matrix_transposedd, x, scale_B);
    }

}

extern "C"
__global__ void rotate_B_double(double *rotated_B, const int n, const double *B, double const * rotation_matrix) {

    int x = blockIdx.x * block_size_x + threadIdx.x;

    if (x < n) {
        rotate_B_point(rotated_B, rotation_matrix, x, B);
    }

}
