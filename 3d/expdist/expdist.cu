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
        
/*
 * Host part for calling GPUExpDist from the CPU 
 */

#include <stdint.h>
#include <cuda_runtime.h>
#include <math.h>

#include "expdist.h"

//tuned for Nvidia K40
#ifndef block_size_x //if not using kernel tuner
#define block_size_x 32
#define block_size_y 4
#define tile_size_x 2
#define tile_size_y 4
#define use_shared_mem 1

#endif
#define reduce_block_size 256


#define INTERMEDIATE_BUFFER_MULTIPLE 20

//#include "kernels.cuh"
#include "kernels.cu"

GPUExpDist::GPUExpDist(int n, int argdim) {
    //allocate GPU memory for size max_n
    max_n = n;
    dim = argdim;
    int elems = max_n * dim;

    //pseudo load balancing across available GPUs
    int count;
    cudaGetDeviceCount(&count);
    int rand_int = (int)rand();
    int id = rand_int % count;
    cudaSetDevice(id);
    //printf("DEBUG Expdist: count returned %d, rand returned %d, using GPU %d\n", count, rand_int, id);

    scale_A_dim = 1;
    scale_B_dim = 1;
    if (dim == 3) {
        scale_A_dim = 2;
        scale_B_dim = 9;
    }

    cudaError_t err;

    err = cudaMalloc((void **)&d_A, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_B, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    /*
    err = cudaMalloc((void **)&d_B_temp, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    */
    err = cudaMalloc((void **)&d_scale_A, scale_A_dim*max_n*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_scale_B, scale_B_dim*max_n*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_cross_term, INTERMEDIATE_BUFFER_MULTIPLE*max_n*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    if (dim == 3) {
        err = cudaMalloc((void **)&d_scale_B_temp, scale_A_dim*max_n*sizeof(double));
        if (err != cudaSuccess) {
            fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }

    err = cudaGetSymbolAddress((void**)&ptrto_rotation_matrixd, rotation_matrixd);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaGetSymbolAddress: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaGetSymbolAddress((void**)&ptrto_rotation_matrix_transposedd, rotation_matrix_transposedd);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaGetSymbolAddress: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaEventCreate(&event);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaEventCreate: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaStreamCreate(&stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaStreamCreate: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaStreamCreate(&stream_b);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaStreamCreate: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    cudaDeviceSynchronize();
} 

GPUExpDist::~GPUExpDist() {
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_scale_A);
    cudaFree(d_scale_B);
    cudaFree(d_cross_term);
    cudaEventDestroy(event);
    cudaStreamDestroy(stream);
    cudaStreamDestroy(stream_b);
} 

double GPUExpDist::compute(const double *A, const double *B, int m, int n, const double *scale_A, const double *scale_B) {
  return GPUExpDist::compute(A, B, m, n, scale_A, scale_B, (const double *)NULL, 0);
}

double GPUExpDist::compute(const double *A, const double *B, int m, int n, const double *scale_A, const double *scale_B, const double *rotation_matrix, int use_prerotated_scale_B) {

    double cost;
    cudaError_t err;

    //move data to the GPU
    err = cudaMemcpyAsync(d_A, A, m*dim*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_A: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_B, B, n*dim*sizeof(double), cudaMemcpyHostToDevice, stream); //should be d_B_temp when rotating B on GPU
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_B: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_scale_A, scale_A, scale_A_dim*m*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_scale_A: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    //prepare scale_B if needed
    if (dim == 3 && use_prerotated_scale_B == 0) {
        err = cudaMemcpyAsync(d_scale_B_temp, scale_B, scale_A_dim*n*sizeof(double), cudaMemcpyHostToDevice, stream_b);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error in cudaMemcpyAsync d_scale_B: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        //compute tranpose of rotation matrix and send both to constant memory
        double transposed_rotation_matrix[9];
        transpose_matrix<double, 9, 3>(transposed_rotation_matrix, reinterpret_cast<const double(&)[9]>(*rotation_matrix));
        
        err = cudaMemcpyToSymbolAsync(rotation_matrixd, rotation_matrix, 9*sizeof(double), 0, cudaMemcpyHostToDevice, stream_b);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error in cudaMemcpyToSymbolAsync rotation_matrix: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        err = cudaMemcpyToSymbolAsync(rotation_matrix_transposedd, transposed_rotation_matrix, 9*sizeof(double), 0, cudaMemcpyHostToDevice, stream_b);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error in cudaMemcpyToSymbolAsync rotation_matrix_transposed: %s\n", cudaGetErrorString(err));
            exit(1);
        }

        //call rotate scales kernel
        dim3 threads(block_size_x, 1, 1);
        dim3 grid((int) ceil(n / (float)(block_size_x)), 1, 1);
        rotate_scales_double<<<grid, threads, 0, stream_b>>>(d_scale_B, n, d_scale_B_temp);
        cudaEventRecord(event, stream_b);


    }
    else {
        err = cudaMemcpyAsync(d_scale_B, scale_B, scale_B_dim*n*sizeof(double), cudaMemcpyHostToDevice, stream);
        if (err != cudaSuccess) {
            fprintf(stderr, "Error in cudaMemcpyAsync d_scale_B: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }

    if (dim == 3) {
        cudaStreamWaitEvent(stream, event, 0);
    }
    /*
    if (dim == 3) {
        //also rotate B, after the recorded event because
        //transfer of B (in stream) and rotation matrix
        //(in stream_b) needs to be completed
        //
        //yes, at the moment only for the 3D case the B coordinates are rotated on the GPU, the 2D version assumes
        //that the coordinates have been rotated before data is transferred to the GPU
        dim3 threads(block_size_x, 1, 1);
        dim3 grid((int) ceil(n / (float)(block_size_x)), 1, 1);
        rotate_B_double<<<grid, threads, 0, stream>>>(d_B, n, d_B_temp, ptrto_rotation_matrixd);

    }
    */

    //compute number of thread blocks that would be used by the ExpDist kernel for this m and n
    int nblocks = ((int) ceil(m / (float)(block_size_x*tile_size_x)) * (int) ceil(n / (float)(block_size_y*tile_size_y)));

    //setup kernel execution parameters
    dim3 threads(block_size_x, block_size_y, 1);
    dim3 grid(1, 1, 1); //to be overwritten

    //check if the number of thread blocks does not exceed the allocated space
    //if it does, run the ExpDist_column kernel that uses fewer thread blocks
    if (nblocks < INTERMEDIATE_BUFFER_MULTIPLE*max_n) {
        //setup kernel execution parameters
        grid.x = (int) ceilf(m / (float)(block_size_x * tile_size_x));
        grid.y = (int) ceilf(n / (float)(block_size_y * tile_size_y));
    
        //call the first kernel
        if (dim == 2) {
            ExpDist<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 
        } else {
            ExpDist3D<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 
        }

    } else {

        fprintf(stderr, "Error in Expdist GPU: number of blocks exceeds buffer for intermediate results\n");
        exit(1);

        /*
        //setup kernel execution parameters
        grid.x = (int) ceilf(m / (float)(block_size_x * tile_size_x));
    
        //call the first kernel
        if (dim == 2) {
            ExpDist_column<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 
        } else {
            ExpDist_column3D<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 
        }
        */

    }

    //call the second kernel
    dim3 threads2(reduce_block_size, 1, 1);
    dim3 grid2(1, 1, 1);
    reduce_cross_term<<<grid2, threads2, 0, stream>>>(d_cross_term, d_cross_term, m, n, grid.x*grid.y);

    err = cudaMemcpyAsync(&cost, d_cross_term, 1*sizeof(double), cudaMemcpyDeviceToHost, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyDeviceToHost cross_term: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    return cost;
}



extern "C"
float test_GPUExpDistHost(double *cost, const double* A, const double* B,
            int m, int n, int dim, const double *scale_A, const double *scale_B, int max_n, const double *rotation_matrix, int use_prerotated_scale_B) {

    GPUExpDist gpu_expdist(1000000, dim);

    if (dim == 2) {
        *cost = gpu_expdist.compute(A, B, m, n, scale_A, scale_B);
    } else {
        *cost = gpu_expdist.compute(A, B, m, n, scale_A, scale_B, rotation_matrix, use_prerotated_scale_B);
    }

    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in test_GPUExpDistHost: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    return 0.0;
}


