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

#ifndef EXPDIST_FUNCTIONS_CUH
#define EXPDIST_FUNCTIONS_CUH

/* 
 * This file contains the shared functionality between the CPU and GPU version of
 * the Expdist 3D implementation
 */
 
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE 
#endif // __CUDACC__

#include "matrix_functions.cuh"
#include <math.h>

template <typename T>
HOST_DEVICE
void rotate_scale(T *rotated_scales, const T *rotation_matrix, const T *transposed_rotation_matrix, const int i, const T *scale_B) {

    //construct matrix sigma with uncertainties on diagonal
    T sigma[9];
    zero_matrix(sigma);
    sigma[0] = scale_B[i*2+0];   // 0 1 2
    sigma[4] = scale_B[i*2+0];   // 3 4 5
    sigma[8] = scale_B[i*2+1];   // 6 7 8
            
    T rotation_matrix2[9];
    zero_matrix(rotation_matrix2);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            rotation_matrix2[i*3+j] = rotation_matrix[i+j*3]; 
        }
    }

    T transposed_rotation_matrix2[9];
    zero_matrix(transposed_rotation_matrix2);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            transposed_rotation_matrix2[i*3+j] = transposed_rotation_matrix[i+j*3]; 
        }
    }

    //multiply sigma with transposed rotation matrix
    T temp[9];
    multiply_matrix<T, 9, 3>(temp, sigma, reinterpret_cast<const T(&)[9]>(*transposed_rotation_matrix2));

    //multiply with rotation matrix, reuse sigma to store result
    multiply_matrix<T, 9, 3>(sigma, reinterpret_cast<const T(&)[9]>(*rotation_matrix2), temp);

    //store result
    for (int ii=0; ii<9; ii++) {
        rotated_scales[i*9+ii] = sigma[ii];
    }

}

template <typename T>
HOST_DEVICE
void rotate_B_point(T *rotated_B, const T *rotation_matrix, const int i, const T *B) {

    multiply_matrix_vector<T, 3>(reinterpret_cast<T(&)[3]>(*(rotated_B+i*3)),
                                 reinterpret_cast<const T(&)[9]>(*rotation_matrix),
                                 reinterpret_cast<const T(&)[3]>(*(B+i*3)));
}


template <typename T, int dim>
HOST_DEVICE
T compute_expdist_3D(const T (&A)[dim], const T (&B)[dim], const T (&Sigma_i)[9], const T (&Sigma_j)[9]) {

    T cross_term = 0;

    //determine (x_t,i - M(x_m,j)) and store as dist_ij
    T dist_ij[dim];
    for (int d=0; d<dim; d++) {
        dist_ij[d] = A[d] - B[d];
    }

    T temp_matrix[9];
    //determine (Sigma_{t,i} + R*Sigma_{m,j}*R^{T})^-1   (rotations to Sigma_j already applied)
    add_matrix(temp_matrix, Sigma_i, Sigma_j);

    //obtain inverted matrix
    T Inverted[9];
    invert_matrix(Inverted, temp_matrix);
    /*
    printf("input:\n");
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            printf("%f ", Sigma_i[ii*3+jj]);
        }
        printf("\n");
    }
    printf("inverted:\n");
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            printf("%f ", Inverted[ii*3+jj]);
        }
        printf("\n");
    }
    */

    //multiply with dist_ij
    T temp[dim];
    multiply_matrix_vector(temp, Inverted, dist_ij);

    //dot product to obtain scalar value
    T exponent = 0;
    dot_product(exponent, dist_ij, temp);
    //printf("%f, %f, [%f, %f, %f]\n", -1.0*exponent, exp(-1.0*exponent), temp[0], temp[1], temp[2]);

    //cross_term += exp(-dist_ij / (scale_A[i] + scale_B[j]) );
    T norm = 0; 
    norm = sqrt(determinant(temp_matrix));
    cross_term += exp(-1.0*exponent)/norm;

    return cross_term;

}

#endif // !EXPDIST_FUNCTIONS
