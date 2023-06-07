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

#include "expdist_functions.cuh"

double expdist(const double *A, const double *B, int m, int n, int dim, const double *scale_A, const double *scale_B);


/*
* The following function is a full 3D implementation of the Bhattacharya distance
* scale_A is an array with 2 values per localization
* scale_B contains the pre-rotated matrix of uncertainties for B
*/
template <typename T>
T expdist3D(const T *A, const T *B, const int m, const int n, const T *scale_A, const T *scale_B) {
    int i, j;
    T cross_term = 0.0;
    const int dim = 3;

    T pA[dim];
    T pB[dim];

    for (i = 0; i<m; i++) {

        //prefetch point Ai
        for (int d = 0; d<dim; d++) {
            int id = i + d * m;
            pA[d] = A[id];
        }

        //assume sigma in x and y are equal and scale_A only stores 2 values per localization
        T Sigma_i[9];
        zero_matrix(Sigma_i);
        Sigma_i[0] = scale_A[i * 2 + 0];   // 0 1 2
        Sigma_i[4] = scale_A[i * 2 + 0];   // 3 4 5
        Sigma_i[8] = scale_A[i * 2 + 1];   // 6 7 8

        for (j = 0; j<n; j++) {

            //prefetch point Bj
            for (int d = 0; d<dim; d++) {
                int jd = j + d * n;
                pB[d] = B[jd];
            }

            //assume sigma_j has been rotated properly beforehand so that it can be used directly
            T Sigma_j[9];
            load_matrix(Sigma_j, scale_B, j);

            cross_term += compute_expdist_3D<T, 3>(pA, pB, Sigma_i, Sigma_j);
        }
    }

    return cross_term;
}


/*
* This function rotates the uncertainties for the 3D bhattacharya distance
*
* It is assumed that the scale_B array contains 2 uncertainty values per localization in
* the particle. One value for the uncertainty in X and Y and one for the uncertainty in Z (depth).
*
* The output array rotated_scales contains a 3x3 matrix for each localization.
*/
template <typename T>
void rotate_scales(T *rotated_scales, const T *rotation_matrix, const int n, const T *scale_B) {

    T transposed_rotation_matrix[9];
    transpose_matrix<T, 9, 3>(transposed_rotation_matrix, reinterpret_cast<const T(&)[9]>(*rotation_matrix));

    for (int i = 0; i<n; i++) {

        rotate_scale(rotated_scales, rotation_matrix, transposed_rotation_matrix, i, scale_B);
    }

}

/*
* This function rotates the coordinates of the localizations in the B particle for the 3D bhattacharya distance
*
* The output array rotated_B contains the x,y,z coordinates of each localization.
*/
template <typename T>
void rotate_B(T *rotated_B, const T *rotation_matrix, const int n, const T *B) {

    for (int i = 0; i<n; i++) {
        rotate_B_point(rotated_B, rotation_matrix, i, B);
    }

}


