# (C) Copyright 2018-2020      
# Faculty of Applied Sciences
# Delft University of Technology
# Ben van Werkhoven, November 2020.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

from nose.tools import nottest

import numpy
import numpy as np

from numpy.linalg import inv

from kernel_tuner import run_kernel

from test_utils import get_kernel_path, generate_inputs, call_reference_function
from test_matrix_functions import generate_wrapper


def test_expdist_ref():

    size = numpy.int32(100)
    ndim = numpy.int32(2)
    cost, A, B, scale_A, scale_B = generate_inputs(size, size, ndim, 1)

    arguments = [cost, A, B, size, size, ndim, scale_A, scale_B]

    kernel_string = get_kernel_path('expdist')+'expdist_c.cu'

    cost = call_reference_function(size, size, ndim, A, B, scale_A, scale_B, cost)

    print("cost")
    print(cost)

    print("A")
    print(A)
    print("B")
    print(B)
    print("scale_A")
    print(scale_A)
    print("scale_B")
    print(scale_B)

    assert 4000.0 < cost and cost < 6000.0



@nottest
def bhatdist_python_reference(cost, A, B, m, n, ndim, scale_A, scale_B):

    cross_term = 0.0;

    for i in range(m):
        pA = np.array([A[i], A[i+m], A[i+2*m]])

        #construct Sigma_i
        Sigma_i = np.zeros((3,3), dtype=np.float64)
        Sigma_i[0,0] = scale_A[i*2+0]
        Sigma_i[1,1] = scale_A[i*2+0]
        Sigma_i[2,2] = scale_A[i*2+1]

        #print(Sigma_i)

        for j in range(n):
            pB = np.array([B[j], B[j+n], B[j+2*n]])

            #determine (x_t,i - M(x_m,j))
            dist = pA - pB

            #read in (R x Sigma_j x R^T)
            Sigma_j = np.zeros((3,3), dtype=np.float64)
            Sigma_j[0,:] = scale_B[j*9:j*9+3]
            Sigma_j[1,:] = scale_B[j*9+3:j*9+6]
            Sigma_j[2,:] = scale_B[j*9+6:j*9+9]

            #obtain inverted matrix from addition of the two covariance matrices
            inputmatrix = Sigma_i+Sigma_j

            Sigma_inv = inv(inputmatrix)
            #print(inputmatrix)
            #print(Sigma_inv)


            #multiply inverted matrix with vector (x_t,i - M(x_m,j)), obtain vector temp
            temp = Sigma_inv.dot(dist)

            #multiply vector (x_t,i - M(x_m,j)) with vector temp, obtain scalar exponent
            exponent = -1.0*dist.dot(temp)

            #print(exponent, np.exp(exponent), temp)

            #calculate exponent function and add to cross_term
            cross_term += np.exp(exponent)

    return cross_term;



@nottest
def rotate_scales_python(rotated_scales, rotation_matrix, scale_B):

    R = rotation_matrix.reshape(3,3)
    RT = R.T

    for i in range(len(scale_B)//2):

        sigma = np.zeros((3,3), dtype=np.float64)
        sigma[0,0] = scale_B[i*2+0]
        sigma[1,1] = scale_B[i*2+0]
        sigma[2,2] = scale_B[i*2+1]
        print(sigma)

        temp = sigma.dot(RT)

        rotated_scales[i] = R.dot(temp)

    return rotated_scales


def test_rotate_scales():

    #define a rotation matrix
    theta_x = np.deg2rad(32)
    Rx = np.array([[np.cos(theta_x), -np.sin(theta_x), 0],
                   [np.sin(theta_x), np.cos(theta_x), 0],
                   [0, 0, 1]])

    theta_y = np.deg2rad(64)
    Ry = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
                   [0, 1, 0],
                   [-np.sin(theta_y), 0, np.cos(theta_y)]])
    R = Ry.dot(Rx)

    rotated_scales = np.zeros((2,3,3), dtype=np.float64)
    scale_B = np.array([0.75, 0.9, 0.02, 0.65], dtype=np.float64)

    n = np.int32(2)
    args = [rotated_scales.flatten(), R, n, scale_B]

    kernel_string = generate_wrapper("rotate_scales", "expdist_ref.h", args)
    cp = ['--std=c++11', '-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"]

    print(kernel_string)

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")


    expected = rotate_scales_python(np.zeros_like(rotated_scales), R, scale_B)

    print(answer[0].reshape(2,3,3))
    print(expected)


    assert np.allclose(answer[0], expected.ravel())


def test_expdist_ref3D():

    size = numpy.int32(100)
    ndim = numpy.int32(3)
    cost, A, B, scale_A, scale_B = generate_inputs(size, size, ndim, 1)

    arguments = [cost, A, B, size, size, ndim, scale_A, scale_B]

    #with open(get_kernel_path('expdist')+'expdist_c.cu', 'r') as f:
    #    kernel_string = f.read()
    kernel_string = get_kernel_path('expdist')+'expdist_c.cu'

    answer = run_kernel("call_expdist3D_double", kernel_string, size, arguments, {},
               lang="C", compiler_options=['-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"], compiler="nvcc")

    print("cost")
    print(answer[0])
    cost = answer[0]

    print("cost computed by Python reference func:")
    python_cost = bhatdist_python_reference(*arguments)
    print(python_cost)

    print("A")
    print(A)
    print("B")
    print(B)
    print("scale_A")
    print(scale_A)
    print("scale_B")
    print(scale_B)

    assert np.isclose(answer[0], python_cost)

    assert 2000 < cost and cost < 4000
