#!/usr/bin/env python
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
from kernel_tuner import run_kernel

from test_expdist_ref import bhatdist_python_reference
from test_utils import get_kernel_path, generate_inputs, call_reference_function
from test_matrix_functions import generate_wrapper

compiler_options = ['-I'+get_kernel_path('expdist')]



@nottest
def test_against_reference(cost, A, B, scale_A, scale_B, m, n, ndim, nblocks, params):

    #call the GPU function
    with open(get_kernel_path('expdist')+'kernels.cu', 'r') as f:
        kernel_string = f.read()

    #mimic Hamid's testcase
    #A = numpy.ones_like(A)
    #B = 2.0*numpy.ones_like(B)
    #scale_A = 0.1*numpy.ones_like(scale_A)
    #scale_B = 0.1*numpy.ones_like(scale_B)

    print(A, B, scale_A, scale_B)

    nblocks = numpy.int32( numpy.ceil(m / float(params["block_size_x"]*params["tile_size_x"])) *
                           numpy.ceil(n / float(params["block_size_y"]*params["tile_size_y"])) )

    arguments = [A, B, m, n, scale_A, scale_B, cost]

    grid_div_x = ["block_size_x", "tile_size_x"]
    grid_div_y = ["block_size_y", "tile_size_y"]

    if ndim == 2:
        answer = run_kernel("ExpDist", kernel_string, (m, n), arguments, params,
                   compiler_options=compiler_options, grid_div_x=grid_div_x, grid_div_y=grid_div_y)
    else:
        answer = run_kernel("ExpDist3D", kernel_string, (m, n), arguments, params,
                   compiler_options=compiler_options, grid_div_x=grid_div_x, grid_div_y=grid_div_y)

    #collect the results from the first kernel
    cross_term = answer[6]
    print("intermediate cross_term")
    print(cross_term)
    print(numpy.sum(cross_term))

    #call the second kernel to reduce the per thread block cross terms to a single value
    out = numpy.zeros(1).astype(numpy.float64)

    arguments = [out, cross_term, m, n, nblocks]
    answer = run_kernel("reduce_cross_term", kernel_string, 1, arguments, {"block_size_x": 128},
               compiler_options=compiler_options, grid_div_x=[])

    #call the reference function
    ref_cost = call_reference_function(m, n, ndim, A, B, scale_A, scale_B, cost)

    #final cross term
    cost = answer[0][0]

    if ndim == 3:
        refpy_cost = bhatdist_python_reference(cost, A, B, m, n, ndim, scale_A, scale_B)
        print("reference python")
        print(refpy_cost)

    print("reference")
    print(ref_cost)
    print("answer")
    print(cost)

    print("reference")
    print("%30.20e" % ref_cost)
    print("answer")
    print("%30.20e" % cost)
    if ndim == 3:
        print("reference python")
        print("%30.20e" % refpy_cost)

    assert numpy.isclose(ref_cost, cost, atol=1e-5)

@nottest
def test_expdist_kernel(dim=2):

    ndim = numpy.int32(dim)
    m = numpy.int32(103)
    n = numpy.int32(59)

    #block_size_x=32, block_size_y=4, tile_size_x=2, tile_size_y=1, use_shared_mem=1

    params = dict()
    params["block_size_x"] = 32
    params["block_size_y"] = 4
    params["tile_size_x"] = 2
    params["tile_size_y"] = 1
    params["use_shared_mem"] = 1

    nblocks = numpy.int32( numpy.ceil(m / float(params["block_size_x"]*params["tile_size_x"])) *
                           numpy.ceil(n / float(params["block_size_y"]*params["tile_size_y"])) )

    cost, A, B, scale_A, scale_B = generate_inputs(m, n, ndim, nblocks)

    test_against_reference(cost, A, B, scale_A, scale_B, m, n, ndim, nblocks, params)


def test_expdist_kernel2D():
    test_expdist_kernel(dim=2)

def test_expdist_kernel3D():
    test_expdist_kernel(dim=3)




@nottest
def test_expdist_kernel_column(dim=2):

    #setup test input
    allocation_size = int(3000)
    ndim = numpy.int32(dim)
    size = numpy.int32(2000)

    params = dict()
    params["block_size_x"] = 32
    params["block_size_y"] = 4
    params["tile_size_x"] = 2
    params["tile_size_y"] = 4
    params["use_shared_mem"] = 1

    nblocks = numpy.int32( numpy.ceil(size / float(params["block_size_x"]*params["tile_size_x"])) )

    cost, A, B, scale_A, scale_B = generate_inputs(allocation_size, size, ndim, nblocks)

    #call the reference function
    ref_cost = call_reference_function(size, size, ndim, A, B, scale_A, scale_B, cost)

    #call the GPU function
    with open(get_kernel_path('expdist')+'kernels.cu', 'r') as f:
        kernel_string = f.read()

    arguments = [A, B, size, size, scale_A, scale_B, cost]

    grid_div_x = ["block_size_x", "tile_size_x"]

    if ndim == 2:
        answer = run_kernel("ExpDist_column", kernel_string, size, arguments, params,
                   compiler_options=compiler_options, grid_div_x=grid_div_x)
    else:
        answer = run_kernel("ExpDist_column3D", kernel_string, size, arguments, params,
                   compiler_options=compiler_options, grid_div_x=grid_div_x)

    #collect the results from the first kernel
    cross_term = answer[6]
    print("intermediate cross_term")
    print(cross_term)

    #call the second kernel to reduce the per thread block cross terms to a single value
    out = numpy.zeros(1).astype(numpy.float64)

    arguments = [out, cross_term, size, size, nblocks]
    answer = run_kernel("reduce_cross_term", kernel_string, 1, arguments, {"block_size_x": 128},
               compiler_options=compiler_options, grid_div_x=[])

    #final cross term
    cost = answer[0][0]

    print("reference")
    print(ref_cost)
    print("answer")
    print(cost)

    print("reference")
    print("%30.20e" % ref_cost)
    print("answer")
    print("%30.20e" % cost)

    assert numpy.isclose(ref_cost, cost, atol=1e-5)

@nottest
def test_expdist_kernel_column2D():
    test_expdist_kernel_column(dim=2)

@nottest
def test_expdist_kernel_column3D():
    test_expdist_kernel_column(dim=3)

@nottest
def test_hostfunction(dim=2):

    #setup test input
    ndim = numpy.int32(dim)

    m = numpy.int32(2003)
    n = numpy.int32(1009)
    nblocks = numpy.int32(numpy.ceil(m / (32*2)) * numpy.ceil(n / (4*4)))

    cost, A, B, scale_A, scale_B = generate_inputs(m, n, ndim, nblocks)
    #host function will do the rotation, so we need to supply the scales indirectly
    scale_B = numpy.absolute(0.01*numpy.random.randn(n*2).astype(numpy.float64))
    rotation_matrix = numpy.eye(3).astype(numpy.float64).flatten()

    #mimic Hamid's testcase
    A = numpy.ones_like(A)
    B = 2.0*numpy.ones_like(B)
    scale_A = 0.1*numpy.ones_like(scale_A)
    scale_B = 0.1*numpy.ones_like(scale_B)

    #call the reference function
    #with open(get_kernel_path('expdist')+'expdist_c.cu', 'r') as f:
    #    kernel_string = f.read()
    kernel_string = get_kernel_path('expdist')+'expdist_c.cu'

    f = "call_expdist"
    scale_B_rot = scale_B
    if ndim == 3:
        #first call the rotate scales kernel
        rotated_scales = numpy.zeros(n*9).astype(numpy.float64)
        args = [rotated_scales, rotation_matrix, n, scale_B]
        answer = run_kernel("call_rotate_scales_double", kernel_string, 1, args, {},
                lang="C", compiler_options=['-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"], compiler='nvcc')
        scale_B_rot = answer[0]
        f = "call_expdist3D_double"


    arguments = [cost, A, B, m, n, ndim, scale_A, scale_B_rot]
    answer = run_kernel(f, kernel_string, 1, arguments, {},
               lang="C", compiler_options=['-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"], compiler='nvcc')

    ref_cost = answer[0][0]


    #call the host function
    arguments = [cost, A, B, m, n, ndim, scale_A, scale_B, numpy.int32(100000), rotation_matrix, np.int32(0)]
    #with open(get_kernel_path('expdist')+'expdist.cu', 'r') as f:
    #    kernel_string = f.read()
    kernel_string = get_kernel_path('expdist')+'expdist.cu'
    answer = run_kernel("test_GPUExpDistHost", kernel_string, 1, arguments, {},
               lang="C", compiler_options=compiler_options+['-arch=sm_30'])
    cost = answer[0][0]

    print("reference")
    print(ref_cost)

    print("answer")
    print(cost)

    assert numpy.isclose(ref_cost, cost, atol=1e-5)

def test_hostfunction2D():
    test_hostfunction(dim=2)

def test_hostfunction3D():
    test_hostfunction(dim=3)


@nottest
def test_hostfunction_largeN():

    #setup test input
    allocation_size = numpy.int32(1e6)
    size = numpy.int32(40000)
    ndim = numpy.int32(2)

    params = dict()
    params["block_size_x"] = 32
    params["block_size_y"] = 4
    params["tile_size_x"] = 2
    params["tile_size_y"] = 1
    params["use_shared_mem"] = 1

    #compute nblocks for when using the expdist kernel
    nblocks = numpy.int32( numpy.ceil(size / float(params["block_size_x"]*params["tile_size_x"])) *
                           numpy.ceil(size / float(params["block_size_y"]*params["tile_size_y"])) )

    #ensure that this test actually causes the host code to call the column kernel
    assert nblocks > allocation_size

    #compute the nblocks actually used by the column kernel
    nblocks = numpy.int32(numpy.ceil(size / float(params["block_size_x"] * params["tile_size_x"])))

    #generate input data
    cost, A, B, scale_A, scale_B = generate_inputs(allocation_size, allocation_size, ndim, nblocks)

    #call the ExpDist_column kernel directly for reference
    arguments = [A, B, size, size, scale_A, scale_B, cost]
    grid_div_x = ["block_size_x", "tile_size_x"]
    with open(get_kernel_path('expdist')+'kernels.cu', 'r') as f:
        kernel_string = f.read()
    answer = run_kernel("ExpDist_column", kernel_string, size, arguments, params,
               compiler_options=compiler_options, grid_div_x=grid_div_x)
    ref_cost = numpy.sum(answer[6])

    #call the host function
    rot_matrix = numpy.eye(3).astype(numpy.float64)
    arguments = [cost, A, B, size, size, ndim, scale_A, scale_B, allocation_size, rot_matrix, np.int32(0)]
    with open(get_kernel_path('expdist')+'expdist.cu', 'r') as f:
        kernel_string = f.read()
    answer = run_kernel("test_GPUExpDistHost", kernel_string, size, arguments, {},
               lang="C", compiler_options=compiler_options+['-arch=sm_30'])
    cost = answer[0][0]

    print("reference")
    print(ref_cost)

    print("answer")
    print(cost)

    assert numpy.isclose(ref_cost, cost, atol=1e-5)






def test_rotate_scales():

    n = np.int32(200)

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

    rotated_scales = np.zeros((n,3,3), dtype=np.float64)

    scale_B = numpy.absolute(0.001*numpy.random.randn(n*2).astype(numpy.float64))

    args = [rotated_scales.flatten(), R, n, scale_B]

    kernel_string = generate_wrapper("rotate_scales", "expdist_ref.h", args)
    cp = ['--std=c++11', '-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"]

    print(kernel_string)

    expected = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    # call the GPU function
    with open(get_kernel_path('expdist')+'expdist.cu', 'r') as f:
        kernel_string = f.read()
    cmem_args = {"rotation_matrixd": R.flatten(), "rotation_matrix_transposedd": R.T.flatten()}

    params = dict(block_size_x = 128)
    arguments = [rotated_scales.flatten(), n, scale_B]
    answer = run_kernel("rotate_scales_double", kernel_string, n, arguments,
                        params, lang="CUDA", compiler_options=cp, cmem_args=cmem_args)

    print("reference")
    print(expected[0])
    print("answer")
    print(answer[0])

    assert np.allclose(answer[0], expected[0])




if __name__ == "__main__":
    test_expdist_kernel(dim=3)
