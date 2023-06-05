# (C) Copyright 2018-2020      
# Faculty of Applied Sciences
# Delft University of Technology
# Ben van Werkhoven, November 2020.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

import os

import numpy
from nose.tools import nottest
from kernel_tuner import run_kernel

@nottest
def get_kernel_path(dir=None):
    """ get path to the kernels as a string """
    path = "/".join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
    return path+'/' + (dir + '/' if dir else '')


@nottest
def generate_inputs(m, n, ndim=2, nblocks=1):
    A = numpy.random.randn(m*ndim).astype(numpy.float64)
    B = numpy.random.randn(n*ndim).astype(numpy.float64)
    if ndim==3:
        scale_A = numpy.absolute(2+0.01*numpy.random.randn(m*2).astype(numpy.float64))
        scale_B = numpy.absolute(2+0.01*numpy.random.randn(n*9).astype(numpy.float64))
    else:
        scale_A = numpy.absolute(2+0.01*numpy.random.randn(m).astype(numpy.float64))
        scale_B = numpy.absolute(2+0.01*numpy.random.randn(n).astype(numpy.float64))
    cost = numpy.zeros((nblocks)).astype(numpy.float64)
    return cost, A, B, scale_A, scale_B



@nottest
def call_reference_function(m, n, ndim, A, B, scale_A, scale_B, cost):
    arguments = [cost, A, B, m, n, ndim, scale_A, scale_B]
    with open(get_kernel_path('expdist')+'expdist_c.cu', 'r') as f:
        kernel_string = f.read()
    f = "call_expdist"
    if ndim == 3:
        f = "call_expdist3D_double"
    answer = run_kernel(f, kernel_string, 1, arguments, {},
               lang="C", compiler_options=['-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"], compiler='nvcc')

    ref_cost = answer[0][0]
    return ref_cost



