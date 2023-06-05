# (C) Copyright 2018-2020      
# Faculty of Applied Sciences
# Delft University of Technology
# Ben van Werkhoven, November 2020.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

import numpy

from kernel_tuner import run_kernel

from test_utils import get_kernel_path

def test_gausstransform_ref():

    size = numpy.int32(2000)
    ndim = numpy.int32(2)
    A = numpy.random.randn(size*ndim).astype(numpy.float64)
    B = numpy.random.randn(size*ndim).astype(numpy.float64)
    scale = numpy.float64(10.0)
    grad = numpy.zeros(size*ndim).astype(numpy.float64)
    cost = numpy.zeros((1)).astype(numpy.float64)

    arguments = [cost, A, B, size, size, ndim, scale, grad]

    with open(get_kernel_path('gausstransform')+'gausstransform_c.cpp', 'r') as f:
        kernel_string = f.read()

    answer = run_kernel("call_GaussTransform", kernel_string, size, arguments, {},
               lang="C", compiler_options=['-I'+get_kernel_path('gausstransform')])

    cost = answer[0]
    print(cost)

    assert 1.0 > cost and cost > 0.0

    gradient = answer[7]
    print(gradient)

    #TODO: somehow assert that the gradient results are sensible
