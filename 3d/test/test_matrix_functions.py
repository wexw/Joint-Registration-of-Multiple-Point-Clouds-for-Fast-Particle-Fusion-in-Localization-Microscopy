# (C) Copyright 2018-2020      
# Faculty of Applied Sciences
# Delft University of Technology
# Ben van Werkhoven, November 2020.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

import numpy as np

from kernel_tuner import run_kernel

from test_utils import get_kernel_path

# trying to think of a general function that could wrap C++ functions
# allowing them to be called by Kernel Tuner


def generate_wrapper(function_name, filename, args, convert_to_array=None, template_parameters=None):
    """ Generate a wrapper to call C++ functions from Python

        This function allows Kernel Tuner to call templated simple-typed C++ functions.

        There is support to convert function arguments from plain pointers
        to array references, if this is needed there should be a True value in
        convert_to_array in the location corresponding to the location in the args array.

    """


    type_map = {"int8": "char",
                "int16": "short",
                "int32": "int",
                "float32": "float",
                "float64": "double"}

    def type_str(arg):
        typestring = type_map[str(arg.dtype)]
        if isinstance(arg, np.ndarray):
            typestring += " *"
        return typestring + " "

    signature = ",".join([type_str(arg) + "arg" + str(i) for i, arg in enumerate(args)])

    if not convert_to_array:
        call_args = ",".join(["arg" + str(i) for i in range(len(args))])
    else:
        call_args = []
        for i, arg in enumerate(args):
            if convert_to_array[i]:
                if np.prod(arg.shape) > 1:
                    #convert pointer to a reference to an array
                    arg_shape = "".join("[%d]" % i for i in arg.shape)
                    arg_str = "*reinterpret_cast<" + type_map[str(arg.dtype)] + "(*)" + arg_shape + ">(arg" + str(i) + ")"
                else:
                    #a reference is accepted rather than a pointer, just dereference
                    arg_str = "*arg" + str(i)
                call_args.append(arg_str)
                #call_args = ",".join(["*reinterpret_cast<double(*)[9]>(arg" + str(i) + ")" for i in range(len(args))])
            else:
                call_args.append("arg" + str(i))
        call_args = ",".join(call_args)

    if template_parameters != None:
        template_parameters = "<" + template_parameters + ">"
    else:
        template_parameters = ""

    return """
    #include "%s"

    extern "C"
    float call_function(%s) {

        %s%s(%s);

        return 0.0f;
    }""" % (filename, signature, function_name, template_parameters, call_args)


filename = "matrix_functions.cuh"
cp = ['--std=c++11', '-I'+get_kernel_path('expdist'), "-Wno-deprecated-gpu-targets"]

def assert_invert_matrix(matrix):

    function_name = "invert_matrix"

    matrix = np.random.randn(9).astype(np.float64)
    inverted_matrix = np.zeros(9, dtype=np.float64)

    args = [inverted_matrix, matrix]

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=[True for _ in args])

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = np.linalg.inv(matrix.reshape(3,3))

    assert np.allclose(answer[0].reshape(3,3), expected, atol=1e-6)

def test_multiply_matrix():

    function_name = "multiply_matrix"

    a = np.random.randn(9).astype(np.float64)
    b = np.random.randn(9).astype(np.float64)
    c = np.zeros_like(a)

    #args = [c, a, b, np.int32(3)]
    args = [c, a, b]
    convert = [True for _ in args]
    #convert[-1] = False

    template_parameters = "double, 9, 3"

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=convert,
                        template_parameters=template_parameters)

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = a.reshape(3,3).dot(b.reshape(3,3))

    print("answer")
    print(answer[0].reshape(3,3))

    print("expected")
    print(expected)

    assert np.allclose(answer[0].reshape(3,3), expected, atol=1e-6)


def test_multiply_matrix_vector():

    function_name = "multiply_matrix_vector"

    a = np.random.randn(9).astype(np.float64)
    b = np.random.randn(3).astype(np.float64)
    c = np.zeros_like(b)

    args = [c, a, b]
    convert = [True for _ in args]

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=convert)

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = a.reshape(3,3).dot(b)

    print("answer")
    print(answer[0])

    print("expected")
    print(expected)

    assert np.allclose(answer[0], expected, atol=1e-6)



def test_dot_product():

    function_name = "dot_product"

    a = np.random.randn(3).astype(np.float64)
    b = np.random.randn(3).astype(np.float64)
    c = np.zeros((1), dtype=np.float64)

    args = [c, a, b]
    convert = [True, True, True]

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=convert)

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = a.dot(b)

    print("answer")
    print(answer[0])

    print("expected")
    print(expected)

    assert np.allclose(answer[0], expected, atol=1e-6)



def test_add_matrix():

    function_name = "add_matrix"

    a = np.random.randn(9).astype(np.float64)
    b = np.random.randn(9).astype(np.float64)
    c = np.zeros_like(a)

    args = [c, a, b]
    convert = [True for _ in args]

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=convert)

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = a + b

    print("answer")
    print(answer[0])

    print("expected")
    print(expected)

    assert np.allclose(answer[0], expected, atol=1e-6)




def test_transpose_matrix():

    function_name = "transpose_matrix"

    a = np.random.randn(9).astype(np.float64)
    b = np.zeros_like(a)

    args = [b, a]
    convert = [True for _ in args]

    kernel_string = generate_wrapper(function_name, filename, args, convert_to_array=convert, template_parameters="double, 9, 3")

    answer = run_kernel("call_function", kernel_string, 1, args, {},
               lang="C", compiler_options=cp, compiler="nvcc")

    expected = a.reshape(3,3).T.flatten()

    print("answer")
    print(answer[0])

    print("expected")
    print(expected)

    assert np.allclose(answer[0], expected, atol=1e-6)

