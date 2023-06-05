# Software for 3D Fast Particle Fusion in Localization Microscopy.



## Requirements

- [MATLAB](https://nl.mathworks.com/products/matlab.html) (version>2018a)
    - distrib_computing_toolbox  
    - statistics_toolbox  
- [CMake](https://cmake.org/) (version>3.14.3)
- [CUDA toolkit](https://developer.nvidia.com/cuda-downloads) (version>8.0) 
- [CUB library](https://nvlabs.github.io/cub/) (version>1.8.0)
- [The DIPImage toolbox](http://www.diplib.org) (version>2.9)
- [GNU C compiler](https://gcc.gnu.org/) (version>5.5.0, for Linux)
- [Viusal Studio](https://visualstudio.microsoft.com/downloads/) (Visual Studio 15 2017, for Windows)

## Installation and usage on Linux

### Get the sources

The Git repository uses submodules. Include them in a _git clone_ action using the _--recursive_ option.
```bash

git clone https://github.com/wexw/Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy.git
cd Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy/3d

````
### Compile the code
In the following

- BUILD_DIRECTORY is the directory where the project will be built
- SOURCE_DIRECTORY is the root directory of the sources
- CUB_DIRECTORY is the root directory of the downloaded [CUB library](https://nvlabs.github.io/cub/) sources
- MATLAB_DIRECTORY is the root of MATLAB installation directory (e.g. /usr/local/MATLAB/R2019a)

Use the following commands to build the necessary libraries for this software:

```bash

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=gcc-5 -DCUB_ROOT_DIR=CUB_DIRECTORY SOURCE_DIRECTORY
make
````
### Use the code
Next, we need to locate the built libraries for MATLAB:
```bash

cd ..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:MATLAB_DIRECTORY/runtime/glnxa64:MATLAB_DIRECTORY/bin/glnxa64:MATLAB_DIRECTORY/sys/os/glnxa64:MATLAB_DIRECTORY/sys/opengl/lib/glnxa64:BUILD_DIRECTORY/mex
``` 
Then, run MATLAB and the demo script
```bash
matlab
>> Demo_JCC_3D.m
```

before you run it make sure that you have CUDA installed:
```
pip install -r requirements.txt
```

The tests can be run using ``nose``, for example by typing the following in 
the top-level or test directory:
```
nosetests -v
```
