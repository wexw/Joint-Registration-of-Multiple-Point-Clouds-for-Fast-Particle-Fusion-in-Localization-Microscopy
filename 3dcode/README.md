# Software for Template-Free 3D Particle Fusion in Localization Microscopy.

This software implements a template-free particle fusion algorithm based on 
an all-to-all registration, which provides robustness against individual 
mis-registrations and underlabeling. The method does not assume any prior
knowledge about the structure to be reconstructed (template-free) and directly
works on localization data not pixelated images.

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

git clone --single-branch --branch develop https://github.com/imphys/smlm_datafusion3d.git --recursive
cd smlm_datafusion3d/
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
>> demo1.m
```
## Installation and usage on Windows

### Get the sources

The Git repository uses submodules. Include them in a _git clone_ action using the _--recursive_ option.
```bash

git clone --single-branch --branch develop git@github.com:berndrieger/alltoall3D.git --recursive
````
### Compile the code
In the following

- BUILD_DIRECTORY is the directory where the project will be built
- SOURCE_DIRECTORY is the root directory of the sources
- CUB_DIRECTORY is the root directory of the downloaded [CUB library](https://nvlabs.github.io/cub/) sources

Use the following commands from Command Prompt (type `cmd` in Run) to build the necessary libraries for this software:

```bash

mkdir build
cd build
cmake -G "Visual Studio 15 2017 Win64" -DCUB_ROOT_DIR=CUB_DIRECTORY SOURCE_DIRECTORY
make
````

Then open MS Visual Studio, by double clicking on ALL_BUILD.vcxproj and build all the targets (press `F7`). You might need to manually build the targets expdist and gausstransform before building the other targets.

### Use the code

Open MATLAB and comment/uncomment the lines which add mex files to MATLAB path. Then, run the demo script
```bash
>> demo1.m
````

## Troubleshooting

- Matlab mex headers not found  
The Makefile tries to find the directories in which MATLAB was installed on 
your system. If this fails, you can manually insert the path to your MATLAB 
installation (ending with `/extern/include`) inside the Makefile. 

- CUDA headers not found  
The Makefile also tries to automatically find the directories with headers 
and libraries needed to compile the CUDA codes. If this fails, these can as well be 
inserted at the top of the Makefile.

- <cub/cub.cuh> not found  
The GPU code has only one external dependency, which is the CUB library. You 
can download it from: https://nvlabs.github.io/cub/index.html. The easiest 
way to install CUB is to add the top-level directory of where you've 
unpacked the CUB source codes to your ``$CPATH`` environment variable. For 
example, if you've unzipped the CUB sources into a directory called 
``/home/username/cub-version.number``, you can use 
``export CPATH=$CPATH:/home/username/cub-version.number/:`` to install CUB. In this way the 
nvcc compiler is able to find the CUB headers.

- Program tries to run GPU code when no GPU is present  
Note that the mex files for the GPU code will be produced by `make` if your 
machine has `nvcc`. Once the mex files for the GPU code have been produced, 
the MATLAB code will prefer to use the GPU functions instead of the CPU 
functions. If you have no GPU available but did compile the mex files for 
the GPU code, you will get errors and MATLAB will exit. To disable the use 
of the GPU code type `make clean` and use `make cpu` instead of `make` or 
`make all`.

## Developer instructions

The testing and tuning scripts for the GPU code have been written in Python, 
using [Kernel Tuner](https://github.com/benvanwerkhoven/kernel_tuner). This 
section provides information on how to setup a development environment. Note 
that these steps are only needed if you are interested in modifying the CUDA 
and C++ codes.

### Python 3

The tests for the GPU code and several of the C functions are written in 
Python, to run these a Python 3 installation is required. The easiest way to 
get this is using [Miniconda](https://conda.io/miniconda.html).

On Linux systems one could type the following commands to download and 
install Python 3 using Miniconda:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

All the required Python packages can be installed using the following command,
before you run it make sure that you have CUDA installed:
```
pip install -r requirements.txt
```

The tests can be run using ``nose``, for example by typing the following in 
the top-level or test directory:
```
nosetests -v
```

## Further questions

For any other questions regarding this software, you can  
- Search [issues section](https://github.com/berndrieger/alltoall3D/issues) or open a new topic there.
- Contact the authors:  
[(Hamidreza Heydarian)](https://github.com/hrheydarian) <H.Heydarian@tudelft.nl>   
[(Ben van Werkhoven)](https://github.com/benvanwerkhoven) <b.vanwerkhoven@esciencecenter.nl>  
[(Bernd Rieger)](https://github.com/berndrieger) <b.rieger@tudelft.nl>  

## Acknowledgement
Some files have been reused and adapted from the following sources:  

- [GMM registration](https://github.com/bing-jian/gmmreg)    
	[1] Jian, B. & Vemuri, B. C. Robust point set registration using Gaussian 
    mixture models. IEEE PAMI 33, 16331645 (2011).
- [Lie-algebraic averaging](http://www.ee.iisc.ac.in/labs/cvl/research/rotaveraging/)  
    [2] Govindu, V. Lie-algebraic averaging for globally consistent motion estimation. 
    In Proc. IEEE Conf. on Computer Vision and Pattern Recognition (2004).  
    [3] Chatterjee, A. Geometric calibration and shape refinement for 3D reconstruction
    PhD thesis. Indian Institute of Science (2015).
- [l1-magic optimization toolbox](https://statweb.stanford.edu/~candes/software/l1magic/)    
- [Natural-Order Filename Sort](https://nl.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort) 
