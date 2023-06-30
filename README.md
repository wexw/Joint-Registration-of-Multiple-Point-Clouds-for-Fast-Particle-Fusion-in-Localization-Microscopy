# Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy
Matlab code for "Joint registration of multiple point clouds for fast particle fusion in localization microscopy." Bioinformatics 38.12 (2022): 3281-3287.

The main code is written in MATLAB and some of the compute-intensive kernels have been written in CUDA and C++.

## Requirements

* [MATLAB](https://www.mathworks.com/products/matlab.html) (version>2019b)
* [Dipimage](http://www.diplib.org/) (version>2.9)
* [Cuda](https://developer.nvidia.com/cuda-downloads) (version>11.2)
* [GCC](https://gcc.gnu.org/) (version>5.5.0)
* [CUB](https://nvlabs.github.io/cub/) (version>1.8.0)
* [Cmake](https://cmake.org/) (version>3.14.3)


Running the demo script in Matlab requires the installation of the following toolboxes in Matlab:
* Image Processing Toolbox
* Parallel Computing Toolbox
*	Statistics and Machine Learning Toolbox
* Matlab Compiler


## Operating system

The recommended OS to run this package is Linux. For windows instalation please see the documentation to install requirements on this [link](https://github.com/imphys/smlm_datafusion3d).

## Installation of 2D version

1. Download the software or clone it using the code below: 

``` git clone https://github.com/wexw/Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy.git ```

2. Install the requirements. If you plan to use GPU, the correct version of the Cuda and CUB must be installed. 


Use the following commands to build the necessary libraries for this software:

```bash
cd Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy/2d

cmake3 .
make
make install
````

Please ensure the correct installation of the CUB library. If CMake cannot locate the CUB library, you should manually add the location of the CUB library in your CMakeLists.txt file, like so:
````bash
include_directories(../../cub-1.8.0)
````
If you do not have a CUDA capable GPU you could install the software without GPU acceleration. Do note that the code will be orders of magnitude slower. Use the following commands to install the CPU-only version:
```bash

cmake -DUSE_GPU=OFF .
make
make install
```

CMake locates the code's dependencies and generates a Makefile. Make compiles the mex files and necessary shared libraries. Make install copies the mex files to the right directory for use in Matlab, it does not require priviledges to run.
For furthere information please see [here](https://github.com/imphys/smlm_datafusion3d) for the instruction.

4.  To use this package it is required to have Mex files (`mex_expdist`) ready on your device. These are automatically generated after compilation of the c files. Each package (2D,3D) requires its own MEX files to be used in Matlab. It is needed to include the <compiled libraries (`.so`) path> in `LD_LIBRARY_PATH`. Afterwards you should add the <location of compiled mex-files (`.mexa64`)> to Matlab path. This can be done using `addpath` command in Matlab. Please note, once you need to use other package (2D,3D) you need to reassign new Mex files to the linux path. For adding the mex files into linux path you should put your software directory to `<Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy>` in the codes below.

   ```bash
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy>/2d
   ```

5. Start MATLAB and set `USE_GPU_GAUSSTRANSFORM=0' and `USE_GPU_EXPDIST=0' to determine if you want to use the CPU or set `USE_GPU_GAUSSTRANSFORM=1' and `USE_GPU_EXPDIST=1' to use the GPU.



4. Start MATLAB and run demo codes.

### Example Usage
Examples of how to use the code on experimental data is shown in the MATLAB script DemoJCC.m for the 2D case.

## Installation of 3D version

### Compile the code
In the following

- <SOURCE_DIRECTORY> is the root directory of the sources
- <CUB_DIRECTORY> is the root directory of the downloaded [CUB library](https://nvlabs.github.io/cub/) sources


Use the following commands to build the necessary libraries for this software:

```bash
cd <SOURCE_DIRECTORY>/3d
mkdir build
cd build
cmake -DCUB_ROOT_DIR=CUB_DIRECTORY <SOURCE_DIRECTORY>/3d
make
````
### Use the code
To use this package it is required to have Mex files (mex_expdist) ready on your device. Each package (2D,3D) reuires its own MEX files to be used in Matlab. On the command line it is needed to include the path where the compiled libraries (.so) are located in LD_LIBRARY_PATH and in Matlab the path where the compiled mex-files (.mexa64)  are located needed to be added with addpath. Please note, once you need to use other package (2D,3D) you need to reassign new Mex files to the linux path. For adding the mex files into linux path you should put your software directory to Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy> in the codes below.
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy/2d
````
Next, we need to locate the built libraries for MATLAB:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy/3d/build/mex/
``` 
Then, run MATLAB and set `USE_GPU_GAUSSTRANSFORM=0' and `USE_GPU_EXPDIST=0' to determine if you want to use the CPU or set `USE_GPU_GAUSSTRANSFORM=1' and `USE_GPU_EXPDIST=1' to use the GPU. 

### Example Usage
Examples of how to use the code on simulated data is shown in the MATLAB script Demo_JCC_3D.m for the 3D case.

## Troubleshooting

1. After runing cmake. getting error message that Matlab was not found. This is because in the installation of Matlab, it is necessary to create symbolic links to MATLAB scripts. This can be done by checking the box ‘Create symbolic links to MATLAB scripts in:’ during installation. Or after installation, it can be done with the command:
```bash 
sudo ln -s $MATLAB/bin/matlab /usr/local/bin/matlab 
``` 
where $MATLAB is the location where matlab is installed.

2. The LD_LIBRARY_PATH is only set inside the current terminal. There fore it is recommended to open matlab in the same terminal other wise .so files are not found.
Another solution is to change the settings from my PC to not reset the LD_LIBRARY_PATH.



## Acknowledgements
Thanks to Ronald Ligteringen for his help in testing and compiling the code.
Thanks to Sobhan Haghparast for his efforts in improving the mex files.
Thanks to Isabel Droste for her careful testing and valuable suggestions for improving the code.


We reused and adapted some files from the following sources:

[1] Haghparast, Sobhan, et al. <S.Haghparast@tudelft.nl>  Continuous heterogeneity detection https://gitlab.tudelft.nl/imphys/ci/chd

[2] Heydarian, Hamidreza, etThanks to Isabel Droste for her careful testing and valuable suggestions for improving the code. al. "Template-free 2D particle fusion in localization microscopy." Nature methods 15.10 (2018): 781-784.

[3] Heydarian, Hamidreza, et al. "3D particle averaging and detection of macromolecular symmetry in localization microscopy." Nature Communications 12.1 (2021): 2847.

[4] Huijben, Teun APM, et al. "Detecting structural heterogeneity in single-molecule localization microscopy data." Nature communications 12.1 (2021): 3791.

## Further questions
For any other questions regarding this software, you can
Search issues section or open a new topic there.
Contact the authors: 
(Wenxiu Wang) w.wang-21@tudelft.nl
(Hamidreza Heydarian) H.Heydarian@tudelft.nl
(Teun Huijben) teunhuijben@hotmail.com
(Sjoerd Stallinga) s.stallinga@tudelft.nl
(Bernd Rieger) b.rieger@tudelft.nl
