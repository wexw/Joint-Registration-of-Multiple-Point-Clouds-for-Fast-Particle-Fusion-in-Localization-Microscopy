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

## Operating system

The recommended OS to run this package is linux. For windows instalation please see the documentation to install requirements on this [link](https://github.com/imphys/smlm_datafusion3d).

## Installation of 2D version

1. Download the software or clone it using the code below: 

``` git clone https://github.com/wexw/Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy.git ```

2. Install the requirements. If you plan to use GPU, the correct version of the Cuda and CUB must be installed. 


Use the following commands to build the necessary libraries for this software:

```bash
cd Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy/2dCode
module load cuda/8.0
module load matlab/2019b

cmake3 .
make
make install
````


If you do not have a CUDA capable GPU you could install the software without GPU acceleration. Do note that the code will be orders of magnitude slower. Use the following commands to install the CPU-only version:
```bash

cmake -DUSE_GPU=OFF .
make
make install
```

CMake locates the code's dependencies and generates a Makefile. Make compiles the mex files and necessary shared libraries. Make install copies the mex files to the right directory for use in Matlab, it does not require priviledges to run.
For furthere information please see [here](https://github.com/imphys/smlm_datafusion3d) for the instruction.

3. To use this package it is required to have Mex files (mex_expdist) ready on your device. Each package (2D,3D) reuires its own MEX files to be used in Matlab. On the command line it is needed to include the path where the compiled libraries (.so) are located in LD_LIBRARY_PATH and in Matlab the path where the compiled mex-files (.mexa64)  are located needed to be added with addpath. Please note, once you need to use other package (2D,3D) you need to reassign new Mex files to the linux path. For adding the mex files into linux path you should put your software directory to <FPF_DIR> in the codes below.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<FPF_DIR>/2dCode

```

4. Start MATLAB and run demo codes.
## Example Usage
Examples of how to use the code on experimental and simulated data is shown in the MATLAB script DemoJCC.m for 2D cases and DemoJCC_3D.m for 3D cases.

## Acknowledgements

We reused and adapted some files from the following sources:

[1] Huijben, Teun APM, et al. "Detecting structural heterogeneity in single-molecule localization microscopy data." Nature communications 12.1 (2021): 3791.

[2] Heydarian, Hamidreza, et al. "Template-free 2D particle fusion in localization microscopy." Nature methods 15.10 (2018): 781-784.

[3] Heydarian, Hamidreza, et al. "3D particle averaging and detection of macromolecular symmetry in localization microscopy." Nature Communications 12.1 (2021): 2847.

[4] Haghparast, Sobhan, et al. <S.Haghparast@tudelft.nl>  Continuous heterogeneity detection https://gitlab.tudelft.nl/imphys/ci/chd
## Further questions
For any other questions regarding this software, you can
Search issues section or open a new topic there.
Contact the authors: 
(Wenxiu Wang) w.wang-21@tudelft.nl
(Hamidreza Heydarian) H.Heydarian@tudelft.nl
(Teun Huijben) teunhuijben@hotmail.com
(Sjoerd Stallinga) s.stallinga@tudelft.nl
(Bernd Rieger) b.rieger@tudelft.nl
