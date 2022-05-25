# Joint-Registration-of-Multiple-Point-Clouds-for-Fast-Particle-Fusion-in-Localization-Microscopy
Matlab code for Joint Registration of Multiple Point Clouds for Fast Particle Fusion in Localization Microscopy.
This method is fast, avoids symmetry artifacts, applies to 2D and 3D datasets, and reconstructs poor data with a limited number of particles, a low density of labelling and a large localization uncertainty.

#MATLAB
The code is written in MATLAB, and tested to work in MATLAB R2018-R2020.

The DIPImage toolbox for MATLAB is required, please see http://www.diplib.org for installation instructions.


2D Mex file: mex_expdist for 2D code in Detecting Structural Heterogeneity in Single-Molecule Localization Microscopy Data is required for 2D Classification, please see https://github.com/imphys/smlm_classification2d
3D Mex file:mex_expdist for 3D code in Software for Template-Free 3D Particle Fusion in Localization Microscopy is required for 3D classification, please see https://github.com/imphys/smlm_datafusion3d

#Example Usage
Examples of how to use the code on experimental and simulated data is shown in the MATLAB script DemoJCC.m for 2D cases and DemoJCC_3D.m for 3D cases.

#Further questions
For any other questions regarding this software, you can
Search issues section or open a new topic there.
Contact the authors: 
(Wenxiu Wang) w.wang-21@tudelft.nl
(Hamidreza Heydarian) H.Heydarian@tudelft.nl
(Teun Huijben) teunhuijben@hotmail.com
(Sjoerd Stallinga) s.stallinga@tudelft.nl
(Bernd Rieger) b.rieger@tudelft.nl
