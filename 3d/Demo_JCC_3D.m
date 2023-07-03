%3D JRMPC+Classify+Connect
clc
clear
%close all

path_mex_matlab1 = genpath('build/mex/');
path_function=genpath('Functions_3D');

addpath(path_mex_matlab1)
addpath(path_function)


% CPU/GPU settings (CPU = 0, GPU = 1)

USE_GPU_EXPDIST =1;
if USE_GPU_EXPDIST==0
    disp(['Please set the number of cores used here: Core_numer=']);
    Core_number=40;%%%You can set this number based on your hardware
end
CCD_pixelsize=130;
cluster_mode=2;          %1: H,2:MDS
Repeattimes=2;           %Number of independent JCC progress for a single NUP model
sizeview=120; %radius of Field of interest
limit=0;
load('SimulatedNup107With65%DOL.mat')
% particles has 1xN cell, with 1X1struct in each, the unit in the cells are nm

Particles=particles;

Prapare_Input_particles
M=size(V1,1);
raster_spacing=30;       % expected resolution [nm]
%Step1: Joint registration
%Prepare parameters for JRMCP
disp(['Start to estimate the number of Gaussian components in GMM']);


% K = 18;
K=K_estimate(V2,meanlocs*2,30,raster_spacing, CCD_pixelsize);%Estimate number of Gaussian models
disp(['The estimated number of Gaussian components in GMM is ' num2str(K)]);
disp(['You can also set the K with prior knowledge']);
Sa=sqrt(1000)./Cscale;% Size of initial Gaussian components
nc=2; %number of clustersminClustSize=M/(nc+1);  %The minumum number of  the particle
dimention=30;%round(N/6)*10; If classification doesn't work, try different dimention values
minClustSize=M/(nc+1);
Run_JRMPC3D


%Step2A: Classification:Dissimilarity Matrix Calculation
Run_DissimilarityMatrixCalculation3D


%Step2B: Classification
Run_Classification3D


%Step3 Connection
Run_Connection    %no need to change anything inside

time_All=Time_JRMPC+time_Dissimilarity+Time_Connect+time_Classification; % Total Computational Time

% effectiveParticlesNumber=size(TVPick,1);%number of particles included in the final results.
% data2c=cell2mat(TVPick')';

%filenameTVpick2=[ 'UncorrectedNPCK=' num2str(K)    'S=' num2str(Sa) 'nc=' num2str(nc) 'Size=' num2str(effectiveParticlesNumber) 'dimention=' num2str(dimention) '.avi' ];
%Draw the Plot by renderprojection
Draw3DAviProjection

%%%You can also Draw the Plot by density plot
%USE_GPU_GAUSSTRANSFORM = 0;
%dip_initialise  %%diplib sould be loaded to use density plot function
%Draw3DFigure%
%%USE_GPU_GAUSSTRANSFORM is used here
