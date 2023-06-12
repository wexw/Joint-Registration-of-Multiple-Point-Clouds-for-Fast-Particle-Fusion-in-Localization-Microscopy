%New 3D JRMPC+Classify+Connect
clc
clear
%close all

path_mex_matlab1 = genpath('build/mex/');
path_mex_matlab2 = genpath('build/mex/');
path_matlab = genpath('MATLAB');
path_function=genpath('Functions_3D');

addpath(path_mex_matlab1)
addpath(path_mex_matlab2)
addpath(path_matlab)
addpath(path_function)


% CPU/GPU settings (CPU = 0, GPU = 1)
USE_GPU_GAUSSTRANSFORM = 1;
USE_GPU_EXPDIST = 1;

cluster_mode=2;          %1: H,2:MDS
Repeattimes=2;           %Number of independent JCC progress for a single NUP model
sizeview=120; %radius of Field of interest
limit=0;
load('SimulatedNup107With65%DOL.mat')
% particles has 1xN cell, with 1X1struct in each, the unit in the cells are nm

Particles=particles;



%dol=0.4
N=size(Particles,2);

for i=1:N % make all the particles around (0,0,0)
    Particles{1,i}.points(:,1:3)= Particles{1,i}.points(:,1:3)-mean(Particles{1,i}.points(:,1:3));
    max_nm(i)=max(max(Particles{1,i}.points(:,1:3)));
end
Cscale=round(max(max_nm));
for i=1:N % to make all the particles inside the unit box
    Particles{1,i}.points(:,1:3)= Particles{1,i}.points(:,1:3)./Cscale-mean(Particles{1,i}.points(:,1:3)./Cscale);
    Particles{1,i}.sigma(:,1:2)=   Particles{1,i}.sigma(:,1:2)./Cscale;
end
minLo=10;
maxLo=10000;%maximumlocalization in one particle
%Prepare the input for the JRMPC
[~,V1,meanlocs,~,~] = ProcessedExperimentalDataPraticle(Particles,1,N,minLo,maxLo);
[Particles1,V2,~,~,~] = ProcessedExperimentalDataPraticle(Particles,Cscale,N,minLo,maxLo);


M=size(V1,1);
raster_spacing=40;       % expected resolution [nm]
K = 18;
%K=K_estimate(V2,meanlocs,20,0.5*raster_spacing, CCD_pixelsize);%Estimate number of Gaussian models

Sa=sqrt(1000)./Cscale;% Size of initial Gaussian components


%Step1: Joint registration
Run_JRMPC3D




Run_DissimilarityMatrixCalculation3D


nc=2; %number of clusters
cluster_mode=2; %2: MDS clustering, 1: Hierarchical agglomerative clustering approach as alternative to MDS
minClustSize=M/(nc+1);  %The minumum number of  the particle
%Step3: Classification

dimention=30;%round(N/6)*10; If classification doesn't work, try different dimention values
Run_Classification3D


%Step4: Connection
Run_Connection    %no need to change anything inside

time_All=Time_JRMPC+time_Dissimilarity+Time_Connect+time_Classification; % Total Computational Time

% effectiveParticlesNumber=size(TVPick,1);%number of particles included in the final results.
% data2c=cell2mat(TVPick')';

%filenameTVpick2=[ 'UncorrectedNPCK=' num2str(K)    'S=' num2str(Sa) 'nc=' num2str(nc) 'Size=' num2str(effectiveParticlesNumber) 'dimention=' num2str(dimention) '.avi' ];

Draw3DAviProjection
% Draw3DFigure

