clc
clear
dip_initialise % Initialize Diplib
%%Load data


path_function=genpath('Functions/')
addpath(path_function)

USE_GPU_GAUSSTRANSFORM = 1;
USE_GPU_EXPDIST = 1;

load('tud_442particles_dol50.mat')
CCD_pixelsize=130;     %CCD_pixelsize for the data
Repeattimes=2;           %Number of initializations of the JRMPC
raster_spacing=5;        
ParticleNumber=size(particles,2);%Number of input particles
M=ParticleNumber;

[Particles1,V1,meanlocs] =PrepareInputParticles(particles,CCD_pixelsize);
%Particles1:Particles with uncertainties and coordinates
%V1:particles' coordinates
%mean locs: average number of localizations in the particles



%Step1: Joint registration
Run_JRMPC  %We can change Xin and JRMPC parameters inside

%Step2: Dissimilarity Matrix Calculation
Run_DissimilarityMatrixCalculation %no need to change anything inside


%Two parameters defined for the classification
nc=2; %number of clusters
cluster_mode=2; %2: MDS clustering, 1: Hierarchical agglomerative clustering approach as alternative to MDS
minClustSize=M/(nc+1);  %The minumum number of  the particle 
%Step3: Classification
Run_Classification %We can set different nc,minClustSize before Run_Classification

%Step4: Connection
Run_Connection    %no need to change anything inside

time_All=time_JRMPC+time_Dissimilarity+time_Connection+time_Classification; % Total Computational Time

%Step5: FRC resolution measurement
%%FRC resolution measurement installation needed
%%Used to calculate FRC resolution See https://github.com/imphys/FRC_resolution
%Run_MainClusterFRC

%Draw the density Plot
XX=[];
pixelnum=300;%number of pixels in the density plot
diameter=100;%[nm]
XX=cell2mat(TVPick')';%final localizations
nef=size(TVPick,1);%number of particles in the final reconstruction
lo=size(XX,1);%number of localizationss in the final reconstruction
%Diplib needed
titlename= ['N=' num2str(nef) 'time=' num2str(time_All) ];

[dip, Z] = visualizeCloud2DW(XX,pixelnum,diameter,titlename);%Draw the density plot
