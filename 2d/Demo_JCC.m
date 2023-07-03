clc
clear
path_function=genpath(['Functions/']);
addpath(path_function)

load('tud_120particles_dol80.mat')

USE_GPU_EXPDIST = 1;

if USE_GPU_EXPDIST==0
    disp(['Please set the number of cores used here: Core_numer=...']);
    Core_number=24;%%%You can set this number based on your hardware 
end

CCD_pixelsize=130;     %CCD_pixelsize for the data
Repeattimes=2;           %Number of initializations of the JRMPC
raster_spacing=5;        
ParticleNumber=size(particles,2);%Number of input particles
M=ParticleNumber;


%Two parameters defined for the classification
nc=2; %number of clusters
cluster_mode=2; %2: MDS clustering, 1: Hierarchical agglomerative clustering approach as alternative to MDS
minClustSize=M/(nc+1);%The minumum number of the particle in a 'good' cluster


[Particles1,V1,meanlocs] =PrepareInputParticles(particles,CCD_pixelsize);
%Particles1:Particles with uncertainties and coordinates
%V1:particles' coordinates
%mean locs: average number of localizations in the particles

%estimate number of GMMcenters K for JRMPC
disp(['Start to estimate the number of Gaussian components in GMM']);
K =K_estimate(V1,meanlocs,20,0.8*raster_spacing,CCD_pixelsize); % Estimate the number of the GMM centers; We can also set the value.
disp(['The estimated number of Gaussian components in GMM is ' num2str(K)]);
disp(['You can also set the K with prior knowledge']);

%Step1: Joint registration
Run_JRMPC  %We can change Xin and JRMPC parameters inside
           %Usually the default setting works well
%Step2A: Classification:Dissimilarity Matrix Calculation
Run_DissimilarityMatrixCalculation %no need to change anything inside

%Step2B: Classification
Run_Classification %We can set different nc,minClustSize before Run_Classification

%Step3: Connection
Run_Connection    %no need to change anything inside

time_All=time_JRMPC+time_Dissimilarity+time_Connection+time_Classification; % Total Computational Time

%Step5: FRC resolution measurement
%%FRC resolution measurement installation needed
%%Used to calculate FRC resolution See https://github.com/imphys/FRC_resolution
%Run_MainClusterFRC

%Draw the Plot by renderprojection
XX=[];
XX=cell2mat(TVPick')';%final localizations
fig=renderprojection2D(XX(:,1),XX(:,2), [-80 80],[-80 80],1,1, 1);
axis off
axis equal
