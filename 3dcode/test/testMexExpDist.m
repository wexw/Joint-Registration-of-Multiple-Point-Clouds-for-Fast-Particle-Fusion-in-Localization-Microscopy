%%
% (C) Copyright 2018-2020      
% Faculty of Applied Sciences
% Delft University of Technology
%
% Hamidreza Heydarian, November 2020.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%    http://www.apache.org/licenses/LICENSE-2.0

%% Initialization
clear all 
close all
clc

path(pathdef)
addpath(genpath('test'))
addpath(genpath('MATLAB'))
addpath(genpath('build/mex/'))      %remove everything from path, except the correct mex files 




%% Uncertainty matrix as mx2:
%this code assumes the new mex-files, where the indexing has been changed
%and the costfunction is normalized with respect to the sigmas. 
%If the old mex-files are used, the values will be different/wrong
%initialize particles (S and M have two locs)

S.points = [0,0,0;1,1,1];
S.sigma = [1 3;1 1];
M.points = [1 -1 1; 1 2 3];
M.sigma = [10 1; 3 4];

%test without rotation or a random rotation
%RM = eye(3); 
RM = rand(3);

%CPU
mex_expdist_cpu(S.points,M.points,correct_uncer(S.sigma),correct_uncer(M.sigma),RM)                

%GPU
mex_expdist(S.points,M.points,correct_uncer(S.sigma),correct_uncer(M.sigma),RM)  

%manual
manualCostFunction(S,M,RM)

%all three values should be the same 

%% Uncertainty matrix as mx9 (for M)
S.points = [0,0,0;1,1,1];
S.sigma = [1 3;1 1];
M.points = [1 -1 1; 1 2 3];
M.sigma = [reshape(diag([10,10,1])*rotx(32)*roty(14)*rotz(11),1,9); reshape(diag([3,3,4])*rotx(123)*roty(110)*rotz(13),1,9)] ;

RM = eye(3);
%RM = rand(3,3);     %RM does nothing when M.sigma is mx9

%CPU
% > not yet implemented
%GPU
mex_expdist(S.points,M.points,correct_uncer(S.sigma),correct_uncer(M.sigma),RM)
%manual
manualCostFunction(S,M,RM)      

%% including pairFitting (on real particles)
load('data/data.mat')
subParticles = cell(1,2);
for i=1:2
    subParticles{1,i}.points = particles{1,i}.coords(:,1:3);
    subParticles{1,i}.sigma = [particles{1,i}.coords(:,5).^2 particles{1,i}.coords(:,10).^2];
end

%GPU
M = subParticles{1}; 
S = subParticles{2}; 
[param,~,cost] = pairFitting3D_parallel(M, S,ones(1,size(S.points,1)), 0.1, 1, 1, 1);

M_trans.points = transform_by_rigid3d(M.points, param);
M_trans.sigma = M.sigma; 
RM = quaternion2rotation(param(1:4));

mex_expdist_cpu(S.points,M_trans.points,correct_uncer(S.sigma),correct_uncer(M_trans.sigma),RM)                
mex_expdist(S.points,M_trans.points,correct_uncer(S.sigma),correct_uncer(M_trans.sigma),RM)
manualCostFunction(S,M_trans,RM)
%these two values must be the same 


