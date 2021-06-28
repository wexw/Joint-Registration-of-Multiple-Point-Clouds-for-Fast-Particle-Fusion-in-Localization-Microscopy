function [Particles1,V1,meanlocs] = PrepareInputParticles(Particles,CCD_pixelsize)
%Input:
%Particles: original experiment data
%N: number of particlesprocessed in this function
%Output
%Particles1: experimental data with suitable format, input of cost function
%V1: scaled points, input of JRMPC
% (C) Copyright 2019               CI Group
%  All rights reserved               Faculty of Applied Physics
%                                              Delft University of Technology
%                                              Lorentzweg 1
%                                             2628 CJ Delft
%                                             The Netherlands
%
% Wenxiu Wang , July 2020.
M=size(Particles,2);

V1={};

for i=1:M
    Particles1{1,i}.sigma=double(Particles{1,i}.sigma);
    Particles1{1,i}.points=double(Particles{1,i}.points);
    V1{i,1}=Particles1{1,i}.points'.*CCD_pixelsize;
end
%------------------Load raw data-----------------
VS = cellfun(@(V) size(V,2),V1,'uniformoutput',false);
VSM=cell2mat(VS);
%hist(VSM)
meanlocs=mean(VSM);
end