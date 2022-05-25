function [Particles1,V1,meanlocs,points,sigma] = ProcessedExperimentalDataPraticle(Particles,CCD_pixelsize,N,minLo,maxLo)
%Input:
%Particles: original experiment data
%N: number of particlesprocessed in this function
%limit: only output particles which has more than limit points
%Output
%Particles1: experimental data with suitable format, input of cost function
%V1: scaled points, input of JRMPC
% (C) Copyright 2019               QI Group
%  All rights reserved               Faculty of Applied Physics
%                                              Delft University of Technology
%                                              Lorentzweg 1
%                                             2628 CJ Delft
%                                             The Netherlands
%
% Wenxiu Wang , July 2020.
ParticleNumber=size(Particles,2);

V0={};
if N>ParticleNumber
    N=ParticleNumber
end
for i=1:N
    Particles0{1,i}.sigma=double(Particles{1,i}.sigma);
    Particles0{1,i}.points=double(Particles{1,i}.points);
    V0{i,1}=Particles0{1,i}.points'.*CCD_pixelsize;
    points0{i,1}=double(Particles{1,i}.points);
    sigma0{i,1}= double(Particles{1,i}.sigma);
end
%------------------Load raw data-----------------
VS = cellfun(@(V) size(V,2),V0,'uniformoutput',false);
VSM=cell2mat(VS);
%hist(VSM)
meanlocs=mean(VSM);
position_large=find(VSM>minLo);
position_small=find(VSM<maxLo);
position_common=intersect( position_large, position_small);
Ns=size(position_common,1);
for i=1:Ns
    p=position_common(i);
    V1{i,1}=Particles0{1,p}.points'.*CCD_pixelsize;
    Particles1{1,i}.sigma=Particles0{1,p}.sigma;
    Particles1{1,i}.points=Particles0{1,p}.points;
    points{i,1}=  points0{p,1};
    sigma{i,1}= sigma0{p,1};
end
end

