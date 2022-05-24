function K = K_estimate(V0,meannlocs,N_set,bandwidth,CCD_pixelsize)
%Input
%V0 is the orginal particles
%meanlocs is the mean localization in the particles
%N_set is the number of particles used to estimate K
%bandwidth is the distance between each two adjecent centers.
%Output
%K is the number of GMM centers estimated for V0
% (C) Copyright 2019               QI Group
%  All rights reserved               Faculty of Applied Physics
%                                              Delft University of Technology
%                                              Lorentzweg 1
%                                             2628 CJ Delft
%                                             The Netherlands
%
% Wenxiu Wang , July 2020.


%set orignal number of K as meanlocs
K_origin=ceil(meannlocs);


%CCD_pixelsize=130;
N=size(V0,1);
%limit number of particles used no larger than N_set
if N>N_set
    Ns=N_set;
    for i=1:Ns
        V1{i,1}=V0{i,1};
    end
else
    V1=V0;
end

Dimension=size (V1{1,1},1);
if Dimension==2
    %For 2D case
az1 = 2*pi*rand(1,K_origin);
el1 = 2*pi*rand(1,K_origin);
Xin1 = [cos(el1); sin(el1)];
Xin1 = Xin1*CCD_pixelsize;
%------------------Generate Xin1-----------------
%------------------JRMPC 2D-----------------
[R ,t,X,S,a] = jrmpc_2D(V1,Xin1,'maxNumIter',500, 'epsilon', 1e-5,'gamma',0.01);   %JRMPC
else 
    %For 3D case
  az= 2*pi*rand(1,K_origin);
el= 2*pi*rand(1,K_origin);
Xin = [cos(az).*cos(el); sin(el); sin(az).*cos(el)];
Xin1 = Xin*CCD_pixelsize;
%------------------Generate Xin1-----------------
%------------------JRMPC 2D-----------------
[R ,t,X,S,a] = jrmpc(V1,Xin1,'maxNumIter',500, 'epsilon', 1e-5,'gamma',0.01);   %JRMPC
end

%Use meanshift to get the distribution of centers
dataPts=X;
[clustCent,~,cluster2dataCell] = MeanShiftCluster(dataPts,bandwidth);

VS = cellfun(@(V) size(V,2),cluster2dataCell,'uniformoutput',false);
VSM=cell2mat(VS);
limit=1;
position_large=find(VSM>limit);
Nl=size(position_large,1);
if Nl>0
    for i=1:Nl
        p=position_large(i);
        K_final(:,i)=clustCent(:,p);
    end
    
    K=size(K_final,2);
else
    K=1;
end