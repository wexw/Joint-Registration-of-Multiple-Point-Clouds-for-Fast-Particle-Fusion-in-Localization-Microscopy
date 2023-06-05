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

if K_origin<100
    K_origin=100;
end
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
    %For 2D case
az1 = 2*pi*rand(1,K_origin);
el1 = 2*pi*rand(1,K_origin);
Xin1 = [cos(az1); sin(el1)];
Xin1 = Xin1*CCD_pixelsize;
%------------------Generate Xin1-----------------
%------------------JRMPC 2D-----------------
[R ,t,X] = jrmpc_2D_origin(V1,Xin1,'maxNumIter',100, 'epsilon', 1e-5,'gamma',0.01);   %JRMPC

%{
Drawmatrix=[];
Drawmatrix=X;
x=[];
y=[];

x=Drawmatrix(1,:);
y=Drawmatrix(2,:);
figure()
plot(x,y,'.');    
xlim([-60 60])
ylim([-60 60])
%}
%Use meanshift to get the distribution of centers
dataPts=X;
[clustCent,~,cluster2dataCell] = MeanShiftCluster(dataPts,bandwidth);
%{
Drawmatrix=[];
Drawmatrix=clustCent;
x=[];
y=[];

x=Drawmatrix(1,:);
y=Drawmatrix(2,:);
figure()
plot(x,y,'.');    
xlim([-60 60])
ylim([-60 60])
%}
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
    %{
    Drawmatrix=[];
Drawmatrix=K_final;
x=[];
y=[];

x=Drawmatrix(1,:);
y=Drawmatrix(2,:);
figure()
plot(x,y,'.');    
xlim([-60 60])
ylim([-60 60])
    %}
    K=size(K_final,2);
else
    K=8;
end