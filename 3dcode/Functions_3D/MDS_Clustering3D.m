%Use Multiple Dimention scaling to classify
%input
%cost_norm: the normalized cost matrix from Cost_2D or Cost_3D
%numClust: the number of clusters after classification
%minClustSize: the minimum number of particles included in the cluster that
%the custer can be seen as an effective cluster
%Output
%cluster:effective cluster
%clusterAll: all cluster including the cluster which has less than
%minClustSize particles
function [clusters,clustersAll]=MDS_Clustering3D(cost_norm,numClust,minClustSize,dimention)
if nargin < 4
   dimention=80;
end
D = cost_norm+cost_norm';       %use costnorm to try your output
D = max(D(:))-D; 
D = D - diag(diag(D)); 

opts=   struct('TolFun',1e-3, 'TolX',1e-3);
mds = mdscale(D,dimention,'Criterion','metricstress','Options',opts);

[c] = kmeans(mds(:,1:end),numClust,'replicates',1000);

for i = 1:max(c)
    clustersAll{i} = find(c==i);
end

clusters = clustersAll( cellfun(@(v) length(v),clustersAll) >= minClustSize);      
% plotting
%figure(); scatter3(mds(:,1),mds(:,2),mds(:,3),'o','filled'); title 'MDS'
%figure(); scatter3(mds(:,1),mds(:,2),mds(:,3),[],c,'o','filled'); title 'clusters'
end