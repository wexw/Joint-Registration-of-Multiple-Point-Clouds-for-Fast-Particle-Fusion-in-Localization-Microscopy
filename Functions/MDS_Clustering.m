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
function [clusters,clustersAll]=MDS_Clustering(cost_norm,numClust,minClustSize)

D = cost_norm+cost_norm';
D = max(D(:))-D+1; 
D = D - diag(diag(D)); 
mds = mdscale(D,30,'Criterion','metricstress');     %you can vary the number of dimension (30) and the criterion, but for me these values always seem to work

% cluster the MDS space with k-means

[c] = kmeans(mds,numClust,'replicates',100);    % I do at least 100 replicates to make it more robust, because k-means depends on the initialization of the means 
for i = 1:max(c)
    clustersAll{i} = find(c==i);
end
clusters = clustersAll( cellfun(@(v) length(v),clustersAll) >= minClustSize);      
% plotting
%figure(); scatter3(mds(:,1),mds(:,2),mds(:,3),'o','filled'); title 'MDS'
%figure(); scatter3(mds(:,1),mds(:,2),mds(:,3),[],c,'o','filled'); title 'clusters'
end