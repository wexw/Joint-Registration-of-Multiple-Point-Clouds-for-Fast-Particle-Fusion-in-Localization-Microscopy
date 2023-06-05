function [clusters,clustersAll]=H_Clustering(cost_norm,numclust,minClustSize)


%% Classification with Tree
mat = cost_norm;
distanceMethod = 'complete';        %single,average and complete are possible, preferably COMPLETE!

%numclust = 8;       %number of clusters to cut the tree in 
%minClustSize = K_c/10%minimum size of a cluster (determine this value by looking at the tree)

K_c=size(cost_norm,1);
threshold = 0;      %not important for now
order = 1:K_c;            %not important for now 
%plotTree = false;     %wether or not to plot the tree
names = MakeLeafNames(order,threshold); 

clustersAll = makeTreeAndCluster(mat, distanceMethod, names, numclust,1);       %all clusters found by Tree
if minClustSize==0
minClustSize = max( cellfun(@(V) length(V),clustersAll))/1.5;
end
clusters = clustersAll( cellfun(@(v) length(v),clustersAll) >= minClustSize);      
%{
clus = clusters;      
for c = 1:length(clus)
        

        TBindingSites2 = TBindingSites( clus{c}' );
        [SumDN,DmatrixN]=particleDistSum(TBindingSites2);
        AverageDN=SumDN/32;

end
%}
end