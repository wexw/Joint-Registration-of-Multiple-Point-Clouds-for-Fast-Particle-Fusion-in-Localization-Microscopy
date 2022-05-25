disp('Start to Run Connection');
tic
%------------------Connect------------------
ClusterViewMatrix=cell2mat(ClusterView1_All);
[maxrow,maxrowposition]=max(ClusterViewMatrix(:,3));
maxR=ClusterViewMatrix(maxrowposition,1);
maxC=ClusterViewMatrix(maxrowposition,2);
clus_initial=clus_All{maxR,1}{1,maxC};                   %%%%%%%%%%%%%main cluster to connect all possible clusters
[TVPick, ReorderParticles, Rpick,tpick] = ConnectParticlesInDifferentClusterWithSigma(limit,Repeattimes,clus_All,clus_initial,maxR,R_All,t_All,TV_All,Particles1);
%------------------Connect------------------

toc
Time_Connect=toc;
disp('Connected finished');