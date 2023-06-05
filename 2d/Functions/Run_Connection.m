disp('Start to Run Connection');
tic
ClusterViewMatrix=cell2mat(ClusterView1_All);
[maxrow,maxrowposition]=max(ClusterViewMatrix(:,3));
maxR=ClusterViewMatrix(maxrowposition,1); %The cluster with most particles from the maxR_th initialization
maxC=ClusterViewMatrix(maxrowposition,2);%The cluster with most particles from the maxC_th cluster in the 
clus_initial=clus_All{maxR,1}{1,maxC}; %Initialize the main cluster
%Connect the Particles in the different Clusters
[TVPick,Rpick,tpick] = ConnectParticlesInDifferentCluster(Repeattimes,clus_All,clus_initial,maxR,R_All,t_All,TV_All);
time_Connection=toc;
disp('Connection finished');