disp('Start to Classify');
tic

parfor J_C1=1:Repeattimes%%%%%%%%%%%%%%
    TV_p=TV_All{ J_C1,1};
    MatrixAfterAll2all_norm=  cost_norm_All{ J_C1,1};
    if cluster_mode==1
        %Hierarchical agglomerative clustering approach as alternative to MDS
        %See Detecting structural heterogeneity in single-molecule localization microscopy data by T.A.P.M. Huijben
        [clus,clusterfull]=H_Clustering(MatrixAfterAll2all_norm,nc,minClustSize); 
    else
        [clus,clusterfull]=MDS_Clustering(MatrixAfterAll2all_norm,nc,minClustSize);% MDS clustering
    end
    
    %Prepare the output of the classification
    if size(clus,2)==0
        ClusterView=[0,0]
        clus=0;
    else
        for ic=1:1:size(clus,2)
            ClusterView(ic,:)=[J_C1,ic,size(clus{1,ic},1)]; %[iteration, cluster sequence, number of particles in the cluster]
            
        end
    end
    
    clus_All{ J_C1,1}=clus;
    ClusterView1_All{ J_C1,1}=ClusterView;
    clusfull_All{ J_C1,1}=clusterfull;
end
%------------------Classify 1st------------------
disp('Classification finished');
time_Classification=toc;