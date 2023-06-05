%----------------Classify all clusters-----------
disp('Start to Classify');
tic
clus_All=cell(Repeattimes,1);
ClusterView1_All=cell(Repeattimes,1);
clusfull_All=cell(Repeattimes,1);
ClusterView=[];
for J_C1=1:Repeattimes%%%%%%%%%%%%%%

    % TV_p=TV_All{ J_C1,1};
    MatrixAfterAll2all_norm= cost_norm_All{ J_C1,1};
    %  [~, MatrixAfterAll2all_norm]=Cost_2D(TV_p,Particles1,CCD_pixelsize);
    clusterfull=[];
    if cluster_mode==1
        [ clus,clusterfull]=H_Clustering(MatrixAfterAll2all_norm,nc,minClustSize);
    else
        [clus,clusterfull]=MDS_Clustering3D(MatrixAfterAll2all_norm,nc,minClustSize,dimention);% First Classify
    end
    clusfull_All{ J_C1,1}=clusterfull;
    clus_All{ J_C1,1}=clus;
end

for J_C1=1:Repeattimes%%%%%%%%%%%%%%%%%%%%
    clusterfull=clusfull_All{ J_C1,1};
    clus= clus_All{ J_C1,1};

    if size(clus,2)==0
        ClusterView=[0,0]
        clus=0;
    elseif size(clus,2)==1
        ClusterView=[J_C1,1,size(cell2mat(clus),1)];
    else
        for ic=1:1:size(clus,2)
            ClusterView(ic,:)=[J_C1,ic,size(clus{1,ic},1)];
        end
    end
    ClusterView1_All{ J_C1,1}=ClusterView;
end


disp('Classification finished');
time_Classification=toc;
