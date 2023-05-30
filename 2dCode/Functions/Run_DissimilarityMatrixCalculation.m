disp('Start to Calculate Dissimilarity Matrix');
tic
%Prepare Output parameters in Parfor loop
clus_All=cell(Repeattimes,1);
ClusterView1_All=cell(Repeattimes,1);
cost_norm_All=cell(Repeattimes,1);
clusfull_All=cell(Repeattimes,1);

if USE_GPU_EXPDIST==1
    disp('Start to calculate Dissimilarity Matrix by GPU');

parfor J_C1=1:Repeattimes%%%%%%%%%%%%%%
    TV_p=TV_All{ J_C1,1};      %Input transport particles from the JRMPC
    % Calculate the cost matrix
    tic
    [~,  cost_norm_All{ J_C1,1}]=Cost_2D_gpu(TV_p,Particles1,CCD_pixelsize); 
toc
end
time_Dissimilarity=toc;
disp(' Dissimilarity Matrix Calculation finished by GPU');
else
      disp('Start to calculate Dissimilarity Matrix by CPU');
parfor J_C1=1:Repeattimes%%%%%%%%%%%%%%
    TV_p=TV_All{ J_C1,1};      %Input transport particles from the JRMPC
    % Calculate the cost matrix
    tic
    [~,  cost_norm_All{ J_C1,1}]=Cost_2D_cpu(TV_p,Particles1,CCD_pixelsize); 
toc
end
time_Dissimilarity=toc;
disp(' Dissimilarity Matrix Calculation finished by CPU');

end
