disp('Start to Calculate Dissimilarity Matrix');
tic
%Prepare Output parameters in Parfor loop
clus_All=cell(Repeattimes,1);
ClusterView1_All=cell(Repeattimes,1);
cost_norm_All=cell(Repeattimes,1);
clusfull_All=cell(Repeattimes,1);
parfor J_C1=1:Repeattimes%%%%%%%%%%%%%%
    TV_p=TV_All{ J_C1,1};      %Input transport particles from the JRMPC
    % Calculate the cost matrix
    [~,  cost_norm_All{ J_C1,1}]=Cost_2D(TV_p,Particles1,CCD_pixelsize); 
end
time_Dissimilarity=toc;
disp(' Dissimilarity Matrix Calculation finished');