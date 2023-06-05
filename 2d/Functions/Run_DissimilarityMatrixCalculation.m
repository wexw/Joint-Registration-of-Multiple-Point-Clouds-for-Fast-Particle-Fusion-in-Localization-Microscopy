disp('Start to Calculate Dissimilarity Matrix');
tic
%Prepare Output parameters in Parfor loop
clus_All=cell(Repeattimes,1);
ClusterView1_All=cell(Repeattimes,1);
cost_norm_All=cell(Repeattimes,1);
clusfull_All=cell(Repeattimes,1);

if USE_GPU_EXPDIST==1
    disp('Start to calculate Dissimilarity Matrix by GPU');
    for J_C1=1:Repeattimes%%%%%%%%%%%%%%
        TV_p=TV_All{ J_C1,1};      %Input transport particles from the JRMPC
        % Calculate the cost matrix
   
        [ ~,cost_norm_All{J_C1,1}]=Cost_2D_gpu(TV_p,Particles1,CCD_pixelsize);
    
    end
    time_Dissimilarity=toc
    disp(' Dissimilarity Matrix Calculation finished by GPU');
else


    tic
    if isempty(gcp('nocreate'))
        parpool('local',40);
    end
    disp('Start to calculate Dissimilarity Matrix by CPU');

    for J_C1=1:Repeattimes%%%%%%%%%%%%%%
        TV_p=TV_All{ J_C1,1};      %Input transport particles from the JRMPC
        % Calculate the cost matrix
        [~,  cost_norm_All{ J_C1,1}]=Cost_2D_cpu(TV_p,Particles1,CCD_pixelsize);
    end

    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if ~isempty(poolobj)
        delete(poolobj);
    end

    time_Dissimilarity=toc
    disp(' Dissimilarity Matrix Calculation finished by CPU');
end
