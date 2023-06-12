%...-----------------Calculate Cost Function-----------

tic
%V2's unit is nm, V2 is the original input particles' localizations
cost_norm_All=cell(Repeattimes,1);
if USE_GPU_EXPDIST ==1
    disp('Start to calculate Dissimilarity Matrix by GPU');
    for J_C1=1:Repeattimes%%%%%%%%%%%%%%
        TV_p=TV_All{ J_C1,1};
        R_p=R_All{ J_C1,1};
        t_p=t_All{ J_C1,1};
        [~, MatrixAfterAll2all_norm]=Cost_3DReal_gpu(V2,TV_p,R_p,t_p,Particles1,Cscale);
        cost_norm_All{J_C1,1}=MatrixAfterAll2all_norm;
    end
    disp(' Dissimilarity Matrix Calculation finished by GPU');
else
    disp('Start to calculate Dissimilarity Matrix by CPU');
    
    if isempty(gcp('nocreate'))
        parpool('local',20);
    end
    for J_C1=1:Repeattimes%%%%%%%%%%%%%%
        TV_p=TV_All{ J_C1,1};
        R_p=R_All{ J_C1,1};
        t_p=t_All{ J_C1,1};
        
        
        [~, MatrixAfterAll2all_norm]=Cost_3DReal_cpu(V2,TV_p,R_p,t_p,Particles1,Cscale);
        
        cost_norm_All{J_C1,1}=MatrixAfterAll2all_norm;
    end
    disp(' Dissimilarity Matrix Calculation finished by CPU');
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if ~isempty(poolobj)
        delete(poolobj);
    end
end
%----------------Calculate Cost Function-----------

time_Dissimilarity=toc;
