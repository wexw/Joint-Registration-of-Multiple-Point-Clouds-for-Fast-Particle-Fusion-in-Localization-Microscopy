   %...-----------------Calculate Cost Function-----------
   disp('Start to Calculate Dissimilarity Matrix 3D');
   tic
%V2's unit is nm, V2 is the original input particles' localizations
   cost_norm_All=cell(Repeattimes,1);
   for J_C1=1:Repeattimes%%%%%%%%%%%%%%
       TV_p=TV_All{ J_C1,1};
       R_p=R_All{ J_C1,1};
       t_p=t_All{ J_C1,1};
       [~, MatrixAfterAll2all_norm]=Cost_3DReal(V2,TV_p,R_p,t_p,Particles1,Cscale);
       cost_norm_All{ J_C1,1}=MatrixAfterAll2all_norm;
   end
   %----------------Calculate Cost Function-----------

   time_Dissimilarity=toc;
   disp(' Dissimilarity Matrix Calculation finished');