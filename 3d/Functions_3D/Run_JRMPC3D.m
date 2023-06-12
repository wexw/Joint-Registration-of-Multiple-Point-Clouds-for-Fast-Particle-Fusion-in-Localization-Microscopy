%Run_JRMPC_3D
disp('Start to Run JRMPC 3D');

tic
%------------------prepare output format for parallel for------------------%
R_All=cell(Repeattimes,1);              %rotation matrix
t_All=cell(Repeattimes,1);               %translation matrix
TV_All=cell(Repeattimes,1);            %Transformed particles
X_All=cell(Repeattimes,1);
X_inAll=cell(Repeattimes,1);
%------------------prepare output format for parallel for------------------%

%The inputs of the JRMPC_3D are of the unit nm/CCD_pixelsize
%So that some values calculated in the JRMPC would not approach
% the maimum/minimum limit of the MATLAB
%The outputs are of the unit of nm by multiply CCD_pixelsize
for JK=1:Repeattimes
    %estimate number of GMMcenters K for JRMPC
    %------------------Generate Xin1-----------------
    az1 = 2*pi*rand(1,K);
    el1 = 2*pi*rand(1,K);
    Xin1 = [cos(az1).*cos(el1); sin(el1); sin(az1).*cos(el1)];
    [R ,t,X,S,a] = jrmpc_Nooutliers(V1,Xin1,'maxNumIter',500, 'epsilon', 1e-5,'S',Sa);%'gamma',0.2);%'S',1000);   %JRMPC

    TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V1,R,t ,'uniformoutput',false); %Transformed V


    R_All{JK,1}=R;                        %all Rotation ma                        az1 = 2*pi*rand(1,K);
    t_All{JK,1}=cellfun(@(t) t.*Cscale,t ,'uniformoutput',false);
    %all translation matrix from PMKth JRMPC
    TV_All{JK,1}=cellfun(@(TV) TV.*Cscale,TV ,'uniformoutput',false);                    %all transformed particles from PMKth JRMPC
%...-----------------Calculate Cost Function-----------


    X_All{JK,1}= X.*Cscale;
    X_inAll{JK,1}=Xin1.*Cscale;
    S_All{JK,1}=S.*Cscale;

    %------------------JRMPC 2D-----------------
end
disp('JRMPC 3D for finished');
toc
Time_JRMPC=toc;
    