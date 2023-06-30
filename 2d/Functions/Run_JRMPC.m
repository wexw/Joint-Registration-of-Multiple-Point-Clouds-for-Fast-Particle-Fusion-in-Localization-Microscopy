disp('Start to Run JRMPC');
tic
%------------------prepare output format for parallel for------------------%
R_All=cell(Repeattimes,1);              %rotation matrix
t_All=cell(Repeattimes,1);               %translation matrix
TV_All=cell(Repeattimes,1);            %Transformed particles
X_All=cell(Repeattimes,1);              %Final GMM centers
S_All=cell(Repeattimes,1);               %Final GMM width
%------------------prepare output format for parallel for------------------%

for JK=1:Repeattimes%%%%%%%%%%%%%%
    disp(['Start to run JRMPC with initialization ' num2str(JK)]);
    %------------------Generate Xin1-----------------
    az1 = 2*pi*rand(1,K);     %  modelname='tud383particlesdol80';
    el1 = 2*pi*rand(1,K);
    Xin1 = [cos(az1); sin(el1)]*CCD_pixelsize;
    %------------------Generate Xin1-----------------
    iterations=300;% Number of iterations for JRMPC
    disp(['...']);

    %------------------JRMPC 2D-----------------
    [R ,t,X,S,a] = jrmpc_2D_Nooutliers(V1,Xin1,'maxNumIter',iterations, 'epsilon', 1e-5,'gamma',0.01);   %JRMPC
    %[R ,t,X,S,a] = jrmpc_2D_Nooutliers(V1,Xin1,'maxNumIter',iterations, 'epsilon', 1e-5,'S',Sa);   %JRMPC

    TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V1,R,t ,'uniformoutput',false); %Transformed V
    %------------------JRMPC 2D-----------------
    R_All{JK,1}=R;                        %all Rotation matrix from PMKth JRMPC
    t_All{JK,1}=t;                        %all translation matrix from PMKth JRMPC
    TV_All{JK,1}=TV;                      %all transformed particles from PMKth JRMPC
    X_All{JK,1}=Xin1;
    S_All{JK,1}=S;
    %------------------JRMPC 2D-----------------
    disp(['JRMPC with initialization ' num2str(JK) ' is finished']);


end
time_JRMPC=toc
disp('JRMPC finished');
