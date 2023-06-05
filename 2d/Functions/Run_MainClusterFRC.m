tic
TVforth_Main=TV_All{maxR,1};   %all Transformed particles from JRMPC of the best cluster
TVOrigin_Main=TVforth_Main(clus_initial,1);%Transformed particles of the best cluster
sizePart(1)=size(TVOrigin_Main,1);
sizePart(3)=size(TVPick,1);
sizePart(2)=sizePart(3)-sizePart(1);
TVCompare_Main=TVPick(sizePart(1)+1:1:sizePart(3),1);
TV_Main{1,1}=cell2mat(TVOrigin_Main');
TV_Main{2,1}=cell2mat(TVCompare_Main');
pixelsize=1;
zoomfactor=1;
SR_pixelsize=pixelsize/zoomfactor;
for part=1:2
    pixelnum=200;
    diameter=90;
    X_Main=[];
    X_Main=TV_Main{part,1}';
    % Visualize aligned super-particles seperately
    titlename= [ num2str(sizePart(part)) 'particles'];
    [dip, Z] = visualizeCloud2DW(X_Main, pixelnum,diameter,titlename,modelname);
    
    ROIradius = 0.5*diameter;
    X_Main = X_Main(find(X_Main(:,1) < ROIradius & X_Main(:,1) > (-ROIradius)),:);
    X_Main = X_Main(find(X_Main(:,2) < ROIradius & X_Main(:,2) > (-ROIradius)),:);
    X_Main=X_Main+ROIradius;  %X unit [nm]
    X_Main=X_Main./SR_pixelsize; %X unit [pixel]
    InImages_Main{part,1}=X_Main;
end
toc
szx=diameter*zoomfactor;%=diameter*superzoom
szy=diameter*zoomfactor;

%   zoomfactor in binlocalizations
%   Pixels per unit distance in the positions list
%   The input has the unit of pixel
in1_Main = binlocalizations( InImages_Main{1,1},szx,szy);
in2_Main = binlocalizations( InImages_Main{2,1},szx,szy);

[fire_value_Main, frc_curve_Main, fire_high_Main, fire_low_Main] = frcresolution(in1_Main,in2_Main);
time_MainFRC=toc;
fprintf('FRC value %2.2f +- %2.2f nm.\n', fire_value_Main*SR_pixelsize, (fire_low_Main-fire_high_Main)/2*SR_pixelsize);
fprintf('FRC value %2.2f +- %2.2f pixels.\n', fire_value_Main, (fire_low_Main-fire_high_Main)/2);