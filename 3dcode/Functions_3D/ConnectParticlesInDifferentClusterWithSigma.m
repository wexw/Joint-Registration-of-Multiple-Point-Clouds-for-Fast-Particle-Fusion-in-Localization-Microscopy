%This is used to connect particles in differnet cluster through the same
%particle in different cluster.
%Input:
%Repeattimes: the progress of JRMPC+Classify is repeated Repeattimes
%clusAll: all good clusters(including Serial number pf particles in each cluster) from Repeattimes JRMPC+Classify
%clus_initial: the cluster with most particles. All the particles will be transferred into this cluster.
%maxp: the clus_initial is comes from the maxp^th JRMPC+Classify
%RAll: all Rotation matrix for R times JRMPC
%tAll: all translation matrix for R times JRMPC
%TVAll: all transformed particles for R times JRMPC
%TVrefinedAll: all refinedtransformed particles for R times JRMPC
%Output:
%TVPick: transfomed particles after connecting
%TVPick: transfomed refined particles after connecting
%RPick: Rotation matrix for transformed particles
%tPick: translation matrix for TVPick
function [TVPick, TVrefinedpick, Rpick,tpick]  = ConnectParticlesInDifferentClusterWithSigma(limit,Repeattimes,clusAll,clus_initial,maxp,RAll,tAll,TVAll,TVrefinedAll)
clus_accumulate=clus_initial;
Rforth=RAll{maxp,1};     %all Rotation matrix from JRMPC of the best cluster
tforth=tAll{maxp,1};     %all translation matrix from JRMPC of the best cluster
TVforth=TVAll{maxp,1};   %all Transformed particles from JRMPC of the best cluster
TVrefinedforth= TVrefinedAll';%all Transformed bindingsites from JRMPC of the best cluster
TVPick=TVforth(clus_initial,1);%Transformed particles of the best cluster
TVrefinedpick=TVrefinedforth(clus_initial,1);%Transformed bindingsites of the best cluster
iterR=0;
ADD=[];
Rpick=cell(100,1);%inicialize final rotation matrix
tpick=cell(100,1);%inicialize final translation matrix
Rpick(clus_initial)=Rforth(clus_initial);
tpick(clus_initial)=tforth(clus_initial);

Secondtransferredparticles=[];
Mfinalsets=[];
%Combine All available clusters together

for i=1:1:Repeattimes
    clusers_Number=size(clusAll{i,1},2);
    Rturn=RAll{i,1}; %all Rotation matrix from ith JRMPC
    tturn=tAll{i,1}; %all translation matrix from ith JRMPC
    TVturn=TVAll{i,1}; %all transformed particles from ith JRMPC
    TVrefinedturn=TVrefinedAll; %all Transformed bindingsitesfrom ith JRMPC
    for j=1:1:clusers_Number
        clus_change=[];
        clus_change=clusAll{i,1}{1,j};           %particle index in cluster j from ith JCC
        [output,common]=Compare2Clusters(clus_initial, clus_change,clus_accumulate,Rturn,tturn,Rforth,tforth,limit); %if common particles exist,output is unique particle which is not included in clus_accumulate.otherwise, output is 0.
        
        if output~=0
            clus_accumulate=[clus_accumulate;output];%update clus_accumulate
            
            for pi=output'                            %Add information of unique particle into final resultes one by one
                Mfinal=TVturn{pi,1};
       
                Mback=Rturn{common,1}'*(Mfinal-tturn{common,1}); %transfer transformed unique particles back to original position of common particle
                %                 TBback=Rturn{common,1}'*(TBfinal-tturn{common,1});%transfer bindingsites of transformed unique particles back to original position of common particle
                Mforth=Rforth{common,1}*Mback+tforth{common,1};  %transfer unique particles to target cluster's position.
                %                 TBforth=Rforth{common,1}*TBback+tforth{common,1};%transfer bindingsites of unique particles to target cluster's position.
                TBforth=TVrefinedAll{1,pi};
                TVPick=[TVPick;Mforth];                          %update transformed particles
                %Secondtransferredparticles= [Secondtransferredparticles,Mforth];
                %Mfinalsets=[Mfinalsets,Mfinal];
                TVrefinedpick=[TVrefinedpick;TBforth];   %update transformed BindingSites
                Rpick{pi,1}=Rforth{common,1}*Rturn{common,1}';
                tpick{pi,1}=tforth{common,1}-Rforth{common,1}*Rturn{common,1}'*tturn{common,1};
            end
            %       [dip, Z] = visualizeCloud2DW(Secondtransferredparticles', 200, 120);
            %[dip, Z] = visualizeCloud2DW(Mfinalsets', 200, 120);TV28_12=TVAll{12,1}{28,1};
            %TVcommon_i=TVAll{i,1}{common,1};
            %TVcommon_maxp=TVAll{maxp,1}{common,1};
            %[dip, Z] = visualizeCloud2DW(TVcommon_i', 200, 120);
            %[dip, Z] = visualizeCloud2DW(TVcommon_maxp', 200, 120);
            
        end
    end
end

end

