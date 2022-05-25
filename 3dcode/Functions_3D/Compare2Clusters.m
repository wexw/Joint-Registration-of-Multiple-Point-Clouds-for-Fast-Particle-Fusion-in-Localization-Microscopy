function [output,common]=Compare2Clusters(clus_fix, clus_change,clus_accumulate,Rturn,tturn,Rforth,tforth,limit)
clus_common=intersect( clus_fix, clus_change);
scc=size(clus_common,1);
if scc>limit
    clus_delete=intersect(clus_accumulate,clus_change);
    for pcp=1:1:size(clus_common,1)
        pc=clus_common(pcp,1);
        Rcommon{pcp,1}=Rforth{pc,1}*Rturn{pc,1}';
        tcommon{pcp,1}=[tforth{pc,1}-Rforth{pc,1}*Rturn{pc,1}'*tturn{pc,1}]';
    end
    
    for i=1:1:size(clus_common,1)
        for j=1:1:size(clus_common,1)
            corrt(i,j)=corr2(tcommon{i,1},tcommon{j,1});
            corrR(i,j)=corr2(Rcommon{i,1},Rcommon{j,1});
        end
        Sumt(i)=sum(corrt(i,:));
        SumR(i)=sum(corrR(i,:));
        SumTrans(i)=     Sumt(i)+SumR(i);
    end
      [selectvalue,select_position]=max(SumTrans);
    %{
    tcommon_matrix=cell2mat(tcommon);
    tcommon_median=median(tcommon_matrix);
    tcommon_minus=tcommon_matrix-tcommon_median;
    
    if size(tcommon_matrix,2)==3
        tcommon_sum=abs(tcommon_minus(:,1))+abs(tcommon_minus(:,2))+abs(tcommon_minus(:,3));
    else
        tcommon_sum=abs(tcommon_minus(:,1))+abs(tcommon_minus(:,2));
    end

    [selectvalue,select_position]=min(tcommon_sum);

    %}
    for i=1:1:size(clus_delete)
        deleteposition=find(clus_change==clus_delete(i));
        clus_change(deleteposition)=[];
    end
    
    if size(clus_change,1)>0
        output=clus_change;
        common=clus_common(select_position);
        %   [r]=find(tcommon_matrix(:,1)==tcommon_median(1,1))
        %  common=clus_common(max_commonposition);
        %common=clus_common(1);
    else
        output=0;
        common=0;
    end
else
    output=0;
    common=0;
end
end

           
 