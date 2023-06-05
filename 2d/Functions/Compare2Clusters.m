function [output,common]=Compare2Clusters(clus_fix, clus_change,clus_accumulate,Rturn,tturn,Rforth,tforth)
clus_common=intersect( clus_fix, clus_change);
scc=size(clus_common);%number of common particles
if scc>0
    for pcp=1:1:size(clus_common,1)
        pc=clus_common(pcp,1);
        Rcommon{pcp,1}=Rforth{pc,1}*Rturn{pc,1}';
        tcommon{pcp,1}=[tforth{pc,1}-Rforth{pc,1}*Rturn{pc,1}'*tturn{pc,1}]';
    end
    Rt = cellfun(@(Rcommon,tcommon) cat(1,Rcommon,tcommon),Rcommon,tcommon ,'uniformoutput',false); %Transformed V
    sizeColoum=numel(Rt{1,1});
    RowRt = cellfun(@(Rt) reshape(Rt,1,sizeColoum),Rt ,'uniformoutput',false); %Transformed V
    
    RowRt_matrix=cell2mat(RowRt);
    RowRt_median=median(RowRt_matrix);
    RowRt_minus=RowRt_matrix-RowRt_median;
    RowRt_sum=sum(abs(RowRt_minus),2);
    
    [minvalue,min_position]=min(RowRt_sum);
    clus_delete=intersect(clus_accumulate,clus_change);
    
    for i=1:1:size(clus_delete)
        deleteposition=find(clus_change==clus_delete(i));
        clus_change(deleteposition)=[];
    end
    
    if size(clus_change,1)>0
        output=clus_change;
        common=clus_common(min_position);
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

           
 