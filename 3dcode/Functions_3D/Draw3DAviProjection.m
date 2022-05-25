  %filter out localizations out of the field of view
  filenameTVpick2='Demo.avi'
  data2c=cell2mat(TVPick')';
% data2c=BottomRing';
  sizeview=200;
idx = find(data2c(:,1) <  sizeview& ...
    data2c(:,1) > -sizeview  & ...
    data2c(:,2) <  sizeview & ...
    data2c(:,2) > -sizeview  & ...
    data2c(:,3) <  sizeview & ...
    data2c(:,3) > -sizeview);
data2c = data2c(idx,:);
% 

  set(0,'defaultfigurecolor','black') 

        v = VideoWriter(filenameTVpick2);
        open(v);
        for beta=0:1:180
       
            renderprojection(data2c(:,1),data2c(:,2), data2c(:,3), 0,beta,[-150 150],[-150 150],1,1, 1,10);
    
            drawnow;
%             set(gcf,'color','black')
            axis off

            
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
        close(v);