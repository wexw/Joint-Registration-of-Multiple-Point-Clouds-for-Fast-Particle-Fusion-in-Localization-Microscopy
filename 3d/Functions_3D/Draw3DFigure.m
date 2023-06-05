disp(['Start to Draw 3D figures']);
%Draw figure

%filter out localizations out of the field of view
idx = find(data2c(:,1) <  sizeview& ...
    data2c(:,1) > -sizeview  & ...
    data2c(:,2) <  sizeview & ...
    data2c(:,2) > -sizeview  & ...
    data2c(:,3) <  sizeview & ...
    data2c(:,3) > -sizeview);
data2c = data2c(idx,:);
  set(0,'defaultfigurecolor','black') 
  densityc=visualizeCloud3D2(data2c,10,0);
figure()
hold on
scatter3(data2c(:,1),data2c(:,2), data2c(:,3),1,densityc,'.');
title(sprintf('TV Number of particles=%d',effectiveParticlesNumber), 'fontweight','bold','fontsize',12,'Color', 'w');
colormap(hot);
colordef black
hold off
%         v = VideoWriter(filename);
%         open(v);
%         for beta=0:180
%        
%             renderprojection(data2c(:,1),data2c(:,2), data2c(:,3), 0,beta,[-100 100],[-120 120],1,1, 1,10);
%     
%             drawnow;
%             
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         end
%         
%         close(v);