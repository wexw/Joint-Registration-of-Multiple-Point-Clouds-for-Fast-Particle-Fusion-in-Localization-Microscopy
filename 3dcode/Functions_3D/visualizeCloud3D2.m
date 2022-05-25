%VISUALCLOUD2D 3D scatter plot of the points with an intensity proportional
%to its local density
% 
% SYNOPSIS:
%   visualizeCloud2D(pointcloud, bins, diameter, angle)
% 
% INPUT
%   pointcloud
%       [x y z; ...] (for example Particle{10}.points)
%   scale
%       the scale parameters of Gauss transform
% 
% OUTPUT
%   dip
%       dip image structure
%   Z
%       image matrix
%
% DEFAULTS:
%   none
% 
% NOTES:
% 
% (C) Copyright 2015               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Diederik Feilzer, July 2015
% Revision: Hamidreza Heydarian, 2018

function [density] = visualizeCloud3D2(data, scale, USE_GPU)

    density = zeros(size(data,1),1);   
    
%     idxTop = find(data(:,3) > 0);
%     for i=1:size(idxTop,1)
% 
%         [density(idxTop(i)),~] = GaussTransform(data(idxTop,:),data(idxTop(i),:), scale);
% 
%     end        
%     
%     idxBottom = find(data(:,3) <= 0);
%     for i=1:size(idxBottom,1)
% 
%         [density(idxBottom(i)),~] = GaussTransform(data(idxBottom,:),data(idxBottom(i),:), scale);
% 
%     end        
    
    % compute the intensity of each point using the Gauss transform
    for i=1:size(data,1)

        [density(i),~] = GaussTransform(data,data(i,:), scale, USE_GPU);

    end    
    
    % uncomment for new figure
%     figure('pos',[100 100 500 500]);

%     cdensity = buildcmap(density);
%     figure,showset(data,cdensity);    

end