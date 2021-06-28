%VISUALCLOUD2D Displays the binned density of the pointcloud in 2D.
% 
% SYNOPSIS:
%   visualizeCloud2D(pointcloud, bins, diameter, angle)
% 
% INPUT
% X
%       [x y z; ...] (for example Particle{10}.points)
%n
%       number of bins in all directions
%   diameter
%       the diameter of the ROI in px

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
% Revision: Hamidreza Heydarian, 2017
%Revision:WenxiuWang,2020

function [dip, Z] = visualizeCloud2DW(X, n,diameter,caption,caption2,caption3,colorscale)

if nargin <4
    caption = '';
end
if nargin <5
    caption2 = '';
end
if nargin <6
    caption3 = '';
end
if nargin<7
    colorscale=250;
end

% rotation matrix
% Rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];

%    % apply rotation to particle
%X = (Rot*(X'))';

% center particle around zero
X(:,1) = X(:,1) - mean(X(:,1),1);
X(:,2) = X(:,2) - mean(X(:,2),1);

% extract the ROI
ROIradius = 0.5*diameter;
X = X(find(X(:,1) < ROIradius & X(:,1) > (-ROIradius)),:);
X = X(find(X(:,2) < ROIradius & X(:,2) > (-ROIradius)),:);

% define the grid
xi = linspace(-ROIradius,ROIradius,n);
yi = linspace(-ROIradius,ROIradius,n);

% discretize the particle coordinates
xr = interp1(xi,1:numel(xi),X(:,1),'nearest');
yr = interp1(yi,1:numel(yi),X(:,2),'nearest');

% binning
Z = accumarray([xr yr],1, [n n]);

% display the image with hot colormap
% dip = dipshow(Z','lin');
dip = dipshow(Z','lin');

dipmapping(dip,[0 max(Z(:))],'COLORMAP',hot(colorscale))

%bar = 10*n/diameter/CCD_pixel; %20nm bar

if caption
    text(n/2,0, caption,'Color','black','FontSize',10,'FontWeight' ...
        ,'bold','HorizontalAlignment','center','VerticalAlignment' ...
        ,'bottom');
end

if caption2
    text(n/2,n+30, caption2,'Color','black','FontSize',10,'FontWeight' ...
        ,'bold','HorizontalAlignment','center','VerticalAlignment' ...
        ,'bottom');
end
if caption3
    text(n/2,n+15, caption3,'Color','black','FontSize',10,'FontWeight' ...
        ,'bold','HorizontalAlignment','center','VerticalAlignment' ...
        ,'bottom');
end

%  rectangle('Position',[n-140,n-70,bar,20],'FaceColor','white')

end

