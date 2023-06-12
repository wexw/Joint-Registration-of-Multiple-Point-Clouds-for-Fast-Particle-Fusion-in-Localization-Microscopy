function img=renderprojection2D(x,y,rangex,rangey,pixelsize,sigma, ploton)
%cmax: maximum intensity in rendering: important to keep the intensity
%scaling constant

rx=rangex(1):pixelsize:rangex(end);
ry=rangey(1):pixelsize:rangey(end);
img=histcounts2(y,x,rx,ry);
cmax=max(max(img));
cmax
if sigma>0 %filter
    h=fspecial('gaussian',ceil(2.5*sigma),sigma);
    img=filter2(h,img);
end

if ploton
   nexttile
    imagesc(img,[0, cmax]);
    colormap hot
    axis equal
end
end
