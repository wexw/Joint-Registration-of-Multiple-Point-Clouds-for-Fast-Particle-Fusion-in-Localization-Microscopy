function [cost,cost_norm]=Cost_2D_gpu(TV,particles,CCD_pixelsize)
%Calculate the cost and normalization of cost between orignal particle(particles) and
%Transformed particles(TV).
%CCD_pixelsize is used to make the TV the same scale with orginal
%particlee
K_c= size(TV,1); 
cost = zeros(K_c,K_c); 
cost_norm = zeros(K_c,K_c); 
TVOrigin = cellfun(@(TV) TV./CCD_pixelsize,TV,'uniformoutput',false);
SigmaSquare=cellfun(@(particles)  (particles.sigma).^2,particles,'uniformoutput',false);
for i = 1:K_c
    for j = (i+1):K_c%%%%%%%%%%%%%   
        %%%%%%%mex_expdist need Classification code.
            cost(i,j) = mex_expdist(TVOrigin{i}',TVOrigin{j}', SigmaSquare{1,i}, SigmaSquare{1,j});    
            cost_norm(i,j) = cost(i,j)/(size(TVOrigin{j}',1)*size(TVOrigin{i}',1));
    end
end
end