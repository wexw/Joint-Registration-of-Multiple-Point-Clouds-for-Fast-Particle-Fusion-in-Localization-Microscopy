function [cost,cost_norm]=Cost_2D_cpu(TV,particles,CCD_pixelsize)
%Calculate the cost and normalization of cost between orignal particle(particles) and
%Transformed particles(TV).
%CCD_pixelsize is used to make the TV the same scale with orginal
%particle
K_c=size(TV,1);
TVOrigin = cellfun(@(TV) TV./CCD_pixelsize,TV,'uniformoutput',false);
SigmaSquare=cellfun(@(particles)  (particles.sigma).^2,particles,'uniformoutput',false);
% Pre-compute size values
size_values = cellfun(@(x) size(x',1), TVOrigin);
% Generate all pairs of indices for the upper triangular part
[ind_i, ind_j] = find(triu(ones(K_c),1));
% Initialize vectors for cost and cost_norm values
cost_values = zeros(size(ind_i));
norm_values = zeros(size(ind_i));
TVOrigin_i = TVOrigin(ind_i);
TVOrigin_j = TVOrigin(ind_j);
SigmaSquare_i = SigmaSquare(1,ind_i);
SigmaSquare_j = SigmaSquare(1,ind_j);
parfor k = 1:length(ind_i)
    cost_values(k,1) = mex_expdist_cpu(TVOrigin_i{k}',TVOrigin_j{k}', SigmaSquare_i{k}, SigmaSquare_j{k});
end
cost = full(sparse(ind_i, ind_j, cost_values, K_c, K_c));
norm_values = cost_values ./ (size_values(ind_j) .* size_values(ind_i));
cost_norm =full(sparse(ind_i, ind_j, norm_values, K_c, K_c));
end
