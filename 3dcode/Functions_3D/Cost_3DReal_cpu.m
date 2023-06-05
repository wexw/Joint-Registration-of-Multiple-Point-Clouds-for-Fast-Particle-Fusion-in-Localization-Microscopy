function [cost,cost_norm]=Cost_3DReal_cpu(V,TV,R,t,particles,CCD_pixels)
    if nargin <6
       CCD_pixels = 130;
    end
    K_c= size(TV,1);
    SigmaSquare1=cellfun(@(particles)  [particles.sigma(:,1).^2.*CCD_pixels^2 particles.sigma(:,2).^2.*CCD_pixels^2],particles,'uniformoutput',false);
    % Pre-compute size values
    size_values = cellfun(@(x) size(x',1), V);

    % Generate all pairs of indices for the upper triangular part
    [ind_i, ind_j] = find(triu(ones(K_c),1));

    % Initialize vectors for cost and cost_norm values
    cost_values = zeros(size(ind_i));
    norm_values = zeros(size(ind_i));

    % Parallel computation of cost and cost_norm
 
% Parallel computation of cost and cost_norm
for k = 1:length(ind_i)
    i = ind_i(k);
    j = ind_j(k);

    Mfinal = TV{j};             %final orientation of M
    Mback = R{i}'*(Mfinal-t{i});%inverse transform M with the transformation of S
    [~,~,procr] = procrustes(Mback',V{j}','scaling',false,'reflection',false);   %find transformation how to transform M(native) to inverse transformed M
    
    cost_values(k) = mex_expdist_cpu(V{i}', (V{j}'*procr.T + procr.c), correct_uncer(SigmaSquare1{1,i}), SigmaSquare1{1,j}, procr.T); 
    norm_values(k) = cost_values(k) / (size_values(j) * size_values(i));
end
    cost = full(sparse(ind_i, ind_j, cost_values, K_c, K_c));
    cost_norm = full(sparse(ind_i, ind_j, norm_values, K_c, K_c));
end
