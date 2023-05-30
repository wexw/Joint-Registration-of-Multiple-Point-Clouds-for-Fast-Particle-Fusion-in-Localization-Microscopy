function [cost,cost_norm]=Cost_3DReal_cpu(V,TV,R,t,particles,CCD_pixels)
    if nargin <6
       CCD_pixels = 130;
    end
K_c= size(TV,1); 

cost = zeros(K_c,K_c); 
cost_norm = zeros(K_c,K_c); 
SigmaSquare1=cellfun(@(particles)  [particles.sigma(:,1).^2.*CCD_pixels^2 particles.sigma(:,2).^2.*CCD_pixels^2],particles,'uniformoutput',false);
for i = 1:K_c
    for j = (i+1):K_c
            S.points = V{i}';           %take S in native orientation
            M.points = V{j}'; 
            
            Mfinal = TV{j};             %final orientation of M
            Mback = R{i}'*(Mfinal-t{i});%inverse transform M with the transformation of S
            [~,~,procr] = procrustes(Mback',V{j}','scaling',false,'reflection',false);   %find transformation how to transform M(native) to inverse transformed M
            M.points = M.points*procr.T + procr.c;  %do transformation
            S.sigma =SigmaSquare1{1,i};
            M.sigma =SigmaSquare1{1,j};
            
            cost(i,j) = mex_expdist_cpu(S.points, M.points, correct_uncer(S.sigma), M.sigma, procr.T);    %for now: rot-matrix set to identity matrix (IS WRONG)
            cost_norm(i,j) = cost(i,j)/(size(M.points,1)*size(S.points,1));
    end
end
