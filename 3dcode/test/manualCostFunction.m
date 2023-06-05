function [costValue] = manualCostFunction(S,M,RM)
%Function to calculate the (sigma) normalized costfunction manually in
%Matlab
% SYNOPSIS:
%       [costvalue] = manualCostFunction(S,M,RM)
%
%   Input: 
%      S,M
%           Matlab struct with 2 fields. The 'points' field is a Nx3
%           matrix with the 3D coordinates of N Localizations. The 'sigma'
%           field is a Nx2 matrix with for each of the N localizations the
%           sigmaXY and sigmaZ localization uncertainty. (particles S and M
%           can have different localizations, but have the same structure) 
%      RM
%           A 3x3 rotation matrix used to transform the uncertainties of
%           particle M. NOTE: the localizations of M already need to be
%           rotated.
%
%   Output:
%      costValue
%           The calculated costFunction value between S and M
%
% (C) Copyright 2018-2020      
% Faculty of Applied Sciences
% Delft University of Technology
%
% Hamidreza Heydarian, November 2020.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%    http://www.apache.org/licenses/LICENSE-2.0

    cross_term = 0; 
    for i = 1:size(S.points,1)                                  %for each loc in S
        for j = 1:size(M.points,1)                              %for each loc in M
            r =(S.points(i,:)'-M.points(j,:)');                 %distance vector     
            Sa = diag([S.sigma(i,1) S.sigma(i,1) S.sigma(i,2)]);%3x3 uncertainty matrix loc i of particle S

            if size(M.sigma,2)==2
                Sb = diag([M.sigma(j,1) M.sigma(j,1) M.sigma(j,2)]);%3x3 uncertainty matrix loc j of particle M
                Sb = RM*Sb*RM';                                 %rotate the uncertainties of M
            elseif size(M.sigma,2)==9
                Sb = reshape(M.sigma(j,:),3,3); 
            else
                disp('error: uncertainties of M have wrong shape')
            end
            
            norm_factor = sqrt(det(Sa+Sb));                     %normalizations factor for a multivariate gaussian
            exponent = r'*inv(Sa+Sb)*r;                         %r'*((Sa+Sb)\r) can be used faster for big matrices
            cross_term = cross_term + exp(-exponent)/norm_factor;
        end
    end
    costValue = cross_term;

end

