% correct_uncer
% Makes sure that the uncertainties have the right format to be used by
% mex_expdist. Note that this function must only be used to reformat the 
% sigmas of the first particle (which is not transformed) and when using 
% the GPU function (mex_expdist) and 
%
%   SYNOPSIS:
%       [sigmas_corrected] =correct_uncer(sigmas_in)
%
%   Input: 
%      sigmas_in
%           a Nx2 array containing all sigma_xy and sigma_z per
%           localization ordered horizontally: 
%           sigmas_in = [Sig1_xy Sig1_z ; 
%                                 Sig2_xy Sig2_z ; 
%                                                     ... ]
%
%   Output:
%       sigmas_corrected
%           corrected format for the sigmas, necessary to use mex_expdist.
%           Now the uncertainties per particle are ordered vertically,
%           starting in the upperleft corner, filling the matrix
%           column-wise
%           sigmas_corrected =[Sig1_xy Sig50_xy ; 
%                                            Sig1_z   Sig50_z ; 
%                                            Sig2_xy Sig51_xy;
%                                                                       ... ]
%
%
% (C) Copyright 2019                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Teun Huijben , Feb 2020.

function [sigmas2] = correct_uncer(sigmas1)

    sigmas2 = reshape(sigmas1',size(sigmas1,1),2); 
    
end




