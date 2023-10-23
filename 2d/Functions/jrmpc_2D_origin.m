function [R,t,X,Q,a,pk,T] = jrmpc_2D_origin(V,X,varargin)

sqe = @(Y,X) sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);

% ======================================================================== %
%                           C H E C K S                                    %
% ======================================================================== %

V = V(:);
M = numel(V);

[dim,K] = size(X);

if dim ~= 2
    error('X must be a 2 x K matrix.');
end

%if mod(numel(varargin),2)
%    error('odd number of optional parameters, Opt parames must be given as string-value pairs.')
%end

for j=1:M
    if size(V{j},1) ~= 2
        error('V must be an M x 1 cell of 2 x .. matrices, V{%d} has %d in dimension 1.',j,size(V{j},1));
    end
end

% ======================================================================== %
%                       V A R A R G I N  P A R S E R                       %
% ======================================================================== %

isSetR = 0;
isSetT = 0;
isSetQ = 0;
isSetMaxNumIter = 0;
isSetInitialPriorsOrGamma = 0;
isSetEpsilon = 0;
isSetUpdatePriors = 0;

for i=1:2:numel(varargin)
    if ~isSetR && strcmpi(varargin{i},'r')
        
        R = varargin{i+1};
        R = R(:);
        isSetR = 1;
        
    elseif ~isSetT && strcmpi(varargin{i},'t')
        
        t = varargin{i+1};
        t = t(:);
        isSetT = 1;
        
    elseif ~isSetQ && strcmpi(varargin{i},'s')
        
        if isscalar(varargin{i+1})
            
            Q = repmat(varargin{i+1},K,1);
        else
            
            Q = varargin{i+1};
        end
                
        isSetQ = 1;
        Q = 1./Q;
        
    elseif ~isSetMaxNumIter && strcmpi(varargin{i},'maxnumiter')
        
        maxNumIter = varargin{i+1};
                
        isSetMaxNumIter = 1;
        
    elseif ~isSetEpsilon && strcmpi(varargin{i},'epsilon')
        
        epsilon = varargin{i+1};
                
        isSetEpsilon = 1;
        
    elseif ~isSetUpdatePriors && strcmpi(varargin{i},'updatepriors')
        
        updatePriors = varargin{i+1}; % don use the flag,it will affect subsequent parse
                
        isSetUpdatePriors = 1;
        
    elseif ~isSetInitialPriorsOrGamma && strcmpi(varargin{i},'gamma')
        
        gamma = varargin{i+1};
                
        pk = repmat(1/(K*(gamma+1)),K,1);
        
        isSetInitialPriorsOrGamma = 1;
        
    elseif ~isSetInitialPriorsOrGamma && strcmpi(varargin{i},'initialpriors')
        
        if isscalar(varargin{i+1})
            
            pk = repmat(varargin{i+1},K,1);
        else
            
            pk = varargin{i+1};
            
        end
        
        gamma = (1-sum(pk))/sum(pk);
        
        isSetInitialPriorsOrGamma = 1;
        
    else
        
        if isSetInitialPriorsOrGamma
            
            error('Only one of the parameters ''initialPriors'' and ''gamma'' must be set.');
        else
            error('uknown option %s, or already set.',varargin{i});
        end
    end
end

% ======================================================================== %
%                   I N I T I A L I Z E   D E F A U L T S                  %
% ======================================================================== %

if ~isSetR
    
    R = repmat({eye(2)},M,1);
    
end

if ~isSetT
    
    t = cellfun(@(V) (-mean(V,2)+mean(X,2)),V,'uniformoutput',false);
    
end

% transformed sets based on inititial R & t (\phi(v) in the paper)
TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);

if ~isSetQ
    
    [minXyZ,maxXyZ] = cellfun(@(x) deal(min(x,[],2),max(x,[],2)),[TV;X],'uniformoutput',false);
    
    minXyZ = min(cat(2,minXyZ{:}),[],2);
    
    maxXyZ = max(cat(2,maxXyZ{:}),[],2);
    
    Q = repmat(1./(sqe(minXyZ,maxXyZ)),K,1);

end

if ~isSetMaxNumIter
    
    maxNumIter = 100;
    
end

if ~isSetEpsilon
    
    epsilon = 1e-6;
    
end

if ~isSetUpdatePriors
    
    updatePriors = 0;
    
end

if ~isSetInitialPriorsOrGamma
    
    gamma = 1/K;
    
    pk = repmat(1/(K+1),K,1);
    
end

% ======================================================================== %
%                                   E  M                                   %
% ======================================================================== %

% if requested, allocate an empty cell in maxNumIter dimension
if nargout > 6
    T = cell(M,2,maxNumIter);
end

% parameter h in the paper (this should be proportional to the volume that
% encompasses all the point sets). Above, we initially translate the sets around
% (0,0,0), and we compute accordingly the initial variances (and precisions)
%. Thus, we compute h in a similar way.
h = 2/mean(Q); 

beta = gamma/(h*(gamma+1));

%keyboard

pk = pk'; % used as a row


for iter = 1:maxNumIter
    
    % POSTERIORS
    
    % sqe (squared differences between TV & X)
    a = cellfun(@(TV) sqe(TV,X),TV,'uniformoutput',false);
    
    % pk*S^-1*exp(-.5/S^2*||.||)
    a = cellfun(@(a) bsxfun(@times,pk.*(Q'.^1),exp(bsxfun(@times,-.5*Q',a))),a,'uniformoutput',false);
    
    % normalize
    a = cellfun(@(a) bsxfun(@rdivide,a,sum(a,2)+beta),a,'uniformoutput',false);    
   
   
    % ------  weighted UMEYAMA ------ 
    
    lambda = cellfun(@(a) sum(a)',a,'uniformoutput',false); % 1 x K rows 
    
    W = cellfun(@(V,a) bsxfun(@times,V*a,Q'),V,a,'uniformoutput',false);
    
    % weights, b
    b = cellfun(@(lambda) lambda.*Q,lambda,'uniformoutput',false);
    
    % mean of W
    mW = cellfun(@(W) sum(W,2),W,'uniformoutput',false);
    
    % mean of X
    mX = cellfun(@(b) X*b,b,'uniformoutput',false);
    
    % sumOfWeights
    sumOfWeights = cellfun(@ (lambda) dot(lambda,Q),lambda,'uniformoutput',false);
    
    % P
    P = cellfun(@(W,sumOfWeights,mW,mX) X*W' - mX*mW'/sumOfWeights, W,sumOfWeights,mW,mX,'uniformoutput',false);
    
    
    % SVD
    [uu,~,vv] = cellfun(@svd,P,'uniformoutput',false);

    % compute R and check reflection
    R = cellfun(@(uu,vv) uu*diag([1 det(uu*vv)])*vv',uu,vv,'uniformoutput',false);
    
    % solve for t
    t = cellfun(@(mW,mX,R,sumOfWeights) (mX-R*mW)/sumOfWeights,mW,mX,R,sumOfWeights,'uniformoutput',false);
    
    
    % populate T
    if nargout > 6
        T(:,1,iter) = R;
        
        T(:,2,iter) = t;
    end
    
    
    % transformed sets
    TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);
    
    
    % UPDATE X

    den = sum(cell2mat(lambda'),2)'; % den is used for S's update as well

    X = cellfun(@(TV,a) TV*a,TV,a,'uniformoutput',false);%
    
    X = sum(cat(3,X{:}),3);%cat: connect N cells which contain 2*K matrix to 2*K*N matrix. sum(cat(3,X{:}),3): add the 3rd line of 2*K*N matrix
    
    X = bsxfun(@rdivide,X,den);
    
    
    % UPDATE S
    
    % denominators for each j
    wnormes = cellfun(@(TV,a) sum(a.*sqe(TV,X)), TV,a, 'uniformoutput',false);
    
    
    Q = transpose(3*den ./ (sum(cat(3,wnormes{:}),3) + 3*den*epsilon));
    
    
    % UPDATE pk
    
    if updatePriors
        
        pk = den / ((gamma+1)*sum(den));
        
    end

end


% return variances
if nargout > 3
    
    Q = 1./Q;
    

end


