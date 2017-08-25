function [res] = doGMRA(X,d)          
GMRAopts.ManifoldDimension = d;
GMRAopts.smallestMetisNet =10;
% whether to use best approximations
GMRAopts.addTangentialCorrections = true;
GMRAopts.errorType = 'relative';

GMRAopts.Predictor ='none';
GMRAopts.coeffs_threshold = 1e-10;
GMRAopts.MaxDim =d;

% whether to sparsify the scaling functions and wavelet bases
GMRAopts.sparsifying = false;
GMRAopts.sparsifying_method = 'ksvd'; % or 'spams'

% whether to split the wavelet bases into a common intersection and
% children-specific parts
GMRAopts.splitting = false;

% METIS parameters
GMRAopts.knn = 30;
GMRAopts.knnAutotune = 20;
GMRAopts.smallestMetisNet = 10;

% whether to output time
GMRAopts.verbose = 0;

% method for shrinking the wavelet coefficients
GMRAopts.shrinkage = 'hard';

% whether to merge the common part of the wavelet subspaces
% associated to the children into the scaling function of the parent.
GMRAopts.mergePsiCapIntoPhi  = false;

%%%%%%%%%%%%%%%%%%%%%%%% SetGMRAParameters

GMRAopts.GWTversion      = 0;
GMRAopts.PartitionType              = 'covertree';
GMRAopts.CoverTreeExpandOpts        = struct('ExtRange','max');
    if isfield(GMRAopts,'ManifoldDimension') && (GMRAopts.ManifoldDimension>0)
        if GMRAopts.ManifoldDimension==1,
            theta = 0.75;
        else
            theta = 1-1/(2*GMRAopts.ManifoldDimension);
        end
        GMRAopts.CoverTreeBuildOpts      = struct( 'theta'     , theta, ...
            'numlevels' , max([1,int32(round(log(size(X,2)/(10*GMRAopts.ManifoldDimension))/log(1/theta)))]), ...
            'minlevel'  , int32(0), 'NTHREADS'  , int32(feature('numcores')), 'BLOCKSIZE' , int32(2048));        
        GMRAopts.CoverTreeTrimOpts       = struct( 'TrimType','Size','TrimValue',int32(2*GMRAopts.ManifoldDimension));
        %GMRAopts.CoverTreeBuildOpts.numlevels = 100; GMRAopts.CoverTreeTrimOpts.TrimValue = 1;                                 % Go deep: useful for seing bias/variance trade-off
    else
        theta = 0.9;
        GMRAopts.CoverTreeBuildOpts      = struct( 'theta',theta,'numlevels',2*max([1,int32(round(0.5*log(1+size(X,2))/log(1/theta)))]),'minlevel',int32(0), ...
            'NTHREADS',int32(feature('numcores')),'BLOCKSIZE',int32(2048));
        GMRAopts.CoverTreeTrimOpts       = struct( 'TrimType','Size','TrimValue',int32(min([size(X,2),10])));
        %        GMRAopts.CoverTreeBuildOpts.NTHREADS = int32(0);
    end
    GMRAopts.CoverTreeOpts.RefineCoverTree = false;
    GMRAopts.CoverTreeOpts.MaxSamples = Inf;

GMRAopts.parallel                   = false;
GMRAopts.ComputeWavelets            = true;
GMRAopts.ConstructGraph             = false;
GMRAopts.threshold0                 = 0.001;
GMRAopts.threshold1                 = 1e-6;
GMRAopts.threshold2                 = 0.001;
GMRAopts.addTangentialCorrections   = true;
GMRAopts.precision                  = 1e-3;
% Parameters for diffusion maps
GMRAopts.ConstructGraph             = false;                                                                                    % Change this to construct graph and diffusion embedding
GMRAopts.knn                        = 20;
GMRAopts.knnAutotune                = 20;
GMRAopts.graphNormalization         = 'beltrami';
GMRAopts.CoverTreeTrimOpts.TrimValue = 5; 

%X_train = X(:,1:size(X,2)-2);
%X_test = X(:,size(X,2)-2:size(X,2));

gMRA = GMRA(X, GMRAopts);
XGWT = FGWT(gMRA, X); 

[Projections, ~] = IGWT(gMRA, XGWT);

leafNodeIdxs   = XGWT.leafNodeIdxs;  % each point corresponds to which leaf, at J level, k_J leaves (1,...k_J)
allLeafs = gMRA.LeafNodes;  % all the leaf nodes index among all nodes, at all levels, sum k_j nodes total 


phatIdx = allLeafs(leafNodeIdxs(size(X,2)-1)); %which leaf 
phatproj = gMRA.ScalBasis{1,phatIdx} * (X(:,size(X,2)-1)- gMRA.Centers{1,phatIdx});

%phatproj = gMRA.ScalBasis{phatIdx (@minus, X(:,size(X,2)-1), gMRA.Centers{phatIdx})}

p0Idx = allLeafs(leafNodeIdxs(size(X,2))); 
p0proj = gMRA.ScalBasis{1,p0Idx} * (X(:,size(X,2))- gMRA.Centers{1,p0Idx});

%if leafNodeIdxs(size(X,2)-1) > leafNodeIdxs(size(X,2))
 %   for i = leafNodeIdxs(size(X,2)):leafNodeIdxs(size(X,2)-1)
  %  end 
%elseif leafNodeIdxs(size(X,2)-1) < leafNodeIdxs(size(X,2))
    
%else leafNodeIdxs(size(X,2)-1) == leafNodeIdxs(size(X,2))
%end


%res = norm(p0proj-phatproj);
J =size(Projections,3);
N= size(Projections,2);
%res = norm(Projections(:,N,J-3)-Projections(:,N-1,J-3) );
res =[];
for j = 1:J
    res(j) = norm(Projections(:,N,j)-Projections(:,N-1,j) );
end

return