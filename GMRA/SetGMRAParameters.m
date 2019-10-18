GMRAopts.GWTversion      = 0;
if Params.PartitionMethods == 0
    GMRAopts.PartitionType      = 'nesdis';
    GMRAopts.smallestMetisNet   = 5;
else
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
end

GMRAopts.parallel                   = false;
GMRAopts.ComputeWavelets            = true;
GMRAopts.ConstructGraph             = false;
GMRAopts.threshold0                 = 0.5;
GMRAopts.threshold1                 = 1e-6;
GMRAopts.addTangentialCorrections   = true;
GMRAopts.precision                  = 0.1;
% Parameters for diffusion maps
GMRAopts.ConstructGraph             = false;                                                                                    % Change this to construct graph and diffusion embedding
GMRAopts.knn                        = 50;
GMRAopts.knnAutotune                = 20;
GMRAopts.graphNormalization         = 'beltrami';



