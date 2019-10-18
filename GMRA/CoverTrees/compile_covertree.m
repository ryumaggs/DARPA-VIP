%% Script for compiling covertree code interface to Matlab
% This has been tested on OS X, and a few versions of Linux (Fedora, Ubuntu, CentOS). Please let me know of trouble you encounter
% with compilation on these or other platforms, and the command line if you succeed in compiling on them.

%%
% Unfortunately different versions (with corresponding different MEX files) need to be compiled for different precisions.
% Float is much faster when (most) GPU's are used by covertree in non-Euclidean distance computations.
FloatOrDouble = {'FLOAT','DOUBLE'};                                         % Compile both FLOAT and DOUBLE precision
Debug = false;
EIGEN_DIR = [];
EIGEN_DIRs = {'/Users/mauro/ownCloudMath/Code/Eigen','/home/mauro/ownCloudMath/Code/Eigen'};
for k = 1:length(EIGEN_DIRs)
    if exist(EIGEN_DIRs{k}), EIGEN_DIR = EIGEN_DIRs{k}; break; end
end
if isempty(EIGEN_DIR), warning('EIGEN directory not found among the ones specified.'); end

%% Compile
for k = 1:length(FloatOrDouble),   
    Options = ['-compatibleArrayDims COPTIMFLAGS=''-Ofast -DNDEBUG''  -D_EIGEN_ -DMEX -D' FloatOrDouble{k} ' -I../ -I/usr/local/include -L/usr/local/lib/ -lfftw3 -lfftw3f -I' EIGEN_DIR];   %-I/usr/local/include/eigen3 
    
    os = system_dependent('getos');
    if ~isempty(findstr('Darwin',os))
        Options = [Options ' -DNAMED_SEMAPHORES -D_ACCELERATE_ON_']; % -Dchar16_t=UINT16_T'];
        Options = [Options ' -D__MACOSX_CORE__ LDFLAGS=''\$LDFLAGS -framework Accelerate'''];
    end
    
    Files = ['../ThreadsWithCounter.C ../IDLList.C ../IDLListNode.C ../Vector.C ../Cover.C ../Point.C ../CoverNode.C ../EnlargeData.C ../Timer.C covertree_MEXint.C ' ...
            '../Distances.C ../ImageUtils.cpp ../convolution_fftw.C ../factorize.C ../TimeUtils.cpp ../CoverForest.C ../findNearest.C ../FindNearestData.C ../findWithin.C ../FindWithinData.C'];
    
    %% Compile covertree
    fprintf('\n Compiling covertree...');
    if strcmpi(FloatOrDouble{k},'DOUBLE')
        outputfilenameaddition = 'D';
    else
        outputfilenameaddition = '';
    end
    if Debug
        Options = ['-g ' Options];
    end
    eval(['mex ' Options ' covertree.C ' Files ' -lpthread -output covertree' outputfilenameaddition]);
    
    fprintf('done.');
    
    %% Compile findnearest
    fprintf('\n Compiling findnearest...');
    eval(['mex ' Options ' findnearestMEX.C ' Files ' ../FindWithinData.C ../findNearest.C -lpthread -output findnearestMEX' outputfilenameaddition']);
    fprintf('done.');  

    %% Compile findwithin
    fprintf('\n Compiling findwithin...');
    eval(['mex ' Options ' findwithinMEX.C ' Files ' ../FindWithinData.C ../findWithin.C -lpthread -output findwithinMEX' outputfilenameaddition]);
    fprintf('done.');
    
end

% Remove object files resulting from compilation
system('rm -rf *.o');

fprintf('\n done.\n');

