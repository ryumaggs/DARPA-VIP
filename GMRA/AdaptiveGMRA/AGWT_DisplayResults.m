function AGWT_DisplayResults( gMRA, XGWT, DataError,displayopts )
%
% function to display results, including
%   1. Error versus scale of GMRA and center approximation
%   2. Error versus partition size of GMRA and adaptive GMRA
%   3. Approximation variation |P_j x - P_{j+1} x| at each point x
%   4. Magnitude of wavelet coefficients
%   5. Dimension of wavelet subspaces
%   6. Wavelet coefficients versus scale
% 
%   gMRA        : GMRA data structure
%   XGWT        : FGWT of data
%   DataError   : structure containing errors of uniform and adaptive GMRA
%   [displayopts]    : data structure of options:
%       [errortype]  : 0 absolute error  1 relative error
%       [errornorm]  : '2' L^2 error 'inf' L^\infty error
%
% (c) Wenjing Liao, Mauro Maggioni


if nargin < 4
    displayopts = struct('errortype',0,'errornorm',2);                     % default: absolute error; L2 norm
end


%% Display L^2 error versus scales/partition size for uniform and adaptive GMRA
fprintf('\n Displaying MSE versus scales/partition size for uniform and adaptive GMRA');
tic
AGWT_DisplayApproxErr(gMRA,DataError,displayopts);
fprintf('done. (%f secs)',toc);

%% Display the approximation variations \|P_j x - P_{j+1} x\| at each point
fprintf('\n Displaying the approximation variations |P_j x - P_{j+1} x| at each point');
AGWT_DisplayDelta(XGWT,DataError.UniformScales)

% Wavelet coefficients
GWT_DisplayCoeffs(gMRA,XGWT);
fprintf('\n\n  ');

if exist('BaillInCell','var')
    if isempty(gMRA.opts.ManifoldDimension)
        AGWT_DisplayBallRadii(gMRA,BallInCell)
    else
        AGWT_DisplayBallRadii(gMRA,BallInCell,gMRA.opts.ManifoldDimension)
    end
end


%% Display linear sets
%AGWT_DisplayLinearSets(gMRA,X_test,XGWT_test)


% %% Save file fro web UI
% % image parameters metadata for visualization
% if ~isempty(fieldnames(imgOpts))
%     imgOpts.imageData = true;
%     imgOpts.Labels = int32(reshape(imgOpts.Labels,length(imgOpts.Labels),1));
%     imgOpts.LabelNames{1} = 'img_cat';
%     imgOpts.LabelDescriptions{1} = 'Image Label';
% else
%     imgOpts.imageData = false;
% end
% 
% imgOpts.LabelDataTypes{1} = 'int32';
% imgOpts.LabelVariableTypes{1} = 'categorical';
% imgOpts.LabelNames{1} = 'no label';
% imgOpts.LabelDescriptions{1} = 'no label';
% 
% if ~isfield(imgOpts,'Labels') || isempty(imgOpts.Labels)
%     imgOpts.Labels = ones(size(X,2),1);
% end
% 
% return
% 
% %matlab_to_hdf5_write(gMRA, imgOpts,[pVisDir pExampleNames{pExampleIdx} '.hdf5'] );
% pause;
% %%
% figure;
% for j= 1:size(Projections,3)
%     plot3(squeeze(Projections(1,:,j)),squeeze(Projections(2,:,j)),squeeze(Projections(3,:,j)),'.');
%     pause;
% end;
% 
% fprintf('\ndone.\n');



return

