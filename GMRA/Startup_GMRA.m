function Startup_GMRA(BaseDirectory)

% Startup_GMRA adds the various directories used by the GMRA code to the current path.  
% Startup_GMRA will assume the current directory is the base directory unless another directory is specified
fprintf('Startup_GMRA.m: setting GMRA paths ... \n');

if nargin==0
    Prefix  = [pwd filesep];
else
    Prefix  = [BaseDirectory filesep];
end;

appendpath([Prefix]);

% Look for Diffusion Geometry
if ~exist('GraphDiffusion.m','file')
    if ~exist(['..' filesep 'DiffusionGeometry/GraphDiffusion.m'],'file')
        fprintf('\nDiffusion Geometry is needed for the GMRA package, but it could not be found. Please install it and/or add it to the Matlab path');
    else
        cd('../DiffusionGeometry');
        Startup_DiffusionGeometry
        cd(Prefix);
    end
end

return

function appendpath(string)

fprintf('\t%s\\ \n', string);
addpath(genpath(string));

return;

