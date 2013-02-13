function csq_deps(varargin)
% csq_deps(.....)
% This function handles all of the proper dependencies, making sure that
% the desired directories are all reachable to allow for desired
% functionality.

% Clearing out your path -- no problems!
restoredefaultpath

% Get current directory
cd_str = cd;
repo_name = 'CSQuantization';
repo_dir = cd_str(1:(strfind(cd_str,repo_name)+length(repo_name)-1));

% Check inputs and get strings to add
for i=1:nargin
    switch varargin{i}
        % Common directories
        case 'common'
            addpath([repo_dir '/common']);
            addpath([repo_dir '/common/image']);
            addpath([repo_dir '/common/csq']);
            addpath([repo_dir '/common/metrics']);
        case 'common-image'
            addpath([repo_dir '/common/image']);
        case 'common-csq'
            addpath([repo_dir '/common/csq'])
        case 'common-metrics'
            addpath([repo_dir '/common/metrics']);
            
        % Toolbox directories
        case 'ssim'
            addpath([repo_dir '/toolboxes/ssim']);
        case 'wavelet'
            addpath([repo_dir '/toolboxes/WaveletSoftware']);
        case 'bpdq'
            addpath([repo_dir '/toolboxes/bpdq_toolbox']);
        case 'bcs-spl'
            addpath([repo_dir '/toolboxes/BCS-SPL-1.5-1']);
        case 'mc-cs-pl'
            addpath([repo_dir '/toolboxes/MC-CS-PL']);
        
        % Top level directories
        case 'biht'
            addpath([repo_dir '/biht']);
        
        otherwise
            return_str = sprintf('Unknown dependency: %s',varargin{i});
            warning('csq_deps:Uknown',return_str);
    end
end

