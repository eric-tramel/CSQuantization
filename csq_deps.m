function csq_deps(varargin)
% csq_deps(.....)
% This function handles all of the proper dependencies, making sure that
% the desired directories are all reachable to allow for desired
% functionality.

% E.W.T. Removing default path restoration. Takes forever and it seems like
% a bad idea overall.
% % Clearing out your path -- no problems!
% restoredefaultpath

% Get current directory
cd_str = cd;
repo_name = 'CSQuantization';
repo_dir = cd_str(1:(strfind(cd_str,repo_name)+length(repo_name)-1));
addpath(repo_dir);

% Check inputs and get strings to add
for i=1:nargin
    switch varargin{i}
        % Common directories
        case 'common'
            addpath([repo_dir '/common']);
            csq_deps('common-image','common-csq','common-metrics','common-threshold','git');
        case 'common-image'
            if (exist('im2col') + exist('col2im')) == 0
                addpath([repo_dir '/common/image/conflicts']);
            end
            addpath([repo_dir '/common/image']);
        case 'common-csq'
            addpath([repo_dir '/common/csq'])
        case 'common-metrics'
            addpath([repo_dir '/common/metrics']);
        case 'common-threshold'
            addpath([repo_dir '/common/threshold']);
          
        % Toolbox directories
        case 'ssim'
            addpath([repo_dir '/toolboxes/ssim']);
        case 'wavelet'
            % Check to see if we have some issues with the expand function, mostly to
            % do with using octave.
            if csq_in_octave()
                addpath([repo_dir '/common/conflicts/expand']);
            end

            addpath([repo_dir '/toolboxes/WaveletSoftware']);
        case 'bpdq'
            addpath([repo_dir '/toolboxes/bpdq_toolbox']);
        case 'bcs-spl'
            addpath([repo_dir '/toolboxes/BCS-SPL-1.5-1']);
        case 'mc-cs-pl'
            addpath([repo_dir '/toolboxes/MC-CS-PL']);
        case 'bcs-spl-dpcm'
            addpath([repo_dir '/toolboxes/BCS-SPL-DPCM-1.0-2']);
        case 'srm'
            addpath([repo_dir '/toolboxes/SRM/Measurements']);
        case 'qamp'
            addpath([repo_dir '/toolboxes/MessagePassingAlgorithms/QuantizedApproximateMessagePassing']);
            addpath([repo_dir '/toolboxes/MessagePassingAlgorithms/QuantizedRelaxedBeliefPropagation']);
        case 'biht-aop'
            addpath([repo_dir '/toolboxes/AOP1BCSv1']);
        case 'inpaint'
            addpath([repo_dir '/toolboxes/Inpaint_nans']);
        case 'git'
            addpath([repo_dir '/toolboxes/git']);
            
        % Top level directories
        case 'biht'
            addpath([repo_dir '/biht']);
        case 'qbcsspl'
            addpath([repo_dir '/qbcsspl']);
        case 'experiments'
            addpath([repo_dir '/experiments']); 
            addpath([repo_dir '/experiments/experiment_parameters/defaults']);  
            addpath([repo_dir '/experiments/experiment_parameters/1bit']);  
        case 'proj'
            addpath([repo_dir '/projections']);

        otherwise
            return_str = sprintf('Unknown dependency: %s',varargin{i});
            warning('csq_deps:Uknown',return_str);
    end
end

