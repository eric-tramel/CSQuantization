% calling_experiment.m
% 
% Demonstrates how to call the experiment modules.
clear
csq_deps('common','experiments');
repo_dir = csq_get_repo_dir();

% Parameter settings
run([repo_dir '/experiments/experiment_parameters/1bit/1bit_srm-full_dwt.m']);

params.trials = 2;

% Call the module
image_ratedistortion_experiment('lena.jpg',...
								linspace(0.05,2,5),...
								[repo_dir '/results/rd/lena_1bit-srm-full_dwt.mat'],...
								@experiment_module_biht2d,...
								params);
