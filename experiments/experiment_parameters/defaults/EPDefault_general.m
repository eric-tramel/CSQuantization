% calling_experiment.m
% 
% Demonstrates how to call the experiment modules.
csq_deps('common');

% General Settings
params.rand_seed = 1;
params.verbose = 1;

% Get the repository directory in case it is needed
repo_dir = csq_get_repo_dir();
