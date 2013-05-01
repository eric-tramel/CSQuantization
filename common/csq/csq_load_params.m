function params = csq_load_params(filename)

% Get to the proper path directory
csq_deps('experiments');

% Get the repository path name
repo_dir = csq_get_repo_dir();

% Load in the parameters
run([repo_dir '/experiments/experiment_parameters/' filename]);