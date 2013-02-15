function repodir = csq_get_repo_dir()
% repodir = csq_get_repo_dir()
% Returns the full path name of the CSQ repository.

cd_str = cd;
repo_name = 'CSQuantization';
repodir = cd_str(1:(strfind(cd_str,repo_name)+length(repo_name)-1));
