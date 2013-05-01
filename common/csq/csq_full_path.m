function full_path = csq_full_path(local_path)

repo_dir = csq_get_repo_dir();

full_path = [repo_dir '/' local_path];