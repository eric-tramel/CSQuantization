function X = csq_load_data(datatype,name)
% X = csq_load_data(datatype,name)
% This function loads in data from the data directory and outputs
% it as X. The particular datatype directory as well as the name
% of the file should be given.

repo_dir = csq_get_repo_dir();

switch datatype
case 'image'
	filename = sprintf('%s/data/image/%s',repo_dir,name);
	if exist(filename,'file')
		X = double(imread(filename));
	else
		return_str = sprintf('%s not found.',filename);
		error('load_data:MissingData',return_str);
	end

otherwise
	return_str = sprintf('Dataset of type "%s" is unknown.',datatype);
	error('load_data:UknownType',return_str);
end
