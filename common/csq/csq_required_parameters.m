function csq_required_parameters(params,varargin)
% csq_required_parameters(params, 'field1', 'field2, ...)
% This function will throw an error if the structure 'params' does not
% contain any one of the fields specified as strings in varagin.

for i=1:length(varargin)
    if ~isfield(params,varargin{i})
        return_str = sprintf('Missing required parameter: %s',varargin{i});
        error('required_parameter:MissingParameter',return_str);
    end
end