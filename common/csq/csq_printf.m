function csq_printf(str_format,varargin)
% csq_printf("string",var1,var2,....)
%
% A wrapper function for printf which ensures that 
% immediate line flushing occurs for Octave as well as 
% Matlab.

fprintf(str_format,varargin{:});

% Force flushing in octave
if csq_in_octave()
	% This is an octave only function. If Matlab gets to
	% here it will crash. 
	fflush(stdout);
end
