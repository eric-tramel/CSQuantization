function smoothing = csq_generate_smoothing(params)
%CSQ_GENERATE_SMOOTHING return function handle for smoothing function.
%
% function smoothing = csq_generate_smoothing(params)
%
% This funciton generates function handles for a smoothing function
% based upon the specified parameters. Different types of smoothing
% require differing parametesr for proper functionality.
%
%	Inputs:
%		params 		- Structure containing function parameters.
%					  Required parameters:
%					  params.imsize	
%					  params.block_based
%					  params.block_dim
%
%		Smoothing specific parameters
%		'weiner' : 	params.window_dim,
%		'deblock':	params.radius
%
%
%	Outputs:
%		smoothing 	- Function handle of the form smoothing(x) where
%					  x is a rasterized image vector.
%
% Written by: Eric W. Tramel, Ph.D.
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%% Check Inputs
csq_required_parameters(params,'imsize','block_based');
csq_required_parameters(params.smoothing,'id');

if params.block_based
	csq_required_parameters(params,'block_dim');
end

%% Default parameter values
if ~isfield(params.smoothing,'radius')
	params.smoothing.radius = 2;
end

if ~isfield(params.smoothing,'window_dim')
	params.smoothing.window_dim = [3 3];
end

%% Generate Function Handles
switch params.smoothing.id
case 'weiner'
	% Just a spelling check ;)
	smoothing = csq_generate_smoothing('wiener',params);

case 'wiener'
	smoothing = @(z) csq_vectorize( wiener2( reshape(z,params.imsize),...
										     params.smoothing.window_dim) );

case 'deblock'
	if params.block_based
	    smoothing = @(z) csq_vectorize( deblocking_filter( reshape(z,params.imsize),...
														   params.block_dim,...
														   params.smoothing.radius));
	else
		return_str = sprintf('Should not use deblocking if image is not block based.');
		error('csq_generate_smoothing:Blocking',return_str);
	end

case 'none'
	smoothing = @(z) z;

otherwise
	return_str = sprintf('Threshold "%s" is unsupported.',params.smoothing.id);
    error('csq_generate_smoothing:UnsupportedSmoothing',return_str);
end
