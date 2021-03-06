% argselectAssign : Assign variables in calling workspace
%
%  function argselectAssign(variable_value_pairs)
%  
%  Inputs : 
%  variable_value_pairs is a cell list of form
%  'variable1',value1,'variable2',value2,...
%  This function assigns variable1=value1 ... etc in the *callers* workspace
%
%  This is used at beginning of function to simulate keyword argument
%  passing. Typical usage is
%
%  argselectAssign(control_params);
%  argselectCheck(control_params,varargin);
%  argselectAssign(varargin);
%
%  where control_params is a cell list of variable,value pairs containing
%  the default parameter values.
% 
% See also argselectCheck
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function argselectAssign(variable_value_pairs)
for j =1:2:length(variable_value_pairs)
  pname=variable_value_pairs{j};
  pval=variable_value_pairs{j+1};
  assignin('caller',pname,pval);
end

% The BPDQ Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%  
% The BPDQ Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License
% along with The BPDQ Toolbox.  If not, see <http://www.gnu.org/licenses/>.
