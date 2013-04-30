function repo_info = csq_get_repo_info()
%CSQ_GET_REPO_INFO returns a structure containing informaiton about the git repo.
%
% function repo_info = csq_get_repo_info()
%
%	Outputs a structure whose fields contain data about the current 
%   repository state.
%
%	repo_info.branch_name:    The name of the active branch.
%   repo_info.branch_version: SHA-1 hash of the current commit.
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

% Make sure that the git package is included
csq_deps('git');

% Get the current branch name
repo_info.branch_name = git('rev-parse','--abbrev-ref','HEAD');

% Get the current SHA-1 for the commit
repo_info.branch_version = git('rev-parse','HEAD');