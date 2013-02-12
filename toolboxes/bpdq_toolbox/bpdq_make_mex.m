% bpdq_make_mex : compiles mex files used in bpdq code
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

%
% This runs mex from within matlab. The source can also be 
% compiled from the command line, using the Makefile
% in the common/mex directory, 
% by setting location of the matlab mex script in the Makefile.
bpdqdir=pwd;
fprintf('Changing to directory mex\n');
cd mex
fprintf('Compiling mex functions proj_lpball_newton_mex and prox_tv_mex\n');
mex bpdq_proj_lpball_mex.cpp LpNormalizer.cpp LpProjector.cpp
mex bpdq_prox_tv_mex.c

cd(bpdqdir);

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
