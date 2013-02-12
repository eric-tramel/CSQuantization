% bpdq_compute_snr : Calculate SNR
%
% snr=bpdq_compute_snr(u,v)
% 
% Note : Signals are assumed to be zero mean, so signal and noise
% power computed without explicitly subtracting off the mean.
%
% Inputs :
% u,v input signals
%
% Outputs:
% snr - calculated SNR (in dB)
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function snr=bpdq_compute_snr(u,v)
  snr=20*log10(norm(u)/norm(u-v));
  
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
  