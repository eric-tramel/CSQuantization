% bpdq_random_fourier_locations : Compute random Fourier locations
%
% [useind,loc]=bpdq_random_fourier_locations(dim,K,varargin)
%
% Choses 2*K (or 2*K-1 if center is chosen) random locations 
% respecting fft symmetry,  i.e. if x of dimension dim is real
% then ifft2( ifftshift( fftshift(fft2(x)) .* loc ) )
% will also be real.
%
% Inputs:
% dim - dimensions of image
% K - # of independent locations desired
%
% Selectable Control Parameters :
% forcecenter - set to 1 to enforce DC coefficient to be chosen
% seed - pseudorandom seed, set to a real number to seed random number
% generator
%
% Outputs:
% useind : indices of Fourier coefficient locations
% loc : binary mask of size dim indicating Fourier coefficient locations
%
% This file is part of BPDQ Toolbox (Basis Pursuit DeQuantizer)
% Copyright (C) 2009, the BPDQ Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

function [useind,loc]=bpdq_random_fourier_locations(dim,K,varargin)
  control_params={'forcecenter',1,...
                 'seed',nan};
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);
  
  if ~isnan(seed)
    rand('seed',seed);
  end
  
  % odim : smallest even dimension less than or equal to dim
  odim=floor((dim-1)/2)*2+1;

  if (2*K>prod(odim))
    error('too many locations wanted')
  end
  
  % pick random locations from odim
  rloc=rand(odim);
  rloc=rloc+flipud(fliplr(rloc));
  if forcecenter==1
    ctr=floor((odim+1)/2);
    rloc(ctr)=-1; % smaller, will be picked first
  end
  
  [v,ind]=sort(rloc(:));
  thresh=v(2*K);
  eind=ind(find(v<=thresh));
  oloc=zeros(odim);
  oloc(eind)=1;
    
  % then pad with zeros to recover size dim
  % If I do this correctly, the center of size odim matrix
  % will be center of size dim matrix
  % and zero padded regions will correspond to highest frequencies
  % that must be real-valued for fft to be real
  
  if(odim(1)~=dim(1))
    % pad on top
    oloc=[zeros(1,size(oloc,2));oloc];
  end
  if(odim(2)~=dim(2))
    oloc=[zeros(size(oloc,1),1),oloc];
    % pad on left
  end
  loc=oloc;
  useind=find(loc);

  
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
