function h = fspecial(type,P1,P2)
%FSPECIAL Special 2-D filters.
%	H = FSPECIAL('gaussian',N,STD) returns a rotational symmetric
%	  Gaussian N-by-N filter with mainlobe width STD.
%	H = FSPECIAL('sobel') returns a 3-by-3 Sobel y-derivative filter.
%	H = FSPECIAL('prewitt') returns a 3-by-3 Prewitt y-derivative
%	  filter.
%	H = FSPECIAL('laplacian',ALPHA) returns a 3-by-3 Laplacian
%	  filter.  The parameter ALPHA controls the shape of the
%	  Laplacian and must be between 0 and 1. The default is 1/5.
%	H = FSPECIAL('log',N,STD) returns an N-by-N Laplacian of the
%	  Gaussian filter.
%	H = FSPECIAL('average',N) returns an N-by-N averaging filter,
%	  ONES(N,N)/(N*N).
%	H = FSPECIAL('unsharp',ALPHA) returns a 3-by-3 unsharp
%	  contrast enhancement filter.  The unsharp filter is 
%	  created from the Laplacian filter.  ALPHA controls the 
%	  shape of the Laplacian as for 'laplacian' above.
%
%	X-derivative filters are created by using -H'.  M-by-N
%	filters are returned when N is a 2-by-1 size vector.

%	Clay M. Thompson 11-17-92
%	Copyright (c) 1992 by The MathWorks, Inc.
%	$Revision: 1.9 $  $Date: 1993/08/24 20:24:04 $

if nargin==0, error('Not enough input arguments.'); end
type = [type,'  '];
code = lower(type(1:2));
if nargin>1,
  if ~(all(size(P1)==[1 1]) | all(size(P1)==[1 2])),
     error('The second parameter must be a scalar or a 1-by-2 size vector.');
  end
  if length(P1)==1, siz = [P1 P1]; else siz = P1; end
end

if all(code=='ga'), % Gaussian filter
  if nargin<2, siz = [3 3]; end
  if nargin<3, std = .5; else std = P2; end
  [x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(1)-1)/2:(siz(1)-1)/2);
  h = exp(-(x.*x + y.*y)/(2*pi*std*std));
  h = h/sum(sum(h));

elseif all(code=='so'), % Sobel filter
  h = [1 2 1;0 0 0;-1 -2 -1];

elseif all(code=='pr'), % Prewitt filter
  h = [1 1 1;0 0 0;-1 -1 -1];

elseif all(code=='la'), % Laplacian filter
  if nargin<2, alpha = 1/5; else alpha = P1; end
  alpha = max(0,min(alpha,1));
  h1 = alpha/(alpha+1); h2 = (1-alpha)/(alpha+1);
  h = [h1 h2 h1;h2 -4/(alpha+1) h2;h1 h2 h1];

elseif all(code=='lo'), % Laplacian of Gaussian
  if nargin<2, siz = [5 5]; end
  if nargin<3, std = .5; else std = P2; end
  [x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(1)-1)/2:(siz(1)-1)/2);
  std2 = std*std;
  h1 = exp(-(x.*x + y.*y)/(2*pi*std2));
  h = h1.*(x.*x + y.*y - 2*pi*std2)/(((pi*std2).^2)*sum(sum(h1)));

elseif all(code=='av'), % Smoothing filter
  if nargin<2, siz = [3 3]; end
  h = ones(siz)/prod(siz);

elseif all(code=='un'), % Unsharp filter
  if nargin<2, alpha = 1/5; else alpha = P1; end
  h = [0 0 0;0 1 0;0 0 0] - fspecial('laplacian',alpha);

else
  error('Unknown filter type.');

end