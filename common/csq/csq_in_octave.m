function in = csq_in_octave()
% in = csq_in_octave()
%
% Function returns TRUE if called from Octave and FALSE 
% if called from Matlab.
% This code is based on suggestions made by David Grohmann,
% Tom Holroyd, and David Bateman on <octave-help> 
% (http://octave.1599824.n4.nabble.com/How-to-determine-if-you-are-in-octave-or-matlab-td1624960.html).

 persistent csqtmpinout;

 if isempty(csqtmpinout),
   csqtmpinout = exist('OCTAVE_VERSION','builtin') ~= 0;
 end;
 in = csqtmpinout;

return;
