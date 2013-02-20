function nrms = NRMS(x,y)
% nrms = NRMS(x,y)
%               ^--- Observed
% Calcualte the normalized root meat square between two vectors. 

rms = RMS(x,y);
nrms = rms ./ (max(y(:)) - min(y(:)));