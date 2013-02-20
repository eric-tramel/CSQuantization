function cvrms = CVRMS(x,y)
% cvrms = CVRMS(x,y)
%                 ^--- Observed
% Calcualte the normalized root meat square between two vectors. 

rms = RMS(x,y);
cvrms = rms ./ mean(y(:));