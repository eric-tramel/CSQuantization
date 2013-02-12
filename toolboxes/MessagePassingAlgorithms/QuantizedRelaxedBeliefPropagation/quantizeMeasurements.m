function q = quantizeMeasurements(u, quantizer)
% Quantizes the vector u according to quantizer object. Returns output
% values, left and right boudaries of the cell for each u(i).

% Number of measurements
M = length(u);

% Quantized output and cell boundaries
q = zeros(M, 1);

% For each cell find values in it
for i = 1:quantizer.nLevels
    % Cell boundaries
    left = quantizer.x(i);
    right = quantizer.x(i+1);
    
    % Indices
    ix = (u > left) & (u <= right);
    
    % Bin number
    b = quantizer.binTable(i, 2);
    
    q(ix) = b;
end