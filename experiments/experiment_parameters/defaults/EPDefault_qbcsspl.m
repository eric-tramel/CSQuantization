% Default parameter settings for using Q-BCS-SPL recovery
clear
csq_deps('common','qbcsspl');

% Set Q-BCS-SPL Parameters
params.qbcsspl.maxIter = 400;
params.qbcsspl.tol = 0.0001;
params.qbcsspl.quant = 'dpcm';
params.qbcsspl.meanSubtraction = 0;


