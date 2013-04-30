function image_ratedistortion_experiment(image,bitrates,filename,module,params)
% image_ratedistortion_experiment(image,bitrates,params)
%
% Generates a single RD curve for the given image over the given bitrates
% given in bits per pixel (bpp). Final curve results are saved to specified
% filename.

% Make sure that we have git on the path
csq_deps('git');

% Load in image data
X = csq_load_data('image',image);

% Set some general parameters
params.imsize = size(X);
params.N = size(X,1)*size(X,2);

% Set up result curves
rd_points = length(bitrates);
rates = zeros(rd_points,1);
run_times = zeros(rd_points,1);
iterations = zeros(rd_points,1);
distortion.mse_curve = zeros(rd_points,1);
distortion.snr_curve = zeros(rd_points,1);
distortion.rms_curve = zeros(rd_points,1);
distortion.ssim_curve = zeros(rd_points,1);

for rate_idx=1:rd_points
    params.experiment.target_bitrate = bitrates(rate_idx);
    [XF point_results] = module(X,params);
    
    % Set full raw results
    full_results(rate_idx) = point_results;
    
    % Save data to curves
    iterations(rate_idx) = point_results.iterations;
    run_times(rate_idx) = point_results.run_time;
    rates(rate_idx) = point_results.true_bitrate;
    distortion.mse_curve(rate_idx)  = MSE(X(:),XF(:));
    distortion.snr_curve(rate_idx)  = SNR(X(:),XF(:));
    distortion.psnr_curve(rate_idx)  = PSNR(X(:),XF(:));
    distortion.rms_curve(rate_idx)  = RMS(X(:),XF(:));
    distortion.ssim_curve(rate_idx) = ssim(X,XF);
end




%% Save Results
% Store inputs
inputs.image = image;
inputs.bitrates = bitrates;
inputs.params = params;
inputs.module_name = func2str(module);

% Store repository information
repo.branch_name = git('rev-parse','--abbrev-ref','HEAD');
repo.branch_version = git('rev-parse','HEAD');

save_date = date;


save(filename,'inputs',...                  % Experiment inputs
              'repo',...                    % Code repository information (for reproduction)
              'save_date',...               % Date of experiment completion
              'rates', ...                  % True experiment bitrates
              'run_times',...               % Experiment module run time at each bitrate
              'iterations',...              % Recovery iterations at each bitrate
              'distortion',...              % Distortion metrics at each bitrate
              'full_results');              % Full result dump at each bitrate

