function image_ratedistortion_experiment(image,bitrates,filename,module,params)
% image_ratedistortion_experiment(image,bitrates,params)
%
% Generates a single RD curve for the given image over the given bitrates
% given in bits per pixel (bpp). Final curve results are saved to specified
% filename.
csq_deps('common','ssim');


% Store repository information
%  We are doing this early in case somehow the repository changes while
%  running. This probably shouldn't happen, but lets be safe.
repo = csq_get_repo_info();

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
distortion.psnr_curve = zeros(rd_points,1);
distortion.snr_curve = zeros(rd_points,1);
distortion.rms_curve = zeros(rd_points,1);
distortion.ssim_curve = zeros(rd_points,1);

% Main Loop
%  This experiment iterates over specified target bitrates
for trial = 1:params.trials
  if params.trials > 1 
  % If there are multiple trials being run, we need to make sure to
  % update the random seed. Otherwise, leave the random seed alone.
    csq_reset_seed();
  end

  for rate_idx=1:rd_points
      params.experiment.target_bitrate = bitrates(rate_idx);
      [XF point_results] = module(X,params);
      
      % Set full raw results
      full_results(rate_idx) = point_results;
      
      % Save data to curves
      iterations(rate_idx) = iterations(rate_idx) + point_results.iterations;
      run_times(rate_idx)  = run_times(rate_idx)  + point_results.run_time;
      rates(rate_idx)      = rates(rate_idx)      + point_results.true_bitrate;

      distortion.mse_curve(rate_idx)  = distortion.mse_curve(rate_idx)  + MSE(X(:),XF(:));
      distortion.snr_curve(rate_idx)  = distortion.snr_curve(rate_idx)  + SNR(X(:),XF(:));
      distortion.psnr_curve(rate_idx) = distortion.psnr_curve(rate_idx) + PSNR(X(:),XF(:));
      distortion.rms_curve(rate_idx)  = distortion.rms_curve(rate_idx)  + RMS(X(:),XF(:));
      distortion.ssim_curve(rate_idx) = distortion.ssim_curve(rate_idx) + ssim(X,XF);
  end
end

% Average the results
  iterations = iterations./ params.trials;
  run_times  = run_times ./ params.trials;
  rates      = rates     ./ params.trials;

  distortion.mse_curve  = distortion.mse_curve  ./ params.trials;
  distortion.snr_curve  = distortion.snr_curve  ./ params.trials;
  distortion.psnr_curve = distortion.psnr_curve ./ params.trials;
  distortion.rms_curve  = distortion.rms_curve  ./ params.trials;
  distortion.ssim_curve = distortion.ssim_curve ./ params.trials;



%% Save Results
% Store inputs
inputs.image = image;
inputs.bitrates = bitrates;
inputs.params = params;
inputs.module_name = func2str(module);

save_date = date;

save(filename,'inputs',...                  % Experiment inputs
              'repo',...                    % Code repository information (for reproduction)
              'save_date',...               % Date of experiment completion
              'rates', ...                  % True experiment bitrates
              'run_times',...               % Experiment module run time at each bitrate
              'iterations',...              % Recovery iterations at each bitrate
              'distortion',...              % Distortion metrics at each bitrate
              'full_results');              % Full result dump at each bitrate

