function PSNR_mc = run_simple_GOP_MC(sequence_name, subrate, total_num_frames)

addpath(genpath('C:\Dropbox\work\Common\Sequences'));
addpath(genpath('C:\Dropbox\work\Common\BCS-SPL'));
addpath(genpath('C:\Dropbox\work\Common\WaveletSoftware'));

key_subrate = 0.7;
subrate = 0.1;


isSave = 1;
sequence_name = 'cardiac';
GOP_size = 8;
total_num_frames = 24;

% global frames

if floor(GOP_size/2) ~= ceil(GOP_size/2)
    error('Frame number should be even.')
end
if total_num_frames < GOP_size
    error('Total number of frames should be larger than GOP size')
end
if mod(total_num_frames,GOP_size) ~= 0
    error('Total number of frames should be multiple of GOP size')
end

PSNR_mc = zeros(1,total_num_frames);
PSNR_indep = zeros(1,total_num_frames);
for k = 1:GOP_size:total_num_frames;
    frames = readframes(sequence_name, GOP_size, k);    
    [num_rows num_cols] = size(frames{1});
    % encoding
    [y Phis] = MC_CS_PL_Encoder(frames, key_subrate, subrate);
    
    % decoding
    tic
    [reconstructed_frames reconstructed_frames_indep] ...
        = MC_CS_PL_Decoder_GOP(y,Phis,GOP_size, num_rows, num_cols);
    toc
    disp('< PSNR result Comparison (in dB) >');
    for i = 1 : GOP_size;
        PSNR_mc(k+i-1) = PSNR(frames{i}, reconstructed_frames{i});
        PSNR_indep(k+i-1) = PSNR(frames{i}, reconstructed_frames_indep{i});
        disp(['Frame # '  num2str(k+i-1) '  Independent: ' num2str(PSNR_indep(k+i-1)) ...
            ', MC-CS-PL: ' num2str(PSNR_mc(k+i-1)) ]);
    end

    % save image result of center frame of GOP
    if (k < GOP_size)&&isSave
        file_name = [sequence_name '.00' num2str(round(GOP_size/2)) '_' ...
            num2str(key_subrate) '_' num2str(subrate) '.pgm'];
        imwrite(uint8(reconstructed_frames{5}),file_name);
    end
end
mean(PSNR_mc)
if isSave
    filename = ['results/' sequence_name '_mc_telescope_noiteration_weight' num2str(subrate) '.mat'];
    save(filename, 'PSNR_mc', 'PSNR_indep');
end


function frames = readframes(sequence_name, GOP_size, start_sequence)
frames = cell(1,GOP_size);
index = 1;
for i = start_sequence: start_sequence+GOP_size
    if i <= 10
        file_name = [ sequence_name '.00' num2str(i-1) '.pgm'];
    elseif (i>10) && (i<=100)
        file_name = [ sequence_name '.0' num2str(i-1) '.pgm'];
    else
        file_name = [ sequence_name '.' num2str(i-1) '.pgm'];
    end
    frames{index} = double(imread(file_name));
    index = index +1;
end
