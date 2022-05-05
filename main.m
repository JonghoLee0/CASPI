clear all;




%% Scenes
app = 'LiDAR'; scene = 'Art'; t_bin = 80e-12; N_bin = 1024; N_sim = 10; s_patch = 8; r_intra = 10; r_inter = 0; FWHM = 400e-12; d_min = 1.40; d_max = 2.17; cmap_name = 'custom2';
% app = 'FLIM'; scene = 'autoFluorescence1'; t_bin = 48e-12; N_bin = 256; N_sim = 10; s_patch = 8; r_intra = 10; r_inter = 1; FWHM = 400e-12; t_min = 1.1e-9; t_max = 2.0e-9; cmap_name = 'custom4';




%% Parameters
octave = 0;                                 % 1:octave, 0:matlab
save_interm_int = 0;                        % 1:save all intermediate intensities, 0:don't save them


% filtering
skp = 1;                                    % filtering skip size, 1 means dense filtering
s_intra = 2*r_intra + 1;                    % search size within a frame


% intensity-based filtering
pseudo_int = 1;                             % 1:pseudo intensity, 0:true intensity
IBF = 1;                                    % 1:yes, 0:no
th_IBF3D = 0.8;                             % default:0.8
th_IBF4D = 0.9;                             % default:0.9
max_int_prct = 100;                         % 100: normalize intensity with max


% system
N_cycle = 1000;                             % laser cycle numbers
c = 3e8;                                    % light speed


% misc
d_range = c*N_bin*t_bin/2;                  % measurable depth range
sigma_t = FWHM/(2*sqrt(2*log(2)));          % std for the given FWHM in time domain
sigma_f = 1/(2*pi*sigma_t);                 % std for the given FWHM in freq domain
bin_f = 1/(t_bin*N_bin);                    % bin size in freq domain
N_sig_f = ceil(3*sigma_f/bin_f);            % number of signal bins in freq domain


sigma_bin = sigma_t/t_bin;
bin_max = 10*round(sigma_bin);
psf = normpdf([1:bin_max]', bin_max/2, sigma_bin);



%% Load data
% histogram
load([app, '\Histograms\', scene, '_hst.mat']);


% data size
size_y = size(hst_map_set, 1);
size_x = size(hst_map_set, 2);
N_frame = size(hst_map_set, 4);


% get pseudo intensity from histogram
if pseudo_int
    for frame = 1 : N_frame
    
    
        % histogram
        hst_map = hst_map_set(:, :, :, frame);
        
        
        % intensity
        i_map_hst = sum(hst_map, 3);
        i_map_hst = i_map_hst/prctile(i_map_hst(:), max_int_prct);


        % save
        if save_interm_int
            imwrite(i_map_hst, [app, '\Results\1_i_hst_', num2str(frame) , '.png']);
        end
    end
end



%% Coates Correction
fprintf('\nCoates Correction:\n'); tic

flux_map_set_noisy = zeros(size_y, size_x, N_bin, N_frame);

for frame = 1 : N_frame
    
    
    % histogram
    hst_map = hst_map_set(:, :, :, frame);
    
    
    % coates correction
    flux_map_noisy = coates_grp(hst_map, N_cycle);
    flux_map_set_noisy(:, :, :, frame) = flux_map_noisy;
    
    
    % get pseudo intensity from coates correction
    if pseudo_int
        
        i_map_coates = sum(flux_map_noisy, 3);
        i_map_coates = i_map_coates/prctile(i_map_coates(:), max_int_prct);

        
        % save
        if save_interm_int
            imwrite(i_map_coates, [app, '\Results\2_i_coates_', num2str(frame) , '.png']);
        end
    end
end


% compute run-time
run_time1 = toc; fprintf('%e %s\n\n', run_time1, 'sec');



%% 1D FFT
N_bin_half = round((N_bin-1)/2)+1;
FLUX_map_set_noisy1D = fft(flux_map_set_noisy, N_bin, 3);
FLUX_map_set_noisy1D = FLUX_map_set_noisy1D(:, :, 1:N_bin_half, :);



%% Get intensities after 3D filtering for similar patch finding
fprintf('Get intensities:\n'); tic

if pseudo_int

    i_map_set_wiener3D = zeros(size_y, size_x, N_frame);
    
    for frame = 1 : N_frame
        
        fprintf('frame %d: ', frame);


        % pseudo intensity from 1D frequency filtering
        FLUX_map_noisy1D = FLUX_map_set_noisy1D(:, :, :, frame);
        flux_map_1D = ifft(FLUX_map_noisy1D(:, :, 1:N_sig_f), N_bin, 3, 'symmetric');
        flux_map_1D(flux_map_1D < 0) = 0;
        i_map_1D = sum(flux_map_1D, 3);
        i_map_1D = i_map_1D/prctile(i_map_1D(:), max_int_prct);
        
        if save_interm_int
            imwrite(i_map_1D, [app, '\Results\3_i_1D_', num2str(frame) , '.png']);
        end

        
        % 3D Hard-thresholding
        fprintf('3D HT: ');
        [flux_map_ht3D, MSE_map, IBFmask3D] = hard_threshold_3D(FLUX_map_noisy1D, i_map_1D, IBF, th_IBF3D, size_y, size_x, N_bin, s_patch, N_sig_f, skp, octave);
        
        
        % pseudo intensity from 3D HT
        i_map_ht3D = sum(flux_map_ht3D, 3);
        i_map_ht3D = i_map_ht3D/prctile(i_map_ht3D(:), max_int_prct);
        

        % save
        if save_interm_int
            imwrite(i_map_ht3D, [app, '\Results\4_i_ht3D_', num2str(frame) , '.png']);
            imwrite(IBFmask3D, [app, '\Results\IBFmask3D_', num2str(frame), '.png']);
        end
        
                
        % 3D Wiener filtering
        fprintf('3D Wiener:\n');
        flux_map_wiener3D = wiener_filter_3D(flux_map_ht3D, FLUX_map_noisy1D, MSE_map, size_y, size_x, N_bin, N_bin_half, s_patch, N_sig_f, skp, octave);


        % pseudo intensity from 3D Wiener
        i_map_wiener3D = sum(flux_map_wiener3D, 3);
        i_map_wiener3D = i_map_wiener3D/prctile(i_map_wiener3D(:), max_int_prct);
        i_map_set_wiener3D(:, :, frame) = i_map_wiener3D;
        
        
        % save
        if save_interm_int
            imwrite(i_map_wiener3D, [app, '\Results\5_i_wiener3D_', num2str(frame) , '.png']);
        end
    end
else
    
    i_map_set_wiener3D = i_map_set;
    
end

% compute run-time
run_time2 = toc; fprintf('%e %s\n\n', run_time2, 'sec');



%% Find similar patches
fprintf('Find similar patches:\n'); tic
[x_sim_map_set, y_sim_map_set, frame_sim_map_set] = find_similar_patches(i_map_set_wiener3D, size_x, size_y, N_frame, N_sim, s_patch, r_intra, r_inter);


% compute run-time
run_time3 = toc; fprintf('%e %s\n\n', run_time3, 'sec');



%% 4D Hard-thresholding
fprintf('4D HT:\n'); tic


% hard-thresholding
[flux_map_set_ht4D, MSE_map_set, IBFmask4D_set] = hard_threshold_4D(FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, i_map_set_wiener3D, IBF, th_IBF4D, size_y, size_x, N_bin, N_bin_half, N_sim, N_frame, s_patch, N_sig_f, skp, octave);


% IBF masks
% for frame = 1 : N_frame
%     
%     imwrite(IBFmask4D_set(:,:,frame), [app, '\Results\IBFmask4D_', num2str(frame), '.png']);
%     
% end


% get pseudo intensities
if pseudo_int

    for frame = 1 : N_frame


        % flux
        flux_map_ht4D = flux_map_set_ht4D(:,:,:,frame);


        % pseudo intensity map from 4D HT
        i_map_ht4D = sum(flux_map_ht4D, 3);
        i_map_ht4D = i_map_ht4D/prctile(i_map_ht4D(:), max_int_prct);

        if save_interm_int
            imwrite(i_map_ht4D, [app, '\Results\6_i_ht4D_', num2str(frame) , '.png']);
        end
    end
end


% compute run-time
run_time4 = toc; fprintf('%e %s\n\n', run_time4, 'sec');



%% 4D Wiener filtering
fprintf('4D Wiener:\n'); tic


% wiener filtering
flux_map_set_wiener4D = wiener_filter_4D(flux_map_set_ht4D, FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, MSE_map_set, size_y, size_x, N_bin, N_bin_half, N_sim, N_frame, s_patch, N_sig_f, skp, octave);


% get pseudo intensities
if pseudo_int

    for frame = 1 : N_frame
        

        % flux
        flux_map_wiener4D = flux_map_set_wiener4D(:,:,:,frame);


        % pseudo intensity map from 4D HT
        i_map_wiener4D = sum(flux_map_wiener4D, 3);
        i_map_wiener4D = i_map_wiener4D/prctile(i_map_wiener4D(:), max_int_prct);

        if save_interm_int
            imwrite(i_map_wiener4D, [app, '\Results\7_i_wiener4D_', num2str(frame) , '.png']);
        end
    end
end


% compute run-time
run_time5 = toc; fprintf('%e %s\n', run_time5, 'sec');


% total run-time
total_run_time = run_time1 + run_time2 + run_time3 + run_time4 + run_time5;
fprintf('\n%s %e %s\n\n', 'Total run-time: ', total_run_time/3600, 'hrs');



%% Save filtered photon flux
save([app, '\Results\flux_hat.mat'], 'flux_map_set_wiener4D');



%% LiDAR results
if strcmp(app, 'LiDAR')
    
    
    % load ground-truth for error computation
    load([app, '\Groundtruth\', scene, '.mat'])
    load(['Color_maps/', cmap_name, '.mat']);

    
    % Conventional: matched filtering (MF)
    fprintf('Matched filtering:\n')
    method = 2;         % 1:max, 2:matched
    mask_set = ~isnan(d_map_set);
    d_matched_map_set = get_d_hat_map_set(flux_map_set_noisy, mask_set, psf, c, t_bin, cmap, d_min, d_max, method, app, 'matched');
    compute_metric_set(d_map_set, d_matched_map_set, mask_set)
    
      
    % CASPI + MF
    fprintf('CASPI + Matched filtering:\n')
    method = 2;         % 1:max, 2:matched
    mask_set = ~isnan(d_map_set);
    flux_map_set = flux_map_set_wiener4D;
    d_hat_map_set = get_d_hat_map_set(flux_map_set, mask_set, psf, c, t_bin, cmap, d_min, d_max, method, app, 'wiener4d');
    compute_metric_set(d_map_set, d_hat_map_set, mask_set)

end



%% FLIM results
if strcmp(app, 'FLIM')
    
    max_int_prct = 95;
    int_exp = 1.7;
    mask_th = 0.0;
    N_sample = 35;

    load(['Color_maps/', cmap_name, '.mat']);
    cmap = flipud(cmap);
    
    for frame = 1 : N_frame
        
        
        % Conventional: linear fitting to log-transformed data
        hst_map = hst_map_set(:, :, :, frame);
        i_map = sum(hst_map, 3);
        i_map = i_map/prctile(i_map(:), max_int_prct);
        i_map_cnv = i_map.^int_exp;
        mask_cnv = i_map_cnv > mask_th;
        
        
        % linear fitting after spatial binning
        [t_map, t_map_color, ti_map_color] = make_t_map(hst_map, i_map_cnv, mask_cnv, N_sample, t_bin, cmap, t_min, t_max, 3);
        t_map(t_map < t_min) = t_min;
        t_map(t_map > t_max) = t_max;
        
        t_map(~mask_cnv) = nan;
        t_map(t_map == 0) = nan;
        
        imwrite(ti_map_color, ['FLIM\Results\binning_linearfit_', num2str(frame), '.png'])
        
        
        
        % CASPI
        flux_map = flux_map_set_wiener4D(:, :, :, frame);
        i_map = sum(flux_map, 3);
        i_map = i_map/prctile(i_map(:), max_int_prct);
        i_map_csp = i_map.^int_exp;
        mask_csp = i_map_csp > mask_th;
        
        [t_map, t_map_color, ti_map_color] = make_t_map(flux_map, i_map_csp, mask_csp, N_sample, t_bin, cmap, t_min, t_max, 0);
        t_map(t_map < t_min) = t_min;
        t_map(t_map > t_max) = t_max;
        
        t_map(~mask_csp) = nan;
        t_map(t_map == 0) = nan;

        imwrite(ti_map_color, ['FLIM\Results\caspi_linearfit', num2str(frame), '.png'])
        
    end
end