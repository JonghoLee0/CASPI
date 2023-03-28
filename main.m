% Minimum working code for the reference:
%
% "CASPI: Colloaborative Photon Processing for Active Single-Photon Imaging"
% Jongho Lee, Atul Ingle, JV Chacko, KW Eliceiri, and Mohit Gupta
% Nature Communications, 2023

clear all;



%% Scenes
% application = 'LiDAR';                      % CASPI application
% scene = 'art';                              % scene name
% delta_t = 80e-12;                           % time bin size (s)
% R_inter = 0;                                % search size over cubes; S_inter = 2*R_inter + 1
% FWHM = 400e-12;                             % full-width at half-maximum of Gaussian laser pulse (s)


application = 'FLIM';                       % CASPI application
scene = 'autoFluorescence1';                % scene name
delta_t = 48e-12;                           % time bin size (s)
R_inter = 1;                                % search size over cubes; S_inter = 2*R_inter + 1
FWHM = 400e-12;                             % full-width at half-maximum of Gaussian laser pulse (s)



%% Parameters
pseudo_int = 1;                             % 1:pseudo intensity, 0:true intensity

skp = 1;                                    % pixel skip size (e.g., 1:dense processing, 2:processing for every 2nd pixel (both row and column directions))
N_sim = 10;                                 % number of similar cubelets
C_size = 8;                                 % photon cubelet size (C_size = C_x = C_y)
R_intra = 10;                               % search size within cube; S_intra = 2*R_intra + 1

R_th_LC = 0.8;                              % threshold to determine between thresholding & guided photon processing (3D photon transient cubelet)
R_th_NLC = 0.9;                             % threshold to determine between thresholding & guided photon processing (4D photon transient set)

N_cycle = 1000;                             % number of laser cycles
c = 3e8;                                    % light speed (m/s)
int_prct = 99;                              % 100: normalize intensity with maximum value

sigma_t = FWHM/(2*sqrt(2*log(2)));          % standard deviation of Gaussian laser pulse in time domain
sigma_f = 1/(2*pi*sigma_t);                 % standard deviation of Gaussian laser pulse in freq domain



%% Load data
fprintf('Loading data:\n\n');

load(['Data\', application, '\', scene]);


% photon data size
N_y = size(hst_map_set, 1);                 % number of rows of photon transient cube
N_x = size(hst_map_set, 2);                 % number of columns of photon transient cube
N_t = size(hst_map_set, 3);                 % number of time bins of photon transient
N_c = size(hst_map_set, 4);                 % number of photon transient cubes

delta_f = 1/(delta_t*N_t);                  % frequency bin size
N_sig_f = ceil(3*sigma_f/delta_f);          % number of frequency bins of laser pulse



%% Coates correction to reduce pileup
disp('Coates correction:');
tic

flux_noisy_map_set = zeros(N_y, N_x, N_t, N_c);

for frame = 1 : N_c
    
    
    % histogram
    hst_map = hst_map_set(:, :, :, frame);
    
    
    % coates correction
    flux_noisy_map = coates_grp(hst_map, N_cycle);
    flux_noisy_map_set(:, :, :, frame) = flux_noisy_map;
    
end

run_time1 = toc;
fprintf('%e %s\n\n', run_time1, 'sec');



%% Photon processing using local correlations
disp('Photon processing using LC:');
tic


% 1D FFT
N_t_half = round((N_t - 1)/2) + 1;
FLUX_noisy_map_set1D = fft(flux_noisy_map_set, N_t, 3);
FLUX_noisy_map_set1D = FLUX_noisy_map_set1D(:, :, 1:N_t_half, :);

if pseudo_int

    pseudo_i_map_set = zeros(N_y, N_x, N_c);
    
    for frame = 1 : N_c
        
        fprintf(['frame ', num2str(frame), ': ']);


        % pseudo intensity from 1D processing result
        FLUX_noisy_map1D = FLUX_noisy_map_set1D(:, :, :, frame);
        flux_map1D = ifft(FLUX_noisy_map1D(:, :, 1:N_sig_f), N_t, 3, 'symmetric');
        flux_map1D(flux_map1D < 0) = 0;
        i_map1D = sum(flux_map1D, 3);
        i_map1D = i_map1D/prctile(i_map1D(:), int_prct);

        
        % initial flux estimation
        fprintf('Initial flux estimation: ');
        [flux_map_init_LC, MSE_map, GP_mask] = estimate_init_flux_LC(FLUX_noisy_map1D, i_map1D, R_th_LC, N_y, N_x, N_t, C_size, N_sig_f, skp);
        
          
        % Wiener filtering
        fprintf('Wiener filtering: \n');
        flux_map_wiener_LC = estimate_wiener_LC(flux_map_init_LC, FLUX_noisy_map1D, MSE_map, N_y, N_x, N_t, N_t_half, C_size, N_sig_f, skp);


        % pseudo intensity
        pseudo_i_map = sum(flux_map_wiener_LC, 3);
        pseudo_i_map = pseudo_i_map/prctile(pseudo_i_map(:), int_prct);
        pseudo_i_map_set(:, :, frame) = pseudo_i_map;
        
    end
else
    
    pseudo_i_map_set = i_map_set;
    
end

run_time2 = toc;
fprintf('%e %s\n\n', run_time2, 'sec');



%% Find similar patches
disp('Finding similar patches:');
tic

[x_sim_map_set, y_sim_map_set, frame_sim_map_set] = find_similar_patches(pseudo_i_map_set, N_x, N_y, N_c, N_sim, C_size, R_intra, R_inter);

run_time3 = toc;
fprintf('%e %s\n\n', run_time3, 'sec');



%% Photon processing using local and non-local correlations

% initial flux estimation
disp('Initial flux estimation using LC and NLC:');
tic

[flux_map_set_init_NLC, MSE_map_set, GP_mask_set] = estimate_init_flux_NLC(FLUX_noisy_map_set1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, pseudo_i_map_set, R_th_NLC, N_y, N_x, N_t, N_t_half, N_sim, N_c, C_size, N_sig_f, skp);

run_time4 = toc;
fprintf('%e %s\n\n', run_time4, 'sec');




% Wiener filtering
disp('Wiener filtering using LC and NLC:');
tic

flux_map_set_wiener_NLC = estimate_wiener_NLC(flux_map_set_init_NLC, FLUX_noisy_map_set1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, MSE_map_set, N_y, N_x, N_t, N_t_half, N_sim, N_c, C_size, N_sig_f, skp);

run_time5 = toc;
fprintf('%e %s\n', run_time5, 'sec');



% total run time
run_time = run_time1 + run_time2 + run_time3 + run_time4 + run_time5;
fprintf('\n%s %e %s\n\n', 'Total run-time: ', run_time/60, 'min');



%% Save recovered photon fluxes
disp('Saving recovered fluxes:');
flux_recovered_map_set = flux_map_set_wiener_NLC;
save(['Data\', application, '\flux_recovered.mat'], 'flux_recovered_map_set');



%% LiDAR results
if strcmp(application, 'LiDAR')
    
    disp('Saving depth maps:');
    
    
    % get normalized laser pulse
    psf = get_laser_pulse(sigma_t, delta_t);

    
    % load colormap
    load('colormap1.mat');

    
    % depth range
    d_min = 1.40;
    d_max = 2.17;
    
    
    % conventional matched filtering (MF)
    [d_MF, d_color_MF] = get_d_map_set(flux_noisy_map_set, psf, c, delta_t, cmap, d_min, d_max);

    
    % CASPI + MF
    [d_CASPI, d_color_CASPI] = get_d_map_set(flux_recovered_map_set, psf, c, delta_t, cmap, d_min, d_max);
    

    % save depth maps
    for frame = 1 : N_c
        
        % depth maps
        imwrite(d_color_MF{frame}, ['Data\LiDAR\d_MF_frame', num2str(frame), '.png']);
        imwrite(d_color_CASPI{frame}, ['Data\LiDAR\d_CASPI_frame', num2str(frame), '.png']);
        
        
        % intensity images
        i_raw = sum(hst_map_set(:, :, :, frame), 3);
        i_raw = i_raw/prctile(i_raw(:), int_prct);
        i_CASPI = sum(flux_recovered_map_set(:, :, :, frame), 3);
        i_CASPI = i_CASPI/prctile(i_CASPI(:), int_prct);
        
        imwrite(i_raw, ['Data\LiDAR\i_raw_frame', num2str(frame), '.png']);
        imwrite(i_CASPI, ['Data\LiDAR\i_CASPI_frame', num2str(frame), '.png']);
        
    end
    
    
    % compare between photon fluxes
    p_y = 41; p_x = 36;
    hst_raw = squeeze(hst_map_set(p_y, p_x, :, 1));
    flux_recovered = squeeze(flux_recovered_map_set(p_y, p_x, :, 1));
    figure; plot(hst_raw); axis tight; title('Raw histogram')
    figure; plot(flux_recovered); axis tight; title('Recovered by CASPI')
    
    
    % compare between depth maps
    figure; imshow(d_color_MF{1}); title('Raw histogram + MF')
    figure; imshow(d_color_CASPI{1}); title('CASPI + MF')
    
end



%% FLIM results
if strcmp(application, 'FLIM')
    
    disp('Saving lifetime images:');
    
    
    % load colormap
    load('colormap2.mat');
    cmap = flipud(cmap);
    
    
    % lifetime range
    t_min = 1.1e-9;
    t_max = 2.0e-9;

    
    % intensity adjustment
    int_prct = 95;
    int_exp = 1.7;
    
    
    % number of samples for fitting
    N_sample = 35;

    
    % spatial binning + naive linear fitting (to log-transformed data)
    binning_size = 3;   % 3 means 7 x 7 spatial binning
    [t_binning, t_color_binning] = get_t_map_set(hst_map_set, int_prct, int_exp, N_sample, binning_size, delta_t, cmap, t_min, t_max);
    
    
    % CASPI + naive linear fitting
    binning_size = 0;
    [t_CASPI, t_color_CASPI] = get_t_map_set(flux_recovered_map_set, int_prct, int_exp, N_sample, binning_size, delta_t, cmap, t_min, t_max);
    
    
    
    % save lifetime images
    for frame = 1 : N_c
        
        imwrite(t_color_binning{frame}, ['Data\FLIM\binning_frame', num2str(frame), '.png']);
        imwrite(t_color_CASPI{frame}, ['Data\FLIM\CASPI_frame', num2str(frame), '.png']);
    end
    
    
    % compare between photon fluxes
    p_y = 200; p_x = 210;
    hst_raw = squeeze(hst_map_set(p_y, p_x, :, 1));
    flux_recovered = squeeze(flux_recovered_map_set(p_y, p_x, :, 1));
    figure; plot(hst_raw); axis tight; title('Raw histogram')
    figure; plot(flux_recovered); axis tight; title('Recovered by CASPI')
    
    
    % compare between lifetime images
    figure; imshow(t_color_binning{1}); title('Raw histogram + spatial binning + linear fit')
    figure; imshow(t_color_CASPI{1}); title('CASPI + linear fit')
   
end