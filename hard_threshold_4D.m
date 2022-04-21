function [flux_map_set_ht4D, MSE_map_set, IBFmask4D_set] = hard_threshold_4D(FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, i_map_set, IBF, th_IBF, size_y, size_x, N_bin, N_bin_half, N_sim, N_frame, s_patch, N_sig_f, skp, octave)
                                                                            

% buffers
flux_map_set_ht4D = zeros(size_y, size_x, N_bin, N_frame);
weight_sum_map_set = zeros(size_y, size_x, N_frame);
MSE_map_set = nan(size_y, size_x, N_frame);
IBFmask4D_set = zeros(size_y, size_x, N_frame);


if octave
    
    [y_grid, x_grid, t_grid, s_grid] = ndgrid([1:s_patch], [1:s_patch], [1:N_sig_f], [1:round((N_sim-1)/2)+1]);

    y_grid = y_grid(:);
    x_grid = x_grid(:);
    t_grid = t_grid(:);
    s_grid = s_grid(:);

    y_grid2 = mod(s_patch - (y_grid-1), s_patch) + 1;
    x_grid2 = mod(s_patch - (x_grid-1), s_patch) + 1;
    t_grid2 = mod(N_bin - (t_grid-1), N_bin) + 1;
    s_grid2 = mod(N_sim - (s_grid-1), N_sim) + 1;
    
end


%% Loop
for frame = 1 : N_frame
    
    fprintf('frame %d:\n', frame);
    
    for y = 1 : skp : size_y - (s_patch-1)
        for x = 1 : skp :  size_x - (s_patch-1) 



            %% Collect similar patches
            FLUX_patchset_noisy3D = zeros(s_patch, s_patch, N_bin_half, N_sim);

            for sim_idx = 1 : N_sim

                y_sim = y_sim_map_set(y, x, sim_idx, frame);
                x_sim = x_sim_map_set(y, x, sim_idx, frame);
                frame_sim = frame_sim_map_set(y, x, sim_idx, frame);

                FLUX_patch_noisy1D = FLUX_map_set_noisy1D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim);


                % 3D FFT
                FLUX_patchset_noisy3D(:, :, :, sim_idx) = fft2(FLUX_patch_noisy1D);

            end



            %% 4D FFT    
            FLUX_patchset_noisy4D = fft(FLUX_patchset_noisy3D, N_sim, 4);        



            %% Signal & Noise separation
            % signal
            SIG_patchset_noisy4D = FLUX_patchset_noisy4D(:, :, 1:N_sig_f, :);


            % noise
            NOISE_patchset_noisy4D = FLUX_patchset_noisy4D(:, :, N_sig_f+1:end, :);



            %% Noise statistics
            % mean noise & mean squared noise
            ME = mean(abs(NOISE_patchset_noisy4D(:)));
            MSE = mean(abs(NOISE_patchset_noisy4D(:)).^2);

            MSE_map_set(y, x, frame) = MSE;


            % noise threshold
            if N_frame == 1 && N_sim > 1
                noise_sigma = 6;
            else
                noise_sigma = 3;
            end
            
            noise_scale = 1 + noise_sigma*sqrt((gamma(1)/gamma(1.5))^2 - 1);
            noise_th = noise_scale*ME;
            
            
%             %
%             figure; hold on; grid on;
%             plot(abs(FLUX_patchset_noisy4D(:)), 'color', [102 153 255]/255, 'linewidth', 6);
%             yline(noise_th, 'color', [255 0 0]/255, 'lineWidth', 6);
%             axis tight
%             ylim([0 0.4])
%             xticklabels(''); yticklabels('');



            %% Intensity-based filtering
            if IBF      

                % noise-to-signal ratio
                N_power = MSE;
                S_power = mean(abs(SIG_patchset_noisy4D(:)).^2);               % pseudo signal power before hard-thresholding


                if N_power/S_power >= th_IBF

                    IBFmask4D_set(y, x, frame) = 1;


                    % intensity patch
                    i_patch = i_map_set(y : y+(s_patch-1), x : x+(s_patch-1), frame);
                    i_patch = i_patch/sum(i_patch(:));
                    I_patch = abs(fftn(i_patch));


                    % intensity-based filtering
                    SIG_patchset_ht4D = SIG_patchset_noisy4D.*repmat(I_patch, 1, 1, N_sig_f, N_sim);

                else
                    % Hard thresholding
                    [SIG_patchset_ht4D, hard_coeff] = shrink_hard_sig(SIG_patchset_noisy4D, noise_th);
                end        
            else
                % Hard thresholding
                [SIG_patchset_ht4D, hard_coeff] = shrink_hard_sig(SIG_patchset_noisy4D, noise_th);
            end



            %% Inverse 4D FFT
            FLUX_patchset = zeros(s_patch, s_patch, N_bin, N_sim);
            FLUX_patchset(:, :, 1:N_sig_f, :) = SIG_patchset_ht4D;
            
            if octave
                
                FLUX_patchset(:, :, :, round((N_sim-1)/2)+2 : end) = 0;
                FLUX_patchset(sub2ind([s_patch, s_patch, N_bin, N_sim], y_grid2, x_grid2, t_grid2, s_grid2)) =  conj(FLUX_patchset(sub2ind([s_patch, s_patch, N_bin, N_sim], y_grid, x_grid, t_grid, s_grid)));
                flux_patchset = ifftn(FLUX_patchset);
                
            else % matlab
                
                flux_patchset = ifftn(FLUX_patchset, 'symmetric');
                
            end



            %% Aggregate
            weight = 1/(MSE + eps);

            for sim_idx = 1 : N_sim

                y_sim = y_sim_map_set(y, x, sim_idx, frame);
                x_sim = x_sim_map_set(y, x, sim_idx, frame);
                frame_sim = frame_sim_map_set(y, x, sim_idx, frame);

                weight_sum_map_set(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, frame_sim) = weight_sum_map_set(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, frame_sim) + weight;
                flux_map_set_ht4D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim) = flux_map_set_ht4D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim) + weight*flux_patchset(:,:,:, sim_idx);

            end
        end
    end
end



% filtered flux by hard thresholding
weight_sum_map_set = reshape(weight_sum_map_set, size_y, size_x, 1, N_frame);
flux_map_set_ht4D = flux_map_set_ht4D./repmat(weight_sum_map_set, 1, 1, N_bin, 1);
flux_map_set_ht4D(flux_map_set_ht4D < 0) = 0;
