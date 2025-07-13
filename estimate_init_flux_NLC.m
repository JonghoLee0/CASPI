% Estimate initial fluxes using local and non-local correlations


function [flux_map_set, MSE_map_set, GP_mask_set] = estimate_init_flux_NLC(FLUX_noisy_map_set1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, i_map_set, R_th, N_y, N_x, N_t, N_t_half, N_sim, N_c, C_size, N_sig_f, skp)
                                                                            


%% Buffers
flux_map_set = zeros(N_y, N_x, N_t, N_c);
MSE_map_set = nan(N_y, N_x, N_c);
GP_mask_set = zeros(N_y, N_x, N_c);
weight_sum_map_set = zeros(N_y, N_x, N_c);



%% Initial flux estimation using local and non-local correlations
for frame = 1 : N_c
    
    disp(['frame ', num2str(frame), ':']);
    
    for y = 1 : skp : N_y - (C_size-1)
        for x = 1 : skp :  N_x - (C_size-1) 



            %% Collect similar patches
            FLUX_noisy_cubelet_set3D = zeros(C_size, C_size, N_t_half, N_sim);

            for sim_idx = 1 : N_sim

                y_sim = y_sim_map_set(y, x, sim_idx, frame);
                x_sim = x_sim_map_set(y, x, sim_idx, frame);
                frame_sim = frame_sim_map_set(y, x, sim_idx, frame);

                FLUX_noisy_cubelet1D = FLUX_noisy_map_set1D(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim);


                % 3D FFT of cubelet
                FLUX_noisy_cubelet_set3D(:, :, :, sim_idx) = fft2(FLUX_noisy_cubelet1D);

            end


            % 4D FFT of similar cubelet set (4D photon transient set)
            FLUX_noisy_cubelet_set4D = fft(FLUX_noisy_cubelet_set3D, N_sim, 4);        



            %% Signal & Noise separation
            
            % signal components
            SIG_cubelet_set = FLUX_noisy_cubelet_set4D(:, :, 1:N_sig_f, :);


            % noise components
            NOISE_cubelet_set = FLUX_noisy_cubelet_set4D(:, :, N_sig_f+1:end, :);



            %% Noise statistics
            
            % mean noise & mean squared noise
            ME = mean(abs(NOISE_cubelet_set(:)));
            MSE = mean(abs(NOISE_cubelet_set(:)).^2);

            MSE_map_set(y, x, frame) = MSE;


            % noise threshold       
            noise_scale = 1 + 5*sqrt((gamma(1)/gamma(1.5))^2 - 1);
            noise_th = noise_scale*ME;



            %% Adaptive guided photon processing
            
            % noise-to-signal ratio (NSR)
            N_power = MSE;
            S_power = mean(abs(SIG_cubelet_set(:)).^2);               % pseudo signal power before hard-thresholding


            % adaptive filtering
            if N_power/S_power >= R_th

                GP_mask_set(y, x, frame) = 1;


                % intensity patch
                i_patch = i_map_set(y : y+(C_size-1), x : x+(C_size-1), frame);
                i_patch = i_patch/sum(i_patch(:));
                I_patch = abs(fftn(i_patch));


                % guided photon processing
                SIG_cubelet_set_hat = SIG_cubelet_set.*repmat(I_patch, 1, 1, N_sig_f, N_sim);

            else
                
                
                % thresholding
                SIG_cubelet_set_hat = shrink_hard(SIG_cubelet_set, noise_th);
            end        



            %% Inverse 4D FFT
            FLUX_cubelet_set = zeros(C_size, C_size, N_t, N_sim);
            FLUX_cubelet_set(:, :, 1:N_sig_f, :) = SIG_cubelet_set_hat;               
            flux_cubelet_set = ifftn(FLUX_cubelet_set, 'symmetric');



            %% Aggregate
            weight = 1/(MSE + eps);

            for sim_idx = 1 : N_sim

                y_sim = y_sim_map_set(y, x, sim_idx, frame);
                x_sim = x_sim_map_set(y, x, sim_idx, frame);
                frame_sim = frame_sim_map_set(y, x, sim_idx, frame);

                weight_sum_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, frame_sim) = weight_sum_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, frame_sim) + weight;
                flux_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim) = flux_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim) + weight*flux_cubelet_set(:,:,:, sim_idx);

            end
        end
    end
end



% filtered flux by hard thresholding
weight_sum_map_set = reshape(weight_sum_map_set, N_y, N_x, 1, N_c);
flux_map_set = flux_map_set./repmat(weight_sum_map_set, 1, 1, N_t, 1);
flux_map_set(flux_map_set < 0) = 0;
