% Estimate initial fluxes using local correlations


function [flux_map, MSE_map, GP_mask] = estimate_init_flux_LC(FLUX_noisy_map1D, i_map_1D, R_th, N_y, N_x, N_t, C_size, N_sig_f, skp)



%% Buffers
flux_map = zeros(N_y, N_x, N_t);                        % initial flux map
MSE_map = nan(N_y, N_x);                                % mean squared error map
GP_mask = zeros(N_y, N_x);                              % mask for guided photon processing
weight_sum_map = zeros(N_y, N_x);                       % weight sum map





%% Initial flux estimation using local correlations
for y = 1 : skp : N_y - (C_size-1)
    for x = 1 : skp : N_x - (C_size-1)

        
        
        %% Extract cubelet
        FLUX_noisy_cubelet1D = FLUX_noisy_map1D(y : y+(C_size-1), x : x+(C_size-1), :);
        
        
             
        %% 3D FFT of cubelet
        FLUX_noisy_cubelet = fft2(FLUX_noisy_cubelet1D);
        


        %% Signal & noise separation
        
        % signal components
        SIG_cubelet = FLUX_noisy_cubelet(:, :, 1:N_sig_f);
        
        
        % noise components
        NOISE_cubelet = FLUX_noisy_cubelet(:, :, N_sig_f+1 : end);
        
        
                  
        %% Noise statistics
        
        % mean noise & mean squared noise
        ME = mean(abs(NOISE_cubelet(:)));                  % for noise threshold computation
        MSE_map(y, x) = mean(abs(NOISE_cubelet(:)).^2);    % for weight computation and Wiener coefficient computation
        

        % noise threshold
        noise_scale = 1 + 3*sqrt((gamma(1)/gamma(1.5))^2 - 1);
        noise_th = noise_scale*ME;
        
        
                             
        %% Adaptive guided photon processing
        
        % noise-to-signal ratio (NSR)
        N_power = MSE_map(y, x);
        S_power = mean(abs(SIG_cubelet(:)).^2);
        
        
        % adaptive filtering
        if N_power/S_power >= R_th

            GP_mask(y, x) = 1;


            % intensity patch
            i_patch = i_map_1D(y : y+(C_size-1), x : x+(C_size-1));
            i_patch = i_patch/sum(i_patch(:));
            I_patch = abs(fftn(i_patch));


            % guided photon processing
            SIG_cubelet_hat = SIG_cubelet.*repmat(I_patch, 1, 1, N_sig_f);

        else
            
            
            % thresholding
            SIG_cubelet_hat = shrink_hard(SIG_cubelet, noise_th);
        end


                    
        %% Inverse 3D FFT
        FLUX_cubelet = zeros(C_size, C_size, N_t);
        FLUX_cubelet(:, :, 1:N_sig_f) = SIG_cubelet_hat;     
        flux_cubelet = ifftn(FLUX_cubelet, 'symmetric');
            


        %% Aggregate
        % compute weight
        weight = 1/(MSE_map(y, x) + eps);
        weight_sum_map(y : y+(C_size-1), x : x+(C_size-1)) = weight_sum_map(y : y+(C_size-1), x : x+(C_size-1)) + weight;
        
        
        % weighted sum
        flux_map(y : y+(C_size-1), x : x+(C_size-1), :) = flux_map(y : y+(C_size-1), x : x+(C_size-1), :) + weight*flux_cubelet;
        
     
        
    end
end  


% initial flux using local correlations
flux_map = flux_map./weight_sum_map;
flux_map(flux_map < 0) = 0;