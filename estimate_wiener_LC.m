% Wiener filtering using local correlations


function flux_map_wiener = estimate_wiener_LC(flux_map_init, FLUX_noisy_map1D, MSE_map, N_y, N_x, N_t, N_t_half, C_size, N_sig_f, skp)



%% Buffers
flux_map_wiener = zeros(N_y, N_x, N_t);
weight_sum_map = zeros(N_y, N_x);



%% 1D FFT
FLUX_map_init1D = fft(flux_map_init, N_t, 3);
FLUX_map_init1D = FLUX_map_init1D(:, :, 1:N_t_half);



%% Wiener filtering using local correlations
for y = 1 : skp : N_y - (C_size-1)
    for x = 1 : skp :  N_x - (C_size-1)    
        
        
    
        %% 3D FFT
        % initial flux
        FLUX_cubelet_init1D = FLUX_map_init1D(y : y+(C_size-1), x : x+(C_size-1), :);
        FLUX_cubelet_init3D = fft2(FLUX_cubelet_init1D);        
        SIG_cubelet_init = FLUX_cubelet_init3D(:, :, 1:N_sig_f);

        
        % noisy flux
        FLUX_cubelet_noisy1D = FLUX_noisy_map1D(y : y+(C_size-1), x : x+(C_size-1), :);
        FLUX_cubelet_noisy3D = fft2(FLUX_cubelet_noisy1D); 
        SIG_cubelet_noisy = FLUX_cubelet_noisy3D(:, :, 1:N_sig_f);
        
              
              
        %% Wiener shrinkage
        % noise-to-signal ratio
        N_power = MSE_map(y, x);
        S_power = abs(SIG_cubelet_init).^2;  
        S_power(S_power == 0) = eps;

        
        % Wiener coefficient
        wiener_coeff = 1./(1 + N_power./S_power);
        
  
        % Wiener shrinkage
        SIG_cubelet_hat = wiener_coeff.*SIG_cubelet_noisy;
        
        
        
        %% Inverse 3D FFT
        FLUX_cubelet = zeros(C_size, C_size, N_t);
        FLUX_cubelet(:, :, 1:N_sig_f) = SIG_cubelet_hat;
        flux_patch_wiener3D = ifftn(FLUX_cubelet, 'symmetric');
          
        
        
        %% Aggregate
        % compute weight
        weight = 1/(N_power + eps);
        weight_sum_map(y : y+(C_size-1), x : x+(C_size-1)) = weight_sum_map(y : y+(C_size-1), x : x+(C_size-1)) + weight;
        
        
        % weighted sum
        flux_map_wiener(y : y+(C_size-1), x : x+(C_size-1), :) = flux_map_wiener(y : y+(C_size-1), x : x+(C_size-1), :) + weight*flux_patch_wiener3D;


        
    end
end


% refined flux using local correlations
flux_map_wiener = flux_map_wiener./weight_sum_map;
flux_map_wiener(flux_map_wiener < 0) = 0;