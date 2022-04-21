function flux_map_wiener3D = wiener_filter_3D(flux_map_ht3D, FLUX_map_noisy1D, MSE_map3D, size_y, size_x, N_bin, N_bin_half, s_patch, N_sig_f, skp, octave)



%% buffers
flux_map_wiener3D = zeros(size_y, size_x, N_bin);
weight_sum_map = zeros(size_y, size_x);


if octave
    
    [y_grid, x_grid, t_grid] = ndgrid([1:s_patch], [1:s_patch], [1:N_sig_f]);

    y_grid = y_grid(:);
    x_grid = x_grid(:);
    t_grid = t_grid(:);

    y_grid2 = mod(s_patch - (y_grid-1), s_patch) + 1;
    x_grid2 = mod(s_patch - (x_grid-1), s_patch) + 1;
    t_grid2 = mod(N_bin - (t_grid-1), N_bin) + 1;

end



%% 1D FFT of hard-thresholded
FLUX_map_ht1D = fft(flux_map_ht3D, N_bin, 3);
FLUX_map_ht1D = FLUX_map_ht1D(:,:, 1:N_bin_half);



%%
for y = 1 : skp : size_y-(s_patch-1)
    for x = 1 : skp :  size_x-(s_patch-1)    
        
        
    
        %% 3D FFT
        % hard-thresholded
        FLUX_patch_ht1D = FLUX_map_ht1D(y : y+(s_patch-1), x : x+(s_patch-1), :);
        FLUX_patch_ht3D = fft2(FLUX_patch_ht1D);        
        SIG_patch_ht3D = FLUX_patch_ht3D(:, :, 1:N_sig_f);

        
        % noisy
        FLUX_patch_noisy1D = FLUX_map_noisy1D(y : y+(s_patch-1), x : x+(s_patch-1), :);
        FLUX_patch_noisy3D = fft2(FLUX_patch_noisy1D); 
        SIG_patch_noisy3D = FLUX_patch_noisy3D(:, :, 1:N_sig_f);
        
              
              
        %% Wiener shrinkage
        % noise-to-signal ratio
        N_power = MSE_map3D(y, x);
        S_power = abs(SIG_patch_ht3D).^2;  
        S_power(S_power == 0) = eps;

        
        % Wiener coefficient
        wiener_coeff = 1./(1 + N_power./S_power);
        
  
        % Wiener shrinkage
        SIG_patch_wiener3D = wiener_coeff.*SIG_patch_noisy3D;
        
        
        
        %% Inverse 3D FFT
        FLUX_patch_wiener3D = zeros(s_patch, s_patch, N_bin);
        FLUX_patch_wiener3D(:, :, 1:N_sig_f) = SIG_patch_wiener3D;
        
        
        if octave
            
            FLUX_patch_wiener3D(sub2ind([s_patch, s_patch, N_bin], y_grid2, x_grid2, t_grid2)) =  conj(FLUX_patch_wiener3D(sub2ind([s_patch, s_patch, N_bin], y_grid, x_grid, t_grid)));
            flux_patch_wiener3D = ifftn(FLUX_patch_wiener3D);
            
        else % matlab
            
            flux_patch_wiener3D = ifftn(FLUX_patch_wiener3D, 'symmetric');
        
        end
  
        
        
        %% Aggregate
        % compute weight
        weight = 1/(N_power + eps);
        weight_sum_map(y : y+(s_patch-1), x : x+(s_patch-1)) = weight_sum_map(y : y+(s_patch-1), x : x+(s_patch-1)) + weight;
        
        
        % weighted sum
        flux_map_wiener3D(y : y+(s_patch-1), x : x+(s_patch-1), :) = flux_map_wiener3D(y : y+(s_patch-1), x : x+(s_patch-1), :) + weight*flux_patch_wiener3D;


        
    end
end


% weighted-average flux
flux_map_wiener3D = flux_map_wiener3D./weight_sum_map;
flux_map_wiener3D(flux_map_wiener3D < 0) = 0;