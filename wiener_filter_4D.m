function flux_map_set_wiener4D = wiener_filter_4D(flux_map_set_ht4D, FLUX_map_set_noisy1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, MSE_map_set, size_y, size_x, N_bin, N_bin_half, N_sim, N_frame, s_patch, N_sig_f, skp, octave)



%% Buffers
flux_map_set_wiener4D = zeros(size_y, size_x, N_bin, N_frame);
weight_sum_map_set = zeros(size_y, size_x, N_frame);


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



%% 1D FFT of hard-thresholded
FLUX_map_set_ht1D = fft(flux_map_set_ht4D, N_bin, 3);
FLUX_map_set_ht1D = FLUX_map_set_ht1D(:,:, 1:N_bin_half, :);



%% Loop
for frame = 1 : N_frame

	fprintf('frame %d:\n', frame);

	for y = 1 : skp : size_y-(s_patch-1)
		for x = 1 : skp :  size_x-(s_patch-1)         
			
			
			  
			%% Collect similar patches
			FLUX_patchset_noisy3D = zeros(s_patch, s_patch, N_bin_half, N_sim);
			FLUX_patchset_ht3D = zeros(s_patch, s_patch, N_bin_half, N_sim);
	   
			for sim_idx = 1 : N_sim

				y_sim = y_sim_map_set(y, x, sim_idx, frame);
				x_sim = x_sim_map_set(y, x, sim_idx, frame);
				frame_sim = frame_sim_map_set(y, x, sim_idx, frame);
				
				
				% 3D FFT of hard-thresholded
				FLUX_patch_ht1D = FLUX_map_set_ht1D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim);
				FLUX_patchset_ht3D(:, :, :, sim_idx) = fft2(FLUX_patch_ht1D);
				
				
				% 3D FFT of noisy
				FLUX_patch_noisy1D = FLUX_map_set_noisy1D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim);
				FLUX_patchset_noisy3D(:, :, :, sim_idx) = fft2(FLUX_patch_noisy1D);
				
				
			end
			
			
				 
			%% 4D FFT
			% hard-thresholded
			FLUX_patchset_ht4D = fft(FLUX_patchset_ht3D, N_sim, 4);
			SIG_patchset_ht4D = FLUX_patchset_ht4D(:,:,1:N_sig_f,:);
			
			
			% noisy
			FLUX_patchset_noisy4D = fft(FLUX_patchset_noisy3D, N_sim, 4);
			SIG_patchset_noisy4D = FLUX_patchset_noisy4D(:,:,1:N_sig_f,:);
			
			
		
			%% Wiener shrinkage
			% noise-to-signal ratio
			N_power = MSE_map_set(y, x, frame);
			S_power = abs(SIG_patchset_ht4D).^2;  
			S_power(S_power == 0) = eps;

			
			% Wiener coefficient
			wiener_coeff = 1./(1 + N_power./S_power);
			
	  
			% Wiener shrinkage
			SIG_patch_wiener4D = wiener_coeff.*SIG_patchset_noisy4D;
			
			
			
			%% Inverse 4D FFT
			FLUX_wiener_patchset = zeros(s_patch, s_patch, N_bin, N_sim);
			FLUX_wiener_patchset(:, :, 1:N_sig_f, :) = SIG_patch_wiener4D;
            
            if octave
                
                FLUX_wiener_patchset(:, :, :, round((N_sim-1)/2)+2 : end) = 0;
                FLUX_wiener_patchset(sub2ind([s_patch, s_patch, N_bin, N_sim], y_grid2, x_grid2, t_grid2, s_grid2)) =  conj(FLUX_wiener_patchset(sub2ind([s_patch, s_patch, N_bin, N_sim], y_grid, x_grid, t_grid, s_grid)));
                flux_wiener_patchset = ifftn(FLUX_wiener_patchset);
                
            else % matlab
                
                flux_wiener_patchset = ifftn(FLUX_wiener_patchset, 'symmetric');
                
            end
			
	  
			
			%% Aggregate
			weight = 1/(N_power + eps);
			
			for sim_idx = 1 : N_sim

				y_sim = y_sim_map_set(y, x, sim_idx, frame);
				x_sim = x_sim_map_set(y, x, sim_idx, frame);
				frame_sim = frame_sim_map_set(y, x, sim_idx, frame);	
				
				weight_sum_map_set(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, frame_sim) = weight_sum_map_set(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, frame_sim) + weight;
				flux_map_set_wiener4D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim) = flux_map_set_wiener4D(y_sim : y_sim+s_patch-1, x_sim : x_sim+s_patch-1, :, frame_sim) + weight*flux_wiener_patchset(:,:,:,sim_idx);
				
			end 
		end
	end
end


% weighted-average flux
weight_sum_map_set = reshape(weight_sum_map_set, size_y, size_x, 1, N_frame);
flux_map_set_wiener4D = flux_map_set_wiener4D./repmat(weight_sum_map_set, 1, 1, N_bin, 1);
flux_map_set_wiener4D(flux_map_set_wiener4D < 0) = 0;