% Wiener filtering using local and non-local correlations


function flux_map_set_wiener = estimate_wiener_NLC(flux_map_set_init_NLC, FLUX_noisy_map_set1D, x_sim_map_set, y_sim_map_set, frame_sim_map_set, MSE_map_set, N_y, N_x, N_t, N_t_half, N_sim, N_c, C_size, N_sig_f, skp)



%% Buffers
flux_map_set_wiener = zeros(N_y, N_x, N_t, N_c);
weight_sum_map_set = zeros(N_y, N_x, N_c);



%% 1D FFT of hard-thresholded
FLUX_map_set_init1D = fft(flux_map_set_init_NLC, N_t, 3);
FLUX_map_set_init1D = FLUX_map_set_init1D(:,:, 1:N_t_half, :);



%% Wiener filtering using local and non-local correlations
for frame = 1 : N_c

	disp(['frame ', num2str(frame), ':']);

	for y = 1 : skp : N_y-(C_size-1)
		for x = 1 : skp :  N_x-(C_size-1)         
			
			
			  
			%% Collect similar patches
			FLUX_cubelet_noisy3D = zeros(C_size, C_size, N_t_half, N_sim);
			FLUX_cubelet_init3D = zeros(C_size, C_size, N_t_half, N_sim);
	   
			for sim_idx = 1 : N_sim

				y_sim = y_sim_map_set(y, x, sim_idx, frame);
				x_sim = x_sim_map_set(y, x, sim_idx, frame);
				frame_sim = frame_sim_map_set(y, x, sim_idx, frame);
				
				
				% 3D FFT of initial flux
				FLUX_cubelet_init1D = FLUX_map_set_init1D(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim);
				FLUX_cubelet_init3D(:, :, :, sim_idx) = fft2(FLUX_cubelet_init1D);
				
				
				% 3D FFT of noisy flux
				FLUX_cubelet_noisy1D = FLUX_noisy_map_set1D(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim);
				FLUX_cubelet_noisy3D(:, :, :, sim_idx) = fft2(FLUX_cubelet_noisy1D);
				
			end
			
			
				 
			% 4D FFT of initial flux
			FLUX_cubelet_init4D = fft(FLUX_cubelet_init3D, N_sim, 4);
			SIG_cubelet_init4D = FLUX_cubelet_init4D(:,:,1:N_sig_f,:);
			
			
			% 4D FFT of noisy flux
			FLUX_cubelet_noisy4D = fft(FLUX_cubelet_noisy3D, N_sim, 4);
			SIG_cubelet_noisy4D = FLUX_cubelet_noisy4D(:,:,1:N_sig_f,:);
			
			
		
			%% Wiener shrinkage
			% noise-to-signal ratio
			N_power = MSE_map_set(y, x, frame);
			S_power = abs(SIG_cubelet_init4D).^2;  
			S_power(S_power == 0) = eps;

			
			% Wiener coefficient
			wiener_coeff = 1./(1 + N_power./S_power);
			
	  
			% Wiener shrinkage
			SIG_patch_wiener4D = wiener_coeff.*SIG_cubelet_noisy4D;
			
			
			
			%% Inverse 4D FFT
			FLUX_wiener_patchset = zeros(C_size, C_size, N_t, N_sim);
			FLUX_wiener_patchset(:, :, 1:N_sig_f, :) = SIG_patch_wiener4D;  
            flux_wiener_patchset = ifftn(FLUX_wiener_patchset, 'symmetric');
			
	  
			
			%% Aggregate
			weight = 1/(N_power + eps);
			
			for sim_idx = 1 : N_sim

				y_sim = y_sim_map_set(y, x, sim_idx, frame);
				x_sim = x_sim_map_set(y, x, sim_idx, frame);
				frame_sim = frame_sim_map_set(y, x, sim_idx, frame);	
				
				weight_sum_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, frame_sim) = weight_sum_map_set(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, frame_sim) + weight;
				flux_map_set_wiener(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim) = flux_map_set_wiener(y_sim : y_sim+C_size-1, x_sim : x_sim+C_size-1, :, frame_sim) + weight*flux_wiener_patchset(:,:,:,sim_idx);
				
			end 
		end
	end
end


% refined flux using local and non-local correlations
weight_sum_map_set = reshape(weight_sum_map_set, N_y, N_x, 1, N_c);
flux_map_set_wiener = flux_map_set_wiener./repmat(weight_sum_map_set, 1, 1, N_t, 1);
flux_map_set_wiener(flux_map_set_wiener < 0) = 0;