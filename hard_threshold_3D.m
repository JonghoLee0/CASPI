function [flux_map_ht3D, MSE_map, IBFmask3D] = hard_threshold_3D(FLUX_map_noisy1D, i_map_1D, IBF, IBF_th, size_y, size_x, N_bin, s_patch, N_sig_f, skp, octave)



%% buffers
MSE_map = nan(size_y, size_x);             % map of mean squared noise
weight_sum_map = zeros(size_y, size_x);                     % map of weight sum
flux_map_ht3D = zeros(size_y, size_x, N_bin);               % recovered flux map
IBFmask3D = zeros(size_y, size_x);


if octave
    
    [y_grid, x_grid, t_grid] = ndgrid([1:s_patch], [1:s_patch], [1:N_sig_f]);

    y_grid = y_grid(:);
    x_grid = x_grid(:);
    t_grid = t_grid(:);

    y_grid2 = mod(s_patch - (y_grid-1), s_patch) + 1;
    x_grid2 = mod(s_patch - (x_grid-1), s_patch) + 1;
    t_grid2 = mod(N_bin - (t_grid-1), N_bin) + 1;
    
end



%%
% highSNR:(y,x)=(58,70), medimumSNR:(138,99), lowSNR:(28,103)
for y = 1 : skp : size_y - (s_patch-1)
    for x = 1 : skp : size_x - (s_patch-1)

        
        
        %% Extract a patch
        FLUX_patch_noisy1D = FLUX_map_noisy1D(y : y+(s_patch-1), x : x+(s_patch-1), :);
        
        
             
        %% 3D FFT
        FLUX_patch_noisy3D = fft2(FLUX_patch_noisy1D);
        


        %% Signal & Noise separation
        % signal
        SIG_patch_noisy3D = FLUX_patch_noisy3D(:, :, 1:N_sig_f);
        
        
        % noise
        NOISE_patch_noisy3D = FLUX_patch_noisy3D(:, :, N_sig_f+1 : end );
        
        
                  
        %% Noise statistics
        % mean noise & mean squared noise
        ME = mean(abs(NOISE_patch_noisy3D(:)));                  % for noise threshold computation
        MSE_map(y, x) = mean(abs(NOISE_patch_noisy3D(:)).^2);    % for weight computation and Wiener coefficient computation
        

        % noise threshold
        noise_scale = 1 + 3*sqrt((gamma(1)/gamma(1.5))^2 - 1);
        noise_th = noise_scale*ME;
        
        
%         % 
%         figure; hold on; grid on;
%         plot(abs(FLUX_patch_noisy3D(:)), 'color', [102 153 255]/255, 'linewidth', 6);
%         yline(noise_th, 'color', [255 0 0]/255, 'lineWidth', 6);
%         axis tight
%         ylim([0 0.4])
%         xticklabels(''); yticklabels('');

        
        
                             
        %% Intensity-based filtering or Hard thresholding
        if IBF
            
            % noise-to-signal ratio
            N_power = MSE_map(y, x);
            S_power = mean(abs(SIG_patch_noisy3D(:)).^2);                           % pseudo signal power before hard-thresholding

            
            if N_power/S_power >= IBF_th

                IBFmask3D(y, x) = 1;


                % intensity patch
                i_patch = i_map_1D(y : y+(s_patch-1), x : x+(s_patch-1));
                i_patch = i_patch/sum(i_patch(:));
                I_patch = abs(fftn(i_patch));

         
                % intensity-based filtering
                SIG_patch_ht3D = SIG_patch_noisy3D.*repmat(I_patch, 1, 1, N_sig_f);
                
            else
                % Hard thresholding
                [SIG_patch_ht3D, coeff_hard] = shrink_hard_sig(SIG_patch_noisy3D, noise_th);
            end
        else 
            % Hard thresholding
            [SIG_patch_ht3D, coeff_hard] = shrink_hard_sig(SIG_patch_noisy3D, noise_th);   
        end

        
              
        %% Inverse 3D FFT
        FLUX_patch = zeros(s_patch, s_patch, N_bin);
        FLUX_patch(:, :, 1:N_sig_f) = SIG_patch_ht3D;
        
        
        if octave
            
            FLUX_patch(sub2ind([s_patch, s_patch, N_bin], y_grid2, x_grid2, t_grid2)) =  conj(FLUX_patch(sub2ind([s_patch, s_patch, N_bin], y_grid, x_grid, t_grid)));
            flux_patch = ifftn(FLUX_patch);
            
        else % matlab
        
            flux_patch = ifftn(FLUX_patch, 'symmetric');
            
        end



        %% Aggregate
        % compute weight
        weight = 1/(MSE_map(y, x) + eps);
        weight_sum_map(y : y+(s_patch-1), x : x+(s_patch-1)) = weight_sum_map(y : y+(s_patch-1), x : x+(s_patch-1)) + weight;
        
        
        % weighted sum
        flux_map_ht3D(y : y+(s_patch-1), x : x+(s_patch-1), :) = flux_map_ht3D(y : y+(s_patch-1), x : x+(s_patch-1), :) + weight*flux_patch;
        
     
    end
end  


% filtered flux by hard thresholding
flux_map_ht3D = flux_map_ht3D./weight_sum_map;
flux_map_ht3D(flux_map_ht3D < 0) = 0;