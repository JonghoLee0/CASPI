function [t_map, t_map_color, t_i_map_color, mask_hat] = make_t_map(flux_noisy_map, pseudo_i_map, mask, N_sample, t_bin, cmap, t_min, t_max, bin_size)



%% parameters
size_y = size(flux_noisy_map, 1);
size_x = size(flux_noisy_map, 2);
N_bin = size(flux_noisy_map, 3);

t_grid = [1 : N_bin]'*t_bin;
t_map = zeros(size_y, size_x);
mask_hat = mask;



%% find fitting domain
% find a location of the maximum flux
flux_noisy_sum = sum(flux_noisy_map, 3);
[max_val, max_idx] = max(flux_noisy_sum(:));
[y_idx, x_idx] = ind2sub([size_y, size_x], max_idx);


% extract a 11 x 11 histogram patch
y_start = max(y_idx-5, 1);
y_end = min(y_idx+5, size_y);
x_start = max(x_idx-5, 1);
x_end = min(x_idx+5, size_x);

flux_noisy_patch =flux_noisy_map(y_start:y_end, x_start:x_end, :);


% set a fitting domain from the initial 1/3 samples of the averaged flux
flux_noisy_avg = mean(flux_noisy_patch, 1);
flux_noisy_avg = mean(flux_noisy_avg, 2);
flux_noisy_avg = squeeze(flux_noisy_avg);
[max_val, start_idx] = max(flux_noisy_avg(1:round(N_bin/2)));
end_idx = min(start_idx + N_sample - 1, N_bin);



%% Loop
for y = 1 : size_y
    for x = 1 : size_y
        
        if mask(y, x)
            
            y_start = max(y-bin_size, 1);
            y_end = min(y+bin_size, size_y);
            x_start = max(x-bin_size, 1);
            x_end = min(x+bin_size, size_x);
            
            flux_noisy = flux_noisy_map(y_start:y_end, x_start:x_end, :);
            flux_noisy = sum(flux_noisy, 1);
            flux_noisy = sum(flux_noisy, 2);
            flux_noisy = squeeze(flux_noisy);

            time = t_grid(start_idx : end_idx);
            data = flux_noisy(start_idx : end_idx);
            

            % extract non-zero elements
            nonZeroIdx = find(data>0);
            
            if numel(nonZeroIdx) > 1
                time = time(nonZeroIdx);
                data = data(nonZeroIdx);
                
                
                % linear fitting with log data
                data = log(data);
                p = polyfit(time, data, 1);


                % save lifetime
%                 t_hat = -1/p(1)-0.07e-9;  % calibration for brain
                t_hat = -1/p(1);
                t_map(y, x) = max(t_hat, 0);
                
            else
                
                mask_hat(y, x) = 0;
                
            end
        end
    end
end


% color-code lifetime map
t_map_color = color_code(t_map, cmap, t_min, t_max);
t_i_map_color = t_map_color.*repmat(pseudo_i_map.*mask, 1, 1, 3);