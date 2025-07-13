% Get a lifetime image by naive linear fitting to log-transformed photon data


function [t_map, t_map_color] = get_t_map_linear(flux_map, i_map, mask, N_sample, delta_t, cmap, t_min, t_max, binning_size)

N_y = size(flux_map, 1);
N_x = size(flux_map, 2);
N_t = size(flux_map, 3);

t_grid = [1 : N_t]'*delta_t;
t_map = zeros(N_y, N_x);


% position of max photon count
flux_map_sum = sum(flux_map, 3);
[max_val, max_idx] = max(flux_map_sum(:));
[y_highSNR, x_highSNR] = ind2sub([N_y, N_x], max_idx);


% extract a 11 x 11 cubelet
y_start = max(y_highSNR - 5, 1);
y_end = min(y_highSNR + 5, N_y);
x_start = max(x_highSNR - 5, 1);
x_end = min(x_highSNR + 5, N_x);

cubelet_highSNR = flux_map(y_start:y_end, x_start:x_end, :);


% set a fitting domain
flux_highSNR = mean(cubelet_highSNR, 1);
flux_highSNR = mean(flux_highSNR, 2);
flux_highSNR = squeeze(flux_highSNR);
[max_val, fit_start] = max(flux_highSNR(1:round(N_t/2)));
fit_end = min(fit_start+N_sample-1, N_t);



for y = 1 : N_y
    for x = 1 : N_y
        
        if mask(y, x)
            
            
            % spatial binning if necessary
            y_start = max(y-binning_size, 1);
            y_end = min(y+binning_size, N_y);
            x_start = max(x-binning_size, 1);
            x_end = min(x+binning_size, N_x);
            
            flux = flux_map(y_start:y_end, x_start:x_end, :);
            flux = sum(flux, 1);
            flux = sum(flux, 2);
            flux = squeeze(flux);

            time = t_grid(fit_start : fit_end);
            data = flux(fit_start : fit_end);
            

            % extract non-zero elements
            nonzero_idx = find(data>0);
            
            
            if numel(nonzero_idx) > 1
                
                time = time(nonzero_idx);
                data = data(nonzero_idx);
                
                
                % naive linear fitting to log-transformed photon data
                data = log(data);
                p = polyfit(time, data, 1);


                % save lifetime
                t_hat = -1/p(1);
                
                if t_hat < t_min
                    t_hat = t_min;
                end
                
                if t_hat > t_max
                    t_hat = t_max;
                end
                
                t_map(y, x) = t_hat;
                
            end
        end
    end
end

t_map_color = color_code(t_map, cmap, t_min, t_max);
t_map_color = t_map_color.*repmat(i_map.*mask, 1, 1, 3);