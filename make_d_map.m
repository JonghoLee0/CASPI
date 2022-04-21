function [d_map, d_map_color] = make_d_map(flux_map, mask, psf, c, bin_t, cmap, d_min, d_max, method)


switch method
    
    
    % max
    case 1
        
        [dummy, tof_idx_map] = max(flux_map, [], 3);
        tof_hat_map = tof_idx_map*bin_t;
        d_map = c*tof_hat_map/2;
        d_map(~mask) = nan;
          
        if d_max < 20
            d_map_color = color_code(d_map, cmap, d_min, d_max);
        else
            d_map_color = color_code(log(d_map), cmap, log(d_min), log(d_max));
        end
        
        
    % matched filter
    case 2
        
        bin_N = size(flux_map, 3);
        
        flux_map_sum = sum(flux_map, 3);
        flux_map_norm = flux_map./repmat(flux_map_sum , 1, 1, bin_N);
        
        correl = imfilter(flux_map_norm, reshape(psf, 1, 1, size(psf, 1)), 'corr');
        
        [dummy, tof_idx_map] = max(correl, [], 3);

        tof_hat_map = tof_idx_map*bin_t;
        d_map = c*tof_hat_map/2;
        d_map(~mask) = nan;
        
        if d_max < 20
            d_map_color = color_code(d_map, cmap, d_min, d_max);
        else
            d_map_color = color_code(log(d_map), cmap, log(d_min), log(d_max));
        end
         
end