% Get depth map from flux map by matched filtering


function [d_map, d_map_color] = get_d_map_MF(flux_map, mask, psf, c, delta_t, cmap, d_min, d_max)
    
bin_t = size(flux_map, 3);


% normalize flux
flux_map_sum = sum(flux_map, 3);
flux_map_norm = flux_map./repmat(flux_map_sum , 1, 1, bin_t);


% get correlation
correl = imfilter(flux_map_norm, reshape(psf, 1, 1, size(psf, 1)), 'corr');


% get index for max correlation
[dummy, tof_idx_map] = max(correl, [], 3);


% compute depth
tof_hat_map = tof_idx_map*delta_t;
d_map = c*tof_hat_map/2;
d_map(~mask) = nan;
d_map_color = color_code(d_map, cmap, d_min, d_max);