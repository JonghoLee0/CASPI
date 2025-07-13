% Get depth maps from photon fluxes by matched filtering


function [d_map_set, d_color_map_set] = get_d_map_set(flux_map_set, psf, c, delta_t, cmap, d_min, d_max)

N_c = size(flux_map_set, 4);

d_map_set = cell(N_c, 1);
d_color_map_set = cell(N_c, 1);

for frame = 1 : N_c
    
    flux_map = flux_map_set(:, :, :, frame);
    mask = ~isnan(sum(flux_map, 3));
    
    
    % get depth map by matched filtering
    [d_map, d_map_color] = get_d_map_MF(flux_map, mask, psf, c, delta_t, cmap, d_min, d_max);
    
    d_map_set{frame} = d_map;
    d_color_map_set{frame} = d_map_color;
    
end