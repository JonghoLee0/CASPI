% Get lifetime images from photon fluxes by naive linear filtering to log-transformed data


function [t_map_set, t_color_map_set] = get_t_map_set(flux_map_set, int_prct, int_exp, N_sample, binning_size, delta_t, cmap, t_min, t_max)

N_c = size(flux_map_set, 4);

t_map_set = cell(N_c, 1);
t_color_map_set = cell(N_c, 1);

for frame = 1 : N_c
    
    flux_map = flux_map_set(:, :, :, frame);
    mask = ~isnan(sum(flux_map, 3));
    
    
    % intensity
    i_map = sum(flux_map, 3);
    i_map = i_map/prctile(i_map(:), int_prct);
    i_map = i_map.^int_exp;
    
    
    % get lifetime images by naive linear fitting to log-transformed data
    [t_map, t_color_map] = get_t_map_linear(flux_map, i_map, mask, N_sample, delta_t, cmap, t_min, t_max, binning_size);
    
    t_map_set{frame} = t_map;
    t_color_map_set{frame} = t_color_map;
    
end