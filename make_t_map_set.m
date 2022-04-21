function [t_map_set, t_map_set_color, ti_map_set_color] = make_t_map_set(flux_map_set, i_map_set, mask_set, N_sample, t_bin, cmap, t_min, t_max)



%% parameters
size_y = size(flux_map_set, 1);
size_x = size(flux_map_set, 2);
N_frame = size(flux_map_set, 4);

t_map_set = nan(size_y, size_x, N_frame);
t_map_set_color = zeros(size_y, size_x, 3, N_frame);
ti_map_set_color = zeros(size_y, size_x, 3, N_frame);



%% Loop
for frame = 1 : N_frame
    
    flux_map = squeeze(flux_map_set(:, :, :, frame));
    i_map = squeeze(i_map_set(:, :, frame));
    mask = squeeze(mask_set(:, :, frame));

    [t_map, t_map_color, ti_map_color] = make_t_map(flux_map, i_map, mask, N_sample, t_bin, cmap, t_min, t_max);
    
    t_map_set(:, :, frame) = t_map;
    t_map_set_color(:, :, :, frame) = t_map_color;
    ti_map_set_color(:, :, :, frame) = ti_map_color;
    
end


