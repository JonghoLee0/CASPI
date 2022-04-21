function d_hat_map_set = get_d_hat_map_set(flux_map_set, mask_set, psf, c, t_bin, cmap, visual_d_min, visual_d_max, method, app, name)

size_y = size(flux_map_set, 1);
size_x = size(flux_map_set, 2);
N_frame = size(flux_map_set, 4);
d_hat_map_set = zeros(size_y, size_x, N_frame);

for frame = 1 : N_frame
    
    flux_map = flux_map_set(:, :, :, frame);
    mask = mask_set(:, :, frame);
    
    
    % depth map 
    [d_hat_map, d_hat_map_color] = make_d_map(flux_map, mask, psf, c, t_bin, cmap, visual_d_min, visual_d_max, method);
    
    
    % save monochrome depth
    d_hat_map_set(:, :, frame) = d_hat_map;
    
    
    % write color depth
    file_name = [app, '\Results\d_color_', num2str(frame), '_', name, '.png'];
    imwrite(d_hat_map_color, file_name);

    
end