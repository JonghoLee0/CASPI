function compute_metric_set(d_map_set, d_hat_map_set, mask_set)

N_frame = size(d_map_set, 3);

fprintf('Error metric:\n')
for frame = 1 : N_frame
    
    fprintf('frame %d:\n', frame);
    
    d_map = d_map_set(:, :, frame);
    d_hat_map = d_hat_map_set(:, :, frame);
    mask = mask_set(:, :, frame);
    
    compute_metric(d_map(mask), d_hat_map(mask));
    
end