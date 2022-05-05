function i_map_adjust = adjust_intensity(src_i_map, tgt_i_map)


src_mean = mean(src_i_map(:));
src_std = std(src_i_map(:));

tgt_mean = mean(tgt_i_map(:));
tgt_std = std(tgt_i_map(:));


i_map_adjust = (src_i_map - src_mean)/src_std*tgt_std + tgt_mean;
i_map_adjust(i_map_adjust < 0) = 0;
i_map_adjust(i_map_adjust > 1) = 1;


