% data: M x N array
% cmap: color map
% data_colored: M x N array


function data_colored = color_code(data, cmap, data_min, data_max)

y_size = size(data, 1);
x_size = size(data, 2);

data = data(:);
data(data<data_min) = data_min;
data(data>data_max) = data_max;

x = [data_min : (data_max - data_min)/(size(cmap,1) - 1) : data_max]';

R = cmap(:, 1);
G = cmap(:, 2);
B = cmap(:, 3);

data_R = interp1(x, R, data);
data_G = interp1(x, G, data);
data_B = interp1(x, B, data);

data_colored = nan(y_size, x_size, 3);

tmp = reshape(data_R, y_size, x_size);
tmp(isnan(tmp)) = 1;
data_colored(:, :, 1) = tmp;

tmp = reshape(data_G, y_size, x_size);
tmp(isnan(tmp)) = 1;
data_colored(:, :, 2) = tmp;

tmp = reshape(data_B, y_size, x_size);
tmp(isnan(tmp)) = 1;
data_colored(:, :, 3) = tmp;