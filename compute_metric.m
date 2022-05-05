function compute_metric(d_map, d_hat_map)


d_error = abs(d_map - d_hat_map);


% 0.2% inlier
inlier_idx = d_error < d_map*0.002;
inlier_error = d_error(inlier_idx);
inlierRate1 = sum(inlier_idx(:))/numel(d_error);
RMSE1 = sqrt(mean(inlier_error.^2));


% 0.5% inlier
inlier_idx = d_error < d_map*0.005;
inlier_error = d_error(inlier_idx);
inlierRate2 = sum(inlier_idx(:))/numel(d_error);
RMSE2 = sqrt(mean(inlier_error.^2));


% 1% inlier
inlier_idx = d_error < d_map*0.01;
inlier_error = d_error(inlier_idx);
inlierRate3 = sum(inlier_idx(:))/numel(d_error);
RMSE3 = sqrt(mean(inlier_error.^2));


% % 5% inlier
% inlier_idx = d_error < d_map*0.05;
% inlier_error = d_error(inlier_idx);
% inlierRate4 = sum(inlier_idx(:))/numel(d_error);
% RMSE4 = sqrt(mean(inlier_error.^2));


% print
fprintf('%s%2.4f%s%2d%s\n', '0.2% RMSE: ', RMSE1, '(',  round(inlierRate1*100),'%)');
fprintf('%s%2.4f%s%2d%s\n', '0.5% RMSE: ', RMSE2, '(',  round(inlierRate2*100),'%)');
fprintf('%s%2.4f%s%2d%s\n\n', '1.0% RMSE: ', RMSE3, '(',  round(inlierRate3*100),'%)');
% fprintf('%s%2.4f%s%2d%s\n\n', '5.0% RMSE: ', RMSE4, '(',  round(inlierRate4*100),'%)');