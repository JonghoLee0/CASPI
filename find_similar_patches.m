% Find similar patches


function [x_sim_map_set, y_sim_map_set, frame_sim_map_set] = find_similar_patches(pseudo_i_map_set, N_x, N_y, N_c, N_sim, C_size, R_intra, R_inter)


%% Buffers
x_sim_map_set = zeros(N_y, N_x, N_sim, N_c);
y_sim_map_set = zeros(N_y, N_x, N_sim, N_c);
frame_sim_map_set = zeros(N_y, N_x, N_sim, N_c);

S_intra = 2*R_intra + 1;  % search size within cube
S_inter = 2*R_inter + 1;  % search size over cubes



%% Find similar patches
for frame = 1 : N_c
    
    disp(['frame ', num2str(frame), ':']);
    
    for y = 1 : N_y - (C_size-1)
        for x = 1 : N_x - (C_size-1)


            % Reference patch
            patch = pseudo_i_map_set(y : y + (C_size-1), x : x + (C_size-1), frame);
            R = patch(:);


            % set search region
            frame_start = max(frame - R_inter, 1);
            frame_end = min(frame + R_inter, N_c);

            y_start = max(y - R_intra, 1);
            y_end = min(y + R_intra, N_y - (C_size-1));

            x_start = max(x - R_intra, 1);
            x_end = min(x + R_intra, N_x - (C_size-1));
            
            
            % buffers
            x_buffer = nan(S_intra*S_intra*S_inter, 1);
            y_buffer = nan(S_intra*S_intra*S_inter, 1);
            frame_buffer = nan(S_intra*S_intra*S_inter, 1);
            dist_buffer = nan(S_intra*S_intra*S_inter, 1);

            
            % compute patch distance over search volume
            for frame_srch = frame_start : frame_end
                for y_srch = y_start : y_end
                    for x_srch = x_start : x_end


                        % target patch
                        patch = pseudo_i_map_set(y_srch : y_srch+(C_size-1), x_srch : x_srch+(C_size-1), frame_srch);
                        T = patch(:);


                        % distance
                        dist = sum((R - T).^2);


                        % save
                        ind_srch = sub2ind([S_intra, S_intra, S_inter], (R_intra+1)+(y_srch-y), (R_intra+1)+(x_srch-x), (R_inter+1)+(frame_srch-frame));
                        
                        dist_buffer(ind_srch, 1) = dist;
                        x_buffer(ind_srch, 1) = x_srch;
                        y_buffer(ind_srch, 1) = y_srch;
                        frame_buffer(ind_srch, 1) = frame_srch;


                    end
                end
            end
            
            
            % sort based on spatial distance
            xyDist_buffer = (x - x_buffer).^2 + (y - y_buffer).^2;
            [sorted_xyDist, sorted_xyDistIdx] = sort(xyDist_buffer, 'ascend');

            x_buffer = x_buffer(sorted_xyDistIdx);
            y_buffer = y_buffer(sorted_xyDistIdx);
            frame_buffer = frame_buffer(sorted_xyDistIdx);
            dist_buffer = dist_buffer(sorted_xyDistIdx);


            % sort based on frame distance
            frameDist_buffer = (frame - frame_buffer).^2;
            [sorted_frameDist, sorted_frameDistIdx] = sort(frameDist_buffer, 'ascend');

            x_buffer = x_buffer(sorted_frameDistIdx);
            y_buffer = y_buffer(sorted_frameDistIdx);
            frame_buffer = frame_buffer(sorted_frameDistIdx);
            dist_buffer = dist_buffer(sorted_frameDistIdx);


            % sort based on patch distance
            [sorted_dist, sorted_idx] = sort(dist_buffer, 'ascend');
            sorted_idx = sorted_idx(1:N_sim);

            x_sim_map_set(y, x, :, frame) = x_buffer(sorted_idx);
            y_sim_map_set(y, x, :, frame) = y_buffer(sorted_idx);
            frame_sim_map_set(y, x, :, frame) = frame_buffer(sorted_idx);

        end
    end
end
