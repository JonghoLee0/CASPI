function [x_sim_map_set, y_sim_map_set, frame_sim_map_set] = find_similar_patches(pseudo_i_map_set, size_x, size_y, N_frame, N_sim, s_patch, r_intra, r_inter)


% buffers for similar patches
x_sim_map_set = zeros(size_y, size_x, N_sim, N_frame);
y_sim_map_set = zeros(size_y, size_x, N_sim, N_frame);
frame_sim_map_set = zeros(size_y, size_x, N_sim, N_frame);

s_intra = 2*r_intra + 1;  % intra-frame search size
s_inter = 2*r_inter + 1;  % inter-frame search size



%% Loop
for frame = 1 : N_frame
    
    fprintf('frame %d:\n', frame);
    
    for y = 1 : size_y - (s_patch-1)
        for x = 1 : size_x - (s_patch-1)



            %% Reference patch
            patch = pseudo_i_map_set(y : y + (s_patch-1), x : x + (s_patch-1), frame);
            R = patch(:);



            %% Target patches
            % set search region
            frame_start = max(frame - r_inter, 1);
            frame_end = min(frame + r_inter, N_frame);

            y_start = max(y - r_intra, 1);
            y_end = min(y + r_intra, size_y - (s_patch-1));

            x_start = max(x - r_intra, 1);
            x_end = min(x + r_intra, size_x - (s_patch-1));
            
            
            % buffers for temporary positions
            x_buffer = nan(s_intra*s_intra*s_inter, 1);      % absolute x position
            y_buffer = nan(s_intra*s_intra*s_inter, 1);      % absolute y position
            frame_buffer = nan(s_intra*s_intra*s_inter, 1);  % frame position
            dist_buffer = nan(s_intra*s_intra*s_inter, 1);   % distance between patches

            for frame_srch = frame_start : frame_end
                for y_srch = y_start : y_end
                    for x_srch = x_start : x_end


                        % target patch
                        patch = pseudo_i_map_set(y_srch : y_srch+(s_patch-1), x_srch : x_srch+(s_patch-1), frame_srch);
                        T = patch(:);


                        % distance
                        dist = sum((R - T).^2);


                        % save
                        ind_srch = sub2ind([s_intra, s_intra, s_inter], (r_intra+1)+(y_srch-y), (r_intra+1)+(x_srch-x), (r_inter+1)+(frame_srch-frame)); % linear search index for current position

                        dist_buffer(ind_srch, 1) = dist;
                        x_buffer(ind_srch, 1) = x_srch;
                        y_buffer(ind_srch, 1) = y_srch;
                        frame_buffer(ind_srch, 1) = frame_srch;


                    end
                end
            end


            
            %% find simiar patches
            
            
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
