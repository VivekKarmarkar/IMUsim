%% SET PREFERENCES
activity_type_walking = true;
load_uvm_data_bool = false;

%% LOAD DATA
if load_uvm_data_bool
    load('WalkingFP_RawData_reprocessed.mat');
end
load('marker_names_saved.mat');
load('body_segment_input_data.mat');

%% PROCESSING
num_activities = length(activity_reprocessed);

activity_start_idx = 2;
activity_end_idx = num_activities;
if activity_type_walking
    activity_start_idx = 5;
    activity_end_idx = 14;
end

side_list = ["right", "left"];
segment_list = strings(1,2);
for k=1:length(body_segment_input)
    segment_list(1,k) = body_segment_input(k).name;
end

technical_frame_markers_rigidBody_info = struct;
all_markers_info = struct;
num_segments = length(segment_list);
num_sides = length(side_list);

for j=activity_start_idx:activity_end_idx
    current_activity_idx = j;
    current_activity_name = strrep(activity_reprocessed(current_activity_idx).name, " ", "_");
    for s=1:25
        current_subject_idx = s;
        
        for segment_idx=1:num_segments
            current_segment = segment_list(segment_idx);
            
            for side_idx=1:num_sides
                current_side = side_list(side_idx);

                input_struct_val = struct;
                input_struct_val.marker_names_all_list_val = marker_names_all_list;
                input_struct_val.current_activity_idx_val = current_activity_idx;
                input_struct_val.current_activity_name_val = current_activity_name;
                input_struct_val.current_subject_idx_val = current_subject_idx;
                input_struct_val.current_segment_val = current_segment;
                input_struct_val.current_side_val = current_side;

                [current_technical_markers_val, current_markers_val] = generate_technical_markers_rigidBody(input_struct_val, activity_reprocessed);
                technical_frame_markers_rigidBody_info.(current_activity_name)(current_subject_idx).(current_segment).(current_side) = strings(3,1);
                all_markers_info.(current_activity_name)(current_subject_idx).(current_segment).(current_side) = strings(1,1);
                if ~isempty(current_technical_markers_val)
                    technical_frame_markers_rigidBody_info.(current_activity_name)(current_subject_idx).(current_segment).(current_side) = current_technical_markers_val;
                    all_markers_info.(current_activity_name)(current_subject_idx).(current_segment).(current_side) = current_markers_val;
                end
                
                D = [s, current_segment, current_side, "Finished Processing"];
                disp(D);
            end
        end
    end
end
save('technical_frame_markers_rigidBody_info_data', 'technical_frame_markers_rigidBody_info', '-v7.3');
save('all_markers_info_data', 'all_markers_info', '-v7.3');

function [current_technical_markers, marker_names_rb_list] = generate_technical_markers_rigidBody(input_struct, activity)
    current_activity_idx = input_struct.current_activity_idx_val;
    current_activity_name = input_struct.current_activity_name_val;
    current_subject_idx = input_struct.current_subject_idx_val;
    current_segment = input_struct.current_segment_val;
    current_side = input_struct.current_side_val;
    marker_names_all_list = input_struct.marker_names_all_list_val;

    default_trial_idx = 1;
    static_cal_file_idx = 22;
    
    sideBound = 160;

    marker_names_all = struct;
    marker_names_all.femoral.left = marker_names_all_list([10:15, 55:56]);
    marker_names_all.femoral.right = marker_names_all_list([31:36, 63:64]);
    marker_names_all.tibial.left = marker_names_all_list([19:22, 57:59]);
    marker_names_all.tibial.right = marker_names_all_list([40:44, 65:67]);

    rigidBody_markers = struct;
    activity_idx_list = [static_cal_file_idx, current_activity_idx];
    marker_names_rb_list = marker_names_all.(current_segment).(current_side);
    for j=2:length(activity_idx_list)
        current_iterate_activity_idx = activity_idx_list(j);
        if ~(all(activity(current_iterate_activity_idx).VICON(current_subject_idx, default_trial_idx).markers.present == 0))
            for k=1:length(marker_names_rb_list)
                current_marker_name = marker_names_rb_list(k);
                current_marker_idx = marker_names_all_list == current_marker_name;
                if ~activity(current_iterate_activity_idx).VICON(current_subject_idx, default_trial_idx).markers.present(current_marker_idx)
                    marker_names_rb_list(k)="str_rem";
                end
            end
        end
    end
    marker_names_rb_list = marker_names_rb_list(marker_names_rb_list~="str_rem");
    rigidBody_markers.(current_activity_name)(current_subject_idx).(current_segment).(current_side) = marker_names_rb_list;

    time_vector = activity(current_activity_idx).VICON(current_subject_idx, default_trial_idx).movement_data.time;
    num_measurements = length(time_vector);

    all_marker_info = struct;

    current_rigidBody_markers = rigidBody_markers.(current_activity_name)(current_subject_idx).(current_segment).(current_side);
    combinations_idx = nchoosek(1:length(current_rigidBody_markers), 3);

    num_combinations = length(combinations_idx);
    technical_marker_array = strings(num_combinations, 3);
    sdSum_list = nan(num_combinations, 1);
    maxSide_list = nan(num_combinations, 1);
    for j=1:num_combinations

        current_triangle_idx = combinations_idx(j,:);
        current_triangle_markers = current_rigidBody_markers(current_triangle_idx);
        current_triangle_marker_info = generate_triangle_marker_info(activity, current_triangle_markers, current_activity_idx, current_subject_idx);
        current_triangle_info = generate_triangle_info(current_triangle_marker_info);
        current_triangle_info_validated = triangle_sanity_check(current_triangle_info, current_triangle_marker_info, num_measurements);

        current_sd_sum = 0;
        current_mean_side_list = nan(3,1);
        for c=1:3
            current_mean_side = mean(current_triangle_info_validated(c).side, 'omitnan');
            current_angle_vector = current_triangle_info_validated(c).angleDeg;
            current_sd = std(current_angle_vector, 'omitnan');
            current_mean_side_list(c) = current_mean_side;
            current_sd_sum = current_sd_sum + current_sd;
        end
        current_max_side = max(current_mean_side_list);
        D = [j, current_sd_sum, current_max_side];
        
        if current_max_side < 150
            disp("triangle sides within upper bound")
            disp(D);
        end

        all_marker_info(j).names = current_triangle_markers;
        all_marker_info(j).sdSum = current_sd_sum;

        technical_marker_array(j,:) = current_triangle_markers;
        sdSum_list(j) = current_sd_sum;
        maxSide_list(j) = current_max_side;

    end
    
    sideBound_idx = maxSide_list < sideBound;
    sdSum_bound_list = sdSum_list(sideBound_idx);
    technical_marker_bound_array = technical_marker_array(sideBound_idx, :);
    
    sdMin_idx = sdSum_bound_list == min(sdSum_bound_list);
    current_technical_markers = transpose(technical_marker_bound_array(sdMin_idx, :));
end

function triangle_info_validated = triangle_sanity_check(triangle_info, triangle_marker_info, num_measurements)
    for n=1:num_measurements
        r1_now = [triangle_marker_info(1).x; triangle_marker_info(1).y; triangle_marker_info(1).z];
        r2_now = [triangle_marker_info(2).x; triangle_marker_info(2).y; triangle_marker_info(2).z];
        r3_now = [triangle_marker_info(3).x; triangle_marker_info(3).y; triangle_marker_info(3).z];
        ptCloud = pointCloud(horzcat(r1_now,r2_now,r3_now));
        pts_collinear = true;
        if ~any(any(isnan(ptCloud.Location)))
            ptCloud_tp = pointCloud(transpose(ptCloud.Location));
            [~,plane_inlierIndices,~] = pcfitplane(ptCloud_tp,1);
            triangle_inequality_satisfied_bool = apply_triangle_inequality(ptCloud);
            if or(~isempty(plane_inlierIndices), triangle_inequality_satisfied_bool)
                pts_collinear = false;
            end
        end

        if ~pts_collinear
            triangle_sum = triangle_info(1).angleDeg(n) + triangle_info(2).angleDeg(n) + triangle_info(3).angleDeg(n);
            if (abs(triangle_sum - 180) > tol)
                for c=1:3
                    triangle_info(c).angleDeg(n) = nan;
                    triangle_info(c).angleRad(n) = nan;
                    triangle_info(c).side(n) = nan;
                    disp("Sum of Angles is not 180 degrees")
                end
            end
        end
        
    end
    triangle_info_validated = triangle_info;
    
end

function triangle_info = generate_triangle_info(marker_info)
    triangle_idx = [1, 2, 3];
    triangle_info = struct;
    triangle_adjacent_idx = circshift(triangle_idx, 2);
    triangle_opposite_idx = circshift(triangle_idx, 4);
    for k=1:3
        current_idx = triangle_idx(k);
        adjacent_idx = triangle_adjacent_idx(k);
        opposite_idx = triangle_opposite_idx(k);
        current_side_x_val = marker_info(current_idx).x - marker_info(adjacent_idx).x;
        current_side_y_val = marker_info(current_idx).y - marker_info(adjacent_idx).y;
        current_side_z_val = marker_info(current_idx).z - marker_info(adjacent_idx).z;
        current_side_val = sqrt(current_side_x_val.^2 + current_side_y_val.^2 + current_side_z_val.^2);
        current_side_vector = horzcat(current_side_x_val, current_side_y_val, current_side_z_val);
        adjacent_side_x_val = marker_info(current_idx).x - marker_info(opposite_idx).x;
        adjacent_side_y_val = marker_info(current_idx).y - marker_info(opposite_idx).y;
        adjacent_side_z_val = marker_info(current_idx).z - marker_info(opposite_idx).z;
        adjacent_side_vector = horzcat(adjacent_side_x_val, adjacent_side_y_val, adjacent_side_z_val);
        current_angle_rad = acos(dot(current_side_vector, adjacent_side_vector, 2)./(vecnorm(current_side_vector, 2, 2).*vecnorm(adjacent_side_vector, 2, 2)));
        current_angle_deg = rad2deg(current_angle_rad);
        triangle_info(k).side = current_side_val;
        triangle_info(k).angleRad = current_angle_rad;
        triangle_info(k).angleDeg = current_angle_deg;
    end
end

function marker_info = generate_triangle_marker_info(activity, triangle_markers, activity_idx, subject_idx)
    marker_info = struct;
    for j=1:3
        marker_x = triangle_markers(j) + "_x";
        marker_y = triangle_markers(j) + "_y";
        marker_z = triangle_markers(j) + "_z";
        if isfield(activity(activity_idx).VICON(subject_idx, 1).movement_data, marker_x)
            marker_info(j).x = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_x);
            marker_info(j).y = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_y);
            marker_info(j).z = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_z);
        else
            marker_info(j).x = nan;
            marker_info(j).y = nan;
            marker_info(j).z = nan;
        end
    end
end

function triangle_inequality_satisfied = apply_triangle_inequality(ptCloud)
    triangle_inequality_satisfied = false;
    tol = 1.0e-4;
    if ~any(any(isnan(ptCloud.Location)))
        r1 = ptCloud.Location(:,1);
        r2 = ptCloud.Location(:,2);
        r3 = ptCloud.Location(:,3);
        r_12 = norm(r2-r1);
        r_13 = norm(r3-r1);
        r_23 = norm(r3-r2);
        sum_of_sides = r_12 + r_23;
        third_side = r_13;
        diff_side = abs(sum_of_sides - third_side);
        if diff_side > tol
            triangle_inequality_satisfied = true;
        end
    end
end