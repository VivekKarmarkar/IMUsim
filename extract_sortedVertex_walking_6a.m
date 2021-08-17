%% LOAD HEAVY DATA

clc
clear
close

load('WalkingFP_RawData_reprocessed.mat');

%% LOAD NON HEAVY DATA

load('body_segment_input_data.mat');
load('technical_frame_markers_info_data.mat');

%% GENERATE AND POPULATE DATA CONTAINERS

side_list = ["right", "left"];
segment_list = strings(1,2);
for k=1:length(body_segment_input)
    segment_list(1,k) = body_segment_input(k).name;
end
num_subjects = 25;
activity_start_idx = 5;
activity_end_idx = 14;

input_data_val = struct;
input_data_val.side_list_val = side_list;
input_data_val.segment_list_val = segment_list;
input_data_val.technical_frame_markers_info_val = technical_frame_markers_info;

sorted_vertex_walking_rigidBody_info = struct;

%% PROCESSING AND SAVING DATA

for j=activity_start_idx:activity_end_idx
    current_activity_idx = j;
    input_data_val.activity_idx_val = current_activity_idx;
    current_activity_name = strrep(activity_reprocessed(current_activity_idx).name, " ", "_");
    disp(current_activity_name);
    
    for s=1:num_subjects
        disp(s);
        input_data_val.subject_idx_val = s;
        [local_vectors_output, marker_output] = generate_solidification_data(input_data_val, activity_reprocessed);
        for segment_idx = 1:length(segment_list)
            current_segment = segment_list(segment_idx);
            for side_idx = 1:length(side_list)
                current_side = side_list(side_idx);
                sorted_vertex_walking_rigidBody_info.(current_activity_name)(s).(current_segment).(current_side) = marker_output.(current_segment).(current_side).sorted_vertex_val;
            end
        end
    end
    
end
save('sorted_vertex_walking_rigidBody_info_data', 'sorted_vertex_walking_rigidBody_info', '-v7.3');

function [local_vectors_output, marker_output] = generate_solidification_data(input_data, activity)
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val;
    segment_list = input_data.segment_list_val;
    side_list = input_data.side_list_val;
    for segment_idx=1:length(segment_list)
        current_segment = segment_list(segment_idx);
        for side_idx=1:length(side_list)
            current_side = side_list(side_idx);
            segment_meta_data_struct = struct;
            segment_meta_data_struct.activity_idx_val = activity_idx;
            segment_meta_data_struct.subject_idx_val = subject_idx;
            segment_meta_data_struct.segment_val = current_segment;
            segment_meta_data_struct.side_val = current_side;
            segment_meta_data_struct.technical_frame_markers_info_val = technical_frame_markers_info;
            [current_local_vectors_output, current_marker_output] = triangle_solidification(segment_meta_data_struct, activity);
            local_vectors_output.(current_segment).(current_side) = current_local_vectors_output;
            marker_output.(current_segment).(current_side) = current_marker_output;
        end
    end
end

function [local_vectors_output, marker_output] = triangle_solidification(segment_meta_data, activity)

    activity_idx = segment_meta_data.activity_idx_val;
    subject_idx = segment_meta_data.subject_idx_val;
    segment = segment_meta_data.segment_val;
    side = segment_meta_data.side_val;
    technical_frame_markers_info = segment_meta_data.technical_frame_markers_info_val;
    technical_markers_list = technical_frame_markers_info(subject_idx).(segment).(side);

    num_measurements = length(activity(activity_idx).VICON(subject_idx, 1).movement_data.time);

    triangle_info = struct;
    marker_info = struct;
    for j=1:3
        marker_x = technical_markers_list(j) + "_x";
        marker_y = technical_markers_list(j) + "_y";
        marker_z = technical_markers_list(j) + "_z";
        if isfield(activity(activity_idx).VICON(subject_idx, 1).movement_data, marker_x)
            marker_info(j).x = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_x);
            marker_info(j).y = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_y);
            marker_info(j).z = activity(activity_idx).VICON(subject_idx, 1).movement_data.(marker_z);
        else
            if num_measurements > 1
                marker_info(j).x = nan(num_measurements, 1);
                marker_info(j).y = nan(num_measurements, 1);
                marker_info(j).z = nan(num_measurements, 1);
            else
                marker_info(j).x = nan;
                marker_info(j).y = nan;
                marker_info(j).z = nan;
            end
        end
    end

    triangle_idx = [1, 2, 3];
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

    all_idx = 1:num_measurements;
    nancount_cutoff = 0.3 * length(all_idx);
    tol = 10^-4;
    for n=1:num_measurements
        triangle_sum = triangle_info(1).angleDeg(n) + triangle_info(2).angleDeg(n) + triangle_info(3).angleDeg(n);
        if (abs(triangle_sum - 180) > tol)
            for c=1:3
                triangle_info(c).angleDeg(n) = nan;
                triangle_info(c).angleRad(n) = nan;
                triangle_info(c).side(n) = nan;
                disp("Sum of Angles is not 180 degrees")
            end
        end
        sum_angleDeg = nan(num_measurements,1);
        for m=1:num_measurements
            current_sum = 0;
            for c=1:3
                angleDeg_listVal = triangle_info(c).angleDeg;
                current_angleDeg_listVal = vertcat(angleDeg_listVal(1:n-1), nan, angleDeg_listVal(n+1:end));
                current_angleDeg_meanVal = mean(current_angleDeg_listVal, 'omitnan');
                current_term = (current_angleDeg_listVal(m) - current_angleDeg_meanVal)^2;
                current_sum = current_sum + current_term;
            end
            sum_angleDeg(m) = current_sum;
        end
        max_idx = find(sum_angleDeg == max(sum_angleDeg));
        all_idx(max_idx) = nan;
        for c=1:3
            triangle_info(c).angleDeg(max_idx) = nan;
            triangle_info(c).angleRad(max_idx) = nan;
            triangle_info(c).side(max_idx) = nan;
        end
        current_nancount = length(all_idx(isnan(all_idx)));
        if current_nancount > nancount_cutoff
            break
        end
    end

    angle_deform_list = nan(3,1);
    for c=1:3
        angle_deform_list(c) = std(triangle_info(c).angleDeg, 'omitnan');
    end
    [~, sorted_vertex] = sort(angle_deform_list, 'ascend');

    first_angle_val = mean(triangle_info(sorted_vertex(1)).angleDeg, 'omitnan');
    second_angle_val = mean(triangle_info(sorted_vertex(2)).angleDeg, 'omitnan');
    third_angle_val = 180 - (first_angle_val + second_angle_val);
    first_side_x_val_list = marker_info(sorted_vertex(2)).x - marker_info(sorted_vertex(1)).x;
    first_side_y_val_list = marker_info(sorted_vertex(2)).y - marker_info(sorted_vertex(1)).y;
    first_side_z_val_list = marker_info(sorted_vertex(2)).z - marker_info(sorted_vertex(1)).z;
    first_side_val_list = sqrt(first_side_x_val_list.^2 + first_side_y_val_list.^2 + first_side_z_val_list.^2);
    first_side_val = mean(first_side_val_list, 'omitnan');
    second_side_val = first_side_val * (sin(deg2rad(first_angle_val)) / sin(deg2rad(third_angle_val)));
    third_side_val = first_side_val * (sin(deg2rad(second_angle_val)) / sin(deg2rad(third_angle_val)));

    solid_triangle = struct;

    solid_triangle(1).vertex_idx = sorted_vertex(1);
    solid_triangle(1).angleDeg = first_angle_val;
    solid_triangle(1).angleRad = deg2rad(first_angle_val);
    solid_triangle(1).side = first_side_val;

    solid_triangle(2).vertex_idx = sorted_vertex(2);
    solid_triangle(2).angleDeg = second_angle_val;
    solid_triangle(2).angleRad = deg2rad(second_angle_val);
    solid_triangle(2).side = second_side_val;

    solid_triangle(3).vertex_idx = sorted_vertex(3);
    solid_triangle(3).angleDeg = third_angle_val;
    solid_triangle(3).angleRad = deg2rad(third_angle_val);
    solid_triangle(3).side = third_side_val;

    [success_bool_val, solid_triangle]= triangle_sanity_check(solid_triangle);

    local_vectors = struct;
    
    local_vectors(1).name = technical_markers_list(sorted_vertex(1));
    local_vectors(1).x = solid_triangle(1).side;
    local_vectors(1).y = 0;
    local_vectors(1).z = 0;

    local_vectors(2).name = technical_markers_list(sorted_vertex(2));
    local_vectors(2).xTrue = 0;
    local_vectors(2).x = solid_triangle(1).side - 1;
    local_vectors(2).y = 0;
    local_vectors(2).z = 0;

    local_vectors(3).name = technical_markers_list(sorted_vertex(3));
    local_vectors(3).xTrue = solid_triangle(2).side * cos(solid_triangle(2).angleRad);
    local_vectors(3).yTrue = solid_triangle(2).side * sin(solid_triangle(2).angleRad);
    local_vectors(3).x = solid_triangle(1).side - cos(solid_triangle(1).angleRad);
    local_vectors(3).y = sin(solid_triangle(1).angleRad);
    local_vectors(3).z = 0;
    
    local_vectors(1).vec = [local_vectors(1).x; local_vectors(1).y; local_vectors(1).z];
    local_vectors(2).vec = [local_vectors(2).x; local_vectors(2).y; local_vectors(2).z];
    local_vectors(3).vec = [local_vectors(3).x; local_vectors(3).y; local_vectors(3).z];
    
    local_vectors(1).vecTrue = local_vectors(1).vec;
    local_vectors(2).vecTrue = [local_vectors(2).xTrue; local_vectors(2).y; local_vectors(2).z];
    local_vectors(3).vecTrue = [local_vectors(3).xTrue; local_vectors(3).yTrue; local_vectors(3).z];
    
    local_vectors_mean = struct;
    local_vectors_mean.x = mean([local_vectors(1).x, local_vectors(2).x, local_vectors(3).x], 'omitnan');
    local_vectors_mean.y = mean([local_vectors(1).y, local_vectors(2).y, local_vectors(3).y], 'omitnan');
    local_vectors_mean.z = mean([local_vectors(1).z, local_vectors(2).z, local_vectors(3).z], 'omitnan');
    local_vectors_mean.vec = [local_vectors_mean.x; local_vectors_mean.y; local_vectors_mean.z];

    global_vectors = struct;
    q_local = struct;
    for c=1:3
        global_vectors(c).name = technical_markers_list(sorted_vertex(c));
        global_vectors(c).x = marker_info(sorted_vertex(c)).x;
        global_vectors(c).y = marker_info(sorted_vertex(c)).y;
        global_vectors(c).z = marker_info(sorted_vertex(c)).z;

        q_local(c).x = local_vectors(c).x - local_vectors_mean.x;
        q_local(c).y = local_vectors(c).y - local_vectors_mean.y;
        q_local(c).z = local_vectors(c).z - local_vectors_mean.z;
        q_local(c).vec = [q_local(c).x; q_local(c).y; q_local(c).z];
    end

    local_vectors_output = struct;
    local_vectors_output.local_vectors_val = local_vectors;
    local_vectors_output.local_vectors_mean_val = local_vectors_mean;
    local_vectors_output.q_local_val = q_local;
    local_vectors_output.solid_triangle_val = solid_triangle;

    marker_output = struct;
    marker_output.technical_markers_list_val = technical_markers_list;
    marker_output.marker_info_val = marker_info;
    marker_output.sorted_vertex_val = sorted_vertex;
    marker_output.success_bool_val = success_bool_val;

end

function [success_bool_val, solid_triangle] = triangle_sanity_check(solid_triangle)
    success_bool_val = true;
    tol = 10^-4;
    angleSum = solid_triangle(1).angleDeg + solid_triangle(2).angleDeg + solid_triangle(3).angleDeg;
    angleSum_bool = abs(angleSum - 180) > tol;
    sideSum_bool = solid_triangle(1).side + solid_triangle(2).side < solid_triangle(3).side;
    sineLaw_bool =  abs((solid_triangle(1).side/sin(solid_triangle(3).angleRad)) - (solid_triangle(2).side/sin(solid_triangle(1).angleRad))) > tol;
    maxSide = max([solid_triangle(1).side, solid_triangle(2).side, solid_triangle(3).side]);
    if maxSide >= 150
        disp("Max side greater equals 15 cm");
    end
    if or(or(angleSum_bool , sideSum_bool) , sineLaw_bool)
        for c=1:3
            solid_triangle(c).angleDeg = nan;
            solid_triangle(c).angleRad = nan;
            solid_triangle(c).side = nan;
            success_bool_val = false;
        end
    end
end