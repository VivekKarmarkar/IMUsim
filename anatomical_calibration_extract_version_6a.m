%% LOAD DATA

load_uvm_data_processed_bool = true;
if load_uvm_data_processed_bool
    load('StaticCal_RawData_reprocessed.mat');
end
load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');
load('sorted_vertex_walking_rigidBody_info_data');

%% SWAP DATA

activity = activity_reprocessed;
clearvars activity_reprocessed;

%% SET PREFERENCES

activity_idx = 5;

%% PROCESS DATA

activity_name = strrep(activity(activity_idx).name, " ", "_");

default_trial_idx = 1;
static_cal_idx = 22;

num_subjects = 25;

side_list = ["right", "left"];
num_side = length(side_list);

num_body_segment = length(body_segment_input);
segment_list = strings(1,2);
for k=1:num_body_segment
    segment_list(1,k) = body_segment_input(k).name;
end

input_data_val = struct;
input_data_val.side_list_val = side_list;
input_data_val.segment_list_val = segment_list;
input_data_val.activity_idx_val = static_cal_idx;
input_data_val.technical_frame_markers_info_val = technical_frame_markers_info;
input_data_val.sorted_vertex_old_info_val = sorted_vertex_walking_rigidBody_info.(activity_name);

anatomical_calibration_parameters_processed_rigid_body_info = struct;

for n=1:num_subjects
    disp(n)
    current_subject_idx = n;
    static_cal_empty = all(~activity(static_cal_idx).VICON(n, default_trial_idx).markers.present);
    input_data_val.subject_idx_val = current_subject_idx;
    [local_vectors_output, marker_output] = generate_solidification_data(input_data_val, activity);
    for j=1:num_body_segment
        current_segment = body_segment_input(j).name;
        for k=1:num_side
            current_side = side_list(k);
            current_anatomical_marker_info = struct;
            anatomical_marker_name_list = anatomical_markers_info(current_subject_idx).(current_segment).(current_side).names;
            num_marker = length(anatomical_marker_name_list);
            for m=1:num_marker
                current_anatomical_marker_info(m).name = strings(1,1);
                current_anatomical_marker_info(m).x = nan;
                current_anatomical_marker_info(m).y = nan;
                current_anatomical_marker_info(m).z = nan;
                if ~static_cal_empty
                    current_technical_markers = technical_frame_markers_info(current_subject_idx).(current_segment).(current_side);
                    current_markers = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).markers;
                    current_technical_markers_idx = find(ismember(current_markers.names, current_technical_markers));
                    current_technical_markers_present = all(current_markers.present(current_technical_markers_idx));
                    current_anatomical_markers_idx = find(ismember(current_markers.names, anatomical_marker_name_list));
                    current_anatomical_markers_present = all(current_markers.present(current_anatomical_markers_idx));
                    if current_technical_markers_present && current_anatomical_markers_present
                        current_anatomical_marker_name_val = anatomical_marker_name_list(m);
                        current_anatomical_marker_info(m).name = current_anatomical_marker_name_val;
                        t_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.time;
                        x1_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(1)+'_x');
                        y1_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(1)+'_y');
                        z1_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(1)+'_z');
                        r1_list = horzcat(x1_list,y1_list,z1_list);
                        x2_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(2)+'_x');
                        y2_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(2)+'_y');
                        z2_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(2)+'_z');
                        r2_list = horzcat(x2_list,y2_list,z2_list);
                        x3_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(3)+'_x');
                        y3_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(3)+'_y');
                        z3_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_technical_markers(3)+'_z');
                        r3_list = horzcat(x3_list,y3_list,z3_list);
                        xa_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name_val+'_x');
                        ya_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name_val+'_y');
                        za_list = activity(static_cal_idx).VICON(current_subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name_val+'_z');
                        ra_list = horzcat(xa_list,ya_list,za_list);
                        num_measurements = length(t_list);
                        position_local_x_list = nan(num_measurements,1);
                        position_local_y_list = nan(num_measurements,1);
                        position_local_z_list = nan(num_measurements,1);
                        for t=1:num_measurements
                            t_now = t_list(t);
                            ra = ra_list(t,:)';
                            
                            marker_input_struct = struct;
                            marker_input_struct.technical_markers_list_val = marker_output.(current_segment).(current_side).technical_markers_list_val;
                            marker_input_struct.marker_info_val = marker_output.(current_segment).(current_side).marker_info_val;
                            marker_input_struct.sorted_vertex_val = marker_output.(current_segment).(current_side).sorted_vertex_val;
                            marker_input_struct.success_bool_val = marker_output.(current_segment).(current_side).success_bool_val;
                            marker_input_struct.sorted_vertex_val

                            local_vectors_input_struct = struct;
                            local_vectors_input_struct.local_vectors_val = local_vectors_output.(current_segment).(current_side).local_vectors_val;
                            local_vectors_input_struct.local_vectors_mean_val = local_vectors_output.(current_segment).(current_side).local_vectors_mean_val;
                            local_vectors_input_struct.q_local_val = local_vectors_output.(current_segment).(current_side).q_local_val;
                            local_vectors_input_struct.solid_triangle_val = local_vectors_output.(current_segment).(current_side).solid_triangle_val;

                            [T, R, unit_vector_hat, success_bool, reconstructed_vectors] = pose_estimate_technical_frame_rigid_body(t, marker_input_struct, local_vectors_input_struct);
                            
                            position_local = global2local(T, R, ra);
                            position_local_x_list(t) = position_local(1);
                            position_local_y_list(t) = position_local(2);
                            position_local_z_list(t) = position_local(3);
                        end
                        current_anatomical_marker_info(m).x = mean(position_local_x_list, 'omitnan');
                        current_anatomical_marker_info(m).y = mean(position_local_y_list, 'omitnan');
                        current_anatomical_marker_info(m).z = mean(position_local_z_list, 'omitnan');
                    end
                end
            end
            anatomical_calibration_parameters_processed_rigid_body_info(current_subject_idx).(current_segment).(current_side) = current_anatomical_marker_info;
            disp('----finished processing this round----')
        end
    end
end
save('anatomical_calibration_parameters_processed_rigid_body_info_data', 'anatomical_calibration_parameters_processed_rigid_body_info', '-v7.3');

function [local_vectors_output, marker_output] = generate_solidification_data(input_data, activity)
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val;
    segment_list = input_data.segment_list_val;
    side_list = input_data.side_list_val;
    sorted_vertex_old_info = input_data.sorted_vertex_old_info_val;
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
            segment_meta_data_struct.sorted_vertex_old_info_val = sorted_vertex_old_info;
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
    sorted_vertex_old_info = segment_meta_data.sorted_vertex_old_info_val;
    sorted_vertex_old = sorted_vertex_old_info(subject_idx).(segment).(side);

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
            marker_info(j).x = nan;
            marker_info(j).y = nan;
            marker_info(j).z = nan;
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
    
    if all(sorted_vertex_old == circshift(sorted_vertex, 1))
        
        local_vectors(1).name = technical_markers_list(sorted_vertex(3));
        local_vectors(1).x = solid_triangle(3).side;
        local_vectors(1).y = 0;
        local_vectors(1).z = 0;
        
        local_vectors(2).name = technical_markers_list(sorted_vertex(1));
        local_vectors(2).xTrue = 0;
        local_vectors(2).x = solid_triangle(3).side - 1;
        local_vectors(2).y = 0;
        local_vectors(2).z = 0;
        
        local_vectors(3).name = technical_markers_list(sorted_vertex(2));
        local_vectors(3).xTrue = solid_triangle(1).side * cos(solid_triangle(1).angleRad);
        local_vectors(3).yTrue = solid_triangle(1).side * sin(solid_triangle(1).angleRad);
        local_vectors(3).x = solid_triangle(3).side - cos(solid_triangle(3).angleRad);
        local_vectors(3).y = sin(solid_triangle(3).angleRad);
        local_vectors(3).z = 0;
        
    elseif all(sorted_vertex_old == circshift(sorted_vertex, 2))
        
        local_vectors(1).name = technical_markers_list(sorted_vertex(2));
        local_vectors(1).x = solid_triangle(2).side;
        local_vectors(1).y = 0;
        local_vectors(1).z = 0;
        
        local_vectors(2).name = technical_markers_list(sorted_vertex(3));
        local_vectors(2).xTrue = 0;
        local_vectors(2).x = solid_triangle(2).side - 1;
        local_vectors(2).y = 0;
        local_vectors(2).z = 0;
        
        local_vectors(3).name = technical_markers_list(sorted_vertex(1));
        local_vectors(3).xTrue = solid_triangle(3).side * cos(solid_triangle(3).angleRad);
        local_vectors(3).yTrue = solid_triangle(3).side * sin(solid_triangle(3).angleRad);
        local_vectors(3).x = solid_triangle(2).side - cos(solid_triangle(2).angleRad);
        local_vectors(3).y = sin(solid_triangle(2).angleRad);
        local_vectors(3).z = 0;
        
    elseif all(sorted_vertex_old == circshift(sorted_vertex, 3))
        
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
        
    elseif sorted_vertex_old(1) == sorted_vertex(1)
        
        local_vectors(1).name = technical_markers_list(sorted_vertex(1));
        local_vectors(1).x = solid_triangle(3).side;
        local_vectors(1).y = 0;
        local_vectors(1).z = 0;
        
        local_vectors(2).name = technical_markers_list(sorted_vertex(3));
        local_vectors(2).xTrue = 0;
        local_vectors(2).x = solid_triangle(3).side - 1;
        local_vectors(2).y = 0;
        local_vectors(2).z = 0;
        
        local_vectors(3).name = technical_markers_list(sorted_vertex(2));
        local_vectors(3).xTrue = solid_triangle(2).side * cos(solid_triangle(3).angleRad);
        local_vectors(3).yTrue = solid_triangle(2).side * sin(solid_triangle(3).angleRad);
        local_vectors(3).x = solid_triangle(3).side - cos(solid_triangle(1).angleRad);
        local_vectors(3).y = sin(solid_triangle(1).angleRad);
        local_vectors(3).z = 0;
        
    elseif sorted_vertex_old(end) == sorted_vertex(end)
        
        local_vectors(1).name = technical_markers_list(sorted_vertex(2));
        local_vectors(1).x = solid_triangle(1).side;
        local_vectors(1).y = 0;
        local_vectors(1).z = 0;
        
        local_vectors(2).name = technical_markers_list(sorted_vertex(1));
        local_vectors(2).xTrue = 0;
        local_vectors(2).x = solid_triangle(1).side - 1;
        local_vectors(2).y = 0;
        local_vectors(2).z = 0;
        
        local_vectors(3).name = technical_markers_list(sorted_vertex(3));
        local_vectors(3).xTrue = solid_triangle(3).side * cos(solid_triangle(1).angleRad);
        local_vectors(3).yTrue = solid_triangle(3).side * sin(solid_triangle(1).angleRad);
        local_vectors(3).x = solid_triangle(1).side - cos(solid_triangle(2).angleRad);
        local_vectors(3).y = sin(solid_triangle(2).angleRad);
        local_vectors(3).z = 0;
        
    else
        
        local_vectors(1).name = technical_markers_list(sorted_vertex(3));
        local_vectors(1).x = solid_triangle(2).side;
        local_vectors(1).y = 0;
        local_vectors(1).z = 0;
        
        local_vectors(2).name = technical_markers_list(sorted_vertex(2));
        local_vectors(2).xTrue = 0;
        local_vectors(2).x = solid_triangle(2).side - 1;
        local_vectors(2).y = 0;
        local_vectors(2).z = 0;
        
        local_vectors(3).name = technical_markers_list(sorted_vertex(1));
        local_vectors(3).xTrue = solid_triangle(1).side * cos(solid_triangle(2).angleRad);
        local_vectors(3).yTrue = solid_triangle(1).side * sin(solid_triangle(2).angleRad);
        local_vectors(3).x = solid_triangle(2).side - cos(solid_triangle(3).angleRad);
        local_vectors(3).y = sin(solid_triangle(3).angleRad);
        local_vectors(3).z = 0;
    
    end

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
        
        global_vectors(c).name = technical_markers_list(sorted_vertex_old(c));
        global_vectors(c).x = marker_info(sorted_vertex_old(c)).x;
        global_vectors(c).y = marker_info(sorted_vertex_old(c)).y;
        global_vectors(c).z = marker_info(sorted_vertex_old(c)).z;

        q_local(c).x = local_vectors(c).x - local_vectors_mean.x;
        q_local(c).y = local_vectors(c).y - local_vectors_mean.y;
        q_local(c).z = local_vectors(c).z - local_vectors_mean.z;
        q_local(c).vec = [q_local(c).x; q_local(c).y; q_local(c).z];
    end
    
    Local_sideOne = sqrt((local_vectors(1).x - local_vectors(2).x)^2 + (local_vectors(1).y - local_vectors(2).y)^2 + (local_vectors(1).z - local_vectors(2).z)^2);
    Local_sideTwo = sqrt((local_vectors(2).x - local_vectors(3).x)^2 + (local_vectors(2).y - local_vectors(3).y)^2 + (local_vectors(2).z - local_vectors(3).z)^2);
    Local_sideThree = sqrt((local_vectors(1).x - local_vectors(3).x)^2 + (local_vectors(1).y - local_vectors(3).y)^2 + (local_vectors(1).z - local_vectors(3).z)^2);
    disp([Local_sideOne, Local_sideTwo, Local_sideThree]);

    local_vectors_output = struct;
    local_vectors_output.local_vectors_val = local_vectors;
    local_vectors_output.local_vectors_mean_val = local_vectors_mean;
    local_vectors_output.q_local_val = q_local;
    local_vectors_output.solid_triangle_val = solid_triangle;

    marker_output = struct;
    marker_output.technical_markers_list_val = technical_markers_list;
    marker_output.marker_info_val = marker_info;
    marker_output.sorted_vertex_val = sorted_vertex_old;
    marker_output.success_bool_val = success_bool_val;

end

function [T_val, R_val, unit_vectors, success_bool, reconstructed_vectors] = pose_estimate_technical_frame_rigid_body(current_time_idx, marker_input, local_vectors_input)
    
    technical_markers_list = marker_input.technical_markers_list_val;
    marker_info = marker_input.marker_info_val;
    sorted_vertex = marker_input.sorted_vertex_val;
    success_bool = marker_input.success_bool_val;

    local_vectors = local_vectors_input.local_vectors_val;
    local_vectors_mean = local_vectors_input.local_vectors_mean_val;
    q_local = local_vectors_input.q_local_val;
    solid_triangle = local_vectors_input.solid_triangle_val;
    
    current_global_vectors_unit = struct;
    current_global_vectors_true = struct;
    for c=1:3
        current_global_vectors_unit(c).name = technical_markers_list(sorted_vertex(c));
        current_global_vectors_true(c).name = technical_markers_list(sorted_vertex(c));
        if length(marker_info(sorted_vertex(c)).x) > 1
            if c == 2
                s = solid_triangle(1).side;
            elseif c == 3
                s = solid_triangle(3).side;
            else
                s = 1;
            end
            if current_time_idx <= length(marker_info(1).x)
                Lx = abs(marker_info(sorted_vertex(c)).x(current_time_idx) - marker_info(sorted_vertex(1)).x(current_time_idx));
                Ly = abs(marker_info(sorted_vertex(c)).y(current_time_idx) - marker_info(sorted_vertex(1)).y(current_time_idx));
                Lz = abs(marker_info(sorted_vertex(c)).z(current_time_idx) - marker_info(sorted_vertex(1)).z(current_time_idx));
                L = sqrt(Lx^2 + Ly^2 + Lz^2);
                if c==1
                    L = 1;
                end
                NewX_val = (marker_info(sorted_vertex(c)).x(current_time_idx) + (L-1)*marker_info(sorted_vertex(1)).x(current_time_idx))/L;
                NewY_val = (marker_info(sorted_vertex(c)).y(current_time_idx) + (L-1)*marker_info(sorted_vertex(1)).y(current_time_idx))/L;
                NewZ_val = (marker_info(sorted_vertex(c)).z(current_time_idx) + (L-1)*marker_info(sorted_vertex(1)).z(current_time_idx))/L;
                current_global_vectors_unit(c).x = NewX_val;
                current_global_vectors_unit(c).y = NewY_val;
                current_global_vectors_unit(c).z = NewZ_val;
                current_global_vectors_true(c).x = marker_info(sorted_vertex(c)).x(current_time_idx);
                current_global_vectors_true(c).y = marker_info(sorted_vertex(c)).y(current_time_idx);
                current_global_vectors_true(c).z = marker_info(sorted_vertex(c)).z(current_time_idx);
            end
        else
            current_global_vectors_unit(c).x = nan;
            current_global_vectors_unit(c).y = nan;
            current_global_vectors_unit(c).z = nan;
            current_global_vectors_true(c).x = nan;
            current_global_vectors_true(c).y = nan;
            current_global_vectors_true(c).z = nan;
        end
        current_global_vectors_unit(c).vec = [current_global_vectors_unit(c).x; current_global_vectors_unit(c).y; current_global_vectors_unit(c).z];
        current_global_vectors_true(c).vec = [current_global_vectors_true(c).x; current_global_vectors_true(c).y; current_global_vectors_true(c).z];
    end

    current_global_vectors_unit_mean = struct;
    current_global_vectors_unit_mean.x = mean([current_global_vectors_unit(1).x, current_global_vectors_unit(2).x, current_global_vectors_unit(3).x], 'omitnan');
    current_global_vectors_unit_mean.y = mean([current_global_vectors_unit(1).y, current_global_vectors_unit(2).y, current_global_vectors_unit(3).y], 'omitnan');
    current_global_vectors_unit_mean.z = mean([current_global_vectors_unit(1).z, current_global_vectors_unit(2).z, current_global_vectors_unit(3).z], 'omitnan');
    current_global_vectors_unit_mean.vec = [current_global_vectors_unit_mean.x; current_global_vectors_unit_mean.y; current_global_vectors_unit_mean.z];

    q_current_global = struct;
    h_matrix = zeros(3,3);
    for c=1:3
        q_current_global(c).x = current_global_vectors_unit(c).x - current_global_vectors_unit_mean.x;
        q_current_global(c).y = current_global_vectors_unit(c).y - current_global_vectors_unit_mean.y;
        q_current_global(c).z = current_global_vectors_unit(c).z - current_global_vectors_unit_mean.z;
        q_current_global(c).vec = [q_current_global(c).x; q_current_global(c).y; q_current_global(c).z];

        h_matrix = h_matrix + (q_local(c).vec * transpose(q_current_global(c).vec));
    end

    sanityCheck_tol = 10^-8;
    X = nan(3,3);
    if ~any(any(isnan(h_matrix)))
        [U, ~, V] = svd(h_matrix);
        X = V*U';
    end

    if ~all(all(isnan(X)))
        if abs(det(X) - 1) > sanityCheck_tol
            if abs(det(X) + 1) < sanityCheck_tol
                V_new = horzcat(V(:,1:2), -V(:,3));
                X = V_new*U';
            else
                X = nan(3,3);
            end
        end
    end

    R_val = X;
    T_val = current_global_vectors_unit_mean.vec - R_val*local_vectors_mean.vec;
    
    if or(all(all(isnan(R_val))), all(isnan(T_val)))
        success_bool = false;
    end

    reconstructed_vectors = struct;
    for c=1:3
        reconstructed_vectors(c).name = technical_markers_list(sorted_vertex(c));
        
        reconstructed_vectors(c).vec = local2global(T_val, R_val, local_vectors(c).vec);
        reconstructed_vectors(c).res = reconstructed_vectors(c).vec - current_global_vectors_unit(c).vec;
        
        reconstructed_vectors(c).vecTrue = local2global(T_val, R_val, local_vectors(c).vecTrue);
        reconstructed_vectors(c).resTrue = reconstructed_vectors(c).vecTrue - current_global_vectors_true(c).vec;
        
        reconstructed_vectors(c).x = reconstructed_vectors(c).vec(1);
        reconstructed_vectors(c).y = reconstructed_vectors(c).vec(2);
        reconstructed_vectors(c).z = reconstructed_vectors(c).vec(3);
        
        reconstructed_vectors(c).xTrue = reconstructed_vectors(c).vecTrue(1);
        reconstructed_vectors(c).yTrue = reconstructed_vectors(c).vecTrue(2);
        reconstructed_vectors(c).zTrue = reconstructed_vectors(c).vecTrue(3);
        
        reconstructed_vectors(c).Xres = reconstructed_vectors(c).res(1);
        reconstructed_vectors(c).Yres = reconstructed_vectors(c).res(2);
        reconstructed_vectors(c).Zres = reconstructed_vectors(c).res(3);
        
        reconstructed_vectors(c).XresTrue = reconstructed_vectors(c).resTrue(1);
        reconstructed_vectors(c).YresTrue = reconstructed_vectors(c).resTrue(2);
        reconstructed_vectors(c).ZresTrue = reconstructed_vectors(c).resTrue(3);
    end
    
    i_hat_local = [1;0;0];
    j_hat_local = [0;1;0];
    k_hat_local = [0;0;1];
    
    i_hat = R_val * i_hat_local;
    j_hat = R_val * j_hat_local;
    k_hat = R_val * k_hat_local;
    
    if abs(norm(i_hat) - 1) > sanityCheck_tol
        i_hat = nan(3,1);
    end
    
    if abs(norm(j_hat) - 1) > sanityCheck_tol
        j_hat = nan(3,1);
    end
    
    if abs(norm(k_hat) - 1) > sanityCheck_tol
        k_hat = nan(3,1);
    end
    
    if abs(dot(i_hat, j_hat)) > sanityCheck_tol
        i_hat = nan(3,1);
        j_hat = nan(3,1);
    end
    
    if abs(dot(i_hat, k_hat)) > sanityCheck_tol
        i_hat = nan(3,1);
        k_hat = nan(3,1);
    end
    
    if abs(dot(j_hat, k_hat)) > sanityCheck_tol
        j_hat = nan(3,1);
        k_hat = nan(3,1);
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end

function position_local = global2local(T, R, position_global)
    position_local = nan(3,1);
    tol = 1.0e-10;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_local = R'*(position_global - T);
        end
    end
end

function position_global = local2global(T, R, position_local)
    position_global = nan(3,1);
    tol = 1.0e-10;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_global = T + R*position_local;
        end
    end
end

function [success_bool_val, solid_triangle] = triangle_sanity_check(solid_triangle)
    success_bool_val = true;
    tol = 10^-4;
    angleSum = solid_triangle(1).angleDeg + solid_triangle(2).angleDeg + solid_triangle(3).angleDeg;
    angleSum_bool = abs(angleSum - 180) > tol;
    sideSum_bool = solid_triangle(1).side + solid_triangle(2).side < solid_triangle(3).side;
    sineLaw_bool =  abs((solid_triangle(1).side/sin(solid_triangle(3).angleRad)) - (solid_triangle(2).side/sin(solid_triangle(1).angleRad))) > tol;
    if or(or(angleSum_bool , sideSum_bool) , sineLaw_bool)
        for c=1:3
            solid_triangle(c).angleDeg = nan;
            solid_triangle(c).angleRad = nan;
            solid_triangle(c).side = nan;
            success_bool_val = false;
        end
    end
end