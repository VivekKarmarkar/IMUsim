%% LOAD HEAVY DATA

clc
clear
close

load('UVM_data_regrouped.mat');
load('WalkingFP_RawData_reprocessed.mat');

%% SET PREFERENCES

subject_idx_input = 4;
activity_idx = 5;
all_subjects = true;

activityTypeOverGdWalking = true;

plot_flag.rawData = true;
plot_flag.processedData = true;
plot_flag.rigidBody = true;

plot_flag.knee_angle = true;
plot_flag.pose_components = true;

plot_flag.meanValue = false;

animation_flag = false;
animation_movingView = false;
animation_rawData = false;
animation_rigidBody = true;
if animation_rawData
    animation_rigidBody = false;
end

segment_visualize_bool.technical = true;
segment_visualize_bool.reconstructed = true;

ground_contact_visualization_flag = false;

visualization_flag = true;

%% OVERWRITE FLAGS

overwrite_animation_folder = false;
overwrite_visualization_folder = false;
overwrite_GroundCt_visualization_folder = false;
overwrite_plot_folder.knee_angle = false;
overwrite_plot_folder.pose_components = false;

%% LOAD NON HEAVY DATA and ACCOUNT FOR PREFERENCES

if all_subjects
    one_subject = false;
else
    one_subject = true;
end

subject_idx_val = nan;
if one_subject
    subject_idx_val = subject_idx_input;
end

if or(~isfolder("Movement_Videos"), overwrite_animation_folder)
    mkdir Movement_Videos
end

if or(~isfolder("Movement_MovingPtView_Videos"), overwrite_animation_folder)
    mkdir Movement_MovingPtView_Videos
end

if or(~isfolder("Visualization_Videos"), overwrite_visualization_folder)
    mkdir Visualization_Videos
end

if or(~isfolder("Visualization_MovingPtView_Videos"), overwrite_visualization_folder)
    mkdir Visualization_MovingPtView_Videos
end

if or(~isfolder("GroundCt_Visualization_Videos"), overwrite_GroundCt_visualization_folder)
    mkdir GroundCt_Visualization_Videos
end

if or(~isfolder("GroundCt_Visualization_MovingPtView_Videos"), overwrite_GroundCt_visualization_folder)
    mkdir GroundCt_Visualization_MovingPtView_Videos
end

if or(~isfolder("MVT_Knee_Angle_Plots"), overwrite_plot_folder.knee_angle)
    mkdir MVT_Knee_Angle_Plots
end

if or(~isfolder("MVT_Pose_Components_Plots"), overwrite_plot_folder.pose_components)
    mkdir MVT_Pose_Components_Plots
end

load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');
load('technical_frame_markers_rigidBody_info_data.mat');
load('anatomical_calibration_parameters_info_data.mat');
load('anatomical_calibration_parameters_processed_info_data.mat');
load('anatomical_calibration_parameters_processed_rigid_body_info_data.mat');

side_list = ["right", "left"];
static_cal_idx = 22;
default_trial_idx = 1;

activity_name = strrep(activity_reprocessed(activity_idx).name, " ", "_");
now_technical_frame_markers_rb_info = technical_frame_markers_rigidBody_info.(activity_name);
now_anatomical_calibration_parameters_processed_rigid_body_info = anatomical_calibration_parameters_processed_rigid_body_info.(activity_name);

raw_data_bool = true;
rigid_body_bool = false;
ground_contact_rawData = true;

segment_list = strings(1,2);
anatomical_frame_function_idx_list = nan(1,2);
for k=1:length(body_segment_input)
    segment_list(1,k) = body_segment_input(k).name;
    anatomical_frame_function_idx_list(1,k) = body_segment_input(k).anatomical_frame_function_idx;
end

input_data_val = struct;
input_data_val.side_list_val = side_list;
input_data_val.segment_list_val = segment_list;
input_data_val.anatomical_frame_function_idx_list = anatomical_frame_function_idx_list;
input_data_val.activity_idx_val = activity_idx;
input_data_val.anatomical_markers_info_val = anatomical_markers_info;

input_data_val.technical_frame_markers_info_val.raw = technical_frame_markers_info;
input_data_val.technical_frame_markers_info_val.processed = technical_frame_markers_info;
input_data_val.technical_frame_markers_info_val.rigidBody = now_technical_frame_markers_rb_info;

input_data_val.anatomical_calibration_parameters_val.raw = anatomical_calibration_parameters_info;
input_data_val.anatomical_calibration_parameters_val.processed = anatomical_calibration_parameters_processed_info;
input_data_val.anatomical_calibration_parameters_val.rigidBody = now_anatomical_calibration_parameters_processed_rigid_body_info;

meta_data_val = struct;
meta_data_val.side_list_val = side_list;
meta_data_val.segment_list_val = segment_list;
meta_data_val.activity_idx_val = activity_idx;

animation_settings_struct = struct;
animation_settings_struct.rawData_flag_val = animation_rawData;
animation_settings_struct.rigidBody_flag_val = animation_rigidBody;
animation_settings_struct.animation_flag_val = animation_flag;
animation_settings_struct.movingView_flag_val = animation_movingView;

flag_input_struct = struct;
flag_input_struct.visualization_flag_val = visualization_flag;
flag_input_struct.ground_contact_visualization_flag_val = ground_contact_visualization_flag;
flag_input_struct.segment_visualize_bool_val = segment_visualize_bool;
flag_input_struct.all_subject_flag_val = all_subjects;
flag_input_struct.movingView_flag_val = animation_movingView;
flag_input_struct.rawData_plot_flag_val = plot_flag.rawData;
flag_input_struct.processedData_plot_flag_val = plot_flag.processedData;
flag_input_struct.rigidBody_plot_flag_val = plot_flag.rigidBody;
flag_input_struct.rawData_animation_flag_val = animation_rawData;
flag_input_struct.rigidBody_animation_flag_val = animation_rigidBody;
flag_input_struct.ground_contact_rawData_val = ground_contact_rawData;

pose_estimate_anatomical_frame = array_pose_estimate_anatomical_frame;

%% COMPUTATIONS PART

subject_reference_frame_data_mvt = struct;
num_subjects = 25;

if all_subjects
    for s=1:num_subjects
        disp(s);
        input_data_val.subject_idx_val = s;
        meta_data_val.subject_idx_val = s;
        w = waitbar(0, 'Generating Reference Frame Data...');
        [local_vectors_output, marker_output] = generate_solidification_data(input_data_val, activity_reprocessed);
        input_data_val.local_vectors_output_val = local_vectors_output;
        input_data_val.marker_output_val = marker_output;
        reference_frame_data_current.raw = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame, raw_data_bool, rigid_body_bool);
        reference_frame_data_current.processed = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool, rigid_body_bool);
        reference_frame_data_current.rigidBody = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool, ~rigid_body_bool);
        waitbar(0.2, w, 'Generating Figure Scale Values...');
        fig_scale_val = generate_fig_scale_info(reference_frame_data_current.raw);
        waitbar(0.33, w, 'Regrouping Reference Frame Data...');
        pose_data_current.raw = regroup_pose_data(reference_frame_data_current.raw, segment_list);
        pose_data_current.processed = regroup_pose_data(reference_frame_data_current.processed, segment_list);
        pose_data_current.rigidBody = regroup_pose_data(reference_frame_data_current.rigidBody, segment_list);
        waitbar(0.67, w, 'Estimating Knee Angle...');
        knee_angle_data_current.raw = estimate_knee_angle_gs(reference_frame_data_current.raw);
        knee_angle_data_current.processed = estimate_knee_angle_gs(reference_frame_data_current.processed);
        knee_angle_data_current.rigidBody = estimate_knee_angle_gs(reference_frame_data_current.rigidBody);
        time_vector = extract_time_vector(reference_frame_data_current.raw);
        meta_data_val.time_vector_val = time_vector;
        waitbar(0.8, w, 'Generating Visualizations...');
        if activityTypeOverGdWalking
            if ground_contact_rawData
                meta_data_val.sacral_data_val = activity(activity_idx).VICON(meta_data_val.subject_idx_val, default_trial_idx).movement_data.sacrum_cluster1_x;
                direction_vector_val = estimate_direction(meta_data_val);
                ground_contact_data = extract_ground_contact_data(activity, meta_data_val, direction_vector_val);
            else
                meta_data_val.sacral_data_val = activity_reprocessed(activity_idx).VICON(meta_data_val.subject_idx_val, default_trial_idx).movement_data.sacrum_cluster1_x;
                direction_vector_val = estimate_direction(meta_data_val);
                ground_contact_data = extract_ground_contact_data(activity_reprocessed, meta_data_val, direction_vector_val);
            end
        end
        generate_plot_pose_components(plot_flag, time_vector, pose_data_current, s, segment_list);
        generate_plot_knee_angle(plot_flag, time_vector, knee_angle_data_current, s);
        create_animation(animation_settings_struct, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
        create_ground_contact_visualization(flag_input_struct, reference_frame_data_current, meta_data_val, fig_scale_val, ground_contact_data);
        create_visualization(flag_input_struct, reference_frame_data_current, meta_data_val, fig_scale_val, knee_angle_data_current, ground_contact_data);
        waitbar(0.9, w, 'Storing Data...');
        subject_reference_frame_data_mvt(s).reference_frame_data.raw = reference_frame_data_current.raw;
        subject_reference_frame_data_mvt(s).knee_angle_data.raw = knee_angle_data_current.raw;
        subject_reference_frame_data_mvt(s).reference_frame_data.processed = reference_frame_data_current.processed;
        subject_reference_frame_data_mvt(s).knee_angle_data.processed = knee_angle_data_current.processed;
        subject_reference_frame_data_mvt(s).reference_frame_data.rigidBody = reference_frame_data_current.rigidBody;
        subject_reference_frame_data_mvt(s).knee_angle_data.rigidBody = knee_angle_data_current.rigidBody;
        waitbar(1, w, 'Finishing...');
        close(w);
    end
end

if one_subject
    s = subject_idx_val;
    input_data_val.subject_idx_val = s;
    meta_data_val.subject_idx_val = s;
    w = waitbar(0, 'Generating Reference Frame Data...');
    [local_vectors_output, marker_output] = generate_solidification_data(input_data_val, activity_reprocessed);
    input_data_val.local_vectors_output_val = local_vectors_output;
    input_data_val.marker_output_val = marker_output;
    reference_frame_data_current.raw = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame, raw_data_bool, rigid_body_bool);
    reference_frame_data_current.processed = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool, rigid_body_bool);
    reference_frame_data_current.rigidBody = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool, ~rigid_body_bool);
    waitbar(0.2, w, 'Generating Figure Scale Values...');
    fig_scale_val = generate_fig_scale_info(reference_frame_data_current.raw);
    waitbar(0.33, w, 'Regrouping Reference Frame Data...');
    pose_data_current.raw = regroup_pose_data(reference_frame_data_current.raw, segment_list);
    pose_data_current.processed = regroup_pose_data(reference_frame_data_current.processed, segment_list);
    pose_data_current.rigidBody = regroup_pose_data(reference_frame_data_current.rigidBody, segment_list);
    waitbar(0.67, w, 'Estimating Knee Angle...');
    knee_angle_data_current.raw = estimate_knee_angle_gs(reference_frame_data_current.raw);
    knee_angle_data_current.processed = estimate_knee_angle_gs(reference_frame_data_current.processed);
    knee_angle_data_current.rigidBody = estimate_knee_angle_gs(reference_frame_data_current.rigidBody);
    time_vector = extract_time_vector(reference_frame_data_current.raw);
    meta_data_val.time_vector_val = time_vector;
    waitbar(0.8, w, 'Generating Visualizations...');
    if activityTypeOverGdWalking
        if ground_contact_rawData
            meta_data_val.sacral_data_val = activity(activity_idx).VICON(meta_data_val.subject_idx_val, default_trial_idx).movement_data.sacrum_cluster1_x;
            direction_vector_val = estimate_direction(meta_data_val);
            ground_contact_data = extract_ground_contact_data(activity, meta_data_val, direction_vector_val);
        else
            meta_data_val.sacral_data_val = activity_reprocessed(activity_idx).VICON(meta_data_val.subject_idx_val, default_trial_idx).movement_data.sacrum_cluster1_x;
            direction_vector_val = estimate_direction(meta_data_val);
            ground_contact_data = extract_ground_contact_data(activity_reprocessed, meta_data_val, direction_vector_val);
        end
    end
    generate_plot_pose_components(plot_flag, time_vector, pose_data_current, s, segment_list);
    generate_plot_knee_angle(plot_flag, time_vector, knee_angle_data_current, s);
    create_animation(animation_settings_struct, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
    create_ground_contact_visualization(flag_input_struct, reference_frame_data_current, meta_data_val, fig_scale_val, ground_contact_data);
    create_visualization(flag_input_struct, reference_frame_data_current, meta_data_val, fig_scale_val, knee_angle_data_current, ground_contact_data);
    waitbar(0.9, w, 'Storing Data...');
    subject_reference_frame_data_mvt(s).reference_frame_data.raw = reference_frame_data_current.raw;
    subject_reference_frame_data_mvt(s).knee_angle_data.raw = knee_angle_data_current.raw;
    subject_reference_frame_data_mvt(s).reference_frame_data.processed = reference_frame_data_current.processed;
    subject_reference_frame_data_mvt(s).knee_angle_data.processed = knee_angle_data_current.processed;
    subject_reference_frame_data_mvt(s).reference_frame_data.rigidBody = reference_frame_data_current.rigidBody;
    subject_reference_frame_data_mvt(s).knee_angle_data.rigidBody = knee_angle_data_current.rigidBody;
    waitbar(1, w, 'Finishing...');
    close(w);
end

save('reference_frame_data_mvt_file', 'subject_reference_frame_data_mvt', '-v7.3');

function [r1_now_val, r2_now_val, r3_now_val] = current_technical_marker_position(subject_idx, activity_idx, current_segment_val, current_side_val, t_now_idx, technical_frame_markers_info, activity)
    default_trial_idx = 1;
    r1_now_val = nan(3,1);
    r2_now_val = nan(3,1);
    r3_now_val = nan(3,1);
    technical_markers_absent = all(technical_frame_markers_info(subject_idx).(current_segment_val).(current_side_val) == "");
    all_marker_name = technical_frame_markers_info(subject_idx).(current_segment_val).(current_side_val);
    
    if ~technical_markers_absent
        r_now = struct;
            for k=1:3
                current_marker_name = all_marker_name(k);
                x_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_marker_name + "_x")(t_now_idx);
                y_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_marker_name + "_y")(t_now_idx);
                z_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_marker_name + "_z")(t_now_idx);
                r_now_k = [x_now_k; y_now_k; z_now_k];
                r_now(k).value = r_now_k;
            end
        r1_now_val = r_now(1).value;
        r2_now_val = r_now(2).value;
        r3_now_val = r_now(3).value;
    end
end

function [local_vectors_output, marker_output] = generate_solidification_data(input_data, activity)
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val.rigidBody;
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

function [T,R,unit_vectors,success_bool] = pose_estimate_technical_frame(ptCloud)
    pts_collinear = true;
    success_bool = false;

    if ~any(any(isnan(ptCloud.Location)))
        ptCloud_tp = pointCloud(transpose(ptCloud.Location));
        [~,plane_inlierIndices,~] = pcfitplane(ptCloud_tp,1);
        triangle_inequality_satisfied_bool = apply_triangle_inequality(ptCloud);
        if or(~isempty(plane_inlierIndices), triangle_inequality_satisfied_bool)
            success_bool = true;
            pts_collinear = false;
        end
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];
    
    r1 = ptCloud.Location(:,1);

    T = r1;
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear
        
        r2 = ptCloud.Location(:,2);
        r3 = ptCloud.Location(:,3);

        r21 = r2 - r1;
        r31 = r3 - r1;

        u21 = r21/norm(r21);
        u31 = r31/norm(r31);

        i_hat_dir = u21;
        j_hat_dir = cross(u21, u31);
        k_hat_dir = cross(i_hat_dir, j_hat_dir);
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    if ~any(any(isnan(R_calc)))
        same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
        det_R_bool = abs(det(R_calc)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            R = R_calc;
        end
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
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

function technical_frame_plotter(origin_local, unit_vector_local)
    tf_axis_size = 50;
    T = origin_local;
    i_hat = tf_axis_size*unit_vector_local.x;
    j_hat = tf_axis_size*unit_vector_local.y;
    k_hat = tf_axis_size*unit_vector_local.z;
    plot3([T(1) T(1)+i_hat(1)],[-T(3) -T(3)-i_hat(3)],[T(2) T(2)+i_hat(2)], 'r', 'LineWidth',1)
    plot3([T(1) T(1)+j_hat(1)],[-T(3) -T(3)-j_hat(3)],[T(2) T(2)+j_hat(2)], 'g', 'LineWidth',1)
    plot3([T(1) T(1)+k_hat(1)],[-T(3) -T(3)-k_hat(3)],[T(2) T(2)+k_hat(2)], 'b', 'LineWidth',1)
end

function anatomical_frame_plotter(origin_local, unit_vector_local)
    af_axis_size = 60;
    T = origin_local;
    i_hat = af_axis_size*unit_vector_local.x;
    j_hat = af_axis_size*unit_vector_local.y;
    k_hat = af_axis_size*unit_vector_local.z;
    plot3([T(1) T(1)+i_hat(1)],[-T(3) -T(3)-i_hat(3)],[T(2) T(2)+i_hat(2)], 'r', 'LineWidth',1.5)
    plot3([T(1) T(1)+j_hat(1)],[-T(3) -T(3)-j_hat(3)],[T(2) T(2)+j_hat(2)], 'g', 'LineWidth',1.5)
    plot3([T(1) T(1)+k_hat(1)],[-T(3) -T(3)-k_hat(3)],[T(2) T(2)+k_hat(2)], 'b', 'LineWidth',1.5)
end

function global_frame_plotter
    axis_size = 100;
    i = [axis_size;0;0];
    j = [0;axis_size;0];
    k = [0;0;axis_size];
    plot3([0 i(1)],[0 -i(3)],[0 i(2)], 'r', 'LineWidth',2)
    hold on
    plot3([0 j(1)],[0 -j(3)],[0 j(2)], 'g', 'LineWidth',2)
    plot3([0 k(1)],[0 -k(3)],[0 k(2)], 'b', 'LineWidth',2)
end

function pretty_plotter(fig_scale, subject_idx)
    grid

    description_axes = zeros(3, 1);
    description_axes(1) = plot(NaN,NaN,'-r');
    description_axes(2) = plot(NaN,NaN,'-g');
    description_axes(3) = plot(NaN,NaN,'-b');
    description_axes(4) = plot(NaN,NaN,'-m');
    description_axes(5) = plot(NaN,NaN,'-y');
    legend(description_axes, 'x: red','y: green','z: blue', 'technical segment: magenta', 'reconstructed segment: yellow');
    if ~isnan(fig_scale.x_min) && ~isnan(fig_scale.x_max)
        xlim([fig_scale.x_min, fig_scale.x_max])
    end
    if ~isnan(fig_scale.y_min) && ~isnan(fig_scale.y_max)
        ylim([fig_scale.y_min, fig_scale.y_max])
    end
    if ~isnan(fig_scale.z_min) && ~isnan(fig_scale.z_max)
        zlim([fig_scale.z_min, fig_scale.z_max])
    end
    if all(~isnan(fig_scale.view_vector(subject_idx,:)))
        view(fig_scale.view_vector(subject_idx,:));
    end
end

function reference_frame_data = reference_frame_data_generate(input_data, activity, pose_estimate_anatomical_frame, raw_data_bool_val, rigid_body_bool_val)
    reference_frame_data = struct;
    default_trial_idx = 1;
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    anatomical_markers_info = input_data.anatomical_markers_info_val;
    if raw_data_bool_val
        anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.raw;
        technical_frame_markers_info = input_data.technical_frame_markers_info_val.raw;
    else
        if rigid_body_bool_val
            anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.rigidBody;
            technical_frame_markers_info = input_data.technical_frame_markers_info_val.rigidBody;
        else
            anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.processed;
            technical_frame_markers_info = input_data.technical_frame_markers_info_val.processed;
        end
    end
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    anatomical_frame_function_idx_list = input_data.anatomical_frame_function_idx_list;
    t = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.time;
    num_measurements = length(t);
    num_sides = length(side_list);
    num_segments = length(segment_list);
    for n=1:num_measurements
        nan_bool_val = false;
        technical_frames = struct;
        anatomical_frames = struct;
        segment_visualize = struct;
        segment_visualize.technical = struct;
        segment_visualize.reconstructed = struct;
        for j=1:num_segments
            current_segment = segment_list(j);
            current_anatomical_frame_function_idx = anatomical_frame_function_idx_list(j);
            for m=1:num_sides
                current_side = side_list(m);
                
                if rigid_body_bool_val
                    local_vectors_output = input_data.local_vectors_output_val;
                    marker_output = input_data.marker_output_val;
                    
                    marker_input_struct = struct;
                    marker_input_struct.technical_markers_list_val = marker_output.(current_segment).(current_side).technical_markers_list_val;
                    marker_input_struct.marker_info_val = marker_output.(current_segment).(current_side).marker_info_val;
                    marker_input_struct.sorted_vertex_val = marker_output.(current_segment).(current_side).sorted_vertex_val;
                    marker_input_struct.success_bool_val = marker_output.(current_segment).(current_side).success_bool_val;

                    local_vectors_input_struct = struct;
                    local_vectors_input_struct.local_vectors_val = local_vectors_output.(current_segment).(current_side).local_vectors_val;
                    local_vectors_input_struct.local_vectors_mean_val = local_vectors_output.(current_segment).(current_side).local_vectors_mean_val;
                    local_vectors_input_struct.q_local_val = local_vectors_output.(current_segment).(current_side).q_local_val;
                    local_vectors_input_struct.solid_triangle_val = local_vectors_output.(current_segment).(current_side).solid_triangle_val;

                    [T, R, unit_vector_hat, success_bool, reconstructed_vectors] = pose_estimate_technical_frame_rigid_body(n, marker_input_struct, local_vectors_input_struct);
                    r1_now = reconstructed_vectors(1).vecTrue;
                    r2_now = reconstructed_vectors(2).vecTrue;
                    r3_now = reconstructed_vectors(3).vecTrue;
                    current_ptCloud_val = pointCloud(horzcat(r1_now,r2_now,r3_now));
                else
                    [r1_now, r2_now, r3_now] = current_technical_marker_position(subject_idx, activity_idx, current_segment, current_side, n, technical_frame_markers_info, activity);
                    current_ptCloud_val = pointCloud(horzcat(r1_now,r2_now,r3_now));
                    [T, R, unit_vector_hat, success_bool] = pose_estimate_technical_frame(current_ptCloud_val);
                end
                technical_frames.(current_segment).(current_side).ptCloud_val = current_ptCloud_val;
                
                if ~success_bool
                    input_data_nan = any(any(isnan([r1_now, r2_now, r3_now])));
                    if input_data_nan
                        disp("input data has at least one NaN value")
                        nan_bool_val = true;
                    else
                        disp("NO NaN values------------------------")
                    end
                end
                technical_frames.(current_segment).(current_side).T_val = T;
                technical_frames.(current_segment).(current_side).R_val = R;
                technical_frames.(current_segment).(current_side).unit_vector_hat_val = unit_vector_hat;
                trial_pt = struct;
                trial_pt.technical = struct;
                X = nan(3,1);
                Y = nan(3,1);
                Z = nan(3,1);
                V = nan(3,3);
                trial_pt.technical(1).global = r1_now;
                trial_pt.technical(2).global = r2_now;
                trial_pt.technical(3).global = r3_now;
                num_trial_pt.technical = length(trial_pt.technical);
                for k=1:num_trial_pt.technical
                    X(k) = trial_pt.technical(k).global(1);
                    Y(k) = trial_pt.technical(k).global(2);
                    Z(k) = trial_pt.technical(k).global(3);
                    V(:,k) = [X(k); Y(k); Z(k)];
                end
                segment_visualize.technical.(current_segment).(current_side).X_val = X;
                segment_visualize.technical.(current_segment).(current_side).Y_val = Y;
                segment_visualize.technical.(current_segment).(current_side).Z_val = Z;
                current_anatomical_calibration_parameters = anatomical_calibration_parameters(subject_idx).(current_segment).(current_side);
                trial_pt.anatomical = struct;
                Xr = nan(3,1);
                Yr = nan(3,1);
                Zr = nan(3,1);
                Vr = nan(3,3);
                current_reconstructed_position = struct;
                current_reconstructed_position.side = current_side;
                num_trial_pt.anatomical = length(anatomical_markers_info(subject_idx).(current_segment).(current_side).names);
                current_anatomical_marker_name_list = anatomical_markers_info(subject_idx).(current_segment).(current_side).names;
                for k=1:num_trial_pt.anatomical
                    current_anatomical_marker_name = current_anatomical_marker_name_list(k);
                    xa_now_local = current_anatomical_calibration_parameters(k).x;
                    ya_now_local = current_anatomical_calibration_parameters(k).y;
                    za_now_local = current_anatomical_calibration_parameters(k).z;
                    ra_now_local = [xa_now_local; ya_now_local; za_now_local];
                    trial_pt.anatomical(k).local = ra_now_local;
                    trial_pt.anatomical(k).reconstructed = local2global(T, R, ra_now_local);
                    Xr(k) = trial_pt.anatomical(k).reconstructed(1);
                    Yr(k) = trial_pt.anatomical(k).reconstructed(2);
                    Zr(k) = trial_pt.anatomical(k).reconstructed(3);
                    Vr(:,k) = [Xr(k); Yr(k); Zr(k)];
                    current_reconstructed_position.(current_anatomical_marker_name) = Vr(:,k);
                end
                segment_visualize.reconstructed.(current_segment).(current_side).X_val = Xr;
                segment_visualize.reconstructed.(current_segment).(current_side).Y_val = Yr;
                segment_visualize.reconstructed.(current_segment).(current_side).Z_val = Zr;
                [T_af, R_af, unit_vector_hat_af] = pose_estimate_anatomical_frame{current_anatomical_frame_function_idx}(current_reconstructed_position);
                anatomical_frames.(current_segment).(current_side).T_val = T_af;
                anatomical_frames.(current_segment).(current_side).R_val = R_af;
                anatomical_frames.(current_segment).(current_side).unit_vector_hat_val = unit_vector_hat_af;
            end
        end
        reference_frame_data(n).nan_bool = nan_bool_val;
        reference_frame_data(n).time = t(n);
        reference_frame_data(n).technical_frames = technical_frames;
        reference_frame_data(n).anatomical_frames = anatomical_frames;
        reference_frame_data(n).segment_visualize = segment_visualize;
    end
    disp(sum([reference_frame_data(:).nan_bool]))
end

function create_animation(animation_settings, segment_visualize_bool, reference_frame_data_all, meta_data, fig_scale, all_subject_flag)
    
    rawData_flag = animation_settings.rawData_flag_val;
    rigidBody_flag = animation_settings.rigidBody_flag_val;
    animation_flag = animation_settings.animation_flag_val;
    movingView_flag = animation_settings.movingView_flag_val;
    
    if rawData_flag
        reference_frame_data = reference_frame_data_all.raw;
    else
        if rigidBody_flag
            reference_frame_data = reference_frame_data_all.rigidBody;
        else
            reference_frame_data = reference_frame_data_all.processed;
        end
    end
    
    if animation_flag
        
        time_vector_val = meta_data.time_vector_val;
        
        subject_idx = meta_data.subject_idx_val;
        side_list = meta_data.side_list_val;
        segment_list = meta_data.segment_list_val;
        
        num_measurements = length(reference_frame_data);
        num_sides = length(side_list);
        num_segments = length(segment_list);
        
        original_view_vector = fig_scale.view_vector(subject_idx,:);
        original_azimuth = original_view_vector(1);
        original_elevation = original_view_vector(2);
        dTheta = 0;
        if ~isnan(num_measurements)
            if num_measurements > 1
                dTheta = rad2deg((2*pi)/(num_measurements-1));
            end
        end
        
        folder_name_start = "Movement_Videos\movement_subNo_";
        if movingView_flag
            folder_name_start = "Movement_MovingPtView_Videos\movement_subNo_";
        end
        
        my_writer = VideoWriter(folder_name_start + string(subject_idx), 'MPEG-4');
        my_writer.FrameRate = 50;
        open(my_writer);
        
        figh = figure;
        figh.WindowState = 'maximized';

        for k=1:num_measurements
            clf

            technical_frames = reference_frame_data(k).technical_frames;
            anatomical_frames = reference_frame_data(k).anatomical_frames;
            segment_visualize = reference_frame_data(k).segment_visualize;

            global_frame_plotter
            for j=1:num_segments
                current_segment = segment_list(j);
                for m=1:num_sides
                    current_side = side_list(m);
                    if segment_visualize_bool.technical
                        fill3(segment_visualize.technical.(current_segment).(current_side).X_val, -segment_visualize.technical.(current_segment).(current_side).Z_val, segment_visualize.technical.(current_segment).(current_side).Y_val,'m','FaceAlpha',.2)
                        technical_frame_plotter(technical_frames.(current_segment).(current_side).T_val, technical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                    if segment_visualize_bool.reconstructed
                        fill3(segment_visualize.reconstructed.(current_segment).(current_side).X_val, -segment_visualize.reconstructed.(current_segment).(current_side).Z_val, segment_visualize.reconstructed.(current_segment).(current_side).Y_val,'y','FaceAlpha',.2)
                        anatomical_frame_plotter(anatomical_frames.(current_segment).(current_side).T_val, anatomical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                end
            end
            if movingView_flag
                new_azimuth = original_azimuth + (k-1)*dTheta;
                fig_scale.view_vector(subject_idx,:) = [new_azimuth original_elevation];
            end
            
            if ~all(isnan(time_vector_val))
                pretty_plotter(fig_scale, subject_idx)
            end
           
            current_frame = getframe(figh);
            writeVideo(my_writer, current_frame);
        end
        close(my_writer);
        
        if all_subject_flag
            close(figh);
        end

    end

end

function fig_scale_info = generate_fig_scale_info(reference_frame_data_current_val)
    plot_edge_width = 500;
    ground_height = 0;
    ceiling_height = 1200;
    num_measurements = length(reference_frame_data_current_val);
    tibial_x_right_list = nan(num_measurements,1);
    tibial_x_left_list = nan(num_measurements,1);
    tibial_y_right_list = nan(num_measurements,1);
    tibial_y_left_list = nan(num_measurements,1);
    tibial_z_right_list = nan(num_measurements,1);
    tibial_z_left_list = nan(num_measurements,1);
    femoral_x_right_list = nan(num_measurements,1);
    femoral_x_left_list = nan(num_measurements,1);
    femoral_y_right_list = nan(num_measurements,1);
    femoral_y_left_list = nan(num_measurements,1);
    femoral_z_right_list = nan(num_measurements,1);
    femoral_z_left_list = nan(num_measurements,1);
    for n=1:num_measurements
        current_tibial_pos_right = reference_frame_data_current_val(n).anatomical_frames.tibial.right.T_val;
        current_tibial_pos_left = reference_frame_data_current_val(n).anatomical_frames.tibial.left.T_val;
        current_femoral_pos_right = reference_frame_data_current_val(n).anatomical_frames.femoral.right.T_val;
        current_femoral_pos_left = reference_frame_data_current_val(n).anatomical_frames.femoral.left.T_val;
        tibial_x_right_list(n) = current_tibial_pos_right(1);
        tibial_x_left_list(n) = current_tibial_pos_left(1);
        tibial_y_right_list(n) = current_tibial_pos_right(2);
        tibial_y_left_list(n) = current_tibial_pos_left(2);
        tibial_z_right_list(n) = current_tibial_pos_right(3);
        tibial_z_left_list(n) = current_tibial_pos_left(3);
        femoral_x_right_list(n) = current_femoral_pos_right(1);
        femoral_x_left_list(n) = current_femoral_pos_left(1);
        femoral_y_right_list(n) = current_femoral_pos_right(2);
        femoral_y_left_list(n) = current_femoral_pos_left(2);
        femoral_z_right_list(n) = current_femoral_pos_right(3);
        femoral_z_left_list(n) = current_femoral_pos_left(3);
    end
    
    x_max_limit = max([max(femoral_x_right_list), max(tibial_x_right_list), max(femoral_x_left_list), max(tibial_x_left_list)]) + plot_edge_width;
    y_max_limit = ceiling_height;
    z_max_limit = max([max(femoral_z_right_list), max(tibial_z_right_list), max(femoral_z_left_list), max(tibial_z_left_list)]) + plot_edge_width;

    x_min_limit = min([min(femoral_x_right_list), min(tibial_x_right_list), min(femoral_x_left_list), min(tibial_x_left_list)]) - plot_edge_width;
    y_min_limit = ground_height;
    z_min_limit = min([min(femoral_z_right_list), min(tibial_z_right_list), min(femoral_z_left_list), min(tibial_z_left_list)]) - plot_edge_width;

    fig_scale_info = struct;
    
    fig_scale_info.x_max = x_max_limit;
    fig_scale_info.y_max = -z_min_limit;
    fig_scale_info.z_max = y_max_limit;

    fig_scale_info.x_min = x_min_limit;
    fig_scale_info.y_min = -z_max_limit;
    fig_scale_info.z_min = y_min_limit;
    
    view_vector_rep = [-37.5 30];
    view_vector = repmat(view_vector_rep,25,1);
    
    fig_scale_info.view_vector = view_vector;
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

function knee_angle = estimate_knee_angle_cappozzo(reference_frame_data_struct)
    knee_angle = struct;
    side_list = ["right", "left"];
    num_sides = length(side_list);
    num_measurements = length(reference_frame_data_struct);
    for j=1:num_measurements
        for k=1:num_sides
            current_side = side_list(k);
            Rp = reference_frame_data_struct(1).anatomical_frames.femoral.(current_side).R_val;
            Rd = reference_frame_data_struct(1).anatomical_frames.tibial.(current_side).R_val;
            Rj = Rp'*Rd;
            alpha_val = rad2deg(asin(Rj(3,2)));
            beta_val = rad2deg(asin(-Rj(3,1)/cos(alpha_val)));
            gamma_val = rad2deg(asin(-Rj(1,2)/cos(alpha_val)));
            knee_angle(j).(current_side).alpha = alpha_val;
            knee_angle(j).(current_side).beta = beta_val;
            knee_angle(j).(current_side).gamma = gamma_val;
        end
    end
end

function knee_angle = estimate_knee_angle_gs(reference_frame_data_struct)
    knee_angle = struct;
    side_list = ["right", "left"];
    num_sides = length(side_list);
    num_measurements = length(reference_frame_data_struct);
    for m=1:num_sides
        current_side = side_list(m);
        flexion_list = nan(num_measurements,1);
        adduction_list = nan(num_measurements,1);
        externalRotation_list = nan(num_measurements,1);
        for n=1:num_measurements
            
            i = reference_frame_data_struct(n).anatomical_frames.tibial.(current_side).unit_vector_hat_val.x;
            j = reference_frame_data_struct(n).anatomical_frames.tibial.(current_side).unit_vector_hat_val.y;
            k = reference_frame_data_struct(n).anatomical_frames.tibial.(current_side).unit_vector_hat_val.z;

            I = reference_frame_data_struct(n).anatomical_frames.femoral.(current_side).unit_vector_hat_val.x;
            J = reference_frame_data_struct(n).anatomical_frames.femoral.(current_side).unit_vector_hat_val.y;
            K = reference_frame_data_struct(n).anatomical_frames.femoral.(current_side).unit_vector_hat_val.z;
            
            I_gs = K;
            J_gs = I;
            K_gs = J;
            
            i_gs = k;
            j_gs = i;
            k_gs = j;

            e1 = I_gs;
            e3 = k_gs;
            e2 = cross(e3,e1);

            e1_hat = e1/norm(e1);
            e2_hat = e2/norm(e2);
            e3_hat = e3/norm(e3);

            alpha = rad2deg(asin(-dot(e2_hat,K_gs)));
            beta = rad2deg(acos(dot(I_gs,k_gs)));
            gamma = rad2deg(asin(dot(e2_hat,i_gs)));
            
            current_flexion = alpha;
            if current_side == "right"
                current_adduction = beta - 90;
                current_externalRotation = -gamma;
            else
                current_adduction = 90 - beta;
                current_externalRotation = gamma;
            end
            
            flexion_list(n) = current_flexion;
            adduction_list(n) = current_adduction;
            externalRotation_list(n) = current_externalRotation;
        end
        knee_angle.(current_side).flexion = flexion_list;
        knee_angle.(current_side).adduction = adduction_list;
        knee_angle.(current_side).externalRotation = externalRotation_list;
    end
end

function generate_plot_knee_angle(plot_flag, time_vector_val, knee_angle_struct, subject_idx)
    if plot_flag.knee_angle
        if ~all(isnan(time_vector_val))
            fn = string(fieldnames(knee_angle_struct.raw.right));
            num_subplot = length(fn);
            plt_width = 1;
            for j=1:num_subplot
                current_fn = fn(j);
                right_qt.raw = knee_angle_struct.raw.right.(current_fn);
                left_qt.raw = knee_angle_struct.raw.left.(current_fn);
                right_qt.processed = knee_angle_struct.processed.right.(current_fn);
                left_qt.processed = knee_angle_struct.processed.left.(current_fn);
                right_qt.rigidBody = knee_angle_struct.rigidBody.right.(current_fn);
                left_qt.rigidBody = knee_angle_struct.rigidBody.left.(current_fn);
                if(~isempty(right_qt.raw) && ~isempty(left_qt.raw))
                    right_qt.raw_mean = round(mean(right_qt.raw, 'omitnan'),3);
                    left_qt.raw_mean = round(mean(left_qt.raw, 'omitnan'),3);
                    right_qt.raw_std = round(std(right_qt.raw, 'omitnan'),3);
                    left_qt.raw_std = round(std(left_qt.raw, 'omitnan'),3);
                    max_ht_qt = max([max(right_qt.raw), max(left_qt.raw), max(right_qt.processed), max(left_qt.processed), max(right_qt.rigidBody), max(left_qt.rigidBody)], [], 'omitnan') + plt_width;
                    min_ht_qt = min([min(right_qt.raw), min(left_qt.raw), min(right_qt.processed), min(left_qt.processed), min(right_qt.rigidBody), min(left_qt.rigidBody)], [], 'omitnan') - plt_width;
                    title_subplot = current_fn + " [right raw (avg,sd) = (" + string(right_qt.raw_mean) + "," + string(right_qt.raw_std) + "), left raw (avg,sd) = (" + string(left_qt.raw_mean) + "," + string(left_qt.raw_std) + ") ]";
                    subplot(num_subplot,1,j)
                    if plot_flag.rawData
                        plot(time_vector_val, right_qt.raw, 'r-.', 'DisplayName', 'right raw')
                        hold on
                        plot(time_vector_val, left_qt.raw, 'b-.', 'DisplayName', 'left raw')
                        if plot_flag.meanValue
                            if ~isnan(right_qt.raw_mean)
                                yline(right_qt.raw_mean, 'r--', 'DisplayName', 'right raw avg')
                            end
                            if ~isnan(left_qt.raw_mean)
                                yline(left_qt.raw_mean, 'b--', 'DisplayName', 'left raw avg')
                            end
                            plt_limit = [min_ht_qt, max_ht_qt];
                            if all(~isnan(plt_limit))
                                ylim(plt_limit)
                            end
                        end
                        title(title_subplot)
                    end
                    if plot_flag.processedData
                        plot(time_vector_val, right_qt.processed, 'r', 'DisplayName', 'right processed')
                        hold on
                        plot(time_vector_val, left_qt.processed, 'b', 'DisplayName', 'left processed')
                    end
                    if plot_flag.rigidBody
                        plot(time_vector_val, right_qt.rigidBody, 'Color', '#8B0000', 'DisplayName', 'right rigidBody')
                        hold on
                        plot(time_vector_val, left_qt.rigidBody, 'Color', '#8A2BE2', 'DisplayName', 'left rigidBody')
                    end
                    grid
                    legend
                    xlabel("time (s)")
                    ylabel("angle (degrees)")
                end
            end
            suptitle("Knee Angle Plots")
            fig = gcf;
            fig.WindowState = 'maximized';
            saveas(fig,"MVT_Knee_Angle_Plots\knee_angle_subNo_" + string(subject_idx), "jpg");
            close(fig);
        end
    end
end

function generate_plot_pose_components(plot_flag, time_vector_val, pose_struct, subject_idx, segment_list)
    if plot_flag.pose_components
        if ~all(isnan(time_vector_val))
            num_segments = length(segment_list);
            fn = ["T", "i", "j", "k"];
            side_list = ["left", "right"];
            row_subplot = length(fn);
            col_subplot = length(side_list);
            for s=1:num_segments
                current_segment = segment_list(s);
                for n=1:row_subplot
                    for m=1:col_subplot
                        current_side = side_list(m);
                        current_fn = fn(n);
                        current_vector.raw = pose_struct.raw.(current_segment).(current_side).(current_fn);
                        current_vector.processed = pose_struct.processed.(current_segment).(current_side).(current_fn);
                        current_vector.rigidBody = pose_struct.rigidBody.(current_segment).(current_side).(current_fn);
                        x_qt.raw = current_vector.raw.x;
                        y_qt.raw = current_vector.raw.y;
                        z_qt.raw = current_vector.raw.z;
                        x_qt.processed = current_vector.processed.x;
                        y_qt.processed = current_vector.processed.y;
                        z_qt.processed = current_vector.processed.z;
                        x_qt.rigidBody = current_vector.rigidBody.x;
                        y_qt.rigidBody = current_vector.rigidBody.y;
                        z_qt.rigidBody = current_vector.rigidBody.z;
                        if(~isempty(x_qt.raw) && (~isempty(y_qt.raw) && ~isempty(z_qt.raw)))
                            x_qt.raw_mean = round(mean(x_qt.raw, 'omitnan'),3);
                            y_qt.raw_mean = round(mean(y_qt.raw, 'omitnan'),3);
                            z_qt.raw_mean = round(mean(z_qt.raw, 'omitnan'),3);
                            x_qt.raw_std = round(std(x_qt.raw, 'omitnan'),3);
                            y_qt.raw_std = round(std(y_qt.raw, 'omitnan'),3);
                            z_qt.raw_std = round(std(z_qt.raw, 'omitnan'),3);
                            max_ht_qt = max([max(x_qt.raw), max(y_qt.raw), max(z_qt.raw)], [], 'omitnan');
                            min_ht_qt = min([min(x_qt.raw), min(y_qt.raw), min(z_qt.raw)], [], 'omitnan');
                            plt_width_top = 0.4 * abs(max_ht_qt);
                            plt_width_bottom = 0.4 * abs(min_ht_qt);
                            subplot(row_subplot,col_subplot,(n-1)*col_subplot + m)
                            if plot_flag.rawData
                                plot(time_vector_val, x_qt.raw, 'r-.', 'DisplayName', 'x raw')
                                hold on
                                plot(time_vector_val, y_qt.raw, 'g-.', 'DisplayName', 'y raw')
                                plot(time_vector_val, z_qt.raw, 'b-.', 'DisplayName', 'z raw')
                                if plot_flag.meanValue
                                    if ~isnan(x_qt.raw_mean)
                                        yline(x_qt.raw_mean, 'r--', 'DisplayName', 'x raw avg')
                                    end
                                    if ~isnan(y_qt.raw_mean)
                                        yline(y_qt.raw_mean, 'g--', 'DisplayName', 'y raw avg')
                                    end
                                    if ~isnan(z_qt.raw_mean)
                                        yline(z_qt.raw_mean, 'b--', 'DisplayName', 'z raw avg')
                                    end
                                end
                                title_x = " [x raw (avg,sd) = (" + string(x_qt.raw_mean) + "," + string(x_qt.raw_std) + ")";
                                title_y = ", y raw (avg,sd) = (" + string(y_qt.raw_mean) + "," + string(y_qt.raw_std) + ")";
                                title_z = ", z raw (avg,sd) = (" + string(z_qt.raw_mean) + "," + string(z_qt.raw_std) + ")]";
                                title_subplot = current_fn + title_x + title_y + title_z;
                                title(title_subplot)
                            end
                            if plot_flag.processedData
                                plot(time_vector_val, x_qt.processed, 'r', 'DisplayName', 'x processed')
                                hold on
                                plot(time_vector_val, y_qt.processed, 'g', 'DisplayName', 'y processed')
                                plot(time_vector_val, z_qt.processed, 'b', 'DisplayName', 'z processed')
                            end
                            if plot_flag.rigidBody
                                plot(time_vector_val, x_qt.rigidBody, 'r', 'DisplayName', 'x rigidBody')
                                hold on
                                plot(time_vector_val, y_qt.rigidBody, 'g', 'DisplayName', 'y rigidBody')
                                plot(time_vector_val, z_qt.rigidBody, 'b', 'DisplayName', 'z rigidBody')
                            end
                            plt_lower = min_ht_qt - plt_width_bottom;
                            plt_upper = max_ht_qt + plt_width_top;
                            plt_limits = [plt_lower, plt_upper];
                            if all(~isnan(plt_limits))
                                ylim(plt_limits)
                            end
                            grid
                            legend
                            xlabel("time (s)")
                            ylabel("length (mm)")
                        end
                    end
                end
                suptitle("Pose Components Plots: " + current_segment)
                fig = gcf;
                fig.WindowState = 'maximized';
                saveas(fig,"MVT_Pose_Components_Plots\pose_components_" + current_segment + "_subNo_" + string(subject_idx), "jpg");
                close(fig);
            end
        end
    end
end

function time_vector_list = extract_time_vector(reference_frame_data_struct)
    num_measurements = length(reference_frame_data_struct);
    time_vector_list = nan(num_measurements,1);
    for k=1:num_measurements
        time_vector_list(k) = reference_frame_data_struct(k).time;
    end
end

function pose_struct = regroup_pose_data(reference_frame_data_struct, segment_list)
    fn = ["T", "i", "j", "k"];
    side_list = ["left", "right"];
    fn_pose = ["T_val", "unit_vector_hat_val.x", "unit_vector_hat_val.y", "unit_vector_hat_val.z"];
    num_vectors = length(fn_pose);
    num_sides = length(side_list);
    num_segments = length(segment_list);
    num_measurements = length(reference_frame_data_struct);
    pose_struct = struct;
    for m=1:num_vectors
        current_fn = fn(m);
        current_fn_pose = fn_pose(m);
        for s=1:num_segments
            current_segment = segment_list(s);
            for p=1:num_sides
                current_side = side_list(p);
                current_vector_x = nan(num_measurements,1);
                current_vector_y = nan(num_measurements,1);
                current_vector_z = nan(num_measurements,1);
                for n=1:num_measurements
                    if current_fn == "T"
                        vector_from_struct = reference_frame_data_struct(n).anatomical_frames.(current_segment).(current_side).(current_fn_pose);
                    else
                        fn_pose_split = strsplit(current_fn_pose, ".");
                        vector_from_struct = reference_frame_data_struct(n).anatomical_frames.(current_segment).(current_side).(fn_pose_split(1)).(fn_pose_split(2));
                    end
                    current_vector_x(n) = vector_from_struct(1);
                    current_vector_y(n) = vector_from_struct(2);
                    current_vector_z(n) = vector_from_struct(3);
                end
                pose_struct.(current_segment).(current_side).(current_fn).x = current_vector_x;
                pose_struct.(current_segment).(current_side).(current_fn).y = current_vector_y;
                pose_struct.(current_segment).(current_side).(current_fn).z = current_vector_z;
            end
        end
    end
end

function create_visualization(flag_input, reference_frame_data_all, meta_data, fig_scale, knee_angle_struct, ground_contact_data)

    visualization_flag = flag_input.visualization_flag_val;
    segment_visualize_bool = flag_input.segment_visualize_bool_val;
    all_subject_flag = flag_input.all_subject_flag_val;
    movingView_flag = flag_input.movingView_flag_val;
    
    rawData_plot_flag = flag_input.rawData_plot_flag_val;
    processedData_plot_flag = flag_input.processedData_plot_flag_val;
    rigidBody_plot_flag = flag_input.rigidBody_plot_flag_val;
    
    rawData_flag = flag_input.rawData_animation_flag_val;
    rigidBody_flag = flag_input.rigidBody_animation_flag_val;
    
    if rawData_flag
        reference_frame_data = reference_frame_data_all.raw;
    else
        if rigidBody_flag
            reference_frame_data = reference_frame_data_all.rigidBody;
        else
            reference_frame_data = reference_frame_data_all.processed;
        end
    end

    if visualization_flag
        
        subject_idx = meta_data.subject_idx_val;
        side_list = meta_data.side_list_val;
        segment_list = meta_data.segment_list_val;
        
        time_vector_val = meta_data.time_vector_val;
        time_start = nan;
        time_end = nan;
        time_window = 0.65;
        
        num_measurements = length(reference_frame_data);
        num_sides = length(side_list);
        num_segments = length(segment_list);
        
        original_view_vector = fig_scale.view_vector(subject_idx,:);
        original_azimuth = original_view_vector(1);
        original_elevation = original_view_vector(2);
        dTheta = 0;
        if ~isnan(num_measurements)
            if num_measurements > 1
                time_start = time_vector_val(1);
                time_end = time_vector_val(end);
                dTheta = rad2deg((2*pi)/(num_measurements-1));
            end
        end
        
        folder_name_start = "Visualization_Videos\movement_subNo_";
        if movingView_flag
            folder_name_start = "Visualization_MovingPtView_Videos\movement_subNo_";
        end
        
        my_writer = VideoWriter(folder_name_start + string(subject_idx), 'MPEG-4');
        my_writer.FrameRate = 50;
        open(my_writer);
        
        figh = figure;
        figh.WindowState = 'maximized';
        for k=1:num_measurements
            if ~all(isnan(time_vector_val))
                current_time_val = time_vector_val(k);
                fn = string(fieldnames(knee_angle_struct.raw.right));
                num_subplot = length(fn);
                plt_width = 1;
                for j=1:num_subplot
                    current_fn = fn(j);
                    
                    right_qt.raw = knee_angle_struct.raw.right.(current_fn);
                    left_qt.raw = knee_angle_struct.raw.left.(current_fn);
                    
                    right_qt.processed = knee_angle_struct.processed.right.(current_fn);
                    left_qt.processed = knee_angle_struct.processed.left.(current_fn);
                    
                    right_qt.rigidBody = knee_angle_struct.rigidBody.right.(current_fn);
                    left_qt.rigidBody = knee_angle_struct.rigidBody.left.(current_fn);
                    
                    if(~isempty(right_qt.raw) && ~isempty(left_qt.raw))
                        max_ht_qt = max([max(right_qt.raw), max(left_qt.raw), max(right_qt.processed), max(left_qt.processed), max(right_qt.rigidBody), max(left_qt.rigidBody)], [], 'omitnan') + plt_width;
                        min_ht_qt = min([min(right_qt.raw), min(left_qt.raw), min(right_qt.processed), min(left_qt.processed), min(right_qt.rigidBody), min(left_qt.rigidBody)], [], 'omitnan') - plt_width;
                        title_subplot = current_fn;
                        subplot(2,3,j);
                        grid
                        cla;
                        
                        if rawData_plot_flag
                            plot(time_vector_val, right_qt.raw, 'r--')
                            hold on;
                            plot(time_vector_val, left_qt.raw, 'b--')
                        end
                        
                        if processedData_plot_flag
                            plot(time_vector_val, right_qt.processed, 'r')
                            hold on;
                            plot(time_vector_val, left_qt.processed, 'b')
                        end
                        
                        if rigidBody_plot_flag
                            plot(time_vector_val, right_qt.rigidBody, 'Color', '#8B0000')
                            hold on;
                            plot(time_vector_val, left_qt.rigidBody, 'Color', '#8A2BE2')
                        end
                        
                        if ~isnan(time_vector_val(k))
                            xline(time_vector_val(k),'LineWidth', 1.2);
                        end
                        
                        time_min_lt = max(time_start, current_time_val-time_window);
                        time_max_lt = min(time_end, current_time_val+time_window);
                        time_limit = [time_min_lt, time_max_lt];
                        plt_limit = [min_ht_qt, max_ht_qt];
                        if all(~isnan(time_limit))
                            xlim(time_limit);
                        end
                        if all(~isnan(plt_limit))
                            ylim(plt_limit);
                        end
                        
                        RTO = ground_contact_data.right_TO.time;
                        RHS = ground_contact_data.right_HS.time;
                        num_RHS = length(RHS);
                        num_RTO = length(RTO);
                        for m1=1:num_RHS
                            xline(RHS(m1),'LineWidth', 1.2, 'Color', 'magenta', 'Linestyle', '-');
                        end
                        for m2=1:num_RTO
                            xline(RTO(m2),'LineWidth', 1.2, 'Color', 'magenta', 'Linestyle', '--');
                        end
                        LTO = ground_contact_data.left_TO.time;
                        LHS = ground_contact_data.left_HS.time;
                        num_LHS = length(LHS);
                        num_LTO = length(LTO);
                        for n1=1:num_LHS
                            xline(LHS(n1),'LineWidth', 1.2, 'Color', 'green', 'Linestyle', '-');
                        end
                        for n2=1:num_LTO
                            xline(LTO(n2),'LineWidth', 1.2, 'Color', 'green', 'Linestyle', '--');
                        end
                        grid;
                        xlabel("time (s)");
                        ylabel("angle (degrees)");
                        title(title_subplot);
                        if j==1
                            description_axes = zeros(10, 1);
                            description_axes(1) = plot(NaN,NaN,'--r');
                            description_axes(2) = plot(NaN,NaN,'r');
                            description_axes(3) = plot(NaN,NaN,'Color','#8B0000');
                            description_axes(4) = plot(NaN,NaN,'--b');
                            description_axes(5) = plot(NaN,NaN,'b');
                            description_axes(6) = plot(NaN,NaN,'Color','#8A2BE2');
                            description_axes(7) = plot(NaN,NaN,'-m');
                            description_axes(8) = plot(NaN,NaN,'--m');
                            description_axes(9) = plot(NaN,NaN,'-g');
                            description_axes(10) = plot(NaN,NaN,'--g');
                            legend(description_axes, 'RightLegRaw: RedDash', 'RightLegProcessed: Red', 'RightLegRigidBody: Brown', 'LeftLegRaw: BlueDash', 'LeftLegProcessed: Blue', 'LeftLegRigidBody: Voilet', 'RightHS: Magenta', 'RightTO: MagentaDash', 'LeftHS: Green', 'LeftTO: GreenDash');
                        end
                    end
                end
            end
            
            subplot(2,3,[4,5,6]);
            cla;
            
            technical_frames = reference_frame_data(k).technical_frames;
            anatomical_frames = reference_frame_data(k).anatomical_frames;
            segment_visualize = reference_frame_data(k).segment_visualize;

            global_frame_plotter
            for j=1:num_segments
                current_segment = segment_list(j);
                for m=1:num_sides
                    current_side = side_list(m);
                    if segment_visualize_bool.technical
                        fill3(segment_visualize.technical.(current_segment).(current_side).X_val, -segment_visualize.technical.(current_segment).(current_side).Z_val, segment_visualize.technical.(current_segment).(current_side).Y_val,'m','FaceAlpha',.2)
                        technical_frame_plotter(technical_frames.(current_segment).(current_side).T_val, technical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                    if segment_visualize_bool.reconstructed
                        fill3(segment_visualize.reconstructed.(current_segment).(current_side).X_val, -segment_visualize.reconstructed.(current_segment).(current_side).Z_val, segment_visualize.reconstructed.(current_segment).(current_side).Y_val,'y','FaceAlpha',.2)
                        anatomical_frame_plotter(anatomical_frames.(current_segment).(current_side).T_val, anatomical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                end
            end
            if movingView_flag
                new_azimuth = original_azimuth + (k-1)*dTheta;
                fig_scale.view_vector(subject_idx,:) = [new_azimuth original_elevation];
            end
            if ~all(isnan(time_vector_val))
                pretty_plotter(fig_scale, subject_idx)
            end
            figh = gcf;
            current_frame = getframe(figh);
            writeVideo(my_writer, current_frame);
        end
        close(my_writer);
        
        if all_subject_flag
            close(figh);
        end

    end

end

function create_ground_contact_visualization(flag_input, reference_frame_data_all, meta_data, fig_scale, ground_contact_data)

    ground_contact_rawData = flag_input.ground_contact_rawData_val;
    ground_contact_visualization_flag = flag_input.ground_contact_visualization_flag_val;
    segment_visualize_bool = flag_input.segment_visualize_bool_val;
    all_subject_flag = flag_input.all_subject_flag_val;
    movingView_flag = flag_input.movingView_flag_val;
    
    if ground_contact_rawData
        reference_frame_data = reference_frame_data_all.raw;
    else
        reference_frame_data = reference_frame_data_all.raw;
    end

    if ground_contact_visualization_flag
        
        subject_idx = meta_data.subject_idx_val;
        side_list = meta_data.side_list_val;
        segment_list = meta_data.segment_list_val;
        
        time_vector_val = meta_data.time_vector_val;
        time_start = nan;
        time_end = nan;
        time_window = 0.65;
        
        num_measurements = length(reference_frame_data);
        num_sides = length(side_list);
        num_segments = length(segment_list);
        
        original_view_vector = fig_scale.view_vector(subject_idx,:);
        original_azimuth = original_view_vector(1);
        original_elevation = original_view_vector(2);
        dTheta = 0;
        if ~isnan(num_measurements)
            if num_measurements > 1
                time_start = time_vector_val(1);
                time_end = time_vector_val(end);
                dTheta = rad2deg((2*pi)/(num_measurements-1));
            end
        end
        
        folder_name_start = "GroundCt_Visualization_Videos\movement_subNo_";
        if movingView_flag
            folder_name_start = "GroundCt_Visualization_MovingPtView_Videos\movement_subNo_";
        end
        
        my_writer = VideoWriter(folder_name_start + string(subject_idx), 'MPEG-4');
        my_writer.FrameRate = 50;
        open(my_writer);
        
        figh = figure;
        figh.WindowState = 'maximized';
        for k=1:num_measurements
            
            if ~all(isnan(time_vector_val))
                current_time_val = time_vector_val(k);
                subplot(2,2,1);
                grid;
                cla;
                plot(time_vector_val, ground_contact_data.right_heel, 'r', 'DisplayName', 'heel');
                hold on;
                plot(time_vector_val, ground_contact_data.right_foot, 'b', 'DisplayName', 'foot');
                xline(time_vector_val(k),'LineWidth', 1.2, 'DisplayName', 'time');
                RTO = ground_contact_data.right_TO.time;
                RHS = ground_contact_data.right_HS.time;
                num_RHS = length(RHS);
                num_RTO = length(RTO);
                for m1=1:num_RHS
                    xline(RHS(m1),'LineWidth', 1.2, 'Color', 'black', 'Linestyle', '-');
                end
                for m2=1:num_RTO
                    xline(RTO(m2),'LineWidth', 1.2, 'Color', 'black', 'Linestyle', '--');
                end
                grid;
                xlabel('time (s)');
                ylabel('displacement (mm)');
                title('right leg');
                description_axes = zeros(4, 1);
                description_axes(1) = plot(NaN,NaN,'-r');
                description_axes(2) = plot(NaN,NaN,'-b');
                description_axes(3) = plot(NaN,NaN,'-black');
                description_axes(4) = plot(NaN,NaN,'--black');
                legend(description_axes, 'Heel: Red', 'Foot: Blue', 'HS: Black', 'TO: BlackDash');
                time_min_lt = max(time_start, current_time_val-time_window);
                time_max_lt = min(time_end, current_time_val+time_window);
                xlim([time_min_lt, time_max_lt]);
                
                subplot(2,2,2);
                grid;
                cla;
                plot(time_vector_val, ground_contact_data.left_heel, 'r', 'DisplayName', 'heel');
                hold on;
                plot(time_vector_val, ground_contact_data.left_foot, 'b', 'DisplayName', 'foot');
                xline(time_vector_val(k),'LineWidth', 1.2, 'DisplayName', 'time');
                LTO = ground_contact_data.left_TO.time;
                LHS = ground_contact_data.left_HS.time;
                num_LHS = length(LHS);
                num_LTO = length(LTO);
                for n1=1:num_LHS
                    xline(LHS(n1),'LineWidth', 1.2, 'Color', 'black', 'Linestyle', '-');
                end
                for n2=1:num_LTO
                    xline(LTO(n2),'LineWidth', 1.2, 'Color', 'black', 'Linestyle', '--');
                end
                grid;
                xlabel('time (s)');
                ylabel('displacement (mm)');
                title('left leg');
                time_min_lt = max(time_start, current_time_val-time_window);
                time_max_lt = min(time_end, current_time_val+time_window);
                xlim([time_min_lt, time_max_lt]);
            end
            
            subplot(2,2,[3,4]);
            cla;
            
            technical_frames = reference_frame_data(k).technical_frames;
            anatomical_frames = reference_frame_data(k).anatomical_frames;
            segment_visualize = reference_frame_data(k).segment_visualize;

            global_frame_plotter
            for j=1:num_segments
                current_segment = segment_list(j);
                for m=1:num_sides
                    current_side = side_list(m);
                    if segment_visualize_bool.technical
                        fill3(segment_visualize.technical.(current_segment).(current_side).X_val, -segment_visualize.technical.(current_segment).(current_side).Z_val, segment_visualize.technical.(current_segment).(current_side).Y_val,'m','FaceAlpha',.2)
                        technical_frame_plotter(technical_frames.(current_segment).(current_side).T_val, technical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                    if segment_visualize_bool.reconstructed
                        fill3(segment_visualize.reconstructed.(current_segment).(current_side).X_val, -segment_visualize.reconstructed.(current_segment).(current_side).Z_val, segment_visualize.reconstructed.(current_segment).(current_side).Y_val,'y','FaceAlpha',.2)
                        anatomical_frame_plotter(anatomical_frames.(current_segment).(current_side).T_val, anatomical_frames.(current_segment).(current_side).unit_vector_hat_val)
                    end
                end
            end
            if movingView_flag
                new_azimuth = original_azimuth + (k-1)*dTheta;
                fig_scale.view_vector(subject_idx,:) = [new_azimuth original_elevation];
            end
            if ~all(isnan(time_vector_val))
                pretty_plotter(fig_scale, subject_idx)
            end
            figh = gcf;
            current_frame = getframe(figh);
            writeVideo(my_writer, current_frame);
        end
        close(my_writer);
        
        if all_subject_flag
            close(figh);
        end

    end

end

function ground_contact_struct = extract_ground_contact_data(activity, meta_data, direction_vector)
    activity_idx = meta_data.activity_idx_val;
    subject_idx = meta_data.subject_idx_val;
    time_vector_val = meta_data.time_vector_val;
    default_trial_idx = 1;
    sacrum = meta_data.sacral_data_val;
   
    ground_contact_struct = struct;
    
    ground_contact_struct.right_heel = direction_vector.*(activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.right_heel_x - sacrum);
    ground_contact_struct.left_heel = direction_vector.*(activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.left_heel_x - sacrum);
    ground_contact_struct.right_foot = direction_vector.*(activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.right_metatarsal2_x - sacrum);
    ground_contact_struct.left_foot = direction_vector.*(activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.left_metatarsal2_x - sacrum);
    
    fn_list = fieldnames(ground_contact_struct);
    fn_num = length(fn_list);
    for j=1:fn_num
        fn_current = string(fn_list{j});
        fn_val = ground_contact_struct.(fn_current);
        fn_val_noNaN = fn_val(~isnan(fn_val));
        fn_val_avg = mean(fn_val_noNaN);
        if contains(fn_current,"heel")
            pks = nan;
            locs = [];
            if length(fn_val) > 3
                [pks,locs] = findpeaks(fn_val);
                if ~isempty(locs)
                    for n=1:length(locs)
                        if ~(direction_vector(locs(n)) == direction_vector(locs(n)+1))
                            locs(n) = nan;
                            pks(n) = nan;
                        end
                        start_idx = max(1, locs(n)-5);
                        end_idx = min(length(fn_val), locs(n)+5);
                        fn_val_now = fn_val(start_idx:end_idx);
                        time_val_now = time_vector_val(start_idx:end_idx);
                        sgf_now = sgolayfilt(fn_val_now,4,7);
                        [pks_filt, ~] = findpeaks(sgf_now);
                        if isempty(pks_filt)
                            locs(n) = nan;
                            pks(n) = nan;
                        end
                        if pks(n) < fn_val_avg
                            locs(n) = nan;
                            pks(n) = nan;
                        end
                    end
                end
            end
            pks = pks(~isnan(pks));
            locs = locs(~isnan(locs));
            fn_peak = pks;
            fn_time = time_vector_val(locs);
            if contains(fn_current,"right")
                ground_contact_struct.right_HS.time = fn_time;
                ground_contact_struct.right_HS.val = fn_peak;
            else
                ground_contact_struct.left_HS.time = fn_time;
                ground_contact_struct.left_HS.val = fn_peak;
            end
        end
        
        if contains(fn_current,"foot")
            valley = nan;
            locs = [];
            if length(fn_val) > 3
                [valley,locs] = findpeaks(-fn_val);
                if ~isempty(locs)
                    for n=1:length(locs)
                        if ~(direction_vector(locs(n)) == direction_vector(locs(n)+1))
                            locs(n) = nan;
                            valley(n) = nan;
                        end
                        start_idx = max(1, locs(n)-5);
                        end_idx = min(length(fn_val), locs(n)+5);
                        fn_val_now = fn_val(start_idx:end_idx);
                        time_val_now = time_vector_val(start_idx:end_idx);
                        sgf_now = sgolayfilt(fn_val_now,4,7);
                        [valley_filt, ~] = findpeaks(-sgf_now);
                        if isempty(valley_filt)
                            locs(n) = nan;
                            valley(n) = nan;
                        end
                        if -valley(n) > fn_val_avg
                            locs(n) = nan;
                            valley(n) = nan;
                        end
                    end
                end
            end
            valley = valley(~isnan(valley));
            locs = locs(~isnan(locs));
            fn_valley = -valley;
            fn_time = time_vector_val(locs);
            if contains(fn_current,"right")
                    ground_contact_struct.right_TO.time = fn_time;
                    ground_contact_struct.right_TO.val = fn_valley;
            else
                ground_contact_struct.left_TO.time = fn_time;
                ground_contact_struct.left_TO.val = fn_valley;
            end
        end
    end
    
end

function direction_vector = estimate_direction(meta_data)
    plot_check_flag = false;
    if plot_check_flag
        overwrite_dirVec_folder = false;
        if or(~isfolder("DirectionVector_Plots"), overwrite_dirVec_folder)
            mkdir DirectionVector_Plots
        end
        subject_idx = meta_data.subject_idx_val;
    end
    
    data = meta_data.sacral_data_val;
    time_vector = meta_data.time_vector_val;
    
    num_measurements = length(time_vector);
    direction_vector = nan(num_measurements, 1);
    if num_measurements <= 1
        return
    end
    data_NaN_idx = isnan(data);
    data_Not_NaN_idx = ~data_NaN_idx;
    numeric_data_Not_NaN_idx = double(data_Not_NaN_idx);
    jump_data = ischange(numeric_data_Not_NaN_idx);
    [~, locs] = findpeaks(double(jump_data));
    jump_ct = length(locs);
    
    start_idx_list = [];
    end_idx_list = [];
    
    if ~isnan(data(1))
    start_idx_list = [start_idx_list; 1];
    end
    
    for j=1:jump_ct
        current_idx = locs(j);
        jump_direction = numeric_data_Not_NaN_idx(current_idx);
        if jump_direction == 1
            start_idx_list = [start_idx_list; current_idx];
        else
            end_idx_list = [end_idx_list; current_idx-1];
        end
    end
    
    if ~isnan(data(end))
        end_idx_list = [end_idx_list; length(data)];
    end
    
    num_segments = length(start_idx_list);
    
    for m=1:num_segments
        f = find(isnan(data(start_idx_list(m):end_idx_list(m))));
        if ~isempty(f)
            new_end = start_idx_list(m)+f(1)-2;
            end_idx_list(m) = new_end;
        end
    end
    
    for k=1:num_segments
        segment_pt_peaks_idx = [];
        segment_pt_valleys_idx = [];
        current_start_idx = start_idx_list(k);
        current_end_idx = end_idx_list(k);
        current_diff_idx = current_end_idx - current_start_idx;
        current_list_idx = linspace(current_start_idx, current_end_idx, current_diff_idx+1);
        current_data = data(current_start_idx:current_end_idx);
        current_time = time_vector(current_start_idx:current_end_idx);
        data_locs = [];
        if length(current_time) > 3
            [~, data_locs] = findpeaks(current_data);
        end
        time_pks = current_time(data_locs);
        current_pks_idx = nan;
        if ~isempty(time_pks)
            for j=1:length(time_pks)
                if or(current_data(data_locs(j)) < current_data(1), current_data(data_locs(j)) < current_data(end))
                    time_pks(j) = nan;
                end
                if ~isnan(time_pks(j))
                    current_pks_idx = current_list_idx(data_locs(j));
                    N = 5;
                    fit_start_idx = max(current_start_idx, current_pks_idx - N);
                    fit_end_idx = min(current_end_idx, current_pks_idx + N);
                    time_fit = time_vector(fit_start_idx:fit_end_idx);
                    data_fit = data(fit_start_idx:fit_end_idx);
                    [~,S] = polyfit(time_fit, data_fit, 2);
                    if ~isnan(S.normr)
                        if abs(S.normr) > 2.0
                            time_pks(j) = nan;
                        end
                    end
                end
                if ~isnan(time_pks(j))
                    segment_pt_peaks_idx = [segment_pt_peaks_idx; current_pks_idx];
                end
            end
        end
        
        data_vall_locs = [];
        if length(current_data) > 3
            [~, data_vall_locs] = findpeaks(-current_data);
        end
        time_vall = current_time(data_vall_locs);
        current_vall_idx = nan;
        if ~isempty(time_vall)
            for j=1:length(time_vall)
                if or(current_data(data_vall_locs(j)) > current_data(1), current_data(data_vall_locs(j)) > current_data(end))
                    time_vall(j) = nan;
                end
                if ~isnan(time_vall(j))
                    current_vall_idx = current_list_idx(data_vall_locs(j));
                    N = 5;
                    fit_start_idx = max(current_start_idx, current_vall_idx - N);
                    fit_end_idx = min(current_end_idx, current_vall_idx + N);
                    time_fit = time_vector(fit_start_idx:fit_end_idx);
                    data_fit = data(fit_start_idx:fit_end_idx);
                    [~,S] = polyfit(time_fit, data_fit, 2);
                    if ~isnan(S.normr)
                        if abs(S.normr) > 2.0
                            time_vall(j) = nan;
                        end
                    end
                end
                if ~isnan(time_vall(j))
                    segment_pt_valleys_idx = [segment_pt_valleys_idx; current_vall_idx];
                end
            end
        end
        
        if isempty(segment_pt_peaks_idx) && isempty(segment_pt_valleys_idx)
            if ~isempty(current_data)
                if current_data(1) < current_data(end)
                    direction_vector(current_start_idx:current_end_idx) = 1;
                else
                    direction_vector(current_start_idx:current_end_idx) = -1;
                end
            end
        else
            positive_velocity = false;
            negative_velocity = false;
            positive_switch = 0;
            negative_switch = 0;
            for n=current_start_idx:current_end_idx
                
                if ismember(n, segment_pt_peaks_idx)
                    positive_velocity = false;
                    negative_velocity = true;
                    negative_switch = negative_switch + 1;
                end
                
                if ismember(n, segment_pt_valleys_idx)
                    positive_velocity = true;
                    negative_velocity = false;
                    positive_switch = positive_switch + 1;
                end
                
                if positive_velocity
                    direction_vector(n) = 1;
                end
                
                if negative_velocity
                    direction_vector(n) = -1;
                end
                
                if positive_switch == 1 && all(isnan(direction_vector(current_start_idx:n-1)))
                    direction_vector(current_start_idx:n-1) = -1;
                end
                
                if negative_switch == 1 && all(isnan(direction_vector(current_start_idx:n-1)))
                    direction_vector(current_start_idx:n-1) = 1;
                end
                
            end
         
        end
    end
    
    if plot_check_flag
        scaling_factor = max(abs(max(data)), abs(min(data)));
        scaled_direction_vector = scaling_factor * direction_vector;
        figure
        plot(time_vector, data, 'r', 'DisplayName', 'data')
        hold on
        plot(time_vector, scaled_direction_vector, 'b', 'DisplayName', 'direction')
        yline(0, 'black', 'DisplayName', 'x-axis')
        grid
        ylim([-1.2*scaling_factor, 1.2*scaling_factor])
        legend
        figh = gcf;
        figh.WindowState = 'maximized';
        saveas(figh,"DirectionVector_Plots\DirectionVectorPlot_subNo_" + string(subject_idx), "jpg");
        close(figh);
    end
    
end