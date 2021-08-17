%% USER INPUT

activity_start_idx = 5;
activity_end_idx = 5;
%activity_end_idx = 14;

%% LOAD HEAVY DATA

load('WalkingFP_RawData_reprocessed.mat');

%% LOAD NON-HEAVY DATA

load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');
load('anatomical_calibration_parameters_processed_info_data.mat');

%% GENERATE AND POPULATE DATA CONTAINERS

side_list = ["right", "left"];
segment_list = strings(1,2);
for k=1:length(body_segment_input)
    segment_list(1,k) = body_segment_input(k).name;
end
num_subjects = 25;

input_data_val = struct;
input_data_val.side_list_val = side_list;
input_data_val.segment_list_val = segment_list;
input_data_val.technical_frame_markers_info_val.processed = technical_frame_markers_info;
input_data_val.anatomical_markers_info_val = anatomical_markers_info;
input_data_val.anatomical_calibration_parameters_val.processed = anatomical_calibration_parameters_processed_info;

%% PROCESSING

for j=activity_start_idx:activity_end_idx
    activity_idx = j;
    activity_name = strrep(activity_reprocessed(activity_idx).name, " ", "_");
    input_data_val.activity_idx_val = activity_idx;
    for s=1:num_subjects
        input_data_val.subject_idx_val = s;
        rigidBody_anatomical_markers = reconstructed_data_anatomical_markers_generate(input_data_val, activity_reprocessed);
        [rigidBody_anatomical_markers_rgp, technical_frame_markers_rigidBody] = regrouped_data_generate(rigidBody_anatomical_markers, input_data_val);
        rigidBody_anatomical_markers_info.(activity_name)(s) = rigidBody_anatomical_markers_rgp;
        technical_frame_markers_rigidBody_info.(activity_name)(s) = technical_frame_markers_rigidBody;
        disp(s);
    end
end
save('rigidBody_anatomical_markers_info_data', 'rigidBody_anatomical_markers_info', '-v7.3');
save('technical_frame_markers_rigidBody_info_data', 'technical_frame_markers_rigidBody_info', '-v7.3');

function [regrouped_data, technical_frame_markers_rigidBody] = regrouped_data_generate(reconstructed_data, input_data)
    regrouped_data = struct;
    technical_frame_markers_rigidBody = struct;
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    num_measurements = length(reconstructed_data);
    for segment_idx=1:length(segment_list)
        current_segment = segment_list(segment_idx);
        for side_idx=1:length(side_list)
            current_side = side_list(side_idx);
            x_list = nan(num_measurements, 1);
            y_list = nan(num_measurements, 1);
            z_list = nan(num_measurements, 1);
            num_markers = length(reconstructed_data(1).anatomical_markers.(current_segment).(current_side));
            name_list = strings(3,1);
            for m=1:num_markers
                first_anatomical_struct = reconstructed_data(1).anatomical_markers.(current_segment).(current_side)(m);
                current_name = first_anatomical_struct.name;
                if ~contains(current_name, "fibular")
                    if current_segment == "tibial"
                        k = m-1;
                    else
                        k = m;
                    end
                    name_list(k,1) = current_name;
                    regrouped_data.(current_segment).(current_side)(k).name = current_name;
                    for n=1:num_measurements
                        current_anatomical_struct = reconstructed_data(n).anatomical_markers.(current_segment).(current_side)(m);
                        x_list(n) = current_anatomical_struct.x;
                        y_list(n) = current_anatomical_struct.y;
                        z_list(n) = current_anatomical_struct.z;
                    end
                    regrouped_data.(current_segment).(current_side)(k).x = x_list;
                    regrouped_data.(current_segment).(current_side)(k).x = x_list;
                    regrouped_data.(current_segment).(current_side)(k).y = y_list;
                    regrouped_data.(current_segment).(current_side)(k).z = z_list;
                end
            end
            technical_frame_markers_rigidBody.(current_segment).(current_side) = name_list;
        end
    end
end

function reconstructed_data = reconstructed_data_anatomical_markers_generate(input_data, activity)
    reconstructed_data = struct;
    default_trial_idx = 1;
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    anatomical_markers_info = input_data.anatomical_markers_info_val;
    anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.processed;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val.processed;
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    t = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.time;
    num_measurements = length(t);
    num_sides = length(side_list);
    num_segments = length(segment_list);
    for n=1:num_measurements
        technical_frames = struct;
        anatomical_markers = struct;
        for j=1:num_segments
            current_segment = segment_list(j);
            for m=1:num_sides
                current_side = side_list(m);
            
                [r1_now, r2_now, r3_now] = current_technical_marker_position(subject_idx, activity_idx, current_segment, current_side, n, technical_frame_markers_info, activity);
                current_ptCloud_val = pointCloud(horzcat(r1_now,r2_now,r3_now));
                [T, R, unit_vector_hat, success_bool] = pose_estimate_technical_frame(current_ptCloud_val);
                technical_frames.(current_segment).(current_side).ptCloud_val = current_ptCloud_val;
                
                if ~success_bool
                    input_data_nan = any(any(isnan([r1_now, r2_now, r3_now])));
                    if input_data_nan
                        disp("input data has at least one NaN value")
                    else
                        disp("NO NaN values------------------------")
                    end
                end
                technical_frames.(current_segment).(current_side).T_val = T;
                technical_frames.(current_segment).(current_side).R_val = R;
                technical_frames.(current_segment).(current_side).unit_vector_hat_val = unit_vector_hat;
                
                trial_pt = struct;
                current_anatomical_calibration_parameters = anatomical_calibration_parameters(subject_idx).(current_segment).(current_side);
                trial_pt.anatomical = struct;
                Xr = nan(3,1);
                Yr = nan(3,1);
                Zr = nan(3,1);
                current_reconstructed_position = struct;
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
                    current_reconstructed_position(k).name = current_anatomical_marker_name;
                    current_reconstructed_position(k).x = Xr(k);
                    current_reconstructed_position(k).y = Yr(k);
                    current_reconstructed_position(k).z = Zr(k);
                end
                anatomical_markers.(current_segment).(current_side) = current_reconstructed_position;
            end
        end
        reconstructed_data(n).anatomical_markers = anatomical_markers;
    end
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