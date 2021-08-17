%% LOAD DATA

load_uvm_data_processed_bool = true;
if load_uvm_data_processed_bool
    load('StaticCal_RawData_reprocessed.mat');
end
load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');

%% SWAP DATA

activity = activity_reprocessed;
clearvars activity_reprocessed;

%% PROCESS DATA

default_trial_idx = 1;
static_cal_idx = 22;
side_list = ["right", "left"];
num_side = length(side_list);
num_body_segment = length(body_segment_input);
num_subjects = 25;

anatomical_calibration_parameters_processed_info = struct;

for n=1:num_subjects
    disp(n)
    current_subject_idx = n;
    static_cal_empty = all(~activity(static_cal_idx).VICON(n, default_trial_idx).markers.present);
    for j=1:num_body_segment
        current_segment_name = body_segment_input(j).name;
        for k=1:num_side
            current_side = side_list(k);
            current_anatomical_marker_info = struct;
            anatomical_marker_name_list = anatomical_markers_info(current_subject_idx).(current_segment_name).(current_side).names;
            num_marker = length(anatomical_marker_name_list);
            for m=1:num_marker
                current_anatomical_marker_info(m).name = strings(1,1);
                current_anatomical_marker_info(m).x = nan;
                current_anatomical_marker_info(m).y = nan;
                current_anatomical_marker_info(m).z = nan;
                if ~static_cal_empty
                    current_technical_markers = technical_frame_markers_info(current_subject_idx).(current_segment_name).(current_side);
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
                            r1 = r1_list(t,:)';
                            r2 = r2_list(t,:)';
                            r3 = r3_list(t,:)';
                            ra = ra_list(t,:)';
                            ptCloud_val = pointCloud(horzcat(r1,r2,r3));
                            [T, R, unit_vector_hat] = pose_estimate_technical_frame(ptCloud_val);
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
            anatomical_calibration_parameters_processed_info(current_subject_idx).(current_segment_name).(current_side) = current_anatomical_marker_info;
            disp('----finished processing this round----')
        end
    end
end
save('anatomical_calibration_parameters_processed_info_data', 'anatomical_calibration_parameters_processed_info', '-v7.3');

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