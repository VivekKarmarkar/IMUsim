%% LOAD HEAVY DATA

clc
clear
close

load('UVM_data_regrouped.mat');
load('StaticCal_RawData_reprocessed.mat');

%% SET PREFERENCES

subject_idx_input = 10;
all_subjects = true;

animation_flag = true;
animation_movingView = true;
animation_rawData = false;
if animation_rawData
    animation_processedData = false;
else
    animation_processedData = true;
end

animation_settings_struct = struct;
animation_settings_struct.rawData_flag_val = animation_rawData;
animation_settings_struct.animation_flag_val = animation_flag;
animation_settings_struct.movingView_flag_val = animation_movingView;

segment_visualize_bool.technical = true;
segment_visualize_bool.anatomical = true;
segment_visualize_bool.reconstructed = true;

plot_flag.rawData = true;
plot_flag.processedData = true;

plot_flag.knee_angle = true;
plot_flag.pose_components = true;

overwrite_animation_folder = false;
overwrite_plot_folder.knee_angle = false;
overwrite_plot_folder.pose_components = false;

%% LOAD NON HEAVY DATA AND ACCOUNT FOR PREFERENCES

if or(~isfolder("Anatomical_Calibration_Videos"), overwrite_animation_folder)
    mkdir Anatomical_Calibration_Videos
end

if or(~isfolder("Anatomical_Calibration_MovingPtView_Videos"), overwrite_animation_folder)
    mkdir Anatomical_Calibration_MovingPtView_Videos
end

if or(~isfolder("Knee_Angle_Plots"), overwrite_plot_folder.knee_angle)
    mkdir Knee_Angle_Plots
end

if or(~isfolder("Pose_Components_Plots"), overwrite_plot_folder.pose_components)
    mkdir Pose_Components_Plots
end

if all_subjects
    one_subject = false;
else
    one_subject = true;
end

subject_idx_val = nan;
if one_subject
    subject_idx_val = subject_idx_input;
end

raw_data_bool = true;

load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');
load('anatomical_calibration_parameters_info_data.mat');
load('anatomical_calibration_parameters_processed_info_data.mat');

static_cal_idx = 22;

side_list = ["right", "left"];

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
input_data_val.activity_idx_val = static_cal_idx;
input_data_val.technical_frame_markers_info_val = technical_frame_markers_info;
input_data_val.anatomical_markers_info_val = anatomical_markers_info;
input_data_val.anatomical_calibration_parameters_val.raw = anatomical_calibration_parameters_info;
input_data_val.anatomical_calibration_parameters_val.processed = anatomical_calibration_parameters_processed_info;

meta_data_val = struct;
meta_data_val.side_list_val = side_list;
meta_data_val.segment_list_val = segment_list;

fig_scale_val = generate_fig_scale_info;

pose_estimate_anatomical_frame = array_pose_estimate_anatomical_frame;

%% COMPUTATIONS PART

subject_reference_frame_data_ac = struct;
num_subjects = 25;

if all_subjects
    for s=1:num_subjects
        input_data_val.subject_idx_val = s;
        meta_data_val.subject_idx_val = s;
        w = waitbar(0, 'Generating Reference Frame Data...');
        reference_frame_data_current.raw = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame, raw_data_bool);
        reference_frame_data_current.processed = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool);
        waitbar(0.33, w, 'Regrouping Reference Frame Data...');
        pose_data_current.raw = regroup_pose_data(reference_frame_data_current.raw, segment_list);
        pose_data_current.processed = regroup_pose_data(reference_frame_data_current.processed, segment_list);
        waitbar(0.67, w, 'Estimating Knee Angle...');
        knee_angle_data_current.raw = estimate_knee_angle_gs(reference_frame_data_current.raw);
        knee_angle_data_current.processed = estimate_knee_angle_gs(reference_frame_data_current.processed);
        time_vector = extract_time_vector(reference_frame_data_current.raw);
        waitbar(0.8, w, 'Generating Visualizations...');
        generate_plot_pose_components(plot_flag, time_vector, pose_data_current, s, segment_list);
        generate_plot_knee_angle(plot_flag, time_vector, knee_angle_data_current, s);
        create_animation(animation_settings_struct, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
        waitbar(0.9, w, 'Storing Data...');
        subject_reference_frame_data_ac(s).reference_frame_data.raw = reference_frame_data_current.raw;
        subject_reference_frame_data_ac(s).reference_frame_data.processed = reference_frame_data_current.processed;
        subject_reference_frame_data_ac(s).knee_angle_data.raw = knee_angle_data_current.raw;
        subject_reference_frame_data_ac(s).knee_angle_data.processed = knee_angle_data_current.processed;
        subject_reference_frame_data_ac(s).error_values.raw = invalidate_anatomical_calibration(reference_frame_data_current.raw, input_data_val, raw_data_bool);
        subject_reference_frame_data_ac(s).error_values.processed = invalidate_anatomical_calibration(reference_frame_data_current.processed, input_data_val, ~raw_data_bool);
        waitbar(1, w, 'Finishing...');
        close(w);
    end
end

if one_subject
    s = subject_idx_val;
    input_data_val.subject_idx_val = s;
    meta_data_val.subject_idx_val = s;
    w = waitbar(0, 'Generating Reference Frame Data...');
    reference_frame_data_current.raw = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame, raw_data_bool);
    reference_frame_data_current.processed = reference_frame_data_generate(input_data_val, activity_reprocessed, pose_estimate_anatomical_frame, ~raw_data_bool);
    waitbar(0.33, w, 'Regrouping Reference Frame Data...');
    pose_data_current.raw = regroup_pose_data(reference_frame_data_current.raw, segment_list);
    pose_data_current.processed = regroup_pose_data(reference_frame_data_current.processed, segment_list);
    waitbar(0.67, w, 'Estimating Knee Angle...');
    knee_angle_data_current.raw = estimate_knee_angle_gs(reference_frame_data_current.raw);
    knee_angle_data_current.processed = estimate_knee_angle_gs(reference_frame_data_current.processed);
    time_vector = extract_time_vector(reference_frame_data_current.raw);
    waitbar(0.8, w, 'Generating Visualizations...');
    generate_plot_pose_components(plot_flag, time_vector, pose_data_current, s, segment_list);
    generate_plot_knee_angle(plot_flag, time_vector, knee_angle_data_current, s);
    create_animation(animation_settings_struct, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
    waitbar(0.9, w, 'Storing Data...');
    subject_reference_frame_data_ac(s).reference_frame_data.raw = reference_frame_data_current.raw;
    subject_reference_frame_data_ac(s).reference_frame_data.processed = reference_frame_data_current.processed;
    subject_reference_frame_data_ac(s).knee_angle_data.raw = knee_angle_data_current.raw;
    subject_reference_frame_data_ac(s).knee_angle_data.processed = knee_angle_data_current.processed;
    subject_reference_frame_data_ac(s).error_values.raw = invalidate_anatomical_calibration(reference_frame_data_current.raw, input_data_val, raw_data_bool);
    subject_reference_frame_data_ac(s).error_values.processed = invalidate_anatomical_calibration(reference_frame_data_current.processed, input_data_val, ~raw_data_bool);
    waitbar(1, w, 'Finishing...');
    close(w);
end

save('reference_frame_data_ac_file', 'subject_reference_frame_data_ac', '-v7.3');

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

function ra_now_val = current_anatomical_marker_position(subject_idx, activity_idx, current_segment_val, current_side_val, t_now_idx, anatomical_markers_info, activity)
    default_trial_idx = 1;
    anatomical_markers_absent = anatomical_markers_info(subject_idx).(current_segment_val).(current_side_val).no_data_bool;
    anatomical_marker_name = anatomical_markers_info(subject_idx).(current_segment_val).(current_side_val).names;
    
    ra_now_val = struct;
    num_anatomical_markers = length(anatomical_marker_name);
    for k=1:num_anatomical_markers
        current_anatomical_marker_name = anatomical_marker_name(k);
        ra_now_val(k).value = nan(3,1);
        if ~anatomical_markers_absent
            xa_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name + "_x")(t_now_idx);
            ya_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name + "_y")(t_now_idx);
            za_now_k = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.(current_anatomical_marker_name + "_z")(t_now_idx);
            ra_now_k = [xa_now_k; ya_now_k; za_now_k];
            ra_now_val(k).value = ra_now_k;
            ra_now_val(k).name = current_anatomical_marker_name;
        end
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
    description_axes(5) = plot(NaN,NaN,'-c');
    description_axes(6) = plot(NaN,NaN,'-y');
    legend(description_axes, 'x: red','y: green','z: blue', 'technical segment: magenta', 'anatomical segment: cyan', 'reconstructed segment: yellow');
    xlim([fig_scale.x_min, fig_scale.x_max])
    ylim([fig_scale.y_min, fig_scale.y_max])
    zlim([fig_scale.z_min, fig_scale.z_max])
    view(fig_scale.view_vector(subject_idx,:));
end

function reference_frame_data = reference_frame_data_generate(input_data, activity, pose_estimate_anatomical_frame, raw_data_bool_val)
    reference_frame_data = struct;
    default_trial_idx = 1;
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val;
    anatomical_markers_info = input_data.anatomical_markers_info_val;
    if raw_data_bool_val
        anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.raw;
    else
        anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val.processed;
    end
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    anatomical_frame_function_idx_list = input_data.anatomical_frame_function_idx_list;
    t = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.time;
    num_measurements = length(t);
    num_sides = length(side_list);
    num_segments = length(segment_list);
    for n=1:num_measurements
        technical_frames = struct;
        anatomical_frames = struct;
        segment_visualize = struct;
        segment_visualize.technical = struct;
        segment_visualize.anatomical = struct;
        segment_visualize.reconstructed = struct;
        for j=1:num_segments
            current_segment = segment_list(j);
            current_anatomical_frame_function_idx = anatomical_frame_function_idx_list(j);
            for m=1:num_sides
                current_side = side_list(m);
                [r1_now, r2_now, r3_now] = current_technical_marker_position(subject_idx, activity_idx, current_segment, current_side, n, technical_frame_markers_info, activity);
                technical_frames.(current_segment).(current_side).ptCloud_val = pointCloud(horzcat(r1_now,r2_now,r3_now));
                [T, R, unit_vector_hat, success_bool] = pose_estimate_technical_frame(technical_frames.(current_segment).(current_side).ptCloud_val);
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
                ra_now = current_anatomical_marker_position(subject_idx, activity_idx, current_segment, current_side, n, anatomical_markers_info, activity);
                current_anatomical_calibration_parameters = anatomical_calibration_parameters(subject_idx).(current_segment).(current_side);
                trial_pt.anatomical = struct;
                Xa = nan(3,1);
                Ya = nan(3,1);
                Za = nan(3,1);
                Va = nan(3,3);
                Xr = nan(3,1);
                Yr = nan(3,1);
                Zr = nan(3,1);
                Vr = nan(3,3);
                current_reconstructed_position = struct;
                current_reconstructed_position.side = current_side;
                Xe = nan(3,1);
                Ye = nan(3,1);
                Ze = nan(3,1);
                Ve = nan(3,3);
                num_trial_pt.anatomical = length(anatomical_markers_info(subject_idx).(current_segment).(current_side).names);
                current_anatomical_marker_name_list = anatomical_markers_info(subject_idx).(current_segment).(current_side).names;
                for k=1:num_trial_pt.anatomical
                    current_anatomical_marker_name = current_anatomical_marker_name_list(k);
                    trial_pt.anatomical(k).name = current_anatomical_marker_name;
                    trial_pt.anatomical(k).global = ra_now(k).value;
                    Xa(k) = trial_pt.anatomical(k).global(1);
                    Ya(k) = trial_pt.anatomical(k).global(2);
                    Za(k) = trial_pt.anatomical(k).global(3);
                    Va(:,k) = [Xa(k); Ya(k); Za(k)];
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
                    Xe(k) = abs(Xa(k)-Xr(k));
                    Ye(k) = abs(Ya(k)-Yr(k));
                    Ze(k) = abs(Za(k)-Zr(k));
                    Ve(:,k) = [Xe(k); Ye(k); Ze(k)];
                    current_reconstructed_position.(current_anatomical_marker_name) = Vr(:,k);
                end
                segment_visualize.anatomical.(current_segment).(current_side).X_val = Xa;
                segment_visualize.anatomical.(current_segment).(current_side).Y_val = Ya;
                segment_visualize.anatomical.(current_segment).(current_side).Z_val = Za;
                segment_visualize.reconstructed.(current_segment).(current_side).X_val = Xr;
                segment_visualize.reconstructed.(current_segment).(current_side).Y_val = Yr;
                segment_visualize.reconstructed.(current_segment).(current_side).Z_val = Zr;
                segment_visualize.error.(current_segment).(current_side).X_val = mean(Xe);
                segment_visualize.error.(current_segment).(current_side).Y_val = mean(Ye);
                segment_visualize.error.(current_segment).(current_side).Z_val = mean(Ze);
                [T_af, R_af, unit_vector_hat_af] = pose_estimate_anatomical_frame{current_anatomical_frame_function_idx}(current_reconstructed_position);
                anatomical_frames.(current_segment).(current_side).T_val = T_af;
                anatomical_frames.(current_segment).(current_side).R_val = R_af;
                anatomical_frames.(current_segment).(current_side).unit_vector_hat_val = unit_vector_hat_af;
            end
        end
        reference_frame_data(n).time = t(n);
        reference_frame_data(n).technical_frames = technical_frames;
        reference_frame_data(n).anatomical_frames = anatomical_frames;
        reference_frame_data(n).segment_visualize = segment_visualize;
    end
end

function create_animation(animation_settings, segment_visualize_bool, reference_frame_data_all, meta_data, fig_scale, all_subject_flag)

    rawData_flag = animation_settings.rawData_flag_val;
    animation_flag = animation_settings.animation_flag_val;
    movingView_flag = animation_settings.movingView_flag_val;

    if rawData_flag
        reference_frame_data = reference_frame_data_all.raw;
    else
        reference_frame_data = reference_frame_data_all.processed;
    end

    if animation_flag
        
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
        
        folder_name_start = "Anatomical_Calibration_Videos\anatomical_calibration_subNo_";
        if movingView_flag
            folder_name_start = "Anatomical_Calibration_MovingPtView_Videos\anatomical_calibration_subNo_";
        end
        my_writer = VideoWriter(folder_name_start + string(subject_idx), 'MPEG-4');
        my_writer.FrameRate = 20;
        open(my_writer);
        
        figh = figure('units','normalized','outerposition',[0 0 1 1]);

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
                    if segment_visualize_bool.anatomical
                        fill3(segment_visualize.anatomical.(current_segment).(current_side).X_val, -segment_visualize.anatomical.(current_segment).(current_side).Z_val, segment_visualize.anatomical.(current_segment).(current_side).Y_val,'c','FaceAlpha',.2)
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
            pretty_plotter(fig_scale, subject_idx)
            current_frame = getframe(figh);
            writeVideo(my_writer, current_frame);
        end
        close(my_writer);
        
        if all_subject_flag
            close(figh);
        end

    end

end

function fig_scale_info = generate_fig_scale_info
    fig_scale_info = struct;
    
    fig_scale_info.x_max = 1000;
    fig_scale_info.y_max = 600;
    fig_scale_info.z_max = 1200;

    fig_scale_info.x_min = -400;
    fig_scale_info.y_min = -200;
    fig_scale_info.z_min = 0;
    
    view_vector_rep = [-109 7];
    view_vector = repmat(view_vector_rep,25,1);
    view_vector(2,:) = [-133 12];
    view_vector(6,:) = [-112 28];
    view_vector(13,:) = [-118 24];
    view_vector(18,:) = [-151 -2];
    
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

function error_values = invalidate_anatomical_calibration(reference_frame_data, input_data, raw_data_bool_val)
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    subject_idx = input_data.subject_idx_val;
    if raw_data_bool_val
        anatomical_calibration_parameters_info = input_data.anatomical_calibration_parameters_val.raw;
    else
        anatomical_calibration_parameters_info = input_data.anatomical_calibration_parameters_val.processed;
    end
    error_values = struct;
    error_tol = 2.0;
    for m=1:length(segment_list)
        current_segment = segment_list(m);
        error_values.(current_segment) = struct;
        for n=1:length(side_list)
            current_side = side_list(n);
            error_values.(current_segment).(current_side) = struct;
            x = nan(length(reference_frame_data),1);
            y = nan(length(reference_frame_data),1);
            z = nan(length(reference_frame_data),1);
            t = nan(length(reference_frame_data),1);
            for k=1:length(reference_frame_data)
                x(k) = reference_frame_data(k).segment_visualize.error.(current_segment).(current_side).X_val;
                y(k) = reference_frame_data(k).segment_visualize.error.(current_segment).(current_side).Y_val;
                z(k) = reference_frame_data(k).segment_visualize.error.(current_segment).(current_side).Z_val;
                t(k) = reference_frame_data(k).time;
            end
            error_values.(current_segment).(current_side).x = x;
            error_values.(current_segment).(current_side).y = y;
            error_values.(current_segment).(current_side).z = z;
            error_values.(current_segment).(current_side).t = t;
            error_representative = max([mean(x), mean(y), mean(z)]);
            current_anatomical_parameters_info = anatomical_calibration_parameters_info(subject_idx).(current_segment).(current_side);
            if error_representative > error_tol
                invalidate_anatomical_parameters_msg = ["subNo: ", string(subject_idx), current_segment, current_side, "error(mm): ", string(round(error_representative,3))];
                disp("--------FAILED: Anatomical Calibration--------")
                disp(invalidate_anatomical_parameters_msg)
                for m1=1:length(current_anatomical_parameters_info)
                    current_anatomical_parameters_info(m1).x = nan;
                    current_anatomical_parameters_info(m1).y = nan;
                    current_anatomical_parameters_info(m1).z = nan;
                end
            end
            anatomical_calibration_parameters_info(subject_idx).(current_segment).(current_side) = current_anatomical_parameters_info;
        end
    end
    anatomical_calibration_parameters_info_validated = anatomical_calibration_parameters_info;
    save('anatomical_calibration_parameters_info_validated_data', 'anatomical_calibration_parameters_info_validated', '-v7.3');
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
                if(~isempty(right_qt.raw) && ~isempty(left_qt.raw))
                    right_qt.raw_mean = round(mean(right_qt.raw, 'omitnan'),3);
                    left_qt.raw_mean = round(mean(left_qt.raw, 'omitnan'),3);
                    right_qt.raw_std = round(std(right_qt.raw, 'omitnan'),3);
                    left_qt.raw_std = round(std(left_qt.raw, 'omitnan'),3);
                    max_ht_qt = max([max(right_qt.raw), max(left_qt.raw)], [], 'omitnan') + plt_width;
                    min_ht_qt = min([min(right_qt.raw), min(left_qt.raw)], [], 'omitnan') - plt_width;
                    title_subplot = current_fn + " [right raw (avg,sd) = (" + string(right_qt.raw_mean) + "," + string(right_qt.raw_std) + "), left raw (avg,sd) = (" + string(left_qt.raw_mean) + "," + string(left_qt.raw_std) + ") ]";
                    subplot(num_subplot,1,j)
                    if plot_flag.rawData
                        plot(time_vector_val, right_qt.raw, 'r-.', 'DisplayName', 'right raw')
                        hold on
                        if ~isnan(right_qt.raw_mean)
                            yline(right_qt.raw_mean, 'r--', 'DisplayName', 'right raw avg')
                        end
                        plot(time_vector_val, left_qt.raw, 'b-.', 'DisplayName', 'left raw')
                        if ~isnan(left_qt.raw_mean)
                            yline(left_qt.raw_mean, 'b--', 'DisplayName', 'left raw avg')
                        end
                        plt_limit = [min_ht_qt, max_ht_qt];
                        if all(~isnan(plt_limit))
                            ylim(plt_limit)
                        end
                        title(title_subplot)
                    end
                    if plot_flag.processedData
                        plot(time_vector_val, right_qt.processed, 'r', 'DisplayName', 'right processed')
                        hold on
                        plot(time_vector_val, left_qt.processed, 'b', 'DisplayName', 'left processed')
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
            saveas(fig,"Knee_Angle_Plots\knee_angle_subNo_" + string(subject_idx), "jpg");
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
                        x_qt.raw = current_vector.raw.x;
                        y_qt.raw = current_vector.raw.y;
                        z_qt.raw = current_vector.raw.z;
                        x_qt.processed = current_vector.processed.x;
                        y_qt.processed = current_vector.processed.y;
                        z_qt.processed = current_vector.processed.z;
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
                                if ~isnan(x_qt.raw_mean)
                                    yline(x_qt.raw_mean, 'r--', 'DisplayName', 'x raw avg')   
                                end
                                if ~isnan(y_qt.raw_mean)
                                    yline(y_qt.raw_mean, 'g--', 'DisplayName', 'y raw avg')
                                end
                                if ~isnan(z_qt.raw_mean)
                                    yline(z_qt.raw_mean, 'b--', 'DisplayName', 'z raw avg')
                                end
                                plot(time_vector_val, y_qt.raw, 'g-.', 'DisplayName', 'y raw')
                                plot(time_vector_val, z_qt.raw, 'b-.', 'DisplayName', 'z raw')
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
                saveas(fig,"Pose_Components_Plots\pose_components_" + current_segment + "_subNo_" + string(subject_idx), "jpg");
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
                current_vector.raw_x = nan(num_measurements,1);
                current_vector.raw_y = nan(num_measurements,1);
                current_vector.raw_z = nan(num_measurements,1);
                for n=1:num_measurements
                    if current_fn == "T"
                        vector_from_struct = reference_frame_data_struct(n).anatomical_frames.(current_segment).(current_side).(current_fn_pose);
                    else
                        fn_pose_split = strsplit(current_fn_pose, ".");
                        vector_from_struct = reference_frame_data_struct(n).anatomical_frames.(current_segment).(current_side).(fn_pose_split(1)).(fn_pose_split(2));
                    end
                    current_vector.raw_x(n) = vector_from_struct(1);
                    current_vector.raw_y(n) = vector_from_struct(2);
                    current_vector.raw_z(n) = vector_from_struct(3);
                end
                pose_struct.(current_segment).(current_side).(current_fn).x = current_vector.raw_x;
                pose_struct.(current_segment).(current_side).(current_fn).y = current_vector.raw_y;
                pose_struct.(current_segment).(current_side).(current_fn).z = current_vector.raw_z;
            end
        end
    end
end