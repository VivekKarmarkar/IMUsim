%% LOAD HEAVY DATA

clc
clear
close

load('UVM_data_regrouped.mat');

%% LOAD NON HEAVY DATA and SET PREFERENCES

animation_flag = true;
overwrite_animation_folder = false;
segment_visualize_bool.technical = true;
segment_visualize_bool.anatomical = true;
segment_visualize_bool.reconstructed = true;

all_subjects = false;
subject_idx_input = 25;

if or(~isfolder("Anatomical_Calibration_Videos"), overwrite_animation_folder)
    mkdir Anatomical_Calibration_Videos
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

load('body_segment_input_data.mat');
load('anatomical_markers_info_data.mat');
load('technical_frame_markers_info_data.mat');
load('anatomical_calibration_parameters_info_data.mat');

static_cal_idx = 22;

side_list = ["right", "left"];

segment_list = strings(1,2);
for k=1:length(body_segment_input)
    segment_list(1,k) = body_segment_input(k).name;
end

input_data_val = struct;
input_data_val.side_list_val = side_list;
input_data_val.segment_list_val = segment_list;
input_data_val.activity_idx_val = static_cal_idx;
input_data_val.technical_frame_markers_info_val = technical_frame_markers_info;
input_data_val.anatomical_markers_info_val = anatomical_markers_info;
input_data_val.anatomical_calibration_parameters_val = anatomical_calibration_parameters_info;

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
        reference_frame_data_current = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame);
        create_animation(animation_flag, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
        subject_reference_frame_data_ac(s).reference_frame_data = reference_frame_data_current;
        subject_reference_frame_data_ac(s).error_values = invalidate_anatomical_calibration(reference_frame_data_current, input_data_val);
    end
end

if one_subject
    s = subject_idx_val;
    input_data_val.subject_idx_val = s;
    meta_data_val.subject_idx_val = s;
    reference_frame_data_current = reference_frame_data_generate(input_data_val, activity, pose_estimate_anatomical_frame);
    create_animation(animation_flag, segment_visualize_bool, reference_frame_data_current, meta_data_val, fig_scale_val, all_subjects);
    subject_reference_frame_data_ac(s).reference_frame_data = reference_frame_data_current;
    subject_reference_frame_data_ac(s).error_values = invalidate_anatomical_calibration(reference_frame_data_current, input_data_val);
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
        end
    end
end

function [T,R,unit_vectors,success_bool] = pose_estimate_technical_frame(ptCloud)
    pts_collinear = true;
    success_bool = false;

    if ~any(any(isnan(ptCloud.Location)))
        [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
        if ~isempty(plane_inlierIndices)
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
    T = origin_local;
    i_hat = unit_vector_local.x;
    j_hat = unit_vector_local.y;
    k_hat = unit_vector_local.z;
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

function reference_frame_data = reference_frame_data_generate(input_data, activity, pose_estimate_anatomical_frame)
    reference_frame_data = struct;
    default_trial_idx = 1;
    subject_idx = input_data.subject_idx_val;
    activity_idx = input_data.activity_idx_val;
    technical_frame_markers_info = input_data.technical_frame_markers_info_val;
    anatomical_markers_info = input_data.anatomical_markers_info_val;
    anatomical_calibration_parameters = input_data.anatomical_calibration_parameters_val;
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    t = activity(activity_idx).VICON(subject_idx, default_trial_idx).movement_data.time;
    num_measurements = length(t);
    num_sides = length(side_list);
    num_segments = length(segment_list);
    for n=1:num_measurements
        technical_frames = struct;
        segment_visualize = struct;
        segment_visualize.technical = struct;
        segment_visualize.anatomical = struct;
        segment_visualize.reconstructed = struct;
        for j=1:num_segments
            current_segment = segment_list(j);
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
                Xe = nan(3,1);
                Ye = nan(3,1);
                Ze = nan(3,1);
                Ve = nan(3,3);
                num_trial_pt.anatomical = length(ra_now);
                for k=1:num_trial_pt.anatomical
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
            end
        end
        reference_frame_data(n).time = t(n);
        reference_frame_data(n).technical_frames = technical_frames;
        reference_frame_data(n).segment_visualize = segment_visualize;
    end
end

function create_animation(animation_flag, segment_visualize_bool, reference_frame_data, meta_data, fig_scale, all_subject_flag)

    if animation_flag
        
        subject_idx = meta_data.subject_idx_val;
        side_list = meta_data.side_list_val;
        segment_list = meta_data.segment_list_val;
        
        num_measurements = length(reference_frame_data);
        num_sides = length(side_list);
        num_segments = length(segment_list);
        
        figh = figure('units','normalized','outerposition',[0 0 1 1]);

        for k=1:num_measurements
            clf

            technical_frames = reference_frame_data(k).technical_frames;
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
                    end
                end
            end
            pretty_plotter(fig_scale, subject_idx)
            movie_vector(k) = getframe(figh);
        end

        my_writer = VideoWriter("Anatomical_Calibration_Videos\anatomical_calibration_subNo_" + string(subject_idx), 'MPEG-4');
        my_writer.FrameRate = 20;

        open(my_writer);
        writeVideo(my_writer, movie_vector);
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

function error_values = invalidate_anatomical_calibration(reference_frame_data, input_data)
    side_list = input_data.side_list_val;
    segment_list = input_data.segment_list_val;
    subject_idx = input_data.subject_idx_val;
    anatomical_calibration_parameters_info = input_data.anatomical_calibration_parameters_val;
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