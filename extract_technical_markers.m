save_data_bool = true;

load_uvm_data_bool = false;
if load_uvm_data_bool
    load('UVM_data_regrouped.mat');
end
load('marker_names_saved.mat')
load('body_segment_input_data.mat')

anatomical_markers_info = struct;
technical_frame_markers_info = struct;

[anatomical_markers_info, technical_frame_markers_info] = marker_info(activity, marker_names_all_list, body_segment_input(1), anatomical_markers_info, technical_frame_markers_info);
[anatomical_markers_info, technical_frame_markers_info] = marker_info(activity, marker_names_all_list, body_segment_input(2), anatomical_markers_info, technical_frame_markers_info);

if save_data_bool
    save('anatomical_markers_info_data', 'anatomical_markers_info', '-v7.3');
    save('technical_frame_markers_info_data', 'technical_frame_markers_info', '-v7.3');
end

function [anatomical_markers_info, technical_frame_markers_info] = marker_info(activity, marker_names_all_list, body_segment_input, anatomical_markers_info, technical_frame_markers_info)

    body_segment_name = body_segment_input.name;
    body_segment_marker_names_anatomical = body_segment_input.marker_names_anatomical;
    body_segment_marker_names_global = body_segment_input.marker_names_global;

    global_technical_frame_markers = struct;

    side_list = ["right", "left"];
    static_cal_file_idx = 22;

    body_segment_empty_global = strings(3,1);
    body_segment_marker_names_global_left = body_segment_empty_global;
    body_segment_marker_names_global_right = body_segment_empty_global;
    for k=1:length(body_segment_empty_global)
        body_segment_marker_names_global_left(k) = "left_" + body_segment_marker_names_global(k);
        body_segment_marker_names_global_right(k) = "right_" + body_segment_marker_names_global(k);
    end

    for m=1:25
        body_segment = body_segment_name;
        current_subject_am_info = subject_am_info(m, activity, marker_names_all_list, body_segment, body_segment_marker_names_anatomical);
        if ~(all(activity(static_cal_file_idx).VICON(m, 1).markers.present == 0))
            global_technical_frame_markers(m).(body_segment).left = body_segment_marker_names_global_left;
            global_technical_frame_markers(m).(body_segment).right = body_segment_marker_names_global_right;
        else
            global_technical_frame_markers(m).(body_segment).left = body_segment_empty_global;
            global_technical_frame_markers(m).(body_segment).right = body_segment_empty_global;
        end

        for k2=1:2
            side = side_list(k2);
            tf_min_bool_val = current_subject_am_info.(body_segment).(side).technical_frame_minimum_bool;
            af_bool_val = current_subject_am_info.(body_segment).(side).anatomical_frame_bool;
            anatomical_markers_info(m).(body_segment).(side) = current_subject_am_info.(body_segment).(side);
            if ~af_bool_val
                if tf_min_bool_val
                    technical_frame_markers_info(m).(body_segment).(side) = current_subject_am_info.(body_segment).(side).technical_markers;
                else
                    technical_frame_markers_info(m).(body_segment).(side) = global_technical_frame_markers(m).(body_segment).(side);
                end
            else
                technical_frame_markers_info(m).(body_segment).(side) = strings(3,1);
            end
        end

    end

end

function anatomical_markers = subject_am_info(subject_idx, activity, marker_names_all_list, body_segment, body_segment_markers)
    default_trial_idx = 1;
    marker_names_common_list = marker_names_all_list;
    sum_null_activities = 0;
    for j=2:length(activity)
        if ~(all(activity(j).VICON(subject_idx, default_trial_idx).markers.present == 0))
            for k=1:length(marker_names_all_list)
                if ~activity(j).VICON(subject_idx, default_trial_idx).markers.present(k)
                    marker_names_common_list(k)="str_rem";
                end
            end
        else
            sum_null_activities = sum_null_activities + 1;
        end
    end
    marker_names_common_list = marker_names_common_list(marker_names_common_list~="str_rem");
    marker_names_common_non_imu_list = marker_names_common_list(~contains(marker_names_common_list, "cluster"));
    no_data_bool_val = false;
    
    if sum_null_activities == length(activity)-1
        marker_names_common_non_imu_list = strings(1,1);
        no_data_bool_val = true;
    end
    
    body_segment_empty = strings(size(body_segment_markers));
    body_segment_left_markers = body_segment_empty;
    body_segment_right_markers = body_segment_empty;
    for k=1:length(body_segment_empty)
        body_segment_left_markers(k) = "left_" + body_segment_markers(k);
        body_segment_right_markers(k) = "right_" + body_segment_markers(k);
    end

    anatomical_markers = struct;
    anatomical_markers.no_data_bool = no_data_bool_val;
    anatomical_markers.(body_segment).left = anatomical_markers_segment(body_segment_left_markers, marker_names_common_non_imu_list);
    anatomical_markers.(body_segment).right = anatomical_markers_segment(body_segment_right_markers, marker_names_common_non_imu_list);

end

function segment_markers = anatomical_markers_segment(marker_names, marker_common)
    length_mn = length(marker_names);
    false_mn = false(length_mn);
    marker_present = false_mn(:,1);
    for j=1:length(marker_common)
        marker_name_val = marker_common(j);
        if ismember(marker_name_val, marker_names)
            marker_name_idx = marker_names == marker_name_val;
            marker_present(marker_name_idx) = true;
        end
    end
    no_data_bool_val = false;
    if all(marker_common == "")
        no_data_bool_val = true;
    end
    segment_markers = struct;
    segment_markers.names = marker_names;
    segment_markers.present = marker_present;
    segment_markers.table_visualization = table(marker_names, marker_present);
    segment_markers.anatomical_frame_bool = all(marker_present);
    segment_markers.technical_frame_bool = sum(marker_present) >= 3;
    segment_markers.technical_frame_minimum_bool = sum(marker_present) == 3;
    segment_markers.technical_markers = marker_names(marker_present);
    segment_markers.no_data_bool = no_data_bool_val;
end