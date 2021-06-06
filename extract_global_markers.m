load_uvm_data_bool = false;
if load_uvm_data_bool
    load('UVM_data_regrouped.mat');
end

load('marker_names_saved.mat')
default_trial_idx = 1;
static_cal_file_idx = 22;
marker_names_common_global_list = marker_names_all_list;
for n=1:25
    subject_idx = n;
    for j=2:length(activity)
        if ~(all(activity(j).VICON(subject_idx, default_trial_idx).markers.present == 0))
            for k=1:length(marker_names_all_list)
                if ~activity(j).VICON(subject_idx, default_trial_idx).markers.present(k)
                    marker_names_common_global_list(k)="str_rem";
                end
            end
        end
    end
end
marker_names_common_global_list = marker_names_common_global_list(marker_names_common_global_list~="str_rem");
marker_names_common_global_non_imu_list = marker_names_common_global_list(~contains(marker_names_common_global_list, "cluster"));