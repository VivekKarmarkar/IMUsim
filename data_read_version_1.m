% Start fresh by clearing the workspace and cleaning the command window
clear
clc
close all


% Setting up the commands for calling Python scripts functions saved in working directory
pyExec = 'C:\Anaconda3\';
pyRoot = fileparts(pyExec);
p = getenv('PATH');
p = strsplit(p, ';');
addToPath = {
   pyRoot
   fullfile(pyRoot, 'Library', 'mingw-w64', 'bin')
   fullfile(pyRoot, 'Library', 'usr', 'bin')
   fullfile(pyRoot, 'Library', 'bin')
   fullfile(pyRoot, 'Scripts')
   fullfile(pyRoot, 'bin')
};
p = [addToPath(:); p(:)];
p = unique(p, 'stable');
p = strjoin(p, ';');
setenv('PATH', p);


% Setting up folder paths and some constants
glob_path_folder = "C:\Users\vkarmarkar\Documents\Healthy Subjects v3\Healthy Subjects v3";
location_val = "Lab";
session_str_val = "\Session_1\";
VICON_file_ext_val = ".trc";
VICON_name_str_val = "\VICON\";
activity_str_del_list = ["VICON", "MLK", "2019", "accel", "gyro", "elec", ".h5", "ReadMe", "gaitReport"];


% Load exhaustive list of marker names
load('marker_names_saved.mat');
length_mn = length(marker_names_all_list);
false_mn = false(length_mn);
marker_names_false_all_list = false_mn(:,1);


% Extract meta data info from example file
example_file_path = 'C:\Users\vkarmarkar\Documents\Healthy Subjects v3\Healthy Subjects v3\S0025\Session_1\Lab\VICON\Air Squats.trc';
example_file_python_data_trc = py.trc_read.regroup_data(example_file_path);
example_file_meta_data = example_file_python_data_trc(1);
example_file_meta_data_dict = example_file_meta_data{1,1};
example_file_meta_data_struct = struct(example_file_meta_data_dict);
meta_data_fieldnames = fieldnames(example_file_meta_data_struct);
meta_data_default_struct = struct;
for i=1:length(meta_data_fieldnames)
    meta_data_default_struct.(string(meta_data_fieldnames{i,1})) = string;
end


% Initialize struct for tabular data
tabular_data_standardized_default_struct = struct;
tabular_data_standardized_default_struct.("frames") = nan;
tabular_data_standardized_default_struct.("time") = nan;
for j=1:length(marker_names_all_list)
    marker_name_val = marker_names_all_list(j);
    marker_name_x_val = marker_name_val + "_x";
    marker_name_y_val = marker_name_val + "_y";
    marker_name_z_val = marker_name_val + "_z";
    tabular_data_standardized_default_struct.(marker_name_x_val)=nan;
    tabular_data_standardized_default_struct.(marker_name_y_val)=nan;
    tabular_data_standardized_default_struct.(marker_name_z_val)=nan;
end


% Extracting subject and location information
glob_path_folder_val = glob_path_folder;
subject_idx_list = folder_content_glob(glob_path_folder);
location_all_list = folder_content_loc(glob_path_folder, subject_idx_list, "");


% Calculating maximum number of subjects
subject_max_str = subject_idx_list(end);
subject_max_str_list = strsplit(subject_max_str, "S");
subject_max_val = str2num(subject_max_str_list(end));


% Initializing list for sensor and activity information
sensor_name_all_list = string.empty;
sensor_placement_all_list = string.empty;
sensor_unknown_all_list = string.empty;
activity_description_all_list = string.empty;
activity_name_all_list = string.empty;
activity_file_ext_all_list = string.empty;


% Extracting sensor and activity information
for k=1:length(location_all_list)
    location_str = "\" + location_all_list(k);
    sensor_name_location_list = folder_content_loc(glob_path_folder, subject_idx_list, location_str);
    for j=1:length(sensor_name_location_list)
        sensor_name_bool = ~(sensor_name_location_list(j) == "APDM");
        sensor_name_str = location_str + "\" + sensor_name_location_list(j);
        sensor_next_level_sensor_name_list = folder_content_loc(glob_path_folder, subject_idx_list, sensor_name_str);
        sensor_last_level_bool = isempty(sensor_next_level_sensor_name_list);
        activity_sensor_name_list = string.empty;
        activity_str = "";
        activity_str_ini = glob_path_folder + "\" + subject_idx_list(k) + "\Session_1" + sensor_name_str;
        if sensor_last_level_bool
            activity_str = activity_str_ini;
            activity_sensor_name_struct = dir(activity_str);
            activity_sensor_name_cell = struct2cell(activity_sensor_name_struct);
            activity_sensor_name_list = string(activity_sensor_name_cell(2,:)) + "\" + string(activity_sensor_name_cell(1,:));
            activity_description_all_list = [activity_description_all_list, activity_sensor_name_list];
        else
            for m=1:length(sensor_next_level_sensor_name_list)
                activity_str = activity_str_ini + "\" + sensor_next_level_sensor_name_list(m);
                activity_sensor_name_struct = dir(activity_str);
            activity_sensor_name_cell = struct2cell(activity_sensor_name_struct);
            activity_sensor_name_list = string(activity_sensor_name_cell(2,:)) + "\" + string(activity_sensor_name_cell(1,:));
            activity_description_all_list = [activity_description_all_list, activity_sensor_name_list];
            end
        end
        if sensor_name_bool
            sensor_placement_all_list = [sensor_placement_all_list; sensor_next_level_sensor_name_list];
        else
            sensor_unknown_all_list = [sensor_unknown_all_list; sensor_next_level_sensor_name_list];
        end
    end
    sensor_name_all_list = [sensor_name_all_list; sensor_name_location_list];
end


% Processing activity information
for k=1:length(activity_description_all_list)
    [path, name, ext] = fileparts(activity_description_all_list(k));
    name = strrep(name,'_', ' ');
    for k2=1:length(activity_str_del_list)
        [match, no_match] = regexp(name,activity_str_del_list(k2),'match','split');
        if ~isempty(match)
            name = "";
        end
    end
    name = strtrim(name);
    activity_name_all_list = [activity_name_all_list; name];
    activity_file_ext_all_list = [activity_file_ext_all_list; ext];
end


% List formating: Removal of duplicates and blanks
sensor_name_duplicates = tabulate(sensor_name_all_list);
sensor_name_all_list = string(sensor_name_duplicates(:,1));
sensor_name_del_idx = cellfun(@isempty, sensor_name_all_list) == 0;
sensor_name_all_list = sensor_name_all_list(sensor_name_del_idx);

sensor_placement_duplicates = tabulate(sensor_placement_all_list);
sensor_placement_all_list = string(sensor_placement_duplicates(:,1));
sensor_placement_del_idx = cellfun(@isempty, sensor_placement_all_list) == 0;
sensor_placement_all_list = sensor_placement_all_list(sensor_placement_del_idx);

sensor_unknown_duplicates = tabulate(sensor_unknown_all_list);
sensor_unknown_all_list = string(sensor_unknown_duplicates(:,1));
sensor_unknown_del_idx = cellfun(@isempty, sensor_unknown_all_list) == 0;
sensor_unknown_all_list = sensor_unknown_all_list(sensor_unknown_del_idx);

activity_name_duplicates = tabulate(activity_name_all_list);
activity_name_all_list = string(activity_name_duplicates(:,1));
activity_name_del_idx = cellfun(@isempty, activity_name_all_list) == 0;
activity_name_all_list = activity_name_all_list(activity_name_del_idx);

activity_file_ext_duplicates = tabulate(activity_file_ext_all_list);
activity_file_ext_all_list = string(activity_file_ext_duplicates(:,1));
activity_file_ext_del_idx = cellfun(@isempty, activity_file_ext_all_list) == 0;
activity_file_ext_all_list = activity_file_ext_all_list(activity_file_ext_del_idx);


% Defining the bounds of some parameter values
n_max_activity = length(activity_name_all_list);
n_max_subject = subject_max_val;
n_max_trial = 2;


% Defining the structures
activity = struct;
activity(n_max_activity).name = string;
activity(n_max_activity).IMU = struct;
activity(n_max_activity).VICON = struct;


% Populating the structures
for j=1:n_max_activity
    activity_name_idx = j;
    activity_name_input = activity_name_all_list(j);
    disp(activity_name_input)
    
    activity(activity_name_idx).name = activity_name_input;
    activity(activity_name_idx).VICON = struct;
    activity(activity_name_idx).IMU = struct;
    
    activity(activity_name_idx).VICON(n_max_subject, n_max_trial).path = string;
    activity(activity_name_idx).VICON(n_max_subject, n_max_trial).file_ext = string;
    activity(activity_name_idx).VICON(n_max_subject, n_max_trial).meta_data = struct;
    activity(activity_name_idx).VICON(n_max_subject, n_max_trial).tabular_data = struct;
    activity(activity_name_idx).VICON(n_max_subject, n_max_trial).markers = struct;
    
    activity(activity_name_idx).IMU(n_max_subject, n_max_trial).path = string;
    activity(activity_name_idx).IMU(n_max_subject, n_max_trial).file_ext = string;
    activity(activity_name_idx).IMU(n_max_subject, n_max_trial).frame_rate = string;
    
    for n=1:n_max_subject
        subject_idx_input = n;
        
        for m=1:n_max_trial
            subject_idx_val = "\S" + num2str(subject_idx_input,'%04d');
            activity_name_val = activity_name_input;
            trial_number_val = m;
            
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).path = string;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).file_ext = string;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).meta_data = struct;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).tabular_data = struct;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers = struct;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.names = marker_names_all_list;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.present = marker_names_false_all_list;
            marker_present = activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.present;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.table = table(marker_names_all_list, marker_present);
            activity(activity_name_idx).IMU(subject_idx_input, trial_number_val).path = string;
            activity(activity_name_idx).IMU(subject_idx_input, trial_number_val).file_ext = string;
            activity(activity_name_idx).IMU(subject_idx_input, trial_number_val).frame_rate = string;
            
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).meta_data = meta_data_default_struct;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).tabular_data = tabular_data_standardized_default_struct;
            
            tabular_data_standardized_struct = tabular_data_standardized_default_struct;
            VICON_path_val = glob_path_folder_val + subject_idx_val + session_str_val + location_val + VICON_name_str_val + activity_name_val + VICON_file_ext_val;
            if isfile(VICON_path_val)
                activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).path = VICON_path_val;
                activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).file_ext = VICON_file_ext_val;
                python_data_trc = py.trc_read.regroup_data(VICON_path_val);
                meta_data = python_data_trc(1);
                meta_data_dict = meta_data{1,1};
                meta_data_struct = struct(meta_data_dict);
                tabular_data = python_data_trc(2);
                tabular_data_dict = tabular_data{1,1};
                tabular_data_struct = struct(tabular_data_dict);
                marker_names_data = python_data_trc(3);
                marker_names_dict = marker_names_data{1,1};
                
                marker_names_cell = cell(marker_names_dict);
                for k=1:length(marker_names_cell)
                    marker_name_val = string(marker_names_cell{1,k});
                    if ismember(marker_name_val, marker_names_all_list)
                        marker_name_idx = find(marker_names_all_list == marker_name_val);
                        activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.present(marker_name_idx) = true;
                    end
                    marker_names_cell{1,k} = marker_name_val;
                end
                marker_present = activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.present;
                activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).markers.table = table(marker_names_all_list, marker_present);
                marker_names_string = string(marker_names_cell);
                marker_names = transpose(marker_names_string);
                
                meta_data_fn = fieldnames(meta_data_struct);
                for k=1:length(meta_data_fn)
                    fn_cell = meta_data_fn(k);
                    fn = fn_cell{1,1};
                    fn_val_py_type = meta_data_struct.(fn);
                    fn_val_str = string(fn_val_py_type);
                    meta_data_struct.(fn) = fn_val_str;
                end
                
                tabular_data_fn = fieldnames(tabular_data_struct);
                for k=1:length(tabular_data_fn)
                    fn_cell = tabular_data_fn(k);
                    fn = fn_cell{1,1};
                    fn_val_py_type = tabular_data_struct.(fn);
                    fn_val_cell = cell(fn_val_py_type);
                    fn_val_array = cellfun(@double, fn_val_cell);
                    tabular_data_struct.(fn) = fn_val_array;
                end
                
                tabular_data_standardized_struct_size = length(tabular_data_struct.("frames"));
                tabular_data_standardized_struct.("frames") = tabular_data_struct.("frames");
                tabular_data_standardized_struct.("time") = tabular_data_struct.("time");
                
                for i=1:length(marker_names_all_list)
                    marker_name_val = marker_names_all_list(i);
                    marker_name_x_val = marker_name_val + "_x";
                    marker_name_y_val = marker_name_val + "_y";
                    marker_name_z_val = marker_name_val + "_z";

                    if ismember(marker_name_x_val, tabular_data_fn)
                        tabular_data_standardized_struct.(marker_name_x_val) = tabular_data_struct.(marker_name_x_val);
                    else
                        tabular_data_standardized_struct.(marker_name_x_val) = nan(1,tabular_data_standardized_struct_size);
                    end

                    if ismember(marker_name_y_val, tabular_data_fn)
                        tabular_data_standardized_struct.(marker_name_y_val) = tabular_data_struct.(marker_name_y_val);
                    else
                        tabular_data_standardized_struct.(marker_name_y_val) = nan(1,tabular_data_standardized_struct_size);
                    end

                    if ismember(marker_name_z_val, tabular_data_fn)
                        tabular_data_standardized_struct.(marker_name_z_val) = tabular_data_struct.(marker_name_z_val);
                    else
                        tabular_data_standardized_struct.(marker_name_z_val) = nan(1,tabular_data_standardized_struct_size);
                    end
                end
                
                activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).meta_data = meta_data_struct;
                activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).tabular_data = tabular_data_standardized_struct;
            end
            tabular_data_standardized_struct=structfun(@transpose,tabular_data_standardized_struct,'UniformOutput',false);
            tabular_data_table = struct2table(tabular_data_standardized_struct);
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).tabular_data = tabular_data_standardized_struct;
            activity(activity_name_idx).VICON(subject_idx_input, trial_number_val).tabular_data.table = tabular_data_table;
        end
    end
end


% Functions required for extracting folder content information
function f = folder_content_glob(path_folder)
    subfolder = dir(path_folder);
    num_subfolder = size(subfolder);
    num_subfolder = num_subfolder(1);
    name_subfolder_list = strings(num_subfolder, 1);
    path_subfolder_list = strings(num_subfolder, 1);
    name_subfolder_expression = 'S[00]';
    for j=1:num_subfolder
        name_subfolder_val = subfolder(j).name;
        name_subfolder_str = convertCharsToStrings(name_subfolder_val);
        name_subfolder_bool = ~isempty(regexp(name_subfolder_str,name_subfolder_expression, 'once'));
        path_subfolder_bool = subfolder(j).isdir;
        if path_subfolder_bool && name_subfolder_bool
            name_subfolder_list(j) = name_subfolder_str;
            path_subfolder_str = path_folder + "\" + name_subfolder_str;
            path_subfolder_list(j) = path_subfolder_str;
        end
    end
    idx = cellfun(@isempty, name_subfolder_list) == 0;
    name_subfolder_list = name_subfolder_list(idx);
    f = name_subfolder_list;
end

function g = folder_content_loc(path_folder, subject_idx_list, level_str)
    level_str_all_list = string.empty;
    for k=1:length(subject_idx_list)
        path_folder_new = path_folder + "\" + subject_idx_list(k) + "\Session_1" + level_str;
        subfolder = dir(path_folder_new);
        num_subfolder = size(subfolder);
        num_subfolder = num_subfolder(1);
        path_subfolder_list = strings(num_subfolder, 1);
        for j=1:num_subfolder
            name_subfolder_val = subfolder(j).name;
            name_subfolder_str = convertCharsToStrings(name_subfolder_val);
            name_subfolder_bool = true;
            path_subfolder_bool = subfolder(j).isdir;
            if path_subfolder_bool && name_subfolder_bool
                if ~ismember(name_subfolder_str, level_str_all_list)
                    level_str_all_list = [level_str_all_list; name_subfolder_str];
                end
                path_subfolder_str = path_folder_new + "\" + name_subfolder_str;
                path_subfolder_list(j) = path_subfolder_str;
            end
        end
    end
    level_str_all_list = level_str_all_list(~(level_str_all_list=='.'));
    level_str_all_list = level_str_all_list(~(level_str_all_list=='..'));
    g = level_str_all_list;
end

