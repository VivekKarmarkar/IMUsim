save_data_bool = true;

body_segment_input = struct;

body_segment_input(1).name = "femoral";
body_segment_input(1).marker_names_anatomical = ["lat_femoral_epicondyle";"med_femoral_epicondyle";"trochanter"];
body_segment_input(1).marker_names_global = ["lat_femoral_epicondyle";"anterior_thigh15";"lateral_thigh25"];

body_segment_input(2).name = "tibial";
body_segment_input(2).marker_names_anatomical = ["fibular_head";"tibial_tuberosity";"lateral_malleolus";"medial_malleolus"];
body_segment_input(2).marker_names_global = ["fibular_head";"lateral_malleolus";"anterior_shank30"];

if save_data_bool
    save('body_segment_input_data', 'body_segment_input', '-v7.3');
end