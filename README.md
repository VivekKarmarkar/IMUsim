# IMUsim

Simulation of IMU data using optical motion capture data

The read_data_version_1 reads VICON trc data files into MATLAB struct via python function
The read_data_version_2 has support for OneDrive and local computer drive

Presentation_1 has information regarding data type conventions and standardization. It also has information about data size standardization.

Make sure to load marker_names_all_list while using this program. It can be regenerated through the read_data code if needed and be saved for use.

Set glob_path_folder to point to the correct drive / computer / folder_path

Make sure python trc_read module is in the same folder as MATLAB

Save the activity_struct and reuse it in a fresh code without worrying about the read_data code if desired
