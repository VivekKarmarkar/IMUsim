Simulation of IMU data using optical motion capture data

The read_data_version_1 reads VICON trc data files into MATLAB struct via python function
The read_data_version_2 has support for OneDrive and local computer drive
The read_data_version_3 has support for reading IMU h5 data files

Presentation_1 has information regarding data type conventions and standardization. It also has information about data size standardization.
Presentation_5 has updated information regarding data type and size standardization and convention. It also contains information about the nomenclature map from the VICON marker names to the IMU sensor location names.

Make sure to load marker_names_all_list while using this program. It can be regenerated through the read_data code if needed and be saved for use.

The marker_names_all_list also has a text file included if required

Set glob_path_folder to point to the correct drive / computer / folder_path

Make sure python trc_read module is in the same folder as MATLAB

Make sure IMU_read, IMU_unpack_components functions is in the same folder as MATLAB

Special attention should be paid to the difference between referencing / calling VICON and IMU data files which has been included in Presentation_5

Save the activity_struct and reuse it in a fresh code without worrying about the read_data code if desired
