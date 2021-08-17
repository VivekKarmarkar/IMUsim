Simulation of IMU data using optical motion capture data

The read_data_version_1 reads VICON trc data files into MATLAB struct via python function
The read_data_version_2 has support for OneDrive and local computer drive
The read_data_version_3 has support for reading IMU h5 data files
The read_data_version_4 handles erroneous files which exist but have no data

Presentation_1 has information regarding data type conventions and standardization. It also has information about data size standardization.
Presentation_5 has updated information regarding data type and size standardization and convention. It also contains information about the nomenclature map from the VICON marker names to the IMU sensor location names.

Make sure to load marker_names_all_list while using this program. It can be regenerated through the read_data code if needed and be saved for use.

The marker_names_all_list also has a text file included if required

Set glob_path_folder to point to the correct drive / computer / folder_path

Make sure python trc_read module is in the same folder as MATLAB

Make sure IMU_read, IMU_unpack_components functions is in the same folder as MATLAB

Special attention should be paid to the difference between referencing / calling VICON and IMU data files which has been included in Presentation_5

Save the activity_struct and reuse it in a fresh code without worrying about the read_data code if desired

Anatomical_Frames_implementation_version_1 details below: (Version_1 has basic pose estimator, noise present and uses GT marker)
The extract_global_markers extracts the common marker set and anatomical marker set that can be used to determine the technical marker set
The extract_technical_markers extracts the markers that can be used to construct the technical marker frame for each subject
The body_segment_input_info has information regarding the inputs that are required repeatedly for various body segments in other codes
The anatomical_calibration_extract extracts the anatomical calibration parameters from the static calibration files for all subjects
The anatomical_calibration_validate creates the animations and runs invalidation functions to check if the anatomical_calibration_extraction was successful
The movement_analysis_reference_frame creates animations that can be used to visualize the technical marker frames and reconstructed anatomical segments for a given movement trial

Note that anatomical frames are yet to be constructed. Version_1 only reconstructs the anatomical segment

Make sure that the following data has been loaded while using Version_1:
body_segment_input_data
marker_names_saved
UVM_data_regrouped
anatomical_markers_info_data
anatomical_calibration_parameters_info_data
technical_frame_markers_info_data

Anatomical_Frames_implementation_version_2 details below: (Differences from version_1)
PointCloud object transposed for collinearity test as this is conceptually correct although results are the same
Collinearity test has been made more robust by including a triangle inequality based check as well
Anatomical frame pose estimator for femoral and tibial segments have been included

Note that version_2 uses the Greater Trochanter (GT) marker for construction of the femoral reference frame

Anatomical_Frames_implementation_version_2 details below: (Differences from version_1,2)
Knee Angle plot functions have been included
Pose Component plot functions have been included

Anatomical_Frames_implementation_version_4 details below:
Memory issue fixed for animations
Rotating viewpoint functionality introduced
Appropriate scales have been applied to animation
Algorithm to estimate current movement direction
Gait event detection using direction estimation, sacral and leg data
Visualization tool incorporating animation and knee angle plots annoted with gait events and current time
Knee angle plots embedded in the visualization tool have a zoom window around current time

Anatomical_Frames_implementation_version_5_details below:
Raw data processed using cubic spline fit followed by 4th order low-pass Butterworth filter
Reprocessing of data to remove outliers
Anatomical Calibration performed on reporcessed data
Anatomical Calibration validated
Movement Analysis performed after validating Anatomical Calibration
All Force Plate Walk trials reprocessed
Anatomical Calibration and Movement Analysis results only generated for Force Plate Walk 1
Code refactored to accomodate simultaneous processing of Raw and Processed data
