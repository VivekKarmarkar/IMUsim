# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 16:55:22 2021

Python script to read and regroup VICON data files with 'trc' extension

Function in this module must be called from a MATLAB scipt
MATLAB script and this module should be stored in the same folder

@author: vkarmarkar
"""

import os
from trc import TRCData
mocap_data = TRCData()

path = r'C:\Users\vkarmarkar\Documents\Healthy Subjects v3\Healthy Subjects v3\S0012\Session_1\Lab\VICON\45 Static Lean.trc'


def regroup_data(trc_path_input: str):
    meta_data = dict()
    tabular_data = dict()
    trc_path_str = str(trc_path_input)
    trc_path_input = trc_path_input.encode('unicode_escape') # Converting to a raw string
    trc_path_exists = os.path.exists(trc_path_str)
    if trc_path_exists:
        mocap_data.load(trc_path_str)
        frames = mocap_data['Frame#']
        time = mocap_data['Time']
        marker = mocap_data['Markers']
        size_data = len(frames)
        size_marker = len(marker)
        marker_name_x = [str(marker[idx_x]) for idx_x in range(size_marker)]
        marker_name_y = [str(marker[idx_y]) for idx_y in range(size_marker)]
        marker_name_z = [str(marker[idx_z]) for idx_z in range(size_marker)]
        marker_name = marker_name_x + marker_name_y + marker_name_z
        tabular_data.fromkeys(marker_name)
        tabular_data['frames'] = frames
        tabular_data['time'] = time
        for j in range(size_marker):
            marker_name = str(marker[j])
            marker_name_x = marker_name + '_x'
            marker_name_y = marker_name + '_y'
            marker_name_z = marker_name + '_z'
            marker_data = mocap_data[marker_name]
            marker_data_x = []
            marker_data_y = []
            marker_data_z = []
            for k in range(size_data):
                marker_data_x.append(float(marker_data[k][0]))
                marker_data_y.append(float(marker_data[k][1]))
                marker_data_z.append(float(marker_data[k][2]))
            tabular_data[marker_name_x] = marker_data_x
            tabular_data[marker_name_y] = marker_data_y
            tabular_data[marker_name_z] = marker_data_z
        meta_data['Path_File_Type'] = str(mocap_data['PathFileType'])
        meta_data['Data_Format'] = str(mocap_data['DataFormat'])
        meta_data['Filename'] = str(mocap_data['FileName'])
        meta_data['Data_Rate'] = str(mocap_data['DataRate'])
        meta_data['Camera_Rate'] = str(mocap_data['CameraRate'])
        meta_data['Number_Of_Frames'] = str(mocap_data['NumFrames'])
        meta_data['Number_Of_Markers'] = str(mocap_data['NumMarkers'])
        meta_data['Units'] = str(mocap_data['Units'])
        meta_data['Orig_Data_Rate'] = str(mocap_data['OrigDataRate'])
        meta_data['Orig_Data_Start_Frame'] = str(mocap_data['OrigDataStartFrame'])
        meta_data['Orig_Num_Frames'] = str(mocap_data['OrigNumFrames'])
        marker_data = marker
    else:
        print("WARNING: trc_path does not exist")
    # convert to appropriate MATLAb data type in the MATLAB script
    return meta_data, tabular_data, marker_data
