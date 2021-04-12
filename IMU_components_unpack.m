function [sample_IMU_struct_renamed] = IMU_components_unpack(sample_IMU_struct)
sample_IMU_struct_renamed = sample_IMU_struct;
sample_IMU_struct_renamed_fieldnames = fieldnames(sample_IMU_struct_renamed);
sample_IMU_struct_renamed_substruct_fieldnames = fieldnames(sample_IMU_struct_renamed.sacrum);
for j=1:length(sample_IMU_struct_renamed_fieldnames)
    IMU_location_val = string(sample_IMU_struct_renamed_fieldnames(j));
    if and(IMU_location_val ~= "button", IMU_location_val ~= "time")
        for k=1:length(sample_IMU_struct_renamed_substruct_fieldnames)
           IMU_quantity_val = string(sample_IMU_struct_renamed_substruct_fieldnames(k));
           current_data = sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val);
           current_size = size(current_data);
           current_size_col = current_size(2);
           if current_size_col == 3
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_x") = current_data(:,1);
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_y") = current_data(:,2);
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_z") = current_data(:,3);
               sample_IMU_struct_renamed.(IMU_location_val) = rmfield(sample_IMU_struct_renamed.(IMU_location_val),IMU_quantity_val);
           end
           if current_size_col == 4
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_real") = current_data(:,1);
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_x") = current_data(:,2);
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_y") = current_data(:,3);
               sample_IMU_struct_renamed.(IMU_location_val).(IMU_quantity_val+"_z") = current_data(:,4);
               sample_IMU_struct_renamed.(IMU_location_val) = rmfield(sample_IMU_struct_renamed.(IMU_location_val),IMU_quantity_val);
           end
        end
    end
end
end

