function [IMU] = IMU_read(filepath)
filename = convertStringsToChars(filepath);

data_info1 = h5info(filename,'/Sensors');
data_info2 = h5info(filename,'/Processed');

num_sensors = length(data_info1.Groups);

button = h5read(filename,'/Annotations');
IMU.button.push_time = button.Time;
IMU.button.sensor = button.SensorID;
IMU.button.label = cellstr(button.Annotation');

sensor_names = cell(num_sensors,1); 
vectorlengths = zeros(num_sensors,1);

sensor_data_test = data_info1.Groups(1).Name;
sensor_data_size = length(h5read(filename,[sensor_data_test,'/Time']));
sensor_data_empty = sensor_data_size == 0;
if sensor_data_empty
    IMU = struct;
    return
end


for ii = 1:num_sensors
    sensor_data_loc1 = data_info1.Groups(ii).Name;
    sensor_data_loc2 = data_info2.Groups(ii).Name;
    
    sensor_info = h5info(filename,[sensor_data_loc1,'/Configuration']);
    sensor_label = sensor_info.Attributes(1).Value;
    sensor_label = strrep(sensor_label,' ','_');
    
    sensor_names(ii) = {sensor_label};
    IMU.(sensor_label).a = h5read(filename,[sensor_data_loc1,'/Accelerometer'])';
    IMU.(sensor_label).w = h5read(filename,[sensor_data_loc1,'/Gyroscope'])';
    IMU.(sensor_label).m = h5read(filename,[sensor_data_loc1,'/Magnetometer'])';
    IMU.(sensor_label).q = h5read(filename,[sensor_data_loc2,'/Orientation'])';
    IMU.(sensor_label).b = h5read(filename,[sensor_data_loc1,'/Barometer']);
    IMU.(sensor_label).temp = h5read(filename,[sensor_data_loc1,'/Temperature']);
    IMU.(sensor_label).time = h5read(filename,[sensor_data_loc1,'/Time']);
    
    vectorlengths(ii) = length(IMU.(sensor_label).time);
end

[maxlength, indmax] = max(vectorlengths);

for jj = 1:num_sensors
    if vectorlengths(jj) == maxlength
        
        if length(IMU.(char(sensor_names(jj))).a) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'a');
        end
        if length(IMU.(char(sensor_names(jj))).w) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'w');
        end
        if length(IMU.(char(sensor_names(jj))).m) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'m');
        end
        if length(IMU.(char(sensor_names(jj))).b) == length(IMU.(char(sensor_names(jj))).time)
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'b');
        end
        if length(IMU.(char(sensor_names(jj))).temp) == length(IMU.(char(sensor_names(jj))).time)      
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'temp');
        end
        
    elseif vectorlengths(jj) ~= maxlength
        index = 1:maxlength;
        
        if length(IMU.(char(sensor_names(jj))).a) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).a = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).a,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'a');
        end
        if length(IMU.(char(sensor_names(jj))).w) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).w = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).w,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'w');
        end
        if length(IMU.(char(sensor_names(jj))).m) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).m = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).m,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'m');
        end
        if length(IMU.(char(sensor_names(jj))).b) == length(IMU.(char(sensor_names(jj))).time)
            IMU.(char(sensor_names(jj))).b = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).b,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'b');
        end
        if length(IMU.(char(sensor_names(jj))).temp) == length(IMU.(char(sensor_names(jj))).time)      
            IMU.(char(sensor_names(jj))).temp = interp1(IMU.(char(sensor_names(jj))).time,IMU.(char(sensor_names(jj))).temp,IMU.(char(sensor_names(indmax))).time,'linear','extrap');
        else
            IMU.(char(sensor_names(jj))) = rmfield(IMU.(char(sensor_names(jj))),'temp');
        end
          
        [match, ~] = ismember(IMU.(char(sensor_names(indmax))).time,IMU.(char(sensor_names(jj))).time);
        index_createALL = index(~match);
        
        new_q = zeros(length(IMU.(char(sensor_names(indmax))).time),4);
        new_q(match,:) = IMU.(char(sensor_names(jj))).q;

        for qq = 1:length(index_createALL)
            index_create = index_createALL(qq);     
            t1 = IMU.(char(sensor_names(indmax))).time(index_create - 1);
            t2 = IMU.(char(sensor_names(indmax))).time(index_create + 1);
            tc = IMU.(char(sensor_names(indmax))).time(index_create);

            t2 = t2 - t1;
            tc = tc - t1;
            tc = tc/t2;

            q1 = IMU.(char(sensor_names(jj))).q(index_create - 1,:);
            q2 = IMU.(char(sensor_names(jj))).q(index_create,:);
            qc = slerp(q1, q2, tc);

            new_q(index_create,:) = qc;

        end
        
        IMU.(char(sensor_names(jj))).q = new_q;      
    end
end

IMU.button.push_index = zeros(length(IMU.button.push_time),1);

for qq = 1:length(IMU.button.push_time)
    findmatch = IMU.(char(sensor_names(1))).time == IMU.button.push_time(qq);
    if sum(findmatch) == 1
        IMU.button.push_index(qq) = find(IMU.(char(sensor_names(1))).time == IMU.button.push_time(qq));
    else
        [~,IMU.button.push_index(qq)] = min(abs(double(IMU.(char(sensor_names(1))).time) - double(IMU.button.push_time(qq))));
    end
end

IMU.time = double((IMU.(char(sensor_names(1))).time - IMU.(char(sensor_names(1))).time(1)))/(10^6);
IMU.button.push_time = double((IMU.button.push_time - IMU.(char(sensor_names(1))).time(1)))/(10^6);

for ii = 1:num_sensors
    IMU.(char(sensor_names(ii))) = rmfield(IMU.(char(sensor_names(ii))),'time');
end

end

