%% LOAD HEAVY DATA
load('StaticCal_RawData_reprocessed.mat');

%% SWAP DATA

activity = activity_reprocessed;
clearvars activity_reprocessed;

%% PROCESS DATA

overwrite_ProcessDataPlots_folder = false;
if or(~isfolder("ProcessedData_Plots"), overwrite_ProcessDataPlots_folder)
    mkdir ProcessedData_Plots
end

overwrite_OutlierPlots_folder = false;
if or(~isfolder("OutlierPlots"), overwrite_OutlierPlots_folder)
    mkdir OutlierPlots
end

load('marker_names_saved.mat')

gap_plot_flag = false;
processed_data_plot_flag = false;
OutlierPlot_flag = false;
activity_type_walking = true;

activity_idx_input = 5;

activity_start_idx = 2;
activity_end_idx = length(activity);
if activity_type_walking
    activity_start_idx = 5;
    activity_end_idx = 14;
end
static_cal_idx = 22;

one_activity = false;
if one_activity
    activity_start_idx = activity_idx_input;
    activity_end_idx = activity_idx_input;
end

filt_order = 4;
fc = 6;

coord_list = ["_x", "_y", "_z"];
error_threshold = [200.0, 200.0, 200.0];
marker_threshold = [4.0, 2.5, 2.5];

processing_inputs_struct = struct;
processing_inputs_struct.activity_start_idx_val = activity_start_idx;
processing_inputs_struct.activity_end_idx_val = activity_end_idx;
processing_inputs_struct.gap_plot_flag_val = gap_plot_flag;
processing_inputs_struct.processed_data_plot_flag_val = processed_data_plot_flag;
processing_inputs_struct.filt_order_val = filt_order;
processing_inputs_struct.fc_val = fc;
processing_inputs_struct.coord_list_val = coord_list;
processing_inputs_struct.error_threshold_val = error_threshold;
processing_inputs_struct.marker_threshold_val = marker_threshold;
processing_inputs_struct.marker_names_all_list_val = marker_names_all_list;
processing_inputs_struct.OutlierPlot_flag_val = OutlierPlot_flag;
if activity_type_walking
    activity_processed = save_processed_data(activity, processing_inputs_struct);
    processing_inputs_struct.remove_idx_val = manually_remove_outliers;
    [outliers_info, activity_reprocessed] = generate_processing_outliers(activity, activity_processed, processing_inputs_struct);
end

%% SAVE PROCESSED DATA

save('Walking_RawData_processed', 'activity_processed', '-v7.3');
save('Walking_RawData_reprocessed', 'activity_reprocessed', '-v7.3');

%% SANITY CHECK FOR PROCESSING

num_activities = length(activity);
for j=2:num_activities
    sum_residual = 0;
    sum_residual_old = 0;
    for s=1:25
        sum_residual_old = sum_residual_old + double(isfield(activity(j).VICON(s,1).markers, 'residual'));
        sum_residual = sum_residual + double(isfield(activity_reprocessed(j).VICON(s,1).markers, 'residual'));
    end
    output_disp = [activity(j).name sum_residual sum_residual_old];
    disp(output_disp);
end

%%

function [outliers_info, activity_reprocessed] = generate_processing_outliers(activity, activity_processed, processing_inputs)
    activity_reprocessed = activity_processed;
    marker_names_all_list = processing_inputs.marker_names_all_list_val;
    activity_start_idx = processing_inputs.activity_start_idx_val;
    activity_end_idx = processing_inputs.activity_end_idx_val;
    coord_list = processing_inputs.coord_list_val;
    marker_threshold = processing_inputs.marker_threshold_val;
    OutlierPlot_flag = processing_inputs.OutlierPlot_flag_val;
    filt_order = processing_inputs.filt_order_val;
    fc = processing_inputs.fc_val;
    error_threshold = processing_inputs.error_threshold_val;
    gap_plot_flag = processing_inputs.gap_plot_flag_val;
    remove_idx = processing_inputs.remove_idx_val;
    outlier_idx = 0;
    outliers_info = struct;
    max_gap = 10;
    for j=activity_start_idx:activity_end_idx
        for s=1:25
            for m=1:70
                Marker_name = marker_names_all_list(m);
                for c=1:3
                    Marker_reference = Marker_name + coord_list(c);
                    Time = activity(j).VICON(s, 1).movement_data.time;
                    RawData = activity(j).VICON(s, 1).movement_data.(Marker_reference);
                    ProcessedData = activity_processed(j).VICON(s, 1).movement_data.(Marker_reference);
                    if ~all(isnan(RawData))
                        Residual = activity_processed(j).VICON(s, 1).markers.residual(m, c);
                        if Residual > marker_threshold(c)
                            MissingData = true;
                            if all(~isnan(RawData))
                                MissingData = false;
                            end
                            outlier_idx = outlier_idx + 1;
                            outliers_info(outlier_idx).idx = outlier_idx;
                            outliers_info(outlier_idx).activity_idx = j;
                            outliers_info(outlier_idx).subject_idx = s;
                            outliers_info(outlier_idx).marker_idx = m;
                            outliers_info(outlier_idx).coordinate_idx = c;
                            outliers_info(outlier_idx).marker_name = Marker_name;
                            outliers_info(outlier_idx).coordinate = coord_list(c);
                            outliers_info(outlier_idx).RawData = RawData;
                            outliers_info(outlier_idx).ProcessedData = ProcessedData;
                            outliers_info(outlier_idx).MissingData = MissingData;
                            outliers_info(outlier_idx).Time = Time;
                            current_marker_name = Marker_name;
                            
                            ReprocessedData = nan(length(Time), 1);
                            Residual_current = nan;
                            
                            ModfiedData = RawData;
                            remove_idx_list = remove_idx(outlier_idx).listVal;
                            num_outliers = length(remove_idx_list);
                            for n=1:num_outliers
                                current_remove_idx = remove_idx_list(n);
                                ModfiedData(current_remove_idx) = nan;
                            end
                            RawData = ModfiedData;
                            
                            RawDataNAN = find(isnan(RawData));
                            if ~isempty(RawDataNAN)
                                End_idx = RawDataNAN(1);
                                Start_idx = RawDataNAN(end);
                                if length(RawData) - End_idx + 1 <= 1.5*max_gap
                                    RawData_current = RawData(1:End_idx);
                                    Time_current = Time(1:End_idx);
                                    [ProcessedData_current, Residual_current] = generate_processed_data(RawData_current, Time_current, filt_order, fc, error_threshold);
                                    ReprocessedData(1:End_idx) = ProcessedData_current;
                                end
                                if Start_idx <= 1.5*max_gap
                                    RawData_current = RawData(Start_idx:end);
                                    Time_current = Time(Start_idx:end);
                                    [ProcessedData_current, Residual_current] = generate_processed_data(RawData_current, Time_current, filt_order, fc, error_threshold);
                                    ReprocessedData(Start_idx:end) = ProcessedData_current;
                                end
                            end
                            
                            if all(isnan(ReprocessedData))
                                current_ErrorThreshold = error_threshold(c);
                                current_data = RawData;
                                current_data_copy = current_data;
                                [segmentStartIdx, segmentStopIdx] = gap_detector(s, current_data, gap_plot_flag);
                                [manual_segmentStartIdx, manual_segmentStopIdx, manual_segmentBoolVal] = manually_create_segment(j, s, m, c, processing_inputs);
                                if manual_segmentBoolVal
                                    segmentStartIdx = manual_segmentStartIdx;
                                    segmentStopIdx = manual_segmentStopIdx;
                                end
                                num_segments = length(segmentStartIdx);
                                ResidualSegment = nan(num_segments,1);
                                for ns=1:num_segments
                                    current_segment_start_idx = segmentStartIdx(ns);
                                    current_segment_stop_idx = segmentStopIdx(ns);
                                    current_RawData = current_data(current_segment_start_idx: current_segment_stop_idx);
                                    current_Time = Time(current_segment_start_idx: current_segment_stop_idx);
                                    current_RawDataNAN = find(isnan(current_RawData));
                                    current_ReprocessedData = nan(length(current_Time), 1);
                                    if ~isempty(current_RawDataNAN)
                                        current_End_idx = current_RawDataNAN(1);
                                        current_Start_idx = current_RawDataNAN(end);
                                        if length(current_RawData) - current_End_idx + 1 <= 1.5*max_gap
                                            RawData_trim = current_RawData(1:current_End_idx);
                                            Time_trim = current_Time(1:current_End_idx);
                                            [current_ProcessedData, current_Residual] = generate_processed_data(RawData_trim, Time_trim, filt_order, fc, current_ErrorThreshold);
                                            current_ReprocessedData(1:current_End_idx) = current_ProcessedData;
                                        end
                                        if current_Start_idx <= 1.5*max_gap
                                            RawData_trim = current_RawData(current_Start_idx:end);
                                            Time_trim = current_Time(current_Start_idx:end);
                                            [current_ProcessedData, current_Residual] = generate_processed_data(RawData_trim, Time_trim, filt_order, fc, current_ErrorThreshold);
                                            current_ReprocessedData(current_Start_idx:end) = current_ProcessedData;
                                        end
                                    end
                                    if all(isnan(current_ReprocessedData))
                                        [current_ProcessedData, current_Residual] = generate_processed_data(current_RawData, current_Time, filt_order, fc, current_ErrorThreshold);
                                        current_ReprocessedData = current_ProcessedData;
                                    end
                                    current_data_copy(current_segment_start_idx: current_segment_stop_idx) = current_ReprocessedData;
                                    ResidualSegment(ns) = current_Residual;
                                end
                                ResidualSegment_wonan = ResidualSegment(~isnan(ResidualSegment));
                                ResidualSegment_avg = mean(ResidualSegment_wonan);
                                dataPresent_idx = find(~isnan(current_data));
                                current_data_copy(1: dataPresent_idx(1)) = nan;
                                current_data_copy(dataPresent_idx(end): end) = nan;
                                ReprocessedData = current_data_copy;
                                Residual_current = ResidualSegment_avg;
                            end
                            
                            [manual_DataStartIdx, manual_DataStopIdx] = manually_chop_data(j, s, m, c, processing_inputs);
                            if ~isnan(manual_DataStartIdx)
                                ReprocessedData(1:manual_DataStartIdx) = nan;
                            end
                            if ~isnan(manual_DataStopIdx)
                                ReprocessedData(manual_DataStopIdx:end) = nan;
                            end

                            outliers_info(outlier_idx).ReprocessedData = ReprocessedData;
                            outliers_info(outlier_idx).Residual = Residual_current;
                            outliers_info(outlier_idx).OldResidual = Residual;
                            if OutlierPlot_flag
                                figure
                                broken_str = strsplit(coord_list(c));
                                marker_name_whiteSpace = strrep(current_marker_name, "_", " ");
                                plot(Time, activity_processed(j).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)), 'b', 'DisplayName','Processed Data');
                                hold on
                                plot(Time, activity(j).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)), 'r--', 'DisplayName','Raw Data');
                                if isnan(Residual_current)
                                    plotTitle = "PltNo " + string(outlier_idx) + " subNo " + string(s) + " " + marker_name_whiteSpace + broken_str(end) + " , res val = " + string(round(Residual,3)) + " mm";
                                else
                                    plotTitle = "PltNo " + string(outlier_idx) + " subNo " + string(s) + " " + marker_name_whiteSpace + broken_str(end) + " , res val = " + string(round(Residual,3)) + " mm" + " , res val new = " + string(round(Residual_current,3));
                                    plot(Time, ReprocessedData, 'g', 'DisplayName','Reprocessed Data');
                                end
                                grid
                                xlabel("Time (s)");
                                ylabel("Position (mm)");
                                legend;
                                fig = gcf;
                                fig.WindowState = 'maximized';
                                title(plotTitle);
                                figfile_name = "OutlierPlots\OutlierPlot_PltNo_" + string(outlier_idx) + "_subNo_" + string(s) + "_" + current_marker_name + coord_list(c) + "_ActivityIdx_" + string(j);
                                figfile_name_old = "ProcessedData_Plots\ProcessedDataPlot_subNo_" + string(s) + "_" + current_marker_name + "_" + coord_list(c);
                                saveas(fig, figfile_name, "jpg");
                                saveas(fig, figfile_name_old, "jpg");
                                close(fig);
                            end
                            if Residual_current < Residual
                                activity_reprocessed(j).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)) = ReprocessedData;
                                activity_reprocessed(j).VICON(s, 1).markers.residual(m, c) = Residual_current;
                            end
                        end
                    end
                    
                end
            end
            
        end
        
    end
end

function activity_new = save_processed_data(activity, processing_inputs)
    activity_new = activity;

    activity_start_idx = processing_inputs.activity_start_idx_val;
    activity_end_idx = processing_inputs.activity_end_idx_val;
    gap_plot_flag = processing_inputs.gap_plot_flag_val;
    processed_data_plot_flag = processing_inputs.processed_data_plot_flag_val;
    filt_order = processing_inputs.filt_order_val;
    fc = processing_inputs.fc_val;
    coord_list = processing_inputs.coord_list_val;
    error_threshold = processing_inputs.error_threshold_val;
    
    for k=activity_start_idx:activity_end_idx
        for s=1:25
            marker_name_list = activity(k).VICON(s, 1).markers.names;
            marker_present_list = activity(k).VICON(s, 1).markers.present;
            num_markers = length(marker_name_list);
            residual_markers = nan(num_markers,3);
            for n=1:num_markers
                Dval = [s n];
                disp(Dval);
                current_marker_name = marker_name_list(n);
                current_marker_present = marker_present_list(n);
                if current_marker_present
                    Time = activity(k).VICON(s, 1).movement_data.time;
                    Residual = nan(3,1);
                    for c=1:length(coord_list)
                        current_ErrorThreshold = error_threshold(c);
                        current_data = activity(k).VICON(s, 1).movement_data.(current_marker_name + coord_list(c));
                        current_data_copy = current_data;
                        [segmentStartIdx, segmentStopIdx] = gap_detector(s, current_data, gap_plot_flag);
                        num_segments = length(segmentStartIdx);
                        ResidualSegment = nan(num_segments,1);
                        for j=1:num_segments
                            current_segment_start_idx = segmentStartIdx(j);
                            current_segment_stop_idx = segmentStopIdx(j);
                            current_RawData = current_data(current_segment_start_idx: current_segment_stop_idx);
                            current_Time = Time(current_segment_start_idx: current_segment_stop_idx);
                            [current_ProcessedData, current_Residual] = generate_processed_data(current_RawData, current_Time, filt_order, fc, current_ErrorThreshold);
                            current_data_copy(current_segment_start_idx: current_segment_stop_idx) = current_ProcessedData;
                            ResidualSegment(j) = current_Residual;
                        end
                        ResidualSegment_wonan = ResidualSegment(~isnan(ResidualSegment));
                        ResidualSegment_avg = mean(ResidualSegment_wonan);
                        Residual(c) = ResidualSegment_avg;
                        dataPresent_idx = find(~isnan(current_data));
                        current_data_copy(1: dataPresent_idx(1)) = nan;
                        current_data_copy(dataPresent_idx(end): end) = nan;
                        activity_new(k).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)) = current_data_copy;
                        if processed_data_plot_flag
                            figure
                            broken_str = strsplit(coord_list(c));
                            marker_name_whiteSpace = strrep(current_marker_name, "_", " ");
                            plotTitle = "subNo " + string(s) + " " + marker_name_whiteSpace + broken_str(end) + " , res val = " + string(round(Residual(c),3)) + " mm";
                            plot(Time, activity_new(k).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)), 'b', 'DisplayName','Processed Data');
                            hold on
                            plot(Time, activity(k).VICON(s, 1).movement_data.(current_marker_name + coord_list(c)), 'r--', 'DisplayName','Raw Data');
                            grid
                            xlabel("Time (s)");
                            ylabel("Position (mm)");
                            legend;
                            fig = gcf;
                            fig.WindowState = 'maximized';
                            title(plotTitle);
                            figfile_name = "ProcessedData_Plots\ProcessedDataPlot_subNo_" + string(s) + "_" + current_marker_name + "_" + coord_list(c);
                            saveas(fig, figfile_name, "jpg");
                            close(fig);
                        end
                    end
                    residual_markers(n,:) = transpose(Residual);
                end
            end
            activity_new(k).VICON(s, 1).markers.residual = residual_markers;
        end
    end

end

function [segmentStartIdx, segmentStopIdx] = gap_detector(subject_idx, data, gap_plot_flag)

    nanStartIdx = find(diff(isnan([0;data;0]))==1);
    nanStopIdx = find(diff(isnan([0;data;0]))==-1);
    
    gapMax = 10;
    segmentMin = 13;

    gapStartIdx = nanStartIdx;
    gapStopIdx = nanStopIdx-1;
    
    gapStartIdx_copy = gapStartIdx;
    gapStopIdx_copy = gapStopIdx;
    
    gapNum = nan;
    if ~isempty(gapStartIdx)
        gapNum = length(gapStartIdx);
    end
    if ~isnan(gapNum)
        for k=1:gapNum
            gapSizeNow = gapStopIdx(k) - gapStartIdx(k) + 1;
            if gapSizeNow < gapMax
                gapStartIdx_copy(k) = nan;
                gapStopIdx_copy(k) = nan;
            end
        end
    end
    gapStartIdx_copy = gapStartIdx_copy(~isnan(gapStartIdx_copy));
    gapStopIdx_copy = gapStopIdx_copy(~isnan(gapStopIdx_copy));
    
    idx_NotNaN = find(~isnan(data));
    
    segmentStartIdx = [];
    if length(data) > 1
        if or(~isnan(data(1)), idx_NotNaN(1) <= gapMax)
            segmentStartIdx = [segmentStartIdx; 1];
        end
    end
    segmentStartIdx = [segmentStartIdx; gapStopIdx_copy + 1];
    
    segmentStopIdx = gapStartIdx_copy - 1;
    if length(data) > 1
        if or(~isnan(data(end)), length(data) - idx_NotNaN(end) + 1 < gapMax)
            segmentStopIdx = [segmentStopIdx; length(data)];
        else
            segmentStopIdx = [segmentStopIdx; idx_NotNaN(end)];
        end
    end
    
    if ~isempty(segmentStopIdx)
        if segmentStopIdx(1) == 0
            segmentStopIdx(1) = [];
        end
    end
    
    segmentStartIdx_copy = segmentStartIdx;
    segmentStopIdx_copy = segmentStopIdx;
    segmentNum = nan;
    if ~isempty(segmentStartIdx)
        segmentNum = length(segmentStartIdx);
    end
    if ~isnan(segmentNum)
        for k=1:segmentNum
            segmentSizeNow = segmentStopIdx(k) - segmentStartIdx(k) + 1;
            if segmentSizeNow < segmentMin
                segmentStartIdx_copy(k) = nan;
                segmentStopIdx_copy(k) = nan;
            end
        end
    end
    segmentStartIdx_copy = segmentStartIdx_copy(~isnan(segmentStartIdx_copy));
    segmentStopIdx_copy = segmentStopIdx_copy(~isnan(segmentStopIdx_copy));
    
    segmentStartIdx = segmentStartIdx_copy;
    segmentStopIdx = segmentStopIdx_copy;

    if gap_plot_flag

        for m=1:length(segmentStartIdx)
            x_min_lt = max(1, segmentStartIdx(m)-20);
            x_max_lt = min(segmentStopIdx(m)+20, length(data));
            ax = cla(); 
            plot(~isnan(data), '-k', 'LineWidth', 2)
            xlim([x_min_lt, x_max_lt])
            ylim([-.2, 1.5])
            xlabel('index')
            ylabel('segment')
            set(ax, 'YTick', [0,1], 'YTickLabel', {'False', 'True'})
            hold on
            arrayfun(@(start)xline(start, 'b-','Start'), segmentStartIdx)
            arrayfun(@(stop)xline(stop, 'c-','Stop'), segmentStopIdx)
            color = [0, 0, 1];
            for k=1:length(segmentStartIdx)
                shade_x = [segmentStartIdx(k), segmentStartIdx(k), segmentStopIdx(k), segmentStopIdx(k)];
                shade_y = [0, 1, 1, 0];
                a = fill(shade_x, shade_y, color);
                a.FaceAlpha = 0.1;
            end
            fig = gcf;
            fig.WindowState = 'maximized';
            figfile_name = "GapDetect_Plots\GapDetectPlot_subNo_" + string(subject_idx) + "_segNo_" + string(m);
            saveas(fig, figfile_name, "jpg");
            close(fig);
        end
        
    end

end

function [ProcessedData, Residual] = generate_processed_data(RawData, Time, filt_order, fc, error_threshold)
    num_measurements = length(Time);
    ProcessedData = nan(num_measurements, 1);
    Residual = nan;

    if ~all(isnan(RawData))
        pp = spline(Time, RawData);
        InterpolatedData = ppval(pp, Time);
        fs = 100;
        [b,a] = butter(filt_order,fc/(fs/2),'low');
        DataOut_zeroshift = filtfilt(b,a,InterpolatedData);
        res_vector = DataOut_zeroshift - RawData;
        res_vector_wonan = res_vector(~isnan(res_vector));
        error_val = std(res_vector_wonan);
        if error_val < error_threshold
            ProcessedData = DataOut_zeroshift;
            Residual = error_val;
        end
    end

end

function remove_idx = manually_remove_outliers
    total_outliers = 236;
    remove_idx = struct;
    
    for k=1:total_outliers
        remove_idx(k).idxVal = k;
        remove_idx(k).listVal = [];
    end
    
    remove_idx(21).listVal = [127];
    remove_idx(28).listVal = [6];
    remove_idx(29).listVal = [6];
    remove_idx(30).listVal = [294];
    remove_idx(31).listVal = [294];
    remove_idx(32).listVal = [294];
    remove_idx(33).listVal = [408];
    remove_idx(51).listVal = 9:13;
    remove_idx(52).listVal = 9:13;
    remove_idx(56).listVal = [279];
    remove_idx(57).listVal = [279];
    remove_idx(62).listVal = 269:271;
    remove_idx(63).listVal = 269:271;
    remove_idx(78).listVal = 272:288;
    remove_idx(79).listVal = [272];
    remove_idx(80).listVal = 262:292;
    remove_idx(93).listVal = 1:18;
    remove_idx(96).listVal = [197];
    remove_idx(100).listVal = 1:13;
    remove_idx(124).listVal = [281];
    remove_idx(125).listVal = [281];
    remove_idx(126).listVal = [281];
    remove_idx(127).listVal = [40];
    remove_idx(133).listVal = [248];
    remove_idx(149).listVal = 293:307;
    remove_idx(153).listVal = 1:6;
    remove_idx(154).listVal = 1:6;
    remove_idx(160).listVal = 305:312;
    remove_idx(161).listVal = 305:312;
    remove_idx(162).listVal = 305:312;
    remove_idx(175).listVal = 325:327;
    remove_idx(176).listVal = 224:226;
    remove_idx(188).listVal = 306:314;
    remove_idx(189).listVal = 306:314;
    remove_idx(190).listVal = 283:288;
    remove_idx(192).listVal = 350:355;
    remove_idx(193).listVal = 350:355;
    remove_idx(199).listVal = 319:326;
    remove_idx(200).listVal = 326:362;
    remove_idx(204).listVal = 332:346;
    remove_idx(233).listVal = [10];
    remove_idx(234).listVal = [10];
    remove_idx(236).listVal = 369:380;
    
end

function [manual_DataStartIdx, manual_DataStopIdx] = manually_chop_data(j_val, s_val, m_val, c_val, processing_inputs)
    activity_start_idx = processing_inputs.activity_start_idx_val;
    activity_end_idx = processing_inputs.activity_end_idx_val;
    manual_data = struct;
    for j=activity_start_idx:activity_end_idx
        for s=1:25
            for m=1:70
                for c=1:3
                    manual_data(j,s,m,c).BoolVal = false;
                    manual_data(j,s,m,c).StartIdxVal = nan;
                    manual_data(j,s,m,c).StopIdxVal = nan;
                end
            end
        end
    end
    
    manual_data(5,1,3,3).StartIdxVal = nan;
    manual_data(5,1,3,3).StopIdxVal = 291;
    
    manual_data(6,14,12,3).StartIdxVal = nan;
    manual_data(6,14,12,3).StopIdxVal = 247;
    
    manual_data(7,1,25,3).StartIdxVal = nan;
    manual_data(7,1,25,3).StopIdxVal = 282;
    
    manual_data(8,2,23,1).StartIdxVal = nan;
    manual_data(8,2,23,1).StopIdxVal = 264;
    
    manual_data(10,12,7,3).StartIdxVal = 22;
    manual_data(10,12,7,3).StopIdxVal = nan;
    
    manual_data(11,6,17,3).StartIdxVal = 12;
    manual_data(11,6,17,3).StopIdxVal = nan;
    
    manual_data(11,14,5,2).StartIdxVal = 12;
    manual_data(11,14,5,2).StopIdxVal = nan;
    
    manual_data(11,14,17,2).StartIdxVal = 16;
    manual_data(11,14,17,2).StopIdxVal = nan;
    
    manual_data(11,17,25,3).StartIdxVal = nan;
    manual_data(11,17,25,3).StopIdxVal = 278;
    
    manual_data(11,20,25,2).StartIdxVal = nan;
    manual_data(11,20,25,2).StopIdxVal = 344;
    
    manual_data(12,8,39,2).StartIdxVal = 12;
    manual_data(12,8,39,2).StopIdxVal = nan;
    
    manual_data(12,8,39,3).StartIdxVal = 17;
    manual_data(12,8,39,3).StopIdxVal = nan;
    
    manual_data(12,10,11,3).StartIdxVal = nan;
    manual_data(12,10,11,3).StopIdxVal = 212;
    
    manual_data(13,1,3,3).StartIdxVal = nan;
    manual_data(13,1,3,3).StopIdxVal = 308;
    
    manual_data(13,19,8,3).StartIdxVal = 17;
    manual_data(13,19,8,3).StopIdxVal = nan;
    
    manual_DataStartIdx = manual_data(j_val, s_val, m_val, c_val).StartIdxVal;
    manual_DataStopIdx = manual_data(j_val, s_val, m_val, c_val).StopIdxVal;
end

function [manual_segmentStartIdx, manual_segmentStopIdx, manual_segmentBoolVal] = manually_create_segment(j_val, s_val, m_val, c_val, processing_inputs)
    activity_start_idx = processing_inputs.activity_start_idx_val;
    activity_end_idx = processing_inputs.activity_end_idx_val;
    manual_segment = struct;
    for j=activity_start_idx:activity_end_idx
        for s=1:25
            for m=1:70
                for c=1:3
                    manual_segment(j,s,m,c).BoolVal = false;
                    manual_segment(j,s,m,c).StartIdxVal = nan;
                    manual_segment(j,s,m,c).StopIdxVal = nan;
                end
            end
        end
    end
    
    manual_segment(12,13,20,2).BoolVal = true;
    manual_segment(12,13,20,2).StartIdxVal = [1;200];
    manual_segment(12,13,20,2).StopIdxVal = [192;261];
    
    manual_segment(12,13,20,3).BoolVal = true;
    manual_segment(12,13,20,3).StartIdxVal = [1;200];
    manual_segment(12,13,20,3).StopIdxVal = [192;261];
    
    manual_segment(12,18,23,3).BoolVal = true;
    manual_segment(12,18,23,3).StartIdxVal = [1;219];
    manual_segment(12,18,23,3).StopIdxVal = [209;282];
    
    manual_segmentStartIdx = manual_segment(j_val, s_val, m_val, c_val).StartIdxVal;
    manual_segmentStopIdx = manual_segment(j_val, s_val, m_val, c_val).StopIdxVal;
    manual_segmentBoolVal = manual_segment(j_val, s_val, m_val, c_val).BoolVal;
end