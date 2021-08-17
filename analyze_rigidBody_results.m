%% LOAD DATA

load('reference_frame_data_ac_file.mat');
load('reference_frame_data_mvt_file.mat');
load('UVM_data_regrouped.mat');

%% SET PREFERENCES

activity_idx = 5;
overwrite_plot_folder = false;

%% CREATE DATA CONTAINERS

num_subjects = 25;
plt_width = 5;

flexionMax_list = nan(num_subjects, 1);
adductionMax_list = nan(num_subjects, 1);
externalRotationMax_list = nan(num_subjects, 1);

flexionRange_list = nan(num_subjects, 1);
adductionRange_list = nan(num_subjects, 1);
externalRotationRange_list = nan(num_subjects, 1);

flexionMeanVal_list = nan(num_subjects, 1);
adductionMeanVal_list = nan(num_subjects, 1);
externalRotationMeanVal_list = nan(num_subjects, 1);

flexionMean_list = nan(num_subjects, 1);
adductionMean_list = nan(num_subjects, 1);
externalRotationMean_list = nan(num_subjects, 1);

noData_idx = [20];

if or(~isfolder("MVT_RigidBody_Error_Plots"), overwrite_plot_folder)
    mkdir MVT_RigidBody_Error_Plots
end

%% POPULATE DATA CONTAINERS

for j=1:num_subjects
    
    flexionMean_list(j) = mean(subject_reference_frame_data_ac(j).knee_angle_data.processed.right.flexion, 'omitnan');
    adductionMean_list(j) = mean(subject_reference_frame_data_ac(j).knee_angle_data.processed.right.adduction, 'omitnan');
    externalRotationMean_list(j) = mean(subject_reference_frame_data_ac(j).knee_angle_data.processed.right.externalRotation, 'omitnan');
    
    if ismember(j, noData_idx)
        flexionMax_list(j) = nan;
        adductionMax_list(j) = nan;
        externalRotationMax_list(j) = nan;
    else
        flexionMax_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion);
        adductionMax_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.adduction);
        externalRotationMax_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.externalRotation);
        
        flexionRange_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion) - min(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion);
        adductionRange_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.adduction) - min(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.adduction);
        externalRotationRange_list(j) = max(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.externalRotation) - min(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.externalRotation);
        
        flexionMeanVal_list(j) = mean(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion, 'omitnan');
        adductionMeanVal_list(j) = mean(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.adduction, 'omitnan');
        externalRotationMeanVal_list(j) = mean(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.externalRotation, 'omitnan');
    end
    
end

%% CORRELATION WITH GREATER TROCHANTER

flexionX = [flexionMean_list,adductionMean_list,externalRotationMean_list];
flexionMaxLM = fitlm(flexionX,flexionMax_list);
flexionRangeLM = fitlm(flexionX,flexionRange_list);
flexionMeanValLM = fitlm(flexionX,flexionMeanVal_list);
flexionMaxGT_bool = any(flexionMaxLM.Coefficients.pValue(2:end) < 0.05);
if flexionMaxGT_bool
    disp("-----Max Flexion related to offset caused by GT marker-----")
end
flexionRangeGT_bool = any(flexionRangeLM.Coefficients.pValue(2:end) < 0.05);
if flexionRangeGT_bool
    disp("-----Range Flexion related to offset caused by GT marker-----")
end
flexionMeanValGT_bool = any(flexionMeanValLM.Coefficients.pValue(2:end) < 0.05);
if flexionMeanValGT_bool
    disp("-----MeanVal Flexion related to offset caused by GT marker-----")
end

adductionX = [flexionMean_list,adductionMean_list,externalRotationMean_list];
adductionMaxLM = fitlm(adductionX,adductionMax_list);
adductionRangeLM = fitlm(adductionX,adductionRange_list);
adductionMeanValLM = fitlm(adductionX,adductionMeanVal_list);
adductionMaxGT_bool = any(adductionMaxLM.Coefficients.pValue(2:end) < 0.05);
if adductionMaxGT_bool
    disp("-----Max Adduction related to offset caused by GT marker-----")
end
adductionRangeGT_bool = any(adductionRangeLM.Coefficients.pValue(2:end) < 0.05);
if adductionRangeGT_bool
    disp("-----Range Adduction related to offset caused by GT marker-----")
end
adductionMeanValGT_bool = any(adductionMeanValLM.Coefficients.pValue(2:end) < 0.05);
if flexionMeanValGT_bool
    disp("-----MeanVal Adduction related to offset caused by GT marker-----")
end

adductionMaxLM_predCoeff = fitlm(adductionMean_list,adductionMax_list).Coefficients.Estimate;
adductionMaxLM_pred = adductionMaxLM_predCoeff(1) + adductionMaxLM_predCoeff(2)*adductionMean_list;
plot(adductionMean_list, adductionMaxLM_pred, 'r')
hold on
scatter(adductionMean_list,adductionMax_list, 'b')
grid
title("Correlation of RigidBody Results with offset due to GT")
xlabel("Mean Adduction Static Calibration (degrees)")
ylabel("Max Adduction Walking Trial (degrees)")
fig = gcf;
fig.WindowState = 'maximized';
saveas(fig,"MVT_RigidBody_Error_Plots\GT_correlation_plot", "jpg");
close(fig);

externalRotationX = [flexionMean_list,adductionMean_list,externalRotationMean_list];
externalRotationMaxLM = fitlm(externalRotationX,externalRotationMax_list);
externalRotationRangeLM = fitlm(externalRotationX,externalRotationRange_list);
externalRotationMeanValLM = fitlm(externalRotationX,externalRotationMeanVal_list);
externalRotationMaxGT_bool = any(externalRotationMaxLM.Coefficients.pValue(2:end) < 0.05);
if externalRotationMaxGT_bool
    disp("-----Max ExternalRotation related to offset caused by GT marker-----")
end
externalRotationRangeGT_bool = any(externalRotationRangeLM.Coefficients.pValue(2:end) < 0.05);
if externalRotationRangeGT_bool
    disp("-----Range ExternalRotation related to offset caused by GT marker-----")
end
externalRotationMeanValGT_bool = any(externalRotationMeanValLM.Coefficients.pValue(2:end) < 0.05);
if externalRotationMeanValGT_bool
    disp("-----MeanVal ExternalRotation related to offset caused by GT marker-----")
end

%% SYSTEMATIC PERIODIC KNEE ANGLE ERROR PLOTS

for j=1:num_subjects
    
    if ~ismember(j, noData_idx)
        time = activity(activity_idx).VICON(j, 1).movement_data.time;

        flexionError = subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion - subject_reference_frame_data_mvt(j).knee_angle_data.processed.right.flexion;
        adductionError = subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.adduction - subject_reference_frame_data_mvt(j).knee_angle_data.processed.right.adduction;
        externalRotationError = subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.externalRotation - subject_reference_frame_data_mvt(j).knee_angle_data.processed.right.externalRotation;
        
        if length(flexionError) > 5
            [pks_FlexionError, locs_FlexionError] = findpeaks(flexionError);
            [pks_FlexionRB, locs_FlexionRB] = findpeaks(subject_reference_frame_data_mvt(j).knee_angle_data.rigidBody.right.flexion);
        else
            pks_FlexionError = [];
            locs_FlexionError = [];
            pks_FlexionRB = [];
            locs_FlexionRB = [];
        end
        

        subplot(3,1,1)
        plot(time, flexionError, 'b', 'DisplayName', 'right flexionError');
        hold on
        scatter(time(locs_FlexionError), pks_FlexionError, 'r', 'DisplayName', 'right flexionErrorPeaks');
        arrayfun(@(FlexPeak)xline(FlexPeak, 'black-'), time(locs_FlexionRB));
        max_ht = max(flexionError) + plt_width;
        min_ht = min(flexionError) - plt_width;
        plt_limit = [min_ht, max_ht];
        if all(~isnan(plt_limit))
            ylim(plt_limit);
        end
        description_axes = zeros(3, 1);
        description_axes(1) = plot(NaN,NaN,'-b');
        description_axes(2) = plot(NaN,NaN,'-r');
        description_axes(3) = plot(NaN,NaN,'-black');
        legend(description_axes, 'FlexionError: blue','FlexionErrorPeaks: red','FlexionRigidBodyPeakLines: black');
        grid
        xlabel("time (s)")
        ylabel("angle (degrees)")
        title("Right Knee Flexion Error")

        subplot(3,1,2)
        plot(time, adductionError, 'b', 'DisplayName', 'right adductionError')
        max_ht = max(adductionError) + plt_width;
        min_ht = min(adductionError) - plt_width;
        plt_limit = [min_ht, max_ht];
        if all(~isnan(plt_limit))
            ylim(plt_limit);
        end
        grid
        legend
        xlabel("time (s)")
        ylabel("angle (degrees)")
        title("Right Knee Adduction Error")

        subplot(3,1,3)
        plot(time, externalRotationError, 'b', 'DisplayName', 'right externalRotationError')
        max_ht = max(externalRotationError) + plt_width;
        min_ht = min(externalRotationError) - plt_width;
        plt_limit = [min_ht, max_ht];
        if all(~isnan(plt_limit))
            ylim(plt_limit);
        end
        grid
        legend
        xlabel("time (s)")
        ylabel("angle (degrees)")
        title("Right Knee ExternalRotation Error")

        suptitle("Rigid Body Error Plots")
        fig = gcf;
        fig.WindowState = 'maximized';
        saveas(fig,"MVT_RigidBody_Error_Plots\RigidBody_error_subNo_" + string(j), "jpg");
        close(fig);
        
    end
end