
clear all;
close all;
clc

movements = {'ngait5','mtpgait3','tmfast1','bouncy7'};
thresholds=[1 2];
femfaces_values = [188];

weight_combinations = [
    20, 40, 25
    ]; 
MsAc_weights_to_try = 10;  % Assuming this is constant
trial_number = 0;

for i = 1:length(movements)
    nametrial_id = movements{i};
    for j = 1:length(thresholds)
        for useReducedPolynomials = [0 1]
            err_poly = thresholds(j);
            for femfaces = femfaces_values
                % Set options
                Options.nfacesFem = femfaces;
                Options.useReducedPolynomials = useReducedPolynomials;
                Options.err_poly = err_poly;

                % Create output folder
                outputFolder = fullfile('Results', nametrial_id);
                if ~exist(outputFolder, 'dir')
                    mkdir(outputFolder);
                end

                % Loop through defined weight combinations
                for k = 1:size(weight_combinations, 1)
                    trial_number = trial_number + 1;

                    W.Qs  = weight_combinations(k, 1);
                    W.KCF = weight_combinations(k, 2);
                    W.GRF = weight_combinations(k, 3);
                    W.a   = MsAc_weights_to_try;

                    W.Qdots = 10;
                    W.GRM = 10;
                    W.ID_act = 0;
                    W.minPelvisRes = 0.2;
                    W.u = 0.03;
                    W.u_qd2dot = 0.003;
                    W.u_qd2dot_kneesecdof=50;
                    W.u_vA = 0.52;

                    % Unique save name
                    savename_suffix = sprintf('_Qs%d_KCF%d_a%d_GRF%d_T%d', ...
                        W.Qs, W.KCF, W.a, W.GRF, trial_number);

                    % Run simulation
                    Results_3D = TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);

                    % Save simulation output
                    save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
                    save(save_filename, 'Results_3D');

                    % Extract and compute RMSE & R2
                    exp_hip_flex = Results_3D.NMesh_40.Qs_toTrack(:,10);
                    sim_hip_flex = Results_3D.Simulated.Qs_opt(:,10);
                    exp_hip_add  = Results_3D.NMesh_40.Qs_toTrack(:,11);
                    sim_hip_add  = Results_3D.Simulated.Qs_opt(:,11);
                    exp_hip_rot  = Results_3D.NMesh_40.Qs_toTrack(:,12);
                    sim_hip_rot  = Results_3D.Simulated.Qs_opt(:,12);
                    exp_knee_flex = Results_3D.NMesh_40.Qs_toTrack(:,14);
                    sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);

                    rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
                    r2   = @(y, yhat) 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);

                    results(trial_number).TrialID = trial_number;
                    results(trial_number).W_Qs = W.Qs;
                    results(trial_number).W_KCF = W.KCF;
                    results(trial_number).W_GRF = W.GRF;
                    results(trial_number).W_a = W.a;
                    results(trial_number).RMSE_HipFlex = rmse(exp_hip_flex, sim_hip_flex);
                    results(trial_number).R2_HipFlex   = r2(exp_hip_flex, sim_hip_flex);
                    results(trial_number).RMSE_HipAdd  = rmse(exp_hip_add, sim_hip_add);
                    results(trial_number).R2_HipAdd    = r2(exp_hip_add, sim_hip_add);
                    results(trial_number).RMSE_HipRot  = rmse(exp_hip_rot, sim_hip_rot);
                    results(trial_number).R2_HipRot    = r2(exp_hip_rot, sim_hip_rot);
                    results(trial_number).RMSE_Knee    = rmse(exp_knee_flex, sim_knee_flex);
                    results(trial_number).R2_Knee      = r2(exp_knee_flex, sim_knee_flex);
                    results(trial_number).SaveName     = save_filename;
                end
            end
        end
    end
end

% Convert results to table and sort
ResultsTable = struct2table(results);
disp(ResultsTable);
sortedTable = sortrows(ResultsTable, {'RMSE_HipFlex','RMSE_HipAdd','RMSE_HipRot','RMSE_Knee'}, 'ascend');
