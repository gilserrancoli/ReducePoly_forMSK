% % List of movements
% clear all;
% close all;
% clc
% % movements = { 'ngait_og1', 'bouncy4', 'mtpgait3', 'ngait_tm_fast1', 'ngait_og5', 'bouncy7','mtpgait9','ngait_tm_set1'};
% movements = { 'ngait_og1'};
% 
% % Threshold values
% % thresholds = [2, 1, 0.5];
% thresholds=[2];
% femfaces_values = [171];
% 
% % for i = 1:length(movements)
% %     nametrial_id = movements{i};
% %     for j = 1:length(thresholds) 
% %         for useReducedPolynomials = [0]
% %             err_poly = thresholds(j);
% %                 for femfaces = femfaces_values  % Iterate over femfaces values (171 and 258)
% % 
% %                 % Create Options structure and assign femfaces value
% %                 Options.nfacesFem = femfaces;
% % 
% %                 TrackSim_3D_GC_v2(nametrial_id,useReducedPolynomials,err_poly,Options);
% %                 end
% % 
% %         end
% %     end
% % end
% for i = 1:length(movements)
%     nametrial_id = movements{i};
%     for j = 1:length(thresholds) 
%         for useReducedPolynomials = [1]
% 
%             err_poly = thresholds(j); % Keep err_poly assignment before the femfaces loop
% 
%             for femfaces = femfaces_values  % Iterate over femfaces values (171 and 258)
% 
%                 % Create Options structure and assign femfaces value
%                 Options.nfacesFem = femfaces;  %  Add femfaces to Options struct
%                 Options.useReducedPolynomials = useReducedPolynomials;
%                 Options.err_poly = err_poly;
% 
%                 % Call the function with Options struct
%                 TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options);
% 
%             end
%         end
%     end
% end
clear all;
close all;
clc

movements = {'mtpgait9'};
thresholds=[1 2];
femfaces_values = [188];
% 
% Define ranges of weights to try
% Qs_weights_to_try  = [50,40,30];
% KCF_weights_to_try = [20,25,30];
% GRF_weights_to_try = [20,25,30];
% MsAc_weights_to_try =10;
% % Initialize result storage
% results = [];
% trial_number = 0;
% 
% 
% for i = 1:length(movements)
%     nametrial_id = movements{i};
%     for j = 1:length(thresholds) 
%         for useReducedPolynomials = [1]
%             err_poly = thresholds(j);
%             for femfaces = femfaces_values
%                 % Set options
%                 Options.nfacesFem = femfaces;
%                 Options.useReducedPolynomials = useReducedPolynomials;
%                 Options.err_poly = err_poly;
% 
%                 % Create output folder
%                 outputFolder = fullfile('Results', nametrial_id);
%                 if ~exist(outputFolder, 'dir')
%                     mkdir(outputFolder);
%                 end
% 
% 
%                 % Tuning loop
%                 for w1 = 1:length(Qs_weights_to_try)
%                     for w2 = 1:length(KCF_weights_to_try)
%                         for w3 = 1:length(GRF_weights_to_try)
%                             for w4 = 1:length(MsAc_weights_to_try)
%                                 trial_number = trial_number + 1;
%                                 % Check skip condition
%                                 if Qs_weights_to_try(w1) == 20 && ...
%                                    KCF_weights_to_try(w2) == 30 && ...
%                                    GRF_weights_to_try(w3) == 30
%                                     continue; % Skip this combination
%                                 end
% 
%                                 % Set Weights
%                                 W.Qs = Qs_weights_to_try(w1);
%                                 W.KCF = KCF_weights_to_try(w2);
%                                 W.GRF = GRF_weights_to_try(w3);
%                                 W.a = MsAc_weights_to_try(w4);
%                                 W.Qdots = 10;
%                                 W.GRM = 10;
%                                 W.ID_act = 0;
%                                 W.minPelvisRes = 2;
%                                 W.u = 0.03;
%                                 W.u_qd2dot = 0.003;
%                                 W.u_vA = 0.52;
% 
%                                 % Unique save name
%                                 savename_suffix = sprintf('_Qs%d_KCF%d_a%d_GRF%d_T%d', ...
%                                     W.Qs, W.KCF, W.a, W.GRF, trial_number);
% 
%                                 % Run simulation
%                                 Results_3D = TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);
% 
%                                 % Save simulation output
%                                 save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
%                                 save(save_filename, 'Results_3D');
% 
%                                 % Extract tracking signals
%                                 exp_hip_flex = Results_3D.NMesh_50.Qs_toTrack(:,10);
%                                 sim_hip_flex = Results_3D.Simulated.Qs_opt(:,10);
%                                 exp_hip_add = Results_3D.NMesh_50.Qs_toTrack(:,11);
%                                 sim_hip_add = Results_3D.Simulated.Qs_opt(:,11);
%                                 exp_hip_rot = Results_3D.NMesh_50.Qs_toTrack(:,12);
%                                 sim_hip_rot = Results_3D.Simulated.Qs_opt(:,12);
%                                 exp_knee_flex = Results_3D.NMesh_50.Qs_toTrack(:,14);
%                                 sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);
% 
%                                 % Calculate RMSE and R2
%                                 rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
%                                 r2 = @(y, yhat) 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
% 
%                                 % Store metrics
%                                 results(trial_number).TrialID = trial_number;
%                                 results(trial_number).W_Qs = W.Qs;
%                                 results(trial_number).W_KCF = W.KCF;
%                                 results(trial_number).W_GRF = W.GRF;
%                                 results(trial_number).W_a = W.a;
%                                 results(trial_number).RMSE_HipFlex = rmse(exp_hip_flex, sim_hip_flex);
%                                 results(trial_number).R2_HipFlex = r2(exp_hip_flex, sim_hip_flex);
%                                 results(trial_number).RMSE_HipAdd = rmse(exp_hip_add, sim_hip_add);
%                                 results(trial_number).R2_HipAdd = r2(exp_hip_add, sim_hip_add);
%                                 results(trial_number).RMSE_HipRot = rmse(exp_hip_rot, sim_hip_rot);
%                                 results(trial_number).R2_HipRot = r2(exp_hip_rot, sim_hip_rot);
%                                 results(trial_number).RMSE_Knee = rmse(exp_knee_flex, sim_knee_flex);
%                                 results(trial_number).R2_Knee = r2(exp_knee_flex, sim_knee_flex);
%                                 results(trial_number).SaveName = save_filename;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % After all simulations
% ResultsTable = struct2table(results);
% disp(ResultsTable);
% sortedTable = sortrows(ResultsTable, {'RMSE_HipFlex','RMSE_HipAdd','RMSE_HipRot','RMSE_Knee'}, 'ascend');
% % disp('Top 50 Best Trials:');
% % disp(sortedTable(1:50,:));
% Define specific weight combinations
% weight_combinations = [
%     20, 30, 25;
%     30, 30, 40;
%     30, 35, 40;
%     40, 25, 30
% ];

weight_combinations = [
    20, 40, 25
    ]; %20, 50, 25
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
