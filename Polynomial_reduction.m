% This function takes the full polynomials that parameterize lMT and moment
% arms and reduce the number of coefficients systematically, based on the
% statistical significance.

%Authors: Mohanad Harba and Gil Serrancolí
%Date: 01/10/2024
%

clearvars -except threshold_r threshold_t;
close all;
clc;
tic;
current_folder = pwd;
load MuscleData_subject_GC.mat

% Load the appropriate file based on the threshold
if (threshold_r == 0.0005) && (threshold_t==0.005)
    load('MuscleInfo_full_05.mat');
elseif (threshold_r == 0.001) && (threshold_t==0.01)
    load('MuscleInfo_full_1.mat');
elseif (threshold_r == 0.002) && (threshold_t==0.02)
    load('MuscleInfo_full_2.mat');
else
    error('Invalid threshold value. Please select 0.0005 mm / 0.005, 0.001 mm / 0.01, or 0.002 mm / 0.02');
end

load muscle_spanning_joint_INFO_subject_GC.mat

muscle_sel = 1:50;

muscles_with_low_error = [];
muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO <= 0.0001 & muscle_spanning_joint_INFO >= -0.0001) = 0;
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO ~= 0) = 1;
muscles_satisfy = 0;
muscles_not_satisfy = 0;

% Initialize total counters for all muscles
total_constant_count = 0;
total_linear_count = 0;
total_quadratic_single_count = 0;
total_quadratic_interaction_count = 0;
total_cubic_single_count = 0;
total_cubic_interaction_count = 0;
total_biquadratic_single_count = 0;
total_biquadratic_interaction_count = 0;

% Counter for total coefficients across all muscles
total_coefficients_count = 0;

% Initialize counters for dropped terms across all muscles
total_dropped_constant_count = 0;
total_dropped_linear_count = 0;
total_dropped_quadratic_only_count = 0;
total_dropped_quadratic_mixed_count = 0;
total_dropped_cubic_only_count = 0;
total_dropped_cubic_mixed_count = 0;
total_dropped_biquadratic_only_count = 0;
total_dropped_biquadratic_mixed_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the process for each muscle
for m_nr = muscle_sel
    muscle_index = muscle_sel(m_nr);
    musle_info_summary.muscle_number = muscle_index;
    order = MuscleInfo.muscle(muscle_index).order;
    Total_coefficients = length(MuscleInfo.muscle(muscle_index).coeff{order});
    dofs = MuscleInfo.muscle(muscle_index).DOF;
    muscle_name = MuscleInfo.muscle(muscle_index).m_name;
    muscle_info_summary.DoF = dofs;
    muscle_info_summary.M_name = muscle_name;
    muscle_info_summary.Total_coeff = Total_coefficients;

    index_dof_crossing = find(muscle_spanning_joint_INFO(muscle_index, :) == 1);
    [mat, diff_mat_q] = n_art_mat_9_GC(MuscleData.q(:, index_dof_crossing), order);
    coeff = MuscleInfo.muscle(muscle_index).coeff{order};

    lMT_reco{muscle_index, order} = ...
        mat * MuscleInfo.muscle(muscle_index).coeff{1, order}; 
        %lMT reconstructed from full polynomials
    lMT_matrix = lMT_reco{muscle_index, order};
    has_negative_value = any(lMT_matrix(:) < 0); 
        %to help to indentify if there was an issue with the reconstruction

    for n = 1:length(index_dof_crossing)
        dM_reco{muscle_index, order}(:, n) = ...
            (-squeeze(diff_mat_q(:, :, n))) * MuscleInfo.muscle(muscle_index).coeff{1, order}; 
            %moment arms reconstructed from full polynomials
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computation of the standard error (SE)
    % See the computation of SE from Dowd et al. for instance:
    % https://doi.org/10.1111/1475-6773.12122

    X = mat;
    Y = MuscleData.lMT(:, muscle_index); % Y contains random lMT and dM data

    for j = 1:size(diff_mat_q, 3)
        X = [X; -diff_mat_q(:, :, j)];
        Y = [Y; reshape(MuscleData.dM(:, muscle_index, index_dof_crossing(j)), [], 1)];
    end
    Y_reco = [lMT_reco{muscle_index, order}; dM_reco{muscle_index, order}(:)]; 
    % Y_reco contains lMT and dM data reconstructed

    % Calculate the size of the input matrix X. X can be seen as the 
    % jacobian matrix of the polynomial Y_reco with respect each polynomial 
    % coefficient  
    [m, p] = size(X);  % Number of frames in X values and coefficients
    Jacobian = X;  % Jacobian matrix is the same as X here

    % Calculate the transpose of the Jacobian matrix
    Jacobian_transpose = Jacobian';

    % Compute Y_hat
    y_hat = X * coeff; %y_hat is the same as Y_reco here
    SSE = sum((Y - y_hat).^2);

    jacobian_1 = pinv(Jacobian_transpose * Jacobian);
    SE = sqrt(diag(SSE / (length(Y) - length(coeff)) * jacobian_1)); 
    % SE stands for standard error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the t-values for each coefficient
    t_values = coeff ./ SE;

    % Calculate the p-values for each coefficient using the t-distribution
    df = length(X) - length(coeff);
    p_values = 2 * (1 - tcdf(abs(t_values), df));

    threshold_p = 0.05;

    % Order the coefficients by their p-value
    [p_values_sorted, index] = sort(p_values);
    coeff_sorted = coeff(index);

    % Add coefficients one by one (in the order of p-value) until the
    % threshold is satisfied
    i = 1;
    found = false;
    while (~found) && (i <= length(coeff))
        lMT = MuscleData.lMT(:, muscle_index);
        dM = squeeze(MuscleData.dM(:, m_nr, find(muscle_spanning_joint_INFO(m_nr, :))));

        % Update selected coefficients
        selected_index = index(1:i);

        npoints = size(diff_mat_q, 1);
        ndofs = length(dofs);  % Number of degrees of freedom
        diff_mat_q_all = zeros(npoints * ndofs, size(diff_mat_q, 2));

        for d = 1:ndofs
            start_idx = npoints * (d - 1) + 1;
            end_idx = npoints * d;
            diff_mat_q_all(start_idx:end_idx, :) = -diff_mat_q(:, :, d);
        end

        coeff_aux = [mat(:, selected_index); diff_mat_q_all(:, selected_index)] \ [lMT; dM(:)];
            %coefficients that best fit lMT and moment arms using only the 
            % subset of selected coefficients

        % Update the full coefficient vector with the selected ones
        coeff_2 = zeros(size(mat, 2), 1);
        coeff_2(selected_index) = coeff_aux;

        % Recalculate lMT and dM with the updated selected coefficients
        lMT_reco_2{muscle_index, order} = mat * coeff_2;
        for n = 1:length(index_dof_crossing)
            dM_reco_2{muscle_index, order}(:, n) = (-squeeze(diff_mat_q(:, :, n))) * coeff_2;
        end

        % Calculate RMSE_lMT and RMSE_dM with the updated selected coefficients
        Y_reco_2 = [lMT_reco_2{muscle_index, order}; dM_reco_2{muscle_index, order}(:)];
        RMSE_lMT = sqrt(mean((Y(1:10000) - Y_reco_2(1:10000)).^2));
        RMSE_dM = zeros(1, length(dofs));
        for dof = 1:length(dofs)
            start_index = dof * 10000 + 1;
            end_index = (dof +1)* 10000;
            RMSE_dM(dof) = sqrt(mean((Y(start_index:end_index) - Y_reco_2(start_index:end_index)).^2));
        end

        % Debug: Print intermediate results
        disp(['Iteration: ', num2str(i)]);
        disp(['Selected Coefficients: ', num2str(coeff_2')]);
        disp(['RMSE_lMT: ', num2str(RMSE_lMT)]);
        disp(['RMSE_dM: ', num2str(RMSE_dM)]);

        % Define the threshold for each dof
        %the first corresponds to the threshold for lMT and the others for 
        % moment arms
        threshold = ones(1, length(dofs) + 1) * threshold_r; 
            
        for k = 1:length(dofs)
            if any(strcmp(dofs{k}, {'knee_tx', 'knee_ty', 'knee_tz'}))
                threshold(k + 1) = threshold_t;
            end
        end

        % Check if the RMSE values satisfy the threshold conditions
        if RMSE_lMT <= threshold_r && all(RMSE_dM <= threshold(2:end))
            % If both conditions are satisfied, count this muscle
            muscles_satisfy = muscles_satisfy + 1; %SHOULD WE REMOVE THIS?
            muscles_with_low_error = [muscles_with_low_error, m_nr];  % Optional: store muscle index
        end

        % Check if the error is below the threshold
        if all([RMSE_lMT, RMSE_dM] <= threshold)
            found = true; % Stop adding coefficients if the error is below the threshold
        else
            i = i + 1; % Increment i if the error is not below the threshold
        end
    end

    % Compute Y_hat
    y_hat = X * coeff_2;
    SSE = sum((Y - y_hat).^2);

    SE = sqrt(diag(SSE / (length(Y) - length(coeff)) * jacobian_1));

    % Calculate the t-values for each coefficient
    t_values = coeff ./ SE;

    % Calculate the p-values for each coefficient using the t-distribution
    df = length(X) - length(coeff);
    p_values = 2 * (1 - tcdf(abs(t_values), df));


    % Get the indices of the selected coefficients
    selected_indices = find(coeff_2 ~= 0);

    coeff_p_greater_05 = sum(p_values(selected_indices) > 0.05);
    
    % Display how many selected coefficients have a p-value > 0.05
    disp(['Muscle ', num2str(muscle_index), ': ', num2str(coeff_p_greater_05), ' coefficients have P-value > 0.05']);
    % Count the number of selected coefficients
    nSelected_coeff = length(selected_indices);
    % Display the number of selected coefficients with p-value > 0.05
    disp(['Muscle ', muscle_name, ' needs ', num2str(nSelected_coeff), ' coefficients to satify the threshold accuracy']);

    % Display the standard errors
    disp(['Muscle ', num2str(muscle_index), ' Order ', num2str(order), ' Standard Errors:']);
    disp(num2str(SE(selected_indices)'));

    % Display the p-values
    disp(['Muscle ', num2str(muscle_index), ' Order ', num2str(order), ' P-values:']);
    disp(num2str(p_values(selected_indices)'));

    % Display the RMSE of lMT and dM for each degree of freedom
    disp('RMSE_lMT:');
    disp(num2str(RMSE_lMT));
    disp('RMSE_dM:');
    disp(num2str(RMSE_dM));
    disp('lMT_negative_value:');
    disp(num2str(has_negative_value'));
    disp('Number of selected coefficients:');
    disp(num2str(length(selected_indices)));
    disp('Total Number of coefficients:');
    disp(num2str(length(coeff_2)'));

    MuscleInfo.muscle(muscle_index).coeff{order} = coeff_2;
    percentage_selected = nSelected_coeff / Total_coefficients * 100;
    muscle_info_summary.nSelected_coeff = nSelected_coeff;
    muscle_info_summary.Coeff_percentage = percentage_selected;

    %% Analyse the terms of the reduced polynomial

    terms=Identify_string_select_coefficients(coeff_2,order,dofs);
    disp(terms);

    % -------------------------------
    % Analyze selected coefficients by category
    % -------------------------------

    total_terms = length(terms);
    constant_count = 0;
    linear_count = 0;
    quadratic_only_count = 0;
    quadratic_mixed_count = 0;
    cubic_only_count = 0;
    cubic_mixed_count = 0;
    biquadratic_only_count = 0;
    biquadratic_mixed_count = 0;
    quintic_count=0;
    sextic_count=0;

    % Ensure terms are valid strings
    if ~iscell(terms)
        terms = cellstr(terms); % Convert to cell array of strings if needed
    end
    
    % Loop over terms to analyze categories
    for term_idx = 1:length(terms)
        term = terms{term_idx};
    
        % Check if term is valid
        if ~ischar(term) && ~isstring(term)
            warning(['Invalid term format at index ', num2str(term_idx)]);
            continue; % Skip this term if it's not a valid string or char
        end
    
        % Ensure term is a character vector
        term = char(term);
    
        % Split the term into components
        tokens = strsplit(term, '*'); % Split the term by '*'
        powers = zeros(1, length(tokens)); % Initialize powers array
    
        % Extract powers from tokens
        for token_idx = 1:length(tokens)
            token = strtrim(tokens{token_idx}); % Clean up whitespace
            if contains(token, '^') % Ensure token has a power
                parts = strsplit(token, '^');
                if length(parts) == 2
                    powers(token_idx) = str2double(parts{2}); % Extract the power
                end
            end
        end
    
        % Classify terms based on power sums and individual power distributions
        power_sum = sum(powers);
        if power_sum == 0
            constant_count = constant_count + 1;
        elseif power_sum == 1 && nnz(powers == 1) == 1
            linear_count = linear_count + 1;
        elseif power_sum == 2
            if nnz(powers == 2) == 1 && nnz(powers == 1) == 0
                quadratic_only_count = quadratic_only_count + 1;
            elseif nnz(powers == 1) == 2
                quadratic_mixed_count = quadratic_mixed_count + 1;
            end
        elseif power_sum == 3
            if nnz(powers == 3) == 1 && nnz(powers > 0) == 1
                cubic_only_count = cubic_only_count + 1;
            else
                cubic_mixed_count = cubic_mixed_count + 1;
            end
            elseif power_sum == 4
            if nnz(powers== 4) == 1
               biquadratic_only_count = biquadratic_only_count + 1; % biquadratic_only only term
            else
               biquadratic_mixed_count = biquadratic_mixed_count + 1; % Cubic mixed term  
            end
            
            %to be updated for polynomials with higher order than 4...
        end
        
    end

    quadratic_count=quadratic_only_count+quadratic_mixed_count;
    cubic_count=cubic_only_count+cubic_mixed_count;
    biquadratic_count=biquadratic_only_count+biquadratic_mixed_count;

    % Update global counters for all muscles
    total_constant_count = total_constant_count + constant_count;
    total_linear_count = total_linear_count + linear_count;
    total_quadratic_single_count = total_quadratic_single_count + quadratic_only_count;
    total_quadratic_interaction_count = total_quadratic_interaction_count + quadratic_mixed_count;
    total_cubic_single_count = total_cubic_single_count + cubic_only_count;
    total_cubic_interaction_count = total_cubic_interaction_count + cubic_mixed_count;
    total_biquadratic_single_count = total_biquadratic_single_count + biquadratic_only_count;
    total_biquadratic_interaction_count = total_biquadratic_interaction_count + biquadratic_mixed_count;
    
    % Update total coefficients count
    total_coefficients_count = total_coefficients_count + Total_coefficients; 
    % Total coefficients for the current muscle (full polynomial)

    % Calculate percentages
    constant_percentage_ar = (constant_count / total_terms) * 100;
    linear_percentage_ar = (linear_count / total_terms) * 100;
    quadratic_only_percentage = (quadratic_only_count / total_terms) * 100;
    quadratic_mixed_percentage = (quadratic_mixed_count / total_terms) * 100;
    cubic_only_percentage = (cubic_only_count / total_terms) * 100;
    cubic_mixed_percentage = (cubic_mixed_count / total_terms) * 100;
    biquadratic_only_percentage = (biquadratic_only_count / total_terms) * 100;
    biquadratic_mixed_percentage = (biquadratic_mixed_count / total_terms) * 100;

    % -------------------------------
    % Dropped coefficients
    % -------------------------------
    % Generate all terms (ensure this includes all coefficients)
    dropped_terms=Identify_string_dropped_coefficients(coeff_2,order,dofs);

    % Analyze dropped terms
    % Initialize counters
    dropped_constant_count = 0;
    dropped_linear_count = 0;
    dropped_quadratic_only_count = 0;
    dropped_quadratic_mixed_count = 0;
    dropped_cubic_only_count = 0;
    dropped_cubic_mixed_count = 0;
    dropped_biquadratic_only_count = 0;
    dropped_biquadratic_mixed_count = 0;

    
    % Analyze each dropped term
    for i = 1:length(dropped_terms)
        term = dropped_terms{i};
        tokens = strsplit(term, ' * '); % Split term into individual components
        powers = zeros(1, length(tokens));
        
        % Extract powers from tokens
        for j = 1:length(tokens)
            if contains(tokens{j}, '^')
                power = sscanf(tokens{j}, 'q%d^%d');
                powers(j) = power(2); % Extract the power
            end
        end
        
        % Classify the term
        power_sum = sum(powers);
        if power_sum == 0
            dropped_constant_count = dropped_constant_count + 1; % Constant term
        elseif power_sum == 1
            dropped_linear_count = dropped_linear_count + 1; % Linear term
        elseif power_sum == 2
            if nnz(powers == 2) == 1
                dropped_quadratic_only_count = dropped_quadratic_only_count + 1; % Squared only term
            else
                dropped_quadratic_mixed_count = dropped_quadratic_mixed_count + 1; % Squared mixed term
            end
        elseif power_sum == 3
            if nnz(powers == 3) == 1
                dropped_cubic_only_count = dropped_cubic_only_count + 1; % Cubic only term
            else
                dropped_cubic_mixed_count = dropped_cubic_mixed_count + 1; % Cubic mixed term
            end
        elseif power_sum == 4
            if nnz(powers== 4) == 1
               dropped_biquadratic_only_count = dropped_biquadratic_only_count + 1; % biquadratic_only only term
            else
              dropped_biquadratic_mixed_count = dropped_biquadratic_mixed_count + 1; % Cubic mixed term  
            end
        end
    end
    
    % Accumulate dropped term counts for all muscles
    total_dropped_constant_count = total_dropped_constant_count + dropped_constant_count;
    total_dropped_linear_count = total_dropped_linear_count + dropped_linear_count;
    total_dropped_quadratic_only_count = total_dropped_quadratic_only_count + dropped_quadratic_only_count;
    total_dropped_quadratic_mixed_count = total_dropped_quadratic_mixed_count + dropped_quadratic_mixed_count;
    total_dropped_cubic_only_count = total_dropped_cubic_only_count + dropped_cubic_only_count;
    total_dropped_cubic_mixed_count = total_dropped_cubic_mixed_count + dropped_cubic_mixed_count;
    total_dropped_biquadratic_only_count = total_dropped_biquadratic_only_count + dropped_biquadratic_only_count;
    total_dropped_biquadratic_mixed_count = total_dropped_biquadratic_mixed_count + dropped_biquadratic_mixed_count;

    total_dropped = total_dropped_constant_count +total_dropped_linear_count + total_dropped_quadratic_only_count + total_dropped_quadratic_mixed_count + total_dropped_cubic_only_count + total_dropped_cubic_mixed_count +total_dropped_biquadratic_only_count + total_dropped_biquadratic_mixed_count;

    % -------------------------------
    % Display results
    % -------------------------------
    fprintf('Muscle: %s\n', muscle_name);
    fprintf('Constant Terms: %.2f%%\n', constant_percentage_ar);
    fprintf('Linear Terms: %.2f%%\n', linear_percentage_ar);
    fprintf('Only Squared Terms: %.2f%%\n', quadratic_only_percentage);
    fprintf('Mixed Squared Terms: %.2f%%\n', quadratic_mixed_percentage);
    fprintf('Only Cubic Terms: %.2f%%\n', cubic_only_percentage);
    fprintf('Mixed Cubic Terms: %.2f%%\n', cubic_mixed_percentage);
    fprintf('Only Biquadratic Terms: %.2f%%\n', biquadratic_only_percentage);
    fprintf('Mixed Biquadratic Terms: %.2f%%\n', biquadratic_mixed_percentage);
   
    % Calculate the total count of terms
    total_count = constant_count + linear_count + quadratic_count + cubic_count + biquadratic_count + quintic_count + sextic_count;

    

    % Calculate the percentages
    linear_percentage = (linear_count / total_count) * 100;
    squared_percentage = (quadratic_count / total_count) * 100;
    cubic_percentage = (cubic_count / total_count) * 100;
    biquadratic_percentage = (biquadratic_count / total_count) * 100;
    quintic_percentage = (quintic_count / total_count) * 100;
    sextic_percentage = (sextic_count / total_count) * 100;
    constant_percentage = (constant_count / total_count) * 100;
   
    muscle_info_summary.coeff_p_less_05 = length(find(p_values(selected_indices) <= 0.05));
    muscle_info_summary.selected_coeff_p_greater_05 = nSelected_coeff - muscle_info_summary.coeff_p_less_05;
    muscle_info_summary.linear_count = linear_count;
    muscle_info_summary.quadratic_count = quadratic_count;
    muscle_info_summary.cubic_count = cubic_count;
    muscle_info_summary.cubic_only_count=cubic_only_count;
    muscle_info_summary.cubic_mixed_count=cubic_mixed_count;
    muscle_info_summary.biquadratic_count = biquadratic_count;
    muscle_info.summary.biquadratic_only_count= biquadratic_only_count;
    muscle_info.summary.biquadratic_mixed_count=biquadratic_mixed_count;

    muscle_info_summary.linear = linear_percentage;
    muscle_info_summary.squared = squared_percentage;
    muscle_info_summary.cubic = cubic_percentage;
    muscle_info_summary.biquadratic = biquadratic_percentage;
    muscle_info_summary.quintic = quintic_percentage;
    muscle_info_summary.sextic = sextic_percentage;

    muscle_info_summary.RMSE_lMT = RMSE_lMT;
    muscle_info_summary.RMSE_dM = RMSE_dM;

    output_structure.muscle(m_nr) = muscle_info_summary;
    disp(['Muscle: ', muscle_name]);

end

%%
% Calculate final percentages for all muscles together
final_constant_percentage_ar = (total_constant_count / total_coefficients_count) * 100;
final_linear_percentage_ar = (total_linear_count / total_coefficients_count) * 100;
final_quadratic_only_percentage = (total_quadratic_single_count / total_coefficients_count) * 100;
final_quadratic_mixed_percentage = (total_quadratic_interaction_count / total_coefficients_count) * 100;
final_cubic_only_percentage = (total_cubic_single_count / total_coefficients_count) * 100;
final_cubic_mixed_percentage = (total_cubic_interaction_count / total_coefficients_count) * 100;
final_biquadratic_only_percentage = (total_biquadratic_single_count / total_coefficients_count) * 100;
final_biquadratic_mixed_percentage = (total_biquadratic_interaction_count / total_coefficients_count) * 100;

% Display final percentages
fprintf('Final Percentage of constant terms: %.2f%%\n', final_constant_percentage_ar);
fprintf('Final Percentage of linear terms: %.2f%%\n', final_linear_percentage_ar);
fprintf('Final Percentage of squared-only terms: %.2f%%\n', final_quadratic_only_percentage);
fprintf('Final Percentage of squared-mixed terms: %.2f%%\n', final_quadratic_mixed_percentage);
fprintf('Final Percentage of cubic-only terms: %.2f%%\n', final_cubic_only_percentage);
fprintf('Final Percentage of cubic-mixed terms: %.2f%%\n', final_cubic_mixed_percentage);
fprintf('Final Percentage of biquadratic-only terms: %.2f%%\n', final_biquadratic_only_percentage);
fprintf('Final Percentage of biquadratic-mixed terms: %.2f%%\n', final_biquadratic_mixed_percentage);
% Display the total number of muscles that satisfy the conditions
disp(['Total muscles that satisfy the error conditions: ', num2str(muscles_satisfy)]);
%%
% Calculate total dropped terms
total_dropped_terms = total_dropped_constant_count + total_dropped_linear_count + ...
                      total_dropped_quadratic_only_count + total_dropped_quadratic_mixed_count + ...
                      total_dropped_cubic_only_count + total_dropped_cubic_mixed_count + total_dropped_biquadratic_mixed_count + total_dropped_biquadratic_mixed_count;

% Calculate percentages for each category
percentage_constant_dropped = (total_dropped_constant_count / total_dropped_terms) * 100;
percentage_linear_dropped = (total_dropped_linear_count / total_dropped_terms) * 100;
percentage_quadratic_only_dropped = (total_dropped_quadratic_only_count / total_dropped_terms) * 100;
percentage_quadratic_mixed_dropped = (total_dropped_quadratic_mixed_count / total_dropped_terms) * 100;
percentage_cubic_only_dropped = (total_dropped_cubic_only_count / total_dropped_terms) * 100;
percentage_cubic_mixed_dropped = (total_dropped_cubic_mixed_count / total_dropped_terms) * 100;
percentage_biquadratic_only_dropped = (total_dropped_biquadratic_only_count / total_dropped_terms) * 100;
percentage_biquadratic_mixed_dropped = (total_dropped_biquadratic_mixed_count / total_dropped_terms) * 100;

% Display totals and percentages
fprintf('Total Dropped Terms Across All Muscles:\n');
fprintf('Dropped Constant Terms: %d (%.2f%%)\n', total_dropped_constant_count, percentage_constant_dropped);
fprintf('Dropped Linear Terms: %d (%.2f%%)\n', total_dropped_linear_count, percentage_linear_dropped);
fprintf('Dropped Squared-Only Terms: %d (%.2f%%)\n', total_dropped_quadratic_only_count, percentage_quadratic_only_dropped);
fprintf('Dropped Squared-Mixed Terms: %d (%.2f%%)\n', total_dropped_quadratic_mixed_count, percentage_quadratic_mixed_dropped);
fprintf('Dropped Cubic-Only Terms: %d (%.2f%%)\n', total_dropped_cubic_only_count, percentage_cubic_only_dropped);
fprintf('Dropped Cubic-Mixed Terms: %d (%.2f%%)\n', total_dropped_cubic_mixed_count, percentage_cubic_mixed_dropped);
fprintf('Dropped Biquadratic-Only Terms: %d (%.2f%%)\n', total_dropped_biquadratic_only_count, percentage_biquadratic_only_dropped);
fprintf('Dropped Biquadratic-Mixed Terms: %d (%.2f%%)\n', total_dropped_biquadratic_mixed_count, percentage_biquadratic_mixed_dropped);

% Save MuscleInfo with a specific name based on the threshold
if (threshold_r == 0.0005) && (threshold_t==0.005)
    name_MuscleInfo='MuscleInfo_full_05_ropt.mat';
elseif (threshold_r == 0.001) && (threshold_t==0.01)
    name_MuscleInfo='MuscleInfo_full_1_ropt.mat';
elseif (threshold_r == 0.002) && (threshold_t==0.02)
    name_MuscleInfo='MuscleInfo_full_2_ropt.mat';
else
    disp(['these thresholds have not been used for this study, but ' ...
        'feel free to incorporate them']);
    keyboard;
end

% Save the output_structure to a MAT file
save(name_MuscleInfo,'MuscleInfo');

% Initialize sums for each term
sum_linear = sum([output_structure.muscle.Linear]);
sum_squared = sum([output_structure.muscle.squared]);
sum_cubic = sum([output_structure.muscle.cubic]);
sum_biquadratic = sum([output_structure.muscle.biquadratic]);
sum_quintic = sum([output_structure.muscle.quintic]);
sum_sextic = sum([output_structure.muscle.sextic]);

avg_linear = sum_linear / 50;
avg_squared = sum_squared / 50;
avg_cubic = sum_cubic / 50;
avg_biquadratic = sum_biquadratic / 50;
avg_quintic = sum_quintic / 50;
avg_sextic = sum_sextic / 50;
avg_constant = 100 - sum([avg_linear, avg_squared, avg_cubic, avg_biquadratic, avg_quintic, avg_sextic]);

selected_coefficientss = 0;
total_coefficientss = 0;
for y = 1:50
    selected_coefficientss = selected_coefficientss + output_structure.muscle(y).Selected_coeff;
    total_coefficientss = total_coefficientss + output_structure.muscle(y).Total_coeff;
end

toc;  

function terms=Identify_string_select_coefficients(selected_coefficients,order,dofs)

    selected_coefficients_bool = logical(selected_coefficients); % convert to logical array
    terms = cell(1, sum(selected_coefficients_bool)); % initialize the terms cell array %% NEED TO DOUBLE CHECK, THE ORDER MATTERS
    
    coeff_nr = 1;
    term_idx = 1; % keep track of the current term index
    for n_q1 = 0:order
        if length(dofs) < 2
            n_q2s = 0;
        else
            n_q2s = 0:order - n_q1;
        end
        for n_q2 = n_q2s
            if length(dofs) < 3
                n_q3s = 0;
            else
                n_q3s = 0:order - n_q1 - n_q2;
            end
            for n_q3 = n_q3s
                if length(dofs) < 4
                    n_q4s = 0;
                else
                    n_q4s = 0:order - n_q1 - n_q2 - n_q3;
                end
                for n_q4 = n_q4s
                    if length(dofs) < 5
                        n_q5s = 0;
                    else
                        n_q5s = 0:order - n_q1 - n_q2 - n_q3 - n_q4;
                    end
                    for n_q5 = n_q5s
                        if length(dofs) < 6
                            n_q6s = 0;
                        else
                            n_q6s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5;
                        end
                        for n_q6 = n_q6s
                            if length(dofs) < 7
                                n_q7s = 0;
                            else
                                n_q7s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6;
                            end
                            for n_q7 = n_q7s
                                if length(dofs) < 8
                                    n_q8s = 0;
                                else
                                    n_q8s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7;
                                end
                                for n_q8 = n_q8s
                                    if length(dofs) < 9
                                        n_q9s = 0;
                                    else
                                        n_q9s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7 - n_q8;
                                    end
                                    for n_q9 = n_q9s
                                        if selected_coefficients_bool(coeff_nr) % if this coefficient is selected
                                            % store the term as a string in the terms cell array
                                            terms{term_idx} = sprintf('q1^%d * q2^%d * q3^%d * q4^%d * q5^%d * q6^%d * q7^%d * q8^%d * q9^%d', n_q1, n_q2, n_q3, n_q4, n_q5, n_q6, n_q7, n_q8, n_q9);
                                            term_idx = term_idx + 1; % increment the term index
                                        end
                                        coeff_nr = coeff_nr + 1; % increment the coefficient number
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


function dropped_terms=Identify_string_dropped_coefficients(selected_coefficients,order,dofs)

    selected_coefficients_bool = logical(selected_coefficients); % convert to logical array
    terms_all = cell(1, length(selected_coefficients_bool)); % Initialize
    coeff_nr = 1;
    term_idx = 1; % Counter for all terms

    dropped_terms = {}; % To store dropped terms
    selected_terms = {}; % To store selected terms
    for n_q1 = 0:order
        if length(dofs) < 2
            n_q2s = 0;
        else
            n_q2s = 0:order - n_q1;
        end
        for n_q2 = n_q2s
            if length(dofs) < 3
                n_q3s = 0;
            else
                n_q3s = 0:order - n_q1 - n_q2;
            end
            for n_q3 = n_q3s
                if length(dofs) < 4
                    n_q4s = 0;
                else
                    n_q4s = 0:order - n_q1 - n_q2 - n_q3;
                end
                for n_q4 = n_q4s
                    if length(dofs) < 5
                        n_q5s = 0;
                    else
                        n_q5s = 0:order - n_q1 - n_q2 - n_q3 - n_q4;
                    end
                    for n_q5 = n_q5s
                        if length(dofs) < 6
                            n_q6s = 0;
                        else
                            n_q6s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5;
                        end
                        for n_q6 = n_q6s
                            if length(dofs) < 7
                                n_q7s = 0;
                            else
                                n_q7s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6;
                            end
                            for n_q7 = n_q7s
                                if length(dofs) < 8
                                    n_q8s = 0;
                                else
                                    n_q8s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7;
                                end
                                for n_q8 = n_q8s
                                    if length(dofs) < 9
                                        n_q9s = 0;
                                    else
                                        n_q9s = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7 - n_q8;
                                    end
                                    for n_q9 = n_q9s
                                        % Construct the term string
                                    term_str = sprintf('q1^%d * q2^%d * q3^%d * q4^%d * q5^%d * q6^%d * q7^%d * q8^%d * q9^%d', ...
                                        n_q1, n_q2, n_q3, n_q4, n_q5, n_q6, n_q7, n_q8, n_q9);
                                    
                                    % Store in terms_all
                                    terms_all{term_idx} = term_str;
                                    
                                    % Separate dropped terms
                                    if ~selected_coefficients_bool(coeff_nr)
                                        dropped_terms{end+1} = term_str; % Add to dropped terms
                                    end
                                    
                                    term_idx = term_idx + 1; % Increment term index
                                    coeff_nr = coeff_nr + 1; % Increment coefficient counter
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

