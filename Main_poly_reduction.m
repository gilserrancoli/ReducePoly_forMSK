% %%
clear all;
close all;
clc;
tic;
current_folder = pwd;
load MuscleData_subject_GC.mat

% Define the threshold in mm
threshold_mm = 2; % Change this value to 0.5, 1, or 2 mm as required

% Load the appropriate file based on the threshold
if threshold_mm == 0.5
    load('MuscleInfo_full_05.mat');
elseif threshold_mm == 1
    load('MuscleInfo_full_1.mat');
elseif threshold_mm == 2
    load('MuscleInfo_full_2.mat');
else
    error('Invalid threshold value. Please select 0.5, 1, or 2 mm.');
end

muscle_sel = 1:50;



muscles_with_low_error = [];
muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO <= 0.0001 & muscle_spanning_joint_INFO >= -0.0001) = 0;
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO ~= 0) = 1;
muscles_satisfy = 0;
muscles_not_satisfy = 0;

se_all_muscles = cell(length(muscle_sel), 1);
% Initialize counters for all muscles
total_constant_count = 0;
total_linear_count = 0;
total_squared_only_count = 0;
total_squared_mixed_count = 0;
total_cubic_only_count = 0;
total_cubic_mixed_count = 0;
total_biquadratic_only_count = 0;
total_biquadratic_mixed_count = 0;

% Counter for total coefficients across all muscles
total_coefficients_count = 0;

%%
% Initialize totals for dropped terms across all muscles
total_dropped_constant_count = 0;
total_dropped_linear_count = 0;
total_dropped_squared_only_count = 0;
total_dropped_squared_mixed_count = 0;
total_dropped_cubic_only_count = 0;
total_dropped_cubic_mixed_count = 0;
total_dropped_biquadratic_only_count = 0;
total_dropped_biquadratic_mixed_count = 0;
%%
for m_nr = muscle_sel
    muscle_index = muscle_sel(m_nr);
    % muscle_index = muscle_sel;


    musle_info_summary.muscle_number = muscle_index;
    orders = MuscleInfo.muscle(muscle_index).order;
    Total_coefficients = length(MuscleInfo.muscle(muscle_index).coeff{orders});
    dofs = MuscleInfo.muscle(muscle_index).DOF;
    muscle_name = MuscleInfo.muscle(muscle_index).m_name;
    muscle_info_summary.DoF = dofs;
    muscle_info_summary.M_name = muscle_name;
    muscle_info_summary.Total_coeff = Total_coefficients;

    order = orders;
    index_dof_crossing_withknee = find(muscle_spanning_joint_INFO(muscle_index, :) == 1);
    [mat, diff_mat_q] = n_art_mat_9_GC(MuscleData.q(:, index_dof_crossing_withknee), order);
    coeff = MuscleInfo.muscle(muscle_index).coeff{order};

    lMT_reco{muscle_index, order} = mat * MuscleInfo.muscle(muscle_index).coeff{1, order};
    lMT_matrix = lMT_reco{muscle_index, order};
    has_negative_value = any(lMT_matrix(:) < 0);

    for n = 1:length(index_dof_crossing_withknee)
        dM_reco{muscle_index, order}(:, n) = (-squeeze(diff_mat_q(:, :, n))) * MuscleInfo.muscle(muscle_index).coeff{1, order};
    end

    X = mat;
    Y = MuscleData.lMT(:, muscle_index);

    for j = 1:size(diff_mat_q, 3)
        X = [X; -diff_mat_q(:, :, j)];
        Y = [Y; reshape(MuscleData.dM(:, muscle_index, index_dof_crossing_withknee(j)), [], 1)];
    end
    Y_reco = [lMT_reco{muscle_index, order}; dM_reco{muscle_index, order}(:)];

    % Calculate the Jacobian matrix
    [m, p] = size(X);  % Number of frames in X values and coefficients
    Jacobian = X;  % Jacobian matrix is the same as X here

    % Calculate the transpose of the Jacobian matrix
    Jacobian_transpose = Jacobian';

    % Compute Y_hat
    y_hat = X * coeff;
    SSE = sum((Y - y_hat).^2);

    jacobian_1 = pinv(Jacobian_transpose * Jacobian);
    SE = sqrt(diag(SSE / (length(Y) - length(coeff)) * jacobian_1));

    % Calculate the t-values for each coefficient
    t_values = coeff ./ SE;

    % Calculate the p-values for each coefficient using the t-distribution
    df = length(X) - length(coeff);
    p_values = 2 * (1 - tcdf(abs(t_values), df));



    threshold_p = 0.05;

    selected_indices = find(p_values <= threshold_p);

    % Order the coefficients by their p-value
    [p_values_sorted, index] = sort(p_values);
    coeff_sorted = coeff(index);

    % Initialize the selected coefficients with zeros
    i = 1;
    found = false;
    while (~found) && (i <= length(coeff))
        lMT = MuscleData.lMT(:, muscle_index);
        dM = squeeze(MuscleData.dM(:, m_nr, find(muscle_spanning_joint_INFO(m_nr, :))));

        % Update selected coefficients
        selected_index = index(1:i);

        % Create selected_mat and selected_diff_mat_q_all using only the selected coefficients
        % selected_indices = find(selected_coefficients ~= 0);

        npoints = size(diff_mat_q, 1);
        ndofs = length(dofs);  % Number of degrees of freedom
        diff_mat_q_all = zeros(npoints * ndofs, size(diff_mat_q, 2));

        for d = 1:ndofs
            start_idx = npoints * (d - 1) + 1;
            end_idx = npoints * d;
            diff_mat_q_all(start_idx:end_idx, :) = -diff_mat_q(:, :, d);
        end

        coeff_aux = [mat(:, selected_index); diff_mat_q_all(:, selected_index)] \ [lMT; dM(:)];

        % Update the full coefficient vector with the selected ones
        coeff_2 = zeros(size(mat, 2), 1);
        coeff_2(selected_index) = coeff_aux;

        % Recalculate lMT and dM with the updated selected coefficients
        lMT_reco_2{muscle_index, order} = mat * coeff_2;
        for n = 1:length(index_dof_crossing_withknee)
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

        % Define the threshold array
        if threshold_mm == 0.5
            threshold_value = 0.0005;
            knee_threshold = 0.005;
        elseif threshold_mm == 1
            threshold_value = 0.001;
            knee_threshold = 0.01;
        elseif threshold_mm == 2
            threshold_value = 0.002;
            knee_threshold = 0.02;
        end

        % Define the threshold array
        threshold = ones(1, length(dofs) + 1) * threshold_value;
        for k = 1:length(dofs)
            if any(strcmp(dofs{k}, {'knee_tx', 'knee_ty', 'knee_tz'}))
                threshold(k + 1) = knee_threshold;
            end
        end


            % After calculating RMSE_lMT and RMSE_dM
        % Check if the RMSE values satisfy the threshold conditions
        if RMSE_lMT <= threshold_value && all(RMSE_dM <= threshold(2:end))
            % If both conditions are satisfied, count this muscle
            muscles_satisfy = muscles_satisfy + 1;
            muscles_with_low_error = [muscles_with_low_error, m_nr];  % Optional: store muscle index
        end



        % Check if the error is below the threshold
        if all([RMSE_lMT, RMSE_dM] <= threshold)
            found = true; % Stop adding coefficients if the error is below the threshold
        else
            i = i + 1; % Increment i if the error is not below the threshold
        end
    end

    %Recompute p-values

    Jacobian = X;  % Jacobian matrix is the same as X here

    % Calculate the transpose of the Jacobian matrix
    Jacobian_transpose = Jacobian';

    % Compute Y_hat
    y_hat = X * coeff;
    SSE = sum((Y - y_hat).^2);

    jacobian_1 = pinv(Jacobian_transpose * Jacobian);
    SE = sqrt(diag(SSE / (length(Y) - length(coeff)) * jacobian_1));

    % Calculate the t-values for each coefficient
    t_values = coeff ./ SE;

    % Calculate the p-values for each coefficient using the t-distribution
    df = length(X) - length(coeff);
    p_values = 2 * (1 - tcdf(abs(t_values), df));


    % Get the indices of the selected coefficients
    selected_indices = find(coeff_2 ~= 0);

    coeff_p_greater_05 = sum(p_values > 0.05);
    
    % Display how many selected coefficients have a p-value > 0.05
    disp(['Muscle ', num2str(muscle_index), ': ', num2str(coeff_p_greater_05), ' coefficients have P-value > 0.05']);
    % Count the number of selected coefficients
    selected_coeff_p_greater_005 = sum(p_values > 0.05);
    % Display the number of selected coefficients with p-value > 0.05
    disp(['Muscle ', muscle_name, ' has ', num2str(selected_coeff_p_greater_005), ' coefficients with p-value > 0.05 out of ', num2str(length(find(MuscleInfo.muscle(muscle_index).coeff{order} ~= 0))), ' selected coefficients.']);


    % Display the standard errors
    disp(['Muscle ', num2str(muscle_index), ' Order ', num2str(order), ' Standard Errors:']);
    disp(num2str(SE'));

    % Display the p-values
    disp(['Muscle ', num2str(muscle_index), ' Order ', num2str(order), ' P-values:']);
    disp(num2str(p_values'));

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
    Selected_coeff = length(find(MuscleInfo.muscle(muscle_index).coeff{order} ~= 0));
    percentage_selected = Selected_coeff / Total_coefficients * 100;
    muscle_info_summary.Selected_coeff = Selected_coeff;
    muscle_info_summary.Coeff_percentage = percentage_selected;

    selected_coefficients=coeff_2;
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
    disp(terms);
    counts = containers.Map('KeyType', 'char', 'ValueType', 'any');
    total_count = 0;
    dof_names_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    poly_order_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:length(terms)
        term = terms{i};
        tokens = split(term, ' * ');

        for j = 1:length(tokens)
            if tokens{j}(end) ~= '0' % Ignore the DOF with power 0
                if counts.isKey(tokens{j})
                    counts(tokens{j}) = counts(tokens{j}) + 1;
                else
                    counts(tokens{j}) = 1;
                    dof_index = str2num(tokens{j}(2:end-2)); % get the index of the DOF
                    dof_names_map(tokens{j}) = dofs{dof_index}; % map the term with its corresponding DOF
                end
            end
        end

        total_count = total_count + 1; % increment total count after processing each term
    end
    % Display the mapping of each selected coefficient to the joint (DoF)
    
    disp('Displaying selected terms with non-zero DoF powers:');
    keys = counts.keys;
        for k = 1:length(keys)
            term = keys{k};
            disp(['Term: ', term, ', DoF: ', dof_names_map(term)]);
        end
    
    % Initialize counters for each DoF
    num_dofs = 9;
    dof_count = zeros(1, num_dofs); % Array to store counts for q1 to q9
    
    for term_idx = 1:length(terms)
        % Get the tokens from the term string to analyze the powers and DoFs
        tokens_1 = strsplit(terms{term_idx}, '*');
    
        % Initialize a string to collect non-zero powered DoF information for this term
        dof_info = '';
        include_term = false;  % Flag to track if this term has any non-zero powers
    
        % Loop over each token to gather the power for each DoF
        for token_idx = 1:length(tokens_1)
            token_1 = strtrim(tokens_1{token_idx});  % Remove any extra whitespace
            % Extract DoF identifier (e.g., q1) and power
            parts = strsplit(token_1, '^');
    
            if length(parts) == 2
                dof_1 = parts{1};  % e.g., 'q1'
                power = str2double(parts{2});  % Convert power to a number
    
                % Only include DoF if power is greater than zero
                if power > 0
                    include_term = true;  % Mark this term to be included
                    % Append to the DoF information string
                    dof_info = [dof_info, dof_1, '^', num2str(power), ' * '];
                    
                    % Update the counter for the corresponding DoF
                    dof_idx = str2double(dof_1(2)); % Extract the index from 'q#'
                    dof_count(dof_idx) = dof_count(dof_idx) + 1;
                end
            end
        end
    
        % Display the term only if it has non-zero powers, removing the trailing ' * ' from dof_info
        if include_term
            dof_info = dof_info(1:end-3);  % Remove the trailing ' * ' for cleaner output
            disp(['Term ', num2str(term_idx), ': ', dof_info]);
        end
    end
    
    % Display the count of each DoF
    for i = 1:num_dofs
        disp(['q', num2str(i), ' count: ', num2str(dof_count(i))]);
    end
    % Initialize the count maps for each power of DOFs
    counts_1 = containers.Map('KeyType', 'char', 'ValueType', 'any');
    counts_2 = containers.Map('KeyType', 'char', 'ValueType', 'any');
    counts_3 = containers.Map('KeyType', 'char', 'ValueType', 'any');

    total_count_1 = 0;  % total count of linear DOFs
    total_count_2 = 0;  % total count of squared DOFs
    total_count_3 = 0;  % total count of cubic DOFs
    total_count_4 = 0;  % total count of biquadratic DOFs
    total_count_5 = 0;  % total count of quintic DOFs
    total_count_6 = 0;  % total count of sextic DOFs

    % Initialize counters for linear, squared, and cubic terms
    linear_count = 0;
    squared_count = 0;
    cubic_count = 0;
    biquadratic_count = 0;
    quintic_count = 0;
    sextic_count = 0;
    constant_count = 0;

    % Loop over the terms
    for i = 1:length(terms)
        term = terms{i};
        tokens = strsplit(term, '*');

        % Initialize counters for q^1, q^2, and q^3
        q1_count = 0;
        q2_count = 0;
        q3_count = 0;
        q4_count = 0;

        % Count the powers
        powers = zeros(1, length(tokens));
        for j = 1:length(tokens)
            split_comp = strsplit(tokens{j}, '^');
            powers(j) = str2double(split_comp{2});
        end
        power_sum = sum(powers);

        % Check if the term is linear, quadratic, or cubic
        if power_sum == 1 && nnz(powers == 1) == 1
            linear_count = linear_count + 1;
        elseif power_sum == 2 && (nnz(powers == 1) == 2 || nnz(powers == 2) == 1)
            squared_count = squared_count + 1;
        elseif power_sum == 3 && (nnz(powers == 1) == 3 || nnz(powers == 2) == 1 && nnz(powers == 1) == 1 || nnz(powers == 3) == 1)
            cubic_count = cubic_count + 1;
        elseif power_sum == 4 && (nnz(powers == 1) == 4 || nnz(powers == 2) == 2 || nnz(powers == 4) == 1 || (nnz(powers == 2) == 1 && nnz(powers == 1) == 2) || nnz(powers == 3) == 1 && nnz(powers == 1) == 1)
            biquadratic_count = biquadratic_count + 1;
        elseif power_sum == 5 && (nnz(powers == 1) == 5 || nnz(powers == 2) == 1 && nnz(powers == 1) == 3 || nnz(powers == 2) == 2 && nnz(powers == 1) == 1 || nnz(powers == 5) == 1 || nnz(powers == 3) == 1 && nnz(powers == 1) == 2 || nnz(powers == 3) == 1 && nnz(powers == 2) == 1)
            quintic_count = quintic_count + 1;
        elseif power_sum == 6 && (nnz(powers == 1) == 6 || nnz(powers == 2) == 3 || nnz(powers == 6) == 1 || nnz(powers == 3) == 2 || (nnz(powers == 4) == 1 && nnz(powers == 2) == 1) || (nnz(powers == 2) == 1 && nnz(powers == 1) == 4) || (nnz(powers == 3) == 1 && nnz(powers == 1) == 3) || (nnz(powers == 4) == 1 && nnz(powers == 1) == 2) || (nnz(powers == 5) == 1 && nnz(powers == 1) == 1))
            sextic_count = sextic_count + 1;
        else
            constant_count = constant_count + 1;
        end
    end
%%
    % -------------------------------
    % Analyze term categories
    % -------------------------------
    total_terms = length(terms);
    constant_count = 0;
    linear_count = 0;
    squared_only_count = 0;
    squared_mixed_count = 0;
    cubic_only_count = 0;
    cubic_mixed_count = 0;
    biquadratic_only_count = 0;
    biquadratic_mixed_count = 0;

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
                squared_only_count = squared_only_count + 1;
            elseif nnz(powers == 1) == 2
                squared_mixed_count = squared_mixed_count + 1;
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
        end
        
    end
    % Update global counters for all muscles
    total_constant_count = total_constant_count + constant_count;
    total_linear_count = total_linear_count + linear_count;
    total_squared_only_count = total_squared_only_count + squared_only_count;
    total_squared_mixed_count = total_squared_mixed_count + squared_mixed_count;
    total_cubic_only_count = total_cubic_only_count + cubic_only_count;
    total_cubic_mixed_count = total_cubic_mixed_count + cubic_mixed_count;
    total_biquadratic_only_count = total_biquadratic_only_count + biquadratic_only_count;
    total_biquadratic_mixed_count = total_biquadratic_mixed_count + biquadratic_mixed_count;
    
    % Update total coefficients count
    total_coefficients_count = total_coefficients_count + Total_coefficients; % Total coefficients for the current muscle

    % Calculate percentages
    constant_percentage_ar = (constant_count / total_terms) * 100;
    linear_percentage_ar = (linear_count / total_terms) * 100;
    squared_only_percentage = (squared_only_count / total_terms) * 100;
    squared_mixed_percentage = (squared_mixed_count / total_terms) * 100;
    cubic_only_percentage = (cubic_only_count / total_terms) * 100;
    cubic_mixed_percentage = (cubic_mixed_count / total_terms) * 100;
    biquadratic_only_percentage = (biquadratic_only_count / total_terms) * 100;
    biquadratic_mixed_percentage = (biquadratic_mixed_count / total_terms) * 100;

%%
    % -------------------------------
    % Dropped coefficients
    % -------------------------------
    % Generate all terms (ensure this includes all coefficients)
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
    % Analyze dropped terms
    % Initialize counters
    dropped_constant_count = 0;
    dropped_linear_count = 0;
    dropped_squared_only_count = 0;
    dropped_squared_mixed_count = 0;
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
                dropped_squared_only_count = dropped_squared_only_count + 1; % Squared only term
            else
                dropped_squared_mixed_count = dropped_squared_mixed_count + 1; % Squared mixed term
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
    total_dropped_squared_only_count = total_dropped_squared_only_count + dropped_squared_only_count;
    total_dropped_squared_mixed_count = total_dropped_squared_mixed_count + dropped_squared_mixed_count;
    total_dropped_cubic_only_count = total_dropped_cubic_only_count + dropped_cubic_only_count;
    total_dropped_cubic_mixed_count = total_dropped_cubic_mixed_count + dropped_cubic_mixed_count;
    total_dropped_biquadratic_only_count = total_dropped_biquadratic_only_count + dropped_biquadratic_only_count;
    total_dropped_biquadratic_mixed_count = total_dropped_biquadratic_mixed_count + dropped_biquadratic_mixed_count;

    total_dropped = total_dropped_constant_count +total_dropped_linear_count + total_dropped_squared_only_count + total_dropped_squared_mixed_count + total_dropped_cubic_only_count + total_dropped_cubic_mixed_count +total_dropped_biquadratic_only_count + total_dropped_biquadratic_mixed_count;
    % % Calculate percentages
    % constant_dropped_percentage = (total_dropped_constant_count / total_dropped) * 100;
    % linear_dropped_percentage_ar = (total_dropped_linear_count / total_dropped) * 100;
    % squared_only_dropped_percentage = (total_dropped_squared_only_count / total_dropped) * 100;
    % squared_mixed_dropped_percentage = (total_dropped_squared_mixed_count / total_dropped) * 100;
    % cubic_only_dropped_percentage = (total_dropped_cubic_only_count / total_dropped) * 100;
    % cubic_mixed_dropped_percentage = (total_dropped_cubic_mixed_count / total_dropped) * 100;

    
%%
    % -------------------------------
    % Display results
    % -------------------------------
    fprintf('Muscle: %s\n', muscle_name);
    fprintf('Constant Terms: %.2f%%\n', constant_percentage_ar);
    fprintf('Linear Terms: %.2f%%\n', linear_percentage_ar);
    fprintf('Only Squared Terms: %.2f%%\n', squared_only_percentage);
    fprintf('Mixed Squared Terms: %.2f%%\n', squared_mixed_percentage);
    fprintf('Only Cubic Terms: %.2f%%\n', cubic_only_percentage);
    fprintf('Mixed Cubic Terms: %.2f%%\n', cubic_mixed_percentage);
    fprintf('Only Biquadratic Terms: %.2f%%\n', biquadratic_only_percentage);
    fprintf('Mixed Biquadratic Terms: %.2f%%\n', biquadratic_mixed_percentage);
    % Display the results
    % fprintf('\nPercentage of Dropped Terms Across All Muscles:\n');
    % fprintf('Constant Dropped Terms: %.2f%%\n', constant_dropped_percentage);
    % fprintf('Linear Dropped Terms: %.2f%%\n', linear_dropped_percentage_ar);
    % fprintf('Squared-Only Dropped Terms: %.2f%%\n', squared_only_dropped_percentage);
    % fprintf('Squared-Mixed Dropped Terms: %.2f%%\n', squared_mixed_dropped_percentage);
    % fprintf('Cubic-Only Dropped Terms: %.2f%%\n', cubic_only_dropped_percentage);
    % fprintf('Cubic-Mixed Dropped Terms: %.2f%%\n', cubic_mixed_dropped_percentage);



    % Calculate the total count of terms
    total_count = constant_count + linear_count + squared_count + cubic_count + biquadratic_count + quintic_count + sextic_count;

    

    % Calculate the percentages
    linear_percentage = (linear_count / total_count) * 100;
    squared_percentage = (squared_count / total_count) * 100;
    cubic_percentage = (cubic_count / total_count) * 100;
    biquadratic_percentage = (biquadratic_count / total_count) * 100;
    quintic_percentage = (quintic_count / total_count) * 100;
    sextic_percentage = (sextic_count / total_count) * 100;
    constant_percentage = (constant_count / total_count) * 100;
    % 
    % % Display the percentages
    % fprintf('Percentage of linear terms: %.2f%%\n', linear_percentage);
    % fprintf('Percentage of squared terms: %.2f%%\n', squared_percentage);
    % fprintf('Percentage of cubic terms: %.2f%%\n', cubic_percentage);
    % fprintf('Percentage of biquadratic terms: %.2f%%\n', biquadratic_percentage);
    % fprintf('Percentage of quintic terms: %.2f%%\n', quintic_percentage);
    % fprintf('Percentage of sextic terms: %.2f%%\n', sextic_percentage);
    % fprintf('Percentage of constant terms: %.2f%%\n', constant_percentage);

    muscle_info_summary.coeff_p_less_05 = length(find(p_values <= 0.05));
    muscle_info_summary.selected_coeff_p_greater_05 = Selected_coeff - muscle_info_summary.coeff_p_less_05;
    muscle_info_summary.selected_coeff_p_greater_005 = coeff_p_greater_05 ;
    muscle_info_summary.linear_count = linear_count;
    muscle_info_summary.squared_count = squared_count;
    muscle_info_summary.cubic_count = cubic_count;
    muscle_info_summary.biquadratic_count = biquadratic_count;

    muscle_info_summary.Linear = linear_percentage;
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
final_squared_only_percentage = (total_squared_only_count / total_coefficients_count) * 100;
final_squared_mixed_percentage = (total_squared_mixed_count / total_coefficients_count) * 100;
final_cubic_only_percentage = (total_cubic_only_count / total_coefficients_count) * 100;
final_cubic_mixed_percentage = (total_cubic_mixed_count / total_coefficients_count) * 100;
final_biquadratic_only_percentage = (total_biquadratic_only_count / total_coefficients_count) * 100;
final_biquadratic_mixed_percentage = (total_biquadratic_mixed_count / total_coefficients_count) * 100;

% Display final percentages
fprintf('Final Percentage of constant terms: %.2f%%\n', final_constant_percentage_ar);
fprintf('Final Percentage of linear terms: %.2f%%\n', final_linear_percentage_ar);
fprintf('Final Percentage of squared-only terms: %.2f%%\n', final_squared_only_percentage);
fprintf('Final Percentage of squared-mixed terms: %.2f%%\n', final_squared_mixed_percentage);
fprintf('Final Percentage of cubic-only terms: %.2f%%\n', final_cubic_only_percentage);
fprintf('Final Percentage of cubic-mixed terms: %.2f%%\n', final_cubic_mixed_percentage);
fprintf('Final Percentage of biquadratic-only terms: %.2f%%\n', final_biquadratic_only_percentage);
fprintf('Final Percentage of biquadratic-mixed terms: %.2f%%\n', final_biquadratic_mixed_percentage);
% Display the total number of muscles that satisfy the conditions
disp(['Total muscles that satisfy the error conditions: ', num2str(muscles_satisfy)]);
%%
% Calculate total dropped terms
total_dropped_terms = total_dropped_constant_count + total_dropped_linear_count + ...
                      total_dropped_squared_only_count + total_dropped_squared_mixed_count + ...
                      total_dropped_cubic_only_count + total_dropped_cubic_mixed_count + total_dropped_biquadratic_mixed_count + total_dropped_biquadratic_mixed_count;

% Calculate percentages for each category
percentage_constant_dropped = (total_dropped_constant_count / total_dropped_terms) * 100;
percentage_linear_dropped = (total_dropped_linear_count / total_dropped_terms) * 100;
percentage_squared_only_dropped = (total_dropped_squared_only_count / total_dropped_terms) * 100;
percentage_squared_mixed_dropped = (total_dropped_squared_mixed_count / total_dropped_terms) * 100;
percentage_cubic_only_dropped = (total_dropped_cubic_only_count / total_dropped_terms) * 100;
percentage_cubic_mixed_dropped = (total_dropped_cubic_mixed_count / total_dropped_terms) * 100;
percentage_biquadratic_only_dropped = (total_dropped_biquadratic_only_count / total_dropped_terms) * 100;
percentage_biquadratic_mixed_dropped = (total_dropped_biquadratic_mixed_count / total_dropped_terms) * 100;

% Display totals and percentages
fprintf('Total Dropped Terms Across All Muscles:\n');
fprintf('Dropped Constant Terms: %d (%.2f%%)\n', total_dropped_constant_count, percentage_constant_dropped);
fprintf('Dropped Linear Terms: %d (%.2f%%)\n', total_dropped_linear_count, percentage_linear_dropped);
fprintf('Dropped Squared-Only Terms: %d (%.2f%%)\n', total_dropped_squared_only_count, percentage_squared_only_dropped);
fprintf('Dropped Squared-Mixed Terms: %d (%.2f%%)\n', total_dropped_squared_mixed_count, percentage_squared_mixed_dropped);
fprintf('Dropped Cubic-Only Terms: %d (%.2f%%)\n', total_dropped_cubic_only_count, percentage_cubic_only_dropped);
fprintf('Dropped Cubic-Mixed Terms: %d (%.2f%%)\n', total_dropped_cubic_mixed_count, percentage_cubic_mixed_dropped);
fprintf('Dropped Biquadratic-Only Terms: %d (%.2f%%)\n', total_dropped_biquadratic_only_count, percentage_biquadratic_only_dropped);
fprintf('Dropped Biquadratic-Mixed Terms: %d (%.2f%%)\n', total_dropped_biquadratic_mixed_count, percentage_biquadratic_mixed_dropped);
%%
% Display the results
% fprintf('\nTotal muscles that satisfy the conditions: %d\n', muscles_satisfy);
% fprintf('Total muscles that do NOT satisfy the conditions: %d\n', muscles_not_satisfy);

% Save MuscleInfo with a specific name based on the threshold
if threshold_mm == 0.5
    % save('MuscleInfo_red_05_ropt_1.mat', 'MuscleInfo');
    % output_filename = sprintf('output_structure_05.mat');

elseif threshold_mm == 1
    % save('MuscleInfo_red_1_ropt_1.mat', 'MuscleInfo');
    % output_filename = sprintf('Output_structure_1.mat');
    % save('output_structure_1.mat', 'output_structure');

elseif threshold_mm == 2
    % save('MuscleInfo_red_2_ropt_1.mat', 'MuscleInfo');
    output_filename = sprintf('output_structure_2.mat');
    % save('output_structure_2.mat', 'output_structure');a

    % output_filename = sprintf('2utput_structure_2.mat');
end

% Save the output_structure to a MAT file
% output_filename = sprintf('Output_structure_1.mat');
% save(output_filename, 'output_structure');

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
%%
%%
% clear all;
% close all;
% clc;
% tic;
% current_folder = pwd;
% load MuscleData_subject_GC.mat
% 
% % Define the threshold in mm
% threshold_mm = 2; % Change this value to 0.5, 1, or 2 mm as required
% 
% % Load the appropriate file based on the threshold
% if threshold_mm == 0.5
%     load('MuscleInfo_full_05.mat');
% elseif threshold_mm == 1
%     load('MuscleInfo_full_1.mat');
% elseif threshold_mm == 2
%     load('MuscleInfo_full_2.mat');
% else
%     error('Invalid threshold value. Please select 0.5, 1, or 2 mm.');
% end
% 
% muscle_sel = 1:50;
% muscles_with_low_error = [];
% muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
% muscle_spanning_joint_INFO(muscle_spanning_joint_INFO <= 0.0001 & muscle_spanning_joint_INFO >= -0.0001) = 0;
% muscle_spanning_joint_INFO(muscle_spanning_joint_INFO ~= 0) = 1;
% muscles_satisfy = 0;
% muscles_not_satisfy = 0;
% muscles_satisfy_p_value = 0;
% 
% results_table = {};
% 
% for m_nr = muscle_sel
%     muscle_index = muscle_sel(m_nr);
%     musle_info_summary.muscle_number = muscle_index;
%     orders = MuscleInfo.muscle(muscle_index).order;
%     Total_coefficients = length(MuscleInfo.muscle(muscle_index).coeff{orders});
%     dofs = MuscleInfo.muscle(muscle_index).DOF;
%     muscle_name = MuscleInfo.muscle(muscle_index).m_name;
%     muscle_info_summary.DoF = dofs;
%     muscle_info_summary.M_name = muscle_name;
%     muscle_info_summary.Total_coeff = Total_coefficients;
% 
%     order = orders;
%     index_dof_crossing_withknee = find(muscle_spanning_joint_INFO(muscle_index, :) == 1);
%     [mat, diff_mat_q] = n_art_mat_9_GC(MuscleData.q(:, index_dof_crossing_withknee), order);
%     coeff = MuscleInfo.muscle(muscle_index).coeff{order};
% 
%     lMT_reco{muscle_index, order} = mat * MuscleInfo.muscle(muscle_index).coeff{1, order};
%     lMT_matrix = lMT_reco{muscle_index, order};
%     has_negative_value = any(lMT_matrix(:) < 0);
% 
%     for n = 1:length(index_dof_crossing_withknee)
%         dM_reco{muscle_index, order}(:, n) = (-squeeze(diff_mat_q(:, :, n))) * MuscleInfo.muscle(muscle_index).coeff{1, order};
%     end
% 
%     X = mat;
%     Y = MuscleData.lMT(:, muscle_index);
% 
%     for j = 1:size(diff_mat_q, 3)
%         X = [X; -diff_mat_q(:, :, j)];
%         Y = [Y; reshape(MuscleData.dM(:, muscle_index, index_dof_crossing_withknee(j)), [], 1)];
%     end
%     Y_reco = [lMT_reco{muscle_index, order}; dM_reco{muscle_index, order}(:)];
% 
%     % Calculate the Jacobian matrix
%     [m, p] = size(X);  % Number of frames in X values and coefficients
%     Jacobian = X;  % Jacobian matrix is the same as X here
% 
%     % Calculate the transpose of the Jacobian matrix
%     Jacobian_transpose = Jacobian';
% 
%     % Compute Y_hat
%     y_hat = X * coeff;
%     SSE = sum((Y - y_hat).^2);
% 
%     jacobian_1 = pinv(Jacobian_transpose * Jacobian);
%     SE = sqrt(diag(SSE / (length(Y) - length(coeff)) * jacobian_1));
% 
%     % Calculate the t-values for each coefficient
%     t_values = coeff ./ SE;
% 
%     % Calculate the p-values for each coefficient using the t-distribution
%     df = length(X) - length(coeff);
%     p_values = 2 * (1 - tcdf(abs(t_values), df));
% 
%     % Filter coefficients with p-values <= 0.05
%     p_values_filtered_indices = find(p_values <= 0.05);
%     num_selected_p_values_lessthan_005 = length(p_values_filtered_indices);
%     muscle_info_summary.num_selected_p_values_lessthan_005 = num_selected_p_values_lessthan_005;
% 
%     % Calculate the number of selected coefficients with p-value > 0.05 from the selected coefficients only
%     selected_coeff_indices = find(coeff ~= 0);
%     selected_p_values_greater_005_indices = intersect(selected_coeff_indices, find(p_values > 0.05));
%     num_selected_p_values_greater_005 = length(selected_p_values_greater_005_indices);
%     muscle_info_summary.num_selected_p_values_greater_005 = num_selected_p_values_greater_005;
% 
%     if isempty(p_values_filtered_indices)
%         continue; % Skip if no coefficients satisfy the p-value condition
%     end
% 
%     % Recompute coefficients using only the filtered indices
%     X_filtered = X(:, p_values_filtered_indices);
%     coeff_filtered = X_filtered \ Y;
% 
%     % Calculate the reconstructed values for lMT and dM with filtered coefficients
%     lMT_reco_filtered = X_filtered * coeff_filtered;
% 
%     % Calculate RMSE for lMT and dM
%     RMSE_lMT_filtered = sqrt(mean((Y(1:10000) - lMT_reco_filtered(1:10000)).^2));
% 
%     RMSE_dM_filtered = zeros(1, length(dofs));
%     for dof = 1:length(dofs)
%         start_index = dof * 10000 + 1;
%         end_index = (dof +1) * 10000;
%         RMSE_dM_filtered(dof) = sqrt(mean((Y(start_index:end_index) - lMT_reco_filtered(start_index:end_index)).^2));
%     end
% 
%     % Define the threshold array
%     if threshold_mm == 0.5
%         threshold_value = 0.0005;
%         knee_threshold = 0.005;
%     elseif threshold_mm == 1
%         threshold_value = 0.001;
%         knee_threshold = 0.01;
%     elseif threshold_mm == 2
%         threshold_value = 0.002;
%         knee_threshold = 0.02;
%     end
% 
%     threshold = ones(1, length(dofs) + 1) * threshold_value;
%     for k = 1:length(dofs)
%         if any(strcmp(dofs{k}, {'knee_tx', 'knee_ty', 'knee_tz'}))
%             threshold(k + 1) = knee_threshold;
%         end
%     end
% 
%     % Check if the error conditions are met
%     if RMSE_lMT_filtered <= threshold_value && all(RMSE_dM_filtered <= threshold(2:end))
%         muscles_satisfy_p_value = muscles_satisfy_p_value + 1;
%     end
% 
%     % Append results to the table
%     results_table = [results_table; {muscle_index, muscle_name, num_selected_p_values_lessthan_005, num_selected_p_values_greater_005}];
% 
%     output_structure.muscle(m_nr) = muscle_info_summary;
% end
% 
% % Convert results table to proper format
% results_table = cell2table(results_table, 'VariableNames', {'Muscle_Index', 'Muscle_Name', 'Num_Selected_p_Values_<=_0.05', 'Num_Selected_p_Values_>_0.05'});
% 
% % Display the results table
% disp(results_table);
% 
% toc;
% cubicpercentage = 0;
% linearpercentage = 0;
% qadpercentage=0;
% for k = 1:50
%     cubicpercentage = cubicpercentage + output_structure.muscle(k).cubic;
%     linearpercentage = linearpercentage +output_structure.muscle(k).Linear;
%     qadpercentage = qadpercentage +output_structure.muscle(k).squared;
% end
% 
% cubicpercentage = cubicpercentage/50;
% linearpercentage = linearpercentage/50;
% qadpercentage = qadpercentage/50;
% disp(['cubicpercentage: ', cubicpercentage]);
% disp(['linearpercentage: ', linearpercentage]);
% disp(['qadpercentage: ', qadpercentage]);
% 
% cubicpercentage_1 = 0;
% linear_1 = 0;
% qadpercentage_1=0;
% for k = 1:50
%     cubicpercentage_1 = cubicpercentage_1 + ouptput_1_old.output_structure.muscle(k).cubic;
%     linear_1 = linear_1 +ouptput_1_old.output_structure.muscle(k).Linear;
%     qadpercentage_1 = qadpercentage_1 +ouptput_1_old.output_structure.muscle(k).squared;
% end
% 
% cubicpercentage_1 = cubicpercentage_1/50;
% linear_1 = linear_1/50;
% qadpercentage_1 = qadpercentage_1/50;
