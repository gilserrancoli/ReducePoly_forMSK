% This function computes the polynomials to approximate muscle-tendon 
% lengths, velocities and moment arms. 
%
% Authors: Original code from Wouter Aerts, adapted by Antoine Falisse to
% work with AD, and by Gil Serrancolí to work with muscles spanning 9 DoFs
%
% Date: 01/10/2024
%
function [muscle_spanning_joint_INFO,MuscleInfo] =  PolynomialFit(MuscleData)
  
%% Construct the polynomials for the moment arms and muscle length


    muscle_sel=1:50; %all lower-limb and lumbar muscles

    muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO<=0.0001 & ...
        muscle_spanning_joint_INFO>=-0.0001) = 0;
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO~=0) = 1;

    q_all = MuscleData.q;

    max_order = 9; % maximal polynomial order
    threshold_r = 0.002; % maximal deviation from reference OpenSim data (in m) for rotational degrees of freedom
    threshold_t = 0.02; % maximal deviation from reference OpenSim data (in []) for translational degrees of freedom
    nr_samples = length(q_all(:,1));

    lMT_all_error = zeros(length(muscle_sel), 1);
    DM_all_error = zeros(length(muscle_sel), length(q_all(1,:)));
    order_all = zeros(length(muscle_sel), 1);

    for m_nr=1:length(muscle_sel)
        fprintf([MuscleData.muscle_names_forMA{m_nr} '\n']); % Print the current muscle
        muscle_index = muscle_sel(m_nr);

        index_dof_crossing = ...
            find(muscle_spanning_joint_INFO(muscle_index,:)==1);


        nr_dof_crossing = length(index_dof_crossing);


        lMT = MuscleData.lMT(:,muscle_index);
        dM = zeros(nr_samples, nr_dof_crossing);
        for dof_nr = 1:nr_dof_crossing
            dM(:,dof_nr) = ...
                MuscleData.dM(:,muscle_index,index_dof_crossing(dof_nr));
        end

        criterion_full_filled = 0;
        order = 3;
        while criterion_full_filled==0
            fprintf('Order: %d\n', order); % Print the current order
            
            [mat,diff_mat_q] = n_art_mat_9(q_all(:,index_dof_crossing), order);
                     
            nr_coeffs = length(mat(1,:));

            diff_mat_q_all = zeros(nr_samples*nr_dof_crossing, nr_coeffs);
            for dof_nr = 1:nr_dof_crossing
                diff_mat_q_all(nr_samples*(dof_nr-1)+1:nr_samples*dof_nr,:)= ...
                    -squeeze(diff_mat_q(:,:,dof_nr));
            end

            coeff{order}=[mat ; diff_mat_q_all]\[lMT; dM(:)];
            dM_recon = zeros(nr_samples, nr_dof_crossing);         

            for dof_nr = 1:nr_dof_crossing
                dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff{order};
            end
            lMT_recon=mat*coeff{order};

            lMT_error_rms{order} = sqrt(mean((lMT - lMT_recon).^2));
            dm_error_rms{order} = sqrt(mean((dM - dM_recon).^2));
            
            % Print the RMS errors
            fprintf('lMT RMS error: %.4f\n', lMT_error_rms{order});
            fprintf('dM RMS error: %.4f\n', max(dm_error_rms{order}));

            % RMS errors should be lower than threshold for muscle-tendon
            % lengths and moment arms. Otherwise, we increase the
            % polynomial order.
            if any(contains(MuscleData.dof_names(index_dof_crossing),{'knee_tx','knee_ty','knee_tz'}))
                t_pos=find(contains(MuscleData.dof_names(index_dof_crossing),{'knee_tx','knee_ty','knee_tz'}));
                r_pos=1:length(index_dof_crossing);
                r_pos(t_pos)=[];
                criterion_full_filled = ...
                    lMT_error_rms{order}<=threshold_r &...
                    max(dm_error_rms{order}(r_pos))<=threshold_r &...
                    max(dm_error_rms{order}(t_pos))<=threshold_t;
            else
                criterion_full_filled = ...
                    lMT_error_rms{order}<=threshold_r & max(dm_error_rms{order})<=threshold_r;
            end
            if order==max_order
                criterion_full_filled = 1;
            end
            if criterion_full_filled==0
                order = order+1;
            end
            fprintf('\n');
        end

        MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
        MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names_forMA{muscle_index};
        MuscleInfo.muscle(m_nr).coeff = coeff;
        MuscleInfo.muscle(m_nr).order = order;
        MuscleInfo.muscle(m_nr).lMT_error_rms = lMT_error_rms;
        MuscleInfo.muscle(m_nr).dm_error_rms = dm_error_rms;

        lMT_all_error(m_nr) = lMT_error_rms{order};
        DM_all_error(m_nr, index_dof_crossing) = dm_error_rms{order};
        order_all(m_nr) = order;      

        clear coeff lMT_error_rms dm_error_rms;

    end

    figure();
    hold on;
    plot(lMT_all_error)
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold_r, threshold_r], 'r', 'linewidth', 2)
    title('RMS error on the approximated muscle-tendon length')
    ylabel('RMS error (m)')

    figure();
    hold on;
    plot(max(DM_all_error, [], 2))
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold_t, threshold_t], 'r', 'linewidth', 2); %though this is the threshold for translational dofs
    title('maximal RMS error on the approximated muscle moment arm')
    ylabel('RMS error [m] / []'); 

    figure();
    hold on;
    plot(order_all)
    ylim([0 max_order+1])
    xlimits = get(gca, 'XLim');
    plot(xlimits, [max_order, max_order], 'r', 'linewidth', 2)
    title('Order of the polynomial approximation')
    ylabel('Order')

nematicend
