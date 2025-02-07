% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts, adapted to be used with 
% CasADi from Antoine Falisse, and adapted to use muscles spanning joints 
% with up to 9 DoFs by Gil Serrancolí.
%
% Date: 01/10/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% User inputs
runPolynomialfit = 1;
saveQdot = 0;
savePolynomials = 1;

%% Extract time and angles from dummy motion

pathmain = pwd;
name_dummymotion = 'random_motion_10s.mot';
path_dummymotion = [pathmain,'\MuscleAnalysis\dummy_motion\'];
path_resultsMA = [pathmain,'\MuscleAnalysis\resultsMA\'];

dummy_motion = importdata([path_dummymotion,name_dummymotion]);

% Define the desired order of column names
new_order = {'time', 'hip_flexion', 'hip_adduction', 'hip_rotation', ...
             'knee_adduction', 'knee_rotation', 'knee_flexion', ...
             'knee_tx', 'knee_ty', 'knee_tz', 'ankle_angle', ...
             'subtalar_angle', 'lumbar_extension', 'lumbar_bending', ...
             'lumbar_rotation', 'mtp_angle'};

% Reorder the columns and their corresponding values
[~, old_order] = ismember(new_order, dummy_motion.colheaders);
dummy_motion.data = dummy_motion.data(:, old_order);
dummy_motion.colheaders = dummy_motion.colheaders(old_order);

% % Display the updated column names
% dummy_motion.colheaders

q = dummy_motion.data(:,2:end);
q(:, 1:6) = q(:, 1:6) * pi / 180;
q(:, 10:15) = q(:, 10:15) * pi / 180;

%% Import data
% lMT
lMT = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_Length.sto']);
% hip flexion r
MA.hip.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_flexion.sto']);
% hip adduction r
MA.hip.add = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_adduction.sto']);
% hip rotation r
MA.hip.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_rotation.sto']);
% knee adduction r 
MA.knee.add = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_adduction.sto']);
% knee rotation r 
MA.knee.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_rotation.sto']);
% knee flexion r 
MA.knee.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_flexion.sto']);
% knee tx r 
MA.knee.tx = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_tx.sto']);
% knee ty r 
MA.knee.ty = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_ty.sto']);
% knee tz r 
MA.knee.tz = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_tz.sto']);
% ankle flexion r
MA.ankle.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_ankle_angle.sto']);
% subtalar r
MA.sub = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_subtalar_angle.sto']);
% lumbar extension
MA.trunk.ext = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_extension.sto']);
% lumbar bending
MA.trunk.ben = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_bending.sto']);
% lumbar rotation
MA.trunk.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_rotation.sto']);
% % Metatarsophalangeal joint (mtp angle)
MA.mtp = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_mtp_angle.sto']);

%% Organize MuscleData
if runPolynomialfit
    MuscleData.dof_names = dummy_motion.colheaders(2:end); %right leg deofs and lumbar dofs 
    
    %note that the model has 4 muscles ending with _withreallength
    %this was because in order to compute the knee moment arm, the muscle 
    % needed to be inserted at the tibia, however, in order to be compute 
    % the length, the muscle needs to be inserted at the patella

    muscleNames_forLengths = {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem_withreallength','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint_withreallength','vaslat_withreallength',...
        'vasmed_withreallength','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'}; 
    muscleNames_forMA = {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint','vaslat',...
        'vasmed','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'};
    MuscleData.muscle_names_forLengths = muscleNames_forLengths;
    MuscleData.muscle_names_forMA = muscleNames_forMA;
    for m = 1:length(muscleNames_forLengths)
        MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.colheaders,muscleNames_forLengths{m}));       % lMT    
        MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % hip_flex
        MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % hip_add
        MuscleData.dM(:,m,3)    = MA.hip.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % hip_rot
        MuscleData.dM(:,m,4)    = MA.knee.add.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % knee add
        MuscleData.dM(:,m,5)    = MA.knee.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % knee rot
        MuscleData.dM(:,m,6)    = MA.knee.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % knee flex
        MuscleData.dM(:,m,7)    = MA.knee.tx.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee tx
        MuscleData.dM(:,m,8)    = MA.knee.ty.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee ty
        MuscleData.dM(:,m,9)    = MA.knee.tz.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee tz
        MuscleData.dM(:,m,10)    = MA.ankle.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));  % ankle
        MuscleData.dM(:,m,11)    = MA.sub.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));         % sub 
        MuscleData.dM(:,m,12)    = MA.trunk.ext.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk ext
        MuscleData.dM(:,m,13)    = MA.trunk.ben.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk ben
        MuscleData.dM(:,m,14)    = MA.trunk.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk rot
        MuscleData.dM(:,m,15)    = MA.mtp.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % mtp angle
    end
    MuscleData.q = q;
end

%% Call PolynomialFit
if runPolynomialfit
    [muscle_spanning_joint_INFO,MuscleInfo] = PolynomialFit(MuscleData);

    if savePolynomials
        save MuscleData_subject_GC MuscleData
        save muscle_spanning_joint_INFO_subject_GC muscle_spanning_joint_INFO
        save MuscleInfo
    end
end


%% Example of CasADi function to create the full polynomial
%f_lMT_vMT_dM is a function taking lower-limb kinematics as input (qin,
% qdotin) and provides lMT,vMT,dM as output

import casadi.*
% Order mobilities: hip_flex, hip_add, hip_rot, knee_angle, ankle-angle, 
load muscle_spanning_joint_INFO_subject_GC.mat
NMuscle = length(MuscleInfo.muscle);
q_leg_trunk = 14;
qin     = SX.sym('qin',1,q_leg_trunk);
qdotin  = SX.sym('qdotin',1,q_leg_trunk);
lMT     = SX(NMuscle,1);
vMT     = SX(NMuscle,1);
dM      = SX(NMuscle,q_leg_trunk);
for i=1:NMuscle  
    index_dof_crossing = find(muscle_spanning_joint_INFO(i,:)==1);
    order                        = MuscleInfo.muscle(i).order;
    %here we consider that mtp is locked (at 0º)
    if  any(index_dof_crossing(index_dof_crossing==15))
        [mat,diff_mat_q]             = n_art_mat_9_cas_SX([qin(1,index_dof_crossing(1:end-1)) 0],order);
    else
        [mat,diff_mat_q]             = n_art_mat_9_cas_SX(qin(1,index_dof_crossing),order);
    end
    lMT(i,1)                     = mat*MuscleInfo.muscle(i).coeff{order};
    vMT(i,1)                     = 0;
    dM(i,1:q_leg_trunk)          = 0;
    nr_dof_crossing              = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        if index_dof_crossing(dof_nr)==15
        else
            dM(i,index_dof_crossing(dof_nr)) = (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff{order};
            vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*qdotin(1,index_dof_crossing(dof_nr)));
        end
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

%% Check results
load MuscleData_subject_GC.mat
lMT_out_r = zeros(size(q,1),NMuscle);
vMT_out_r = zeros(size(q,1),NMuscle);
dM_out_r = zeros(size(q,1),NMuscle,q_leg_trunk);
for i = 1:size(q,1)
    MuscleData.qdot=zeros(size(MuscleData.q)); %since only lMT and dM are 
    % evaluated here
    [out1_r,out2_r,out3_r] = f_lMT_vMT_dM(MuscleData.q(i,1:end-1),MuscleData.qdot(i,1:end-1));
    lMT_out_r(i,:) = full(out1_r);
    vMT_out_r(i,:) = full(out2_r);
    for j=1:size(dM_out_r,3)
        dM_out_r(i,:,j) = full(out3_r(:,j));
    end
end
 