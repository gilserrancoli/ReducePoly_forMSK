% This script provides an inital guess for the design variables. The guess is
% data-informed (DI). We use experimental data to provide an initial guess of
% the joint variables but set constant values to the muscle and arm variables.
%
% Author: Gil Serrancolí, 
% modified and adapted from Antoine Falisse
% Date: 31/07/2024
% 
function guess = getGuess_3D_DI_GC_v2(Qs,nq,N,NMuscle,jointi,scaling,d,KCF,GRF,F1,F2_in,f_contactForce,Options,mode)

if Options.KCFasinputstoExternalFunction
    F2_skeletal=F2_in{1};
    F2=F2_in{2};
else
    F2=F2_in;
end

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs.allinterpfilt));
Qs_spline.data(:,1) = Qs.allinterpfilt(1:end,1);
Qdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdots_spline.data(:,1) = Qs.allinterpfilt(1:end,1);
Qdotdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdotdots_spline.data(:,1) = Qs.allinterpfilt(1:end,1);
for i = 2:size(Qs.allfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1),Qs.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allinterpfilt(1:end,1),1);
    [Qs_spline_col.data(:,i),Qdots_spline_col.data(:,i),...
        Qdotdots_spline_col.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allinterpfilt_col(1:end,1),1);
end

%% Adapt for the predictive case
if contains(mode,'pred')
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))=Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))-Qs_spline.data(1,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))=Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))-Qs_spline.data(1,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
    Qs_spline_col.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))=Qs_spline_col.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))-Qs_spline_col.data(1,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
    Qs_spline_col.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))=Qs_spline_col.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))-Qs_spline_col.data(1,strcmp(Qs.colheaders(1,:),'pelvis_tz'));

end

% Filter the accelerations
order = 4;
cutoff_low = 10;
fs=1/mean(diff(Qs_spline.data(:,1)));
if cutoff_low./(0.5*fs)>1
    cutoff_low=0.5*fs-1e-3;
end
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qdotdots_spline.data(:,2:end) = filtfilt(af,bf,Qdotdots_spline.data(:,2:end));  

%% Qs: data-informed
% Pelvis tilt
guess.Qs(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rz'));
% Pelvis list
guess.Qs(:,jointi.pelvis.list) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rx'));
% Pelvis rotation
guess.Qs(:,jointi.pelvis.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ry'));
% Pelvis_tx
guess.Qs(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qs(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qs(:,jointi.pelvis.tz) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qs(:,jointi.hip_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qs(:,jointi.hip_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion'));
% Hip adduction
guess.Qs(:,jointi.hip_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qs(:,jointi.hip_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction'));
% Hip rotation
guess.Qs(:,jointi.hip_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qs(:,jointi.hip_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation'));
% Knee adduction
guess.Qs(:,jointi.knee_add.r) = -0.0380;
%Knee rotation
guess.Qs(:,jointi.knee_rot.r) = 0;
% Knee flexion
guess.Qs(:,jointi.knee_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion_l'));
guess.Qs(:,jointi.knee_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion'));
% Knee tx
guess.Qs(:,jointi.knee_tx.r) = 0;
% Knee ty
guess.Qs(:,jointi.knee_ty.r) = 0.042;
% Knee tz
guess.Qs(:,jointi.knee_tz.r) = 0;
% Ankle
guess.Qs(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qs(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle'));
% Subtalar
guess.Qs(:,jointi.subt.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qs(:,jointi.subt.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle'));
% Trunk extension
guess.Qs(:,jointi.trunk.ext) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qs(:,jointi.trunk.ben) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qs(:,jointi.trunk.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qs(:,jointi.sh_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qs(:,jointi.sh_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qs(:,jointi.sh_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qs(:,jointi.sh_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qs(:,jointi.sh_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qs(:,jointi.sh_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qs(:,jointi.elb.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qs(:,jointi.elb.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));
%Knee ty preoptimized
if Options.preoptimizekneedofs==1
    fprintf('Preoptimize knee superior-inferior translation ...\n')
    guess.Qs(:,jointi.knee_ty.r)=OptimizeForkneety(guess.Qs,jointi,KCF.allinterpfilt,N,F2,Options);
    fprintf('knee superior-inferior translation pre-optimized \n')
elseif (Options.preoptimizekneedofs==2)||(Options.preoptimizekneedofs==3)
    fprintf(['Preoptimize knee superior-inferior translation and ' ...
    'abduction-adduction ...\n'])
    guess.Qs(:,[jointi.knee_add.r jointi.knee_ty.r])=OptimizeForkneeaddty(guess.Qs,jointi,KCF.allinterpfilt,N,F2,Options);
    fprintf(['Knee superior-inferior translation and ' ...
    'abduction-adduction pre-optimized \n'])
elseif Options.preoptimizekneedofs==0;
    %with this, we are considering knee_ty=0.042 and knee_add=-0.038;
else
    keyboard;
end



%% Qdots: data-informed
% Pelvis tilt
guess.Qdots(:,jointi.pelvis.tilt) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rz'));
% Pelvis list
guess.Qdots(:,jointi.pelvis.list) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rx'));
% Pelvis rotation
guess.Qdots(:,jointi.pelvis.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ry'));
% Pelvis_tx
guess.Qdots(:,jointi.pelvis.tx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdots(:,jointi.pelvis.ty) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdots(:,jointi.pelvis.tz) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdots(:,jointi.hip_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdots(:,jointi.hip_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion'));
% Hip adduction
guess.Qdots(:,jointi.hip_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdots(:,jointi.hip_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction'));
% Hip rotation
guess.Qdots(:,jointi.hip_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdots(:,jointi.hip_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation'));
%Knee adduction
guess.Qdots(:,jointi.knee_add.r) = 0;
%Knee rotation
guess.Qdots(:,jointi.knee_rot.r) = 0;
% Knee flexion
guess.Qdots(:,jointi.knee_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion_l'));
guess.Qdots(:,jointi.knee_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion'));
%Knee tx
guess.Qdots(:,jointi.knee_tx.r) = 0;
%Knee ty
guess.Qdots(:,jointi.knee_ty.r) = 0;
%Knee tz
guess.Qdots(:,jointi.knee_tz.r) = 0;
% Ankle
guess.Qdots(:,jointi.ankle.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdots(:,jointi.ankle.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle'));
% Subtalar
guess.Qdots(:,jointi.subt.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdots(:,jointi.subt.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle'));
% Trunk extension
guess.Qdots(:,jointi.trunk.ext) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdots(:,jointi.trunk.ben) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdots(:,jointi.trunk.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdots(:,jointi.sh_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdots(:,jointi.sh_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdots(:,jointi.sh_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdots(:,jointi.sh_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdots(:,jointi.sh_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdots(:,jointi.sh_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdots(:,jointi.elb.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdots(:,jointi.elb.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));

%Optimize pelvis tilt and ty, and their velocities
fprintf(['Preoptimize pelvis tilt and ty and their velocities \n']);
if Options.KCFasinputstoExternalFunction
    [pelvis_rtilt,pelvis_ty,pelvis_rtilt_v,pelvis_ty_v]=OptimizeForPelvis_ty_rtilt(guess.Qs,guess.Qdots,Qs.colheaders,jointi,GRF,N,F1,F2_skeletal,f_contactForce,Options); % pre optimization of pelvis_ty and pelvis_tilt does not improve the convergence
else
    [pelvis_rtilt,pelvis_ty,pelvis_rtilt_v,pelvis_ty_v]=OptimizeForPelvis_ty_rtilt(guess.Qs,guess.Qdots,Qs.colheaders,jointi,GRF,N,F1,F2,f_contactForce,Options); % pre optimization of pelvis_ty and pelvis_tilt does not improve the convergence
end
fprintf(['pelvis tilt and ty and their velocities preoptimized\n']);
guess.Qs(:,jointi.pelvis.tilt)=pelvis_rtilt;
guess.Qs(:,jointi.pelvis.ty)=pelvis_ty; 
guess.Qdots(:,jointi.pelvis.tilt)=pelvis_rtilt_v;
guess.Qdots(:,jointi.pelvis.ty)=pelvis_ty_v; 

%% Qs and Qdots are intertwined
guess.QsQdots = zeros(N+1,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots: data-informed
% Pelvis tilt
guess.Qdotdots(:,jointi.pelvis.tilt) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rz'));
% Pelvis list
guess.Qdotdots(:,jointi.pelvis.list) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ry'));
% Pelvis rotation
guess.Qdotdots(:,jointi.pelvis.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rx'));
% Pelvis_tx
guess.Qdotdots(:,jointi.pelvis.tx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdotdots(:,jointi.pelvis.ty) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdotdots(:,jointi.pelvis.tz) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdotdots(:,jointi.hip_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdotdots(:,jointi.hip_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion'));
% Hip adduction
guess.Qdotdots(:,jointi.hip_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdotdots(:,jointi.hip_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction'));
% Hip rotation
guess.Qdotdots(:,jointi.hip_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdotdots(:,jointi.hip_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation'));
% Knee adduction
guess.Qdotdots(:,jointi.knee_add.r) = 0;
% Knee rotation
guess.Qdotdots(:,jointi.knee_rot.r) = 0;
% Knee flexion
guess.Qdotdots(:,jointi.knee_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion_l'));
guess.Qdotdots(:,jointi.knee_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flexion'));
% Knee tx
guess.Qdotdots(:,jointi.knee_tx.r) = 0;
% Knee ty
guess.Qdotdots(:,jointi.knee_ty.r) = 0;
% Knee tz
guess.Qdotdots(:,jointi.knee_tz.r) = 0;
% Ankle
guess.Qdotdots(:,jointi.ankle.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdotdots(:,jointi.ankle.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle'));
% Subtalar
guess.Qdotdots(:,jointi.subt.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdotdots(:,jointi.subt.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle'));
% Trunk extension
guess.Qdotdots(:,jointi.trunk.ext) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdotdots(:,jointi.trunk.ben) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdotdots(:,jointi.trunk.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdotdots(:,jointi.sh_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdotdots(:,jointi.sh_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdotdots(:,jointi.sh_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdotdots(:,jointi.sh_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdotdots(:,jointi.sh_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdotdots(:,jointi.sh_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdotdots(:,jointi.elb.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdotdots(:,jointi.elb.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));

%% Muscle variables
guess.a = 0.1*ones(N+1,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N+1,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm activations
guess.a_a = 0.1*ones(N,nq.arms);
guess.e_a = 0.1*ones(N,nq.arms);

%% Parameters contact model
% Original values
B_locSphere_s1_r    = [0.00190115788407966, -0.00382630379623308];
B_locSphere_s2_r    = [0.148386399942063, -0.028713422052654];
B_locSphere_s3_r    = [0.133001170607051, 0.0516362473449566];
B_locSphere_s4_r    = [0.06, -0.0187603084619177];    
B_locSphere_s5_r    = [0.0662346661991635, 0.0263641606741698];
B_locSphere_s6_r    = [0.045, 0.0618569567549652];
IG_rad              = 0.032*ones(1,6); 
guess.params = [B_locSphere_s1_r,B_locSphere_s2_r,...
    B_locSphere_s3_r,B_locSphere_s4_r,B_locSphere_s5_r,B_locSphere_s6_r,...
    IG_rad];

%% KCF residuals
guess.KCFresiduals=zeros(1,5)./scaling.KCF_res;

%% Joint moment residuals
guess.JRM_res=zeros(1,21)./scaling.JRM_res;

%% Scaling
%all except ty
guess.QsQdots(:,[1:jointi.knee_tx.r jointi.knee_tz.r:nq.all]*2-1)...
    = guess.QsQdots(:,[1:jointi.knee_tx.r jointi.knee_tz.r:nq.all]*2-1)./repmat(scaling.QsQdots([1:jointi.knee_tx.r jointi.knee_tz.r:nq.all]*2-1),N+1,1);
guess.QsQdots(:,jointi.knee_ty.r*2-1) = (guess.QsQdots(:,jointi.knee_ty.r*2-1) - scaling.knee_ty.a)/scaling.knee_ty.b;
guess.QsQdots(:,2:2:end)=guess.QsQdots(:,2:2:end)./repmat(scaling.QsQdots(2:2:end),N+1,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N+1,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));
% no need to scale the IG for the arm activations / excitations
guess.params        = guess.params.*scaling.params.v + scaling.params.r;

%% Collocation points
guess.a_col = zeros(d*N,NMuscle);
guess.FTtilde_col = zeros(d*N,NMuscle);
guess.QsQdots_col = zeros(d*N,2*nq.all);
guess.a_a_col = zeros(d*N,nq.arms);
guess.dFTtilde_col = zeros(d*N,NMuscle);
guess.Qdotdots_col = zeros(d*N,nq.all);
for k=1:N
    guess.a_col((k-1)*d+1:k*d,:) = repmat(guess.a(k,:),d,1); 
    guess.FTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.FTtilde(k,:),d,1);
    % guess.QsQdots_col((k-1)*d+1:k*d,:) = repmat(guess.QsQdots(k,:),d,1);
    guess.a_a_col((k-1)*d+1:k*d,:) = repmat(guess.a_a(k,:),d,1);
    guess.dFTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.dFTtilde(k,:),d,1);
    % guess.Qdotdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdotdots(k,:),d,1);
end
guess.Qs_col=Qs_spline_col.data(:,2:end);
% guess.Qs_col=interp1(Qs.allfilt(:,1),Qs.allfilt,)
if Options.preoptimizekneedofs==1
       fprintf('Preoptimize knee ty at col points ...\n')
       guess.Qs_col(:,jointi.knee_ty.r)=OptimizeForkneety(guess.Qs_col,jointi,KCF.allinterpfilt_col,N,F2,Options); 
       fprintf('knee ty at col points preoptimized...\n')
elseif (Options.preoptimizekneedofs==2)||(Options.preoptimizekneedofs==3)
       fprintf(['Preoptimize knee superior-inferior translation and ' ...
        'abduction-adduction at col points...\n'])
        guess.Qs_col(:,[jointi.knee_add.r jointi.knee_ty.r])=OptimizeForkneeaddty(guess.Qs_col,jointi,KCF.allinterpfilt_col,N,F2,Options);
        fprintf(['Knee superior-inferior translation and ' ...
        'abduction-adduction at col points pre-optimized \n'])
elseif Options.preoptimizekneedofs==0

else
    keyboard;
       
end

guess.Qdots_col=Qdots_spline_col.data(:,2:end);

fprintf(['Preoptimize pelvis tilt and ty and their velocities at coll points \n']);
if Options.KCFasinputstoExternalFunction
    [pelvis_rtilt_col,pelvis_ty_col, pelvis_rtilt_v_col,pelvis_ty_v_col]=OptimizeForPelvis_ty_rtilt(guess.Qs_col,guess.Qdots_col,Qs.colheaders,jointi,GRF,N,F1,F2_skeletal,f_contactForce,Options); % pre optimization of pelvis_ty and pelvis_tilt does not improve the convergence
else
    [pelvis_rtilt_col,pelvis_ty_col, pelvis_rtilt_v_col,pelvis_ty_v_col]=OptimizeForPelvis_ty_rtilt(guess.Qs_col,guess.Qdots_col,Qs.colheaders,jointi,GRF,N,F1,F2,f_contactForce,Options); % pre optimization of pelvis_ty and pelvis_tilt does not improve the convergence
end
fprintf(['pelvis tilt and ty and their velocities at coll points preoptimized \n']);
guess.Qs_col(:,jointi.pelvis.tilt)=pelvis_rtilt_col;
guess.Qs_col(:,jointi.pelvis.ty)=pelvis_ty_col; 
guess.Qs_col=guess.Qs_col./scaling.Qs;
guess.Qs_col(:,jointi.knee_ty.r)=(guess.Qs_col(:,jointi.knee_ty.r)-scaling.knee_ty.a)/scaling.knee_ty.b;

guess.Qdots_col(:,jointi.pelvis.tilt)=pelvis_rtilt_v_col;
guess.Qdots_col(:,jointi.pelvis.ty)=pelvis_ty_v_col; 
guess.Qdots_col=guess.Qdots_col./scaling.Qdots;
guess.Qdotdots_col=Qdotdots_spline_col.data(:,2:end)./scaling.Qdotdots;
guess.Qdotdots_col=guess.Qdotdots_col./scaling.Qdotdots;

guess.QsQdots_col(:,1:2:end) = guess.Qs_col;
guess.QsQdots_col(:,2:2:end) = guess.Qdots_col;

end

function knee_ty=OptimizeForkneety(Q,jointi,KCF,N,F2,Options)

nq=size(Q,2);
for i=1:size(Q,1)
    X=[];
    if Options.KCFasinputstoExternalFunction
        keyboard; %to write it
    else
        for j=1:nq
            X=[X Q(i,j)' 0];
        end
        Tk=full(F2([X,zeros(1,nq),ones(1,54)]));
        KCF_M0=Tk(48);
        KCF_L0=Tk(47);
    end
    
    q0_unsc=0.042;
    q0=(q0_unsc-0.042)/0.007;
    options_optim=[];
    options_optim.Display='none';
    qopt = lsqnonlin(@(x)fun(x,X,F2,jointi,KCF(i,:),nq,Options),q0,-1,1,options_optim);
    
    if Options.KCFasinputstoExternalFunction
        keyboard; %to write it
    else
        newX(i,:)=X;
        newX(i,jointi.knee_ty.r*2-1)=qopt*0.007+0.042;
        out=full(F2([newX(i,:),zeros(1,nq),ones(1,54)]));
        KCF_M(i)=out(48);
        KCF_L(i)=out(47);
    end

    knee_ty(i)=qopt*0.007+0.042;    
    
    debug_fun=fun(qopt,X,F2,jointi,KCF(i,:),nq,Options);
end

end

function out=OptimizeForkneeaddty(Q,jointi,KCF,N,F2,Options)

nq=size(Q,2);
for i=1:size(Q,1)
    i_change=0;

    X=[];
    if Options.KCFasinputstoExternalFunction
        X=zeros(1,6);
        X(1)=Q(i,jointi.knee_flex.r);
        X(5)=0.042;
        Tk=full(F2(X));
        KCF_M0=Tk(1);
        KCF_L0=Tk(2);
    else
        for j=1:nq
            X=[X Q(i,j)' 0];
        end
        X(jointi.knee_add.r*2-1)=0;
        X(jointi.knee_ty.r*2-1)=0.042;
        if F2.nnz_in==102
            Tk=full(F2([X,zeros(1,nq)]));
            KCF_M0=Tk(36);
            KCF_L0=Tk(35);
        else
            Tk=full(F2([X,zeros(1,nq),ones(1,54)]));
            KCF_M0=Tk(48);
            KCF_L0=Tk(47);
        end
    end
    
    q0_unsc=[0 0.042];
    q0(1)=q0_unsc(1)/0.1;
    q0(2)=(q0_unsc(2)-0.042)/0.007;
    options_optim=[];
    options_optim.Display='none';
    options_optim.Algorithm='trust-region-reflective';
    if Options.preoptimizekneedofs==2
        KCF_in=KCF(i,:);
    elseif Options.preoptimizekneedofs==3
        KCF_in=[200 200];
    else
    end

    [qopt,resnorm,residual] = lsqnonlin(@(x)fun2(x,X,F2,jointi,KCF_in,nq,Options),q0,[-1 -1],[1 1],options_optim);
    if Options.KCFasinputstoExternalFunction
        newX(i,:)=zeros(1,6);
        newX(i,1)=Q(i,jointi.knee_flex.r);
        newX(i,2)=qopt(1)*0.1;
        newX(i,5)=qopt(2)*0.007+0.042;
        out=full(F2(newX(i,:)));
        KCF_M(i)=out(1);
        KCF_L(i)=out(2);
    else
        newX(i,:)=X;
        newX(i,jointi.knee_add.r*2-1)=qopt(1)*0.1;
        newX(i,jointi.knee_ty.r*2-1)=qopt(2)*0.007+0.042;
        if F2.nnz_in==102
            out=full(F2([newX(i,:),zeros(1,nq)]));
            KCF_M(i)=out(36);
            KCF_L(i)=out(35);
        else
            out=full(F2([newX(i,:),zeros(1,nq),ones(1,54)]));
            KCF_M(i)=out(48);
            KCF_L(i)=out(47);
        end
    end

    dev=0.3;
    it=0;
    while any(abs(residual)>0.1)&(it<=120)
        if Options.KCFasinputstoExternalFunction
            newX(i,2)=qopt(1)*0.1;
            newX(i,5)=qopt(2)*0.007+0.042;  
            Tk=full(F2(newX(i,:)));
            KCF_M0=Tk(1);
            KCF_L0=Tk(2);
        else
            newX(i,jointi.knee_add.r*2-1)=qopt(1)*0.1;
            newX(i,jointi.knee_ty.r*2-1)=qopt(2)*0.007+0.042;    
            Tk=full(F2([newX(i,:),zeros(1,nq),ones(1,54)]));
            KCF_M0=Tk(48);
            KCF_L0=Tk(47);
        end
        if KCF_M0>KCF_L0
            q0_unsc(1)=q0_unsc(1)-dev*pi/180;
            i_change=-1;
        else
            q0_unsc(1)=q0_unsc(1)+dev*pi/180;
            i_change=1;
        end
        q0(1)=q0_unsc(1)/0.1;
        q0(2)=(q0_unsc(2)-0.042)/0.007;
        [qopt,resnorm,residual] = lsqnonlin(@(x)fun2(x,newX(i,:),F2,jointi,KCF_in,nq,Options),q0,[-1 -1],[1 1],options_optim,Options);
        
        
        it=it+1;
        if (it==30)
            options_optim.Algorithm='interior-point';
            q0=[0 0];
            q0_unsc(1)=q0(1)*0.1;
            q0_unsc(2)=q0(2)*0.007+0.042;    
        elseif it==60
            options_optim.Algorithm='trust-region-reflective';
            q0_unsc(1)=newX(i-1,jointi.knee_add.r*2-1);
            q0_unsc(2)=newX(i-1,jointi.knee_ty.r*2-1);
            q0(1)=q0_unsc(1)/0.1;
            q0(2)=(q0_unsc(2)-0.042)/0.007;
        elseif it==80
            options_optim.Algorithm='interior-point';
            q0_unsc(1)=newX(i-1,jointi.knee_add.r*2-1);
            q0_unsc(2)=newX(i-1,jointi.knee_ty.r*2-1);
            q0(1)=q0_unsc(1)/0.1;
            q0(2)=(q0_unsc(2)-0.042)/0.007;
        elseif it>100
            if KCF_M0<1 && KCF_L0<1
                q0_unsc(2)=q0_unsc(2)-0.01;
                q0(1)=q0_unsc(1)/0.1;
                q0(2)=(q0_unsc(2)-0.042)/0.007;
            end
        end
    end
    if Options.KCFasinputstoExternalFunction
        newX(i,2)=qopt(1)*0.1;
        newX(i,5)=qopt(2)*0.007+0.042;
        out=full(F2(newX(i,:)));
        KCF_M(i)=out(1);
        KCF_L(i)=out(2);
    else
        newX(i,jointi.knee_add.r*2-1)=qopt(1)*0.1;
        newX(i,jointi.knee_ty.r*2-1)=qopt(2)*0.007+0.042;    
        if F2.nnz_in==102
            out=full(F2([newX(i,:),zeros(1,nq)]));
            KCF_M(i)=out(36);
            KCF_L(i)=out(35);
        else
            out=full(F2([newX(i,:),zeros(1,nq),ones(1,54)]));
            KCF_M(i)=out(48);
            KCF_L(i)=out(47);
        end
    end

    knee_add(i)=qopt(1)*0.1;
    knee_ty(i)=qopt(2)*0.007+0.042;    
    
end

out=[knee_add' knee_ty'];

end

function f=fun(x,X,F2,jointi,KCFexp,nq,Options)
    if Options.KCFasinputstoExternalFunction
        X=zeros(1,6);
        X(5)=x*0.007+0.042;
        Tk=full(F2(X));
        KCF_M=Tk(1);
        KCF_L=Tk(2);
    else
        X(jointi.knee_ty.r*2-1)=x*0.007+0.042;
        if F2.nnz_in==102
            Tk=full(F2([X,zeros(1,nq)]));
            KCF_M=Tk(36);
            KCF_L=Tk(35);
        else
            Tk=full(F2([X,zeros(1,nq),ones(1,54)]));
            KCF_M=Tk(48);
            KCF_L=Tk(47);
        end
    end
    f=(KCFexp-[KCF_M KCF_L])/1000;
end

function f=fun2(x,X,F2,jointi,KCFexp,nq,Options)
    if Options.KCFasinputstoExternalFunction
        X(2)=x(1)*0.1;
        X(5)=x(2)*0.007+0.042;
        Tk=full(F2(X));
        KCF_M=Tk(1);
        KCF_L=Tk(2);
    else
        X(jointi.knee_add.r*2-1)=x(1)*0.1;
        X(jointi.knee_ty.r*2-1)=x(2)*0.007+0.042;
        if F2.nnz_in==102
            Tk=full(F2([X,zeros(1,nq)]));
            KCF_M=Tk(36);
            KCF_L=Tk(35);
        else
            Tk=full(F2([X,zeros(1,nq),ones(1,54)]));
            KCF_M=Tk(48);
            KCF_L=Tk(47);
        end
    end
    f=(KCFexp-[KCF_M KCF_L])/100;
end

function [pelvis_rtilt,pelvis_ty,pelvis_rtilt_v,pelvis_ty_v]=OptimizeForPelvis_ty_rtilt(Q,Qdots,colheaders,jointi,GRF,N,F1,F2_skeletal,f_contactForce,Options)

lb=[-0.2 0.7 -2 -2];
ub=[0.4 1.2 2 2];
nq=size(Q,2);
npoints=size(Q,1);
if npoints==N+1
    GRFexp=GRF.val.allinterp;
else
    GRFexp=GRF.val.allinterp_col;
end

for i=1:npoints
    X=[];

    for j=1:nq
        X=[X Q(i,j)' Qdots(i,j)'];
    end
    
    q0=[Q(i,[jointi.pelvis.tilt  jointi.pelvis.ty]) Qdots(i,[jointi.pelvis.tilt  jointi.pelvis.ty])];
    options_optim.Display='none';
    [qopt,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = lsqnonlin(@(x)fun_preoptGRF(x,X,F1,F2_skeletal,f_contactForce,jointi,GRFexp(i,:),nq,Q(i,:),Qdots(i,:),Options),q0,lb,ub,options_optim);

    pelvis_rtilt(i,1)=qopt(1);
    pelvis_ty(i,1)=qopt(2);
    pelvis_rtilt_v(i,1)=qopt(3);
    pelvis_ty_v(i,1)=qopt(4);

    

end


end

function f=fun_preoptGRF(x,X,F1,F2_skeletal,f_contactForce,jointi,GRF_exp,nq,Qsexp,Qdotsexp,Options)

X(:,jointi.pelvis.tilt*2-1)=x(1);
X(:,jointi.pelvis.ty*2-1)=x(2);
X(:,jointi.pelvis.tilt*2)=x(3);
X(:,jointi.pelvis.ty*2)=x(4);

% Sphere parameters
% Fixed parameters
dissipation = 2; %check this value
normal = [0,1,0];
transitionVelo = 0.2;
staticFriction = 0.8;
dynamicFriction = 0.8;
viscousFriction = 0.5; 
stif = 1000000;  
loc.s1.y    = -0.021859;
loc.s2.y    = -0.021859;
loc.s3.y    = -0.021859;
loc.s4.y    = -0.0214476;
loc.s5.y    = -0.021859;
loc.s6.y    = -0.0214476;
B_locSphere_s1_r    = [0.00190115788407966, -0.00382630379623308];
B_locSphere_s2_r    = [0.148386399942063, -0.028713422052654];
B_locSphere_s3_r    = [0.133001170607051, 0.0516362473449566];
B_locSphere_s4_r    = [0.06, -0.0187603084619177];    
B_locSphere_s5_r    = [0.0662346661991635, 0.0263641606741698];
B_locSphere_s6_r    = [0.045, 0.0618569567549652];
IG_rad              = 0.032*ones(1,6); 
guess.params = [B_locSphere_s1_r,B_locSphere_s2_r,...
    B_locSphere_s3_r,B_locSphere_s4_r,B_locSphere_s5_r,B_locSphere_s6_r,...
    IG_rad];
% External function: F1
calcn.l.omega = 1:3;
calcn.l.v_lin = 4:6;
calcn.r.omega = 7:9;
calcn.r.v_lin = 10:12;
calcn.l.pos = 13:15;
calcn.r.pos = 16:18;
calcn.TR.R.l = 19:27;
calcn.TR.T.l = 28:30;
calcn.TR.R.r = 31:39;
calcn.TR.T.r = 40:42;
toes.l.omega = 43:45;
toes.l.v_lin = 46:48;
toes.r.omega = 49:51;
toes.r.v_lin = 52:54;
toes.l.pos = 55:57;
toes.r.pos = 58:60;
toes.TR.R.l = 61:69;
toes.TR.T.l = 70:72;
toes.TR.R.r = 73:81;
toes.TR.T.r = 82:84;
% Indices variable parameters
loci.s1.r.x = 1;
loci.s1.r.z = 2;
loci.s2.r.x = 3;
loci.s2.r.z = 4;
loci.s3.r.x = 5;
loci.s3.r.z = 6;
loci.s4.r.x = 7;
loci.s4.r.z = 8;
loci.s5.r.x = 9;
loci.s5.r.z = 10;
loci.s6.r.x = 11;
loci.s6.r.z = 12;
radi.s1 = 13;
radi.s2 = 14;
radi.s3 = 15;
radi.s4 = 16;
radi.s5 = 17;
radi.s6 = 18;
paramsCM_nsc=guess.params;
locSphere_s1_l = [paramsCM_nsc(loci.s1.r.x),...
    loc.s1.y,-paramsCM_nsc(loci.s1.r.z)]';    
locSphere_s1_r = [paramsCM_nsc(loci.s1.r.x),...
    loc.s1.y,paramsCM_nsc(loci.s1.r.z)]';        
locSphere_s2_l = [paramsCM_nsc(loci.s2.r.x),...
    loc.s2.y,-paramsCM_nsc(loci.s2.r.z)]';
locSphere_s2_r = [paramsCM_nsc(loci.s2.r.x),...
    loc.s2.y,paramsCM_nsc(loci.s2.r.z)]';        
locSphere_s3_l = [paramsCM_nsc(loci.s3.r.x),...
    loc.s3.y,-paramsCM_nsc(loci.s3.r.z)]';
locSphere_s3_r = [paramsCM_nsc(loci.s3.r.x),...
    loc.s3.y,paramsCM_nsc(loci.s3.r.z)]';        
locSphere_s4_l = [paramsCM_nsc(loci.s4.r.x),...
    loc.s4.y,-paramsCM_nsc(loci.s4.r.z)]';
locSphere_s4_r = [paramsCM_nsc(loci.s4.r.x),...
    loc.s4.y,paramsCM_nsc(loci.s4.r.z)]';          
locSphere_s5_l = [paramsCM_nsc(loci.s5.r.x),...
    loc.s5.y,-paramsCM_nsc(loci.s5.r.z)]';
locSphere_s5_r = [paramsCM_nsc(loci.s5.r.x),...
    loc.s5.y,paramsCM_nsc(loci.s5.r.z)]';          
locSphere_s6_l = [paramsCM_nsc(loci.s6.r.x),...
    loc.s6.y,-paramsCM_nsc(loci.s6.r.z)]';
locSphere_s6_r = [paramsCM_nsc(loci.s6.r.x),...
    loc.s6.y,paramsCM_nsc(loci.s6.r.z)]';          
% Compute contact forces

    [outk1] = F1(X); 

    force_s1_l = f_contactForce(stif,paramsCM_nsc(radi.s1),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s1_l,...
        outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
        outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
    force_s2_l = f_contactForce(stif,paramsCM_nsc(radi.s2),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s2_l,...
        outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
        outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
    force_s3_l = f_contactForce(stif,paramsCM_nsc(radi.s3),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s3_l,...
        outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
        outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
    force_s4_l = f_contactForce(stif,paramsCM_nsc(radi.s4),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s4_l,...
        outk1(toes.l.pos),outk1(toes.l.v_lin),...
        outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l));        
    force_s5_l = f_contactForce(stif,paramsCM_nsc(radi.s5),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s5_l,...
        outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
        outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
    force_s6_l = f_contactForce(stif,paramsCM_nsc(radi.s6),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s6_l,...
        outk1(toes.l.pos),outk1(toes.l.v_lin),...
        outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l));        
    force_s1_r = f_contactForce(stif,paramsCM_nsc(radi.s1),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s1_r,...
        outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
        outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
    force_s2_r = f_contactForce(stif,paramsCM_nsc(radi.s2),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s2_r,...
        outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
        outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
    force_s3_r = f_contactForce(stif,paramsCM_nsc(radi.s3),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s3_r,...
        outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
        outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
    force_s4_r = f_contactForce(stif,paramsCM_nsc(radi.s4),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s4_r,...
        outk1(toes.r.pos),outk1(toes.r.v_lin),...
        outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r));        
    force_s5_r = f_contactForce(stif,paramsCM_nsc(radi.s5),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s5_r,...
        outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
        outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
    force_s6_r = f_contactForce(stif,paramsCM_nsc(radi.s6),...
        dissipation,normal,transitionVelo,staticFriction,...
        dynamicFriction,viscousFriction,locSphere_s6_r,...
        outk1(toes.r.pos),outk1(toes.r.v_lin),...
        outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r));
    in_F2 = [force_s1_l,force_s2_l,force_s3_l,force_s4_l,...
        force_s5_l,force_s6_l,force_s1_r,force_s2_r,force_s3_r,...
        force_s4_r,force_s5_r,force_s6_r,...       
        paramsCM_nsc(loci.s1.r.x),paramsCM_nsc(loci.s1.r.z),...
        paramsCM_nsc(loci.s2.r.x),paramsCM_nsc(loci.s2.r.z),...
        paramsCM_nsc(loci.s3.r.x),paramsCM_nsc(loci.s3.r.z),...
        paramsCM_nsc(loci.s4.r.x),paramsCM_nsc(loci.s4.r.z),...
        paramsCM_nsc(loci.s5.r.x),paramsCM_nsc(loci.s5.r.z),...
        paramsCM_nsc(loci.s6.r.x),paramsCM_nsc(loci.s6.r.z),...
        paramsCM_nsc(radi.s1:radi.s6)];
    % Relevant Indices for F2, GRFs
    GRFi.r              = 35:37;
    GRFi.l              = 38:40;
    GRFi.all            = [GRFi.r,GRFi.l];
     % The second external function (F2) returns joint torques, GRFs,
            % and GRMs based on joint states, controls, contact forces, and
            % several parameters of the contact models.
    if Options.KCFasinputstoExternalFunction
        [Tj] = F2_skeletal([X';zeros(nq,1);in_F2';zeros(6,1)]);
    else
        [Tj] = F2_skeletal([X';zeros(nq,1);in_F2']);
    end
    GRF_v_r=full(Tj(GRFi.r(2)));
    GRF_v_l=full(Tj(GRFi.l(2)));

    f=[GRF_v_r-GRF_exp(3) GRF_v_l-GRF_exp(6)];
    f=[f x(1)-Qsexp(jointi.pelvis.tilt) x(2)-Qsexp(jointi.pelvis.ty*2-1) x(3)-Qdotsexp(jointi.pelvis.tilt) x(4)-Qdotsexp(jointi.pelvis.ty*2-1)];

end