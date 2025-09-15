% This script contains several CasADi-based functions that are
% used when solving the OCPs
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
import casadi.*

if Options.usedCompiled_lMT_vMT_External
else
    %% Polynomial approximation
    pathpolynomial = [pathRepo,'\Polynomials_GC'];
    addpath(genpath(pathpolynomial));
    muscle_spanning_info_m = muscle_spanning_joint_INFO(musi_pol,:);
    MuscleInfo_m.muscle  = MuscleInfo.muscle(musi_pol);                  
    qin     = SX.sym('qin',1,nq.legright);
    qdotin  = SX.sym('qdotin',1,nq.legright);
    lMT     = SX(NMuscle_pol,1);
    vMT     = SX(NMuscle_pol,1);
    dM      = SX(NMuscle_pol,nq.legright);
    for i=1:NMuscle_pol      
    %     index_dof_crossing_withseckneedofs  = find(muscle_spanning_info_m(i,:)==1);
    %     internaldof_knee=any([index_dof_crossing_withseckneedofs==4 ...
    %             index_dof_crossing_withseckneedofs==5 index_dof_crossing_withseckneedofs==7 ...
    %             index_dof_crossing_withseckneedofs==8 index_dof_crossing_withseckneedofs==9]);
    %     
    %     index_dof_crossing_withoutseckneedofs=index_dof_crossing_withseckneedofs;
    %     index_dof_crossing_withoutseckneedofs(index_dof_crossing_withoutseckneedofs==4)=[];
    %     index_dof_crossing_withoutseckneedofs(index_dof_crossing_withoutseckneedofs==5)=[];
    %     index_dof_crossing_withoutseckneedofs(index_dof_crossing_withoutseckneedofs==7)=[];
    %     index_dof_crossing_withoutseckneedofs(index_dof_crossing_withoutseckneedofs==8)=[];
    %     index_dof_crossing_withoutseckneedofs(index_dof_crossing_withoutseckneedofs==9)=[];
        index_dof_crossing_withknee = ...
            find(muscle_spanning_info_m(i,:)==1);    
        index_dof_crossing_withoutknee = ...
            find(muscle_spanning_info_m(i,[1:3 6 10:end])==1);
        order = MuscleInfo_m.muscle(i).order;

        if any(index_dof_crossing_withknee==15)
            [mat,diff_mat_q]    = n_art_mat_9_GC_cas_SX([qin(1,index_dof_crossing_withknee(1:end-1)) 0],...
                    order);
        else
            [mat,diff_mat_q]    = n_art_mat_9_GC_cas_SX(qin(1,index_dof_crossing_withknee),...
                order);
        end


    %     [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing_withoutknee),...
    %         order);
        I_nonzeros=find(MuscleInfo_m.muscle(i).coeff{order});
        lMT(i,1)            = mat(I_nonzeros)*MuscleInfo_m.muscle(i).coeff{order}(I_nonzeros);
        vMT(i,1)            = 0;
        dM(i,1:nq.legright)      = 0;
        nr_dof_crossing     = length(index_dof_crossing_withknee); 
        for dof_nr = 1:nr_dof_crossing
            if index_dof_crossing_withknee(dof_nr)==15
            else
                dM(i,index_dof_crossing_withknee(dof_nr)) = ...
                    (-(diff_mat_q(I_nonzeros,dof_nr)))'*MuscleInfo_m.muscle(i).coeff{order}(I_nonzeros);
                vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing_withknee(dof_nr))*...
                    qdotin(1,index_dof_crossing_withknee(dof_nr)));
            end
        end 
    end
    f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});
end

%% Normalized sum of squared values
% Function for 2 elements 
etemp2 = SX.sym('etemp2',2);
Jtemp2 = 0;
for i=1:length(etemp2)
    Jtemp2 = Jtemp2 + etemp2(i).^2;
end
Jtemp2 = Jtemp2/2;
f_J2 = Function('f_J2',{etemp2},{Jtemp2});
% Function for 5 elements 
etemp5 = SX.sym('etemp5',5);
Jtemp5 = 0;
for i=1:length(etemp5)
    Jtemp5 = Jtemp5 + etemp5(i).^2;
end
Jtemp5 = Jtemp5/5;
f_J5 = Function('f_J5',{etemp5},{Jtemp5});
% Function for 6 elements 
etemp6 = SX.sym('etemp6',6);
Jtemp6 = 0;
for i=1:length(etemp6)
    Jtemp6 = Jtemp6 + etemp6(i).^2;
end
Jtemp6 = Jtemp6/6;
f_J6 = Function('f_J6',{etemp6},{Jtemp6});
% Function for 92 elements 
etemp92 = SX.sym('etemp92',92);
Jtemp92 = 0;
for i=1:length(etemp92)
    Jtemp92 = Jtemp92 + etemp92(i).^2;
end
Jtemp92 = Jtemp92/92;
f_J92 = Function('f_J92',{etemp92},{Jtemp92});
% Function for 94 elements 
etemp94 = SX.sym('etemp94',94);
Jtemp94 = 0;
for i=1:length(etemp94)
    Jtemp94 = Jtemp94 + etemp94(i).^2;
end
Jtemp94 = Jtemp94/94;
f_J94 = Function('f_J94',{etemp94},{Jtemp94});
% Function for 21 elements
etemp21 = SX.sym('etemp21',21);
Jtemp21 = 0;
for i=1:length(etemp21)
    Jtemp21 = Jtemp21 + etemp21(i).^2;
end
Jtemp21 = Jtemp21/21;
f_J21 = Function('f_J21',{etemp21},{Jtemp21});
% Function for 22 elements
etemp22 = SX.sym('etemp22',22);
Jtemp22 = 0;
for i=1:length(etemp22)
    Jtemp22 = Jtemp22 + etemp22(i).^2;
end
Jtemp22 = Jtemp22/22;
f_J22 = Function('f_J22',{etemp22},{Jtemp22});
% Function for 28 elements
etemp28 = SX.sym('etemp28',28);
Jtemp28 = 0;
for i=1:length(etemp28)
    Jtemp28 = Jtemp28 + etemp28(i).^2;
end
Jtemp28 = Jtemp28/28;
f_J28 = Function('f_J28',{etemp28},{Jtemp28});
% Function for 33 elements
etemp33 = SX.sym('etemp33',33);
Jtemp33 = 0;
for i=1:length(etemp33)
    Jtemp33 = Jtemp33 + etemp33(i).^2;
end
Jtemp33 = Jtemp33/33;
f_J33 = Function('f_J33',{etemp33},{Jtemp33});
% Function for 23 elements
etemp23 = SX.sym('etemp23',23);
Jtemp23 = 0;
for i=1:length(etemp23)
    Jtemp23 = Jtemp23 + etemp23(i).^2;
end
Jtemp23 = Jtemp23/23;
f_J23 = Function('f_J23',{etemp23},{Jtemp23});
% Function for 29 elements
etemp29 = SX.sym('etemp29',29);
Jtemp29 = 0;
for i=1:length(etemp29)
    Jtemp29 = Jtemp29 + etemp29(i).^2;
end
Jtemp29 = Jtemp29/29;
f_J29 = Function('f_J29',{etemp29},{Jtemp29});
% Function for 34 elements
etemp34 = SX.sym('etemp34',34);
Jtemp34 = 0;
for i=1:length(etemp34)
    Jtemp34 = Jtemp34 + etemp34(i).^2;
end
Jtemp34 = Jtemp34/34;
f_J34 = Function('f_J34',{etemp34},{Jtemp34});

%% Sum of products 
% Function for 27 elements 
ma_temp27 = SX.sym('ma_temp27',27);
ft_temp27 = SX.sym('ft_temp27',27);
J_sptemp27 = 0;
for i=1:length(ma_temp27)
    J_sptemp27 = J_sptemp27 + ma_temp27(i,1)*ft_temp27(i,1);    
end
f_T27 = Function('f_T27',{ma_temp27,ft_temp27},{J_sptemp27});
% Function for 28 elements 
ma_temp28 = SX.sym('ma_temp28',28);
ft_temp28 = SX.sym('ft_temp27',28);
J_sptemp28 = 0;
for i=1:length(ma_temp28)
    J_sptemp28 = J_sptemp28 + ma_temp28(i,1)*ft_temp28(i,1);    
end
f_T28 = Function('f_T28',{ma_temp28,ft_temp28},{J_sptemp28});
% Function for 13 elements 
ma_temp13 = SX.sym('ma_temp13',13);
ft_temp13 = SX.sym('ft_temp13',13);
J_sptemp13 = 0;
for i=1:length(ma_temp13)
    J_sptemp13 = J_sptemp13 + ma_temp13(i,1)*ft_temp13(i,1);    
end
f_T13 = Function('f_T13',{ma_temp13,ft_temp13},{J_sptemp13});
% Function for 12 elements 
ma_temp12 = SX.sym('ma_temp12',12);
ft_temp12 = SX.sym('ft_temp12',12);
J_sptemp12 = 0;
for i=1:length(ma_temp12)
    J_sptemp12 = J_sptemp12 + ma_temp12(i,1)*ft_temp12(i,1);    
end
f_T12 = Function('f_T12',{ma_temp12,ft_temp12},{J_sptemp12});
% Function for 6 elements 
ma_temp6 = SX.sym('ma_temp6',6);
ft_temp6 = SX.sym('ft_temp6',6);
J_sptemp6 = 0;
for i=1:length(ma_temp6)
    J_sptemp6 = J_sptemp6 + ma_temp6(i,1)*ft_temp6(i,1);    
end
f_T6 = Function('f_T6',{ma_temp6,ft_temp6},{J_sptemp6});

%% Arm activation dynamics
e_a = SX.sym('e_a',nq.arms); % arm excitations
a_a = SX.sym('a_a',nq.arms); % arm activations
dadt = ArmActivationDynamics(e_a,a_a);
f_ArmActivationDynamics = ...
    Function('f_ArmActivationDynamics',{e_a,a_a},{dadt});

%% Muscle contraction dynamics
pathmusclemodel = [pathRepo,'\MuscleModel'];
addpath(genpath(pathmusclemodel));
% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde',NMuscle); % Normalized tendon forces
a           = SX.sym('a',NMuscle); % Muscle activations
dFTtilde    = SX.sym('dFTtilde',NMuscle); % Time derivative tendon forces
lMT         = SX.sym('lMT',NMuscle); % Muscle-tendon lengths
vMT         = SX.sym('vMT',NMuscle); % Muscle-tendon velocities
% tension_SX  = SX.sym('tension',NMuscle); % Tensions not used here
Hilldiff    = SX(NMuscle,1); % Hill-equilibrium
FT          = SX(NMuscle,1); % Tendon forces
Fce         = SX(NMuscle,1); % Contractile element forces
Fiso        = SX(NMuscle,1); % Normalized forces from force-length curve
vMmax       = SX(NMuscle,1); % Maximum contraction velocities
% massM       = SX(NMuscle,1); % Muscle mass
% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam
for m = 1:NMuscle
    [Hilldiff(m),FT(m),Fce(m),Fiso(m),vMmax(m)] = ...
        ForceEquilibrium_FtildeState_GC(a(m),FTtilde(m),dFTtilde(m),...
        lMT(m),vMT(m),MTparameters_m(:,m),Fvparam,Fpparam,Faparam);
        %tension_SX(m));
end
f_forceEquilibrium_FtildeState = ...
    Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
    lMT,vMT},{Hilldiff,FT,Fce,Fiso,vMmax});

%% Passive joint torques
K_pass      = SX.sym('K_pass',4);
theta_pass  = SX.sym('theta_pass',2);
qin_pass    = SX.sym('qin_pass',1);
qdotin_pass = SX.sym('qdotin_pass',1);
% theta_pass 1 and 2 are inverted on purpose.
Tau_pass = K_pass(1,1)*exp(K_pass(2,1)*(qin_pass-theta_pass(2,1))) + ...
    K_pass(3,1)*exp(K_pass(4,1)*(qin_pass-theta_pass(1,1))) ...
    - 0.001*qdotin_pass;
f_PassiveMoments = Function('f_PassiveMoments',{K_pass,theta_pass,...
    qin_pass,qdotin_pass},{Tau_pass});

%% Passive joint torques for knee internal dofs
%Knee add and rot
% Tau_pass_kneeintmom=-30*((qin_pass/(5*pi/180)).^5);
Tau_pass_kneeintmom=-30*((qin_pass/(7.5*pi/180)).^5);
F_pass_txtz=-200*((qin_pass/0.02).^5); 
if Options.dampingInKneeSec==1
    Tau_pass_kneeintmom=Tau_pass_kneeintmom-qdotin_pass*0.001;
    F_pass_txtz=F_pass_txtz-qdotin_pass*1000;
    f_passiveMoments_kneeintmom = Function('f_passiveMoments_kneeintmom',{qin_pass,qdotin_pass},{Tau_pass_kneeintmom});
    f_passiveForce_kneeintf = Function('f_passiveForce_kneeintf', {qin_pass,qdotin_pass},{F_pass_txtz});
else
    f_passiveMoments_kneeintmom = Function('f_passiveMoments_kneeintmom',{qin_pass},{Tau_pass_kneeintmom});
    f_passiveForce_kneeintf = Function('f_passiveForce_kneeintf', {qin_pass},{F_pass_txtz});
end

%% Unscaling
x_bar18      = SX.sym('x_bar18_SX',18,1); 
sc_v18       = SX.sym('sc_v18_SX',18,1);
sc_r18       = SX.sym('sc_r18',18,1);
x18          = (x_bar18 - sc_r18)./sc_v18;
f_nsc18      = Function('f_nsc18',{x_bar18,sc_v18,sc_r18},{x18});

%% Contact forces
stiffnessSX         = SX.sym('stiffnessSX',1);
radiusSX            = SX.sym('radiusSX',1);
dissipationSX       = SX.sym('dissipationSX',1);
normalSX            = SX.sym('normalSX',1,3);
transitionVelocitySX= SX.sym('transitionVelocitySX',1);
staticFrictionSX    = SX.sym('staticFrictionSX',1);
dynamicFrictionSX   = SX.sym('dynamicFrictionSX',1);
viscousFrictionSX   = SX.sym('viscousFrictionSX',1);
spherePosSX         = SX.sym('spherePosSX',3);
orFramePosSX        = SX.sym('orFramePosSX',3);
v_linSX             = SX.sym('v_linSX',3);
omegaSX             = SX.sym('omegaSX',3);
RotSX               = SX.sym('RotSX',9);
TrSX                = SX.sym('TrSX',3);
% Hunt-Crossley contact model
forceSX = HCContactModel(stiffnessSX,radiusSX,dissipationSX,...
    normalSX,transitionVelocitySX,staticFrictionSX,...
    dynamicFrictionSX,viscousFrictionSX,spherePosSX,orFramePosSX,...
    v_linSX,omegaSX,RotSX,TrSX);
f_contactForce = Function('f_contactForce',{stiffnessSX,radiusSX,...
    dissipationSX,normalSX,transitionVelocitySX,staticFrictionSX,...
    dynamicFrictionSX,viscousFrictionSX,spherePosSX,orFramePosSX,...
    v_linSX,omegaSX,RotSX,TrSX},{forceSX});
