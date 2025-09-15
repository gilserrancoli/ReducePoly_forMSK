function Results_3D = TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix)


%% Muscle-driven simulations including mesh-based contact
    import casadi.*
    
    
    solveProblem = true; %Set true to solve the optimal control problem.
    analyseResults = true; %Set true to analyze the results.
    loadResults = false; %Set true to load the results of the optimization.
    saveResults = true; %Set true to save the results of the optimization.
    checkBoundsIG = false; % visualize guess-bounds
    writeMotionFiles = true; %Set true to write motion files for use in OpenSim GUI
    saveOptimalTrajectories = true; %Set true to save optimal trajectories
    writeIKmotion=true; %Set true to write .mot file
    
    subject = 'subject_GC';
    parallelMode = 'thread';
    NThreads = 8; % Number of threads used in parallel.
    
    pathmain = pwd;
    [pathRepo,~,~] = fileparts(pathmain);
    pathSettings = [pathRepo,'/Settings'];
    addpath(genpath(pathSettings));
    % Load settings
    % settings_3D
    % ww = 2; % Index row in matrix settings - Recorder, Approximated Hessian, mumps, data informed initial guess
    %%

    % List of movements
        nametrial.id    = nametrial_id;%'bouncy7'; % Experimental walking trial to track 
        % (available: ngait_og1, ngait_og5, bouncy4, bouncy7, mtpgait3, 
        % mtpgait9, ngait_tm_fast1, ngait_tm_set1) ngait_og1_old to use the
        % previous data
    
    %Use Tangential force with static and dynamic coefficients of 0.8 and 0.8,
        %or 0.1 and 0.1
        Options.useContactModelwithTangentialForce0808=0;
        Options.useContactModelwithTangentialForce0101=0;
    %Related to joint force and moment residuals
        Options.useKCFresiduals=0;
            Options.useKCFconstraintsasBounds=0;
        Options.useJointResMom=0; % originally 1
            Options.useJRMconstraitsasBounds=1;
        Options.usePelvisResMom=1;
        Options.dampingInKneeSec=1; %use damping at the knee secondary dofs
    %Related to polynomials approximating muscle-tendon lengths
        Options.useReducedPolynomials=useReducedPolynomials;
            Options.useReoptimizedPoly=1;
        Options.usedCompiled_lMT_vMT_External=0;
        Options.err_poly = err_poly; % Define option before running the code
    %Related to mesh-based contact model
        %Choose .dll
        Options.maxsmoothness='MellowMax'; %options: logSum, MellowMax, nosmooth
            Options.kInmaxpen=1e4; %with logsum valid numbers 1e3 for TrackSim_2_kneeCont_kmax1e3.dll, 1e5 TrackSim_2_kneeCont_logSum_kmax1e5_kpress5e4.dll, 5e5 TrackSim_2_kneeCont_kmax1e5.dll, 1e4 for TrackSim_2_kneeCont_logSum_kmax1e4_kpress5e4.dll,5e4 for TrackSim_2_kneeCont_logSum_kmax5e4_kpress5e4.dll 
            % 'nosmooth' for non-smooth version TrackSim_2_kneeCont_maxnosmooth.dll
            Options.kInpress=1e4; %with mellowmax valid for 5e4, 5e5, 1e4
        Options.kInCheckContacts=1e3; %valid values 1e3, 5e3 and 1e4
        Options.nfacesTib=100;
        Options.nfacesFem=Options.nfacesFem;
        Options.rad4Pairs=1;
        Options.numpairs=[]; %1378, 499... old way to search for selected .dll
        Options.setLBKCF=0; %low bounds for KCF
    %Related to IG
        IGi=2; % Data-informed initial guess
        Options.multiple_IG=0   ;
        Options.IGn=1;
        Options.force2createNewIGs=0;
        Options.noiseApplied2NewIGs=1; %percentage before was 1
        Options.preoptimizekneedofs=3; %0 is considering all dofs as 
        % experimental data or 0 values, 1 is preoptimizing knee ty, 2 is 
        % preoptimizing knee adduction and knee ty to match knee contact
        % forces, 3 is preoptimizing knee adduction and ty to start with
        % 200 N of knee contact force at each compartment
    %Related to how external function is handled
        Options.KCFasinputstoExternalFunction=1;
    hessi=1; %approximate Hession hessi = 2 is the exact hessian
    deri=1; %use of AD Recorder
    Options.nametrial=nametrial;
    
    % % Variable parameters
    % deri        = settings(ww,1); % derivative supplier identifier
    % deri=1;
    % hessi       = settings(ww,2); % Hessian identifier
    % linsoli     = settings(ww,3); % linear solver identifier
    % IGi         = settings(ww,4); % initial guess identifier
    
    % Fixed parameters
    % W.Qs        = 125;	% weight joint kinematic, % before 50
    % W.Qdots     = 25;    % weight joint velocities % before 10
    % W.GRF       = 10;	% weight ground reaction forces
    % W.GRM       = 10;    % weight ground reaction moments
    % W.ID_act    = 0;   % weight joint kinetics %before 10
    % W.KCF       = 75;   % weight KCF before was 10
    % W.a         = 5;    % weight muscle activations
    % W.minPelvisRes=2;
    % W.u         = 0.03;
    % W.u_qd2dot  = 0.003;
    % W.u_vA      = 0.52;
    % 
    if Options.useKCFresiduals
        W.KCF_res=1e-1;
    else
        W.KCF_res=0;
    end
    if Options.useJointResMom
        W.JRM_res=1;
    else
        W.JRM_res=0;
    end
    dev_cm.loc  = 25;   % allowed deviation from generic locations
    dev_cm.rad  = 50;   % allowed deviation from generic radii
    tol_ipopt   = 4;    % tolerance (means 1e-(tol_ipopt))
    N           = 40;   % number of mesh intervals


    
    % Derivative supplier
    % if deri == 1
        setup.derivatives = 'AD_Recorder'; % Algorithmic differentiation / Recorder     
    % elseif deri == 2
    %     setup.derivatives = 'AD_ADOLC'; % Algorithmic differentiation / ADOL-C
    % elseif deri == 3
    %     setup.derivatives = 'FD'; % Finite differences
    % end
    % Identifiers for experimental data
    nametrial.ID    = ['ID_',nametrial.id];
    nametrial.GRF   = ['GRF_',nametrial.id];
    nametrial.IK    = ['IK_',nametrial.id];
    nametrial.KCF   = ['KCF_',nametrial.id];
    switch nametrial.id
        case 'ngait_og1'
            time_opt = [2.2292 3.3740]; % right leg cycle [1.680833, 2.775] --> need left GRF at the beginning
        case 'ngait_og1_old'
            time_opt = [2.2292 3.3740]; % right leg cycle [1.680833, 2.775] --> need left GRF at the beginning
        case 'ngait_og5'
            time_opt = [2.7158 3.919]; 
        case 'bouncy4'
            time_opt = [3.5433 4.825];
        case 'bouncy7'
            time_opt = [0.7175 1.965];
        case 'mtpgait3'
            time_opt = [2.1508 3.324];
        case 'mtpgait9'
            time_opt = [2.2583 3.482];
        case 'ngait_tm_fast1'
            time_opt = [4.7067 5.7983]; %one cycle apparently neat
        case 'ngait_tm_set1'
            time_opt = [2.2000 3.5033];
    end  
    % Available linear solvers
    linear_solvers = {'mumps','ma27','ma57','ma77','ma86','ma97'}; 
    if Options.useReducedPolynomials
        if Options.useReoptimizedPoly
            poly=['red' num2str(Options.err_poly) '_reopt'];
        else
            poly=['red' num2str(Options.err_poly)];
        end
    else
        poly=['full' num2str(Options.err_poly)];
    end
    % The filename used to save the results depends on the settings 
    if ~exist('savename_suffix', 'var')
        savename_suffix = '';
    end
    savename = ['_', nametrial.id, '_', num2str(Options.IGn), '_poly', poly, savename_suffix];
    savename2 = ['_', nametrial.id, '_', num2str(Options.IGn), '_poly', poly];

    %% Load external functions
    % The external function performs inverse dynamics through the
    % OpenSim/Simbody C++ API. This external function is compiled as a dll from
    % which we create a Function instance using CasADi in MATLAB. 
    % We use different external functions. A first external function extracts 
    % several parameters of the bodies to which the contact spheres are attached.
    % The contact forces are then computed in MATLAB and are inputs of the
    % second external function in which the skeleton dynamics is described. The
    % motivation for this decoupling is to limit the number of times we need to
    % build the model. By defining the contact model in MATLAB, we only need to
    % build the model once per external function, whereas keeping the contact
    % model in the external function would require re-building the model during
    % the optimization.
    pathExternalFunctions = [pathRepo,'/ExternalFunctions_GC'];
    % Loading external functions. 
    % cd(pathExternalFunctions);
    if exist(pathExternalFunctions, 'dir')
       cd(pathExternalFunctions);
    else
        warning('Cannot change directory: %s does not exist.', pathExternalFunctions);
        % Optionally: create the folder if it’s supposed to exist
        % mkdir(pathExternalFunctions);
    end
    switch setup.derivatives
        case 'AD_Recorder'   
            if hessi == 1
                F1 = external('F','TrackSim_1_kneeCont.dll'); 
                if Options.useContactModelwithTangentialForce0101
                    F2 = external('F','TrackSim_2_kneeCont_withTanForce_us01ud01.dll'); 
                elseif Options.useContactModelwithTangentialForce0808
                    F2 = external('F','TrackSim_2_kneeCont_withTanForce_us08ud08.dll'); 
                else
                    if strcmp(Options.maxsmoothness,'logSum')
                        if Options.kInmaxpen==1e3
                            F2 = external('F','TrackSim_2_kneeCont_kmax1e3.dll'); 
                        elseif Options.kInmaxpen==5e3
                            F2 = external('F','TrackSim_2_kneeCont_kmax5e3.dll'); 
                        elseif Options.kInmaxpen==1e4
                            F2 = external('F','TrackSim_2_kneeCont_logSum_kmax1e4_kpress5e4.dll'); 
                        elseif Options.kInmaxpen==5e4
                            F2 = external('F','TrackSim_2_kneeCont_logSum_kmax5e4_kpress5e4.dll');
                        elseif Options.kInmaxpen==1e5
                            F2 = external('F','TrackSim_2_kneeCont_logSum_kmax1e5_kpress5e4.dll');
                        else
                            keyboard; %no function found
                        end
                    elseif strcmp(Options.maxsmoothness,'MellowMax')
                        if  Options.kInmaxpen==1e5
                            F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e5_kpress5e4.dll');
                        elseif Options.kInmaxpen==1e4
                            if (Options.nfacesTib==49)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x171_rad1.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x171_rad1_debug.dll');
                                end
                            elseif (Options.nfacesTib==49)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==0.5)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x171_rad05.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x171_rad05_debug.dll');
                                end
%%

                            elseif (Options.nfacesTib==49)&&(Options.nfacesFem==258)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x258_rad1.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x258_rad1_debug.dll');
                                end
                            elseif (Options.nfacesTib==49)&&(Options.nfacesFem==258)&&(Options.rad4Pairs==0.5)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x258_rad05.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x258_rad05_debug.dll');
                                end                             
%% 
                            elseif (Options.nfacesTib==75)&&(Options.nfacesFem==258)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x258_rad1.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x258_rad1_debug.dll');
                                end 
                            elseif (Options.nfacesTib==75)&&(Options.nfacesFem==258)&&(Options.rad4Pairs==0.5)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x258_rad05.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x258_rad05_debug.dll');
                                end
                            elseif (Options.nfacesTib==75)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==0.5)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x171_rad05.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x171_rad05_debug.dll');
                                end
                            elseif (Options.nfacesTib==75)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x171_rad1.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_75x171_rad1_debug.dll');
                                end
                            elseif (Options.nfacesTib==100)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==0.5)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x171_rad05.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x171_rad05_debug.dll');
                                end
                            elseif (Options.nfacesTib==100)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x171_rad1.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x171_rad1_debug.dll');
                                end
                            elseif (Options.nfacesTib==100)&&(Options.nfacesFem==188)&&(Options.rad4Pairs==1)
                                if (Options.kInpress==1e4)&&(Options.kInCheckContacts==1e3)
                                    if Options.KCFasinputstoExternalFunction
                                        F2_skeletal=external('F','TrackSim_2_kneeCont_KCFasinput.dll');
                                        F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress1e4_checkContact1e3_100x188_rad1.dll');
                                        F2_debug = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress1e4_checkContact1e3_100x188_rad1_forDebug.dll');
                                    end
                                end
                            elseif Options.numpairs==932 %not used anymore
                                if (Options.kInpress==5e4)&&(Options.kInCheckContacts==1e4)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_pairs932.dll');
                                    F2_debug=external('F','..\ExternalFunctions_GC\TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_pairs932_forDebug.dll');
                                elseif (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e4)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e4_pairs932.dll');
                                    F2_debug = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e4_pairs932_forDebug.dll');                  
                                elseif (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_pairs932.dll');
                                    F2_debug = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_pairs932_forDebug.dll');
                                elseif (Options.kInpress==5e4)&&(Options.kInCheckContacts==1e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_checkContact1e3_pairs932.dll');
                                    F2_debug= external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_checkContact1e3_pairs932_forDebug.dll');
                                elseif (Options.kInpress==5e4)&&(Options.kInCheckContacts==5e3)
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_checkContact5e3_pairs932.dll');
                                    F2_debug = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4_checkContact5e3_pairs932_forDebug.dll');
                                else
                                    keyboard;
                                end
                            elseif Options.numpairs==1378
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x188_pairs1378.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_49x188_pairs1378_forDebug.dll');
                                end
                            elseif Options.numpairs==1810
                                if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                                    F2=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x188_pairs1810.dll');
                                    F2_debug=external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e5_checkContact1e3_100x188_pairs1810_forDebug.dll');
                                end
                            else
                                %numpairs=499
                                F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress5e4.dll');
                            end
                        end
                    elseif strcmp(Options.maxsmoothness,'nosmooth')
                        if Options.kInpress==5e5  
                            F2 = external('F','TrackSim_2_kneeCont_hardmax_kpress5e5.dll');
                        else
                            F2 = external('F','TrackSim_2_kneeCont_maxnosmooth.dll');
                        end
                    elseif strcmp(Options.maxsmoothness,'pnorm')
                        if  Options.kInmaxpen==1e4
                            F2 = external('F','TrackSim_2_kneeCont_pnorm.dll');
                        end
                    else
                            keyboard; %no function found
                    end
                end
            elseif hessi == 2
                disp('Memory issue with exact Hessian; case not available')
            end
        case 'AD_ADOLC' 
            disp('ADOL-C cases not available');
        case 'FD'
            %deprecated
            F1 = external('F','TrackSim_1_kneeCont.dll',struct('enable_fd',true,...
                'enable_forward',false,'enable_reverse',false,...
                'enable_jacobian',false,'fd_method','forward'));
            F2 = external('F','TrackSim_2_kneeCont.dll',struct('enable_fd',true,...
                'enable_forward',false,'enable_reverse',false,...
                'enable_jacobian',false,'fd_method','forward'));        
    end
    cd(pathmain);
    
    
    %% Indices external function
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
    % External function: F2
    % First, joint torques. 
    jointi.pelvis.tilt  = 1; 
    jointi.pelvis.list  = 2; 
    jointi.pelvis.rot   = 3; 
    jointi.pelvis.tx    = 4;
    jointi.pelvis.ty    = 5;
    jointi.pelvis.tz    = 6;
    jointi.hip_flex.l   = 7;
    jointi.hip_add.l    = 8;
    jointi.hip_rot.l    = 9;
    jointi.hip_flex.r   = 10;
    jointi.hip_add.r    = 11;
    jointi.hip_rot.r    = 12;
    jointi.knee_flex.l  = 13;
    jointi.knee_flex.r  = 14;
    jointi.knee_add.r   = 15;
    jointi.knee_rot.r   = 16;
    jointi.knee_tx.r    = 17;
    jointi.knee_ty.r    = 18;
    jointi.knee_tz.r    = 19;
    jointi.ankle.l      = 20;
    jointi.ankle.r      = 21;
    jointi.subt.l       = 22;
    jointi.subt.r       = 23;
    jointi.trunk.ext    = 24;
    jointi.trunk.ben    = 25;
    jointi.trunk.rot    = 26;
    jointi.sh_flex.l    = 27;
    jointi.sh_add.l     = 28;
    jointi.sh_rot.l     = 29;
    jointi.sh_flex.r    = 30;
    jointi.sh_add.r     = 31;
    jointi.sh_rot.r     = 32;
    jointi.elb.l        = 33;
    jointi.elb.r        = 34;
    % jointi.prosup.l     = 35; %fixed values
    % jointi.prosup.r     = 36;
    
    % Vectors of indices for later use
    residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all 
    ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
    trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
    armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
    residuals_acti      = [jointi.hip_flex.l:jointi.knee_flex.r jointi.ankle.r:jointi.elb.r]; % all but gr-pelvis and knee
    residual_bptyi      = [jointi.pelvis.tilt:jointi.pelvis.tx,...
        jointi.pelvis.tz:jointi.elb.r]; % all but pelvis_ty
    knee_secDOF         = jointi.knee_add.r:jointi.knee_tz.r;
    % Number of degrees of freedom for later use
    nq.all              = length(residualsi); % all 
    nq.abs              = length(ground_pelvisi); % ground-pelvis
    nq.act              = nq.all-nq.abs;% all but ground-pelvis
    nq.trunk            = length(trunki); % trunk
    nq.arms             = length(armsi); % arms
    nq.leg              = 9; % #joints needed for polynomials
    nq.legright         = 14;
    nq.legleft          = 9;
    Qsi                 = 1:2:2*nq.all; % indices Qs only
    % Second, GRFs
    GRFi.r              = 35:37;
    GRFi.l              = 38:40;
    GRFi.all            = [GRFi.r,GRFi.l];
    nGRF                = length(GRFi.all);
    % Third, GRMs
    GRMi.r              = 41:43;
    GRMi.l              = 44:46;
    GRMi.all            = [GRMi.r,GRMi.l];
    nGRM                = length(GRMi.all);
    KCFi.L              = 47;
    KCFi.M              = 48;
    % Number contact model parameters
    np = 18;  
    
    %% Model info
    body_mass = 62.6;
    body_weight = body_mass*9.81;
    
    %% Collocation scheme
    % We use a pseudospectral direct collocation method, i.e. we use Lagrange
    % polynomials to approximate the state derivatives at the collocation
    % points in each mesh interval. We use d=3 collocation points per mesh
    % interval and Radau collocation points. 
    pathCollocationScheme = [pathRepo,'/CollocationScheme'];
    addpath(genpath(pathCollocationScheme));
    d = 3; % degree of interpolating polynomial
    method = 'radau'; % collocation method
    [tau_root,C,D,B] = CollocationScheme(d,method);
    
    %% Muscle-tendon parameters 
    % Muscles from one leg and from the back
    muscleNames = {'addbrev','addlong','addmagProx','addmagMid',...
            'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
            'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
            'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
            'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
            'recfem','sart','semimem','semiten','soleus','tfl',...
            'tibant','tibpost','vasint','vaslat',...
            'vasmed','ercspn_r','intobl_r','extobl_r',...
            'ercspn_l','intobl_l','extobl_l'};
    % Muscle indices for later use
    pathmusclemodel = [pathRepo,'/MuscleModel'];
    addpath(genpath(pathmusclemodel));    
    % (1:end-3), since we do not want to count twice the back muscles
    musi = MuscleIndices_3D_GC(muscleNames(1:end-3));
    % Total number of muscles
    NMuscle = length(musi)*2;
    % Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
    % fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
    % angles; Row 5: maximal contraction velocities
    load([pathmusclemodel,'/MTparameters_',subject,'.mat']);
    MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
    %Slightly modify lM0 to avoid huge passive forces at the initial guess
    MTparameters_m(2,56)=0.09; %edl l lM0;
    MTparameters_m(2,59)=0.08; %fhl l lM0;
    MTparameters_m(2,60)=0.1; %gaslat l lM0
    MTparameters_m(2,61)=0.1; %gasmed l lM0
    MTparameters_m(2,62)=0.04; %gem l lM0
    MTparameters_m(2,75)=0.07; %perbrev l lM0
    MTparameters_m(2,76)=0.1; %perlong l lM0
    MTparameters_m(2,77)=0.09; %pertert l lM0
    MTparameters_m(2,78)=0.04; %piri l lM0
    MTparameters_m(2,81)=0.1; %recfem l lM0
    MTparameters_m(2,85)=0.1; %soleus l lM0
    MTparameters_m(2,87)=0.1; %tibant l lM0
    MTparameters_m(2,88)=0.06; %tibpost l lM0
    MTparameters_m(2,90)=0.105; %vaslat l lM0
    
    MTparameters_m(2,9)=0.09; %edl r lM0;
    MTparameters_m(2,12)=0.08; %fhl r lM0;
    MTparameters_m(2,13)=0.1; %gaslat r lM0
    MTparameters_m(2,14)=0.1; %gasmed r lM0
    MTparameters_m(2,15)=0.04; %gem r lM0
    MTparameters_m(2,28)=0.07; %perbrev r lM0
    MTparameters_m(2,29)=0.1; %perlong r lM0
    MTparameters_m(2,30)=0.09; %pertert r lM0
    MTparameters_m(2,31)=0.04; %piri r lM0
    MTparameters_m(2,34)=0.1; %recfem r lM0
    MTparameters_m(2,38)=0.1; %soleus r lM0
    MTparameters_m(2,40)=0.1; %tibant r lM0
    MTparameters_m(2,41)=0.06; %tibpost r lM0
    MTparameters_m(2,43)=0.105; %vaslat r lM0
    
    %These are just to check if semimem can operate closer to lMtilde=1
    MTparameters_m(2,36)=0.06; %semimem_l
    MTparameters_m(2,36+47)=0.06; %semimem_r
    MTparameters_m(3,36)=0.345; %semimem_l
    MTparameters_m(3,36+47)=0.345; %semimem_r
    
    MTparameters_m(1,:)=MTparameters_m(1,:)*1.5;
    
    % Indices of the muscles actuating the different joints for later use
    pathpolynomial = [pathRepo,'/Polynomials_GC'];
    addpath(genpath(pathpolynomial));
    tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
    [~,mai] = MomentArmIndices_3D(muscleNames(1:end-3),...
        tl.muscle_spanning_joint_INFO(1:end-3,:));
    
    %% Contact model parameters
    % We optimize the locations of the contact spheres (x and z coordinates) 
    % and the radii. The other parameters are fixed.
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
    
    %% CasADi functions
    % We create several CasADi functions for later use
    pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
    addpath(genpath(pathCasADiFunctions));
    pathContactModel = [pathRepo,'/Contact'];
    addpath(genpath(pathContactModel));
    % We load some variables for the polynomial approximations
    load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
    if Options.useReducedPolynomials
        if Options.useReoptimizedPoly
            switch Options.err_poly 
                case 0.5
                    load([pathpolynomial,'/MuscleInfo_red_05_ropt','.mat']);
                    %MuscleInfo = MuscleInfo_red_05;
                case 1
                    load([pathpolynomial,'/MuscleInfo_red_1_ropt','.mat']);
                    %MuscleInfo = MuscleInfo_red_1;
                case 2
                    load([pathpolynomial,'/MuscleInfo_red_2_ropt','.mat']);
                    %MuscleInfo = MuscleInfo_red_2;
                case 3
                    load([pathpolynomial,'/MuscleInfo_red_3_ropt','.mat']);
                    %MuscleInfo = MuscleInfo_red_3;
            end
        else
            switch Options.err_poly 
                case 0.5
                    load([pathpolynomial,'/MuscleInfo_red_05','.mat']);
                    %MuscleInfo = MuscleInfo_red_05;
                case 1
                    load([pathpolynomial,'/MuscleInfo_red_1','.mat']);
                    %MuscleInfo = MuscleInfo_red_1;
                case 2
                    load([pathpolynomial,'/MuscleInfo_red_2','.mat']);
                    %MuscleInfo = MuscleInfo_red_2;
                case 3
                    load([pathpolynomial,'/MuscleInfo_red_3','.mat']);
                    %MuscleInfo = MuscleInfo_red_3;
            end
        end
    else
        switch Options.err_poly 
            case 0.5
                load([pathpolynomial,'/MuscleInfo_full_05','.mat']);
                %MuscleInfo = MuscleInfo_full_05;
            case 1
                load([pathpolynomial,'/MuscleInfo_full_1','.mat']);
                %MuscleInfo = MuscleInfo_full_1;
            case 2
                load([pathpolynomial,'/MuscleInfo_full_2','.mat']);
                %MuscleInfo = MuscleInfo_full_2;
            case 3
                load([pathpolynomial,'/MuscleInfo_full_3','.mat']);
                %MuscleInfo = MuscleInfo_full_3;
    
        end
    end
    
    % For the polynomials, we want all independent muscles. So we do not need
    % the muscles from both legs, since we assume bilateral symmetry, but want
    % all muscles from the back (indices 48:50).
    musi_pol = [musi,48,49,50];
    NMuscle_pol = NMuscle/2+3;
    CasADiFunctions_3D_GC
    cd(pathExternalFunctions);
    if Options.usedCompiled_lMT_vMT_External
        if Options.useReducedPolynomials
            if Options.useReoptimizedPoly
                switch Options.err_poly 
                    case 0.5 
                        keyboard;
                    case 1
                        f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_dM_red_1.dll');
                    case 2
                        f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_dM_red_2.dll');
                end
            else
                keyboard;
            end
        else
            switch Options.err_poly 
                case 0.5
                    keyboard;
                case 1
                    f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_dM_full_1.dll');
                case 2
                    f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_dM_full_2.dll');
            end
        end
    end
    cd(pathmain);
    
    %% Passive joint torques
    % We extract the parameters for the passive torques of the lower limbs and
    % the trunk
    pathPassiveMoments = [pathRepo,'/PassiveMoments'];
    addpath(genpath(pathPassiveMoments));
    PassiveMomentsData
    
    %% Experimental data
    pathData = [pathRepo,'/OpenSimModel/',subject];
    joints = {'pelvis_rz','pelvis_rx','pelvis_ry','pelvis_tx',...
        'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
        'hip_rotation_l','hip_flexion','hip_adduction','hip_rotation',...
        'knee_flexion_l','knee_flexion','knee_adduction','knee_rotation',...
        'knee_tx','knee_ty','knee_tz','ankle_angle_l','ankle_angle',...
        'subtalar_angle_l','subtalar_angle',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex_l','arm_add_l','arm_rot_l',...
        'arm_flex_r','arm_add_r','arm_rot_r',...
        'elbow_flex_l','elbow_flex_r'};
    pathVariousFunctions = [pathRepo,'/VariousFunctions'];
    addpath(genpath(pathVariousFunctions));
    % Extract joint kinematics
    pathIK = [pathData,'/IK/',nametrial.IK,'.mat'];
    Qs = getIK_GC(pathIK,joints);
    % Extract ground reaction forces and moments
    pathGRF = [pathData,'/GRF/',nametrial.GRF,'.mat'];
    GRF = getGRF_GC(pathGRF);
    % Extract joint kinetics
    pathID = [pathData,'/ID/',nametrial.ID,'.mat'];
    ID = getID_GC(pathID,joints);
    %Extract knee contact forces
    pathKCF=[pathData, '/KCF/',nametrial.KCF,'.mat'];
    KCF=get_KCF(pathKCF);
    % Interpolation experimental data
    sampfreq_kin=1/120;
    time_expi.ID(1) = find((ID.time<(time_opt(1)+sampfreq_kin/2))&(ID.time>=(time_opt(1)-sampfreq_kin/2)));
    time_expi.ID(2) = find((ID.time<(time_opt(2)+sampfreq_kin/2))&(ID.time>=(time_opt(2)-sampfreq_kin/2)));
    time_expi.KCF(1) = find((KCF.time<(time_opt(1)+sampfreq_kin/2))&(KCF.time>=(time_opt(1)-sampfreq_kin/2)));
    time_expi.KCF(2) = find((KCF.time<(time_opt(2)+sampfreq_kin/2))&(KCF.time>=(time_opt(2)-sampfreq_kin/2)));
    sampfreq_GRF=1/1200;
    time_expi.GRF(1) = find((GRF.time<(time_opt(1)+sampfreq_GRF/2))&(GRF.time>=(time_opt(1)-sampfreq_GRF/2)));
    time_expi.GRF(2) = find((GRF.time<(time_opt(2)+sampfreq_GRF/2))&(GRF.time>=(time_opt(2)-sampfreq_GRF/2)));
    % step = (ID.time(time_expi.ID(2))-ID.time(time_expi.ID(1)))/N;
    step=(time_opt(2)-time_opt(1))/N;
    interval = time_opt(1):step:time_opt(2);
    tgrid_col=zeros(N*d,1);
    tgrid_col(1:d:end)=interval(1:end-1)+tau_root(2)*step;
    tgrid_col(2:d:end)=interval(1:end-1)+tau_root(3)*step;
    tgrid_col(3:d:end)=interval(1:end-1)+tau_root(4)*step;
    ID.allinterp = interp1(ID.all(:,1),ID.all,interval);
    ID.allinterp_col = interp1(ID.all(:,1),ID.all,tgrid_col);
    Qs.allinterpfilt = interp1(Qs.allfilt(:,1),Qs.allfilt,interval);
    Qs.allinterpfilt_col = interp1(Qs.allfilt(:,1),Qs.allfilt,tgrid_col');
    KCF.allinterpfilt = interp1(KCF.data(:,1),KCF.data(:,2:3),interval);
    KCF.allinterpfilt_col=interp1(interval,KCF.allinterpfilt,tgrid_col);
    GRF.val.allinterp = interp1(round(GRF.val.all(:,1),4),...
        GRF.val.all,round(interval,4));
    GRF.val.allinterp_col=interp1(interval,GRF.val.allinterp,tgrid_col);
    GRF.MorGF.allinterp = interp1(round(GRF.MorGF.all(:,1),4),...
        GRF.MorGF.all,round(interval,4));
    GRF.MorGF.allinterp_col=interp1(interval,GRF.MorGF.allinterp,tgrid_col);
    GRF.pos.allinterp = interp1(round(GRF.pos.all(:,1),4),...
        GRF.pos.all,round(interval,4));
    GRF.pos.allinterp_col=interp1(interval,GRF.pos.allinterp,tgrid_col);
    GRF.Mcop.allinterp = interp1(round(GRF.Mcop.all(:,1),4),...
        GRF.Mcop.all,round(interval,4));
    GRF.Mcop.allinterp_col=interp1(interval,GRF.Mcop.allinterp,tgrid_col);
    
    %% Bounds
    pathBounds = [pathRepo,'/Bounds'];
    addpath(genpath(pathBounds));
    [bounds,scaling] = getBounds_3D_GC(Qs,NMuscle,nq,jointi,dev_cm,GRF,'');
    
    %% Scale Qs experimental 
    Qs_res_interpfilt_scaled(:,jointi.knee_ty.r)= (Qs.allinterpfilt(:,jointi.knee_ty.r+1)-scaling.knee_ty.a)/scaling.knee_ty.b;
    Qs_res_interpfilt_scaled(:,[1:jointi.knee_tx.r jointi.knee_tz.r:nq.all])=    Qs.allinterpfilt(:, [1:jointi.knee_tx.r jointi.knee_tz.r:nq.all]+1)./scaling.Qs([1:jointi.knee_tx.r jointi.knee_tz.r:end]);          
    Qs_res_interpfilt_scaled_aux=Qs_res_interpfilt_scaled(:,[1:4 6:14 20:34]); %exclude knee internal dofs
    Qs_res_interpfilt_scaled_aux_col=interp1(interval(1:N+1),Qs_res_interpfilt_scaled_aux,tgrid_col);
    
    for i=2:size(Qs.allinterpfilt,2)
        Q_spline(i-1)=spline(Qs.allinterpfilt(:,1),Qs.allinterpfilt(:,i));
        Qdot_spline(i-1)=fnder(Q_spline(i-1),1);
        Qdots.allinterpfilt(:,1)=Qs.allinterpfilt(:,1);
        Qdots.allinterpfilt(:,i)=ppval(Qdot_spline(i-1),Qdots.allinterpfilt(:,1));
    end
    Qdots_res_interpfilt_scaled=Qdots.allinterpfilt(:,2:end)./scaling.Qdots;
    Qdots_res_interpfilt_scaled_aux=Qdots_res_interpfilt_scaled(:,[1:4 6:14 20:34]); %exclude knee internal dofs
    Qdots_res_interpfilt_scaled_aux_col=interp1(interval(1:N+1),Qdots_res_interpfilt_scaled_aux,tgrid_col);
    
    
    %% Initial guess
    pathIG = [pathRepo,'/Guess'];
    addpath(genpath(pathIG));
    if (Options.multiple_IG==0)
    else
        if ~exist([pathIG '\multipleIGs_v2', savename2, '.mat'], 'file')&(Options.multiple_IG==1)
            Options.force2createNewIGs=1; 
        else
            load([pathIG '\multipleIGs_v2', savename2, '.mat']);
            if ~isfield(IGs,nametrial.id)
                Options.force2createNewIGs=1; 
            end
        end
    end
    if (Options.multiple_IG==0)|(Options.force2createNewIGs==1)
        % The initial guess depends on the settings
        if IGi == 1 % Quasi-random initial guess  
            %not used for now
            keyboard;
            guess = getGuess_3D_QR(Qs,nq,N,NMuscle,jointi,time_opt,scaling);
        elseif IGi == 2 % Data-informed initial guess  
            if Options.KCFasinputstoExternalFunction
                guess = getGuess_3D_DI_GC_v2(Qs,nq,N,NMuscle,jointi,scaling,d,KCF,GRF,F1,{F2_skeletal F2},f_contactForce,Options,'track');
            else
                guess = getGuess_3D_DI_GC_v2(Qs,nq,N,NMuscle,jointi,scaling,d,KCF,GRF,F1,F2,f_contactForce,Options,'track');
            end
        elseif IGi == 3 % Data-informed initial guess with muscle information obtained 
            %not used for now
            keyboard;
            % using the muscle redundancy solver 
            % (https://simtk.org/projects/optcntrlmuscle)
            load([pathIG,'/mvar.mat'],'mvar');
            mvarinterp.a = interp1(mvar.t(:,1),mvar.a,interval);
            mvarinterp.FTtilde = interp1(mvar.t(:,1),mvar.FTtilde,interval);
            mvarinterp.vA = interp1(mvar.t(1:end-1,1),mvar.vA,interval);
            mvarinterp.dFTtilde = interp1(mvar.t(1:end-1,1),mvar.dFTtilde,interval);
            guess = getGuess_3D_DIm(Qs,mvarinterp,nq,N,jointi,scaling);
        end
        
    end
    
    if Options.force2createNewIGs==1
        guess_i(1)=guess;
        IGnoise=0.01*Options.noiseApplied2NewIGs;
        for i=1:9
            guess_i(i+1).Qs=guess.Qs+IGnoise*randn([N+1 nq.all]);
            guess_i(i+1).Qdots=guess.Qdots+IGnoise*randn([N+1 nq.all]);
            guess_i(i+1).QsQdots=guess.QsQdots+IGnoise*randn([N+1 nq.all*2]);
            guess_i(i+1).Qdotdots=guess.Qdotdots+IGnoise*randn([N+1 nq.all]);
            guess_i(i+1).a=guess.a+IGnoise*randn([N+1 NMuscle]);
            guess_i(i+1).a(guess_i(i+1).a<1e-4)=1e-4;
            guess_i(i+1).vA=guess.vA+IGnoise*randn([N NMuscle]);
            guess_i(i+1).FTtilde=guess.FTtilde+IGnoise*randn([N+1 NMuscle]);
            guess_i(i+1).FTtilde(guess_i(i+1).FTtilde<1e-4)=1e-4;
            guess_i(i+1).dFTtilde=guess.dFTtilde+IGnoise*randn([N NMuscle]);
            guess_i(i+1).a_a=guess.a_a+IGnoise*randn(size(guess.a_a));
            guess_i(i+1).e_a=guess.e_a+IGnoise*randn(size(guess.e_a));
            guess_i(i+1).params=guess.params+IGnoise*randn(size(guess.params));
            guess_i(i+1).KCFresiduals=guess.KCFresiduals+IGnoise*randn(size(guess.KCFresiduals));
            guess_i(i+1).JRM_res=guess.JRM_res+IGnoise*randn(size(guess.JRM_res));
            guess_i(i+1).a_col=guess.a_col+IGnoise*randn(size(guess.a_col));
            guess_i(i+1).a_col(guess_i(i+1).a_col<1e-4)=1e-4;
            guess_i(i+1).FTtilde_col=guess.FTtilde_col+IGnoise*randn(size(guess.FTtilde_col));
            guess_i(i+1).QsQdots_col=guess.QsQdots_col+IGnoise*randn(size(guess.QsQdots_col));
            guess_i(i+1).a_a_col=guess.a_a_col+IGnoise*randn(size(guess.a_a_col));
        end
        IGs.(nametrial.id).guess_i=guess_i;
        save([pathIG '\multipleIGs_v2',savename2, '.mat'],'IGs');
    end
    if Options.multiple_IG==1
        guess=IGs.(nametrial.id).guess_i(Options.IGn);
    end
    
    % This allows visualizing the initial guess and the bounds
    if checkBoundsIG
        pathPlots = [pathRepo,'/Plots'];
        addpath(genpath(pathPlots));
        plot_BoundsVSInitialGuess_3D_GC
    end
    
    %% Formulate the NLP
    if solveProblem
        % Start with an empty NLP
        % Initialize opti instance.
        % Opti is a collection of CasADi helper classes:
        % https://web.casadi.org/docs/#opti-stack
        opti = casadi.Opti();
        g_names = []; % vector containing names of the constraints. To delete?
        
        % Define static parameters
        % Contact model parameters   
        paramsCM=opti.variable(1,np);
        opti.subject_to(bounds.params.lower < paramsCM < bounds.params.upper);
        opti.set_initial(paramsCM,guess.params');  
            
        % Define states
        % Muscle activations at mesh points
        a = opti.variable(NMuscle,N+1);
        opti.subject_to(bounds.a.lower'*ones(1,N+1) < a < bounds.a.upper'*ones(1,N+1));
        opti.set_initial(a,guess.a');
        % Muscle activations at collocation points
        a_col = opti.variable(NMuscle,d*N);
        opti.subject_to(bounds.a.lower'*ones(1,d*N) < a_col < ...
            bounds.a.upper'*ones(1,d*N));
        opti.set_initial(a_col, guess.a_col');   
        % Muscle-tendon forces at mesh points
        FTtilde = opti.variable(NMuscle,N+1);
        opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1) < FTtilde < ...
            bounds.FTtilde.upper'*ones(1,N+1));
        opti.set_initial(FTtilde, guess.FTtilde');
        % Muscle-tendon forces at collocation points
        FTtilde_col = opti.variable(NMuscle,d*N);
        opti.subject_to(bounds.FTtilde.lower'*ones(1,d*N) < FTtilde_col < ...
            bounds.FTtilde.upper'*ones(1,d*N));
        opti.set_initial(FTtilde_col, guess.FTtilde_col');   
        % Qs and Qdots at mesh points
        X=opti.variable(2*nq.all,N+1);
        opti.subject_to(bounds.QsQdots.lower' < X < bounds.QsQdots.upper');
        opti.set_initial(X,guess.QsQdots');
        % Qs and Qdots at collocation points
        X_col = opti.variable(2*nq.all,d*N);
        opti.subject_to(bounds.QsQdots.lower' < X_col < bounds.QsQdots.upper');
        opti.set_initial(X_col,guess.QsQdots_col');
        % Arm activations at mesh points
        a_a = opti.variable(nq.arms,N);
        opti.subject_to(bounds.a_a.lower'*ones(1,N) < a_a < ...
            bounds.a_a.upper'*ones(1,N));
        opti.set_initial(a_a, guess.a_a');  
        % Arm activations at collocation points
        a_a_col = opti.variable(nq.arms,d*N);
        opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
            bounds.a_a.upper'*ones(1,d*N));
        opti.set_initial(a_a_col, guess.a_a_col');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define controls
        % Time derivative of muscle activations (states) at mesh points
        vA = opti.variable(NMuscle, N);
        opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
            bounds.vA.upper'*ones(1,N));
        opti.set_initial(vA, guess.vA');   
        % Arm excitations
        e_a = opti.variable(nq.arms, N);
        opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
            bounds.e_a.upper'*ones(1,N));
        opti.set_initial(e_a, guess.e_a');
        if Options.useKCFresiduals
    %         % KCF residuals at mesh points
    %         KCF_res=opti.variable(5,N);
    %         opti.subject_to(bounds.KCF_res.lower' < KCF_res < ...
    %             bounds.KCF_res.upper');
    %         opti.set_initial(KCF_res,guess.KCFresiduals'.*ones(1,N));
            % KCF residuals at collocation points
            KCF_res_col=opti.variable(5,d*N);
            opti.subject_to(bounds.KCF_res.lower'*ones(1,d*N) < KCF_res_col < ...
                bounds.KCF_res.upper'*ones(1,d*N));
            opti.set_initial(KCF_res_col,guess.KCFresiduals'.*ones(1,d*N));
        end
        if Options.useJointResMom
            %Joint residual moments at mesh points
    %         JointMom_res = opti.variable(21,N); %all dofs - 5 - arm torque driven dofs
    %         opti.subject_to(bounds.JRM_res.lower' < JointMom_res < ...
    %             bounds.JRM_res.upper');
    %         opti.set_initial(JointMom_res,guess.JRM_res'.*ones(1,N));
            %Joint residual moments at collocation points
            Jointmom_res_col=opti.variable(21,d*N); %all dofs - 5 - arm torque driven dofs
            opti.subject_to(bounds.JRM_res.lower'*ones(1,d*N) < Jointmom_res_col < ...
                bounds.JRM_res.upper'*ones(1,d*N));
        end
        if Options.usePelvisResMom
            Pelvis_res_col=opti.variable(6,d*N);
            opti.subject_to(bounds.Pelvis_res.lower'*ones(1,d*N)< Pelvis_res_col <...
                bounds.Pelvis_res.upper'*ones(1,d*N));
            opti.set_initial(Pelvis_res_col,zeros(6,N*d));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define "slack" controls
        % Time derivative of muscle-tendon forces (states) at collocation points
        dFTtilde_col = opti.variable(NMuscle, d*N);
        opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
            bounds.dFTtilde.upper'*ones(1,d*N));
        opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
        % Time derivative of Qdots (states) at collocation points
        A_col = opti.variable(nq.all, d*N);
        opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
            bounds.Qdotdots.upper'*ones(1,d*N));
        opti.set_initial(A_col, guess.Qdotdots_col');          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel formulation
        % Define CasADi variables for static parameters
        paramsCMk=MX.sym('paramsCM',np); 
        paramsCM_nsc= f_nsc18(paramsCMk,scaling.params.v,scaling.params.r); 
        % Define CasADi variables for states
        ak = MX.sym('ak',NMuscle);
        aj = MX.sym('akmesh',NMuscle,d);
        akj = [ak aj];
        FTtildek = MX.sym('FTtildek',NMuscle); 
        FTtildej = MX.sym('FTtildej',NMuscle,d);
        FTtildekj = [FTtildek FTtildej];
        Xk = MX.sym('Xk',2*nq.all);
        Xj = MX.sym('Xj',2*nq.all,d);
        Xkj = [Xk Xj];
        a_ak = MX.sym('a_ak',nq.arms);
        a_aj = MX.sym('a_akmesh',nq.arms,d);
        a_akj = [a_ak a_aj];
        % Define CasADi variables for controls
        vAk = MX.sym('vAk',NMuscle);
        e_ak = MX.sym('e_ak',nq.arms);
        if Options.useKCFresiduals
            KCF_resk=MX.sym('KCF_resk',5);
            KCF_resj=MX.sym('KCF_resj',5,d);
            KCF_reskj=[KCF_resk KCF_resj];
        else
            KCF_resk=zeros(5,1);
            KCF_resj=zeros(5,d);
            KCF_reskj=[KCF_resk KCF_resj];
        end
        if Options.useJointResMom   
    %         JointMom_resk=MX.sym('JointMom_resk',21);
            JointMom_resj=MX.sym('JointMom_resj',21,d);
    %         JointMom_reskj=[JointMom_resk JointMom_resj];
        end
        % if Options.usePelvisResMom
            Pelvis_res_j=MX.sym('Pelvis_res_j',6,d);
        % end
        % Define CasADi variables for "slack" controls
        dFTtildej = MX.sym('dFTtildej',NMuscle,d);
        Aj = MX.sym('Aj',nq.all,d);   
        %Experimental data to track
        Qs_toTrack_scaled_auxk=MX.sym('Qs_toTrack_scaled_auxk',size(Qs_res_interpfilt_scaled_aux,2)); %resembles Qs_res_interpfilt_scaled_auxk at collocation points
        Qs_toTrack_scaled_auxj=MX.sym('Qs_toTrack_scaled_auxj',size(Qs_res_interpfilt_scaled_aux,2),d); %resembles Qs_res_interpfilt_scaled_auxk at collocation points
        Qs_toTrack_scaled_auxkj=[Qs_toTrack_scaled_auxk Qs_toTrack_scaled_auxj];
        Qdots_toTrack_scaled_auxk=MX.sym('Qdots_toTrack_scaled_auxk',size(Qs_res_interpfilt_scaled_aux,2)); %resembles Qdots_res_interpfilt_scaled_auxk at collocation points
        Qdots_toTrack_scaled_auxj=MX.sym('Qdots_toTrack_scaled_auxj',size(Qs_res_interpfilt_scaled_aux,2),d); %resembles Qdots_res_interpfilt_scaled_auxk at collocation points
        Qdots_toTrack_scaled_auxkj=[Qdots_toTrack_scaled_auxk Qdots_toTrack_scaled_auxj];
        GRF_toTrackk=MX.sym('GRF_toTrackk',6); %resembles GRF.val.allinterp (except first column) at collocation points
        GRF_toTrackj=MX.sym('GRF_toTrackj',6,d); %resembles GRF.val.allinterp (except first column) at collocation points
        GRF_toTrackkj=[GRF_toTrackk GRF_toTrackj];
    %     GRM_toTrackk=MX.sym('GRM_toTrackk',6); %resembles GRF.MorGF.allinterp (except first column) at collocation points
        GRM_toTrackj=MX.sym('GRM_toTrackj',6,d); %resembles GRF.MorGF.allinterp (except first column) at collocation points
    %     GRM_toTrackkj=[GRM_toTrackk GRM_toTrackj];
    %     KCF_toTrackk=MX.sym('KCF_toTrackk',2); %resembles KCF.allinterpfilt at collocation points
        KCF_toTrackj=MX.sym('KCF_toTrackj',2,d); %resembles KCF.allinterpfilt at collocation points
    %     KCF_toTrackkj=[KCF_toTrackk KCF_toTrackj];
    %     ID_toTrackk=MX.sym('ID_toTrackk',length(residuals_acti)); %resembles ID.allinterp(k+1,[(2+nq.abs):(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])
        ID_toTrackj=MX.sym('ID_toTrackj',length(residuals_acti)+6,d); %resembles ID.allinterp(k+1,[(2+nq.abs):(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])
    %     ID_toTrackkj=[ID_toTrackk ID_toTrackj];
        J = 0; % Initialize cost function
        J1=0;
        J1b=0;
        J2=0;
        J3=0;
        J4=0;
        J5=0;
        J6=0;
        J7=0;
        J7b=0;
        J8=0;
        J9=0;
        J10=0;
        J11=0;
        J12=0;
        
        eq_constr = {}; % Initialize equality constraint vector
        ineq_constr1 = {}; % Initialize inequality constraint vector 1
        ineq_constr2 = {}; % Initialize inequality constraint vector 2
        ineq_constr3 = {}; % Initialize inequality constraint vector 3
        ineq_constr4 = {}; % Initialize inequality constraint vector 4
        g_names_coll = {}; % Initialize names of constraints at collocation points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Time step
        h = (time_opt(2)-time_opt(1))/N;
        % Loop over collocation points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unscale variables        
        Xkj_nsc = Xkj.*(scaling.QsQdots'*ones(1,size(Xkj,2)));
        Xkj_nsc(jointi.knee_ty.r*2-1,:) = Xkj(jointi.knee_ty.r*2-1,:).*scaling.knee_ty.b+scaling.knee_ty.a;
        FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
        dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
        Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));  
        vAk_nsc = vAk.*scaling.vA;  
        %just for debugging 
        outk1j=MX(84,d);
        for j=1:d            
            % Call external functions
            % The first external function (F1) returns linear and angular
            % velocities of the calcaneus and toes, positions of the origin of 
            % the calcaneus and toes, and transforms from the calcaneus and 
            % toes to the ground. These variables are used to compute the 
            % contact forces that will be input of the second external function 
            % (F2).
            [outk1] = F1(Xkj_nsc(:,j+1)); 
            %just for debugging
            outk1j(:,j)=outk1;
            % Organize the locations of the contact spheres
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
                paramsCM_nsc(radi.s1:radi.s6)'];
            % The second external function (F2) returns joint torques, GRFs,
            % and GRMs based on joint states, controls, contact forces, and
            % several parameters of the contact models.
            if Options.KCFasinputstoExternalFunction
                q_knee=Xkj_nsc([jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r...
                    jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r]*2-1,j+1);
                out_skeletal=F2(q_knee);
                KCF_Mj=out_skeletal(1);
                KCF_Lj=out_skeletal(2);
                SumForces=out_skeletal(3:5);
                SumMoments=out_skeletal(6:8);
                [Tj] = F2_skeletal([Xkj_nsc(:,j+1);Aj_nsc(:,j);in_F2';SumForces;SumMoments]);
            else
                [Tj] = F2([Xkj_nsc(:,j+1);Aj_nsc(:,j);in_F2']);
                KCF_Mj=Tj(KCFi.M);
                KCF_Lj=Tj(KCFi.L);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Get muscle-tendon lengths, velocities, and moment arms
            % Left leg
            qinj_l = [Xkj_nsc(jointi.hip_flex.l*2-1,j+1),...
                Xkj_nsc(jointi.hip_add.l*2-1,j+1), ...
                Xkj_nsc(jointi.hip_rot.l*2-1,j+1), ...
                0, 0, ...
                Xkj_nsc(jointi.knee_flex.l*2-1,j+1), ...
                0, 0.042, 0, ...
                Xkj_nsc(jointi.ankle.l*2-1,j+1),...
                Xkj_nsc(jointi.subt.l*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ext*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ben*2-1,j+1),...
                Xkj_nsc(jointi.trunk.rot*2-1,j+1)];  
            qdotinj_l = [Xkj_nsc(jointi.hip_flex.l*2,j+1),...
                Xkj_nsc(jointi.hip_add.l*2,j+1),...
                Xkj_nsc(jointi.hip_rot.l*2,j+1),...
                0,0,...
                Xkj_nsc(jointi.knee_flex.l*2,j+1),...
                0,0,0,...
                Xkj_nsc(jointi.ankle.l*2,j+1),...
                Xkj_nsc(jointi.subt.l*2,j+1),...
                Xkj_nsc(jointi.trunk.ext*2,j+1),...
                Xkj_nsc(jointi.trunk.ben*2,j+1),...
                Xkj_nsc(jointi.trunk.rot*2,j+1)]; 
            [lMTj_l,vMTj_l,MAj_l] =  f_lMT_vMT_dM(qinj_l,qdotinj_l);    
            MAj.hip_flex.l   =  MAj_l(mai(1).mus.l',1);
            MAj.hip_add.l    =  MAj_l(mai(2).mus.l',2);
            MAj.hip_rot.l    =  MAj_l(mai(3).mus.l',3);
            MAj.knee_flex.l  =  MAj_l(mai(6).mus.l',6);
            MAj.ankle.l      =  MAj_l(mai(10).mus.l',10);  
            MAj.subt.l       =  MAj_l(mai(11).mus.l',11); 
            % For the back muscles, we want left and right together: left
            % first, right second. In MuscleInfo, we first have the right
            % muscles (44:46) and then the left muscles (47:49). Since the back
            % muscles only depend on back dofs, we do not care if we extract
            % them "from the left or right leg" so here we just picked left.
            MAj.trunk_ext    =  MAj_l([48:50,mai(12).mus.l]',12);
            MAj.trunk_ben    =  MAj_l([48:50,mai(13).mus.l]',13);
            MAj.trunk_rot    =  MAj_l([48:50,mai(14).mus.l]',14);
            % Right leg
            qinj_r = [Xkj_nsc(jointi.hip_flex.r*2-1,j+1),...
                Xkj_nsc(jointi.hip_add.r*2-1,j+1),...
                Xkj_nsc(jointi.hip_rot.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_add.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_rot.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_flex.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_tx.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_ty.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_tz.r*2-1,j+1),...
                Xkj_nsc(jointi.ankle.r*2-1,j+1),...
                Xkj_nsc(jointi.subt.r*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ext*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ben*2-1,j+1),...
                Xkj_nsc(jointi.trunk.rot*2-1,j+1)];  
            qdotinj_r = [Xkj_nsc(jointi.hip_flex.r*2,j+1),...
                Xkj_nsc(jointi.hip_add.r*2,j+1),...
                Xkj_nsc(jointi.hip_rot.r*2,j+1),...
                Xkj_nsc(jointi.knee_add.r*2,j+1),...
                Xkj_nsc(jointi.knee_rot.r*2,j+1),...
                Xkj_nsc(jointi.knee_flex.r*2,j+1),...
                Xkj_nsc(jointi.knee_tx.r*2,j+1),...
                Xkj_nsc(jointi.knee_ty.r*2,j+1),...
                Xkj_nsc(jointi.knee_tz.r*2,j+1),...
                Xkj_nsc(jointi.ankle.r*2,j+1),...
                Xkj_nsc(jointi.subt.r*2,j+1),...
                Xkj_nsc(jointi.trunk.ext*2,j+1),...
                Xkj_nsc(jointi.trunk.ben*2,j+1),...
                Xkj_nsc(jointi.trunk.rot*2,j+1)];      
            [lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qinj_r,qdotinj_r);
             % Here we take the indices from left since the vector is 1:49
            MAj.hip_flex.r   =  MAj_r(mai(1).mus.l',1);
            MAj.hip_add.r    =  MAj_r(mai(2).mus.l',2);
            MAj.hip_rot.r    =  MAj_r(mai(3).mus.l',3);
            MAj.knee_add.r   =  MAj_r(mai(4).mus.l',4);
            MAj.knee_rot.r   =  MAj_r(mai(5).mus.l',5);
            MAj.knee_flex.r  =  MAj_r(mai(6).mus.l',6);
            MAj.knee_tx.r    =  MAj_r(mai(7).mus.l',7);
            MAj.knee_ty.r    =  MAj_r(mai(8).mus.l',8);
            MAj.knee_tz.r    =  MAj_r(mai(9).mus.l',9);
            MAj.ankle.r      =  MAj_r(mai(10).mus.l',10);
            MAj.subt.r       =  MAj_r(mai(11).mus.l',11);
             % Both legs
            % In MuscleInfo, we first have the right back muscles (45:47) and 
            % then the left back muscles (48:50). Here we re-organize so that
            % we have first the left muscles and then the right muscles.
            lMTj_lr     = [lMTj_l([1:44,48:50],1);lMTj_r(1:47,1)];
            vMTj_lr     = [vMTj_l([1:44,48:50],1);vMTj_r(1:47,1)];   
            % Get muscle-tendon forces and derive Hill-equilibrium       
            [Hilldiffj,FTj,~,~,~] =  f_forceEquilibrium_FtildeState(...
                    akj(:,j+1),FTtildekj_nsc(:,j+1),...
                    dFTtildej_nsc(:,j),lMTj_lr,vMTj_lr);
            % Get passive torques
            Tau_passj.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Xkj_nsc(jointi.hip_flex.l*2-1,j+1),...
                Xkj_nsc(jointi.hip_flex.l*2,j+1));
            Tau_passj.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Xkj_nsc(jointi.hip_flex.r*2-1,j+1),...
                Xkj_nsc(jointi.hip_flex.r*2,j+1));
            Tau_passj.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Xkj_nsc(jointi.hip_add.l*2-1,j+1),...
                Xkj_nsc(jointi.hip_add.l*2,j+1));
            Tau_passj.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Xkj_nsc(jointi.hip_add.r*2-1,j+1),...
                Xkj_nsc(jointi.hip_add.r*2,j+1));
            Tau_passj.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Xkj_nsc(jointi.hip_rot.l*2-1,j+1),...
                Xkj_nsc(jointi.hip_rot.l*2,j+1));
            Tau_passj.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Xkj_nsc(jointi.hip_rot.r*2-1,j+1),...
                Xkj_nsc(jointi.hip_rot.r*2,j+1)); 
            if Options.dampingInKneeSec  
                Tau_passj.knee_add.r    = f_passiveMoments_kneeintmom(Xkj_nsc(jointi.knee_add.r*2-1,j+1),Xkj_nsc(jointi.knee_add.r*2,j+1));
                Tau_passj.knee_rot.r    = f_passiveMoments_kneeintmom(Xkj_nsc(jointi.knee_rot.r*2-1,j+1),Xkj_nsc(jointi.knee_rot.r*2,j+1));
            else
                Tau_passj.knee_add.r    = f_passiveMoments_kneeintmom(Xkj_nsc(jointi.knee_add.r*2-1,j+1));
                Tau_passj.knee_rot.r    = f_passiveMoments_kneeintmom(Xkj_nsc(jointi.knee_rot.r*2-1,j+1));
            end
            Tau_passj.knee_flex.l   = f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Xkj_nsc(jointi.knee_flex.l*2-1,j+1),...
                Xkj_nsc(jointi.knee_flex.l*2,j+1));
            if Options.dampingInKneeSec
                Tau_passj.knee_tx.r     = f_passiveForce_kneeintf(Xkj_nsc(jointi.knee_tx.r*2-1,j+1),Xkj_nsc(jointi.knee_tx.r*2,j+1));
                Tau_passj.knee_ty.r     = f_passiveForce_kneeintf(                                0,Xkj_nsc(jointi.knee_ty.r*2,j+1));
                Tau_passj.knee_tz.r     = f_passiveForce_kneeintf(Xkj_nsc(jointi.knee_tz.r*2-1,j+1),Xkj_nsc(jointi.knee_tz.r*2,j+1));
            else
                Tau_passj.knee_tx.r     = f_passiveForce_kneeintf(Xkj_nsc(jointi.knee_tx.r*2-1,j+1));
                Tau_passj.knee_ty.r     = 0;
                Tau_passj.knee_tz.r     = f_passiveForce_kneeintf(Xkj_nsc(jointi.knee_tz.r*2-1,j+1));
            end
            Tau_passj.knee_flex.r        = f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Xkj_nsc(jointi.knee_flex.r*2-1,j+1),...
                Xkj_nsc(jointi.knee_flex.r*2,j+1));
            Tau_passj.ankle.l       = f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Xkj_nsc(jointi.ankle.l*2-1,j+1),...
                Xkj_nsc(jointi.ankle.l*2,j+1));
            Tau_passj.ankle.r       = f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Xkj_nsc(jointi.ankle.r*2-1,j+1),...
                Xkj_nsc(jointi.ankle.r*2,j+1));        
            Tau_passj.subt.l       = f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Xkj_nsc(jointi.subt.l*2-1,j+1),...
                Xkj_nsc(jointi.subt.l*2,j+1));
            Tau_passj.subt.r       = f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Xkj_nsc(jointi.subt.r*2-1,j+1),...
                Xkj_nsc(jointi.subt.r*2,j+1));        
            Tau_passj.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
                theta.pass.trunk.ext,Xkj_nsc(jointi.trunk.ext*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ext*2,j+1));
            Tau_passj.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
                theta.pass.trunk.ben,Xkj_nsc(jointi.trunk.ben*2-1,j+1),...
                Xkj_nsc(jointi.trunk.ben*2,j+1));
            Tau_passj.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
                theta.pass.trunk.rot,Xkj_nsc(jointi.trunk.rot*2-1,j+1),...
                Xkj_nsc(jointi.trunk.rot*2,j+1)); 
            Tau_passj_all = [Tau_passj.hip.flex.l,Tau_passj.hip.flex.r,...
                Tau_passj.hip.add.l,Tau_passj.hip.add.r,...
                Tau_passj.hip.rot.l,Tau_passj.hip.rot.r,...
                Tau_passj.knee_flex.l,Tau_passj.knee_flex.r,Tau_passj.ankle.l,...
                Tau_passj.ankle.r,Tau_passj.subt.l,Tau_passj.subt.r,...
                Tau_passj.trunk.ext,Tau_passj.trunk.ben,...
                Tau_passj.trunk.rot]';
            % Expression for the state derivatives at the collocation points
            Qsp_nsc          = Xkj_nsc(1:2:end,:)*C(:,j+1);
            Qdotsp_nsc       = Xkj_nsc(2:2:end,:)*C(:,j+1);            
            FTtildep_nsc    = FTtildekj_nsc*C(:,j+1);
            ap              = akj*C(:,j+1);
            a_ap            = a_akj*C(:,j+1);
            % Append collocation equations
            % Dynamic constraints are scaled using the same scale
            % factors as the ones used to scale the states
            % Activation dynamics (implicit formulation)
            eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
            g_names_coll = [g_names_coll; repmat({'vA_a_dyn'},NMuscle,1)];
            % Contraction dynamics (implicit formulation)     
            eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./...
                scaling.FTtilde';
            g_names_coll = [g_names_coll; repmat({'dFTtilde_FTtilde_dyn'},NMuscle,1)];
            % Skeleton dynamics (implicit formulation)               
            qdotj_nsc = Xkj_nsc(2:2:end,j+1); % velocity
            eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
            g_names_coll = [g_names_coll; repmat({'qdot_dyn'},length(qdotj_nsc),1)];
            eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./...
                scaling.QsQdots(2:2:end)';
            g_names_coll = [g_names_coll; repmat({'qdotdot_dyn'},length(Qdotsp_nsc),1)];
            % Arm activation dynamics (explicit formulation)   
            dadtj    = f_ArmActivationDynamics(e_ak,a_akj(:,j+1));
            eq_constr{end+1} = (h*dadtj - a_ap)./scaling.a_a;
            g_names_coll = [g_names_coll; repmat({'da_a_a_dyn'},length(a_ap),1)];

            % modify the hip fleion weight
            tracking_indices = [1:4 6:14 20:34];
            q_diff = Xkj(tracking_indices * 2 - 1, j+1) - Qs_toTrack_scaled_auxkj(:, j+1);
            weight_Qs = W.Qs * ones(length(tracking_indices), 1);
            weight_Qs(10) = 2;  % Right hip flexion gets special weight
            % Add contribution to quadrature function
            %J + ...
                   % W.Qs*B(j+1)*(f_J28(Xkj([1:4 6:14 20:34]*2-1,j+1)-...
                    %    Qs_toTrack_scaled_auxkj(:,j+1)))*h

            % Increase weight for secondary DoFs of the right knee
            % Initialize custom acceleration weights
            weight_qddot = 0.003 * ones(size(Aj,1),1);
            weight_qddot(jointi.knee_add.r) = 5;
            weight_qddot(jointi.knee_rot.r) = 5;
            weight_qddot(jointi.knee_tx.r)  = 5;
            weight_qddot(jointi.knee_ty.r)  = 5;
            weight_qddot(jointi.knee_tz.r)  = 5;
            Iallq_exceptkneesecdof=[1:14 20:34];
            Ikneesecdof=15:19;
            

            J = J + B(j+1) * (q_diff' * diag(weight_Qs) * q_diff) * h +...
                    0.1*W.Qs*B(j+1)*(sum(Xkj([jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r]*2-1,j+1).^2))*h+... %to minimize secondary knee dofs add rot, tx and tz
                    W.Qdots*B(j+1)*(f_J28(Xkj([1:4 6:14 20:34]*2,j+1)-...
                        Qdots_toTrack_scaled_auxkj(:,j+1)))*h + ...
                    W.GRF*B(j+1)*(f_J6((Tj(GRFi.all,1)./scaling.GRF')-...
                        GRF_toTrackj(:,j)./scaling.GRF'))*h +...
                    W.GRM*B(j+1)*(f_J6((Tj(GRMi.all,1)./scaling.GRM')-...
                        GRM_toTrackj(:,j)./scaling.GRM'))*h +...
                    W.KCF*B(j+1)*(f_J2(([KCF_Mj; KCF_Lj]./scaling.GRF(2))-...
                        KCF_toTrackj(:,j)./scaling.GRF(2)))*h + ...
                    W.ID_act*B(j+1)*(f_J22((Tj(residuals_acti,1)./scaling.T(1)')-...
                        ID_toTrackj(7:end,j)./scaling.T(1)))*h +...
                    W.a*B(j+1)*(f_J94(akj(:,j+1)))*h + ...
                    W.u_qd2dot*B(j+1)*(f_J29(Aj(Iallq_exceptkneesecdof,j)))*h +...
                    W.u_qd2dot_kneesecdof*B(j+1)*(f_J5(Aj(Ikneesecdof,j)))*h +...
                    W.u_vA*B(j+1)*(f_J94(vAk))*h +...
                    W.u*B(j+1)*(f_J94(dFTtildej(:,j)))*h;  %W.u_qd2dot*B(j+1)*(f_J34(Aj(:,j)))*h
                    % B(j+1) * (Aj(:,j)' * diag(weight_qddot) * Aj(:,j)) * h +...
            % J = J+100*B(j+1)*(Xkj(jointi.hip_flex.r*2-1,j+1)-Qs_toTrack_scaled_auxkj(I_hip_r,j+1))*h;
            if Options.useKCFresiduals
                    J=J+W.KCF_res*B(j+1)*(sum(KCF_resj([1:3 5],j).^2)+10*KCF_resj(4,j).^2)*h;
            end
            if Options.useJointResMom==1
                J=J+W.JRM_res*B(j+1)*(f_J21(JointMom_resj(:,j)))*h; %originally
    %                taken into account
            end
            if Options.usePelvisResMom
                J=J+W.minPelvisRes*B(j+1)*(sum((Pelvis_res_j(1:6,j)).^2))*h;
            else
                J=J+W.minPelvisRes*B(j+1)*(sum((Tj(1:6)/1000).^2))*h;
            end

            %J1 to J12 are for debugging purposes only
            J1=J1+B(j+1) * (q_diff' * diag(weight_Qs) * q_diff) * h + ...
                        0.1 * W.Qs * B(j+1) * (sum(Xkj([jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r]*2-1,j+1).^2)) * h;
            J1b=J1b+W.Qdots*B(j+1)*(f_J28(Xkj([1:4 6:14 20:34]*2,j+1)-...
                        Qdots_toTrack_scaled_auxkj(:,j+1)))*h+...
                        0.1*W.Qdots*B(j+1)*(sum(Xkj([jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r]*2,j+1).^2))*h;
            J2=J2+W.GRF*B(j+1)*(f_J6((Tj(GRFi.all,1)./scaling.GRF')-...
                        GRF_toTrackj(:,j)./scaling.GRF'))*h;
            J3=J3+W.GRM*B(j+1)*(f_J6((Tj(GRMi.all,1)./scaling.GRM')-...
                        GRM_toTrackj(:,j)./scaling.GRM'))*h;
            J4=J4+W.KCF*B(j+1)*(f_J2(([KCF_Mj; KCF_Lj]./scaling.GRF(2))-...
                        KCF_toTrackj(:,j)./scaling.GRF(2)))*h;
            J5=J5+W.ID_act*B(j+1)*(f_J22((Tj(residuals_acti,1)./scaling.T(1)')-...
                        ID_toTrackj(7:end,j)./scaling.T(1)))*h; %be careful since IK and ID indexes do not coincide
            J6=J6+W.a*B(j+1)*(f_J94(akj(:,j+1)))*h;
            J7=J7+W.u_qd2dot*B(j+1)*(f_J29(Aj(Iallq_exceptkneesecdof,j)))*h;
            J7b=J7b+W.u_qd2dot_kneesecdof*B(j+1)*(f_J5(Aj(Ikneesecdof,j)))*h;
            J8=J8+W.u_vA*B(j+1)*(f_J94(vAk))*h;
            J9=J9+W.u*B(j+1)*(f_J94(dFTtildej(:,j)))*h;
            J10=J10+W.KCF_res*B(j+1)*(f_J5(KCF_resj(:,j)))*h;
            if Options.useJointResMom == 1
                J11=J11+W.JRM_res*B(j+1)*(f_J21(JointMom_resk))*h; %ornally
            end 
            if Options.usePelvisResMom
                J12=J12+W.minPelvisRes*B(j+1)*(sum((Pelvis_res_j(1:6,j)).^2))*h;
            else
                J12=J12+W.minPelvisRes*B(j+1)*(sum((Tj(1:6)/1000).^2))*h;
            end
                
            % Add path constraints
            % Pelvis residuals (same as from inverse dynamics)
            % if Options.useJointResMom
            %     eq_constr{end+1} = (ID_toTrackj(1:6,j)-Tj(ground_pelvisi,1)-...
            %         JointMom_resj(1:6,j).*scaling.JRM_res)./scaling.T(1);
            % else
            %     eq_constr{end+1} = (ID_toTrackj(1:6,j)-Tj(ground_pelvisi,1))./scaling.T(1);
            % end
            if Options.usePelvisResMom  
                eq_constr{end+1} = (Pelvis_res_j(1:6,j)*scaling.T(1)-Tj(ground_pelvisi,1))./scaling.T(1);
                g_names_coll = [g_names_coll; repmat({'T_pelvis'},6,1)];
            end
           
            % Muscle-driven joint torques for the lower limbs and the trunk
            % Hip flexion, left
            Ft_hip_flex_l   = FTj(mai(1).mus.l',1);
            T_hip_flex_l    = f_T28(MAj.hip_flex.l,Ft_hip_flex_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
                    Tau_passj.hip.flex.l)-JointMom_resj(jointi.hip_flex.l,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
                    Tau_passj.hip.flex.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_flex_l'},1,1)];
            %Hip flexion, right
            Ft_hip_flex_r   = FTj(mai(1).mus.r',1);
            T_hip_flex_r    = f_T28(MAj.hip_flex.r,Ft_hip_flex_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_flex.r,1)- (T_hip_flex_r + ...
                    Tau_passj.hip.flex.r)-JointMom_resj(jointi.hip_flex.r,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_flex.r,1)- (T_hip_flex_r + ...
                    Tau_passj.hip.flex.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_flex_r'},1,1)];
            % Hip adduction left
            Ft_hip_add_l    = FTj(mai(2).mus.l',1);
            T_hip_add_l     = f_T28(MAj.hip_add.l,Ft_hip_add_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_add.l,1)-(T_hip_add_l + ...
                    Tau_passj.hip.add.l)-JointMom_resj(jointi.hip_add.l,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_add.l,1)-(T_hip_add_l + ...
                    Tau_passj.hip.add.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_add_l'},1,1)];
            % Hip adduction right
            Ft_hip_add_r    = FTj(mai(2).mus.r',1);
            T_hip_add_r     = f_T28(MAj.hip_add.r,Ft_hip_add_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_add.r,1)-(T_hip_add_r + ...
                    Tau_passj.hip.add.r)-JointMom_resj(jointi.hip_add.r,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_add.r,1)-(T_hip_add_r + ...
                    Tau_passj.hip.add.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_add_r'},1,1)];
            %Hip rotation left
            Ft_hip_rot_l    = FTj(mai(3).mus.l',1);
            T_hip_rot_l     = f_T28(MAj.hip_rot.l,Ft_hip_rot_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
                    Tau_passj.hip.rot.l)-JointMom_resj(jointi.hip_rot.l,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
                    Tau_passj.hip.rot.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_rot_l'},1,1)];
            %Hip rotation right
            Ft_hip_rot_r    = FTj(mai(3).mus.r',1);
            T_hip_rot_r     = f_T28(MAj.hip_rot.r,Ft_hip_rot_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
                    Tau_passj.hip.rot.r)-JointMom_resj(jointi.hip_rot.r,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
                    Tau_passj.hip.rot.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_hip_rot_r'},1,1)];
            %Knee add, right
            Ft_knee_add_r   = FTj(mai(4).mus.r',1);
            T_knee_add_r    = f_T13(MAj.knee_add.r,Ft_knee_add_r);
            if Options.useKCFresiduals
                eq_constr{end+1} = Tj(jointi.knee_add.r,1)-(T_knee_add_r + ...
                    Tau_passj.knee_add.r)+KCF_resj(1,j).*scaling.KCF_res(1);
            else
                eq_constr{end+1} = Tj(jointi.knee_add.r,1)-(T_knee_add_r + ...
                   Tau_passj.knee_add.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_add_r'},1,1)];
            %Knee rot, right
            Ft_knee_rot_r   = FTj(mai(5).mus.r',1);
            T_knee_rot_r    = f_T13(MAj.knee_rot.r,Ft_knee_rot_r);
            if Options.useKCFresiduals
                eq_constr{end+1} = Tj(jointi.knee_rot.r,1)-(T_knee_rot_r + ...
                    Tau_passj.knee_rot.r)+KCF_resj(2,j).*scaling.KCF_res(2);
            else
                eq_constr{end+1} = Tj(jointi.knee_rot.r,1)-(T_knee_rot_r + ...
                    Tau_passj.knee_rot.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_rot_r'},1,1)];
            %Knee flexion, left
            Ft_knee_flex_l  = FTj(mai(6).mus.l',1);
            T_knee_flex_l   = f_T13(MAj.knee_flex.l,Ft_knee_flex_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.knee_flex.l,1)-(T_knee_flex_l + ...
                    Tau_passj.knee_flex.l)-JointMom_resj(jointi.knee_flex.l,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.knee_flex.l,1)-(T_knee_flex_l + ...
                    Tau_passj.knee_flex.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_flex_l'},1,1)];
            %Knee flexion, right
            Ft_knee_flex_r= FTj(mai(6).mus.r',1);
            T_knee_flex_r   = f_T13(MAj.knee_flex.r,Ft_knee_flex_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.knee_flex.r,1)-(T_knee_flex_r + ...
                    Tau_passj.knee_flex.r)- JointMom_resj(jointi.knee_flex.l,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.knee_flex.r,1)-(T_knee_flex_r + ...
                    Tau_passj.knee_flex.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_flex_r'},1,1)];
            %Knee, tx
            Ft_knee_tx_r    = FTj(mai(7).mus.r',1);
            T_knee_tx_r     = f_T13(MAj.knee_tx.r,Ft_knee_tx_r);
            if Options.useKCFresiduals
                eq_constr{end+1} = Tj(jointi.knee_tx.r,1) - (T_knee_tx_r + ...
                    Tau_passj.knee_tx.r)+KCF_resj(3,j).*scaling.KCF_res(3);
            else
                eq_constr{end+1} = Tj(jointi.knee_tx.r,1) - (T_knee_tx_r + ...
                    Tau_passj.knee_tx.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_tx_r'},1,1)];
            % Knee ty, right
            Ft_knee_ty_r    = FTj(mai(8).mus.r',1);
            T_knee_ty_r     = f_T13(MAj.knee_ty.r,Ft_knee_ty_r);
            if Options.useKCFresiduals
                eq_constr{end+1} = Tj(jointi.knee_ty.r,1) - (T_knee_ty_r +...
                    Tau_passj.knee_ty.r)+ KCF_resj(4,j).*scaling.KCF_res(4); 
            else
                eq_constr{end+1} = Tj(jointi.knee_ty.r,1) - (T_knee_ty_r +...
                    Tau_passj.knee_ty.r); 
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_ty_r'},1,1)];
            % Knee tz, right
            Ft_knee_tz_r    = FTj(mai(9).mus.r',1);
            T_knee_tz_r     = f_T13(MAj.knee_tz.r,Ft_knee_tz_r);
            if Options.useKCFresiduals
                eq_constr{end+1} = Tj(jointi.knee_tz.r,1)- (T_knee_tz_r + ...
                    Tau_passj.knee_tz.r)+KCF_resj(5,j).*scaling.KCF_res(5);
            else
                eq_constr{end+1} = Tj(jointi.knee_tz.r,1)- (T_knee_tz_r + ...
                    Tau_passj.knee_tz.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_knee_tz_r'},1,1)];
            % Ankle, left
            Ft_ankle_l      = FTj(mai(10).mus.l',1);
            T_ankle_l       = f_T12(MAj.ankle.l,Ft_ankle_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.ankle.l,1)-(T_ankle_l + ...
                    Tau_passj.ankle.l)- JointMom_resj(jointi.ankle.l-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.ankle.l,1)-(T_ankle_l + ...
                    Tau_passj.ankle.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_ankle_l'},1,1)];
            % Ankle, right
            Ft_ankle_r      = FTj(mai(10).mus.r',1);
            T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.ankle.r,1)-(T_ankle_r + ...
                    Tau_passj.ankle.r)- JointMom_resj(jointi.ankle.r-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.ankle.r,1)-(T_ankle_r + ...
                    Tau_passj.ankle.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_ankle_r'},1,1)];
            % Subtalar, left
            Ft_subt_l       = FTj(mai(11).mus.l',1);
            T_subt_l        = f_T12(MAj.subt.l,Ft_subt_l);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.subt.l,1)-(T_subt_l + ...
                    Tau_passj.subt.l) - JointMom_resj(jointi.subt.l-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.subt.l,1)-(T_subt_l + ...
                    Tau_passj.subt.l);
            end
            g_names_coll = [g_names_coll; repmat({'T_subt_l'},1,1)];
            % Subtalar, right
            Ft_subt_r       = FTj(mai(11).mus.r',1);
            T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.subt.r,1)-(T_subt_r + ...
                    Tau_passj.subt.r) - JointMom_resj(jointi.subt.r-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.subt.r,1)-(T_subt_r + ...
                    Tau_passj.subt.r);
            end
            g_names_coll = [g_names_coll; repmat({'T_subt_r'},1,1)];
            % Lumbar extension
            Ft_trunk_ext    = FTj([mai(12).mus.l,mai(12).mus.r]',1);
            T_trunk_ext     = f_T6(MAj.trunk_ext,Ft_trunk_ext);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.trunk.ext,1)-(T_trunk_ext +...
                    Tau_passj.trunk.ext) - JointMom_resj(jointi.trunk.ext-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.trunk.ext,1)-(T_trunk_ext +...
                    Tau_passj.trunk.ext);
            end
            g_names_coll = [g_names_coll; repmat({'T_lumb_ext'},1,1)];
            % Lumbar bending
            Ft_trunk_ben    = FTj([mai(13).mus.l,mai(13).mus.r]',1);
            T_trunk_ben     = f_T6(MAj.trunk_ben,Ft_trunk_ben);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.trunk.ben,1)-(T_trunk_ben + ...
                    Tau_passj.trunk.ben) - JointMom_resj(jointi.trunk.ben-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.trunk.ben,1)-(T_trunk_ben + ...
                    Tau_passj.trunk.ben);
            end
            g_names_coll = [g_names_coll; repmat({'T_lumb_ben'},1,1)];
            % Lumbar rotation
            Ft_trunk_rot    = FTj([mai(14).mus.l,mai(14).mus.r]',1);
            T_trunk_rot     = f_T6(MAj.trunk_rot,Ft_trunk_rot);
            if Options.useJointResMom
                eq_constr{end+1} = Tj(jointi.trunk.rot,1)-(T_trunk_rot +...
                    Tau_passj.trunk.rot) - JointMom_resj(jointi.trunk.rot-5,j)*scaling.JRM_res;
            else
                eq_constr{end+1} = Tj(jointi.trunk.rot,1)-(T_trunk_rot +...
                    Tau_passj.trunk.rot);
            end
            g_names_coll = [g_names_coll; repmat({'T_lumb_rot'},1,1)];
            % Torque-driven joint torques for the arms
            % Arms
            eq_constr{end+1} = Tj(armsi,1)/scaling.ArmTau - a_akj(:,j+1);
            g_names_coll = [g_names_coll; repmat({'T_arms'},8,1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Activation dynamics (implicit formulation)   
            tact = 0.015;
            tdeact = 0.06;
            act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
            act2 = vAk_nsc + ak./(ones(size(ak,1),1)*tact);
            ineq_constr1{end+1} = act1;
            ineq_constr2{end+1} = act2; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Contraction dynamics (implicit formulation)
            eq_constr{end+1} = Hilldiffj;
            g_names_coll = [g_names_coll; repmat({'Hilldiff_eq'},NMuscle,1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %set lower bound for KCFs
            if Options.setLBKCF
                ineq_constr3{end+1}= Tj([KCFi.M]);
                ineq_constr3{end+1}= Tj([KCFi.L])
            end
        end % End loop over collocation points
        eq_constr = vertcat(eq_constr{:});
        ineq_constr1 = vertcat(ineq_constr1{:});
        ineq_constr2 = vertcat(ineq_constr2{:});
        ineq_constr3 = vertcat(ineq_constr3{:});
        if (Options.useJointResMom == 1)&&(Options.useKCFresiduals)
            %deprecated
            f_coll = Function('f_coll',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,KCF_resj,...
                JointMom_resj,Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj},{eq_constr,...
                ineq_constr1,ineq_constr2,J});
            f_coll_map = f_coll.map(N,parallelMode,NThreads);
            [coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2, ...
                Jall] = f_coll_map(paramsCM',a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
                X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
                KCF_res_col, Jointmom_res_col,...
                Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
                Qs_res_interpfilt_scaled_aux_col',GRF.val.allinterp_col(:,2:end)',...
                GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
                ID.allinterp_col(:,[2:(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])');
            %for debugging %%%
            f_coll2=Function('f_coll2',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,KCF_resj,...
                JointMom_resj,Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj},{J1,J2,J3,J4,J5,J6,J7,J7b,J8,J9,J10,J11,J12,Tj,in_F2',force_s1_l',force_s6_r',outk1j,Xk,Xj,Xkj,Xkj_nsc});
        elseif (Options.useJointResMom == 0)&&(Options.useKCFresiduals==1)
            f_coll = Function('f_coll',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,KCF_resj,...
                Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                Qdots_toTrack_scaled_auxk,Qdots_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj},{eq_constr,...
                ineq_constr1,ineq_constr2,J});
            f_coll_map = f_coll.map(N,parallelMode,NThreads);
            [coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2, ...
                Jall] = f_coll_map(paramsCM',a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
                X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
                KCF_res_col, ...
                Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
                Qs_res_interpfilt_scaled_aux_col',...
                Qdots_res_interpfilt_scaled_aux(1:end-1,:)',...
                Qdots_res_interpfilt_scaled_aux_col',...
                GRF.val.allinterp_col(:,2:end)',...
                GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
                ID.allinterp_col(:,[2:(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])');
            %for debugging %%%
            f_coll2=Function('f_coll2',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,KCF_resj,...
                Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                Qdots_toTrack_scaled_auxk,Qdots_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj},{J1,J1b,J2,J3,J4,J5,J6,J7,J7b,J8,J9,J10,J11,J12,Tj,in_F2',force_s1_l',force_s6_r',outk1j,Xk,Xj,Xkj,Xkj_nsc});
            f_coll_map2 = f_coll2.map(N,parallelMode,NThreads);    
            [J1all,J1ball,J2all,J3all,J4all,J5all,J6all,J7all,J7ball,J8all,J9all,J10all,J11all,J12all,Tkall,in_F2all,force_s1_lall,force_s6_rall,outk1all,Xkall,Xjall,Xkjall,Xkj_nscall] = f_coll_map2(paramsCM',a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
            X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
            KCF_res,...
            Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
            Qs_res_interpfilt_scaled_aux_col',Qdots_res_interpfilt_scaled_aux(1:end-1,:)',...
            Qdots_res_interpfilt_scaled_aux_col',GRF.val.allinterp_col(:,2:end)',...
            GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
            ID.allinterp_col(:,[(2+nq.abs):(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])');
        elseif (Options.useJointResMom == 0)&&(Options.useKCFresiduals==0)
            f_coll = Function('f_coll',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,...
                Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                Qdots_toTrack_scaled_auxk,Qdots_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj,Pelvis_res_j},{eq_constr,...
                ineq_constr1,ineq_constr2,ineq_constr3,J});
            f_coll_map = f_coll.map(N,parallelMode,NThreads);
            [coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2, ...
                coll_ineq_constr3,Jall] = f_coll_map(paramsCM',a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
                X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
                Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
                Qs_res_interpfilt_scaled_aux_col',...
                Qdots_res_interpfilt_scaled_aux(1:end-1,:)',...
                Qdots_res_interpfilt_scaled_aux_col',...
                GRF.val.allinterp_col(:,2:end)',...
                GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
                ID.allinterp_col(:,[2:(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])',Pelvis_res_col);
            %for debugging %%%
            f_coll2=Function('f_coll2',{paramsCMk,ak,aj,FTtildek,FTtildej,Xk,Xj,...
                a_ak,a_aj,vAk,e_ak,dFTtildej,Aj,...
                Qs_toTrack_scaled_auxk,Qs_toTrack_scaled_auxj,...
                Qdots_toTrack_scaled_auxk,Qdots_toTrack_scaled_auxj,...
                GRF_toTrackj,GRM_toTrackj,KCF_toTrackj,ID_toTrackj,Pelvis_res_j},{J1,J1b,J2,J3,J4,J5,J6,J7,J7b,J8,J9,J10,J11,J12,Tj,in_F2',force_s1_l',force_s6_r',outk1j,Xk,Xj,Xkj,Xkj_nsc});
            f_coll_map2 = f_coll2.map(N,parallelMode,NThreads);    
            [J1all,J1ball,J2all,J3all,J4all,J5all,J6all,J7all,J7ball,J8all,J9all,J10all,J11all,J12all,Tkall,in_F2all,force_s1_lall,force_s6_rall,outk1all,Xkall,Xjall,Xkjall,Xkj_nscall] = f_coll_map2(paramsCM',a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
            X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
            Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
            Qs_res_interpfilt_scaled_aux_col',Qdots_res_interpfilt_scaled_aux(1:end-1,:)',...
            Qdots_res_interpfilt_scaled_aux_col',GRF.val.allinterp_col(:,2:end)',...
            GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
            ID.allinterp_col(:,[2:(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])',Pelvis_res_col);
        end
    %         f_coll_map2 = f_coll2.map(N,parallelMode,NThreads);
    %     [J1all,J1ball,J2all,J3all,J4all,J5all,J6all,J7all,J8all,J9all,J10all,J11all,Tjall,in_F2all,force_s1_lall,force_s6_rall,outk1all,Xkall,Xjall,Xkjall,Xkj_nscall] = f_coll_map2(repmat(paramsCM',1,N),a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, ...
    %         X(:,1:end-1), X_col, a_a, a_a_col, vA, e_a, dFTtilde_col, A_col,...
    %         KCF_res_col, ...
    %         Qs_res_interpfilt_scaled_aux(1:end-1,:)',...
    %         Qs_res_interpfilt_scaled_aux_col',...
    %         Qdots_res_interpfilt_scaled_aux(1:end-1,:)',...
    %         Qdots_res_interpfilt_scaled_aux_col',...
    %         GRF.val.allinterp_col(:,2:end)',...
    %         GRF.MorGF.allinterp_col(:,2:end)',KCF.allinterpfilt_col',...
    %         ID.allinterp_col(:,[(2+nq.abs):(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])');
        %%%
        opti.subject_to(coll_eq_constr == 0);
        opti.subject_to(coll_ineq_constr1(:) >= 0);
        opti.subject_to(coll_ineq_constr2(:) <= 1/tact);
        if Options.setLBKCF
            opti.subject_to(coll_ineq_constr3(:) >= 0);
        end
        
        % Loop over mesh points
        X_nsc=MX.zeros(size(X));
        X_nsc(1:2:end,:)=X(1:2:end,:).*scaling.Qs';
        X_nsc(jointi.knee_ty.r*2-1,:)=X(jointi.knee_ty.r*2-1,:)*scaling.knee_ty.b+scaling.knee_ty.a;
        X_nsc(2:2:end,:)=X(2:2:end,:).*scaling.Qdots';
        for k=1:N
            % Variables within current mesh interval
            % States   
            akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)];
            FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
            % Xkj = [X(:,k), X_col(:,(k-1)*d+1:k*d)]; do not use the scaled
            % version, otherwise there will be an issue with knee_ty
            Xkj_nsc=[Xkj_nscall(:,(k-1)*(d+1)+1:k*(d+1))];
            a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
            % Add equality constraints (next interval starts with end values of 
            % states from previous interval)
            opti.subject_to(a(:,k+1) == akj*D);
            opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
            % opti.subject_to(X(:,k+1) == Xkj*D); % scaled
            opti.subject_to(X_nsc(:,k+1) == Xkj_nsc*D); % unscaled
        end
        Jall=sum(Jall);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % Create NLP solver
        opti.minimize(Jall);
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.mu_strategy      = 'adaptive';
        options.ipopt.max_iter = 10000;
        options.ipopt.tol = 1*10^(-tol_ipopt);
        % options.ipopt.warm_start_init_point='yes';
        options.ipopt.print_timing_statistics='yes';
        opti.solver('ipopt', options);  
        % Create and save diary
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'/Results'];
        if ~(exist([pathresults,'/',namescript],'dir')==7)
            mkdir(pathresults,namescript);
        end
        if (exist([pathresults,'/',namescript,'/D',savename],'file')==2)
            delete ([pathresults,'/',namescript,'/D',savename])
        end 
        diary([pathresults,'/',namescript,'/D',savename]);  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve problem
        % Opti does not use bounds on variables but constraints. This function
        % adjusts for that.    
        [w_opt,stats,g_opt] = solve_NLPSOL(opti,options);  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        diary off
        % Extract results 
        % Create setup
        setup.tolerance.ipopt = tol_ipopt;
        setup.bounds = bounds;
        setup.scaling = scaling;
        setup.guess = guess;
        % Save results and set    
        save([pathresults,'/',namescript,'/w',savename],'w_opt');
        save([pathresults,'/',namescript,'/stats',savename],'stats');
        
        if analyseResults
            %% Load results
            if loadResults
                p = mfilename('fullpath');
                [~,namescript,~] = fileparts(p);
                pathresults = [pathRepo,'/Results'];
                load([pathresults,'/',namescript,'/w',savename]);
                load([pathresults,'/',namescript,'/g',savename]);
                load([pathresults,'/',namescript,'/s',savename]);
                load([pathresults,'/',namescript,'/stats',savename]);
            end
    %         if Options.useKCFresiduals
    %             if Options.useJointResMom
    %                 NControls = NMuscle+NMuscle+nq.all+nq.arms+5+21;
    %             else
    %                 NControls = NMuscle+NMuscle+nq.all+nq.arms+5;
    %             end
    %         else
    %             if Options.useJointResMom
    %                 NControls = NMuscle+NMuscle+nq.all+nq.arms+21;
    %             else
    %                 NControls = NMuscle+NMuscle+nq.all+nq.arms;
    %             end
    %         end
            NStates = NMuscle+NMuscle+2*nq.all+nq.arms;
            NParameters = np;
            % Here we extract the results and re-organize them for analysis  
            % Static parameters
            paramsCM_opt    = w_opt(1:NParameters);
            starti = NParameters+1;
            % Muscle activations at mesh points
            a_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
            starti = starti + NMuscle*(N+1);
            % Muscle activations at collocation points
            a_col_opt = reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
            starti = starti + NMuscle*(d*N);
            % Tendon force at mesh points
            FTtilde_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
            starti = starti + NMuscle*(N+1);
            % Tendon force at collocation points
            FTtilde_col_opt =reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
            starti = starti + NMuscle*(d*N);
            % Q's and Qdot's at mesh points
            X_opt = reshape(w_opt(starti:starti+2*nq.all*(N+1)-1),2*nq.all,N+1)';
            starti = starti + 2*nq.all*(N+1);
            %Q's and Qdot's at collocation points
            X_col_opt = reshape(w_opt(starti:starti+2*nq.all*(d*N)-1),2*nq.all,d*N)';
            starti = starti + 2*nq.all*(d*N);
            % arm activations at mesh points
            a_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
            starti = starti + nq.arms*N;
            % arm activations at collocation points
            a_a_col_opt = reshape(w_opt(starti:starti+nq.arms*(d*N)-1),nq.arms,d*N)';
            starti = starti + nq.arms*(d*N);
            %Controls
            % Time derivative of muscle activations at mesh points
            vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
            starti = starti + NMuscle*N;
            % Arm excitations at mesh points        
            e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
            starti = starti + nq.arms*N;   
            %KCF residuals at collocation points if any
            if Options.useKCFresiduals
               KCF_res_opt = reshape(w_opt(starti:starti+5*(d*N)-1),5,d*N)';
               starti = starti + 5*(d*N);
            end
            %Joint residual moments at collocation points if any
            if Options.useJointResMom
               JointMom_res_opt = reshape(w_opt(starti:starti+21*(d*N)-1),21,d*N)';
               starti = starti + 21*(d*N);
            end
            if Options.usePelvisResMom
                Pelvis_res_opt = reshape(w_opt(starti:starti+6*(d*N)-1),6,d*N)';
                starti = starti + 6*(d*N);
            end
            %Slack controls
            % Time derivatives of tendon forces at collocation points
            dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
            starti = starti + NMuscle*(d*N);
            % Time derivatives of joint velocities
            qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
            starti = starti + nq.all*(d*N);
    
             %% Unscale results
            % Parameters
            paramsCM_opt_nsc = f_nsc18(paramsCM_opt,scaling.params.v,scaling.params.r);
          
            q_opt=X_opt(:,1:2:end);
            qdot_opt=X_opt(:,2:2:end);
            % Qs (1:N-1)
            q_opt_unsc.rad = q_opt(1:end,:).*repmat(...
            scaling.Qs,size(q_opt(1:end,:),1),1); 
            q_opt_unsc.rad(:,jointi.knee_ty.r) = q_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
            % Convert in degrees
            q_opt_unsc.deg = q_opt_unsc.rad;
            q_opt_unsc.deg(:,[1:3,7:16 20:end]) = q_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
            % Qs (1:N)
            q_opt_unsc_all.rad = q_opt(1:end,:).*repmat(...
                scaling.Qs,size(q_opt(1:end,:),1),1); 
            q_opt_unsc_all.rad(:,jointi.knee_ty.r)=q_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
            % Convert in degrees
            q_opt_unsc_all.deg = q_opt_unsc_all.rad;
            q_opt_unsc_all.deg(:,[1:3,7:16 20:end]) = ...
                q_opt_unsc_all.deg(:,[1:3,7:16 20:end]).*180/pi;
            % Qdots (1:N-1)
            qdot_opt_unsc.rad = qdot_opt(1:end,:).*repmat(...
                scaling.Qdots,size(qdot_opt(1:end,:),1),1);
            % Convert in degrees
            qdot_opt_unsc.deg = qdot_opt_unsc.rad;
            qdot_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
            qdot_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
            
            % Muscle activations
            a_opt_unsc = a_opt(1:end,:).*repmat(...
                scaling.a,size(a_opt(1:end,:),1),size(a_opt,2));
            % Muscle-tendon forces
            FTtilde_opt_unsc = FTtilde_opt(1:end,:).*repmat(...
                scaling.FTtilde,size(FTtilde_opt(1:end,:),1),1);
        
            % State and Controls at collocation points
            % Qs and Qdots
            q_col_opt=X_col_opt(:,1:2:end);
            qdot_col_opt=X_col_opt(:,2:2:end);
            q_col_opt_unsc.rad = q_col_opt(1:end,:).*repmat(...
            scaling.Qs,size(q_col_opt(1:end,:),1),1); 
            q_col_opt_unsc.rad(:,jointi.knee_ty.r) = q_col_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
            qdot_col_opt_unsc.rad = qdot_col_opt(1:end,:).*repmat(...
                scaling.Qdots,size(qdot_col_opt(1:end,:),1),1);
        
            % Time derivative of Qdots
            qdotdot_col_opt_unsc.rad = ...
                qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
            % Convert in degrees
            qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
            qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
                qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
        
            %At mesh points
            % Muscle-tendon forces
            FTtilde_opt_unsc = FTtilde_opt(1:end,:).*repmat(...
                scaling.FTtilde,size(FTtilde_opt(1:end,:),1),1);
            % Arm activations
            a_a_opt_unsc = a_a_opt(1:end,:).*repmat(...
                scaling.a_a,size(a_a_opt(1:end,:),1),size(a_a_opt,2));
            
            % Time derivative of muscle activations (states)
            vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
            tact = 0.015;
            tdeact = 0.06;
            % Get muscle excitations from time derivative of muscle activations
            e_opt_unsc = computeExcitationRaasch(a_opt_unsc(1:N,:),vA_opt_unsc,...
                ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
            %At collocation points
            %Muscle-tendon forces
            FTtilde_col_opt_unsc= FTtilde_col_opt.*repmat(scaling.FTtilde,size(FTtilde_col_opt,1),1);
            % Time derivative of muscle-tendon forces
            dFTtilde_col_opt_unsc = dFTtilde_col_opt.*repmat(...
                scaling.dFTtilde,size(dFTtilde_col_opt,1),size(dFTtilde_col_opt,2));
            a_a_col_opt_unsc=a_a_col_opt.*repmat(scaling.a_a,size(a_a_col_opt,1),size(a_a_col_opt,2));
            % Arm excitations
            e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
                size(e_a_opt,2));
            if Options.useKCFresiduals
                %KCF residuals
                KCF_res_opt_unsc = KCF_res_opt.*repmat(scaling.KCF_res,size(KCF_res_opt,1),1);
            else
                KCF_res_opt_unsc = zeros(size(e_a_opt_unsc,1),5);
                KCF_res_opt=zeros(size(e_a_opt_unsc,1),5);
            end
            if Options.useJointResMom
                JointMom_res_opt_unsc = JointMom_res_opt.*repmat(scaling.JRM_res,N,21);
            else
                JointMom_res_opt_unsc = zeros(N,21);
                JointMom_res_opt = zeros(N,21);
            end
            if Options.usePelvisResMom
                Pelvis_res_opt_unsc=Pelvis_res_opt.*repmat(scaling.T(1),N*d,6);
            else
                Pelvis_res_opt_unsc=zeros(N*d,6);
                Pelvis_res_opt=zeros(N*d,6);
            end
            
            %% Time grid    
            % Mesh points
            tgrid = linspace(time_opt(1),time_opt(end),N+1);
            dtime = zeros(1,d+1);
            for i=1:4
                dtime(i)=tau_root(i)*((time_opt(end)-time_opt(1))/N);
            end
            % Mesh points and collocation points
            tgrid_ext = zeros(1,(d+1)*N+1);
            for i=1:N
                tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
            end
            tgrid_ext(end)=time_opt(end); 
    
            %% Joint torques, ground reaction forces and moments at opt solution
            Xk_Qs_Qdots_opt                 = zeros(N+1,2*size(q_opt_unsc.rad,2));
            Xk_Qs_Qdots_opt(:,1:2:end)      = q_opt_unsc.rad(:,:);
            Xk_Qs_Qdots_opt(:,2:2:end)      = qdot_opt_unsc.rad(:,:);
            Xk_Qdotdots_opt                 = qdotdot_col_opt_unsc.rad(3:3:end,:); %to be removed!! 

            Xk_Qs_Qdots_col_opt             = zeros(N*d,2*size(q_opt_unsc.rad,2));
            Xk_Qs_Qdots_col_opt(:,1:2:end)  = q_col_opt_unsc.rad(1:N*d,:);
            Xk_Qs_Qdots_col_opt(:,2:2:end)  = qdot_col_opt_unsc.rad(1:N*d,:);
            Xk_Qdotdots_col_opt             = qdotdot_col_opt_unsc.rad;  

            for i = 1:N
                if Options.KCFasinputstoExternalFunction
                        [~, GRF_opt_unsc_k(i,:), GRM_opt_unsc_k(i,:), KCF_opt_unsc_k(i,:)]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt(i,:),Xk_Qdotdots_opt(i,:),F1,{F2_skeletal F2},paramsCM_opt_nsc,loc,loci,stif,radi,dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,toes,calcn,deri,f_contactForce,residualsi,GRFi,GRMi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc(i,:),Options,jointi);
                else
                        [~, GRF_opt_unsc_k(i,:), GRM_opt_unsc_k(i,:), KCF_opt_unsc_k(i,:)]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt(i,:),Xk_Qdotdots_opt(i,:),F1,F2,paramsCM_opt_nsc,loc,loci,stif,radi,dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,toes,calcn,deri,f_contactForce,residualsi,GRFi,GRMi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc(i,:),Options,jointi);
                end
                for j=1:d
                    if Options.KCFasinputstoExternalFunction
                        [Tauk_out((i-1)*d+j,:), GRF_opt_unsc((i-1)*d+j,:), GRM_opt_unsc((i-1)*d+j,:), KCF_opt_unsc((i-1)*d+j,:)]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_col_opt((i-1)*d+j,:),Xk_Qdotdots_col_opt((i-1)*d+j,:),...
                            F1,{F2_skeletal F2},paramsCM_opt_nsc,loc,loci,stif,radi,dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,toes,calcn,deri,f_contactForce,residualsi,GRFi,GRMi,KCFi,armsi,scaling,tol_ipopt,a_a_col_opt_unsc((i-1)*d+j,:),Options,jointi);
                    else
                        [Tauk_out((i-1)*d+j,:), GRF_opt_unsc((i-1)*d+j,:), GRM_opt_unsc((i-1)*d+j,:), KCF_opt_unsc((i-1)*d+j,:)]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_col_opt((i-1)*d+j,:),Xk_Qdotdots_col_opt((i-1)*d+j,:),...
                            F1,F2,paramsCM_opt_nsc,loc,loci,stif,radi,dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,toes,calcn,deri,f_contactForce,residualsi,GRFi,GRMi,KCFi,armsi,scaling,tol_ipopt,a_a_col_opt_unsc((i-1)*d+j,:),Options,jointi);
                    end
                end
            end



            %% Muscle data / Hill's model variables
        for i=1:N
            for j=1:d
                % Get muscle-tendon lengths, velocities, and moment arms
                % Left leg
                qin_l = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2-1,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2-1,1), ...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2-1,1), 0,0,...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2-1,1),0,0.042,0 ...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2-1,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2-1,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1,1)];  
                qdotin_l = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2,1),0,0,...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2,1),0,0,0,...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2,1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2,1)];  
                [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);    
                MA.hip_flex.l   =  MA_l(mai(1).mus.l',1);
                MA.hip_add.l    =  MA_l(mai(2).mus.l',2);
                MA.hip_rot.l    =  MA_l(mai(3).mus.l',3);
                MA.knee_flex.l  =  MA_l(mai(6).mus.l',6);
                MA.ankle.l      =  MA_l(mai(10).mus.l',10);  
                MA.subt.l       =  MA_l(mai(11).mus.l',11); 
                % For the back muscles, we want left and right together: left
                % first, right second. In MuscleInfo, we first have the right
                % muscles (45:47) and then the left muscles (48:50). Since the back
                % muscles only depend on back dofs, we do not care if we extract
                % them "from the left or right leg" so here we just picked left.
                MA.trunk_ext    =  MA_l([48:50,mai(12).mus.l]',12);
                MA.trunk_ben    =  MA_l([48:50,mai(13).mus.l]',13);
                MA.trunk_rot    =  MA_l([48:50,mai(14).mus.l]',14);
                % Right leg
                qin_r = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1)];  
                qdotin_r = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2),...
                    Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2)];      
                [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
                % Here we take the indices from left since the vector is 1:47 (is
                % 47?, are right instead of left...)
                MA.hip_flex.r   =  MA_r(mai(1).mus.l',1);
                MA.hip_add.r    =  MA_r(mai(2).mus.l',2);
                MA.hip_rot.r    =  MA_r(mai(3).mus.l',3);
                MA.knee_add.r   =  MA_r(mai(4).mus.l',4);
                MA.knee_rot.r   =  MA_r(mai(5).mus.l',5);
                MA.knee_flex.r  =  MA_r(mai(6).mus.l',6);
                MA.knee_tx.r    =  MA_r(mai(7).mus.l',7);
                MA.knee_ty.r    =  MA_r(mai(8).mus.l',8);
                MA.knee_tz.r    =  MA_r(mai(9).mus.l',9);
                MA.ankle.r      =  MA_r(mai(10).mus.l',10);
                MA.subt.r       =  MA_r(mai(11).mus.l',11);
                MA_mat((i-1)*d+j).MA=MA; %store matrix of MA for later calculations
                % Both legs
                % In MuscleInfo, we first have the right back muscles (45:47) and 
                % then the left back muscles (48:50). Here we re-organize so that
                % we have first the left muscles and then the right muscles.
                lMTk_lr     = [lMTk_l([1:44,48:50],1);lMTk_r(1:47,1)];
                vMTk_lr     = [vMTk_l([1:44,48:50],1);vMTk_r(1:47,1)]; 
        
        [Hilldiffk_opt_aux,FTk_opt_aux,~,~,~] =  f_forceEquilibrium_FtildeState(...
                    a_col_opt((i-1)*d+j,:),FTtilde_col_opt_unsc((i-1)*d+j,:),...
                    dFTtilde_col_opt_unsc((i-1)*d+j,:), lMTk_lr,vMTk_lr);
                Hilldiffk_opt((i-1)*d+j,:)=full(Hilldiffk_opt_aux);
                FTk_opt((i-1)*d+j,:)=full(FTk_opt_aux);
                
        [Hilldiffk_opt_aux2(:,(i-1)*d+j), FTk_opt_aux2(:,(i-1)*d+j), Fce_opt_aux2(:,(i-1)*d+j), Fiso, vMmax, Fpe_opt_aux2(:,(i-1)*d+j), lMtilde_opt_aux2(:,(i-1)*d+j), FMvtilde_aux2(:,(i-1)*d+j), vMtilde_aux2(:,(i-1)*d+j)] = ...
        ForceEquilibrium_FtildeState_GC_test(a_col_opt((i-1)*d+j,:),FTtilde_col_opt_unsc((i-1)*d+j,:),dFTtilde_col_opt_unsc((i-1)*d+j,:),full(lMTk_lr)',full(vMTk_lr)',MTparameters_m,Fvparam,...
            Fpparam,Faparam);    
        
            MAr_all((i-1)*d+j,:,:)=full(MA_r);
            MAl_all((i-1)*d+j,:,:)=full(MA_l);
    
        % Get optimal passive torques
        Tau_pass_opt.hip.flex.l((i-1)*d+j,:)    = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2)));
        Tau_pass_opt.hip.flex.r((i-1)*d+j,:)    = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2)));
        Tau_pass_opt.hip.add.l((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2)));
        Tau_pass_opt.hip.add.r((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2)));
        Tau_pass_opt.hip.rot.l((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2)));
        Tau_pass_opt.hip.rot.r((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2))); 
        if Options.dampingInKneeSec
            Tau_pass_opt.knee_add.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2)));
            Tau_pass_opt.knee_rot.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2)));
        else
            Tau_pass_opt.knee_add.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1)));
            Tau_pass_opt.knee_rot.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1)));
        end
        Tau_pass_opt.knee_flex.l((i-1)*d+j,:)   = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2)));
        if Options.dampingInKneeSec
            Tau_pass_opt.knee_tx.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2)));            
            Tau_pass_opt.knee_ty.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(                                                  0,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2)));            
            Tau_pass_opt.knee_tz.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2)));
        else
            Tau_pass_opt.knee_tx.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1)));            
            Tau_pass_opt.knee_ty.r((i-1)*d+j,:)     = 0;
            Tau_pass_opt.knee_tz.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1)));
        end
        Tau_pass_opt.knee_flex.r((i-1)*d+j,:)        = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2)));
        Tau_pass_opt.ankle.l((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2)));
        Tau_pass_opt.ankle.r((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2)));        
        Tau_pass_opt.subt.l((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2)));
        Tau_pass_opt.subt.r((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2)));        
        Tau_pass_opt.trunk.ext((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2)));
        Tau_pass_opt.trunk.ben((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2)));
        Tau_pass_opt.trunk.rot((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1),...
            Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2))); 
       
        %% Recalculate constraint values for joint moments
    
        % % Pelvis residuals (same as from inverse dynamics)
        if Options.usePelvisResMom
            eq_constr_opt((i-1)*d+j, 1:6) = (Pelvis_res_opt_unsc((i-1)*d+j,:)-...
                    Tauk_out((i-1)*d+j,ground_pelvisi))./scaling.T(1);
        else
            if Options.useJointResMom
                eq_constr_opt((i-1)*d+j, 1:6) = (ID.allinterp_col((i-1)*d+j,[(2+nq.abs):(jointi.knee_flex.r+1) (jointi.ankle.r+1):(nq.all+1)])'-Tauk_out((i-1)*d+j,ground_pelvisi)-...
                    JointMom_res_opt(:,1:6).*scaling.JRM_res)./scaling.T(1);
            else
                eq_constr_opt((i-1)*d+j, 1:6) = (ID.allinterp_col((i-1)*d+j,2:7)-Tauk_out((i-1)*d+j,ground_pelvisi))./scaling.T(1);
            end
        end

        % Hip flexion right
        Ft_hip_flex_r_opt  = FTk_opt((i-1)*d+j,mai(1).mus.r');
        T_hip_flex_r_opt   = full(f_T28(MA.hip_flex.r,Ft_hip_flex_r_opt));
        if Options.useJointResMom
            eq_constr_opt((i-1)*d+j, jointi.hip_flex.r) = Tauk_out((i-1)*d+j, jointi.hip_flex.r)-(T_hip_flex_r_opt + ...
                Tau_pass_opt.hip.flex.r((i-1)*d+j,:))-JointMom_res_opt(:,jointi.knee_flex.r)*scaling.JRM_res;
        else
            eq_constr_opt((i-1)*d+j, jointi.hip_flex.r) = Tauk_out((i-1)*d+j, jointi.hip_flex.r)-(T_hip_flex_r_opt + ...
                Tau_pass_opt.hip.flex.r((i-1)*d+j,:));
        end
        
        % Hip flexion left
        
        Ft_hip_flex_l_opt  = FTk_opt((i-1)*d+j,mai(1).mus.l');
        T_hip_flex_l_opt   = full(f_T28(MA.hip_flex.l,Ft_hip_flex_l_opt));
        
        
        if Options.useJointResMom
            eq_constr_opt((i-1)*d+j, jointi.hip_flex.l) = Tauk_out((i-1)*d+j, jointi.hip_flex.l)-(T_hip_flex_l_opt + ...
                Tau_pass_opt.hip.flex.l((i-1)*d+j,:))-JointMom_res_opt(jointi.knee_flex.l,1)*scaling.JRM_res;
        else
            eq_constr_opt((i-1)*d+j, jointi.hip_flex.l) = Tauk_out((i-1)*d+j, jointi.hip_flex.l)-(T_hip_flex_l_opt + ...
                Tau_pass_opt.hip.flex.l((i-1)*d+j,:));
        end
    
        % 
        % % Hip adduction right
        % 
        % Ft_hip_add_r_opt  = FTk_opt(i,mai(2).mus.r');
        % T_hip_add_r_opt   = full(f_T28(MAk.hip_add.r,Ft_hip_add_r_opt));
        % 
        % 
        % if Options.useJointResMom
        %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
        %         Tau_pass_opt.hip.add.r(i,:))-JointMom_res_opt(jointi.knee_add.r,1)*scaling.JRM_res;
        % else
        %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
        %         Tau_pass_opt.hip.add.r(i,:));
        % end
        % 
        % 
        % 
        % % Hip adduction left
        % 
        % Ft_hip_add_l_opt  = FTk_opt(i,mai(2).mus.l');
        % T_hip_add_l_opt   = full(f_T28(MAk.hip_add.l,Ft_hip_add_l_opt));
        % 
        % 
        % if Options.useJointResMom
        %     eq_constr_opt(i, jointi.hip_add.l) = Tauk_out(i, jointi.hip_add.l)-(T_hip_add_l_opt + ...
        %         Tau_pass_opt.hip.add.l(i,:))-JointMom_res_opt(jointi.knee_add.l,1)*scaling.JRM_res(1);
        % else
        %     eq_constr_opt(i, jointi.hip_add.l) = Tauk_out(i, jointi.hip_add.l)-(T_hip_add_l_opt + ...
        %         Tau_pass_opt.hip.add.l(i,:));
        % end
        % 
        % 
        % 
        % % Hip rotation right
        % 
        % Ft_hip_add_r_opt  = FTk_opt(i,mai(3).mus.r');
        % T_hip_add_r_opt   = full(f_T28(MAk.hip_add.r,Ft_hip_add_r_opt));
        % 
        % 
        % if Options.useJointResMom
        %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
        %         Tau_pass_opt.hip.add.r(i,:))-JointMom_res_opt(jointi.knee_add.r,1)*scaling.JRM_res(2);
        % else
        %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
        %         Tau_pass_opt.hip.add.r(i,:));
        % end
        % 
        % 
        % 
        % % Hip rotation left
        % 
        % Ft_hip_rot_l_opt  = FTk_opt(i,mai(3).mus.l');
        % T_hip_rot_l_opt   = full(f_T28(MAk.hip_rot.l,Ft_hip_rot_l_opt));
        % 
        % 
        % if Options.useJointResMom
        %     eq_constr_opt(i, jointi.hip_rot.l) = Tauk_out(i, jointi.hip_rot.l)-(T_hip_rot_l_opt + ...
        %         Tau_pass_opt.hip.rot.l(i,:))-JointMom_res_opt(jointi.knee_rot.l,1)*scaling.JRM_res;
        % else
        %     eq_constr_opt(i, jointi.hip_rot.l) = Tauk_out(i, jointi.hip_rot.l)-(T_hip_rot_l_opt + ...
        %         Tau_pass_opt.hip.rot.l(i,:));
        % end
    
    
        %Knee add, right
        Ft_knee_add_r_opt   = FTk_opt((i-1)*d+j,mai(4).mus.r');
        T_knee_add_r_opt    = full(f_T13(MA.knee_add.r,Ft_knee_add_r_opt));
        if Options.useKCFresiduals
            eq_constr_opt((i-1)*d+j,jointi.knee_add.r) = Tauk_out((i-1)*d+j,jointi.knee_add.r)-(T_knee_add_r_opt + ...
                Tau_pass_opt.knee_add.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,1).*scaling.KCF_res(1);
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_add.r) = Tauk_out((i-1)*d+j,jointi.knee_add.r,1)-(T_knee_add_r_opt + ...
               Tau_pass_opt.knee_add.r((i-1)*d+j,:));
        end
        %Knee rot, right
        Ft_knee_rot_r_opt   = FTk_opt((i-1)*d+j,mai(5).mus.r');
        T_knee_rot_r_opt    = full(f_T13(MA.knee_rot.r,Ft_knee_rot_r_opt));
        if Options.useKCFresiduals
            eq_constr_opt((i-1)*d+j,jointi.knee_rot.r) = Tauk_out((i-1)*d+j,jointi.knee_rot.r)-(T_knee_rot_r_opt + ...
                Tau_pass_opt.knee_rot.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,2).*scaling.KCF_res(2);
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_rot.r) = Tauk_out((i-1)*d+j,jointi.knee_rot.r)-(T_knee_rot_r_opt + ...
                Tau_pass_opt.knee_rot.r((i-1)*d+j,:));
        end
        %Knee flexion, left
        Ft_knee_flex_l_opt  = FTk_opt((i-1)*d+j,mai(6).mus.l');
        T_knee_flex_l_opt   = full(f_T13(MA.knee_flex.l,Ft_knee_flex_l_opt));
        if Options.useJointResMom
            eq_constr_opt((i-1)*d+j,jointi.knee_flex.l) = Tauk_out((i-1)*d+j,jointi.knee_flex.l)-(T_knee_flex_l_opt + ...
                Tau_pass_opt.knee_flex.l((i-1)*d+j,:))-JointMom_res_opt(jointi.knee_flex.l,1)*scaling.JRM_res;
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_flex.l)= Tauk_out((i-1)*d+j,jointi.knee_flex.l)-(T_knee_flex_l_opt + ...
                Tau_pass_opt.knee_flex.l((i-1)*d+j,:));
        end
        %Knee flexion, right
        Ft_knee_flex_r_opt= FTk_opt((i-1)*d+j,mai(6).mus.r');
        T_knee_flex_r_opt   = full(f_T13(MA.knee_flex.r,Ft_knee_flex_r_opt));
        if Options.useJointResMom
            eq_constr_opt((i-1)*d+j,jointi.knee_flex.r) = Tauk_out((i-1)*d+j,jointi.knee_flex.r)-(T_knee_flex_r_opt + ...
                Tau_pass_opt.knee_flex.r((i-1)*d+j,:))- JointMom_res_opt(jointi.knee_flex.r,1)*scaling.JRM_res;
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_flex.r) = Tauk_out((i-1)*d+j,jointi.knee_flex.r)-(T_knee_flex_r_opt + ...
                Tau_pass_opt.knee_flex.r((i-1)*d+j,:));
        end
        %Knee, tx
        Ft_knee_tx_r_opt    = FTk_opt((i-1)*d+j,mai(7).mus.r');
        T_knee_tx_r_opt     = full(f_T13(MA.knee_tx.r,Ft_knee_tx_r_opt));
        if Options.useKCFresiduals
            eq_constr_opt((i-1)*d+j,jointi.knee_tx.r) = Tauk_out((i-1)*d+j,jointi.knee_tx.r) - (T_knee_tx_r_opt + ...
                Tau_pass_opt.knee_tx.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,3).*scaling.KCF_res(3);
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_tx.r) = Tauk_out((i-1)*d+j,jointi.knee_tx.r) - (T_knee_tx_r_opt + ...
                Tau_pass_opt.knee_tx.r((i-1)*d+j,:));
        end
        % Knee ty, right
        Ft_knee_ty_r_opt    = FTk_opt((i-1)*d+j,mai(8).mus.r');
        T_knee_ty_r_opt((i-1)*d+j,:)     = full(f_T13(MA.knee_ty.r,Ft_knee_ty_r_opt));
        if Options.useKCFresiduals
            eq_constr_opt((i-1)*d+j,jointi.knee_ty.r) = Tauk_out((i-1)*d+j,jointi.knee_ty.r) - (T_knee_ty_r_opt((i-1)*d+j,:)+...
                Tau_pass_opt.knee_ty.r((i-1)*d+j,:))+ KCF_res_opt((i-1)*d+j,4).*scaling.KCF_res(4); 
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_ty.r) = Tauk_out((i-1)*d+j,jointi.knee_ty.r) - (T_knee_ty_r_opt((i-1)*d+j,:)+...
                Tau_pass_opt.knee_ty.r((i-1)*d+j,:)); 
        end
        % Knee tz, right
        Ft_knee_tz_r_opt    = FTk_opt((i-1)*d+j,mai(9).mus.r');
        T_knee_tz_r_opt     = full(f_T13(MA.knee_tz.r,Ft_knee_tz_r_opt));
        if Options.useKCFresiduals
            eq_constr_opt((i-1)*d+j,jointi.knee_tz.r) = Tauk_out((i-1)*d+j,jointi.knee_tz.r)- (T_knee_tz_r_opt + ...
                Tau_pass_opt.knee_tz.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,5).*scaling.KCF_res(5);
        else
            eq_constr_opt((i-1)*d+j,jointi.knee_tz.r) = Tauk_out((i-1)*d+j,jointi.knee_tz.r)- (T_knee_tz_r_opt + ...
                Tau_pass_opt.knee_tz.r((i-1)*d+j,:));
        end
    
    
        %% Recalculate dynamic constraints
        % Skeleton dynamics (implicit formulation)      
    
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Qsp_nsc          = Xkj_nsc(1:2:end,:)*C(:,j+1);
    %         Qdotsp_nsc       = Xkj_nsc(2:2:end,:)*C(:,j+1);    
    %         % Skeleton dynamics (implicit formulation)               
    %         qdotj_nsc = Xkj_nsc(2:2:end,j+1); % velocity
    %         eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
    %         eq_constr{end+1} = (h*Ak_nsc - Qdotsp_nsc)./...
    %                scaling.QsQdots(2:2:end)';
            %%%%%%%%%%%%%%%%%%%%%%%%
    
    
            X_col_opt_nsc = X_col_opt.*(scaling.QsQdots.*ones(size(X_col_opt,1),1));
            X_col_opt_nsc(:,jointi.knee_ty.r*2-1) = X_col_opt(:,jointi.knee_ty.r*2-1).*scaling.knee_ty.b+scaling.knee_ty.a;
            
            q_colkj=[q_opt_unsc.rad(i,:); X_col_opt_nsc((i-1)*d+1:i*d,1:2:end)]; %coordinates for all 4 collocation points of interval i unscaled
            qdot_colkj =[qdot_opt_unsc.rad(i,:); X_col_opt_nsc((i-1)*d+1:i*d,2:2:end)]; % velocity
    
            Qsp_opt_nsc = q_colkj'*C(:,j+1);
            Qdotsp_opt_nsc = qdot_colkj'*C(:,j+1);
            
            FTtildekj_opt_unsc=[FTtilde_opt_unsc(i,:); FTtilde_col_opt_unsc(((i-1)*d+1):((i-1)*d+3),:)];
            FTtildep_opt_unsc=FTtildekj_opt_unsc'*C(:,j+1);
    
            akj_opt_unsc=[a_opt_unsc(i,:); a_col_opt(((i-1)*d+1):((i-1)*d+3),:)];
            ap_opt_unsc=akj_opt_unsc'*C(:,j+1);
    
            eq_constr_dynkin((i-1)*d+j,:)= (h*qdot_colkj(j+1,:) - Qsp_opt_nsc')./scaling.QsQdots(1:2:end);
            
            eq_constr_dynkinvel((i-1)*d+j,:)= (h*qdotdot_col_opt_unsc.rad((i-1)*d+j,:) - Qdotsp_opt_nsc')./scaling.QsQdots(2:2:end);
           
            eq_constr_dynFT((i-1)*d+j,:)=(h*dFTtilde_col_opt_unsc((i-1)*d+j,:)-FTtildep_opt_unsc')./scaling.FTtilde;
    
            eq_constr_act((i-1)*d+j,:)=(h*vA_opt_unsc(i,:)-ap_opt_unsc')./scaling.a;
    
        end
        end
        
        %% Recalculate cost function values at optimal solution
        X_aux = opti.x;
    %     guess=opti.debug.value(X_aux,opti.initial);
        f_cost=Function('f_cost',{X_aux},{J1all,J1ball,J2all,J3all,J4all,J5all,J6all,J7all,J7ball,J8all,J9all,J10all,J11all,J12all});
        [J1all_val,J1ball_val,J2all_val,J3all_val,J4all_val,J5all_val,J6all_val,J7all_val,J7ball_val,J8all_val,J9all_val,J10all_val,J11all_val,J12all_val]=f_cost(w_opt);
        % plot([full(J1all_val') full(J1ball_val') full(J2all_val') full(J3all_val') full(J4all_val') full(J5all_val') full(J6all_val') full(J7all_val') full(J8all_val') full(J9all_val') full(J10all_val') full(J11all_val') full(J12all_val')]);
        JJ=full(J1all_val+J1ball_val+J2all_val+J3all_val+J4all_val+J5all_val+J6all_val+J7all_val+J7ball_val+J8all_val+J9all_val+J10all_val+J11all_val+J12all_val);
        J_vals = [ ...
        full(sum(J1all_val)), ...
        full(sum(J1ball_val)), ...
        full(sum(J2all_val)), ...
        full(sum(J3all_val)), ...
        full(sum(J4all_val)), ...
        full(sum(J5all_val)), ...
        full(sum(J6all_val)), ...
        full(sum(J7all_val)), ...
        full(sum(J7ball_val)), ...
        full(sum(J8all_val)), ...
        full(sum(J9all_val)), ...
        full(sum(J10all_val)), ...
        full(sum(J11all_val)), ...
        full(sum(J12all_val))];
        % Labels for the cost components
        labels = {                          
        'Joint Angle (J1)', ...
        'Tracking coordinates velocity (J1b)', ...
        'Track GRF (J2)', ...
        'Track GRF Moment (J3)', ...
        'Track KCF (J4)', ...
        'Track inverse dynamic (J5)', ...
        'Min muscle activation (J6)', ...
        'Min acceleration all dofs except knee sec(J7)', ...
        'Min acceleration knee sec(J7b)', ...
        'Min derivative of muscle activation (J8)', ...
        'Min derivative of FTtilde (J9)', ...
        'Min residual of KCF(J10)', ...
        'Min residual of moment GRF(J11)', ...
        'Min residual of pelvis (J12)'};
        figure;
        bar(J_vals);
        set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 45);
        ylabel('Cost Contribution');
        title('Cost Function Component Contributions');
        grid on;
        %% Visualization in OpenSim GUI
        % Create .mot file for OpenSim GUI
        if writeIKmotion  
            q_opt_GUI = zeros(N,1+nq.all+2);
            q_opt_GUI(:,1) = tgrid(1:end-1)';
            q_opt_GUI(:,2:nq.all+1)  = q_opt_unsc_all.deg(1:N,:);
            q_opt_GUI(:,end-1:end) = 1.51*180/pi*ones(N,2);% pro_sup (locked)
            JointAngle.labels = {'time','pelvis_rz','pelvis_rx','pelvis_ry','pelvis_tx',...
                'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
                'hip_rotation_l','hip_flexion','hip_adduction','hip_rotation',...
                'knee_flexion_l','knee_flexion','knee_adduction','knee_rotation',...
                'knee_tx','knee_ty','knee_tz','ankle_angle_l','ankle_angle',...
                'subtalar_angle_l','subtalar_angle',...
                'lumbar_extension','lumbar_bending','lumbar_rotation',...
                'arm_flex_l','arm_add_l','arm_rot_l',...
                'arm_flex_r','arm_add_r','arm_rot_r',...
                'elbow_flex_l','elbow_flex_r','pro_sup_l','pro_sup_r'};
            JointAngle.data = q_opt_GUI;
            filenameJointAngles = [pathRepo,'/Results/',namescript,...
                    '/IK',savename,'.mot'];
            write_motionFile(JointAngle, filenameJointAngles)
        end 
        
        %% Recalculate pressures for blender visaulization
        facesTib1=csvread(['../../../contactsKneeProsthesis/facesTib1_' num2str(Options.nfacesTib) '.csv']);
        facesTib2=csvread(['../../../contactsKneeProsthesis/facesTib2_' num2str(Options.nfacesTib) '.csv']);
        tib=stlread(['../../../contactsKneeProsthesis/Tibial Insert simpler_superior_' num2str(Options.nfacesTib) '.stl']);
        [m,Itib]=sort(tib.ConnectivityList(:,1));
        facesTib=tib.ConnectivityList(Itib,:);
        % Generate pressure csv
        cd(pathmain);
        
        if Options.nfacesTib==49
            nfacestib1=26;
            nfacestib2=23;
        elseif Options.nfacesTib==100
            nfacestib1=51;
            nfacestib2=49;
        else
            keyboard;
        end
        ndof=size(q_opt_unsc.rad,2);
        for i=1:size(q_col_opt_unsc.rad,1)
           
            
            atest=zeros(ndof,1); %since pressures do not depend on acceleration
            in_F2=[zeros(36+12,1); ones(6,1)]; %since knee pressures do not depend on GRF on this model
            if Options.KCFasinputstoExternalFunction
                xtest=q_col_opt_unsc.rad(i,[jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r jointi.knee_tx.r jointi.knee_ty.r  jointi.knee_tz.r])';
                Tk=F2(xtest);
                Tk_d=F2_debug(xtest);
                CFMed(i)=full(Tk(1));
                CFLat(i)=full(Tk(2));

                CFMed_d(i)=full(Tk_d(1));
                CFLat_d(i)=full(Tk_d(2));
                P1_d(i,:)=full(Tk_d((8+1):(8+nfacestib1)));
                P2_d(i,:)=full(Tk_d((8+nfacestib1+1):(8+nfacestib1+nfacestib2)));
            else
                xtest(1:2:2*ndof,1)=q_col_opt_unsc.rad(i,:)';
                xtest(2:2:2*ndof,1)=qdot_col_opt_unsc.rad(i,:)';
                [Tk] = F2([xtest;atest;in_F2]); %Tk is a vector with 48 elements containing model joint moments GRF and knee contact forces
                [Tk_d] = F2_debug([xtest;atest;in_F2]);
                % Lateral knee contact force is at index 47, medial contact force is at
                % 48
                CFLat(i)=full(Tk(47));
                CFMed(i)=full(Tk(48));
                
                CFLat_d(i)=full(Tk_d(47));
                CFMed_d(i)=full(Tk_d(48));
                P1_d(i,:)=full(Tk_d(48+1:48+nfacestib1));
                P2_d(i,:)=full(Tk_d(48+26+1:48+nfacestib1+nfacestib2));
            end
            P_d(i,:)=reorderpmat(P1_d(i,:),P2_d(i,:),facesTib1,facesTib2,facesTib);
        end
        
        pmat_blender = P_d;
        writematrix(pmat_blender,[pathresults,'/',namescript,'\pressure2bl_matlab',savename,'.csv']);
    
        %export also kinematics to blender
        fTransf = create_fTransf();
        for i=1:size(q_opt_unsc.rad,1)
            Mtransformtib2=full(fTransf(q_opt_unsc.rad(i,17),...
                q_opt_unsc.rad(i,18),q_opt_unsc.rad(i,19),...
                -q_opt_unsc.rad(i,14),-q_opt_unsc.rad(i,15),...
                q_opt_unsc.rad(i,16)));
            agtheta = asin(Mtransformtib2(3,2));
            agphi = asin(Mtransformtib2(3,1)/(-cos(agtheta)));
            agpsi = asin(Mtransformtib2(1,2)/(-cos(agtheta)));
            kin2blender(i,4:6)= [agtheta agphi agpsi]; %Blender uses rad for Euler rotations
            kin2blender(i,1:3)=Mtransformtib2(1:3,4);
        end
    
    %kin2btransposed=kin2blender;
    kin2btransposed=kin2blender';
    writematrix(kin2btransposed,[pathresults,'/',namescript,'\kin2bl_matlab',savename,'.csv']);


        %% Save results       
        if saveResults
        %         hess_names = {'Approximated','Exact'};
                if (exist([pathresults,'/',namescript,...
                        '/Results_3D_v2',savename,'.mat'],'file')==2) 
                    load([pathresults,'/',namescript,'/Results_3D',savename,'.mat']);
                else
                     Results_3D.Simulated = struct('Qs_opt',[]);
                end
            % Structure results
            Results_3D.Simulated.Qs_opt = q_opt_unsc.deg;
            Results_3D.Simulated.Qdots_opt = qdot_opt_unsc.rad;
            Results_3D.Simulated.Qs_col_opt = q_col_opt_unsc    ;
            Results_3D.Simulated.Qdots_col_opt=qdot_col_opt_unsc.rad;
            Results_3D.Simulated.Acts_opt = a_opt_unsc;
            Results_3D.Simulated.Act_col_opt=a_col_opt;  
            Results_3D.Simulated.e_col_opt=e_opt_unsc;
            Results_3D.Simulated.FTtilde_opt_unsc=FTtildep_opt_unsc;
            Results_3D.Simulated.FTtilde_col_opt_unsc=FTtilde_col_opt_unsc;
            Results_3D.Simulated.dFTtilde_col_opt_unsc=dFTtilde_col_opt_unsc;
            Results_3D.Simulated.Ts_opt = Tauk_out;
            Results_3D.Simulated.GRFs_opt = GRF_opt_unsc;
            Results_3D.Simulated.GRFs_opt_k = GRF_opt_unsc_k;
            Results_3D.Simulated.GRMs_opt = GRM_opt_unsc;   
            Results_3D.Simulated.GRMs_opt_k = GRM_opt_unsc_k;
            Results_3D.Simulated.KCFs_opt = KCF_opt_unsc; 
            Results_3D.Simulated.KCFs_opt_k = KCF_opt_unsc_k; 
            Results_3D.Simulated.Qdotsdots_col_opt = qdotdot_col_opt_unsc.rad;
            Results_3D.Simulated.Cost_Fun.J1a = full(J1all_val);
            Results_3D.Simulated.Cost_Fun.J1b = full(J1ball_val);
            Results_3D.Simulated.Cost_Fun.J2 = full(J2all_val);
            Results_3D.Simulated.Cost_Fun.J3 = full(J3all_val);
            Results_3D.Simulated.Cost_Fun.J4 = full(J4all_val);
            Results_3D.Simulated.Cost_Fun.J5 = full(J5all_val);
            Results_3D.Simulated.Cost_Fun.J6 = full(J6all_val);
            Results_3D.Simulated.Cost_Fun.J7 = full(J7all_val);
            Results_3D.Simulated.Cost_Fun.J7b = full(J7ball_val);
            Results_3D.Simulated.Cost_Fun.J8 = full(J8all_val);
            Results_3D.Simulated.Cost_Fun.J9 = full(J9all_val);
            Results_3D.Simulated.Cost_Fun.J10= full(J10all_val);
            Results_3D.Simulated.Cost_Fun.J11= full(J11all_val);
            Results_3D.Simulated.Cost_Fun.J12= full(J12all_val);
            Results_3D.Simulated.Cost_Fun_val = J_vals;
         
                    Results_3D.Simulated.vA_opt_unsc=vA_opt_unsc;
            if Options.useKCFresiduals
                %KCF residuals
                Results_3D.Simulated.KCF_res_opt_unsc=KCF_res_opt_unsc;
            end
            if Options.usePelvisResMom
                Results_3D.Simulated.PelvisRes=Pelvis_res_opt_unsc;
            end
            Results_3D.Simulated.ParamsCM_opt = full(paramsCM_opt_nsc);     
            Results_3D.Simulated.stats = stats; 
            Results_3D.Simulated.tgrid=tgrid;
            Results_3D.Simulated.tgrid_ext=tgrid_ext;
            Results_3D.Simulated.tgrid_col=tgrid_col;
            Results_3D.Simulated.P1_d=P1_d;
            Results_3D.Simulated.P2_d=P2_d;
            Results_3D.Simulated.P_d=P_d;
            Results_3D.Simulated.FTk_opt = FTk_opt;
            Results_3D.Simulated.Tau_pass_opt = Tau_pass_opt;
            Results_3D.Simulated.Fce_opt_aux2 = Fce_opt_aux2;
            Results_3D.Simulated.Fpe_opt_aux2 = Fpe_opt_aux2;
            Results_3D.Simulated.lMtilde_opt_aux2 = lMtilde_opt_aux2;
            Results_3D.Simulated.FMvtilde_aux2 = FMvtilde_aux2;
            
            % Results_3D.Simulated.Tau_passk = Tau_passk;
           
        
            Results_3D.(['NMesh_',num2str(N)]).Qs_toTrack = Qs.allinterpfilt(:,2:end)*180/pi;
            Results_3D.(['NMesh_',num2str(N)]).Qs_toTrack_col = Qs.allinterpfilt_col(:,2:end)*180/pi;
            Results_3D.(['NMesh_',num2str(N)]).Ts_toTrack = ID.allinterp_col;
            Results_3D.(['NMesh_',num2str(N)]).GRFs_toTrack = GRF.val.allinterp_col;
            Results_3D.(['NMesh_',num2str(N)]).GRMs_toTrack = GRF.MorGF.allinterp_col; 
            Results_3D.(['NMesh_',num2str(N)]).KCF_toTrack = KCF.allinterpfilt_col; 
            Results_3D.(['NMesh_',num2str(N)]).KCF_toTrack50 = KCF.allinterpfilt;
            Results_3D.(['NMesh_',num2str(N)]).time_interval=interval;
            Results_3D.ParamsCM_gen = ...
                full(f_nsc18(guess.params,scaling.params.v,scaling.params.r));
            Results_3D.colheaders.joints = joints;
            Results_3D.colheaders.GRF = {'fore_aft_r','vertical_r','lateral_r',...
                'fore_aft_l','vertical_l','lateral_l'};
            for i = 1:NMuscle/2
                    Results_3D.colheaders.muscles{i} = ...
                        [muscleNames{i}(1:end-2),'_l'];
                    Results_3D.colheaders.muscles{i+NMuscle/2} = ...
                        [muscleNames{i}(1:end-2),'_r'];
            end
            Results_3D.colheaders.paramsCM = {'loc_s1_r_x','loc_s1_r_z',...
                'loc_s2_r_x','loc_s2_r_z','loc_s3_r_x','loc_s3_r_z',...
                'loc_s4_r_x','loc_s4_r_z','loc_s5_r_x','loc_s5_r_z',...
                'loc_s6_r_x','loc_s6_r_z','radius_s1','radius_s2','radius_s3',...
                'radius_s4','radius_s5','radius_s6'};
            Options.W=W;
            Results_3D.Options=Options;
            % Save data
            save([pathresults,'/',namescript,'/Results_3D',savename,'.mat'],'Results_3D');
        end
    
        end
    end
end
function   pmat=reorderpmat(pmat1,pmat2,facesTib1,facesTib2,facesTib)
    %define pmat with the proper order to be exported to Blender
    i1=1;
    i2=1;
    
    %check it for face 36 that it is not contained within pairs?
    for i=1:size(facesTib,1)
        if any(sum(facesTib(i,:)==facesTib1,2)==3)
            pmat(i)=pmat1(i1);
            i1=i1+1;
        elseif any(sum(facesTib(i,:)==facesTib2,2)==3)
            pmat(i)=pmat2(i2);
            i2=i2+1;
        else
            keyboard; %it should not hit this keyboard
        end
    end
    
    
    end

    function [Tauk_out, GRF_opt_unsc, GRM_opt_unsc, KCF_opt_unsc]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt,Xk_Qdotdots_opt,F1,F2_in,paramsCM_opt_nsc,loc,loci,stif,radi,dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,toes,calcn,deri,f_contactForce,residualsi,GRFi,GRMi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc,Options,jointi)

if Options.KCFasinputstoExternalFunction
    F2_skeletal=F2_in{1};
    F2=F2_in{2};
else
    F2=F2_in;
end        
            out1_res = F1(Xk_Qs_Qdots_opt);
            out1_res_opt = full(out1_res);
            % Re-organise the contact sphere locations
            locSphere_s1_l_opt = [paramsCM_opt_nsc(loci.s1.r.x),...
                loc.s1.y,-paramsCM_opt_nsc(loci.s1.r.z)]';
            locSphere_s1_r_opt = [paramsCM_opt_nsc(loci.s1.r.x),...
                loc.s1.y,paramsCM_opt_nsc(loci.s1.r.z)]';
            locSphere_s2_l_opt = [paramsCM_opt_nsc(loci.s2.r.x),...
                loc.s2.y,-paramsCM_opt_nsc(loci.s2.r.z)]';
            locSphere_s2_r_opt = [paramsCM_opt_nsc(loci.s2.r.x),...
                loc.s2.y,paramsCM_opt_nsc(loci.s2.r.z)]';
            locSphere_s3_l_opt = [paramsCM_opt_nsc(loci.s3.r.x),...
                loc.s3.y,-paramsCM_opt_nsc(loci.s3.r.z)]';
            locSphere_s3_r_opt = [paramsCM_opt_nsc(loci.s3.r.x),...
                loc.s3.y,paramsCM_opt_nsc(loci.s3.r.z)]';
            locSphere_s4_l_opt = [paramsCM_opt_nsc(loci.s4.r.x),...
                loc.s4.y,-paramsCM_opt_nsc(loci.s4.r.z)]';
            locSphere_s4_r_opt = [paramsCM_opt_nsc(loci.s4.r.x),...
                loc.s4.y,paramsCM_opt_nsc(loci.s4.r.z)]';
            locSphere_s5_l_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
                loc.s5.y,-paramsCM_opt_nsc(loci.s5.r.z)]';
            locSphere_s5_r_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
                loc.s5.y,paramsCM_opt_nsc(loci.s5.r.z)]';
            locSphere_s6_l_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
                loc.s6.y,-paramsCM_opt_nsc(loci.s6.r.z)]';
            locSphere_s6_r_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
                loc.s6.y,paramsCM_opt_nsc(loci.s6.r.z)]';        
            % Compute contact forces
            force_s1_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_l_opt,...
                out1_res_opt(calcn.l.pos),out1_res_opt(calcn.l.v_lin),...
                out1_res_opt(calcn.l.omega),out1_res_opt(calcn.TR.R.l),...
                out1_res_opt(calcn.TR.T.l));
            force_s1_l_opt_all=full(force_s1_l_opt);
            force_s2_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_l_opt,...
                out1_res_opt(calcn.l.pos),out1_res_opt(calcn.l.v_lin),...
                out1_res_opt(calcn.l.omega),out1_res_opt(calcn.TR.R.l),...
                out1_res_opt(calcn.TR.T.l));
            force_s2_l_opt_all(:,:)=full(force_s2_l_opt);
            force_s3_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s3_l_opt,...
                out1_res_opt(calcn.l.pos),out1_res_opt(calcn.l.v_lin),...
                out1_res_opt(calcn.l.omega),out1_res_opt(calcn.TR.R.l),...
                out1_res_opt(calcn.TR.T.l));
            force_s3_l_opt_all=full(force_s3_l_opt);
            force_s4_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s4_l_opt,...
                out1_res_opt(toes.l.pos),out1_res_opt(toes.l.v_lin),...
                out1_res_opt(toes.l.omega),out1_res_opt(toes.TR.R.l),...
                out1_res_opt(toes.TR.T.l));     
            force_s4_l_opt_all(:,:)=full(force_s4_l_opt);
            force_s5_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s5_l_opt,...
                out1_res_opt(calcn.l.pos),out1_res_opt(calcn.l.v_lin),...
                out1_res_opt(calcn.l.omega),out1_res_opt(calcn.TR.R.l),...
                out1_res_opt(calcn.TR.T.l));
            force_s5_l_opt_all(:,:)=full(force_s5_l_opt);
            force_s6_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s6_l_opt,...
                out1_res_opt(toes.l.pos),out1_res_opt(toes.l.v_lin),...
                out1_res_opt(toes.l.omega),out1_res_opt(toes.TR.R.l),...
                out1_res_opt(toes.TR.T.l));      
            force_s6_l_opt_all=full(force_s6_l_opt);
            force_s1_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_r_opt,...
                out1_res_opt(calcn.r.pos),out1_res_opt(calcn.r.v_lin),...
                out1_res_opt(calcn.r.omega),out1_res_opt(calcn.TR.R.r),...
                out1_res_opt(calcn.TR.T.r));
            force_s1_r_opt_all(:,:)=full(force_s1_r_opt);
            force_s2_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_r_opt,...
                out1_res_opt(calcn.r.pos),out1_res_opt(calcn.r.v_lin),...
                out1_res_opt(calcn.r.omega),out1_res_opt(calcn.TR.R.r),...
                out1_res_opt(calcn.TR.T.r));
            force_s2_r_opt_all(:,:)=full(force_s2_r_opt);
            force_s3_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s3_r_opt,...
                out1_res_opt(calcn.r.pos),out1_res_opt(calcn.r.v_lin),...
                out1_res_opt(calcn.r.omega),out1_res_opt(calcn.TR.R.r),...
                out1_res_opt(calcn.TR.T.r));
            force_s3_r_opt_all(:,:)=full(force_s3_r_opt);
            force_s4_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s4_r_opt,...
                out1_res_opt(toes.r.pos),out1_res_opt(toes.r.v_lin),...
                out1_res_opt(toes.r.omega),out1_res_opt(toes.TR.R.r),...
                out1_res_opt(toes.TR.T.r)); 
            force_s4_r_opt_all=full(force_s4_r_opt);
            force_s5_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s5_r_opt,...
                out1_res_opt(calcn.r.pos),out1_res_opt(calcn.r.v_lin),...
                out1_res_opt(calcn.r.omega),out1_res_opt(calcn.TR.R.r),...
                out1_res_opt(calcn.TR.T.r));
            force_s5_r_opt_all=full(force_s5_r_opt);
            force_s6_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s6_r_opt,...
                out1_res_opt(toes.r.pos),out1_res_opt(toes.r.v_lin),...
                out1_res_opt(toes.r.omega),out1_res_opt(toes.TR.R.r),...
                out1_res_opt(toes.TR.T.r)); 
            force_s6_r_opt_all=full(force_s6_r_opt);
            in_F2_opt = [force_s1_l_opt,force_s2_l_opt,force_s3_l_opt,...
                force_s4_l_opt,force_s5_l_opt,force_s6_l_opt,force_s1_r_opt,...
                force_s2_r_opt,force_s3_r_opt,force_s4_r_opt,force_s5_r_opt,...
                force_s6_r_opt,...
                paramsCM_opt_nsc(loci.s1.r.x),...
                paramsCM_opt_nsc(loci.s1.r.z),...
                paramsCM_opt_nsc(loci.s2.r.x),...
                paramsCM_opt_nsc(loci.s2.r.z),...
                paramsCM_opt_nsc(loci.s3.r.x),...
                paramsCM_opt_nsc(loci.s3.r.z),...
                paramsCM_opt_nsc(loci.s4.r.x),...
                paramsCM_opt_nsc(loci.s4.r.z),...
                paramsCM_opt_nsc(loci.s5.r.x),...
                paramsCM_opt_nsc(loci.s5.r.z),...
                paramsCM_opt_nsc(loci.s6.r.x),...
                paramsCM_opt_nsc(loci.s6.r.z),...
                paramsCM_opt_nsc(radi.s1:radi.s6)'];   
            if Options.KCFasinputstoExternalFunction
                q_knee=Xk_Qs_Qdots_opt([[jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r]*2-1]);
                out=F2(q_knee);
                SumForces=full(out(3:5));
                SumMoments=full(out(6:8));
                out2_res = ...
                        F2_skeletal([Xk_Qs_Qdots_opt';Xk_Qdotdots_opt';in_F2_opt';SumForces;SumMoments]); 
            else
                if deri == 2
                    out2_res = F2(Xk_Qs_Qdots_opt,Xk_Qdotdots_opt,in_F2_opt); 
                else
                    out2_res = ...
                        F2([Xk_Qs_Qdots_opt';Xk_Qdotdots_opt';in_F2_opt']);  
                end
            end
            out2_res_opt(:,:) = full(out2_res); 



            % Optimal joint torques, ground reaction forces and moments
            Tauk_out        = out2_res_opt(residualsi);
            GRF_opt_unsc    = out2_res_opt(GRFi.all);
            GRM_opt_unsc    = out2_res_opt(GRMi.all);
            if Options.KCFasinputstoExternalFunction
                KCF_opt_unsc    = full(out(1:2));
            else
                KCF_opt_unsc    = out2_res_opt([KCFi.M KCFi.L]);
            end
            % assertArmTmax should be 0
            assertArmTmax = max(max(abs(out2_res_opt(armsi)-(a_a_opt_unsc')*...
                scaling.ArmTau))); 
            if assertArmTmax > 1*10^(-tol_ipopt)
                disp('Issue when reconstructing residual forces')
            end 

end
