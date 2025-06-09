%% Settings etc.
clear; restoredefaultpath;

root      = 'D:\NIMADET';
dataDir   = fullfile(root,'Data');
resultDir = fullfile(root,'Results');
spm_dir   = 'D:\spm12'; % SPM

addpath(genpath(fullfile(root,'Analyses')));
addpath(genpath(spm_dir)); 

subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11',...
    'S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24',...
    'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34','S35','S36'};
nSubs    = length(subjects);

% exclusions:
% S04 - no variation in vividness rating
% S08 - failed imagery check 11/16 blocks
% S10 - no variation in vividness rating
% S17 - failed half of imagery checks - should be 75% correct at least 
% S19 - no variation in vividness rating
% S22 - no variation in vividness rating
% S24 - too much movement (> 2x voxel size)
% S25 - no variation in vividness rating
% S29 - only completed two blocks
% S30 - too much movement (> 2x voxel size) 
% [S33 - removed only first run - falling asleep]

which_subjects = [1:3 5:7 9 11:16 18 20 21 23 26:28 31:36];

%% 1. Behaviour
% 1.1 main behavioural effects
cfg         = [];
cfg.root    = root;
cfg.subjects = subjects(which_subjects);
cfg.plotIndividual = false; 
cfg.plotting = true;

[Cr,Dp,FA,H,V,acc] = BehaviourAnalysis(cfg); 

% 1.2 model fitting
cfg         = [];
cfg.root    = root;
cfg.subjects = subjects(which_subjects);
cfg.outDir  = 'Simulations/ModelFit';
cfg.plotGOF = true; % Goodness of Fit
cfg.plotMF  = true; % Model Fit on data

ModelFittingBehaviour(cfg);

%% 2. Pre-processing MRI 

% 2.1 Add details to fil_subject_details then run
fil_subject_details

% 2.2 process fieldmaps
fil_fieldmap_preprocess(which_subjects,subj)

% 2.3 spatial preprocessing of task and functional localizer data 
dir_base    = dataDir;
dir_spm     = spm_dir;
n_sess      = 6; % including func localizers
TR          = 1.45;
nslices     = 72;
slicetiming = 0;
resolEPI    = [2 2 2];

for i = 1:length(which_subjects)    
    
    spatial_preprocessingBatch('dir_base', dir_base,'dir_spm', dir_spm, ...
        'subjects', subjects(which_subjects(i)),'n_sess', n_sess,'TR', TR,...
        'nslices', nslices,'slicetiming', slicetiming,'resolEPI',resolEPI)
end

%% 3. First-level analyses 

% 3.1 Create regressors
cfg         = [];
cfg.root    = dataDir;
cfg.out     = resultDir;
cfg.subjects = subjects(which_subjects);

CreateTaskRegressors(cfg);
CreateBehModRegressors(cfg); % make the pmod regressors based on responses

% 3.2 Estimate first-level GLMs
% 3.2a task conditions
cfg           = [];
cfg.TR        = 1.45;
cfg.root      = root;
cfg.dataDir   = dataDir;
cfg.scanDir   = 'Functional';
cfg.prefix    = 'swuf*';
cfg.subjects  = subjects(which_subjects);
cfg.outDir    = 'FirstLevel/MT_taskfactors';
cfg.regs      = 'Regressors/MT_taskFactors';
cfg.mask      = fullfile(dataDir,'rgm_mask.nii');
cfg.spm_dir   = spm_dir;

FirstLevelGLM(cfg)

% 3.2b behavioural modulators
cfg.outDir    = 'FirstLevel/MT_BehModRegressors';
cfg.regs      = 'Regressors/BehModRegressors';
cfg.spm_dir   = spm_dir;

FirstLevelGLMcat(cfg)

% 3.3 Contrasts
% 3.3a Condition main effects
cfg = [];
cfg.root = resultDir;
cfg.FLdir = 'FirstLevel/MT_taskfactors';
cfg.subjects = subjects(which_subjects);
cfg.spm_dir  = spm_dir;
cfg.tConNames = {'Congruency','Stimulus','Interaction'};
cfg.plotCon   = false;
cfg.tCons{1} = {'cong_','inco_'};  
cfg.tCons{2} = {'abs','pres'};
cfg.tCons{3} = {'congruency x stimulus'};

FirstLevelContrasts(cfg);

% 3.3b Behavioural modulation
cfg.FLdir = 'FirstLevel/MT_BehModRegressors';
cfg.tConNames = {'Vividness','RealityJudgement'};
cfg = rmfield(cfg,'tCons');
cfg.tCons{1} = zeros(46,1); cfg.tCons{4}([3,9,15,21]) = 1; 
cfg.tCons{2} = zeros(46,1); cfg.tCons{5}([5,11,17,23]) = 1; 
FirstLevelContrasts(cfg);

%% 4. Single-trials betas for decoding

% 4.1 Get single trial betas 
% 4.1a From localizer
cfg          = [];
cfg.spm_dir  = spm_dir;
cfg.TR       = 1.45;
cfg.data     = dataDir;
cfg.results  = resultDir;
cfg.subjects = subjects(which_subjects);
cfg.prefix   = 'swuf'; % which data

cfg.type     = 'Localizer';
cfg.sess     = [5 6]; % which functional sessions 
BetasPerTrialPoldrack(cfg)

% 4.1b from task
cfg.type     = 'MT';
cfg.sess     = [1 2 3 4]; % which functional sessions 
BetasPerTrialPoldrack(cfg)

% 4.2 get a t-map for stimulus on versus off to pre-select voxels
cfg          = [];
cfg.spm_dir  = spm_dir;
cfg.TR       = TR;
cfg.results  = resultDir;
cfg.data     = dataDir;
cfg.outDir   = 'FirstLevel/Stimulation';
cfg.sess     = [5 6]; % which functional sessions 
cfg.prefix   = 'swuf'; % which data
cfg.subjects = subjects(which_subjects);

StimulationFirstLevel(cfg)

%% 5. Second-level analyses

% 5.1 Second level one-sample t-test
cfg = [];
cfg.root      = fullfile(root,'Results');
cfg.subjects  = subjects(which_subjects);
cfg.spm_dir   = spm_dir;
cfg.mask      = fullfile(root,'resliced_gm_mask.nii');

% 5.1a Task conditions
cfg.dir       = 'FirstLevel/MT_taskfactors';
for t = 1:3
    cfg.conName   = sprintf('con_%04d.nii',t);
    cfg.outputDir = fullfile(cfg.root,'GroupResults',cfg.dir,cfg.conName);
    SecondLevelOST(cfg)
end

% 5.1b Behavioural modulators
cfg.dir       = 'FirstLevel/MT_BehModRegressors';
for t = 1:2
    cfg.conName   = sprintf('con_%04d.nii',t);
    cfg.outputDir = fullfile(cfg.root,'GroupResults',cfg.dir,cfg.conName);
    SecondLevelOST(cfg)
end

% 5.2 Make sig cluster map
cfg = [];
cfg.root      = fullfile(root,'Results');
cfg.outDir    = 'GroupResults\FirstLevel\MT_BehModRegressorsIC\CongFA';
cfg.prefix    = 'main';
cfg.tmap      = 1;
MakeSigClusterMap(cfg)

%% 6. ROI analyses

% 6.1 Extract RS from model
cfg         = [];
cfg.root    = root;
cfg.subjects = subjects(which_subjects);
cfg.outDir  = 'Simulations/ModelFit';

ExtractModelRS(cfg);

% 6.2 Main effects of condition
cfg = [];
cfg.root = root;
cfg.plot = true;
cfg.subjects = subjects(which_subjects);
cfg.ROI = fullfile(root,'GroupResults','MT_BehModRegressors','con_0002','leftFG.nii');
[~,cfg.ROI_name] = fileparts(cfg.ROI);

cfg.data_dir = 'FirstLevel\MT_taskfactors';
cfg.cond{1}  = 'cong_abs';
cfg.cond{2}  = 'inco_abs';
cfg.cond{3}  = 'cong_pres';
cfg.cond{4}  = 'inco_pres';

ROItaskfactoranalysis(cfg);

% 6.3 Effects of vividness and reality judgement
cfg = [];
cfg.root = root;
cfg.subjects = subjects(which_subjects);
cfg.ROI      = fullfile(root,'GroupResults','MT_BehModRegressors','con_0002','leftFG.nii');
cfg.data_dir = 'FirstLevel/MT_perTrialPoldrack';
cfg.beh_dir  = 'Regressors/Behaviour_matrix';

cfg.cond{1}  = 'R==0 & V<(median(V(R==0)))';
cfg.cond{2}  = 'R==0 & V>(median(V(R==0)))';
cfg.cond{3}  = 'R==1 & V<(median(V(R==1)))';
cfg.cond{4}  = 'R==1 & V>(median(V(R==1)))';

ROIsingletrialanalysis(cfg);

% 6.4 Imagery perception confusions
cfg = [];
cfg.root = root;
cfg.plot = true;
cfg.subjects = subjects(which_subjects);
cfg.ROI = fullfile(root,'GroupResults','MT_BehModRegressors','con_0002','leftFG.nii');
[~,cfg.ROI_name] = fileparts(cfg.ROI);

cfg.data_dir = 'FirstLevel\MT_BehModRegressors';
cfg.cond{1}  = 'cong_absxrealityjudgement'; 
cfg.cond{2}  = 'inco_absxrealityjudgement'; 
ROItaskfactoranalysis(cfg);

%% 7. Decoding analyses
% 7.1 Localizer to task decoding
cfg          = []; 
cfg.root     = root;
cfg.ima_check = true;
cfg.subjects = subjects(which_subjects);

cfg.ROI      = fullfile(cfg.root,'\Results\GroupResults\FirstLevel\MT_taskfactors\con_0001\MainEffectEVC.nii');
cfg.data{1}  = 'FirstLevel/MT_perTrialPoldrack';
cfg.data{2}  = 'FirstLevel/Localizer_perTrialPoldrack';
cfg.stim     = 'FirstLevel/Stimulation/spmT_0001.nii';

cfg.permutation = true;
cfg.plotIndividual = false;
cfg.plotGroup   = true;
cfg.gamma = 0.1;
cfg.outDir   = 'Decoding/LocalizerTask_EVC';
LocalizerTaskDecodingROIs(cfg);

% 7.2 Vividness to reality judgement decoding
cfg          = []; 
cfg.root     = root;
cfg.atlas    = fullfile(cfg.root,'\Results\GroupResults\FirstLevel\MT_BehModRegressors\con_0002\leftFG.nii');
cfg.ima_check = true;
cfg.subjects = subjects(which_subjects);
cfg.dataDir  = 'FirstLevel/MT_perTrialPoldrack';
cfg.outDir   = 'Decoding/VividnessResponse_lFG';
cfg.permutation = true;
cfg.plotIndividual = false;
cfg.plotGroup   = true;
cfg.gamma = 0.1;

VividnessResponseDecodingROIs(cfg)

%% 8. Binary decision network

% 8.1 Reality judgement > vividness
cfg = [];
cfg.root = resultDir;
cfg.FLdir = 'FirstLevel/MT_BehModRegressors';
cfg.subjects = subjects(which_subjects);
cfg.spm_dir  = spm_dir;
cfg.tConNames = {'RJ>Vividness'};
cfg.plotCon   = false;
cfg.tCons{1} = zeros(46,1); cfg.tCons{1}([3,9,15,21]) = -1;
cfg.tCons{1}([5,11,17,23]) = 1; 
FirstLevelContrasts(cfg);

cfg.mask      = fullfile(root,'resliced_gm_mask.nii');
cfg.dir       = 'FirstLevel/MT_BehModRegressors';
cfg.conName   = 'con_0003.nii';
cfg.outputDir = fullfile(cfg.root,'GroupResults',cfg.dir,cfg.conName);
SecondLevelOST(cfg)

% 8.2 Extract RS / threshold relationship model
cfg         = [];
cfg.root    = root;
cfg.subjects = subjects(which_subjects);
cfg.outDir  = 'Simulations/ModelFit';

ModelRSThreshold(cfg);

% 8.3 Functional coupling
% Beta series correlation
cfg = [];
cfg.root     = root;
cfg.subjects = subjects(which_subjects);
cfg.data     = 'FirstLevel/MT_perTrialPoldrack';
FLdir        = fullfile(resultDir,'\GroupResults\FirstLevel\MT_BehModRegressors');

cfg.ROI{1} = fullfile(FLdir,'\con_0002\leftFG.nii');

cfg.ROI{2} = fullfile(FLdir,'\con_0002\leftCaudate.nii');
cfg.ROI{3} = fullfile(FLdir,'\con_0002\leftInsula.nii');
cfg.ROI{4} = fullfile(FLdir,'\con_0002\preSMA.nii');

cfg.outDir = 'GroupResults/FunctionalCoupling/leftFG';
FunctionalCoupling(cfg)

%% 9. Multivariate signatures of RS signal

% 9.1 Decode vividness per searchlight
cfg          = []; 
cfg.root     = root;
cfg.data     = resultDir;
cfg.ima_check = true;
cfg.subjects = subjects(which_subjects);
cfg.dataDir  = '_perTrialPoldrack';
cfg.outDir   = 'Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000';
cfg.radius   = 4;
cfg.func_mask  = fullfile(root,'Results\GroupResults','FirstLevel','GLM_cong','con_0001','mask.nii');

VividnessRegressionDecodingSearchlight(cfg);

% 9.2 Second level one-sample t-test on spearman correlations
cfg = [];
cfg.root      = fullfile(root,'Results');
cfg.subjects  = subjects(which_subjects);
cfg.dir       = 'Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000';
cfg.spm_dir   = spm_dir;
cfg.contrast  = 'accuracy.nii';
cfg.outputDir = fullfile(cfg.root,'GroupResults',cfg.dir,'accuracy');
cfg.mask      = fullfile(root,'resliced_gm_mask.nii');
%cfg.covariate = fullfile(root,'Results','GroupResults','criterion_shift.mat');

SecondLevelOST(cfg)

% 9.3 Do follow up tests on significant searchlights to find RS signatures 
cfg = [];
cfg.root      = fullfile(root,'Results');
cfg.data     = resultDir;
cfg.ima_check = true;
cfg.subjects = subjects(which_subjects);
cfg.dataDir  = '_perTrialPoldrack';
cfg.dir       = 'Revision1_Neuron/Decoding/Vividness_SL_L2_100000';
cfg.tthreshold = 3.450189; % p < 0.001 uncorrected
cfg.mask      = fullfile(root,'resliced_gm_mask.nii');
cfg.radius   = 4;
cfg.func_mask  = fullfile(root,'Results\GroupResults','FirstLevel','GLM_cong','con_0001','mask.nii');

CheckRSpattern(cfg)

% 9.4 Get the significant multivariate RS ROIs
cfg = [];
cfg.root       = root;
cfg.atlas      = 'D:\NIMADET\Analyses\ROI_masks\aal';
cfg.atlasNames = {'Insula_L'};
cfg.tmap       = '\Results\GroupResults\Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000\RS_sig.nii';
cfg.tthresh    = 0.1;
cfg.write      = '\Results\GroupResults\Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000\leftInsula.nii';

GetROIindices(cfg);

%% 10. ROI analyses within mv RS regions
% 10.1 Univariate effects 
% ROI univariate effects
cfg = [];
cfg.root = root;
cfg.subjects = subjects(which_subjects);
cfg.dir = 'D:\NIMADET\Results\GroupResults\Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000';
cfg.ROIs = {'leftInsula','rightInsula','dmPFC'};
cfg.data_dir = 'FirstLevel\MT_BehModRegressorsIC';
cfg.contrast{1}  = 'spmT_0004.nii'; 
cfg.contrast{2}  = 'spmT_0005.nii'; 
cfg.conNames     = {'Vividness','Detection'};
cfg.plotting     = true;

act = ROIcontrastsMultipleROIs(cfg);

cfg = [];
cfg.root = root;
cfg.plot = true;
cfg.subjects = subjects(which_subjects);
cfg.ROI = 'D:\NIMADET\Results\GroupResults\Revision1_Neuron/Decoding/Vividness_SL_L2_100000\midOcc.nii';
[~,cfg.ROI_name] = fileparts(cfg.ROI);

cfg.data_dir = 'FirstLevel\MT_TaskRegressorsGN';
cfg.cond{1}  = 'inco_abs';
cfg.cond{2}  = 'inco_pres';
cfg.cond{3}  = 'cong_abs';
cfg.cond{4}  = 'cong_pres';

act = ROItaskfactoranalysis(cfg);

% 10.2 Predicted RS per ROI
% Look at predicted RS for cong, pres, and high viv, low viv
cfg = [];
cfg.root      = fullfile(root,'Results');
cfg.data      = resultDir;
cfg.ima_check = true;
cfg.subjects  = subjects(which_subjects);
cfg.dataDir   = '_perTrialPoldrack';
cfg.dir       = 'Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000';
cfg.ROIs      = {'leftInsula','rightInsula','dmPFC'};

RSPredictedAct(cfg)

% 10.3 Check if vividness or RJ is better decodeable from these ROIs
cfg          = []; 
cfg.root     = root;
cfg.dir      = 'Revision1_Neuron/Decoding/Vividness_Corrected_SL_L2_100000';
cfg.ROIs     = {'leftInsula','rightInsula','dmPFC'};
cfg.ima_check = true;
cfg.subjects = subjects(which_subjects);
cfg.dataDir  = 'FirstLevel/MT_perTrialPoldrack';
cfg.outDir   = 'Revision1_Neuron/Decoding/VividnessCorrectedVSRJ';
cfg.permutation = true;
cfg.plotIndividual = false;
cfg.plotGroup   = true;
cfg.gamma = 0.1;
cfg.balancetrialsRJViv = false; % balance the trials per RJ/Viv

DecodeVividnessVersusRJDSsamples(cfg)


