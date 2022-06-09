%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% This script is for source reconstruction.
%%% Update Date: 18-March-2022


%% Compatible with  brainstorm (Version: 28-May-2021 or later) and MATLAB2020a

%% Output

% computes the scout time series per subject.

%% Inputs
% % Select the directory where the subject folder is
% dir_subjects = cimec.dir_subjects % 'C:\Users\ASUS\Desktop\4- Thesis\HCP\HCP_Subjects';
% % Select the subject
% subjectid = cimec.subjectid % '111514';
% trialDuration = cimec.trialDuration % 6
% cleanProtocol = cimec.cleanProtocol % false
% fsample = cimec.fsample % 300
% resample = cimec.resample; % true
% inverse_measure = cimec.inverse_measure; % 'dspm2018'

%% Function

function scouts_timeSeries(cimec)

tic
%% ===== SUBJECT INFORMATION (*Be Sure Before Running the Script) =====
% Select the directory where the subject folder is
dir_subjects = cimec.dir_subjects; % 'C:\Users\ASUS\Desktop\4- Thesis\HCP\HCP_Subjects'
% Select the subject
subjectid = cimec.subjectid; % '111514'

trialDuration = cimec.trialDuration; % 6
cleanProtocol = cimec.cleanProtocol; % false

fsample = cimec.fsample; % 300
resample = cimec.resample; % true

inverse_measure = cimec.inverse_measure; % 'dspm2018'

%% ===== FILES TO IMPORT =====

% You have to specify the folder in which the tutorial dataset is unzipped
if isempty(dir_subjects) || ~file_exist(dir_subjects)
    error('The first argument must be the full path to the tutorial dataset folder.');
end

% Build the path of the files to import
extended_AnatDir    = fullfile(dir_subjects, subjectid, 'T1w', subjectid);
HCP_anatomy         = fullfile(dir_subjects, subjectid, 'MEG', 'anatomy');

Run1File   = fullfile(dir_subjects, subjectid, 'MEG', 'Restin', 'cleandata', ...
    [subjectid, '_MEG_3-Restin_cleanData_',  sprintf('%ds.mat', trialDuration)]);
Run2File   = fullfile(dir_subjects, subjectid, 'MEG', 'Restin', 'cleandata', ...
    [subjectid, '_MEG_4-Restin_cleanData_',  sprintf('%ds.mat', trialDuration)]);
Run3File   = fullfile(dir_subjects, subjectid, 'MEG', 'Restin', 'cleandata', ...
    [subjectid, '_MEG_5-Restin_cleanData_',  sprintf('%ds.mat', trialDuration)]);

NoiseFile  = fullfile(dir_subjects, subjectid, 'unprocessed', 'MEG', '1-Rnoise', '4D', 'c,rfDC');

% Output path
OutputDir_scouts = fullfile(dir_subjects, 'Outputs', 'scouts_timeSeries');
OutputDir_fsaverage = fullfile(dir_subjects, 'Outputs', 'fsaverage_anatomy');
mkdir(OutputDir_scouts);
mkdir(OutputDir_fsaverage);

% Check if the folder contains the required files
if ~file_exist(extended_AnatDir) || ~file_exist(Run1File) || ~file_exist(NoiseFile)
    error(['The folder ' dir_subjects ' does not contain the selected subject from the HCP-MEG distribution.']);
end

%% ===== CREATE PROTOCOL =====

% The protocol name has to be a valid folder name (no spaces, no weird characters...)
ProtocolName = ['Thesis_OrhanSoyuhos_', date];

% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end

% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
% Start a new report
bst_report('Start');
% Reset colormaps
bst_colormaps('RestoreDefaults', 'meg');

%% ===== GET THE FIDUCIAL POINTS FROM THE old_anatomy =====

% Process: Import the old_anatomy folder to get fiducials
bst_process('CallProcess', 'process_import_anatomy', [], [], ...
    'subjectname', 'old_anatomy', ...
    'mrifile',     {HCP_anatomy, 'HCPv3'}, ...
    'nvertices',   15000);

ProtocolInfo = bst_get('ProtocolInfo');
dirFolder = bst_fullfile(ProtocolInfo.SUBJECTS, 'old_anatomy/subjectimage_T1w_acpc_dc_restore.mat');
sMri_oldAnatomy = load(dirFolder);

% MRI coordinates of fiducials (for old_anatomy)
NAS = sMri_oldAnatomy.SCS.NAS;
LPA = sMri_oldAnatomy.SCS.LPA;
RPA = sMri_oldAnatomy.SCS.RPA;
AC = sMri_oldAnatomy.NCS.AC;
PC = sMri_oldAnatomy.NCS.PC;
IH = sMri_oldAnatomy.NCS.IH;

% Convert the MRI coordinates to World coordinates (for old_anatomy)
nas = cs_convert(sMri_oldAnatomy, 'mri', 'world', NAS ./ 1000) .* 1000;
lpa = cs_convert(sMri_oldAnatomy, 'mri', 'world', LPA ./ 1000) .* 1000;
rpa = cs_convert(sMri_oldAnatomy, 'mri', 'world', RPA ./ 1000) .* 1000;
ac = cs_convert(sMri_oldAnatomy, 'mri', 'world', AC ./ 1000) .* 1000;
pc = cs_convert(sMri_oldAnatomy, 'mri', 'world', PC ./ 1000) .* 1000;
ih = cs_convert(sMri_oldAnatomy, 'mri', 'world', IH ./ 1000) .* 1000;

%% ===== IMPORT THE DATA FREESURFER ANATOMY (new_anatomy) =====

% find the MRI file to import
MriFile = file_find(extended_AnatDir, 'T1.mgz', [], 1);
% Process: Import MRI (Because we need sMri for the Freesurfer anatomy (new_anatomy))
bst_process('CallProcess', 'process_import_mri', [], [], ...
    'subjectname', 'Tmp_Anat', ...
    'mrifile',     {MriFile, 'MGH'}, ...
    'nas',         [0, 0, 0], ...
    'lpa',         [0, 0, 0], ...
    'rpa',         [0, 0, 0], ...
    'ac',          [0, 0, 0], ...
    'pc',          [0, 0, 0], ...
    'ih',          [0, 0, 0]);

dirFolder = bst_fullfile(ProtocolInfo.SUBJECTS, 'Tmp_Anat/subjectimage_T1.mat');
sMri_newAnatomy = load(dirFolder);

% Convert the World coordinates of fidual points from old_anatomy to the MRI coordinates of new_anatomy
final_NAS = cs_convert(sMri_newAnatomy, 'world', 'mri', nas ./ 1000) .* 1000;
final_LPA = cs_convert(sMri_newAnatomy, 'world', 'mri', lpa ./ 1000) .* 1000;
final_RPA = cs_convert(sMri_newAnatomy, 'world', 'mri', rpa ./ 1000) .* 1000;
final_AC = cs_convert(sMri_newAnatomy, 'world', 'mri', ac ./ 1000) .* 1000;
final_PC = cs_convert(sMri_newAnatomy, 'world', 'mri', pc ./ 1000) .* 1000;
final_IH =cs_convert(sMri_newAnatomy, 'world', 'mri', ih ./ 1000) .* 1000;

% Process: Import FreeSurfer folder
bst_process('CallProcess', 'process_import_anatomy', [], [], ...
    'subjectname', subjectid, ...
    'mrifile',     {extended_AnatDir, 'FreeSurfer'}, ...
    'nvertices',   15000, ...
    'nas', final_NAS, ...
    'lpa', final_LPA, ...
    'rpa', final_RPA, ...
    'ac', final_AC, ...
    'pc', final_PC, ...
    'ih', final_IH);


% Get all the subjects in the protocol
sProtocolSubjects = bst_get('ProtocolSubjects');
% Choose the subject
del1_iSubject = find(strcmpi('Tmp_Anat', {sProtocolSubjects.Subject.Name}));
del2_iSubject = find(strcmpi('old_anatomy', {sProtocolSubjects.Subject.Name}));
% Delete the subject
db_delete_subjects(del1_iSubject);
db_delete_subjects(del2_iSubject);
% Reload the Data
db_reload_database('current');

%% ===== IMPORT the CLEAN and EPOCHED DATA  =====

% Process: Import MEG/EEG: Existing epochs
sFileRun1 = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
    'subjectname',  subjectid, ...
    'condition',    sprintf('3-Restin_cleandata_%ds', trialDuration), ...
    'datafile',     {Run1File, 'FT-TIMELOCK'}, ...
    'iepochs',      [], ...
    'eventtypes',   '', ...
    'createcond',   0, ...
    'channelalign', 1, ...
    'usectfcomp',   1, ...
    'usessp',       1, ...
    'freq',         [], ...
    'baseline',     []);
% Process: DC offset correction: [All file]
sFileRun1 = bst_process('CallProcess', 'process_baseline', sFileRun1, [], ...
    'baseline',    [], ...
    'sensortypes', 'MEG, EEG', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'overwrite',   1);
if resample == true
    % Process: Resample: 300Hz
    sFileRun1 = bst_process('CallProcess', 'process_resample', sFileRun1, [], ...
        'freq',      fsample, ...
        'overwrite', 1);
else
    dirFolder = bst_fullfile(ProtocolInfo.STUDIES, sFileRun1(1).FileName);
    tmp_data = load(dirFolder);
    oldRate = 1 ./ (tmp_data.Time(2)-tmp_data.Time(1));
    fsample = oldRate;
end

% Process: Import MEG/EEG: Existing epochs
sFileRun2 = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
    'subjectname',  subjectid, ...
    'condition',    sprintf('4-Restin_cleandata_%ds', trialDuration), ...
    'datafile',     {Run2File, 'FT-TIMELOCK'}, ...
    'iepochs',      [], ...
    'eventtypes',   '', ...
    'createcond',   0, ...
    'channelalign', 1, ...
    'usectfcomp',   1, ...
    'usessp',       1, ...
    'freq',         [], ...
    'baseline',     []);
% Process: DC offset correction: [All file]
sFileRun2 = bst_process('CallProcess', 'process_baseline', sFileRun2, [], ...
    'baseline',    [], ...
    'sensortypes', 'MEG, EEG', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'overwrite',   1);
if resample == true
    % Process: Resample: 300Hz
    sFileRun2 = bst_process('CallProcess', 'process_resample', sFileRun2, [], ...
        'freq',      fsample, ...
        'overwrite', 1);
end

% Process: Import MEG/EEG: Existing epochs
sFileRun3 = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
    'subjectname',  subjectid, ...
    'condition',    sprintf('5-Restin_cleandata_%ds', trialDuration), ...
    'datafile',     {Run3File, 'FT-TIMELOCK'}, ...
    'iepochs',      [], ...
    'eventtypes',   '', ...
    'createcond',   0, ...
    'channelalign', 1, ...
    'usectfcomp',   1, ...
    'usessp',       1, ...
    'freq',         [], ...
    'baseline',     []);
% Process: DC offset correction: [All file]
sFileRun3 = bst_process('CallProcess', 'process_baseline', sFileRun3, [], ...
    'baseline',    [], ...
    'sensortypes', 'MEG, EEG', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'overwrite',   1);
if resample == true
    % Process: Resample: 300Hz
    sFileRun3 = bst_process('CallProcess', 'process_resample', sFileRun3, [], ...
        'freq',      fsample, ...
        'overwrite', 1);
end

%% !!! check the freq sampling and line freq later!!!
sFileNoise = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',  subjectid, ...
    'datafile',     {NoiseFile, '4D'}, ...
    'channelalign', 1);


sFilesClean = [sFileRun1, sFileRun2, sFileRun3, sFileNoise];
sFilesRest = [sFileRun1, sFileRun2, sFileRun3];

%% ===== SOURCE ESTIMATION =====

% Process: Select file names with tag: task-rest
sFilesNoise = bst_process('CallProcess', 'process_select_tag', sFilesClean, [], ...
    'tag',    '1-Rnoise', ...
    'search', 1, ...  % Search the file names
    'select', 1);  % Select only the files with the tag

% Process: Compute covariance (noise or data)
bst_process('CallProcess', 'process_noisecov', sFilesNoise, [], ...
    'baseline',       [], ...
    'sensortypes',    'MEG', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       1, ...
    'copysubj',       0, ...
    'replacefile',    1);  % Replace

% Process: Compute head model
bst_process('CallProcess', 'process_headmodel', sFilesRest, [], ...
    'sourcespace', 1, ...  % Cortex surface
    'meg',         3);     % Overlapping spheres

% Process: Compute sources [2018]
sSrcRest = bst_process('CallProcess', 'process_inverse_2018', sFilesRest, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', struct(...
    'Comment',        'dSPM: MEG', ...
    'InverseMethod',  'minnorm', ...
    'InverseMeasure', inverse_measure, ...
    'SourceOrient',   {{'fixed'}}, ...
    'Loose',          0.2, ...
    'UseDepth',       1, ...
    'WeightExp',      0.5, ...
    'WeightLimit',    10, ...
    'NoiseMethod',    'reg', ...
    'NoiseReg',       0.1, ...
    'SnrMethod',      'fixed', ...
    'SnrRms',         1e-06, ...
    'SnrFixed',       3, ...
    'ComputeKernel',  1, ...
    'DataTypes',      {{'MEG'}}));

% Reload the Data
db_reload_database('current');

%% ===== IMPORT FSAverage =====

% DO NOT FORGET to put the HCP-MMP1.0-fsaverage inside brainstorm directory

% Add fsaverage subject
[sSubject, iSubject] = db_add_subject('fsaverage');

% Add HCP-MMP1.0-fsaverage template
sTemplates = bst_get('AnatomyDefaults');
iTemplate = find(strcmpi('HCP-MMP1.0-fsaverage', {sTemplates.Name}));
db_set_template(iSubject, sTemplates(iTemplate), 0);

%% ===== PROJECTING SOURCES on FSAverage and SAVING SCOUTS TIME SERIES [360 scouts] as Fieldtrip Structure =====

% Get all the subjects in the protocol
sProtocolSubjects = bst_get('ProtocolSubjects');
% Choose the subject
iSubject = find(strcmpi('fsaverage', {sProtocolSubjects.Subject.Name}));
% Get all the cortex surfaces of selected subject
sCortex = bst_get('SurfaceFileByType', iSubject, 'Cortex', 0);
% Choose the cortex
iCortex = find(strcmpi('cortex_15002V', {sCortex.Comment}));

% Choose the scouts for later
ScoutSel = {'HCP-MMP1.0', {'L_10d_ROI L', 'L_10pp_ROI L', 'L_10r_ROI L', 'L_10v_ROI L', 'L_11l_ROI L', 'L_13l_ROI L', 'L_1_ROI L', 'L_23c_ROI L', 'L_23d_ROI L', 'L_24dd_ROI L', 'L_24dv_ROI L', 'L_25_ROI L', 'L_2_ROI L', 'L_31a_ROI L', 'L_31pd_ROI L', 'L_31pv_ROI L', 'L_33pr_ROI L', 'L_3a_ROI L', 'L_3b_ROI L', 'L_43_ROI L', 'L_44_ROI L', 'L_45_ROI L', 'L_46_ROI L', 'L_47l_ROI L', 'L_47m_ROI L', 'L_47s_ROI L', 'L_4_ROI L', 'L_52_ROI L', 'L_55b_ROI L', 'L_5L_ROI L', 'L_5m_ROI L', 'L_5mv_ROI L', 'L_6a_ROI L', 'L_6d_ROI L', 'L_6ma_ROI L', 'L_6mp_ROI L', 'L_6r_ROI L', 'L_6v_ROI L', 'L_7AL_ROI L', 'L_7Am_ROI L', 'L_7PC_ROI L', 'L_7PL_ROI L', 'L_7Pm_ROI L', 'L_7m_ROI L', 'L_8Ad_ROI L', 'L_8Av_ROI L', 'L_8BL_ROI L', 'L_8BM_ROI L', 'L_8C_ROI L', 'L_9-46d_ROI L', 'L_9a_ROI L', 'L_9m_ROI L', 'L_9p_ROI L', 'L_A1_ROI L', 'L_A4_ROI L', 'L_A5_ROI L', 'L_AAIC_ROI L', 'L_AIP_ROI L', 'L_AVI_ROI L', 'L_DVT_ROI L', 'L_EC_ROI L', 'L_FEF_ROI L', 'L_FFC_ROI L', 'L_FOP1_ROI L', 'L_FOP2_ROI L', 'L_FOP3_ROI L', 'L_FOP4_ROI L', 'L_FOP5_ROI L', 'L_FST_ROI L', 'L_H_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', 'L_IFSa_ROI L', 'L_IFSp_ROI L', 'L_IP0_ROI L', 'L_IP1_ROI L', 'L_IP2_ROI L', 'L_IPS1_ROI L', 'L_Ig_ROI L', 'L_LBelt_ROI L', 'L_LIPd_ROI L', 'L_LIPv_ROI L', 'L_LO1_ROI L', 'L_LO2_ROI L', 'L_LO3_ROI L', 'L_MBelt_ROI L', 'L_MIP_ROI L', 'L_MI_ROI L', 'L_MST_ROI L', 'L_MT_ROI L', 'L_OFC_ROI L', 'L_OP1_ROI L', 'L_OP2-3_ROI L', 'L_OP4_ROI L', 'L_PBelt_ROI L', 'L_PCV_ROI L', 'L_PEF_ROI L', 'L_PF_ROI L', 'L_PFcm_ROI L', 'L_PFm_ROI L', 'L_PFop_ROI L', 'L_PFt_ROI L', 'L_PGi_ROI L', 'L_PGp_ROI L', 'L_PGs_ROI L', 'L_PHA1_ROI L', 'L_PHA2_ROI L', 'L_PHA3_ROI L', 'L_PHT_ROI L', 'L_PH_ROI L', 'L_PIT_ROI L', 'L_PI_ROI L', 'L_POS1_ROI L', 'L_POS2_ROI L', 'L_PSL_ROI L', 'L_PeEc_ROI L', 'L_Pir_ROI L', 'L_PoI1_ROI L', 'L_PoI2_ROI L', 'L_PreS_ROI L', 'L_ProS_ROI L', 'L_RI_ROI L', 'L_RSC_ROI L', 'L_SCEF_ROI L', 'L_SFL_ROI L', 'L_STGa_ROI L', 'L_STSda_ROI L', 'L_STSdp_ROI L', 'L_STSva_ROI L', 'L_STSvp_ROI L', 'L_STV_ROI L', 'L_TA2_ROI L', 'L_TE1a_ROI L', 'L_TE1m_ROI L', 'L_TE1p_ROI L', 'L_TE2a_ROI L', 'L_TE2p_ROI L', 'L_TF_ROI L', 'L_TGd_ROI L', 'L_TGv_ROI L', 'L_TPOJ1_ROI L', 'L_TPOJ2_ROI L', 'L_TPOJ3_ROI L', 'L_V1_ROI L', 'L_V2_ROI L', 'L_V3A_ROI L', 'L_V3B_ROI L', 'L_V3CD_ROI L', 'L_V3_ROI L', 'L_V4_ROI L', 'L_V4t_ROI L', 'L_V6A_ROI L', 'L_V6_ROI L', 'L_V7_ROI L', 'L_V8_ROI L', 'L_VIP_ROI L', 'L_VMV1_ROI L', 'L_VMV2_ROI L', 'L_VMV3_ROI L', 'L_VVC_ROI L', 'L_a10p_ROI L', 'L_a24_ROI L', 'L_a24pr_ROI L', 'L_a32pr_ROI L', 'L_a47r_ROI L', 'L_a9-46v_ROI L', 'L_d23ab_ROI L', 'L_d32_ROI L', 'L_i6-8_ROI L', 'L_p10p_ROI L', 'L_p24_ROI L', 'L_p24pr_ROI L', 'L_p32_ROI L', 'L_p32pr_ROI L', 'L_p47r_ROI L', 'L_p9-46v_ROI L', 'L_pOFC_ROI L', 'L_s32_ROI L', 'L_s6-8_ROI L', 'L_v23ab_ROI L', 'R_10d_ROI R', 'R_10pp_ROI R', 'R_10r_ROI R', 'R_10v_ROI R', 'R_11l_ROI R', 'R_13l_ROI R', 'R_1_ROI R', 'R_23c_ROI R', 'R_23d_ROI R', 'R_24dd_ROI R', 'R_24dv_ROI R', 'R_25_ROI R', 'R_2_ROI R', 'R_31a_ROI R', 'R_31pd_ROI R', 'R_31pv_ROI R', 'R_33pr_ROI R', 'R_3a_ROI R', 'R_3b_ROI R', 'R_43_ROI R', 'R_44_ROI R', 'R_45_ROI R', 'R_46_ROI R', 'R_47l_ROI R', 'R_47m_ROI R', 'R_47s_ROI R', 'R_4_ROI R', 'R_52_ROI R', 'R_55b_ROI R', 'R_5L_ROI R', 'R_5m_ROI R', 'R_5mv_ROI R', 'R_6a_ROI R', 'R_6d_ROI R', 'R_6ma_ROI R', 'R_6mp_ROI R', 'R_6r_ROI R', 'R_6v_ROI R', 'R_7AL_ROI R', 'R_7Am_ROI R', 'R_7PC_ROI R', 'R_7PL_ROI R', 'R_7Pm_ROI R', 'R_7m_ROI R', 'R_8Ad_ROI R', 'R_8Av_ROI R', 'R_8BL_ROI R', 'R_8BM_ROI R', 'R_8C_ROI R', 'R_9-46d_ROI R', 'R_9a_ROI R', 'R_9m_ROI R', 'R_9p_ROI R', 'R_A1_ROI R', 'R_A4_ROI R', 'R_A5_ROI R', 'R_AAIC_ROI R', 'R_AIP_ROI R', 'R_AVI_ROI R', 'R_DVT_ROI R', 'R_EC_ROI R', 'R_FEF_ROI R', 'R_FFC_ROI R', 'R_FOP1_ROI R', 'R_FOP2_ROI R', 'R_FOP3_ROI R', 'R_FOP4_ROI R', 'R_FOP5_ROI R', 'R_FST_ROI R', 'R_H_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R', 'R_IFSa_ROI R', 'R_IFSp_ROI R', 'R_IP0_ROI R', 'R_IP1_ROI R', 'R_IP2_ROI R', 'R_IPS1_ROI R', 'R_Ig_ROI R', 'R_LBelt_ROI R', 'R_LIPd_ROI R', 'R_LIPv_ROI R', 'R_LO1_ROI R', 'R_LO2_ROI R', 'R_LO3_ROI R', 'R_MBelt_ROI R', 'R_MIP_ROI R', 'R_MI_ROI R', 'R_MST_ROI R', 'R_MT_ROI R', 'R_OFC_ROI R', 'R_OP1_ROI R', 'R_OP2-3_ROI R', 'R_OP4_ROI R', 'R_PBelt_ROI R', 'R_PCV_ROI R', 'R_PEF_ROI R', 'R_PF_ROI R', 'R_PFcm_ROI R', 'R_PFm_ROI R', 'R_PFop_ROI R', 'R_PFt_ROI R', 'R_PGi_ROI R', 'R_PGp_ROI R', 'R_PGs_ROI R', 'R_PHA1_ROI R', 'R_PHA2_ROI R', 'R_PHA3_ROI R', 'R_PHT_ROI R', 'R_PH_ROI R', 'R_PIT_ROI R', 'R_PI_ROI R', 'R_POS1_ROI R', 'R_POS2_ROI R', 'R_PSL_ROI R', 'R_PeEc_ROI R', 'R_Pir_ROI R', 'R_PoI1_ROI R', 'R_PoI2_ROI R', 'R_PreS_ROI R', 'R_ProS_ROI R', 'R_RI_ROI R', 'R_RSC_ROI R', 'R_SCEF_ROI R', 'R_SFL_ROI R', 'R_STGa_ROI R', 'R_STSda_ROI R', 'R_STSdp_ROI R', 'R_STSva_ROI R', 'R_STSvp_ROI R', 'R_STV_ROI R', 'R_TA2_ROI R', 'R_TE1a_ROI R', 'R_TE1m_ROI R', 'R_TE1p_ROI R', 'R_TE2a_ROI R', 'R_TE2p_ROI R', 'R_TF_ROI R', 'R_TGd_ROI R', 'R_TGv_ROI R', 'R_TPOJ1_ROI R', 'R_TPOJ2_ROI R', 'R_TPOJ3_ROI R', 'R_V1_ROI R', 'R_V2_ROI R', 'R_V3A_ROI R', 'R_V3B_ROI R', 'R_V3CD_ROI R', 'R_V3_ROI R', 'R_V4_ROI R', 'R_V4t_ROI R', 'R_V6A_ROI R', 'R_V6_ROI R', 'R_V7_ROI R', 'R_V8_ROI R', 'R_VIP_ROI R', 'R_VMV1_ROI R', 'R_VMV2_ROI R', 'R_VMV3_ROI R', 'R_VVC_ROI R', 'R_a10p_ROI R', 'R_a24_ROI R', 'R_a24pr_ROI R', 'R_a32pr_ROI R', 'R_a47r_ROI R', 'R_a9-46v_ROI R', 'R_d23ab_ROI R', 'R_d32_ROI R', 'R_i6-8_ROI R', 'R_p10p_ROI R', 'R_p24_ROI R', 'R_p24pr_ROI R', 'R_p32_ROI R', 'R_p32pr_ROI R', 'R_p47r_ROI R', 'R_p9-46v_ROI R', 'R_pOFC_ROI R', 'R_s32_ROI R', 'R_s6-8_ROI R', 'R_v23ab_ROI R'}};

% empty Fieldtrip data to export later
FT_data = struct;


%%% 3-Restin
% Process: Select Sources in 3-Restin_cleandata_
Sources_1 = bst_process('CallProcess', 'process_select_files_results', [], [], ...
    'subjectname', subjectid, ...
    'condition',   sprintf('3-Restin_cleandata_%ds', trialDuration), ...
    'includebad',  0);

length_s_1 = sum(~cellfun(@isempty,{Sources_1.FileName}));
for idx = 1 : length_s_1
    % Projecting sources on FSAverage
    ResultsFile = cellstr(Sources_1(idx).FileName);
    Projection = bst_project_sources(ResultsFile, sCortex(iCortex).FileName);
    
    % Converting to Fieldtrip structure
    ScoutFunc = 'Max';
    isNorm = 0;
    [ftData, sInput, VertConn] = out_fieldtrip_results( Projection{1}, ScoutSel, ScoutFunc, [], isNorm);
    
    FT_data.trial{idx} = ftData.pow;
    FT_data.trialinfo(idx) = {sInput.Comment};
    
    display(FT_data.trialinfo(idx))
    
    % Process: Delete selected files
    bst_process('CallProcess', 'process_delete', Projection, [], ...
        'target', 1);  % Delete selected files 
end


%%% 4-Restin
% Process: Select Sources in 4-Restin_cleandata_
Sources_2 = bst_process('CallProcess', 'process_select_files_results', [], [], ...
    'subjectname', subjectid, ...
    'condition',   sprintf('4-Restin_cleandata_%ds', trialDuration), ...
    'includebad',  0);

length_s_2 = sum(~cellfun(@isempty,{Sources_2.FileName}));
for idx = 1+length_s_1 : length_s_1+length_s_2
    % Projecting sources on FSAverage
    ResultsFile = cellstr(Sources_2(idx-length_s_1).FileName);
    Projection = bst_project_sources(ResultsFile, sCortex(iCortex).FileName);
    
    % Converting to Fieldtrip structure
    ScoutFunc = 2; % Max
    isNorm = 0;
    [ftData, sInput, VertConn] = out_fieldtrip_results( Projection{1}, ScoutSel, ScoutFunc, [], isNorm);
    
    FT_data.trial{idx} = ftData.pow;
    FT_data.trialinfo(idx) = {sInput.Comment};
    
    display(FT_data.trialinfo(idx))
    
    
    % Process: Delete selected files
    bst_process('CallProcess', 'process_delete', Projection, [], ...
        'target', 1);  % Delete selected files
end


%%% 5-Restin
% Process: Select Sources in 5-Restin_cleandata_
Sources_3 = bst_process('CallProcess', 'process_select_files_results', [], [], ...
    'subjectname', subjectid, ...
    'condition',   sprintf('5-Restin_cleandata_%ds', trialDuration), ...
    'includebad',  0);

length_s_3 = sum(~cellfun(@isempty,{Sources_3.FileName}));
for idx = 1+length_s_1+length_s_2 : length_s_1+length_s_2+length_s_3
    % Projecting sources on FSAverage
    ResultsFile = cellstr(Sources_3(idx-length_s_1-length_s_2).FileName);
    Projection = bst_project_sources(ResultsFile, sCortex(iCortex).FileName);
    
    % Converting to Fieldtrip structure
    ScoutFunc = 2; % Max
    isNorm = 0;
    [ftData, sInput, VertConn] = out_fieldtrip_results( Projection{1}, ScoutSel, ScoutFunc, [], isNorm);
    
    FT_data.trial{idx} = ftData.pow;
    FT_data.trialinfo(idx) = {sInput.Comment};
    
    display(FT_data.trialinfo(idx))
    
    % Process: Delete selected files
    bst_process('CallProcess', 'process_delete', Projection, [], ...
        'target', 1);  % Delete selected files  
end


time = cell(1, length(FT_data.trial));
time(:)={ftData.time};
FT_data.time = time;

FT_data.label = sInput.RowNames;
FT_data.fsample	= fsample;

FT_data.trialDuration = trialDuration;


% Reload the Data
db_reload_database('current');

%% ===== WRITE THE SUBJECT STRUCTURE TO DISK =====

% % Make it memory efficent
% FT_data = ft_struct2single(FT_data);
% fsaverage = ft_struct2single(fsaverage);

save(fullfile(OutputDir_scouts, [subjectid, '_scoutsTimeSeries_', sprintf('%ds.mat', trialDuration)]), 'FT_data', '-v7.3');

filename = fullfile(OutputDir_fsaverage, ['fsaverage_', sprintf('%ds.mat', trialDuration)]);
if ~isfile(filename)
    % Save information related to anatomy
    fsaverage = in_tess_bst( bst_fullfile(ProtocolInfo.SUBJECTS, sCortex(iCortex).FileName), 0);
    save(filename, 'fsaverage', '-v7.3');
end

%% ===== THE END =====
disp(['THE END for ', subjectid, '!'])

if cleanProtocol == true
    gui_brainstorm('DeleteProtocol', ProtocolName);
end

toc
end
