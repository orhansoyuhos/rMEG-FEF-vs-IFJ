%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Update Date: 18-March-2022

clear
close all

%% Part 2: END-to-END script after Preprocessing
%% Compatible with  brainstorm (Version: 28-May-2021 or later), fieldtrip-20210411 and MATLAB2020a

%% Please Select!
% Also, please fill the load_path function inside helper_fun folder for correct directories

seedRegions = {'L_FEF_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', ...
    'R_FEF_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R'}; % seed-based connectivity

trialDuration = 2; % 2, 5, or 10 seconds
connMetric = 'oPEC'; % iCOH, oPEC, dwPLI, or PDC


subjectids = {fill this part};


%% fixed
step2 = true; %% Step 2: scouts_timeSeries
step3 = true; %% Step 3: functional_connectivity
step4 = true; %% Step 4: results_connectivity and group_average

cleanProtocol = true; % meaning delete the protocol you created in Brainstorm
redefine = true; % to fix the floating numbers for cfg.foi
resample = true; fsample = 300;
inverse_measure = 'dspm2018'; % 'dspm2018'
top_k = 5; % save the names of highest functional couplings in Excel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% paths
addpath(genpath(fullfile(pwd, 'helper_fun')));
addpath(genpath(load_path('brainstorm')));
addpath(load_path('fieldtrip-20210411'));
ft_defaults
addpath(genpath(load_path('workingDir_code')));

% directories
dir_subjects = load_path('dir_subjects');
scoutsTimeSeries = load_path('scoutsTimeSeries');
functionalConnectivity = fullfile(load_path('functionalConnectivity'), [connMetric '_' char(string(trialDuration)) 's']);
connectivityResults = fullfile(load_path('connectivityResults'), [connMetric '_' char(string(trialDuration)) 's']);
groupResults = fullfile(load_path('groupResults'), [connMetric '_' char(string(trialDuration)) 's']);


%% Step 2: scouts_timeSeries
if step2 == true
    cimec = struct;
    cimec.dir_subjects = dir_subjects;
    cimec.trialDuration = trialDuration;
    cimec.cleanProtocol = cleanProtocol;
    cimec.fsample = fsample;
    cimec.resample = resample;
    cimec.inverse_measure = inverse_measure;
    
    % Loop over
    f = waitbar(0, 'Looping around...');
    
    for idx=1:length(subjectids)
        
        subjectid = subjectids{idx};
        cimec.subjectid = subjectid;
        
        scouts_timeSeries(cimec);
        
        waitbar(idx/length(subjectids), f, sprintf('Progress: %d %%', floor(idx/length(subjectids)*100)));
    end
    
end
%% Step 3: functional_connectivity
if step3 == true
    cimec = struct;
    cimec.scoutsTimeSeries = scoutsTimeSeries;
    cimec.seedRegions = seedRegions;
    cimec.trialDuration = trialDuration;
    cimec.redefine = redefine;
    cimec.conn_metric = connMetric;
    
    % Loop over
    f = waitbar(0, 'Looping around...');
    
    for idx=1:length(subjectids)
        
        subjectid = subjectids{idx};
        outputdir = fullfile(functionalConnectivity, subjectid);
        mkdir(outputdir);
        
        cimec.subjectid = subjectid;
        cimec.outputdir = outputdir;
        
        functional_connectivity(cimec);
        
        waitbar(idx/length(subjectids), f, sprintf('Progress: %d %%', floor(idx/length(subjectids)*100)));
    end
    
end
%% Step 4A: results_connectivity and group_average_topK_forSelectedROIs (for iCOH and dwPLI; Fieldtrip computes seed-based connectivity for these measures)
if step4 == true
    
    if or(strcmp(connMetric, 'iCOH'), strcmp(connMetric, 'dwPLI'))
        cimec = struct;
        cimec.outputdir = connectivityResults;
        cimec.top_k = top_k;
        cimec.trialDuration = trialDuration;
        
        % Loop over
        f = waitbar(0, 'Looping around...');
        
        for idx=1:length(subjectids)
            
            subjectid = subjectids{idx};
            inputdir = fullfile(functionalConnectivity, subjectid);
            
            cimec.subjectid = subjectid;
            cimec.functionalConnectivity = inputdir;
            
            results_connectivity(cimec);
            
            waitbar(idx/length(subjectids), f, sprintf('Progress: %d %%', floor(idx/length(subjectids)*100)));
            
        end
        
        cimec = struct;
        cimec.subjectids = subjectids;
        cimec.top_k = top_k;
        cimec.trialDuration = trialDuration;
        cimec.connectivityResults = connectivityResults;
        cimec.outputdir = groupResults;
        
        group_average_topK_forSelectedROIs(cimec);
        
%% Step 4B: group_average_topK_wholeBrainConnectivity (oPEC and PDC; Fieldtrip computes the whole-brain connectivity map for these measures)
    elseif or(strcmp(connMetric, 'oPEC'), strcmp(connMetric, 'PDC'))
        cimec = struct;
        cimec.subjectids = subjectids;
        cimec.top_k = top_k;
        cimec.trialDuration = trialDuration;
        cimec.functionalConnectivity = functionalConnectivity;
        cimec.outputdir = groupResults;
        
        group_average_topK_wholeConnectivity(cimec)
    end
    
end
%%
