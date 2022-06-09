%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% Update Date: 18-March-2022

clear
close all

%% Part 1: Preprocessing 
%% Compatible with  megconnectome 3.0, fieldtrip-r10442 and MATLAB2012b

%% Please Select!
subjectids = {fill this part};

trialDuration = 2; % 2, 5, or 10 seconds


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directories & paths
dir_subjects = load_path('dir_subjects');

addpath(genpath(load_path('megconnectome-30')));
addpath(load_path('fieldtrip-r10442'));
ft_defaults
addpath(genpath(load_path('workingDir_code')));

%%
cimec = struct;
cimec.dir_subjects = dir_subjects;
cimec.scanids = {'3-Restin', '4-Restin', '5-Restin'};
cimec.trialDuration = trialDuration; 

%% Step 1: clean_data

f = waitbar(0, 'Looping around...');

for idx=1:length(subjectids)
    
    subjectid = subjectids{idx};
    cimec.subjectid = subjectid;
    cimec.experimentid = [ subjectid, '_', 'MEG'];
    cimec.pipelinedatadir = fullfile(dir_subjects, subjectid, 'MEG', 'Restin');
    
    clean_data(cimec)
    
    waitbar(idx/length(subjectids), f, sprintf('Progress: %d %%', floor(idx/length(subjectids)*100)));
    
end
