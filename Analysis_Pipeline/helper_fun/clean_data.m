%% Project Title: Functional Connectivity Fingerprints of Frontal Eye Field and Inferior Frontal Junction
%%% This script for the preprocessing step.
%%% Update Date: 18-March-2022


%% Compatible with  megconnectome 3.0, fieldtrip-r10442 and MATLAB2012b

%% This function was modified from the parts of megconnectome (3.0) to use my project at CIMEC, University of Trento.

% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)

%% Output

% computes the clean data.


%% Inputs
% %% Info about the data
% subjectid = cimec.subjectid; % e.g. '111514'
% experimentid = cimec.experimentid; % e.g. '111514_MEG'
% scanids = cimec.scanids; % e.g. {'3-Restin', '4-Restin', '5-Restin'}
% dir_subjects = cimec.dir_subjects; % e.g. 'C:\Users\ASUS\Desktop\4- Thesis\HCP\HCP_Subjects'
% pipelinedatadir = cimec.pipelinedatadir; % e.g. 'C:\Users\ASUS\Desktop\Restin'
% trialDuration = cimec.trialDuration; % e.g. 10 (seconds)

%% Function

function clean_data(cimec)

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

%% Info about the data
subjectid = cimec.subjectid; % e.g. '111514'
experimentid = cimec.experimentid; % e.g. '111514_MEG'
scanids = cimec.scanids; % e.g. {'3-Restin', '4-Restin', '5-Restin'}
dir_subjects = cimec.dir_subjects; % e.g. 'C:\Users\ASUS\Desktop\4- Thesis\HCP\HCP_Subjects'


pipelinedatadir = cimec.pipelinedatadir; % e.g. 'C:\Users\ASUS\Desktop\4- Thesis\HCP\HCP_Subjects\175237\MEG\Restin'
% change to the location of the processed data (input and output)
cd(pipelinedatadir)

trialDuration = cimec.trialDuration; % e.g. 10 (seconds)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_baddata = fullfile(pipelinedatadir, 'baddata', '/');
dir_icaclass = fullfile(pipelinedatadir, 'icaclass', '/');
dir_cleandata = fullfile(pipelinedatadir, 'cleandata', '/');
mkdir(dir_cleandata);


for idx=1:length(scanids)
    
    raw_data = fullfile(dir_subjects, subjectid, 'unprocessed', 'MEG', scanids{idx}, '4D', 'c,rfDC');


    resultprefix  = sprintf('%s_%s', experimentid, scanids{idx});
    badchansuffix = 'baddata_badchannels';  % mnemonic for file that contains bad channel and segment info
    badsegmsuffix = 'baddata_badsegments';  % mnemonic for file that contains bad channel and segment info
    icainfosuffix = 'icaclass_vs';          % mnemonic to save results
    
    cfg                       = [];
    cfg.dataset               = raw_data;
    cfg.trialfun              = 'trialfun_Restin';
    cfg.trialdef.trialDuration = trialDuration;
    cfg                       = ft_definetrial(cfg);
    cfg.trl(:,4)              = nan; % hcp_extract_allfromrun assumes a trigger code to be present
    
    cfg.badchanfile           = [dir_baddata, resultprefix, '_', badchansuffix, '.txt'];
    cfg.badsegmfile           = [dir_baddata, resultprefix, '_', badsegmsuffix, '.txt'];
    cfg.icainfofile           = [dir_icaclass, resultprefix, '_', icainfosuffix];
    cfg.montage               = hcp_exgmontage(subjectid, experimentid, scanids{idx});
    cfg.lineFreq              = [60 120];
    cfg.badsegmode            = 'remfull';
    cfg.outputfile            = [dir_cleandata, resultprefix, '_cleanData_', sprintf('%ds', trialDuration)];
    
    hcp_extract_allfromrun(cfg);

end 
end