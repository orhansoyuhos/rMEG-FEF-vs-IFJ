close all
clear

addpath(genpath(fullfile(pwd, 'helper_fun')));

%% MEANING of colors
% y: target region
%
% GREEN  : seed --> y direction is significant. A seed region explains the varience in a target better.
%        ~ sending signal
% MAGENTA : y--> seed direction is significant. A target region explains the varience in a seed better. 
%        ~ receiving signal

%% Select            
statistics.analysis_type = 'ROI';              % 'ROI' 'exploratory'
% 'ROI'         : ROI analysis per hemisphere. FDR correction for 33+3 seed regions. 
% 'exploratory' : exploratory analysis per hemisphere. FDR correction for 180 regions.

duration = '2s'; %fixed                                 % only '2s'
conn_metric = 'pdc';                                    % 'pdc' 
band_name = 'beta';                                     % 'delta' 'theta' 'alpha' 'beta' 'gamma'
seed_target = 'right-right';                            % 'right-right' 'right-left' 'left-left' 'left-right'

%%
% figures
fsaverage_fig = true;
violinplot = false;

% Paired Wilcoxon paired test against zero
statistics.alpha = 0.001;                               % 0.05 0.01 0.001
% FDR correction
statistics.corrected = true;                           % true false
statistics.correction = 'FDR';

% to show the contours of parcels; 
%%keep it false for a faster processing
flag.contours = false; 

% save the figures automatic?
save_fig = false;
dir_fig = load_path('save_fig');

%%
colorlim_value = [-4 4];
%%% fixed
circular_graph = false;
font_size = 15; 
flag.fsaverage = fsaverage_fig; % to show the figure. 
flag.circular_graph = circular_graph;
flag.font_size =  font_size;
flag.save_fig = save_fig;
flag.duration = duration;
flag.dir_fig = dir_fig;
inputfile_conn = fullfile(load_path(['Outputs_' duration '_55subjs_' conn_metric]), ['55Subjects_groupResults_' duration '.mat']);

%% Run the function
cimec.inputfile_conn = inputfile_conn;
cimec.band_name = band_name;
cimec.seed_target = seed_target;
cimec.conn_metric = conn_metric;
cimec.statistics = statistics;
cimec.colorlim_value = colorlim_value;
cimec.flag = flag;

outputs = helper_effectiveConn_fsaverage(cimec);

%% violinplot
if violinplot == true
    helper_violinplot(outputs, conn_metric, band_name, seed_target,  duration)
end