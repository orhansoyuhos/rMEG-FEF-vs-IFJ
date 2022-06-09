band_name = ["delta" "theta" "alpha" "beta" "gamma"];
seed_target = ["left-left" "right-right"];

for ii = 1:5
    for kk = 1:2
        helper_automaticSave_functional(char(band_name(ii)), char(seed_target(kk)))
    end
end


function helper_automaticSave_functional(band_name, seed_target)

addpath(genpath(fullfile(pwd, 'helper_fun')));
% colormapeditor
set(0, 'DefaultFigureColormap', cmap_rbw);

%% Select
statistics.analysis_type = 'ROI';                        % 'ROI' 'exploratory'
% 'ROI'         : ROI analysis per hemisphere. FDR correction for 33 regions. 
% 'exploratory' : exploratory analysis per hemisphere. FDR correction for 180 regions.

duration = '2s';                                        % '2s' '5s' '10s' 
conn_metric = 'powcorrOrtho';                           % 'powcorrOrtho' 'icoh' 'dwPLI'
% band_name = 'delta';                                  % 'delta' 'theta' 'alpha' 'beta' 'gamma' 'collapseFreq'
% seed_target = 'left-left';                            % 'right-right' 'right-left' 'left-left' 'left-right'

% contrast
contrast = true;

% figures
fsaverage_fig = true;
circular_graph = false;
font_size = 15; 

% Wilcoxon signed rank test
statistics.alpha = 0.05;                               % 0.05 0.01 0.001
% FDR correction
statistics.corrected = true;                           % true false
statistics.correction = 'FDR';

% to show the contours of parcels; 
%%keep it false for a faster processing
flag.contours = false; 

% save the figures automatic?
save_fig = true;
dir_fig = load_path('save_fig');

%% fixed
flag.fsaverage = fsaverage_fig; % to show the figure. 
flag.circular_graph = circular_graph;
flag.font_size =  font_size;
flag.save_fig = save_fig;
flag.dir_fig = dir_fig;
flag.duration = duration;
inputfile_conn = fullfile(load_path(['Outputs_' duration '_55subjs_' conn_metric]), ['55Subjects_groupResults_' duration '.mat']);

%% Run the function
cimec.inputfile_conn = inputfile_conn;
cimec.band_name = band_name;
cimec.seed_target = seed_target;
cimec.conn_metric = conn_metric;
cimec.statistics = statistics;
cimec.contrast = contrast;
cimec.flag = flag;

if strcmp(band_name, 'collapseFreq')
    outputs = helper_functionConn_fsaverage_freqCollapsed(cimec);    
else
    outputs = helper_functionConn_fsaverage(cimec);
end

end