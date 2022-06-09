clear
close all

addpath(genpath(fullfile(pwd, 'helper_fun')));

%% Select
analysis_type = 'ROI';                                  % 'ROI' 'exploratory'
% 'ROI'         : ROI analysis per hemisphere. FDR correction for 33 regions.
% 'exploratory' : exploratory analysis per hemisphere. FDR correction for 180 regions.

duration = '2s';                                        % '2s' '5s' '10s'
conn_metric = 'powcorrOrtho';                           % 'powcorrOrtho' 'icoh' 'dwPLI'
band_name = 'collapseFreq';                             % 'delta' 'theta' 'alpha' 'beta' 'gamma' 'collapseFreq'
seed_targets = {'right-right' 'left-left'};             % 'right-right' 'right-left' 'left-left' 'left-right'

% contrast
contrast = true;

% figures
fsaverage_fig = true;
circular_graph = true;

font_size = 18; %18
flag.contours = true;

% save?
flag.save_fig = false;
flag.dir_fig = load_path('save_fig');

%% fixed
flag.script = 'collapseTargets';
flag.seed_targets = seed_targets;

outputs_first = helper_visualizeResults_functionConn_leftORright(analysis_type, contrast, duration, conn_metric, band_name, seed_targets{1});
outputs_second = helper_visualizeResults_functionConn_leftORright(analysis_type, contrast, duration, conn_metric, band_name, seed_targets{2});

alpha = 0.05;
corrected = true;
correction ='FDR';
flag.duration = duration;

%%
seed_target = 'right-right';
idx_IFJa = outputs_first.helper.idx_IFJa;
idx_IFJp = outputs_first.helper.idx_IFJp;
idx_FEF = outputs_first.helper.idx_FEF;
ROIs = outputs_first.helper.ROIs;
label = outputs_first.label;

if strcmp(band_name, 'collapseFreq')
    h_IFJa_all = logical(outputs_first.stat.h_IFJa_all) + logical(outputs_second.stat.h_IFJa_all);
    h_IFJp_all = logical(outputs_first.stat.h_IFJp_all) + logical(outputs_second.stat.h_IFJp_all);
    if contrast == false
        h_FEF_all = logical(outputs_first.stat.h_FEF_all) + logical(outputs_second.stat.h_FEF_all);
    end
else
    h_IFJa_all = logical(outputs_first.stat.h_adj_IFJa) + logical(outputs_second.stat.h_adj_IFJa);
    h_IFJp_all = logical(outputs_first.stat.h_adj_IFJp) + logical(outputs_second.stat.h_adj_IFJp);
    if contrast == false
        h_FEF_all = logical(outputs_first.stat.h_adj_FEF) + logical(outputs_second.stat.h_adj_FEF);
    end
end

tmp_cohspctrm = outputs_first.conn_matrix(181:360,181:360) + outputs_second.conn_matrix(1:180,1:180);
conn_matrix = zeros(360,360);
conn_matrix(181:360,181:360) = tmp_cohspctrm;
conn_matrix(idx_IFJa, ROIs.idx) = conn_matrix(idx_IFJa, ROIs.idx)./h_IFJa_all;
conn_matrix(idx_IFJp, ROIs.idx) = conn_matrix(idx_IFJp, ROIs.idx)./h_IFJp_all;
if contrast == false
    conn_matrix(idx_FEF, ROIs.idx) = conn_matrix(idx_FEF, ROIs.idx)./h_FEF_all;
end
conn_matrix(isinf(conn_matrix)) = 0;
conn_matrix(isnan(conn_matrix)) = 0;
conn_matrix(ROIs.idx, idx_IFJa) = conn_matrix(idx_IFJa, ROIs.idx);
conn_matrix(ROIs.idx, idx_IFJp) = conn_matrix(idx_IFJp, ROIs.idx);
if contrast == false
    conn_matrix(ROIs.idx, idx_FEF) = conn_matrix(idx_FEF, ROIs.idx);
end

conn_mat = struct;
conn_mat.cohspctrm = conn_matrix;
conn_mat.label = outputs_first.label;
conn_mat.brainordinate = outputs_first.helper.brainordinate;

if contrast == true
    colorlim_value = [-4 4];
else
    colorlim_value = [0 6];
end

%% fsaverage
if fsaverage_fig == true

    % to make black the areas outside ROI
    if and(contrast == true, strcmp(analysis_type, 'ROI'))
        cbar = cmap_rbw;
        n = size(cbar, 1);
        up = cbar(1:n/2,:);
        middle = [0.5 0.5 0.5];
        down = cbar(n/2+1:n,:);
        set(0, 'DefaultFigureColormap', [up;middle;down]);
        
        conn_mat.cohspctrm(ROIs.idx(~h_IFJa_all),idx_IFJa) = 0.04;
        conn_mat.cohspctrm(ROIs.idx(~h_IFJp_all),idx_IFJp) = 0.04;
        conn_mat.cohspctrm(idx_IFJa,ROIs.idx(~h_IFJa_all)) = 0.04;
        conn_mat.cohspctrm(idx_IFJp,ROIs.idx(~h_IFJp_all)) = 0.04;
        
    elseif and(contrast == false, strcmp(analysis_type, 'ROI'))
        cbar = hot;
        base1 = [1 1 1];
        base2 = [0.25 0.25 0.25];
        set(0, 'DefaultFigureColormap', [base1; base2; cbar]);
        
        conn_mat.cohspctrm(ROIs.idx(~h_IFJa_all),idx_IFJa) = -0.06;
        conn_mat.cohspctrm(ROIs.idx(~h_IFJp_all),idx_IFJp) = -0.06;
        conn_mat.cohspctrm(ROIs.idx(~h_FEF_all),idx_FEF) = -0.06;
        conn_mat.cohspctrm(idx_IFJa,ROIs.idx(~h_IFJa_all)) = -0.06;
        conn_mat.cohspctrm(idx_IFJp,ROIs.idx(~h_IFJp_all)) = -0.06;
        conn_mat.cohspctrm(idx_FEF,ROIs.idx(~h_FEF_all)) = -0.06;
        
        minv = -0.06;
        maxv = 6;
        colorlim_value = [minv maxv];
    elseif contrast == false
        set(0, 'DefaultFigureColormap', hot);
    end
    
    % in order to run the function
    conn_mat.brainordinate.pos = conn_mat.brainordinate.pos*1000;
    
    %% save the figure
    if flag.save_fig == true
        
        cimec = struct;
        cimec.conn_metric = conn_metric;
        cimec.band_name = band_name;
        cimec.seed_target = seed_target;
        cimec.statistics.analysis_type = analysis_type;
        cimec.statistics.alpha = alpha;
        cimec.statistics.corrected = corrected;
        cimec.statistics.correction = correction;
        cimec.flag = flag;
        
        if seed_target(1) == 'r'
            seed = 'IFJa';
            cimec.seed = seed;
            pos2d_IFJa_R = [52.5323  -69.3638   71.0650; 52.5323   69.4183   71.0650];
            tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_R);
            helper_save_figure(cimec);
            close all;
            
            seed = 'IFJp';
            cimec.seed = seed;
            pos2d_IFJp_R = [45.6139  -69.3638   76.8706 ; 45.6139   69.4183   76.8706];
            tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_R);
            helper_save_figure(cimec);
            close all;
            
            if contrast == false
                seed = 'FEF';
                cimec.seed = seed;
                pos2d_FEF_R = [25.8639  -69.3638   96.3111; 25.8639   69.4183   96.3111];
                tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_R);
                helper_save_figure(cimec);
                close all;
            end
            
        elseif seed_target(1) == 'l'
            seed = 'IFJa';
            cimec.seed = seed;
            pos2d_IFJa_L = [50.5105   69.4183   69.5634 ; 53.1465  -69.3638   68.2860];
            tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_L);
            helper_save_figure(cimec);
            close all;
            
            seed = 'IFJp';
            cimec.seed = seed;
            pos2d_IFJp_L = [43.6193   69.4183   80.0619 ; 43.6193  -69.3638   80.0619];
            tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_L);
            helper_save_figure(cimec);
            close all;
            
            if contrast == false
                seed = 'FEF';
                cimec.seed = seed;
                pos2d_FEF_L = [30.9891   69.4183  105.2148 ; 34.2627  -69.3638   98.8225];
                tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_L);
                helper_save_figure(cimec);
                close all;
            end
        end
    else
        tutorial_nwa_connectivityviewer(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours);
    end
end

%% circularGraph
if circular_graph == true
    results = struct;
    results.stat.analysis_type = analysis_type;
    results.helper.band_name = band_name;
    if strcmp(band_name, 'collapseFreq')
        results.stat.h_IFJa_all = outputs_first.stat.h_IFJa_all + outputs_second.stat.h_IFJa_all;
        results.stat.h_IFJp_all = outputs_first.stat.h_IFJp_all + outputs_second.stat.h_IFJp_all;
        if contrast == false
            results.stat.h_FEF_all = outputs_first.stat.h_FEF_all + outputs_second.stat.h_FEF_all;
        end
    else
        results.stat.h_adj_IFJa = outputs_first.stat.h_adj_IFJa + outputs_second.stat.h_adj_IFJa;
        results.stat.h_adj_IFJp = outputs_first.stat.h_adj_IFJp + outputs_second.stat.h_adj_IFJp;
        if contrast == false
            results.stat.h_adj_FEF = outputs_first.stat.h_adj_FEF + outputs_second.stat.h_adj_FEF;
        end
    end
    
    tmp_data = conn_matrix;
    tmp_data(tmp_data == 1) = NaN;
    tmp_data(idx_IFJa, idx_IFJa) = 0;
    tmp_data(idx_IFJp, idx_IFJp) = 0;
    tmp_data(idx_FEF, idx_FEF) = 0;
    results.conn_matrix = tmp_data;
    
    results.helper.idx_FEF = idx_FEF;
    results.helper.idx_IFJa = idx_IFJa;
    results.helper.idx_IFJp = idx_IFJp;
    results.helper.ROIs = ROIs;
    results.label = label;
    
    %%
    flag.collapseTargets = true;
    flag.contrast = contrast;
    results.cimec = outputs_first.cimec;
    results.cimec.band_name = band_name;
    results.cimec.seed_target = char(string(join(seed_targets)));
    fontSize = font_size;
        
    helper_functionConn_circularGraph(results, fontSize, flag);
end

%% helper function
function outputs = helper_visualizeResults_functionConn_leftORright(analysis_type, contrast, duration, conn_metric, band_name, seed_target)

addpath(genpath(fullfile(pwd, 'helper_fun')));
% colormapeditor
set(0, 'DefaultFigureColormap', cmap_rbw);

statistics.analysis_type = analysis_type;                        % 'ROI' 'exploratory'
% 'ROI'         : ROI analysis per hemisphere. FDR correction for 33 regions.
% 'exploratory' : exploratory analysis per hemisphere. FDR correction for 180 regions.

% figures
fsaverage_fig = false;
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
save_fig = false;
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
