function results = helper_functionConn_fsaverage_freqCollapsed(cimec)

band_name = cimec.band_name;
flag = cimec.flag;
conn_metric = cimec.conn_metric;
contrast = cimec.contrast;
analysis_type = cimec.statistics.analysis_type;
seed_target = cimec.seed_target;
alpha_value = cimec.statistics.alpha;
duration = cimec.flag.duration;
correction = cimec.statistics.correction;
corrected = cimec.statistics.corrected;

%% All freq bands
cimec.flag.circular_graph = false;
cimec.flag.fsaverage = false;
cimec.band_name = 'delta';
delta = helper_functionConn_fsaverage(cimec);
delta.conn_matrix(delta.conn_matrix == 1) = NaN;
cimec.band_name = 'theta';
theta = helper_functionConn_fsaverage(cimec);
theta.conn_matrix(theta.conn_matrix == 1) = NaN;
cimec.band_name = 'alpha';
alpha = helper_functionConn_fsaverage(cimec);
alpha.conn_matrix(alpha.conn_matrix == 1) = NaN;
cimec.band_name = 'beta';
beta = helper_functionConn_fsaverage(cimec);
beta.conn_matrix(beta.conn_matrix == 1) = NaN;
cimec.band_name = 'gamma';
gamma = helper_functionConn_fsaverage(cimec);
gamma.conn_matrix(gamma.conn_matrix == 1) = NaN;

idx_IFJa = delta.helper.idx_IFJa;
idx_IFJp = delta.helper.idx_IFJp;
idx_FEF = delta.helper.idx_FEF;
ROIs.idx = delta.helper.ROIs.idx;

h_IFJa_all = delta.stat.h_adj_IFJa + theta.stat.h_adj_IFJa + alpha.stat.h_adj_IFJa + beta.stat.h_adj_IFJa + gamma.stat.h_adj_IFJa;
h_IFJp_all = delta.stat.h_adj_IFJp + theta.stat.h_adj_IFJp + alpha.stat.h_adj_IFJp + beta.stat.h_adj_IFJp + gamma.stat.h_adj_IFJp;
if contrast == false
    h_FEF_all = delta.stat.h_adj_FEF + theta.stat.h_adj_FEF + alpha.stat.h_adj_FEF + beta.stat.h_adj_FEF + gamma.stat.h_adj_FEF;
end


conn_matrix = zeros(360,360);
tmp_mat = (delta.conn_matrix + theta.conn_matrix + alpha.conn_matrix + beta.conn_matrix + gamma.conn_matrix);
tmp_mat(idx_IFJa, ROIs.idx) = tmp_mat(idx_IFJa, ROIs.idx)./h_IFJa_all;
tmp_mat(idx_IFJp, ROIs.idx) = tmp_mat(idx_IFJp, ROIs.idx)./h_IFJp_all;
if contrast == false
    tmp_mat(idx_FEF, ROIs.idx) = tmp_mat(idx_FEF, ROIs.idx)./h_FEF_all;
end
tmp_mat(isinf(tmp_mat)) = 0;
tmp_mat(isnan(tmp_mat)) = 0;

conn_matrix(:,idx_IFJa) = tmp_mat(idx_IFJa,:);
conn_matrix(:,idx_IFJp) = tmp_mat(idx_IFJp,:);
conn_matrix(idx_IFJa,:) = tmp_mat(idx_IFJa,:);
conn_matrix(idx_IFJp,:) = tmp_mat(idx_IFJp,:);
if contrast == false
    conn_matrix(:,idx_FEF) = tmp_mat(idx_FEF,:);
    conn_matrix(idx_FEF,:) = tmp_mat(idx_FEF,:);
end
conn_matrix(idx_IFJa, idx_IFJa) = 0;
conn_matrix(idx_IFJp, idx_IFJp) = 0;
conn_matrix(idx_FEF, idx_FEF) = 0;


%% Visualize using the helper function
if flag.fsaverage == true
    
    allFreq.cohspctrm = conn_matrix;
    allFreq.brainordinate = delta.helper.brainordinate;
    % in order to run the function
    allFreq.brainordinate.pos = allFreq.brainordinate.pos*1000;
    if contrast == true
        colorlim_value = [-4 4];
    else
        colorlim_value = [0 6];
    end
    
    %% To make black the areas outside ROI
    if and(contrast == true, strcmp(analysis_type, 'ROI'))
        cbar = cmap_rbw;
        n = size(cbar, 1);
        up = cbar(1:n/2,:);
        middle = [0.5 0.5 0.5];
        down = cbar(n/2+1:n,:);
        set(0, 'DefaultFigureColormap', [up;middle;down]);
        
        allFreq.cohspctrm(ROIs.idx(~h_IFJa_all),idx_IFJa) = 0.04;
        allFreq.cohspctrm(ROIs.idx(~h_IFJp_all),idx_IFJp) = 0.04;
        allFreq.cohspctrm(idx_IFJa,ROIs.idx(~h_IFJa_all)) = 0.04;
        allFreq.cohspctrm(idx_IFJp,ROIs.idx(~h_IFJp_all)) = 0.04;
        
    elseif and(contrast == false, strcmp(analysis_type, 'ROI'))
        cbar = hot;
        base1 = [1 1 1];
        base2 = [0.25 0.25 0.25];
        set(0, 'DefaultFigureColormap', [base1; base2; cbar]);
        
        allFreq.cohspctrm(ROIs.idx(~h_IFJa_all),idx_IFJa) = -0.06;
        allFreq.cohspctrm(ROIs.idx(~h_IFJp_all),idx_IFJp) = -0.06;
        allFreq.cohspctrm(ROIs.idx(~h_FEF_all),idx_FEF) = -0.06;
        allFreq.cohspctrm(idx_IFJa,ROIs.idx(~h_IFJa_all)) = -0.06;
        allFreq.cohspctrm(idx_IFJp,ROIs.idx(~h_IFJp_all)) = -0.06;
        allFreq.cohspctrm(idx_FEF,ROIs.idx(~h_FEF_all)) = -0.06;
        
        minv = -0.06;
        maxv = 6;
        colorlim_value = [minv maxv];
    elseif contrast == false
        set(0, 'DefaultFigureColormap', hot);
    end
    
    %% save the figure
    if flag.save_fig == true
        cimec = struct;
        cimec.conn_metric = conn_metric;
        cimec.band_name = band_name;
        cimec.seed_target = seed_target;
        cimec.statistics.analysis_type = analysis_type;
        cimec.statistics.alpha = alpha_value;
        cimec.statistics.corrected = corrected;
        cimec.statistics.correction = correction;
        cimec.flag.duration = duration;
        cimec.flag.dir_fig = flag.dir_fig;
        
        if seed_target(1) == 'r'
            seed = 'IFJa';
            cimec.seed = seed;
            pos2d_IFJa_R = [52.5323  -69.3638   71.0650; 52.5323   69.4183   71.0650];
            tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_R);
            helper_save_figure(cimec);
            close all;
            
            seed = 'IFJp';
            cimec.seed = seed;
            pos2d_IFJp_R = [45.6139  -69.3638   76.8706 ; 45.6139   69.4183   76.8706];
            tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_R);
            helper_save_figure(cimec);
            close all;
            
            if contrast == false
                seed = 'FEF';
                cimec.seed = seed;
                pos2d_FEF_R = [25.8639  -69.3638   96.3111; 25.8639   69.4183   96.3111];
                tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_R);
                helper_save_figure(cimec);
                close all;
            end
            
        elseif seed_target(1) == 'l'
            seed = 'IFJa';
            cimec.seed = seed;
            pos2d_IFJa_L = [50.5105   69.4183   69.5634 ; 53.1465  -69.3638   68.2860];
            tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_L);
            helper_save_figure(cimec);
            close all;
            
            seed = 'IFJp';
            cimec.seed = seed;
            pos2d_IFJp_L = [43.6193   69.4183   80.0619 ; 43.6193  -69.3638   80.0619];
            tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_L);
            helper_save_figure(cimec);
            close all;
            
            if contrast == false
                seed = 'FEF';
                cimec.seed = seed;
                pos2d_FEF_L = [30.9891   69.4183  105.2148 ; 34.2627  -69.3638   98.8225];
                tutorial_nwa_connectivityviewer_save_figures(allFreq, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_L);
                helper_save_figure(cimec);
                close all;
            end
        end
    else
        tutorial_nwa_connectivityviewer(allFreq, 'cohspctrm', colorlim_value, cimec.seed_target, cimec.flag.contours);
    end
end

results.conn_matrix = conn_matrix;
results.contrast = contrast;
results.label = delta.label;
results.helper.brainordinate = delta.helper.brainordinate;
results.helper.idx_FEF = delta.helper.idx_FEF;
results.helper.idx_IFJa = delta.helper.idx_IFJa;
results.helper.idx_IFJp = delta.helper.idx_IFJp;
results.helper.ROIs = delta.helper.ROIs;
results.stat.analysis_type = analysis_type;
results.helper.band_name = band_name;
results.stat.h_IFJa_all = h_IFJa_all;
results.stat.h_IFJp_all = h_IFJp_all;
if exist('h_FEF_all','var'); results.stat.h_FEF_all = h_FEF_all; end

results.perSubject.IFJa = (delta.perSubject.IFJa + theta.perSubject.IFJa + alpha.perSubject.IFJa + beta.perSubject.IFJa + gamma.perSubject.IFJa)/5;
results.perSubject.IFJp = (delta.perSubject.IFJp + theta.perSubject.IFJp + alpha.perSubject.IFJp + beta.perSubject.IFJp + gamma.perSubject.IFJp)/5;
results.perSubject.FEF = (delta.perSubject.FEF + theta.perSubject.FEF + alpha.perSubject.FEF + beta.perSubject.FEF + gamma.perSubject.FEF)/5;

cimec.band_name = band_name;
results.cimec = cimec;

%% circularGraph
if flag.circular_graph == true
    fontSize = flag.font_size;
    flag.collapseTargets = false;
    flag.contrast = contrast; 
        
    helper_functionConn_circularGraph(results, fontSize, flag);
end

end
