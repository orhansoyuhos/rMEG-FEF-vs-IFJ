function results = helper_effectiveConn_fsaverage(cimec)

inputfile_conn = cimec.inputfile_conn;
band_name = cimec.band_name;
seed_target = cimec.seed_target;
conn_metric = cimec.conn_metric;
statistics = cimec.statistics;
colorlim_value = cimec.colorlim_value;
flag = cimec.flag;

%%
addpath(genpath(fullfile(fileparts(fileparts(pwd)), 'helper_fun')));

statistics.contrast = true;
statistics.stat = true;
statistics.f_stats = @signrank;
statistics.f_correction = @fdr_bh;

IFJ_merged = false;
ROIs.do = true; % fixed

colorlim = true;

%%
% fsaverage
inflated = true;
if strcmp(statistics.analysis_type, 'exploratory')
    smoothing = '70'; % percent
else
    smoothing = '100'; % percent
end
inputfile_fs = fullfile(load_path('inputfile_fs'), 'fsaverage.mat');
inputfile_fs_inflated = fullfile(load_path('inputfile_fs'), sprintf('fsaverage_inflated_%s.mat', smoothing));

% seed to target
dorsal_R = {  'R_V6_ROI R', 'R_V6A_ROI R', 'R_V7_ROI R', 'R_IPS1_ROI R', 'R_IP1_ROI R', 'R_MIP_ROI R',  ...
    'R_VIP_ROI R', 'R_LIPd_ROI R', 'R_LIPv_ROI R', 'R_7AL_ROI R', 'R_7PC_ROI R', 'R_7PL_ROI R',    ...
    'R_V3A_ROI R', 'R_V3B_ROI R', 'R_V3CD_ROI R', 'R_LO3_ROI R', 'R_MT_ROI R'};
ventral_R = {'R_V8_ROI R', 'R_VMV1_ROI R', 'R_VMV2_ROI R', 'R_VMV3_ROI R', 'R_PHT_ROI R', 'R_PH_ROI R', ...
    'R_TE1a_ROI R', 'R_TE1m_ROI R', 'R_TE1p_ROI R', 'R_TE2a_ROI R', 'R_TE2p_ROI R',  ...
    'R_TF_ROI R', 'R_TGv_ROI R', 'R_FFC_ROI R', 'R_PIT_ROI R', 'R_VVC_ROI R'};

dorsal_L = {  'L_V6_ROI L', 'L_V6A_ROI L', 'L_V7_ROI L', 'L_IPS1_ROI L', 'L_IP1_ROI L', 'L_MIP_ROI L',  ...
    'L_VIP_ROI L', 'L_LIPd_ROI L', 'L_LIPv_ROI L', 'L_7AL_ROI L', 'L_7PC_ROI L', 'L_7PL_ROI L',    ...
    'L_V3A_ROI L', 'L_V3B_ROI L', 'L_V3CD_ROI L', 'L_LO3_ROI L', 'L_MT_ROI L'};
ventral_L = {'L_V8_ROI L', 'L_VMV1_ROI L', 'L_VMV2_ROI L', 'L_VMV3_ROI L', 'L_PHT_ROI L', 'L_PH_ROI L', ...
    'L_TE1a_ROI L', 'L_TE1m_ROI L', 'L_TE1p_ROI L', 'L_TE2a_ROI L', 'L_TE2p_ROI L',  ...
    'L_TF_ROI L', 'L_TGv_ROI L', 'L_FFC_ROI L', 'L_PIT_ROI L', 'L_VVC_ROI L'};


if strcmp(seed_target, 'right-right')
    between = {'R_IFJa_ROI R', 'R_IFJp_ROI R', 'R_FEF_ROI R'}; % 1-3 & 2-3
    ROIs.name = ['R_FEF_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R', dorsal_R, ventral_R];
elseif strcmp(seed_target, 'right-left')
    between = {'R_IFJa_ROI R', 'R_IFJp_ROI R', 'R_FEF_ROI R'}; % 1-3 & 2-3
    ROIs.name = ['R_FEF_ROI R', 'R_IFJa_ROI R', 'R_IFJp_ROI R', dorsal_L ventral_L];
elseif strcmp(seed_target, 'left-left')
    between = {'L_IFJa_ROI L', 'L_IFJp_ROI L', 'L_FEF_ROI L'}; % 1-3 & 2-3
    ROIs.name = ['L_FEF_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', dorsal_L ventral_L];
elseif strcmp(seed_target, 'left-right')
    between = {'L_IFJa_ROI L', 'L_IFJp_ROI L', 'L_FEF_ROI L'}; % 1-3 & 2-3
    ROIs.name = ['L_FEF_ROI L', 'L_IFJa_ROI L', 'L_IFJp_ROI L', dorsal_R ventral_R];
end


%% Inputs
cimec = struct;
cimec.between = between;
cimec.statistics = statistics;
cimec.band_name = band_name;
cimec.IFJ_merged = IFJ_merged;
cimec.ROIs = ROIs;
cimec.colorlim_value = colorlim_value;
cimec.colorlim = colorlim;
cimec.inputfile_conn = inputfile_conn;
cimec.inputfile_fs = inputfile_fs;
cimec.inflated = inflated;
cimec.inputfile_fs_inflated = inputfile_fs_inflated;
cimec.flag = flag;
cimec.conn_metric = conn_metric;
cimec.seed_target = seed_target;

results = visualize_connectivity(cimec);

%% Function
    function results = visualize_connectivity(cimec)
        
        %% Info about the data
        inputfile_conn = cimec.inputfile_conn; % e.g. '...\HCP_Subjects\Outputs\connectivity_results\175237_connectivityResults__2s';
        inputfile_fs = cimec.inputfile_fs; % e.g. '...\HCP_Subjects\Outputs\fsaverage_anatomy\175237_fsaverage';
        
        contrast = cimec.statistics.contrast; % contrast map or not
        between = cimec.between; % {'R_IFJa_ROI R', 'R_IFJp_ROI R', 'R_FEF_ROI R'}; % 1-2 & 1-3
        
        stat = cimec.statistics.stat;
        f_stats = cimec.statistics.f_stats; % e.g. @signrank
        alpha = cimec.statistics.alpha; % e.g. 0.05
        f_correction = cimec.statistics.f_correction; % correction method e.g. @fdr_bh
        corrected = cimec.statistics.corrected; % e.g. true
        
        IFJ_merged = cimec.IFJ_merged;
        
        ROIs = cimec.ROIs;
        
        band_name = cimec.band_name; % e.g. 'alpha';
        colorlim_value = cimec.colorlim_value; % e.g. [0 0.10];
        colorlim = cimec.colorlim; % e.g. true;
        
        inflated = cimec.inflated; % e.g. true
        inputfile_fs_inflated = cimec.inputfile_fs_inflated;
        
        flag = cimec.flag;
        conn_metric = cimec.conn_metric;
        seed_target = cimec.seed_target;
        
        %% Frequency bands
        band_names  =   {  'delta' ,   'theta' ,   'alpha' ,   'beta'  ,   'gamma'  };
        freq_bands  =   [    1 4   ;     4 8   ;    8 13   ;    13 30  ;    30 100  ];
        
        %% Load data
        % avgFreq_icoh
        load(inputfile_conn, 'group_results');
        conn_results = group_results;
        
        % fsaverage
        load(inputfile_fs, 'fsaverage')
        
        %% Prepare the anatomy
        % Parcellation info (Brainstorm -> Fieldtrip)
        tess_cortex_pial = fsaverage;
        
        if inflated == true
            fs_inflated = load(inputfile_fs_inflated);
            tess_cortex_pial.Vertices = fs_inflated.Vertices;
            tess_cortex_pial.Faces = fs_inflated.Faces;
        end
        
        Scouts = tess_cortex_pial.Atlas(end-1).Scouts;
        Scouts_Vertices = {Scouts.Vertices}';
        Scouts_Label = {Scouts.Label}';
        
        parcellation = zeros(length([tess_cortex_pial.Vertices]), 1);
        for idx = 1:length(Scouts_Vertices)
            per_parcel = Scouts_Vertices{idx};
            
            for ii = 1:length(per_parcel)
                parcellation(per_parcel(ii)) = idx;
            end
        end
        
        brainordinate.parcellation = parcellation;
        brainordinate.parcellationlabel = Scouts_Label;
        brainordinate.pos = tess_cortex_pial.Vertices;
        brainordinate.tri = tess_cortex_pial.Faces;
        brainordinate.brainstructure = tess_cortex_pial.SulciMap+1;
        brainordinate.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'}';
        
        %% Prepare the connectivity matrix
        conn_mat = struct;
        conn_mat.label = Scouts_Label;
        conn_mat.band_names = band_name;
        conn_mat.brainordinate = brainordinate;
        
        idx_band = find(strcmp(band_names, band_name));
        tmp_data = group_results(idx_band).data_perSubject;
        tmp_data(isnan(tmp_data)) = 1;
        
        %% seeds
        idx_IFJa = find(strcmp(conn_mat.label, between{1}));
        idx_IFJp = find(strcmp(conn_mat.label, between{2}));
        idx_FEF = find(strcmp(conn_mat.label, between{3}));
        
        IFJa_12(:, :) = tmp_data(idx_IFJa,:,:);
        IFJp_12(:, :) = tmp_data(idx_IFJp,:,:);
        FEF_12(:, :) = tmp_data(idx_FEF,:,:);
        IFJa_21(:, :) = tmp_data(:,idx_IFJa,:);
        IFJp_21(:, :) = tmp_data(:,idx_IFJp,:);
        FEF_21(:, :) = tmp_data(:,idx_FEF,:);
        
        
        % subtract two directions
        % positive values mean seed region predicts the target
        % negative values mean target is predicted by the seed region
        IFJa(:, :) = IFJa_12 - IFJa_21;
        IFJp(:, :) = IFJp_12 - IFJp_21;
        FEF(:, :) = FEF_12 - FEF_21;
        
        %% ROIs/exploratory
        if strcmp(statistics.analysis_type, 'exploratory')
            all_labels = conn_mat.label;
            if strcmp(between{1}(1), 'R')
                ROIs.name = all_labels(181:360);
            elseif strcmp(between{1}(1), 'L')
                ROIs.name = all_labels(1:180);
            end
        end
        
        nROIs = length(ROIs.name);
        for ii = 1:nROIs
            ROIs.idx(ii) = find(strcmp(conn_mat.label, ROIs.name{ii}));
        end
        
        
        mask = zeros(360,1);
        mask(ROIs.idx) = 1;
        
        %% normality test
        % stat_IFJa = IFJa;
        % stat_IFJa(stat_IFJa(:,1)==1) = NaN;
        % stat_FEF = FEF;
        % stat_FEF(stat_FEF(:,1)==1) = NaN;
        %
        %
        % for ii = 1:nROIs
        %     figure
        % %     histogram(atanh(stat_IFJa(ROIs.idx(ii),:)-stat_FEF(ROIs.idx(ii),:)));
        % %     histogram(stat_IFJa(ROIs.idx(ii),:)-stat_FEF(ROIs.idx(ii),:));
        % %     histogram(stat_IFJa(ROIs.idx(ii),:));
        %     histogram(atanh(stat_IFJa(ROIs.idx(ii),:)));
        % end
        %
        % figure
        % hist(atanh(stat_IFJa(ROIs.idx(33),:)-stat_FEF(ROIs.idx(33),:)));
        % figure
        % hist(stat_IFJa(ROIs.idx(33),:)-stat_FEF(ROIs.idx(33),:));
        
        
        %%
        if contrast == true
            %% Inferential statistics
            %% Wilcoxon signed rank test
            if IFJ_merged == false
                % R_IFJa
                for ii = 1:nROIs
                    pair_1 = IFJa(ROIs.idx(ii),:)';
                    pair_2 = 0;
                    
                    [P, H, ~] = f_stats(pair_1, pair_2, 'Alpha', alpha, 'Tail', 'both');
                    
                    H_Wilcoxon_IFJa(ii) = H;
                    P_Wilcoxon_IFJa(ii) = P;
                end
                
                % R_IFJp
                for ii = 1:nROIs
                    pair_1 = IFJp(ROIs.idx(ii),:)';
                    pair_2 = 0;
                    
                    [P, H, ~] = f_stats(pair_1, pair_2, 'Alpha', alpha, 'Tail', 'both');
                    
                    H_Wilcoxon_IFJp(ii) = H;
                    P_Wilcoxon_IFJp(ii) = P;
                end
                
                % FEF
                for ii = 1:nROIs
                    pair_1 = FEF(ROIs.idx(ii),:)';
                    pair_2 = 0;
                    
                    [P, H, ~] = f_stats(pair_1, pair_2, 'Alpha', alpha, 'Tail', 'both');
                    
                    H_Wilcoxon_FEF(ii) = H;
                    P_Wilcoxon_FEF(ii) = P;
                end
            end
            %% FDR correction
            if corrected == false
                % clear NaNs
                H_Wilcoxon_IFJa(isnan(H_Wilcoxon_IFJa)) = 0;
                H_Wilcoxon_IFJp(isnan(H_Wilcoxon_IFJp)) = 0;
                
                % mask connectivity values
                group_IFJa_contrast = mean(IFJa, 2);
                right_group_IFJa = group_IFJa_contrast(ROIs.idx);
                right_group_IFJa(~H_Wilcoxon_IFJa) = 0;
                group_IFJa = zeros(360,1);
                group_IFJa(ROIs.idx) = right_group_IFJa;
                stats_IFJa = group_IFJa;
                
                group_IFJp_contrast = mean(IFJp, 2);
                right_group_IFJp = group_IFJp_contrast(ROIs.idx);
                right_group_IFJp(~H_Wilcoxon_IFJp) = 0;
                group_IFJp = zeros(360,1);
                group_IFJp(ROIs.idx) = right_group_IFJp;
                stats_IFJp = group_IFJp;
                
                % clear NaNs
                H_Wilcoxon_FEF(isnan(H_Wilcoxon_FEF)) = 0;
                
                % mask connectivity values
                group_FEF_contrast = mean(FEF, 2);
                right_group_FEF = group_FEF_contrast(ROIs.idx);
                right_group_FEF(~H_Wilcoxon_FEF) = 0;
                group_FEF = zeros(360,1);
                group_FEF(ROIs.idx) = right_group_FEF;
                stats_FEF = group_FEF;
                
            elseif corrected == true
                % FDR correction
                P_Wilcoxon_IFJa(isnan(P_Wilcoxon_IFJa)) = 1;
                P_Wilcoxon_IFJp(isnan(P_Wilcoxon_IFJp)) = 1;
                P_Wilcoxon_FEF(isnan(P_Wilcoxon_FEF)) = 1;
                [h_adj_IFJa, ~, ~, p_adj_IFJa] = f_correction(P_Wilcoxon_IFJa(~isnan(P_Wilcoxon_IFJa)), alpha, 'pdep', 'yes');
                [h_adj_IFJp, ~, ~, p_adj_IFJp] = f_correction(P_Wilcoxon_IFJp(~isnan(P_Wilcoxon_IFJp)), alpha, 'pdep', 'yes');
                [h_adj_FEF, ~, ~, p_adj_FEF]   = f_correction(P_Wilcoxon_FEF(~isnan(P_Wilcoxon_FEF)), alpha, 'pdep', 'yes');
                
                
                % z-score
                % for two-tailed: abs(norminv(p/2));
                z_adj_IFJa = abs(norminv(p_adj_IFJa/2));
                z_adj_IFJp = abs(norminv(p_adj_IFJp/2));
                z_adj_FEF = abs(norminv(p_adj_FEF/2));

                
                group_IFJa_contrast = mean(IFJa, 2);
                right_group_IFJa = group_IFJa_contrast(ROIs.idx);
                right_group_IFJa(~h_adj_IFJa) = 0;
                group_IFJa = zeros(360,1);
                group_IFJa(ROIs.idx) = z_adj_IFJa'.*sign(right_group_IFJa);
                stats_IFJa = group_IFJa;
                
                group_IFJp_contrast = mean(IFJp, 2);
                right_group_IFJp = group_IFJp_contrast(ROIs.idx);
                right_group_IFJp(~h_adj_IFJp) = 0;
                group_IFJp = zeros(360,1);
                group_IFJp(ROIs.idx) = z_adj_IFJp'.*sign(right_group_IFJp);
                stats_IFJp = group_IFJp;
                
                group_FEF_contrast = mean(FEF, 2);
                right_group_FEF = group_FEF_contrast(ROIs.idx);
                right_group_FEF(~h_adj_FEF) = 0;
                group_FEF = zeros(360,1);
                group_FEF(ROIs.idx) = z_adj_FEF'.*sign(right_group_FEF);
                stats_FEF = group_FEF;

                %                 % mask connectivity values
                %                 group_IFJa_contrast = mean(IFJa, 2);
                %                 right_group_IFJa = group_IFJa_contrast(ROIs.idx);
                %                 right_group_IFJa(~h_adj_IFJa) = 0;
                %                 group_IFJa = zeros(360,1);
                %                 group_IFJa(ROIs.idx) = right_group_IFJa;
                %                 stats_IFJa = group_IFJa;
                %
                %                 group_IFJp_contrast = mean(IFJp, 2);
                %                 right_group_IFJp = group_IFJp_contrast(ROIs.idx);
                %                 right_group_IFJp(~h_adj_IFJp) = 0;
                %                 group_IFJp = zeros(360,1);
                %                 group_IFJp(ROIs.idx) = right_group_IFJp;
                %                 stats_IFJp = group_IFJp;
                
%                 % mask connectivity values
%                 group_FEF_contrast = mean(FEF, 2);
%                 right_group_FEF = group_FEF_contrast(ROIs.idx);
%                 right_group_FEF(~h_adj_FEF) = 0;
%                 group_FEF = zeros(360,1);
%                 group_FEF(ROIs.idx) = right_group_FEF;
%                 stats_FEF = group_FEF;
            end
            
            %% Final matrix (360x360)
            if IFJ_merged == false
                tmp_data = diag(ones(360,1));
                
                tmp_data(:,idx_IFJa) = stats_IFJa';
                tmp_data(:,idx_IFJp) = stats_IFJp';
                tmp_data(:,idx_FEF) = stats_FEF';
            end
            
        end
        
        %% To make black the areas outside ROI
        cbar = customMap_effectiveConn;
        if strcmp(cimec.statistics.analysis_type, 'ROI')
            n = size(cbar, 1);
            up = cbar(1:n/2,:);
            middle = [0.5 0.5 0.5];
            down = cbar(n/2+1:n,:);
            set(0, 'DefaultFigureColormap', [up;middle;down]);
        else
            set(0, 'DefaultFigureColormap', cbar);
        end
         
        if strcmp(statistics.analysis_type, 'ROI')
            tmp_data(ROIs.idx(~h_adj_IFJa),idx_IFJa) = 0.04;
            tmp_data(ROIs.idx(~h_adj_IFJp),idx_IFJp) = 0.04;
            tmp_data(ROIs.idx(~h_adj_FEF),idx_FEF)   = 0.04;
            tmp_data(idx_IFJa,ROIs.idx(~h_adj_IFJa)) = 0.04;
            tmp_data(idx_IFJp,ROIs.idx(~h_adj_IFJp)) = 0.04;
            tmp_data(idx_FEF,ROIs.idx(~h_adj_FEF))   = 0.04;
        end
 
        
        %%
        conn_mat.cohspctrm = tmp_data;
        
        % output
        results.helper.idx_IFJa = idx_IFJa;
        results.helper.idx_IFJp = idx_IFJp;
        results.helper.idx_FEF = idx_FEF;
        results.helper.ROIs = ROIs;
        results.helper.brainordinate = brainordinate;
        results.conn_matrix = tmp_data;
        results.label = conn_mat.label;
        results.stat.h_adj_IFJa = h_adj_IFJa;
        results.stat.h_adj_IFJp = h_adj_IFJp;
        results.stat.h_adj_FEF = h_adj_FEF;
        
        results.perSubject.IFJa_12 = IFJa_12;
        results.perSubject.IFJp_12 = IFJp_12;
        results.perSubject.FEF_12 = FEF_12;
        
        results.perSubject.IFJa_21 = IFJa_21;
        results.perSubject.IFJp_21 = IFJp_21;
        results.perSubject.FEF_21 = FEF_21;
        
        %% Visualize using the helper function
        if flag.fsaverage == true
            % in order to run the function
            conn_mat.brainordinate.pos = conn_mat.brainordinate.pos*1000;
            
            %% save the figure
            if and(colorlim == true, flag.save_fig == true)

                if seed_target(1) == 'r'
                    seed = 'IFJa';
                    cimec.seed = seed;
                    pos2d_IFJa_R = [52.5323  -69.3638   71.0650; 52.5323   69.4183   71.0650];
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_IFJa_R = [52.5323  -69.3638   71.0650; 52.5323   69.4183   71.0650];
                    else   
                        pos2d_IFJa_R = [49.7454  -69.6590   70.1540; 49.7454   69.1230   70.1540];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_R);
                    helper_save_figure(cimec);
                    close all;
                    
                    seed = 'IFJp';
                    cimec.seed = seed;
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_IFJp_R = [45.6139  -69.3638   76.8706 ; 45.6139   69.4183   76.8706];
                    else   
                        pos2d_IFJp_R = [44.7673  -69.6590   76.5544 ; 44.7673   69.1230   76.5544];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_R);
                    helper_save_figure(cimec);
                    close all;
                    
                    seed = 'FEF';
                    cimec.seed = seed;
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_FEF_R = [30.1309  -69.3638  103.7783 ; 30.1309   69.4183  103.7783];
                    else   
                        pos2d_FEF_R = [30.8997  -69.6590  101.8004 ; 30.8997   69.1230  101.8004];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_R);
                    helper_save_figure(cimec);
                    close all;
                    
                elseif seed_target(1) == 'l'
                    seed = 'IFJa';
                    cimec.seed = seed;
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_IFJa_L = [50.5105   69.4183   69.5634 ; 53.1465  -69.3638   68.2860];
                    else
                        pos2d_IFJa_L = [47.6119   69.1230   68.3761 ; 47.6119  -69.6590   68.3761];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJa_L);
                    helper_save_figure(cimec);
                    close all;
                    
                    seed = 'IFJp';
                    cimec.seed = seed;
                    cimec.seed = seed;
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_IFJp_L = [43.6193   69.4183   80.0619 ; 43.6193  -69.3638   80.0619];
                    else
                        pos2d_IFJp_L = [40.8559   69.1230   76.9099 ; 40.8559  -69.6590   76.9099];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_IFJp_L);
                    helper_save_figure(cimec);
                    close all;
                    
                    seed = 'FEF';
                    cimec.seed = seed;
                    if strcmp(statistics.analysis_type, 'ROI')
                        pos2d_FEF_L = [30.6545   69.4183  105.7922 ; 30.6545  -69.3638  105.7922];
                    else   
                        pos2d_FEF_L = [33.7443   69.1230  101.0893 ; 33.7443  -69.6590  101.0893];
                    end
                    tutorial_nwa_connectivityviewer_save_figures(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours, pos2d_FEF_L);
                    helper_save_figure(cimec);
                    close all;
                end
                
            %% show the figure
            elseif colorlim == true
                tutorial_nwa_connectivityviewer(conn_mat, 'cohspctrm', colorlim_value, seed_target, flag.contours);
            elseif colorlim == false
                tutorial_nwa_connectivityviewer(conn_mat, 'cohspctrm', seed_target, flag.contours);
            end
            
            if contrast == false
                colormap('parula')
            end
        end
        
    end
end
