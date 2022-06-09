function helper_functionConn_circularGraph(results, fontSize, flag)

if flag.contrast == true
    figure;
    contrast_pair = 'FEF-IFJp';
    helper_circularGraph_contrast(results, fontSize, flag, contrast_pair)
    contrast_pair = 'FEF-IFJa';
    helper_circularGraph_contrast(results, fontSize, flag, contrast_pair)
else
    figure;
    seed = 'FEF';
    helper_circularGraph_notContrast(results, fontSize, flag, seed)
    figure;
    seed = 'IFJa';
    helper_circularGraph_notContrast(results, fontSize, flag, seed)
    figure;
    seed = 'IFJp';
    helper_circularGraph_notContrast(results, fontSize, flag, seed)
end

%%
    function helper_circularGraph_contrast(results, fontSize, flag, pair)
        
        tmp_mat = results.conn_matrix;
        idx_FEF = results.helper.idx_FEF;
        idx_IFJa = results.helper.idx_IFJa;
        idx_IFJp = results.helper.idx_IFJp;
        ROIs = results.helper.ROIs;
        analysis_type = results.stat.analysis_type;
        
        band_name = results.cimec.band_name ;
        conn_metric = results.cimec.conn_metric;
        seed_target = results.cimec.seed_target;
        alpha = char(string(results.cimec.statistics.alpha));
        correction = results.cimec.statistics.correction;
        duration = results.cimec.flag.duration;
        
        % https://sashamaps.net/docs/resources/20-colors/
        % https://html-color.codes/blue
        FEF_color_1 = [255,0,0]/255;
        FEF_color_2 = [255,160,122]/255;
        IFJa_color = [0,0,255]/255; %[0, 92, 171]/255;
        IFJp_color = [173,216,230]/255; %[70, 240, 240]/255;
        
        if strcmp(pair, 'FEF-IFJa')
            
            %% ROI analysis
            if strcmp(analysis_type, 'ROI')
                
                idxs = [idx_FEF, idx_IFJa, idx_IFJp, flip(ROIs.idx)];
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                
                % to make FEF results explicit
                data = data*-1;
                
                data(1,data(2,:)<0) = abs(data(2,data(2,:)<0));
                data(2,data(2,:)<0) = 0;
                data(:,1) = data(1,:);
                data(:,2) = data(2,:);
                data(:,3) = 0;
                data(3,:) = 0;
                
                label = results.label(idxs);
                
                if flag.collapseTargets == true
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(3:end-6);
                    end
                else
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(1:end-6);
                    end
                end
                
                % color
                % https://sashamaps.net/docs/resources/20-colors/
                myColorMap = repmat([1 1 1], length(idxs), 1);
                myColorMap(1,:) = FEF_color_1;
                myColorMap(2,:) = IFJa_color;
                
                %% Exploratory analysis
            elseif strcmp(analysis_type, 'exploratory')
                
                figure;
                
                idxs = ROIs.idx;
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                
                % to make FEF results explicit
                data = data*-1;
                
                data(1,data(2,:)<0) = abs(data(2,data(2,:)<0));
                data(2,data(2,:)<0) = 0;
                data(:,1) = data(1,:);
                data(:,2) = data(2,:);
                data(:,3) = 0;
                data(3,:) = 0;
                
                label = results.label(idxs);
                
                if flag.collapseTargets == true
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(3:end-6);
                    end
                else
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(1:end-6);
                    end
                end
                
                % color
                % https://sashamaps.net/docs/resources/20-colors/
                myColorMap = repmat([1 1 1], length(idxs), 1);
                myColorMap(idx_FEF-180,:) = FEF_color_1;
                myColorMap(idx_IFJa-180,:) = IFJa_color;
            end
            
            %% Circular graph for all together
            circularGraph(data', 'Colormap', myColorMap, 'Label', myLabel);
            title('All Frequency Bands')
            set(gcf, 'Position', get(0, 'Screensize'));
            
            fg1 = gca;
            fg2 = get(fg1, 'children');
            
            for ii = 1:length(fg2)
                try
                    tmp_fg = fg2(ii);
                    tmp_fg.FontSize = fontSize;
                catch
                    fprintf('Almost ready...\n', i);
                end
            end
            
            if flag.save_fig == true
                name = [ contrast_pair '_circularGraph_' conn_metric '_' band_name '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
                print(gcf, [flag.dir_fig name '.png'], '-dpng', '-r300');
                close all
            end
            
        elseif strcmp(pair, 'FEF-IFJp')
            %% ROI analysis
            if strcmp(analysis_type, 'ROI')
                
                idxs = [idx_FEF, idx_IFJa, idx_IFJp, flip(ROIs.idx)];
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                
                % to make FEF results explicit
                data = data*-1;
                
                data(1,data(3,:)<0) = abs(data(3,data(3,:)<0));
                data(3,data(3,:)<0) = 0;
                data(:,1) = data(1,:);
                data(:,3) = data(3,:);
                data(:,2) = 0;
                data(2,:) = 0;
                
                label = results.label(idxs);
                
                if flag.collapseTargets == true
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(3:end-6);
                    end
                else
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(1:end-6);
                    end
                end
                
                % color
                % https://sashamaps.net/docs/resources/20-colors/
                myColorMap = repmat([1 1 1], length(idxs), 1);
                myColorMap(1,:) = FEF_color_2;
                myColorMap(3,:) = IFJp_color;
                
                %% Exploratory analysis
            elseif strcmp(analysis_type, 'exploratory')
                
                idxs = ROIs.idx;
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                
                % to make FEF results explicit
                data = data*-1;
                
                data(1,data(3,:)<0) = abs(data(3,data(3,:)<0));
                data(3,data(3,:)<0) = 0;
                data(:,1) = data(1,:);
                data(:,3) = data(3,:);
                data(:,2) = 0;
                data(2,:) = 0;
                
                label = results.label(idxs);
                
                if flag.collapseTargets == true
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(3:end-6);
                    end
                else
                    % labels
                    myLabel = cell(length(label));
                    for i = 1:length(label)
                        myLabel{i} = strrep(label{i},'_','-');
                        myLabel{i} = myLabel{i}(1:end-6);
                    end
                end
                
                % color
                % https://sashamaps.net/docs/resources/20-colors/
                myColorMap = repmat([1 1 1], length(idxs), 1);
                myColorMap(idx_FEF-180,:) = FEF_color_2;
                myColorMap(idx_IFJp-180,:) = IFJp_color;
            end
            %% Circular graph for all together
            circularGraph(data', 'Colormap', myColorMap, 'Label', myLabel);
            title('All Frequency Bands')
            set(gcf, 'Position', get(0, 'Screensize'));
            
            fg1 = gca;
            fg2 = get(fg1, 'children');
            
            for ii = 1:length(fg2)
                try
                    tmp_fg = fg2(ii);
                    tmp_fg.FontSize = fontSize;
                catch
                    fprintf('Almost ready...\n', i);
                end
            end
            
            if flag.save_fig == true
                name = [contrast_pair '_circularGraph_' conn_metric '_' band_name '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
                print(gcf, [flag.dir_fig name '.png'], '-dpng', '-r300');
                close all
            end
            
            
        end
    end

    function helper_circularGraph_notContrast(results, fontSize, flag, seed)
        
        tmp_mat = results.conn_matrix;
        idx_FEF = results.helper.idx_FEF;
        idx_IFJa = results.helper.idx_IFJa;
        idx_IFJp = results.helper.idx_IFJp;
        ROIs = results.helper.ROIs;
        analysis_type = results.stat.analysis_type;
        
        band_name = results.cimec.band_name ;
        conn_metric = results.cimec.conn_metric;
        seed_target = results.cimec.seed_target;
        alpha = char(string(results.cimec.statistics.alpha));
        correction = results.cimec.statistics.correction;
        duration = results.cimec.flag.duration;
        
        % https://sashamaps.net/docs/resources/20-colors/
        % https://html-color.codes/blue
        FEF_color = [255,0,0]/255;
        IFJa_color = [0,0,255]/255; %[0, 92, 171]/255;
        IFJp_color = [173,216,230]/255; %[70, 240, 240]/255;
        
        %% ROI analysis
        if strcmp(analysis_type, 'ROI')
                        
            if strcmp(seed, 'FEF')
                idxs = [idx_FEF, flip(ROIs.idx)];
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                data(:,1) = data(1,:);
            elseif strcmp(seed, 'IFJa')
                idxs = [idx_IFJa, flip(ROIs.idx)];
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                data(:,1) = data(1,:);
            elseif strcmp(seed, 'IFJp')
                idxs = [idx_IFJp, flip(ROIs.idx)];
                data = tmp_mat(idxs, idxs);
                data(isnan(data)) = 0;
                data(:,1) = data(1,:);
            end
            
            label = results.label(idxs);
            
            if flag.collapseTargets == true
                % labels
                myLabel = cell(length(label));
                for i = 1:length(label)
                    myLabel{i} = strrep(label{i},'_','-');
                    myLabel{i} = myLabel{i}(3:end-6);
                end
            else
                % labels
                myLabel = cell(length(label));
                for i = 1:length(label)
                    myLabel{i} = strrep(label{i},'_','-');
                    myLabel{i} = myLabel{i}(1:end-6);
                end
            end
            
            % color
            % https://sashamaps.net/docs/resources/20-colors/
            myColorMap = repmat([1 1 1], length(idxs), 1);
            if strcmp(seed, 'FEF')
                myColorMap(1,:) = FEF_color;
            elseif strcmp(seed, 'IFJa')
                myColorMap(1,:) = IFJa_color;
            elseif strcmp(seed, 'IFJp')
                myColorMap(1,:) = IFJp_color;
            end
        end 

        %% Circular graph for all together
        circularGraph(data', 'Colormap', myColorMap, 'Label', myLabel);
        title('All Frequency Bands')
        set(gcf, 'Position', get(0, 'Screensize'));
        
        fg1 = gca;
        fg2 = get(fg1, 'children');
        
        for ii = 1:length(fg2)
            try
                tmp_fg = fg2(ii);
                tmp_fg.FontSize = fontSize;
            catch
                fprintf('Almost ready...\n', i);
            end
        end
        
        if flag.save_fig == true
            name = [ seed '_circularGraph_' conn_metric '_' band_name '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
            print(gcf, [flag.dir_fig name '.png'], '-dpng', '-r300');
            close all
        end
    end
end
