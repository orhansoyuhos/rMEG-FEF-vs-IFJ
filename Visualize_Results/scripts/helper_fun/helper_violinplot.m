function helper_violinplot(outputs, conn_metric, band_name, seed_target,  duration)

%% Cite: https://github.com/bastibe/Violinplot-Matlab
idx_FEF_R = 242;
idx_IFJa_R = 251;
idx_IFJp_R = 252;

idx_FEF_L = 62;
idx_IFJa_L = 71;
idx_IFJp_L = 72;

if strcmp(band_name, 'delta')
    color = [240 50 230]/255;
elseif strcmp(band_name, 'theta')
    color = [255 225 25]/255;
elseif strcmp(band_name, 'alpha')
    color = [60 180 75]/255;
elseif strcmp(band_name, 'beta')
    color = [245 130 48]/255;
elseif strcmp(band_name, 'gamma')
    color = [128 128 128]/255;
end

if strcmp(band_name, 'gamma')
    yvalues = [0 0.24];
else
    yvalues = [0 0.14];
end

%% 'right-right'
if strcmp(seed_target, 'right-right')
    
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-FEF and R-IFJa'];
    legend = {'R-FEF to R-IFJa' 'R-IFJa to R-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJa_R,:)'; seed_21(idx_IFJa_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-FEF and R-IFJp'];
    legend = {'R-FEF to R-IFJp' 'R-IFJp to R-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_R,:)'; seed_21(idx_IFJp_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    %%
    seed_12 = outputs.perSubject.IFJa_12;
    seed_21 = outputs.perSubject.IFJa_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-IFJa and R-IFJp'];
    legend = {'R-IFJa to R-IFJp' 'R-IFJp to R-IFJa'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_R,:)'; seed_21(idx_IFJp_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    
    %% 'left-left'
elseif strcmp(seed_target, 'left-left')
    
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-FEF and L-IFJa'];
    legend = {'L-FEF to L-IFJa' 'L-IFJa to L-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJa_L,:)'; seed_21(idx_IFJa_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-FEF and L-IFJp'];
    legend = {'L-FEF to L-IFJp' 'L-IFJp to L-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_L,:)'; seed_21(idx_IFJp_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    %%
    seed_12 = outputs.perSubject.IFJa_12;
    seed_21 = outputs.perSubject.IFJa_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-IFJa and L-IFJp'];
    legend = {'L-IFJa to L-IFJp' 'L-IFJp to L-IFJa'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_L,:)'; seed_21(idx_IFJp_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    
    
    %% 'right-left'
elseif strcmp(seed_target, 'right-left')
    
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-FEF and L-IFJa'];
    legend = {'R-FEF to L-IFJa' 'L-IFJa to R-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJa_L,:)'; seed_21(idx_IFJa_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-FEF and L-IFJp'];
    legend = {'R-FEF to L-IFJp' 'L-IFJp to R-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_L,:)'; seed_21(idx_IFJp_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.IFJa_12;
    seed_21 = outputs.perSubject.IFJa_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'R-IFJa and L-IFJp'];
    legend = {'R-IFJa to L-IFJp' 'L-IFJp to R-IFJa'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_L,:)'; seed_21(idx_IFJp_L,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %% 'left-right'
elseif strcmp(seed_target, 'left-right')
    
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-FEF and R-IFJa'];
    legend = {'L-FEF to R-IFJa' 'R-IFJa to L-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJa_R,:)'; seed_21(idx_IFJa_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-FEF and R-IFJp'];
    legend = {'L-FEF to R-IFJp' 'R-IFJp to L-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_R,:)'; seed_21(idx_IFJp_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    %%
    seed_12 = outputs.perSubject.IFJa_12;
    seed_21 = outputs.perSubject.IFJa_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-IFJa and R-IFJp'];
    legend = {'L-IFJa to R-IFJp' 'R-IFJp to L-IFJa'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_R,:)'; seed_21(idx_IFJp_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    %% same regions L and R
    %%
    seed_12 = outputs.perSubject.FEF_12;
    seed_21 = outputs.perSubject.FEF_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-FEF and R-FEF'];
    legend = {'L-FEF to R-FEF' 'R-FEF to L-FEF'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_FEF_R,:)'; seed_21(idx_FEF_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
    %%
    seed_12 = outputs.perSubject.IFJa_12;
    seed_21 = outputs.perSubject.IFJa_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-IFJa and R-IFJa'];
    legend = {'L-IFJa to R-IFJa' 'R-IFJa to L-IFJa'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJa_R,:)'; seed_21(idx_IFJa_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    %%
    seed_12 = outputs.perSubject.IFJp_12;
    seed_21 = outputs.perSubject.IFJp_21;
    
    name = [conn_metric '_' band_name '_' seed_target '_' duration '_' 'L-IFJp and R-IFJp'];
    legend = {'L-IFJp to R-IFJp' 'R-IFJp to L-IFJp'};
    labels = repmat(repmat(legend,55,1), 1, 1);
    labels = labels(:);
    
    values = [seed_12(idx_IFJp_R,:)'; seed_21(idx_IFJp_R,:)'];
    
    figure; vs = violinplot(values, labels, 'ShowMean', true);
    
    ylabel('Partial Directed Coherence'); ylim(yvalues);
    title('Granger Causal Influences')
    vs(1).ViolinColor = color;
    vs(2).ViolinColor = color;
    
    print(gcf,[load_path('save_fig') name '.png'],'-dpng','-r300');
    close all;
    
    
end


end