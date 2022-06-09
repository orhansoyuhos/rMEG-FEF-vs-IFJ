function helper_boxplot(seed1, seed2, seeds, ROIs)

%% Run with MATLAB 2021b and later.

n_seeds = length(seeds);

% labels
label = ROIs.name;
for i = 1:length(label)
    myLabel{i} = strrep(label{i},'_','-');
    myLabel{i} = myLabel{i}(1:end-6);
end

%% dorsal stream
idx_dorsal = 1:17;
dorsal = myLabel(idx_dorsal);
n_targets = length(dorsal);
ROI_FEF = seed1(ROIs.idx(idx_dorsal),:)';
ROI_IFJ = seed2(ROIs.idx(idx_dorsal),:)';

tmp_Seeds = repmat(repmat(seeds,1,1), 1, 55*n_targets)';
Seeds = tmp_Seeds(:);
Targets = repmat(repmat(dorsal',55,1), n_seeds, 1);
Subjects = repmat(repmat(1:55,1,n_targets)', n_seeds, 1);
targetsOrder = dorsal;
Conn_values = [ROI_FEF(:);ROI_IFJ(:)];
results_tab = table(Seeds,Targets,Conn_values);
results_tab.Targets = categorical(results_tab.Targets, targetsOrder);
figure; h_one = boxchart(results_tab.Targets, results_tab.Conn_values, 'GroupByColor', results_tab.Seeds);
legend
title('Dorsam Visual Stream', 'fontweight','bold', 'FontSize', 15)
set(gcf, 'Position', get(0, 'Screensize'));

clear results_tab

h_one(1).SeriesIndex = 2;
h_one(2).SeriesIndex = 1;
h_one(1).WhiskerLineStyle = ':';
h_one(2).WhiskerLineStyle = ':';
h_one(1).MarkerSize = 6;
h_one(2).MarkerSize = 6;


xlabel('ROIs', 'fontweight','bold', 'FontSize', 12)
ylabel('Connectivity Values', 'fontweight','bold', 'FontSize', 12)
ax=gca;
ax.FontSize = 12;

print(gcf,'C:\Users\ASUS\Desktop\dorsal_boxchart.png','-dpng','-r300');


%% ventral stream
idx_ventral = 18:33;
ventral = myLabel(idx_ventral);
n_targets = length(ventral);
ROI_FEF = seed1(ROIs.idx(idx_ventral),:)';
ROI_IFJ = seed2(ROIs.idx(idx_ventral),:)';

tmp_Seeds = repmat(repmat(seeds,1,1), 1, 55*n_targets)';
Seeds = tmp_Seeds(:);
Targets = repmat(repmat(ventral',55,1), n_seeds, 1);
Subjects = repmat(repmat(1:55,1,n_targets)', n_seeds, 1);
targetsOrder = ventral;
Conn_values = [ROI_FEF(:);ROI_IFJ(:)];
results_tab = table(Seeds,Targets,Conn_values);
results_tab.Targets = categorical(results_tab.Targets, targetsOrder);
figure; h_two = boxchart(results_tab.Targets, results_tab.Conn_values, 'GroupByColor', results_tab.Seeds);
legend
title('Ventral Visual Stream', 'fontweight','bold', 'FontSize', 15)
set(gcf, 'Position', get(0, 'Screensize'));

clear results_tab

h_two(1).SeriesIndex = 2;
h_two(2).SeriesIndex = 1;
h_two(1).WhiskerLineStyle = ':';
h_two(2).WhiskerLineStyle = ':';
h_two(1).MarkerSize = 6;
h_two(2).MarkerSize = 6;


xlabel('ROIs', 'fontweight','bold', 'FontSize', 12)
ylabel('Connectivity Values', 'fontweight','bold', 'FontSize', 12)
ax=gca;
ax.FontSize = 12;


print(gcf,'C:\Users\ASUS\Desktop\ventral_boxchart.png','-dpng','-r300');


end