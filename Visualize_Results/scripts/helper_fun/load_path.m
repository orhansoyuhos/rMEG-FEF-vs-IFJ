%% Please set the working directory first before running the scripts.

function path = load_path(directory)

workingDir = 'C:\Users\ASUS\Desktop\Thesis Project\';
fieldTrip_dir = 'C:\Users\ASUS\Desktop\MATLAB\fieldtrip-20210411';
save_fig = 'C:\Users\ASUS\Desktop\Thesis Project\Main Figures\';

% restoredefaultpath % if needed
addpath(fieldTrip_dir);
ft_defaults

if strcmp(directory, 'inputfile_fs')
    path = fullfile(workingDir, 'Visualize results', 'fsaverage'); 
elseif strcmp(directory, 'Outputs_Colclough')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_Colclough');
elseif strcmp(directory, 'Outputs_2s_55subjs_icoh')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_icoh');
elseif strcmp(directory, 'Outputs_2s_55subjs_wPLI')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_wPLI');
elseif strcmp(directory, 'Outputs_2s_55subjs_dwPLI')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_dwPLI');
elseif strcmp(directory, 'Outputs_2s_55subjs_powcorrOrtho')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_powcorrOrtho');
elseif strcmp(directory, 'Outputs_2s_55subjs_granger')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_granger');
elseif strcmp(directory, 'Outputs_2s_55subjs_pdc')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_2s_55subjs_pdc');

elseif strcmp(directory, 'Outputs_5s_55subjs_icoh')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_5s_55subjs_icoh');
elseif strcmp(directory, 'Outputs_5s_55subjs_wPLI')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_5s_55subjs_wPLI');   
elseif strcmp(directory, 'Outputs_5s_55subjs_powcorrOrtho')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_5s_55subjs_powcorrOrtho');    
    
elseif strcmp(directory, 'Outputs_10s_55subjs_icoh')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_10s_55subjs_icoh');
elseif strcmp(directory, 'Outputs_10s_55subjs_powcorrOrtho')
    path = fullfile(workingDir, 'Visualize results', 'results', 'Outputs_10s_55subjs_powcorrOrtho');  
  
elseif strcmp(directory, 'workingDir') 
    path = workingDir;
elseif strcmp(directory, 'helper_fun') 
    path = fullfile(workingDir, 'Visualize results', 'scripts', 'helper_fun');
elseif strcmp(directory, 'save_fig')
    path = save_fig;
    
    
end    
    
    